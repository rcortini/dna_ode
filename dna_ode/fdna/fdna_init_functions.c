/* dna_ode: simulating DNA using the ODE physics engine
 *
 * Copyright (C) 2014, 2015  Ruggero Cortini

 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "config_read.h"
#include "dna_ode_core.h"
#include "dna.h"
#include "trajectory_tools.h"
#include "fdna.h"



/* 
 * INITIALIZATION FUNCTIONS 
 */



/* loads the initial configuration from a file. This file must be constructed as an ordinary trajectory file,
 * consisting in a single frame. You can extract the last frame from an existing trajectory using "trajectory_analyze frame_dump" */
void fdna_load_initial_configuration_from_file (simcontext *simcon, fdna_simcontext *fdna_simcon) {
  unsigned int i;
  const unsigned int nbodies = fdna_simcon->dna->nsegments;
  int frame_load_ret_code;
  trajectory_frame frame;

  /* open the initial configuration file. Note that this function is invoked only after it was checked that this file
   * exists */
  simcon->initial_configuration_file = my_fopen (simcon->simparams->initial_configuration_file_name, "r");

  /* use "load_trajectory_frame" to load all the info from the existing trajectory frame */
  frame_load_ret_code = load_trajectory_frame (nbodies, simcon->initial_configuration_file, &frame, NULL, NULL);

  /* close the initial configuration file */
  fclose (simcon->initial_configuration_file);

  /* check if success */
  if (frame_load_ret_code==LOAD_TRAJ_FRAME_FAIL) {
    err_message ("Error initializing configuration from %s\n", simcon->simparams->initial_configuration_file_name);
    exit (EXIT_FAILURE);
  }

  /* now assign positions and velocities to all bodies */
  for (i=0; i<nbodies; i++) {
    body_set_position (fdna_simcon->dna->seg [i]->b_id, frame.body_i_state [i].r);
    body_set_velocity (fdna_simcon->dna->seg [i]->b_id, frame.body_i_state [i].vel, frame.body_i_state [i].omega);
    body_set_axes (fdna_simcon->dna->seg [i]->b_id, frame.body_i_state [i].u, frame.body_i_state [i].v, frame.body_i_state [i].t);
  }

  /* set simulation initial clock */
  simcon->sim_clock_0 = frame.sim_clock;
}



/* initializes the simulation statistics data structures */
void fdna_init_simulation_context (simcontext *simcon, fdna_simcontext *fdna_simcon) {
  unsigned int i;
  simparameters *simparams = simcon->simparams;

  /* simulation data */
  fdna_simcon->Et = 0.;
  fdna_simcon->Kt = 0.;
  fdna_simcon->Ut = 0.;
  fdna_simcon->ext_Ut = 0.;

  /* SIMULATED OBJECTS CREATION */

  /* creates the DNA */
  fdna_simcon->dna = fdna_create (simcon);
  if (simcon->mechanical_model.f == body_twist_and_torque_nicked_harmonic ||
      simcon->mechanical_model.f == body_twist_and_torque_nicked_kinkable ||
      simcon->mechanical_model.f == body_twist_and_torque_nicked_harmonic4)
    fdna_simcon->dna->nicked = 1;
  else
    fdna_simcon->dna->nicked = 0;

  /* THERMOSTAT INIT */

  /* mE is used to rescale the velocities in the final step of the Langevin-Euler dynamics */
  PAR (mE) = 3*PAR(dna_nsegments)/2.*PAR (beta);

  /* init the body list, which is used by the thermostat */
  simcon->body_list = (sim_body *) malloc (fdna_simcon->dna->nsegments * (sizeof (sim_body)));
  for (i=0; i<fdna_simcon->dna->nsegments; i++) {
    char *name = (char *) malloc (STR_BUFF_SIZE * sizeof (char));
    sprintf (name, "DNA segment %d", i);
    sim_body_init (&simcon->body_list [i], i, name, PAR (sigma_langevin), PAR (ode_step), &fdna_simcon->dna->seg [i]->b_id);
  }

  /* init thermostat */
  init_langevin_thermostat (&simcon->my_thermostat,
      PAR (thermostat_type),
      &fdna_simcon->Kt,
      PAR (sigma_langevin),
      PAR (ode_step),
      PAR (beta),
      PAR (mE));

  /* JOINT SOLUTION WORKSPACE INIT */
  fdna_simcon->ws = fdna_joint_solution_workspace_alloc (fdna_simcon->dna->nsegments);
  fdna_joint_solution_workspace_init (fdna_simcon);

  /* SIM INFO */

  /* print further information to the terminal */
  message ("%i geoms, %ibp\n", dSpaceGetNumGeoms (simcon->space_id), simcon->simparams->dna_nbasepairs_in_segment);

  /* thermostat parameters */
  message ("Thermostat type = %s\n", simcon->simparams->thermostat_type);
  message ("Langevin sigma = %.3f (%.3e s^-1)\n",
      simparams->sigma_langevin,
      simparams->sigma_langevin / simparams->scale_t);
  message ("<E>:%f\n", simcon->simparams->mE);
  print_dna_info (simcon, fdna_simcon->dna->seg [0], &simcon->body_list [0]);
}



/* initializes the simulated objects' initial positions and orientations */
void fdna_init_molecule_positions_orientations (simcontext * simcon, fdna_simcontext * fdna_simcon) {
  dna_segment **dna = fdna_simcon->dna->seg;
  t_vector t0;
 
  /* get tangent to the first cylinder */
  body_axis (t0, fdna_simcon->dna->seg [0]->b_id, 2);

  /* the rest of the DNA molecule is randomly assigned as a worm-like chain segment with
   * zero force. Try repeatedly until the DNA molecule self-intersects. */
  do {
    simcon->numc = 0;
    dna_wlc_create (fdna_simcon->dna->nsegments, simcon->simparams->beta, 1., simcon->simparams->gb, t0, dna);
    dSpaceCollide (simcon->space_id, simcon, &collision_callback);
  } while (simcon->numc != 0);
}



/* initializes the simulated objects' initial velocities */
void fdna_init_molecule_velocities (simcontext * simcon, fdna_simcontext * fdna_simcon) {
  /* shake the DNA */
  dna_shake (fdna_simcon->dna, simcon->body_list, 0, &simcon->mechanical_model);

  /* if user specified that we need to calculate the pair interactions, do so */
  if (simcon->simparams->forces_on) fdna_simcon->Ut = force_loop (&simcon->U, fdna_simcon->dna);

  /* perform collision check, and perform the first step */
  dSpaceCollide (simcon->space_id, simcon, &collision_callback);
  dWorldStep (simcon->world_id, simcon->simparams->ode_step);
  if (simcon->contactjoints)
    dJointGroupEmpty (simcon->contactjoints);
}

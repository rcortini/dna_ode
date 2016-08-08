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
#include "trajectory_tools.h"
#include "dna.h"
#include "cdna.h"



/* 
 * INITIALIZATION FUNCTIONS 
 */



/* loads the initial configuration from a file. This file must be constructed as an ordinary trajectory file,
 * consisting in a single frame. You can extract the last frame from an existing trajectory using "trajectory_analyze frame_dump" */
void cdna_load_initial_configuration_from_file (simcontext *simcon, cdna_simcontext *cdna_simcon) {
  unsigned int i;
  const unsigned int nbodies = cdna_simcon->dna->nsegments;
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
    body_set_position (cdna_simcon->dna->seg [i]->b_id, frame.body_i_state [i].r);
    body_set_velocity (cdna_simcon->dna->seg [i]->b_id, frame.body_i_state [i].vel, frame.body_i_state [i].omega);
    body_set_axes (cdna_simcon->dna->seg [i]->b_id, frame.body_i_state [i].u, frame.body_i_state [i].v, frame.body_i_state [i].t);
  }

  /* set simulation initial clock */
  simcon->sim_clock_0 = frame.sim_clock;
}



void cdna_init_simulation_context (simcontext *simcon, cdna_simcontext *cdna_simcon) {
  unsigned int i, N;
  simparameters *simparams = simcon->simparams;
  config_t * simconfig = PAR (simconfig);

  /* get system-specific parameters */
  cdna_simcon->sigma = mycfg_read_double (simconfig, "sigma");

  /* sets to zero a bunch of stuff */
  cdna_simcon->Wr = 0.;
  cdna_simcon->dTw = 0.;
  cdna_simcon->Lk = 0.;
  cdna_simcon->Et = 0.;
  cdna_simcon->Kt = 0.;
  cdna_simcon->Ut = 0.;
  cdna_simcon->ext_Ut = 0.;

  /* SIMULATED OBJECTS CREATION */

  /* creates the DNA */
  cdna_simcon->dna = cdna_create (simcon);
  if (simcon->mechanical_model.f == body_twist_and_torque_nicked_harmonic ||
      simcon->mechanical_model.f == body_twist_and_torque_nicked_kinkable ||
      simcon->mechanical_model.f == body_twist_and_torque_nicked_harmonic4)
    cdna_simcon->dna->nicked = 1;
  else
    cdna_simcon->dna->nicked = 0;

  cdna_simcon->delta_Lk0 = cdna_simcon->sigma * cdna_simcon->dna->Lk0;
  cdna_simcon->nturns = round (cdna_simcon->delta_Lk0);
  /* we calculated the number of turns that we need to impose to the DNA molecule to have the user-given sigma.
   * In reality, we had to round n to the closest integer, so sigma will be slightly different */
  cdna_simcon->sigma_prime = (t_real) cdna_simcon->nturns / cdna_simcon->dna->Lk0;

  /* LANGEVIN DYNAMICS PARAMETERS */

  /* mE is used to rescale the velocities in the final step of the Langevin-Euler dynamics. TODO: recalculate the number of degrees of freedom */
  PAR (mE) = 3*PAR(dna_nsegments)/2.*PAR (beta);

  /* assign the thermostat type and parameters */
  N = cdna_simcon->dna->nsegments;

  /* init the body list, which is used by the thermostat */
  simcon->body_list = (sim_body *) malloc (cdna_simcon->dna->nsegments * (sizeof (sim_body)));
  for (i=0; i<N; i++) {
    char *name = (char *) malloc (STR_BUFF_SIZE * sizeof (char));
    sprintf (name, "DNA segment %d", i);
    sim_body_init (&simcon->body_list [i], i, name, PAR (sigma_langevin), PAR (ode_step), &cdna_simcon->dna->seg [i]->b_id);
  }

  /* init thermostat */
  init_langevin_thermostat (&simcon->my_thermostat,
      PAR (thermostat_type),
      &cdna_simcon->Kt,
      PAR (sigma_langevin),
      PAR (ode_step),
      PAR (beta),
      PAR (mE));

  /* JOINT SOLUTION WORKSPACE INIT */
  cdna_simcon->ws = ccdna_joint_solution_workspace_alloc (cdna_simcon->dna->nsegments);
  ccdna_joint_solution_workspace_init (cdna_simcon);

  /* SIM INFO */
  message ("Thermostat type = %s\n", simcon->simparams->thermostat_type);
  message ("Langevin sigma = %.3f (%.3e s^-1)\n",
      simparams->sigma_langevin,
      simparams->sigma_langevin / simparams->scale_t);
  message ("<E>:%f\n", simcon->simparams->mE);
  print_dna_info (simcon, cdna_simcon->dna->seg [0], &simcon->body_list [0]);

  message ("Lk0: %f delta_Lk0: %f\n", cdna_simcon->dna->Lk0, cdna_simcon->delta_Lk0);
  message ("sigma = %f, sigma_prime = %f\n", cdna_simcon->sigma, cdna_simcon->sigma_prime);
}



/* initializes the simulated objects' initial positions and orientations */
void cdna_init_molecule_positions_orientations (simcontext * simcon, cdna_simcontext *cdna_simcon) {
  unsigned int i;
  t_vector previous_position, previous_u, previous_v, previous_t;
  t_vector next_position, next_u, next_v, next_t;
  t_vector joint_position, v0;
  (void) simcon;

  /* 
   * these constants are used to calculate the rotation angles that we need to put to impose that
   * 1. the DNA is circular (alpha)
   * 2. the DNA has a constant, imposed twist (twist_angle)
   */
  const t_real alpha = 2*M_PI/((t_real) cdna_simcon->dna->nsegments);
  const t_real twist_angle = 2*M_PI*cdna_simcon->nturns/cdna_simcon->dna->nsegments;

  /* check if the value of the twist angle is greater than pi, in which case we must exit */
  if (twist_angle > M_PI) {
    err_message ("Sigma value too high!\n");
    exit (EXIT_FAILURE);
  }

  /* assign the rotation axis. v0 is used so that every DNA segment is rotated around the y axis, so initially
   * the molecule's centerline is on a plane, and it has zero writhe */
  vector_set (v0, 0., 0., 1.);

  /* assign the position of the first DNA segment */
  /* cycle to assign position and orientations of all the DNA segments, except the first one */
  body_position (next_position, cdna_simcon->dna->seg [0]->b_id);
  body_axis (next_u, cdna_simcon->dna->seg [0]->b_id, 0);
  body_axis (next_v, cdna_simcon->dna->seg [0]->b_id, 1);
  body_axis (next_t, cdna_simcon->dna->seg [0]->b_id, 2);
  for (i=1; i<cdna_simcon->dna->nsegments; i++) {
    t_vector buff_u, buff_v;

    /* get the previous segment's position and orientation */
    dCopyVector3 (previous_position, next_position);
    dCopyVector3 (previous_u, next_u);
    dCopyVector3 (previous_v, next_v);
    dCopyVector3 (previous_t, next_t);

    /* now rotate the axes of the previous segment by an angle given by alpha */
    vector_rotate_around_axis (buff_u, v0, alpha, previous_u);
    vector_rotate_around_axis (buff_v, v0, alpha, previous_v);
    vector_rotate_around_axis (next_t, v0, alpha, previous_t);

    /* adjust the DNA intrinsic twist according to sigma */
    vector_rotate_around_axis (next_u, next_t, twist_angle, buff_u);
    vector_rotate_around_axis (next_v, next_t, twist_angle, buff_v);

    /* now assign the new position of the segment */
    segment_next_joint (joint_position, cdna_simcon->dna->seg [i-1]);
    vector_lincomb (next_position, joint_position, next_t, 1., .5*cdna_simcon->dna->seg [i]->length);
    body_set_position (cdna_simcon->dna->seg [i]->b_id, next_position);

    /* and finally its orientation */
    body_set_axes (cdna_simcon->dna->seg [i]->b_id, next_u, next_v, next_t);
  }
}



/* initializes the simulated objects' initial velocities */
void cdna_init_molecule_velocities (simcontext * simcon, cdna_simcontext *cdna_simcon) {
  /* shake the DNA */
  dna_shake (cdna_simcon->dna, simcon->body_list, 0, &simcon->mechanical_model);

  /* if user specified that we need to calculate the pair interactions, do so */
  if (simcon->simparams->forces_on) cdna_simcon->Ut = force_loop (&simcon->U, cdna_simcon->dna);

  /* perform collision check, and perform the first step */
  dSpaceCollide (simcon->space_id, simcon, &collision_callback);
  dWorldStep (simcon->world_id, simcon->simparams->ode_step);
  if (simcon->contactjoints)
    dJointGroupEmpty (simcon->contactjoints);
}

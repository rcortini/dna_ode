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

#include "dna_ode_core.h"
#include "config_read.h"
#include "trajectory_tools.h"
#include "dna.h"
#include "writhe.h"
#include "sdna.h"



/* 
 * INITIALIZATION FUNCTIONS 
 */



/* loads the initial configuration from a file. This file must be constructed as an ordinary trajectory file,
 * consisting in a single frame. You can extract the last frame from an existing trajectory using "trajectory_analyze frame_dump" */
void sdna_load_initial_configuration_from_file (simcontext *simcon, sdna_simcontext *sdna_simcon) {
  unsigned int i;
  const unsigned int nbodies = sdna_simcon->dna->nsegments+1;
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
  /* ...DNA... */
  for (i=0; i<nbodies-1; i++) {
    body_set_position (sdna_simcon->dna->seg [i]->b_id, frame.body_i_state [i].r);
    body_set_velocity (sdna_simcon->dna->seg [i]->b_id, frame.body_i_state [i].vel, frame.body_i_state [i].omega);
    body_set_axes (sdna_simcon->dna->seg [i]->b_id, frame.body_i_state [i].u, frame.body_i_state [i].v, frame.body_i_state [i].t);
  }

  /* ...and bead */
  body_set_position (sdna_simcon->bead->b_id, frame.body_i_state [nbodies-1].r);
  body_set_velocity (sdna_simcon->bead->b_id, frame.body_i_state [nbodies-1].vel, frame.body_i_state [nbodies-1].omega);
  body_set_axes (sdna_simcon->bead->b_id, frame.body_i_state [nbodies-1].u, frame.body_i_state [nbodies-1].v, frame.body_i_state [nbodies-1].t);

  /* set simulation initial clock */
  simcon->sim_clock_0 = frame.sim_clock;
}


void sdna_init_simulation_context (simcontext *simcon, sdna_simcontext *sdna_simcon) {
  unsigned int i, N;
  t_vector bead_position, bead_R;
  config_t * simconfig = PAR (simconfig);

  /* initial configuration has loop? */
  sdna_simcon->initial_config_has_loop = mycfg_read_int (simconfig, "initial_config_has_loop");
  if (sdna_simcon->initial_config_has_loop) {
    sdna_simcon->loop_nsegments = mycfg_read_int (simconfig, "loop_nsegments");
    sdna_simcon->loop_nturns = mycfg_read_double (simconfig, "loop_nturns");
    sdna_simcon->loop_pitch = mycfg_read_double (simconfig, "loop_pitch");
  }

  /* get the surface parameters */
  sdna_simcon->surface_height_to_reff_ratio = mycfg_read_double (simconfig, "surface_height_to_reff_ratio");

  /* assign parameters of the surface */
  sdna_simcon->surface_height = -1.* sdna_simcon->surface_height_to_reff_ratio * PAR (dna_effective_rad);

  /* SIMULATED OBJECTS CREATION */

  /* creates the surface */
  sdna_simcon->surface = surface_create (simcon->space_id, 0., 0., 1., sdna_simcon->surface_height);

  /* creates the DNA */
  sdna_simcon->dna = sdna_create (simcon);
  N = sdna_simcon->dna->nsegments;
  if (simcon->mechanical_model.f == body_twist_and_torque_nicked_harmonic ||
      simcon->mechanical_model.f == body_twist_and_torque_nicked_kinkable ||
      simcon->mechanical_model.f == body_twist_and_torque_nicked_harmonic4)
    sdna_simcon->dna->nicked = 1;
  else
    sdna_simcon->dna->nicked = 0;

  /* creates the bead */
  sdna_simcon->bead = (bead_ode *) malloc (sizeof (bead_ode));
  bead_init (sdna_simcon->bead,
      N,
      simcon->odepdb_file,
      simcon->world_id,
      simcon->space_id,
      PAR (finite_rotation),
      PAR (gyroscopic_mode),
      PAR (scale_l),
      PAR (scale_e),
      PAR (scale_f),
      simconfig);

  /* sets the bead position */
  body_point_position (bead_position, sdna_simcon->dna->seg [sdna_simcon->dna->nsegments-1]->b_id,
      0., 0., 0.5 * simcon->simparams->dna_segment_length);
  vector_set (bead_R, 0., 0., sdna_simcon->bead->radius);
  vector_inplace_sum (bead_position, bead_R);
  body_set_position (sdna_simcon->bead->b_id, bead_position);

  /* now that the surface and the DNA are created, we can assign the joints */
  sdna_joints_create (simcon, sdna_simcon, sdna_simcon->dna);

  /* assign bead mode */
  if (sdna_simcon->bead->mode==FIXED_EXT) {
    t_vector first_seg_pos;
    joint_anchor (first_seg_pos, sdna_simcon->dna->joint [0]->j_id);

    /* assign trap parameters */
    sdna_simcon->bead->trap_x0 = first_seg_pos[0];
    sdna_simcon->bead->trap_y0 = first_seg_pos[1];
    sdna_simcon->bead->trap_z0 = first_seg_pos[2] + PAR (dna_length)*sdna_simcon->bead->rho + sdna_simcon->bead->radius;
  }

  /* kinetic energy */
  sdna_simcon->Kt = 0.;
  for (i=0; i<N; i++) sdna_simcon->Kt += body_total_kinetic_energy (sdna_simcon->dna->seg [i]->b_id);
  sdna_simcon->Kt += body_total_kinetic_energy (sdna_simcon->bead->b_id);

  /* potential energy */
  sdna_simcon->Ut = 0.;
  if (simcon->simparams->forces_on) sdna_simcon->Ut = force_loop (&simcon->U, sdna_simcon->dna);

  /* total energy */
  sdna_simcon->Et = sdna_simcon->Kt + sdna_simcon->Ut;

  /* JOINT SOLUTION WORKSPACE INIT */
  sdna_simcon->ws = sdna_joint_solution_workspace_alloc (sdna_simcon->dna->nsegments);
  sdna_joint_solution_workspace_init (sdna_simcon);

  /* LANGEVIN DYNAMICS PARAMETERS */

  /* mE is used to rescale the velocities in the final step of the Langevin-Euler global thermostat */
  PAR (mE) = (6*(PAR (dna_nsegments)+1.) - 3*(PAR(dna_nsegments)-1.) - 5. - 3. - 2.)/(2.*PAR(beta));

  /* init the body list, which is used by the thermostat */
  simcon->body_list = (sim_body *) malloc ((N+1) * sizeof (sim_body));
  for (i=0; i<N; i++) {
    char *name = (char *) malloc (STR_BUFF_SIZE * sizeof (char));
    sprintf (name, "DNA segment %u", i);
    sim_body_init (&simcon->body_list [i], i, name, PAR (sigma_langevin), PAR (ode_step), &sdna_simcon->dna->seg [i]->b_id);
    }
  {
    char *name = (char *) malloc (STR_BUFF_SIZE * sizeof (char));
    sprintf (name, "Bead");
    sim_body_init (&simcon->body_list [N], N, name, PAR (sigma_langevin), PAR (ode_step), &sdna_simcon->bead->b_id);
  }

  /* init thermostat */
  init_langevin_thermostat (&simcon->my_thermostat,
      PAR (thermostat_type),
      &sdna_simcon->Kt,
      PAR (sigma_langevin),
      PAR (ode_step),
      PAR (beta),
      PAR (mE));

  /* SIM INFO */
  /* thermostat parameters */
  message ("%i geoms, %ibp\n", dSpaceGetNumGeoms (simcon->space_id), simcon->simparams->dna_nbasepairs_in_segment);
  message ("Thermostat type = %s\n", simcon->simparams->thermostat_type);
  message ("<E>:%f\n", simcon->simparams->mE);
  message ("Ground height: %f nm\n", sdna_simcon->surface_height * simcon->simparams->scale_l * 1e7);

  print_dna_info (simcon, sdna_simcon->dna->seg [0], &simcon->body_list [0]);
  print_bead_info (simcon, sdna_simcon->bead, &simcon->body_list [N]);
}



/* initializes the simulated objects' initial positions and orientations */
void sdna_init_molecule_positions_orientations (simcontext * simcon, sdna_simcontext *sdna_simcon) {
  const unsigned int dna_nsegments = sdna_simcon->dna->nsegments;
  const unsigned int initial_config_has_loop = SDNA_PAR (initial_config_has_loop);
  t_vector next_position, next_u, next_v, next_t, joint_position, z_axis;
  t_real bead_force;
  if (sdna_simcon->bead->mode==FIXED_EXT) {
    t_real rho = sdna_simcon->bead->rho;
    t_real lpb = PAR (dna_lbp);
    bead_force = (rho + (1./((1.-rho)*(1.-rho)) - 1.)/4.)/lpb; /* Marko interpolation formula */
  }
  else {
    bead_force = sdna_simcon->bead->bead_force;
  }

  /* z axis */
   vector_set (z_axis, 0., 0., 1.);

  /* two cases: either the user has specified that the initial configuration has a loop,
   * or it does not have it. */
  if (initial_config_has_loop) {
    /* in this case, we calculate the number of segments that we need to assign to each section
     * of the molecule following a simple method: we divide the molecule into three parts:
     * WLC - LOOP - WLC. The first and last WLC sections have equal length, and the loop is positioned
     * at the center of the molecule */
    unsigned int i, j, loop_nsegments, wlc_1_nsegments, wlc_2_nsegments, checksum, remainder;
    dna_segment **loop, **wlc_1, **wlc_2;
    t_real loop_nturns, loop_pitch, loop_radius;

    /* initialize loop parameters */
    loop_nsegments = SDNA_PAR (loop_nsegments);
    loop_nturns = SDNA_PAR (loop_nturns);
    loop_pitch = SDNA_PAR (loop_pitch);
    loop_radius = (t_real ) loop_nsegments * PAR (dna_segment_length) / (2*M_PI*loop_nturns);

    /* test parameters */
    if (loop_nsegments>dna_nsegments) {
      err_message ("loop_nsegments cannot be greater than dna_nsegments\n");
      exit (EXIT_FAILURE);
    }

    /* calculate the length of the WLC sections */
    wlc_1_nsegments = round ((dna_nsegments-loop_nsegments)/2);
    remainder = (dna_nsegments-loop_nsegments)%2;
    wlc_2_nsegments = wlc_1_nsegments + remainder;

    /* check parameters */
    checksum = wlc_1_nsegments + loop_nsegments + wlc_2_nsegments;
    if (checksum != dna_nsegments) {
      err_message ("Checksum failed.\n");
      exit (EXIT_FAILURE);
    }

    /* WLC 1 init */
    wlc_1 = (dna_segment **) malloc (wlc_1_nsegments*sizeof (dna_segment *));
    j=0;
    for (i=0; i<wlc_1_nsegments; i++) {
      wlc_1 [i] = sdna_simcon->dna->seg [j++];
    }

    /* loop init */
    loop = (dna_segment **) malloc ((loop_nsegments+1)*sizeof (dna_segment *));
    j=wlc_1_nsegments-1;
    for (i=0; i<loop_nsegments+1; i++) {
      loop [i] = sdna_simcon->dna->seg [j++];
    }

    /* WLC 2 init */
    wlc_2 = (dna_segment **) malloc ((wlc_2_nsegments+1)*sizeof (dna_segment *));
    j=wlc_1_nsegments+loop_nsegments-1;
    for (i=0; i<wlc_2_nsegments+1; i++) {
      wlc_2 [i] = sdna_simcon->dna->seg [j++];
    }

    /* cycle until we have no collisions */
    do {
      dna_wlc_create (wlc_1_nsegments, simcon->simparams->beta, bead_force, simcon->simparams->gb, z_axis, wlc_1);
      dna_loop_create (loop_nsegments+1, loop_radius, loop_pitch, 0., loop); /* TODO: find a way of assigning initial total twist */
      dna_wlc_create (wlc_2_nsegments+1, simcon->simparams->beta, bead_force, simcon->simparams->gb, z_axis, wlc_2);

      /* After the DNA segments are assigned, assign the position of the bead */
      segment_next_joint (joint_position, sdna_simcon->dna->seg [dna_nsegments-1]);
      joint_position[2] += sdna_simcon->bead->radius;
      body_set_position (sdna_simcon->bead->b_id, joint_position);

      /* these are the new axes of the bead */
      vector_set (next_u, 1., 0., 0.);
      vector_set (next_v, 0., 1., 0.);
      vector_set (next_t, 0., 0., 1.);
      body_set_axes (sdna_simcon->bead->b_id, next_u, next_v, next_t);

      /* now check if we have collisions */
      simcon->numc = 0;
      dSpaceCollide (simcon->space_id, simcon, &collision_callback);
      if (simcon->contactjoints)
	dJointGroupEmpty (simcon->contactjoints);
    } while (simcon->numc>0);

    /* free the loop, wlc_1 and wlc_2 pointers */
    free (wlc_1);
    free (loop);
    free (wlc_2);
  }
  else {
    /* now, we try to create an initial configuration of the DNA, until we do not have contacts */
    do {
      dna_wlc_create (dna_nsegments, PAR (beta), bead_force, PAR (gb), z_axis, sdna_simcon->dna->seg);

      /* After the DNA segments are assigned, assign the position of the bead */
      body_position (next_position, sdna_simcon->dna->seg [dna_nsegments-1]->b_id);
      body_axis (next_t, sdna_simcon->dna->seg [dna_nsegments-1]->b_id, 2);
      vector_lincomb (joint_position, next_position, next_t, 1., simcon->simparams->dna_segment_length/2.);
      vector_set (next_position, joint_position[0], joint_position[1], joint_position[2] + sdna_simcon->bead->radius);
      body_set_position (sdna_simcon->bead->b_id, next_position);

      /* these are the new axes of the bead */
      vector_set (next_u, 1., 0., 0.);
      vector_set (next_v, 0., 1., 0.);
      vector_set (next_t, 0., 0., 1.);
      body_set_axes (sdna_simcon->bead->b_id, next_u, next_v, next_t);

      /* now check if we have collisions */
      simcon->numc = 0;
      dSpaceCollide (simcon->space_id, simcon, &collision_callback);
      if (simcon->contactjoints)
	dJointGroupEmpty (simcon->contactjoints);
    } while (simcon->numc>0);
  }
}



/* initializes the simulated objects' initial velocities */
void sdna_init_molecule_velocities (simcontext * simcon, sdna_simcontext *sdna_simcon) {
  /* shake the DNA */
  dna_shake (sdna_simcon->dna, simcon->body_list, 0, &simcon->mechanical_model);

  /* single step for the bead */
  langevin_euler_step (&simcon->body_list [sdna_simcon->dna->nsegments]);

  /* if user specified that we need to calculate the pair interactions, do so */
  if (simcon->simparams->forces_on) sdna_simcon->Ut = force_loop (&simcon->U, sdna_simcon->dna);

  /* perform collision check, and perform the first step */
  dSpaceCollide (simcon->space_id, simcon, &collision_callback);
  dWorldStep (simcon->world_id, simcon->simparams->ode_step);
  if (simcon->contactjoints)
    dJointGroupEmpty (simcon->contactjoints);
}



/* initializes the linking number, twist and writhe */
void sdna_init_Lk (sdna_simcontext *sdna_simcon) {
  /* DNA twist */
  sdna_simcon->Tw = dna_twist (sdna_simcon->dna)/(2*M_PI);

  /* DNA writhe */
  sdna_simcon->Wr = writhe_dna (sdna_simcon->dna, 0);

  /* DNA linking number */
  sdna_simcon->Lk = sdna_simcon->Tw + sdna_simcon->Wr;
  sdna_simcon->bead->Lk = &sdna_simcon->Lk;
}

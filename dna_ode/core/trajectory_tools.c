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

#include "trajectory_tools.h"

/* Loads a trajectory from a trajectory file. The trajectory file name is read from
 * a simulation config file. All the information will be stored in the "trajectory"
 * structure, which needs to be alloc'd (and free'd!) OUTSIDE this function */
int load_trajectory (unsigned int nbodies, char *trajectory_file_name, trajectory *traj, int (*frame_operation) (trajectory_frame *, void *), void *frame_operation_parameters) {
  unsigned int vector_size;
  FILE *trajectory_file;

  /* open trajectory file */
  trajectory_file = my_fopen (trajectory_file_name, "r");

  /* check if trajectory file exists */
  if (trajectory_file == NULL) {
    err_message ("Trajectory file %s not found\n", trajectory_file_name);
    return LOAD_TRAJECTORY_FAIL;
  }

  /* Allocate memory to hold the data. We don't know in advance how many
   * frames we will load, so we allocate a number CHUNK_SIZE of frames, and
   * we then reallocate if we need more space. */
  traj->frame = (trajectory_frame *) malloc (CHUNK_SIZE * sizeof (trajectory_frame));

  /* We are ready to scan the file now. */
  traj->nsteps = 0;
  vector_size = CHUNK_SIZE;
  do {
    int frame_load_ret_code;

    /* check if we have enough memory. If not, expand our array */
    if (traj->nsteps>vector_size-1) {
      trajectory_frame * is_null;
      vector_size += CHUNK_SIZE;
      is_null = (trajectory_frame *) realloc (traj->frame, vector_size * sizeof (trajectory_frame));
      if (is_null == NULL) {
	err_message ("No more memory!\n");
	return LOAD_TRAJECTORY_FAIL;
      }
      traj->frame = is_null;
    }

    /* load a trajectory frame from file and increment the step number */
    traj->frame [traj->nsteps].frame_id = traj->nsteps;
    frame_load_ret_code = load_trajectory_frame (nbodies, trajectory_file, &traj->frame [traj->nsteps], frame_operation, frame_operation_parameters);
    traj->nsteps++;

    /* check that we got success return code */
    if (frame_load_ret_code==LOAD_TRAJ_FRAME_FAIL) {
      err_message ("Error reading trajectory frame. Read %d frames.\n", traj->nsteps);
      return LOAD_TRAJECTORY_FAIL;
    }

    if (frame_load_ret_code==LOAD_TRAJ_FRAME_SUCCESS_STOP) {
      fclose (trajectory_file);
      return LOAD_TRAJECTORY_SUCCESS;
    }
  }
  while (!feof (trajectory_file));

  /* close trajectory file and return */
  fclose (trajectory_file);
  return LOAD_TRAJECTORY_SUCCESS;
}



/* Loads a single trajectory frame from a trajectory file. This function assumes that
 * the trajectory file is already open, and has its register in the right position so
 * that it scans starting from the right position. */
int load_trajectory_frame (unsigned int nbodies, FILE *trajectory_file, trajectory_frame *frame, int (*frame_operation) (trajectory_frame *, void *), void *frame_operation_parameters) {
  int retcode;
  unsigned int i, id;
  t_real sim_clock, x, y, z, vx, vy, vz, wx, wy, wz;
  t_real q[4];

  /* alloc memory for the frame */
  frame->nbodies = nbodies;
  frame->body_i_state = (body_state *) malloc (nbodies * sizeof (body_state));
  frame->data = NULL;

  /* check that the system gave the memory */
  if (frame->body_i_state == NULL) {
    err_message ("Out of memory!\n");
    return LOAD_TRAJ_FRAME_FAIL;
  }

  /* we can now read nbodies lines from the file, and assign the trajectory frame data */
  for (i=0; i<nbodies; i++) {
    /* read a single line of the trajectory file and check that it is read correctly */
    retcode = fscanf (trajectory_file, "%lf %u %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
	&sim_clock, &id, &x, &y, &z, &q[0], &q[1], &q[2], &q[3], &vx, &vy, &vz, &wx, &wy, &wz);

    /* check that the line contains the right information */
    if (retcode != 15) {
      err_message ("Error reading trajectory file: line contains incorrect number of data. Expected 15, got %d\n", retcode);
      err_message ("sim_clock = %f body_id = %u\n", sim_clock, id);
      return LOAD_TRAJ_FRAME_FAIL;
    }

    /* the body id that we read from the trajectory file should correspond to i. Otherwise
     * there was a problem either in the positioning of the cursor at the file, or there
     * was a problem writing the trajectory file (ambiguous body id). */
    if (id!=i) {
      err_message ("Error reading trajectory file: invalid body id\n");
      return LOAD_TRAJ_FRAME_FAIL;
    }

    /* now we assign all the body info */
    if (id<nbodies) {
      frame->sim_clock = sim_clock;
      frame->body_i_state [id].id = id;
      vector_set (frame->body_i_state [id].r, x, y, z);
      axis_from_quaternion (frame->body_i_state [id].u, q, 0);
      axis_from_quaternion (frame->body_i_state [id].v, q, 1);
      axis_from_quaternion (frame->body_i_state [id].t, q, 2);
      vector_set (frame->body_i_state [id].vel, vx, vy, vz);
      vector_set (frame->body_i_state [id].omega, wx, wy, wz);
    }
    else {
      err_message ("Error reading trajectory file: body id exceeds nbodies\n");
      return LOAD_TRAJ_FRAME_FAIL;
    }
  }

  /* now invoke the user-defined frame_operation function */
  if (frame_operation != NULL)
    return frame_operation (frame, frame_operation_parameters);

  return LOAD_TRAJ_FRAME_SUCCESS;
}



/* Places a cursor on the n-th frame of a trajectory file. Needs an open file. */
int skip_to_frame_n (unsigned int n, unsigned int nbodies, FILE *trajectory_file) {
  unsigned int i;

  /* check if trajectory file exists */
  if (trajectory_file == NULL) {
    err_message ("Trajectory file not open or non-existent\n");
    return LOAD_TRAJECTORY_FAIL;
  }

  /* places cursor at beginning of stream */
  rewind (trajectory_file);

  /* cycle to the line number n*nbodies */
  i=0;
  while (i != (n*nbodies)) {
    char is_endline;
    is_endline = fgetc (trajectory_file);
    if (is_endline == '\n') i++;
    /* printf ("char = \"%d\", i = %d\n", is_endline); */
    if (feof (trajectory_file)) {
      err_message ("Trajectory read failed. Reached end of file\n");
      return LOAD_TRAJECTORY_FAIL;
    }
  }

  /* close trajectory file and return */
  return LOAD_TRAJECTORY_SUCCESS;
}



void free_trajectory_frame (trajectory_frame *frame) {
  free (frame->body_i_state);
  if (frame->data!=NULL)
    free (frame->data);
}


/* free all the memory associated with a trajectory */
void free_trajectory (trajectory *traj) {
  unsigned int i;
  for (i=0; i<traj->nsteps; i++) {
    free_trajectory_frame (&traj->frame [i]);
    /* free (traj->frame [i].body_i_state);
    if (traj->frame [i].data!=NULL) free (traj->frame [i].data); */
  }
  free (traj->frame);
}



/* OUTPUT UTILITIES */
void print_body_state (FILE *out, trajectory_frame *frame, unsigned int body_id) {
  t_real q [4];
  body_state *b = &frame->body_i_state [body_id];

  /* convert axes to quaternion */
  quaternion_from_axes (q, b->u, b->v, b->t);

  file_message (out, "%.16f %u %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n",
      frame->sim_clock, body_id, b->r[0], b->r[1], b->r[2], q [0], q[1], q[2], q[3],
      b->vel[0], b->vel[1], b->vel[2], b->omega[0], b->omega[1], b->omega[2]);
}


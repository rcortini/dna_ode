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

#ifndef __TRAJECTORY_TOOLS_H__
#define __TRAJECTORY_TOOLS_H__

#include "dna_ode_core.h"
#include "config_read.h"
#include "dna_ode_math.h"

/* useful defines */

#define LOAD_TRAJECTORY_SUCCESS          1001

#define LOAD_TRAJECTORY_FAIL             1002

#define LOAD_TRAJ_FRAME_SUCCESS          1003

#define LOAD_TRAJ_FRAME_FAIL             1004

#define LOAD_TRAJ_FRAME_SUCCESS_STOP     1005



/* this structure contains information on the state of a body */
struct body_state {
  unsigned int id;
  t_vector r;
  t_vector u;
  t_vector v;
  t_vector t;
  t_vector vel;
  t_vector omega;
};
typedef struct body_state body_state;

/* this structure contains all the information on a single frame in a trajectory */
struct trajectory_frame {
  unsigned int frame_id;
  unsigned int nbodies;
  t_real sim_clock;
  body_state *body_i_state;
  void *data;    /* this might become handy as we may calculate additional data that we may want to store */
};
typedef struct trajectory_frame trajectory_frame;

/* this structure contains all the information on a trajectory */
struct trajectory {
  unsigned int nsteps;
  trajectory_frame *frame;
};
typedef struct trajectory trajectory;

int load_trajectory_frame (unsigned int nbodies, FILE *trajectory_file, trajectory_frame *frame, int (*frame_operation) (trajectory_frame *, void *), void *frame_operation_parameters);

int load_trajectory (unsigned int nbodies, char *config_file_name, trajectory *traj, int (*frame_operation) (trajectory_frame *, void *), void *frame_operation_parameters);

int skip_to_frame_n (unsigned int n, unsigned int nbodies, FILE *trajectory_file);

void free_trajectory_frame (trajectory_frame *frame);

void free_trajectory (trajectory *traj);

void print_body_state (FILE *out, trajectory_frame *frame, unsigned int body_id);

#endif

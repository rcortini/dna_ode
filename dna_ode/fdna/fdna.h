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

#ifndef __FDNA_H__
#define __FDNA_H__

#include "dna_ode_core.h"
#include "dna_ode.h"

#define FDNA_PAR(x) (fdna_simcon->x)

/*
 * ---------------------------------------
 *  DATA STRUCTURE DECLARATION
 * ---------------------------------------
 */



/* data structure to contain memory necessary for joint solution */
struct fdna_joint_solution_workspace {
  t_vector *t;
  t_vector *rhs;
  t_vector *lambda;
  t_vector *buffer;
  t_matrix *Hii;
  t_matrix *Hij;
  t_matrix *Dii;
  t_matrix *Dii_inv;
  t_matrix *Lij;
};
typedef struct fdna_joint_solution_workspace fdna_joint_solution_workspace;

struct fdna_simcontext {
  t_real dTw;

  /* energy of the system */
  t_real Et;  /* total energy */
  t_real Ut;  /* total potential energy */
  t_real ext_Ut;  /* total external potential energy */
  t_real Kt;  /* total kinetic energy */

  /* simulated objects */
  dna_ode *dna;

  /* workspace for joints solution */
  fdna_joint_solution_workspace *ws;
};
typedef struct fdna_simcontext fdna_simcontext;



/*
 * OBJECT CREATION (fdna_objects_creation.c)
 */



dna_ode * fdna_create (simcontext *simcon);

void fdna_joints_create (simcontext * simcon, dna_ode * dna);



/*
 * DYNAMICS FUNCTIONS (fdna_dynamics.c)
 */



void fdna_simstep (simcontext * simcon, fdna_simcontext * fdna_simcon);



/*
 * INIT FUNCTIONS PROTOTYPE DECLARATIONS (fdna_init_functions.c)
 */



void fdna_load_initial_configuration_from_file (simcontext *simcon, fdna_simcontext *fdna_simcon);

void fdna_init_simulation_context (simcontext *simcon, fdna_simcontext *fdna_simcon);

void fdna_init_molecule_positions_orientations (simcontext * simcon, fdna_simcontext * fdna_simcon);

void fdna_init_molecule_velocities (simcontext * simcon, fdna_simcontext * fdna_simcon);



/*
 * JOINTS SOLUTION FUNCTIONS (fdna_joints_solution.c)
 */



fdna_joint_solution_workspace * fdna_joint_solution_workspace_alloc (const unsigned int N);

void fdna_joint_solution_workspace_init (fdna_simcontext *fdna_simcon);

void fdna_joint_solution_workspace_free (fdna_joint_solution_workspace *ws);

void fdna_joints_solution (simcontext *simcon, fdna_simcontext *fdna_simcon);



/*
 * FINALIZATION FUNCTIONS PROTOTYPE DECLARATIONS (fdna_fin_functions.c)
 */



void fdna_free_memory (fdna_simcontext *fdna_simcon);



/*
 * FUNCTIONS FOR PRINTING OUTPUT TO FILES (fdna_output_functions.c)
 */



void fdna_string_simulation_data (char * line, fdna_simcontext * fdna_simcon);

void fdna_write_trajectory_file (simcontext * simcon, fdna_simcontext * fdna_simcon);

void fdna_write_data_file (simcontext * simcon, fdna_simcontext * fdna_simcon);

void fdna_output_data (simcontext * simcon, fdna_simcontext * fdna_simcon);



/*
 * FUNCTIONS FOR PRINTING OUTPUT TO TERMINAL (fdna_print_functions.c)
 */



void fdna_print_data (simcontext * simcon, fdna_simcontext * fdna_simcon);



#endif

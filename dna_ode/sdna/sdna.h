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

#ifndef __SDNA_H__
#define __SDNA_H__

#include "dna_ode_core.h"
#include "dna_ode.h"

#define SDNA_PAR(x) (sdna_simcon->x)

/*
 * ---------------------------------------
 *  DATA STRUCTURE DECLARATION
 * ---------------------------------------
 */

/* data structure to contain memory necessary for joint solution */
struct sdna_joint_solution_workspace {
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
typedef struct sdna_joint_solution_workspace sdna_joint_solution_workspace;

struct sdna_simcontext {
  /* initial loop parameters */
  unsigned int initial_config_has_loop;
  unsigned int loop_nsegments;
  t_real loop_nturns;
  t_real loop_pitch;

  /* surface parameters */
  t_real surface_height_to_reff_ratio;
  t_real surface_height;

  /* DNA variables */
  t_real Lk;
  t_real Tw;
  t_real Wr;

  /* energy of the system */
  t_real Et;  /* total energy */
  t_real Ut;  /* total potential energy */
  t_real Kt;  /* total kinetic energy */

  /* simulated objects */
  surface_ode *surface;
  bead_ode *bead;
  dna_ode *dna;
  dBodyID *body_list;

  /* memory for joint solution */
  sdna_joint_solution_workspace *ws;
};
typedef struct sdna_simcontext sdna_simcontext;



/*
 * DYNAMICS FUNCTIONS (sdna_dynamics.c)
 */



void sdna_simstep (simcontext * simcon, sdna_simcontext *sdna_simcon);



/*
 * CREATION OF SIMULATED OBJECTS FUNCTIONS (sdna_objects_creation.c)
 */



dna_ode * sdna_create (simcontext *simcon);

void sdna_joints_create (simcontext * simcon, sdna_simcontext *sdna_simcon, dna_ode * dna);



/*
 * INIT FUNCTIONS PROTOTYPE DECLARATIONS (sdna_init_functions.c)
 */



void sdna_load_initial_configuration_from_file (simcontext *simcon, sdna_simcontext *sdna_simcon);

void sdna_init_simulation_context (simcontext *simcon, sdna_simcontext *sdna_simcon);

void sdna_init_molecule_positions_orientations (simcontext * simcon, sdna_simcontext *sdna_simcon);

void sdna_init_molecule_velocities (simcontext * simcon, sdna_simcontext *sdna_simcon);

void sdna_init_Lk (sdna_simcontext *sdna_simcon);



/*
 * FINALIZATION FUNCTIONS PROTOTYPE DECLARATIONS (sdna_fin_functions.c)
 */



void sdna_free_memory (sdna_simcontext *sdna_simcon);



/*
 * FUNCTIONS FOR PRINTING OUTPUT TO FILES (sdna_output_functions.c)
 */



void sdna_string_simulation_data (char * line, simcontext *simcon, sdna_simcontext *sdna_simcon);

void sdna_write_data_file (simcontext * simcon, sdna_simcontext *sdna_simcon);

void sdna_output_data (simcontext * simcon, sdna_simcontext *sdna_simcon);



/*
 * FUNCTIONS FOR PRINTING OUTPUT TO TERMINAL (sdna_print_functions.c)
 */



void sdna_print_data (simcontext * simcon, sdna_simcontext *sdna_simcon);

void print_joints_data (sdna_simcontext *sdna_simcon);

/*
 * FUNCTIONS FOR THE EXACT SOLUTION OF THE JOINT FORCES (sdna_joints_solution.c)
 */


sdna_joint_solution_workspace * sdna_joint_solution_workspace_alloc (const unsigned int N);

void sdna_joint_solution_workspace_free (sdna_joint_solution_workspace *ws);

void sdna_joint_solution_workspace_init (sdna_simcontext *sdna_simcon);

void sdna_joints_solution (simcontext *simcon, sdna_simcontext *sdna_simcon);


#endif

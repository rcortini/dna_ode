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

#ifndef __CIRCULAR_DNA_H__
#define __CIRCULAR_DNA_H__



#include "dna_ode_core.h"
#include "dna_ode.h"



/*
 * ---------------------------------------
 *  DATA STRUCTURE DECLARATION
 * ---------------------------------------
 */



/* a very convenient macro */
#define CDNA_PAR(x)     cdna_simcon->x



/* data structure to contain memory necessary for joint solution */
struct ccdna_joint_solution_workspace {
  t_vector *t;
  t_vector *rhs;
  t_vector *lambda;
  t_vector *buffer;
  t_matrix *Hii;
  t_matrix *Hij;
  t_matrix *Dii;
  t_matrix *Dii_inv;
  t_matrix *Lij;
  t_matrix *LNj;
};
typedef struct ccdna_joint_solution_workspace ccdna_joint_solution_workspace;



/* system-specific variables and structures */
struct cdna_simcontext {
  /* circular DNA parameters */
  t_real sigma;
  t_real sigma_prime;
  int nturns;
  t_real delta_Lk0;

  /* DNA variables */
  t_real Wr;
  t_real dTw;
  t_real Lk;

  /* energy of the system */
  t_real Et;  /* total energy */
  t_real Ut;  /* total potential energy */
  t_real ext_Ut;  /* total external potential energy */
  t_real Kt;  /* total kinetic energy */

  /* simulated objects */
  dna_ode *dna;

  /* workspace for joint solution */
  ccdna_joint_solution_workspace *ws;
};
typedef struct cdna_simcontext cdna_simcontext;



/*
 * DYNAMICS FUNCTIONS (cdna_dynamics.c)
 */



void cdna_simstep (simcontext * simcon, cdna_simcontext *cdna_simcon);



/*
 * CREATION OF SIMULATED OBJECTS FUNCTIONS (cdna_objects_creation.c)
 */



dna_ode * cdna_create (simcontext *simcon);

void cdna_joints_create (simcontext * simcon, dna_ode * dna);



/*
 * INIT FUNCTIONS PROTOTYPE DECLARATIONS (cdna_init_functions.c)
 */



void cdna_load_initial_configuration_from_file (simcontext *simcon, cdna_simcontext *cdna_simcon);

void cdna_init_simulation_context (simcontext *simcon, cdna_simcontext *cdna_simcon);

void cdna_init_molecule_positions_orientations (simcontext * simcon, cdna_simcontext *cdna_simcon);

void cdna_init_molecule_velocities (simcontext * simcon, cdna_simcontext *cdna_simcon);



/*
 * EXACT SOLUTION FOR JOINTS (cdna_joints_solution.c)
 */



ccdna_joint_solution_workspace * ccdna_joint_solution_workspace_alloc (const unsigned int N);

void ccdna_joint_solution_workspace_init (cdna_simcontext *cdna_simcon);

void ccdna_joint_solution_workspace_free (ccdna_joint_solution_workspace *);

void ccdna_joints_solution (simcontext *simcon, cdna_simcontext *cdna_simcon);



/*
 * FINALIZATION FUNCTIONS PROTOTYPE DECLARATIONS (cdna_fin_functions.c)
 */



void cdna_free_memory (cdna_simcontext *cdna_simcon);



/*
 * FUNCTIONS FOR PRINTING OUTPUT TO FILES (cdna_output_functions.c)
 */


void cdna_string_simulation_data (char * line, cdna_simcontext *cdna_simcon);

void cdna_write_data_file (simcontext * simcon, cdna_simcontext *cdna_simcon);

void cdna_output_data (simcontext * simcon, cdna_simcontext *cdna_simcon);



/*
 * FUNCTIONS FOR PRINTING OUTPUT TO TERMINAL (cdna_print_functions.c)
 */



void cdna_print_data (simcontext * simcon, cdna_simcontext *cdna_simcon);



#endif

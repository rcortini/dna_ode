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

#ifndef __DNA_ODE_CORE_H__
#define __DNA_ODE_CORE_H__



/*
 * ---------------------------------------
 *  HEADERS
 * ---------------------------------------
 */



/* standard library headers */
#include <assert.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <time.h>

/* GSL headers */
#include <gsl/gsl_sf.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_const_num.h>
#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_multifit_nlin.h>

/* ODE headers */
#include <ode/ode.h>
#include <ode/odemath.h>



/*
 * ---------------------------------------
 *  DEBUG SYMBOLS
 * ---------------------------------------
 */



#ifdef DEBUG
  #define DPRINT(...) {printf ( "%s) ", __FUNCTION__ ); printf (""__VA_ARGS__); fflush (stderr);}
#else
  #define DPRINT(...)
#endif



/*
 * ---------------------------------------
 *  PHYSICAL CONSTANTS
 * ---------------------------------------
 */



#define K_B 1.380650400e-16
#define N_AVO GSL_CONST_NUM_AVOGADRO
#define E_CHARGE  4.80320427e-10



/*
 * ---------------------------------------
 *  USEFUL CONSTANTS
 * ---------------------------------------
 */



/* maximum number of contacts in collision_callback */
#define MAX_CONTACTS 64

/* length of a string buffer */
#define STR_BUFF_SIZE 1024        

/* contact mode definition */
#define CONTACT_MODE       (dContactSoftERP | dContactSoftCFM | dContactBounce | dContactMu2);

/* identity of a surface "dBodyID" */
#define SURFACE_BODY_ID   0

/* size of the initial buffer of an unknown size vector */
#define CHUNK_SIZE (int) 100



/*
 * ---------------------------------------
 *  DATA STRUCTURE DECLARATION
 * ---------------------------------------
 */



#define t_real dReal

#define t_vector dVector3

#define t_matrix dMatrix3

#define PAR(x) (simcon->simparams->x)


/* a structure used to contain all body information */
struct sim_body {
  char *name;
  unsigned int odepdb_id;
  dBodyID *b_id;
  t_real trans_friction;
  t_real trans_sigma;
  t_real rot_friction_x;
  t_real rot_friction_y;
  t_real rot_friction_z;
  t_real rot_sigma_x;
  t_real rot_sigma_y;
  t_real rot_sigma_z;
};
typedef struct sim_body sim_body;


/* contains all the information on the DNA segment */
struct dna_segment {
  unsigned int odepdb_id;
  unsigned int id;
  dGeomID g_id;
  dBodyID b_id;
  t_real mass;
  t_real length;
  t_real radius;
  t_matrix inertia;
};
typedef struct dna_segment dna_segment;

/* contains all the information on a joint */
struct joint_ode {
  unsigned int id;
  t_real twist;
  dJointID j_id;
  dJointFeedback jf_id;
  dBodyID *previous_body;
  dBodyID *next_body;
  t_real next_body_size;
};
typedef struct joint_ode joint_ode;

/* this data structure contains all the information on the dna */
struct dna_ode {
  unsigned int nsegments;
  unsigned int njoints;
  unsigned int locked;
  unsigned int nicked;
  dna_segment **seg;
  joint_ode **joint;
  t_real Lk0;
  t_real dTw;
  unsigned int dna_strand_id;
};
typedef struct dna_ode dna_ode;

/* a structure to contain the information on the surface */
struct surface_ode {
  dBodyID b_id; /* a plane does not have a body ID in ODE, but I will set this to zero for consistency with the rest of the code */
  dGeomID g_id;
  t_real a;
  t_real b;
  t_real c;
  t_real d;
};
typedef struct surface_ode surface_ode;



/*
 * ---------------------------------------
 * SYSTEM IDs
 * ---------------------------------------
 */



typedef enum {
  SDNA}
system_id;



#endif

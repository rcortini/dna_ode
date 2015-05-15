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

#include "config.h"
#include "dna_ode_core.h"
#include "dna_ode.h"



/* this function writes the formatted data about a geometry */
void string_geometry_data (char * line, dGeomID g) {
  const t_real *m = dGeomGetPosition (g);
  const t_real *d = dGeomGetRotation (g);
  sprintf (line, "%.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e",
      m [0],
      m [1],
      m [2],
      d [0],
      d [4],
      d [8],
      d [1],
      d [5],
      d [9],
      d [2],
      d [6],
      d [10]);
}



/* this function writes the formatted data about a body */
void string_body_data (char * line, dBodyID b) {
  const t_real * pos = dBodyGetPosition (b);
  const t_real * quat = dBodyGetQuaternion (b);
  const t_real * vel = dBodyGetLinearVel (b);
  const t_real * ang_vel = dBodyGetAngularVel (b);
  sprintf (line, "%.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e",
      pos [0],
      pos [1],
      pos [2],
      quat [0],
      quat [1],
      quat [2],
      quat [3],
      vel [0],
      vel [1],
      vel [2],
      ang_vel [0],
      ang_vel [1],
      ang_vel [2]);
}




/* writes the data to the "geometry" file */
void write_frame (const unsigned int nbodies, t_real sim_clock, FILE *trajectory_file, FILE *restart_file, sim_body *body_list) {
  unsigned int i;
  char line [STR_BUFF_SIZE];

  /* cycle on all simulated bodies */
  for (i=0; i<nbodies; i++) {
    string_body_data (line, *body_list [i].b_id);
    file_message (trajectory_file, "%.5f %u %s\n", sim_clock, body_list [i].odepdb_id, line);
    file_message (restart_file, "%.5f %u %s\n", sim_clock, body_list [i].odepdb_id, line);
  }
}



/* write final information to terminal */
void finish_output (simcontext * simcon) {
  time_t now;
  struct tm tstruct;
  char time_string [STR_BUFF_SIZE];

  /* get the current time and date and convert it to a human-readable string */
  now = time (0);
  tstruct = *localtime (&now);
  strftime (time_string, STR_BUFF_SIZE, "%Y-%m-%d %X", &tstruct);

  /* write total number of steps at the end */
  message ("Steps: %d      Time: %.3f (%.2e s)\n",
      simcon->sim_tick, simcon->sim_clock, simcon->sim_clock*simcon->simparams->scale_t);

  /* finish time */
  message ("Job finished at %s\n", time_string);
} 

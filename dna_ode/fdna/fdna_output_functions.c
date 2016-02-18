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
#include "fdna.h"
#include "config.h"



/*
 * UTILITY FUNCTIONS FOR PRINTING OUTPUT TO FILES
 */



/* TRAJECTORY FILE */



/* SIMULATION DATA FILE */



/* this function writes the data that we need to write about the simulation */
void string_simulation_data (char * line, fdna_simcontext * fdna_simcon) {
  sprintf (line, "%.16e %.16e %.16e %.16e",
      fdna_simcon->Et,
      fdna_simcon->Ut,
      fdna_simcon->ext_Ut,
      fdna_simcon->Kt);
}



/* writes the simulation data */
void fdna_write_data_file (simcontext * simcon, fdna_simcontext * fdna_simcon) {
  char line [STR_BUFF_SIZE];
  string_simulation_data (line, fdna_simcon);
  file_message (simcon->data_file, "%.16f %s\n", simcon->sim_clock, line);
}



/* DATA OUTPUT */



/* once every output_interval number of steps, print output to the output files */
void fdna_output_data (simcontext * simcon, fdna_simcontext * fdna_simcon) {
  /* writes all the data invoking the corresponding utility functions */
  simcon->restart_file = my_fopen (simcon->simparams->restart_file_name, "w");
  write_frame (simcon->nbodies, simcon->sim_clock, simcon->trajectory_file, simcon->restart_file, simcon->body_list);
  fdna_write_data_file (simcon, fdna_simcon);
  fclose (simcon->restart_file);
}

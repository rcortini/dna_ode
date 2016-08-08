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
#include "writhe.h"
#include "cdna.h"
#include "config.h"



/*
 * UTILITY FUNCTIONS FOR PRINTING OUTPUT TO FILES
 */



/* SIMULATION DATA FILE */



/* this function writes the data that we need to write about the simulation */
void cdna_string_simulation_data (char * line, cdna_simcontext *cdna_simcon) {

  /* calculate the writhe of the DNA molecule after a time step */
  cdna_simcon->Wr = writhe_dna (cdna_simcon->dna, 1);
  cdna_simcon->Lk = cdna_simcon->Wr + cdna_simcon->dTw;

  sprintf (line, "%.16e %.16e %.16e %.16e %.16e %.16e %.16e",
      cdna_simcon->Wr,
      cdna_simcon->dTw,
      cdna_simcon->Lk,
      cdna_simcon->Et,
      cdna_simcon->Ut,
      cdna_simcon->ext_Ut,
      cdna_simcon->Kt);
}



/* writes the simulation data */
void cdna_write_data_file (simcontext * simcon, cdna_simcontext *cdna_simcon) {
  char line [STR_BUFF_SIZE];
  cdna_string_simulation_data (line, cdna_simcon);
  file_message (simcon->data_file, "%.16f %s\n", simcon->sim_clock, line);
}



/* DATA OUTPUT */



/* once every output_interval number of steps, print output to the output files */
void cdna_output_data (simcontext * simcon, cdna_simcontext *cdna_simcon) {
  /* writes all the data invoking the corresponding utility functions */
  simcon->restart_file = my_fopen (simcon->simparams->restart_file_name, "w");
  write_frame (simcon->nbodies, simcon->sim_clock, simcon->trajectory_file, simcon->restart_file, simcon->body_list);
  cdna_write_data_file (simcon, cdna_simcon);
  fclose (simcon->restart_file);
}

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
#include "sdna.h"
#include "config.h"



/*
 * UTILITY FUNCTIONS FOR PRINTING OUTPUT TO FILES
 */








/* TRAJECTORY FILE */



/* writes the data to the "geometry" file */
void sdna_write_trajectory_file (simcontext * simcon, sdna_simcontext *sdna_simcon) {
  unsigned int i;
  char line [STR_BUFF_SIZE];

  /* cycle on all the DNA segments */
  for (i=0; i<sdna_simcon->dna->nsegments; i++) {
    string_body_data (line, sdna_simcon->dna->seg [i]->b_id);
    file_message (simcon->trajectory_file, "%.16f %u %s\n", simcon->sim_clock, sdna_simcon->dna->seg [i]->id, line);
  }
  string_body_data (line, sdna_simcon->bead->b_id);
  file_message (simcon->trajectory_file, "%.16f %u %s\n", simcon->sim_clock, i, line);
}



/* SIMULATION DATA FILE */



/* this function writes the data that we need to write about the simulation */
void sdna_string_simulation_data (char * line, simcontext * simcon, sdna_simcontext *sdna_simcon) {
  t_vector last_seg_pos;

  /* get last joint position */
  body_point_position (last_seg_pos, sdna_simcon->dna->seg [sdna_simcon->dna->nsegments-1]->b_id, 0., 0., 0.5 * simcon->simparams->dna_segment_length);

  /* calculate writhe */
  sdna_simcon->Wr = writhe_dna (sdna_simcon->dna, 0);

  sprintf (line, "%.16e %.16e %.16e %.16e %.16e %.16e %.16e",
      last_seg_pos [2] / simcon->simparams->dna_length,
      sdna_simcon->Lk,
      sdna_simcon->Tw,
      sdna_simcon->Wr,
      sdna_simcon->Et,
      sdna_simcon->Ut,
      sdna_simcon->Kt);
}



/* writes the simulation data */
void sdna_write_data_file (simcontext * simcon, sdna_simcontext *sdna_simcon) {
  char line [STR_BUFF_SIZE];
  sdna_string_simulation_data (line, simcon, sdna_simcon);
  file_message (simcon->data_file, "%.16f %s\n", simcon->sim_clock, line);
}



/* DATA OUTPUT */



/* once every output_interval number of steps, print output to the output files */
void sdna_output_data (simcontext * simcon, sdna_simcontext *sdna_simcon) {
  /* writes all the data invoking the corresponding utility functions */
  simcon->restart_file = my_fopen (simcon->simparams->restart_file_name, "w");
  write_frame (simcon->nbodies, simcon->sim_clock, simcon->trajectory_file, simcon->restart_file, simcon->body_list);
  sdna_write_data_file (simcon, sdna_simcon);
  fclose (simcon->restart_file);
}

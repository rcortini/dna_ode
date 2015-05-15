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

#include <getopt.h>
#include "dna_ode_core.h"
#include "dna_ode.h"
#include "dna_ode-analyze.h"



/* usage message */
void print_usage_analyze (const char *program_name) {
  err_message ("Usage: %s <config_file> <analysis_type> <analysis_specific_options>\n", program_name);
  err_message ("\nSupported analysis types:\n");
  err_message ("\tinfo\n");
  err_message ("\tangle_hist <bin_number> <min_angle> <max_angle> <output_file>\n");
  err_message ("\tdistance_hist <bin_number> <max_distance>\n");
  err_message ("\tdistance_map\n");
  err_message ("\tframe_dump <frame_number>\n");
  err_message ("\tbody_trace <segment_to_trace>\n");
  err_message ("\tbead_force\n");
  err_message ("\tlast_segment_position\n");
  err_message ("\tlp_stat <output_file_basename>\n");
}



int dna_ode_analyze (int argc, char *argv [], config_t *simconfig) {
  char *analysis_type;
  int (*analysis) (int, char **, config_t *);

  /* check if we have at least one parameter from command line */
  if (argc<4) {
    print_usage_analyze (DNA_ODE_NAME);
    return 1;
  }

  /* get the string that contains the information on what type of analysis we need to perform */
  analysis_type = argv [3];

  /* now let's understand what is the analysis that we need to do */
  if (strcmp (analysis_type, "info")==0) {
    analysis = ta_info;
  }
  else if (strcmp (analysis_type, "angle_hist")==0) {
    analysis = ta_angle_hist;
  }
  else if (strcmp (analysis_type, "distance_hist")==0) {
    analysis = ta_distance_hist;
  }
  else if (strcmp (analysis_type, "distance_map")==0) {
    analysis = ta_distance_map;
  }
  else if (strcmp (analysis_type, "frame_dump")==0) {
    analysis = ta_frame_dump;
  }
  else if (strcmp (analysis_type, "body_trace")==0) {
    analysis = ta_body_trace;
  }
  else if (strcmp (analysis_type, "bead_force")==0) {
    analysis = ta_bead_force;
  }
  else if (strcmp (analysis_type, "last_segment_position")==0) {
    analysis = ta_last_segment_position;
  }
  else {
    err_message ("Unknown analysis type: %s\n", analysis_type);
    print_usage_analyze (DNA_ODE_NAME);
    return 1;
  }

  /* do the analysis */
  return analysis (argc, argv, simconfig);
}

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
#include "cdna.h"



/*
 * UTILITY FUNCTIONS FOR PRINTING OUTPUT 
 */



/* once every print_interval number of steps, this function is called and output is shown on screen */
void cdna_print_data (simcontext * simcon, cdna_simcontext *cdna_simcon) {
  message ("step %lu     Wr: %.3f  dTw: %.3f  Lk: %.3f        Et: %.2f  Ut: %.2f  ext_Ut: %.2f  Kt: %.2f\n",
      simcon->sim_tick,
      cdna_simcon->Wr,
      cdna_simcon->dTw,
      cdna_simcon->Lk,
      cdna_simcon->Et,
      cdna_simcon->Ut,
      cdna_simcon->ext_Ut,
      cdna_simcon->Kt);
}

/* dna_ode: simulating DNA using the ODE physics engine
 *
 * Copyright (C) 2014, 2015, 2016  Ruggero Cortini

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

#ifndef __CONFIG_READ_H__
#define __CONFIG_READ_H__


#include <libconfig.h>



void * mycfg_safe_read (config_t * conf, char * paramname, unsigned int paramtype);

int mycfg_read_int (config_t * conf, char * paramname);

long mycfg_read_int64 (config_t * conf, char * paramname);

double mycfg_read_double (config_t * conf, char * paramname);

char * mycfg_read_string (config_t * conf, char * paramname);

unsigned int mycfg_read_bool (config_t * conf, char * paramname);



#endif

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

#ifndef __LANGEVIN_H__
#define __LANGEVIN_H__

#include "dna_ode_core.h"
#include "dna_ode_math.h"
#include "ode_body_functions.h"
#include "random.h"
#include "utils.h"

typedef enum {GLOBAL, LOCAL, OFF} thermostat_type;

struct thermostat {
  thermostat_type type;
  void *params;
};
typedef struct thermostat thermostat;

struct global_thermostat_params {
  t_real beta;
  t_real ode_step;
  t_real sigma_langevin;
  t_real mE;
  t_real *Kt;
};

void sim_body_init (sim_body *body,
    unsigned int odepdb_id,
    char *name,
    t_real sigma_langevin,
    t_real ode_step,
    dBodyID *b_id);

void init_langevin_thermostat (thermostat *my_thermostat,
    char *thermostat_type,
    t_real *Kt,
    t_real sigma_langevin,
    t_real ode_step,
    t_real beta,
    t_real mE);

void thermostat_all_bodies (thermostat *my_thermostat, sim_body *body_list, const unsigned int nbodies);

void langevin_euler_step (sim_body *body);

void langevin_euler_global_step (dBodyID b, dReal scale_factor);


#endif

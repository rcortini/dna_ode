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

#ifndef __DNA_H__
#define __DNA_H__

#include "dna_ode_core.h"
#include "dna_ode_math.h"
#include "random.h"
#include "bead.h"
#include "ode_body_functions.h"
#include "mechanical_models.h"
#include "curves.h"
#include "langevin.h"


double gb_function (double x, void *p);

void dna_gb_gt (t_real kuhn_length, t_real segment_length, t_real ltp, t_real tol, t_real *gb, t_real *gt);

t_real dna_twist (dna_ode *dna);

t_real dna_twist_and_torque (dna_ode *dna, mechanical_model_f *);

void dna_shake (dna_ode *dna, sim_body *body_list, unsigned int nstart, mechanical_model_f *m);

void place_dna_segments_on_curve (unsigned int nsegments, dna_segment **seg, curve_f *f, t_real total_twist);

void dna_loop_create (unsigned int nsegments, t_real loop_radius, t_real loop_pitch, t_real total_twist, dna_segment **seg); 

void dna_wlc_create (unsigned int nsegments, t_real beta, t_real force, t_real gb, const t_real * wlc_axis, dna_segment **seg);

void segment_next_joint (t_real *res, dna_segment *seg);

void segment_next_t (t_real *res, dna_segment *seg);

void bead_joint (t_real *res, bead_ode *bead);

void segment_joint_pos_error (t_real *res, dna_segment *previous, dna_segment *next);

void segment_joint_vel_error (t_real *res, dna_segment *previous, dna_segment *next);

void bead_joint_pos_error (t_real *res, dna_segment *last, bead_ode *bead);

void bead_joint_vel_error (t_real *res, dna_segment *last, bead_ode *bead);

void segment_rhs (t_real *rhs,
    dna_segment *prev_seg,
    dna_segment *next_seg,
    const t_real *prev_t,
    const t_real *next_t,
    const t_real ks,
    const t_real kd);

void segment_rhs_ground (t_real *rhs,
    dna_segment *seg,
    const t_real *t,
    const t_real ks,
    const t_real kd);

void segment_rhs_bead (t_real *rhs,
    dna_segment *last_seg,
    bead_ode *bead,
    const t_real *last_t,
    const t_real *bead_t,
    const t_real ks,
    const t_real kd);

void segment_H (t_real *Hii,
    t_real *Hij,
    dna_segment *prev_seg,
    dna_segment *next_seg,
    t_real *prev_T,
    t_real *next_T,
    const t_real cfm,
    const t_real kd);

void bead_H (t_real *Hii,
    dna_segment *last_seg,
    bead_ode *bead,
    t_real *prev_T,
    t_real *next_T,
    const t_real cfm,
    const t_real kd);

#endif

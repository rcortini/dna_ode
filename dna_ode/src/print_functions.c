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
#include "simcontext.h"


char * which_body (const unsigned int nbodies, sim_body *body_list, unsigned int body_id) {
  unsigned int i;
  for (i=0; i<nbodies; i++) {
    if (body_list [i].odepdb_id == body_id)
      return body_list [i].name;
  }
  return "Unknown";
}



/* prints the molecule velocities */
void print_molecule_velocities (const unsigned int nbodies, sim_body *body_list) {
  unsigned int i;
  t_vector vel, ang_vel;

  /* cycle on all the bodies */
  for (i=0; i<nbodies; i++) {
    body_velocity (vel, *body_list [i].b_id);
    body_angular_velocity (ang_vel, *body_list [i].b_id);
    message ("%s velocity:\n", body_list [i].name, i);
    message ("\t");
    print_vector3d (vel);
    message ("\n");
    message ("\t");
    print_vector3d (ang_vel);
    message ("\n");
  }
}



/* print the force that acts on all bodies */
void print_molecule_forces (const unsigned int nbodies, sim_body *body_list) {
  unsigned int i;
  t_vector force;

  /* print the values of the forces */
  for (i=0; i<nbodies; i++) {
    body_force (force, *body_list [i].b_id);
    message ("%s force: ", body_list [i].name, i);
    print_vector3d (force);
    message ("\n");
  }
}



/* print the torques that act on all bodies */
void print_molecule_torques (const unsigned int nbodies, sim_body *body_list) {
  unsigned int i;
  t_vector torque;

  /* print the torques on the DNA segments*/
  for (i=0; i<nbodies; i++) {
    body_torque (torque, *body_list [i].b_id);
    message ("%s torque: ", body_list [i].name, i);
    print_vector3d (torque);
    message ("\n");
  }
}



/* prints the molecule velocities */
void print_molecule_kinetic_energies (const unsigned int nbodies, sim_body *body_list) {
  unsigned int i;
  t_real k;

  /* cycle on all the bodies */
  for (i=0; i<nbodies; i++) {
    k = body_total_kinetic_energy (*body_list [i].b_id);
    message ("%s kinetic energy: %.3f\n", body_list [i].name, k);
  }
}




/* prints all the information on the body */
void print_body_info (sim_body b) {
  dBodyID *b_id = b.b_id;
  t_vector pos, vel, ang_vel, force, torque;

  message ("%s:\n", b.name);

  message ("\tposition:");
  body_position (pos, *b_id);
  print_vector3d (pos);
  message ("\n");

  message ("\tvelocity:");
  body_velocity (vel, *b_id);
  print_vector3d (vel);
  message ("\n");

  message ("\tangular velocity:");
  body_angular_velocity (ang_vel, *b_id);
  print_vector3d (ang_vel);
  message ("\n");

  message ("\tforce:");
  body_force (force, *b_id);
  print_vector3d (force);
  message ("\n");

  message ("\ttorque:");
  body_torque (torque, *b_id);
  print_vector3d (torque);
  message ("\n");
}



void print_bead_info (simcontext *simcon, bead_ode *bead, sim_body *bead_body) {
  t_vector bead_inertia;
  simparameters *simparams = simcon->simparams;

  body_inertia (bead_inertia, bead->b_id);

  message ("BEAD INFO:\n");
  message ("\tmass: %.3e g    (%f)\n", bead->mass * simcon->simparams->scale_m, bead->mass);
  message ("\tradius: %f nm    (%f)\n", bead->radius * simcon->simparams->scale_l * 1e7, bead->radius);
  message ("\tinertia: %f %f %f\n", bead_inertia[0], bead_inertia[1], bead_inertia[2]);
  message ("\ttranslational friction = %.3f (%.3e g s^-1)\n",
      bead_body->trans_friction,
      bead_body->trans_friction * simparams->scale_m / simparams->scale_t);
  message ("\ttranslational sigma = %.3f (%.3e dyn)\n",
      bead_body->trans_sigma,
      bead_body->trans_sigma * simparams->scale_f);
  message ("\trotational sigma = (%.3f %.3f %.3f) (%.3e %.3e %.3e) erg\n",
      bead_body->rot_sigma_x,
      bead_body->rot_sigma_y,
      bead_body->rot_sigma_z,
      bead_body->rot_sigma_x * simparams->scale_e,
      bead_body->rot_sigma_y * simparams->scale_e,
      bead_body->rot_sigma_z * simparams->scale_e);
   message ("\trotational friction = (%.3f %.3f %.3f) (%.3e %.3e %.3e) erg s\n",
      bead_body->rot_friction_x,
      bead_body->rot_friction_y,
      bead_body->rot_friction_z,
      bead_body->rot_friction_x * simparams->scale_e * simparams->scale_t,
      bead_body->rot_friction_y * simparams->scale_e * simparams->scale_t,
      bead_body->rot_friction_z * simparams->scale_e * simparams->scale_t);

   /* mode information */
   message ("\tmode: %i\n", bead->mode);
   if (bead->mode == FIXED_N) {
     message ("\tforce: %.2f pN, n: %f\n", 
	 bead->bead_force * simcon->simparams->scale_f * 1e7,
	 bead->bead_turns);
   }
   else if (bead->mode == FIXED_TORQUE) {
     message ("\tforce: %.2f pN, torque: %.3f\n", 
	 bead->bead_force * simcon->simparams->scale_f * 1e7,
	 bead->bead_torque * simcon->simparams->scale_e * 1e14);
   }
   else if (bead->mode == FIXED_EXT) {
     message ("\trho: %.3f\n", 
	 bead->rho);
     message ("\ttrap_k: %.3f (%.3f pN/nm), trap_x0: %.2f, trap_y0: %.2f, trap_z0: %.2f\n",
	 bead->bead_pos_k, bead->bead_pos_k * simcon->simparams->scale_f / simcon->simparams->scale_l,
	 bead->trap_x0,
	 bead->trap_y0,
	 bead->trap_z0);
   }
}



void print_dna_info (simcontext *simcon, dna_segment *seg, sim_body *seg_body) {
  t_vector seg_inertia;
  simparameters *simparams = simcon->simparams;

  body_inertia (seg_inertia, seg->b_id);

  message ("DNA INFO:\n");
  message ("\tsegment number: %d\n", simparams->dna_nsegments);
  message ("\tlength: %.3f nm (%.2f)\n", simparams->dna_length * simparams->scale_l *1e7, simparams->dna_length);
  message ("DNA SEGMENT INFO:\n");
  message ("\tmass: %.3e g    (%f)\n", seg->mass * simparams->scale_m, seg->mass);
  message ("\tradius: %f nm    (%f)\n", seg->radius * simcon->simparams->scale_l * 1e7, seg->radius);
  message ("\tinertia: %f %f %f\n", seg_inertia[0], seg_inertia[1], seg_inertia[2]);
  message ("\ttranslational friction = %.3f (%.3e g s^-1)\n",
      seg_body->trans_friction,
      seg_body->trans_friction * simparams->scale_m / simparams->scale_t);
  message ("\ttranslational sigma = %.3f (%.3e dyn)\n",
      seg_body->trans_sigma,
      seg_body->trans_sigma * simparams->scale_f);
  message ("\trotational sigma = (%.3f %.3f %.3f) (%.3e %.3e %.3e) erg\n",
      seg_body->rot_sigma_x,
      seg_body->rot_sigma_y,
      seg_body->rot_sigma_z,
      seg_body->rot_sigma_x * simparams->scale_e,
      seg_body->rot_sigma_y * simparams->scale_e,
      seg_body->rot_sigma_z * simparams->scale_e);
   message ("\trotational friction = (%.3f %.3f %.3f) (%.3e %.3e %.3e) erg s\n",
      seg_body->rot_friction_x,
      seg_body->rot_friction_y,
      seg_body->rot_friction_z,
      seg_body->rot_friction_x * simparams->scale_e * simparams->scale_t,
      seg_body->rot_friction_y * simparams->scale_e * simparams->scale_t,
      seg_body->rot_friction_z * simparams->scale_e * simparams->scale_t);
   message ("\tgb: %f, gt: %f\n", simparams->gb, simparams->gt);
}

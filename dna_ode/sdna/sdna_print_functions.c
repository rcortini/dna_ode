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



/*
 * UTILITY FUNCTIONS FOR PRINTING OUTPUT 
 */



/* once every print_interval number of steps, this function is called and output is shown on screen */
void sdna_print_data (simcontext * simcon, sdna_simcontext *sdna_simcon) {
  t_vector last_joint_position;

  /* get last joint position */
  body_point_position (last_joint_position, sdna_simcon->dna->seg [sdna_simcon->dna->nsegments-1]->b_id, 0., 0., 0.5 * simcon->simparams->dna_segment_length);

  sdna_simcon->Wr = writhe_dna (sdna_simcon->dna, 0);

  if (sdna_simcon->bead->mode== FIXED_N) {
    message ("step %lu:    z/L: %.3f   Lk: %.3f    Tw: %.3f   Wr: %.3f  dLk: %.3f    Et: %.2f  Ut: %.2f  Kt: %.2f\n",
	simcon->sim_tick,
	last_joint_position[2]/simcon->simparams->dna_length,
	sdna_simcon->Lk,
	sdna_simcon->Tw,
	sdna_simcon->Wr,
	sdna_simcon->Wr+sdna_simcon->Tw-sdna_simcon->Lk,
	sdna_simcon->Et,
	sdna_simcon->Ut,
	sdna_simcon->Kt);
  }
  else if (sdna_simcon->bead->mode== FIXED_TORQUE) {
    message ("step %lu:    z/L: %.3f   Lk: %.3f    Tw: %.3f   Wr: %.3f     Et: %.2f  Ut: %.2f  Kt: %.2f\n",
	simcon->sim_tick,
	last_joint_position[2]/simcon->simparams->dna_length,
	sdna_simcon->Lk,
	sdna_simcon->Tw,
	sdna_simcon->Wr,
	sdna_simcon->Et,
	sdna_simcon->Ut,
	sdna_simcon->Kt);
  }
  else if (sdna_simcon->bead->mode== FIXED_EXT) {
    message ("step %lu:    z/L: %.3f   Lk: %.3f    Tw: %.3f   Wr: %.3f     Et: %.2f  Ut: %.2f  Kt: %.2f\n",
	simcon->sim_tick,
	last_joint_position[2]/simcon->simparams->dna_length,
	sdna_simcon->Lk,
	sdna_simcon->Tw,
	sdna_simcon->Wr,
	sdna_simcon->Et,
	sdna_simcon->Ut,
	sdna_simcon->Kt);
  }
}



/* prints the data of the joints */
void print_joints_data (sdna_simcontext *sdna_simcon) {
  unsigned int i;
  int joint_type;
  dVector3 joint_position;

  /* prints the values of the position of the DNA joints */
  for (i=0; i<sdna_simcon->dna->njoints; i++) {
    dJointID joint = sdna_simcon->dna->joint [i]->j_id;
    message ("Joint %d: type = ", i);
    joint_type = dJointGetType (joint);
    if (joint_type == dJointTypeBall) {
      message ("ball\n");
      dJointGetBallAnchor (joint, joint_position);
      message ("\tanchor 1 = %.16e %.16e %.16e\n", joint_position [0],
	  joint_position [1],
	  joint_position [2]);
      dJointGetBallAnchor2 (joint, joint_position);
      message ("\tanchor 2 = %.16e %.16e %.16e\n", joint_position [0],
	  joint_position [1],
	  joint_position [2]);
      message ("\tjoint ERP parameter = %f\n", dJointGetBallParam (joint, dParamERP));
    }
    else if (joint_type == dJointTypeHinge) {
      message ("hinge\n");
      dJointGetHingeAnchor (joint, joint_position);
      message ("\tanchor 1 = %.16e %.16e %.16e\n", joint_position [0],
	  joint_position [1],
	  joint_position [2]);
      dJointGetHingeAnchor2 (joint, joint_position);
      message ("\tanchor 2 = %.16e %.16e %.16e\n", joint_position [0],
	  joint_position [1],
	  joint_position [2]);
    }
    else
      message ("unrecognized\n");

    /* prints the joint feedback */
    if (i!=0) {
      dJointFeedback jf_id = sdna_simcon->dna->joint [i]->jf_id;
      message ("\tfeedback (f1) = %.16e %.16e %.16e\n", jf_id.f1 [0], jf_id.f1 [1], jf_id.f1 [2]);
      message ("\tfeedback (t1) = %.16e %.16e %.16e\n", jf_id.t1 [0], jf_id.t1 [1], jf_id.t1 [2]);
      message ("\tfeedback (f2) = %.16e %.16e %.16e\n", jf_id.f2 [0], jf_id.f2 [1], jf_id.f2 [2]);
      message ("\tfeedback (t2) = %.16e %.16e %.16e\n", jf_id.t2 [0], jf_id.t2 [1], jf_id.t2 [2]);
    }
  }
}

#include "dna_ode_core.h"
#include "dna_ode.h"

/* this is the collision callback function, which is used to detect collisions and assign contact joints
 * to the collision points */
void collision_callback (void *data, dGeomID o1, dGeomID o2) {
  unsigned int i;
  dBodyID b1, b2;
  simcontext *simcon = (simcontext *) data;

  /* check that the two objects actually exist */
  assert (o1);
  assert (o2);

  /* get the body identities so we can check that they are not connected */
  b1 = dGeomGetBody (o1);
  b2 = dGeomGetBody (o2);

  /* check if the bodies are connected, and if so return, as there are no collisions between connected bodies */
  if(b1 && b2 && dAreConnected (b1,b2))
    return;
  else {
    /* declare the variables used for the collision detection */
    dContact t_contact [MAX_CONTACTS];
    const unsigned int this_numc = dCollide (o1, o2, MAX_CONTACTS, &(t_contact [0].geom), sizeof (dContact));
   
    /* perform collision detection and increments the global number of contacts*/
    simcon->numc += this_numc;

    /* if there is at least one collision, then cycle on all contacts and create contact joints */
    if (this_numc>0) {
      for (i=0; i<this_numc; i++) {
	dJointID contact_joint;

	/* assign all the parameters to the joint */
	t_contact [i].surface.mode = CONTACT_MODE; /* see the constant CONTACT_MODE in simobjects.h */
	/* friction parameters */
	t_contact [i].surface.mu = simcon->simparams->contact_mu;
	t_contact [i].surface.mu2 = simcon->simparams->contact_mu2;
	/* bounce is the amount of "bounciness". */
	t_contact [i].surface.bounce = simcon->simparams->contact_bounce;
	/* bounce_vel is the minimum incoming velocity to cause a bounce */
	t_contact [i].surface.bounce_vel = simcon->simparams->contact_bounce_vel;
	/* constraint force mixing parameter */
	t_contact [i].surface.soft_cfm = simcon->simparams->cfm_glob;
	/* error reduction parameter */
	t_contact [i].surface.soft_erp = simcon->simparams->erp_glob;

	/* put contact_joint in the contactjoints group so we can delete it before the next step. */
	contact_joint = dJointCreateContact (simcon->world_id, simcon->contactjoints, &t_contact[i]);

	/* attach the joint to the bodies */
	dJointAttach (contact_joint, b1, b2);
      }
    }
  }
}


#ifndef DSRG_PT_H
#define DSRG_PT_H

class molecule;
class dsrg_par;

// State-specific, unrelaxed, semicanonical-only spin-free DSRG-PT2 on a converged CAS
// solve. Builds the state-specific generalized Fock from the selected root's 1-RDM,
// semicanonicalizes, then drives the dsrg_sf_tensors algebra core and the batched
// CCVV/CAVV/CCAV terms. Dispatched from single_point_calc beside CDAS_PT2.
int SA_DSRG_PT2(molecule * M, dsrg_par * dsrg, char * job_name);

#endif

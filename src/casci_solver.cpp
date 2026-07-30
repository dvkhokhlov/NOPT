#include "casci_solver.h"

#include <cstdio>
#include <cstdlib>

#include "common_vars.h"   // out_stream

// Out-of-line anchor: gives casci_solver a key function so its vtable/RTTI is emitted in
// exactly one translation unit instead of every includer.
casci_solver::~casci_solver() {}

// Default: the extended-Koopmans IP/EA blocks have no backend-independent construction, and a
// backend reaching here would leave all four output buffers untouched.
int casci_solver::calc_IPEA_single(double *, double *, double *, double *, int, std::vector<double>) {
    fprintf(out_stream, "ERROR: this CI backend does not implement calc_IPEA_single"
                        " (extended-Koopmans IP/EA blocks)\n");
    exit(EXIT_FAILURE);
}

// Default: both shipped backends snapshot their state set for the dressed re-solve overlap, so
// reaching here is a driver mis-dispatch (our own contract) -- abort loudly.
void casci_solver::snapshot_states(int) {
    fprintf(out_stream, "ERROR: this CI backend cannot snapshot its wavefunction set for the"
                        " dressed re-solve overlap (aldet and DMRG/block2 backends only)\n");
    exit(EXIT_FAILURE);
}

// Default: only the DMRG/block2 backend encodes a dressed general MPO. Any other backend reaching
// here is a driver mis-dispatch (our own contract), so abort loudly naming the supported backend.
void casci_solver::import_dressed_operator(const double*, const double*, const double*, double) {
    fprintf(out_stream, "ERROR: this CI backend does not support dressed-operator import"
                        " (DMRG/block2 backend only)\n");
    exit(EXIT_FAILURE);
}

// Defaults: both shipped backends implement the transition-density read-outs; reaching one
// of these is a driver mis-dispatch (our own contract), so abort loudly, never return silence.
void casci_solver::G_calc_full(double*) {
    fprintf(out_stream, "ERROR: this CI backend does not provide the full transition 2-RDM read-out\n");
    exit(EXIT_FAILURE);
}

void casci_solver::h2caa_overlap(const double*, const double*, int, double*) {
    fprintf(out_stream, "ERROR: this CI backend does not provide the complementary-overlap read-out\n");
    exit(EXIT_FAILURE);
}

#include "casci_solver.h"

#include <cstdio>
#include <cstdlib>

#include "common_vars.h"   // out_stream

// Out-of-line anchor: gives casci_solver a key function so its vtable/RTTI is emitted in
// exactly one translation unit instead of every includer.
casci_solver::~casci_solver() {}

// Default: only the DMRG/block2 backend encodes a dressed general MPO. Any other backend reaching
// here is a driver mis-dispatch (our own contract), so abort loudly naming the supported backend.
void casci_solver::import_dressed_operator(const double*, const double*, const double*, double) {
    fprintf(out_stream, "ERROR: this CI backend does not support dressed-operator import"
                        " (DMRG/block2 backend only)\n");
    exit(EXIT_FAILURE);
}

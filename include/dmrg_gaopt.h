#pragma once
//
// dmrg_gaopt — genetic-algorithm orbital ordering, a port of pyblock2's gaopt driver.
// Deterministic: fixed seeds 1234+task, lexicographically smallest of the distinct results.
// kmat is the row-major n_sites x n_sites coupling matrix.

#include <cstdint>
#include <vector>

std::vector<uint16_t> dmrg_gaopt_order(int n_sites, const std::vector<double> &kmat);

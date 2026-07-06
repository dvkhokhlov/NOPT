c4h6 pi CAS(4,4) -- DMRG-CASSCF vs native CASSCF (state-averaged, cc-pvtz).

Active space: the butadiene pi system, state-averaged over 3 of 7 states.
Two DMRG runs -- a cold-start canonical one and the warm-start PM-localized default --
are checked against the native CASSCF reference; the DMRG runs use m=200. The
state-averaged energy (the variational optimization target) reproduces native to the
print floor (~1e-10) in every run. The individual state energies are not variational at
the SA stationary point, so they match ~1.7e-9 cold but scatter ~1.6e-7 under warm-start
-- a benign, path-dependent effect that cancels in the average. Properties
(dipole/quadrupole/Mulliken) also match. Shared geometry/basis/active-space live in
c4h6.inp; the orbital guess is the RHF set c4h6_RHF.orb.

1a) native CAS-SCF reference
    run_sQM -i ms_cas.inp > ms_cas.log
    compare ms_cas.log with benchmark_cas.log

1b) DMRG-CASSCF (block2 backend, m=200)        [needs a USE_BLOCK2=yes build]
    run_sQM -i ms_dmrg.inp > ms_dmrg.log
    compare ms_dmrg.log with benchmark_dmrg.log

1c) PM-localized DMRG-CASSCF (m=200)           [needs a USE_BLOCK2=yes build]
    run_sQM -i ms_dmrg_loc.inp > ms_dmrg_loc.log
    compare ms_dmrg_loc.log with benchmark_dmrg_loc.log

2) try to open *.out with your favorite visualization tool

All three inputs print the full property set ($prop dipole=1 quad=1 mulliken=1).
The DMRG state and transition 1-RDMs come from the block2 backend and reproduce
the native aldet values; transition-dipole/quadrupole signs are gauge (relative CI
phase), so only the magnitudes are physical and compared.

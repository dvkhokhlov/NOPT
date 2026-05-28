# NOPT

Non-Orthogonal quantum chemistry and Perturbation Theory

Quantum-chemistry package for MRPT calculations and fragment based non-orthogonal methods

# Install
## Dependencies

1) MKL or OpenBLAS

2) LIBINT
- use 3c-ERI and 2c-ERI for RI
- use appropriate l_max (NOPT supports up to l_max=7)
add for LIBINT configuration
```
--enable-eri3 --enable-eri2 --with-max-am=7 --with-eri3-max-am=7 --with-eri2-max-am=7
```

3) LIBGRPP for ECP and GRPP integrals

4) Boost and Eigen (only headers, compilation is not required)

## Compilation

Edit Makefile.vars

Instal standard QM module
```
make run_sQM
```

## Docker image

See [dockerbuild/README.md](dockerbuild/README.md) on how to pull the existing image of the NOPT.
 
Or build a native image using [Dockerfile](dockerbuild/Dockerfile)

# How to Use
```
run_sQM -i input.inp
```
See examples/ for more details

# Testing 


# License

NOPT is distributed under the MIT License - see [LICENSE](LICENSE) for details.

Third-party dependencies retain their own licenses;
see [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md) for full attribution.
# NOPT

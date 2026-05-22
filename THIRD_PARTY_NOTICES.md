# Third-Party Notices

NOPT incorporates or links against the third-party software components listed
below. Each component is the property of its respective copyright holders and
is distributed under its own license. The license of NOPT itself (see
`LICENSE`) does not supersede or alter the terms of these licenses.

This file is provided for compliance and attribution purposes. When NOPT is
redistributed in binary form, this file (or an equivalent notice) must be
included.

---

## 1. Libint (Libint2)

- **Description**: High-performance library for the evaluation of molecular
  integrals of many-body operators over Gaussian functions.
- **Authors / Maintainer**: Edward F. Valeev and contributors.
- **Homepage**: <https://github.com/evaleev/libint>
- **License**: GNU Lesser General Public License v3.0 (LGPL-3.0), distributed
  together with the GNU General Public License v3.0 (GPL-3.0), since LGPL-3.0
  is defined as an additional permission on top of GPL-3.0.
  See `licenses/libint/COPYING` and `licenses/libint/COPYING.LESSER`,
  or <https://www.gnu.org/licenses/lgpl-3.0.html>.
- **Citation**: Valeev, E. F. *Libint: A library for the evaluation of
  molecular integrals of many-body operators over Gaussian functions.*
  <https://github.com/evaleev/libint>
- **Use in NOPT**: Evaluation of standard Gaussian integrals - nuclear
  attraction, electron repulsion, overlap, kinetic energy, dipole, and
  quadrupole integrals.
- **Linkage**: Dynamic (preferred) or static.

> Under LGPL-3.0, end users are entitled to replace the Libint library used by
> NOPT with a modified version of Libint. For statically linked builds, the
> corresponding object files of NOPT will be made available upon request to
> enable such replacement, in accordance with Section 4 of LGPL-3.0.

---

## 2. LIBGRPP

- **Description**: Library for the evaluation of molecular integrals of the
  generalized relativistic pseudopotential (GRPP) operator over Gaussian
  functions.
- **Authors**: A. V. Oleynichenko, A. Zaitsevskii, N. S. Mosyagin,
  A. N. Petrov, E. Eliav, A. V. Titov.
- **Homepage**: <https://github.com/aoleynichenko/libgrpp>
- **License**: MIT License.
  See `licenses/libgrpp/LICENSE`.
- **Citation**: Oleynichenko, A. V.; Zaitsevskii, A.; Mosyagin, N. S.;
  Petrov, A. N.; Eliav, E.; Titov, A. V. *LIBGRPP: A Library for the
  Evaluation of Molecular Integrals of the Generalized Relativistic
  Pseudopotential Operator over Gaussian Functions.* Symmetry **2023**,
  *15* (1), 197. DOI: 10.3390/sym15010197.
- **Use in NOPT**: Optional. When enabled at compile time, provides integrals
  for effective core potentials (ECP) and generalized relativistic
  pseudopotentials (GRPP).

---

## 3. OpenBLAS

- **Description**: Optimized BLAS / LAPACK library.
- **Homepage**: <https://github.com/OpenMathLib/OpenBLAS>
- **License**: BSD 3-Clause License.
  See `licenses/openblas/LICENSE`.
- **Use in NOPT**: Default backend for dense linear-algebra operations
  (matrix multiplication, diagonalization, singular value decomposition,
  Cholesky factorization, etc.).

---

## 4. Intel oneAPI Math Kernel Library (oneMKL) — optional

- **Description**: Vendor-optimized BLAS / LAPACK / FFT library from Intel.
- **Homepage**:
  <https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html>
- **License**: Intel Simplified Software License (or, where applicable, the
  Intel oneAPI Base Toolkit license). Proprietary; redistribution and use are
  subject to Intel's terms and export controls.
- **Use in NOPT**: Optional drop-in replacement for OpenBLAS, selected at
  compile time via a build flag. NOPT is not distributed with MKL; end users
  must obtain MKL separately and comply with its license.

---

## 5. Boost C++ Libraries

- **Description**: Peer-reviewed, free C++ libraries.
- **Homepage**: <https://www.boost.org/>
- **License**: Boost Software License 1.0 (BSL-1.0).
  See `licenses/boost/LICENSE_1_0.txt`.
- **Use in NOPT**: Required transitively by Libint; selected utility
  components may be used directly by NOPT.

---

## 6. Eigen

- **Description**: C++ template library for linear algebra: matrices,
  vectors, numerical solvers, and related algorithms.
- **Homepage**: <https://eigen.tuxfamily.org/>
- **License**: Mozilla Public License 2.0 (MPL-2.0).
  See `licenses/eigen/COPYING.MPL2`.
- **Use in NOPT**: Required transitively by Libint; selected utility
  components may be used directly by NOPT.

> Under MPL-2.0, modifications to Eigen source files must remain under
> MPL-2.0 and be made available. NOPT does not modify Eigen source files;
> Eigen is used as an unmodified header-only dependency.

---

## License Compatibility Summary

NOPT is distributed under the MIT License (see `LICENSE`).
All required dependencies are licensed under terms compatible with MIT:

| Component | License        | Compatibility with MIT         |
| --------- | -------------- | ------------------------------ |
| Libint    | LGPL-3.0       | Compatible (linking permitted) |
| LIBGRPP   | MIT            | Compatible                     |
| OpenBLAS  | BSD-3-Clause   | Compatible                     |
| Boost     | BSL-1.0        | Compatible                     |
| Eigen     | MPL-2.0        | Compatible (file-level)        |
| Intel MKL | Proprietary    | Optional; user obligation      |

---

## How to Obtain Source Code of the Above Components

Each component can be obtained from its homepage listed above. For LGPL-3.0
and MPL-2.0 components, the source code corresponding to the version linked
against any released NOPT binary is preserved and available upon request to
the NOPT maintainers.

## Contact

For licensing questions, please contact the NOPT maintainers: glebov_io@phys.chem.msu.ru, glebov_io@mail.ru

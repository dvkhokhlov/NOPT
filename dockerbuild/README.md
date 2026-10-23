# NOPT — Docker

Containerised build of [NOPT](https://github.com/glebov-io/NOPT) with every
dependency pre-compiled. Mount your working directory, run the binary, get
outputs with your own uid/gid.

```
docker run --rm --shm-size=4g -v "$PWD:/data" pupeza/nopt -i input.inp
```

> **Always pass `--shm-size`** (or `--ipc=host`). Docker gives a container only
> **64 MB** of `/dev/shm`, and the DMRG solver keeps its scratch there — see
> [DMRG and `/dev/shm`](#dmrg-and-devshm). The flag is only a cap: it costs
> nothing for calculations that don't use the space.

## Build configuration

| Item              | Value                                                       |
| ----------------- | ----------------------------------------------------------- |
| Base image        | `ubuntu:24.04` (digest-pinned)                              |
| Platforms         | **linux/amd64** and **linux/arm64**                         |
| Compilers         | GCC 13 (`gcc-13`, `g++-13`, `gfortran-13`)                  |
| Linear algebra    | **OpenBLAS**, OpenMP build — amd64: system `libopenblas0-openmp` (0.3.26); arm64: **0.3.34 built from source** (see below). The pthreads build serializes when called from OpenMP regions; DMRG is 1.3–2.1× faster with the OpenMP one |
| BLAS parallelism  | **non-parallel BLAS inside OMP-parallel loops** (`BLAS_PAR=nopt`) — single-threaded matrix multiplies, threading driven by NOPT itself |
| ECP / GRPP        | [libgrpp](https://github.com/aoleynichenko/libgrpp) @ commit `6e63e88f` (statically linked into the binary) |
| Two-electron ints | [libint](https://github.com/evaleev/libint) @ commit `5fb07b48` (release **v2.13.1**, shared library; built with `BUILD_TESTING=OFF`) |
| DMRG (block2)     | [block2-preview](https://github.com/block-hczhai/block2-preview) @ commit `688dc708` (release **v0.5.4**, shared `libblock2.so`; OpenBLAS, `LARGE_BOND`) — enables `cisolver=dmrg`; `-dmrg` variants only |
| libint config     | `LIBINT2_MAX_AM=7`, ERI3/ERI2 engines enabled at 0-th derivative order (energies only) |
| NOPT angular max  | `L_MAX=7`, `RI_L_MAX=7`; `GTO_MAX` = 20, 27 (default) or 30 per variant |
| NOPT itself       | upstream HEAD at build time; every published variant is built from the same commit (image label `org.nopt.ref`) |

The runtime stage ships only the shared-library closure of `run_sQM`:
`libstdc++6`, `libgcc-s1`, `libgomp1`, `libgfortran5`, `libopenblas0-openmp`,
plus the libint (and, in `-dmrg` variants, block2) shared libs, the NOPT
basis-set library, and the license files under `/usr/share/doc/nopt`.

## Variants and tags

Every image comes in six variants — three `GTO_MAX` values, each with and
without the DMRG backend:

| Variant      | `GTO_MAX` | DMRG (`cisolver=dmrg`) |
| ------------ | --------- | ---------------------- |
| `gto20`      | 20        | no                     |
| `gto20-dmrg` | 20        | yes                    |
| `gto27`      | 27        | no                     |
| **`gto27-dmrg`** = `latest` | **27** | **yes** |
| `gto30`      | 30        | no                     |
| `gto30-dmrg` | 30        | yes                    |

**`GTO_MAX`** is the largest number of primitives in one contracted shell that
the integral engine (libint `max_nprim`) accepts. A basis with a longer
contraction aborts with
`Assertion ... "the max number of primitive pairs exceeded" failed`.
In the bundled basis library:

| Needs `GTO_MAX` ≥ | Basis sets |
| ----------------- | ---------- |
| 21 | cc-pVQZ / cc-pwCVQZ family (incl. aug-, heavy-aug-, -DK) for Sc–Zn; Roos ANO-DZ/TZ |
| 22 | cc-pV5Z for Ca, Ga–Kr; 7ZaPa for Ar |
| 27 | cc-pV5Z / cc-pwCV5Z family for Sc–Zn — the longest contraction in `NOPT_lib` |

The cost is memory: libint keeps a per-thread buffer of `GTO_MAX⁴` primitive
quartets **regardless of the basis actually used**. Measured on the same job
(4ene, cc-pVDZ, 16 threads) — energies and run time are identical:

| | `GTO_MAX=20` | `GTO_MAX=30` |
| --- | --- | --- |
| with RI | 476 MB | 807 MB |
| without RI | 2.6 GB | 12.2 GB |

(≈ 150 MB / 510 MB / 780 MB per thread for 20 / 27 / 30 without RI.) Use
`gto20` on small machines with many threads and no RI, `gto27` (the default)
to run every bundled basis, `gto30` for custom basis sets with up to 30
primitives.

**DMRG** variants link [block2](https://github.com/block-hczhai/block2-preview),
which is **GPL-3.0** — see [License](#license). The variants without it carry
only permissively/LGPL-licensed components.

### Pre-built images

Just pull `latest` (or any `<variant>`) — it automatically gives you the
fastest build your CPU can run, on both x86-64 and ARM:

```
docker pull pupeza/nopt                 # = gto27-dmrg
docker pull pupeza/nopt:gto20           # a specific variant, same auto-selection
docker run --rm --shm-size=4g -v "$PWD:/data" pupeza/nopt -i input.inp
```

To pin one CPU level (reproducible benchmarks, air-gapped mirrors), use
`<variant>-<arch>`, e.g. `pupeza/nopt:gto27-dmrg-x86-64-v3`, where `<arch>` is
one of `x86-64-v1` (any x86-64), `x86-64-v2` (SSE4.2), `x86-64-v3` (AVX2),
`x86-64-v4` (AVX-512), `arm64-v8-a` (any AArch64), `arm64-v8.2-a`,
`arm64-v9-a` (SVE), `apple-m1` (Apple M1–M3), `apple-m4` (M4 and newer).
The bare `<arch>` tags are the same images as `gto27-dmrg-<arch>`.
How the multi-arch tags are put together: [`publish/README.md`](publish/README.md).

## Usage

```
docker run --rm --shm-size=4g -v "$PWD:/data" pupeza/nopt -i input.inp
```

The entrypoint inspects the ownership of `/data` and drops privileges to that
uid/gid via `gosu` before execing `run_sQM`, so any files written into the
mount inherit the caller's ownership instead of root's.

Override the detected uid/gid with `-e HOST_UID=$(id -u) -e HOST_GID=$(id -g)`
when the mount is root-owned (typical on shared HPC scratch).

Everything an input references — including orbital files loaded with
`_MO INP=<file>.orb` — must be in the mounted directory.

### DMRG and `/dev/shm`

The DMRG backend (`cisolver=dmrg`, `-dmrg` variants) keeps its scratch —
renormalized operators and MPS tensors — in **`/dev/shm`** by default, and
Docker caps a container's `/dev/shm` at **64 MB**. That is enough only for toy
active spaces: CAS(8,8) at `m=200` fits, CAS(16,16) at `m=300` already aborts
with

```
terminate called after throwing an instance of 'std::runtime_error'
  what():  DataFrame::save_data on '/dev/shm/nopt_dmrg_.../F0.PART.DMRG.LEFT.6' failed.
```

Pick one:

- **`--shm-size=<N>g`** — raise the container's limit. `/dev/shm` is RAM
  (tmpfs), counted against the container's memory. The need grows with the
  active space and `m`: CAS(16,16) at `m=300` (`examples/dmrg/8ene`) peaks at
  ~150 MB of `/dev/shm` (~5 GB for the whole container), so `2g` is plenty for
  that size; raise it for larger active spaces / bond dimensions.
  **Docker Desktop (macOS/Windows):** the limit must fit inside the VM's
  memory (Settings → Resources; default is often 8 GB), together with the
  calculation itself — `docker info --format '{{.MemTotal}}'` shows it.
- **`--ipc=host`** — use the host's `/dev/shm` (usually half of the host RAM).
- **`save_dir=/data`** in `$DMRG` — put the scratch on the mounted disk instead
  of RAM: no size limit, slower I/O.

```
docker run --rm --shm-size=4g -v "$PWD:/data" pupeza/nopt -i dmrg.inp
```

See `examples/dmrg/DMRG_README.md` for the `$DMRG` keywords.

### Apple Silicon (arm64)

For Apple M-series use the named `apple-m1` target rather than `native`: it pins
the ARMv8.4-A ISA floor shared by every Apple Silicon chip (M1–M4) with a wide
big-core tune (`neoverse-n1`), so the result is portable across all of them
(`apple-m4` raises the floor to ARMv8.7-A for M4 and newer). A named target is
also the correct choice when **cross-building on an x86 host via emulation** —
there `native` would detect the QEMU CPU, not the Apple core.

```
# native on the Mac itself (fast):
docker build --build-arg ARCH=apple-m1 -t nopt:apple -f Dockerfile .

# or cross-build from an x86 host via QEMU (slow; install binfmt emulation once):
docker run --privileged --rm tonistiigi/binfmt --install arm64
docker buildx build --platform linux/arm64 --build-arg ARCH=apple-m1 \
    -t nopt:apple --load -f Dockerfile .
```

GCC 13 has no Apple core model (`-mcpu=apple-m1` is rejected); bump `GCC_SERIES`
to a release that adds it and switch the flag for exact tuning.

## Building from source

The published images are deliberately portable (they exclude any instructions
newer than the chosen microarchitecture level). On the host where you actually
run calculations you can do better.

### OpenBLAS on arm64

Ubuntu 24.04's OpenBLAS (0.3.26) is fine on x86-64, but on arm64 under Docker
Desktop (Apple Silicon) block2's many small concurrent GEMMs serialize inside it:
DMRG ran 2–28× slower than with OpenBLAS 0.3.32/0.3.34, and CDAS-PT2 on a DMRG
reference took ~900 s instead of ~45 s on an M5. arm64 images therefore build
OpenBLAS `OPENBLAS_VERSION` (OpenMP, `DYNAMIC_ARCH`) from source into
`/opt/openblas/lib`, which comes first on `LD_LIBRARY_PATH` (same soname as the
package). On x86-64 (Zen 5) the newer releases bring nothing (0.92–1.04×), so
amd64 keeps the package.

### Building under Docker Desktop (macOS / Windows)

The build runs inside Docker Desktop's VM, whose memory (Settings → Resources,
often 8 GB) is all the compilers get. NOPT's heaviest sources need ~2.3 GB of
RAM per compiler process, so the Dockerfile bounds `make -j` by the VM's
available memory (`JOBS=auto`). If a compile is still killed
(`g++: fatal error: Killed signal terminated program cc1plus` /
`ResourceExhausted ... cannot allocate memory`), give the VM more memory or
force fewer jobs: `--build-arg JOBS=1`.

### Native (best performance)

```
docker build --build-arg ARCH=native -t nopt:native -f Dockerfile .
docker run --rm --shm-size=4g -v "$PWD:/data" nopt:native -i input.inp
```

`ARCH=native` is the default, so the `--build-arg` is optional. It compiles with
`-march=native -mtune=native`, letting GCC use every instruction the build CPU
supports — the resulting image is **not portable** between CPU models, so build
it on the machine you'll run on (or one with the same CPU family). Add
`--build-arg GTO_MAX=20` / `--build-arg DMRG=false` to pick a variant.

### Build arguments

| Arg              | Default                                       | Notes                                            |
| ---------------- | --------------------------------------------- | ------------------------------------------------ |
| `ARCH`           | `native`                                      | `native`; amd64: `x86-64-v1..v4`; arm64: `arm64-v8-a`, `arm64-v8.2-a`, `arm64-v9-a`, `apple-m1` (alias `apple-silicon`), `apple-m4` |
| `GTO_MAX`        | `27`                                          | max primitives per contracted shell (see [Variants](#variants-and-tags)) |
| `DMRG`           | `true`                                        | build block2 and link NOPT with `USE_BLOCK2=yes`; `false` = no DMRG |
| `NOPT_GIT_REF`   | *(empty = upstream HEAD)*                     | NOPT commit/branch to build                      |
| `OPENBLAS_FROM_SOURCE` | `auto`                                   | `auto` = build OpenBLAS from source on arm64 only; `yes` / `no` force it |
| `OPENBLAS_VERSION` | `0.3.34`                                    | OpenBLAS release built from source |
| `JOBS`           | `auto`                                        | parallel compile jobs; `auto` = bounded by available RAM (NOPT ~2.5 GB, block2 ~3 GB, libint ~0.5 GB per job) and `nproc` |
| `LIBINT_MAX_AM`  | `7`                                           | libint angular-momentum cap                      |
| `GCC_SERIES`     | `13`                                          | major version of `gcc-N` / `g++-N` / `gfortran-N` |
| `FIXED`          | `true`                                        | `true` = pinned Ubuntu digest + pinned libint/libgrpp/block2 commits; `false` = floating Ubuntu tag + upstream HEADs |
| `UBUNTU_BASE`    | `ubuntu:24.04@sha256:c4a8d5503dfb…`           | full base image ref                              |

Targeting an arm64 microarch builds for that ISA; pick the one matching your
CPU. `BLAS_DIR` is derived from the compiler triplet, so the same Dockerfile
links OpenBLAS correctly on both x86-64 and AArch64.

### Publishing (maintainers)

The scripts that build and push every variant live in [`publish/`](publish/README.md).

## License

NOPT is MIT-licensed — see [LICENSE](../LICENSE). Third-party dependencies retain
their own licenses; see [THIRD_PARTY_NOTICES.md](../THIRD_PARTY_NOTICES.md). The
images carry these files under `/usr/share/doc/nopt/`.

The **`-dmrg` variants** additionally link **block2, licensed under GPL-3.0**
(its license ships as `/opt/block2/lib/LICENSE.block2`). Redistributing those
images or the `run_sQM` binary from them is subject to the GPL-3.0 terms for
the combined work; block2's source is available at
<https://github.com/block-hczhai/block2-preview> (the pinned commit is listed
above). The variants without `-dmrg` contain no GPL-only components.

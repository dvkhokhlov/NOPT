# NOPT — Docker

Containerised build of [NOPT](https://github.com/glebov-io/NOPT) with every
dependency pre-compiled. Mount your working directory, run the binary, get
outputs with your own uid/gid.

## Build configuration

| Item              | Value                                                       |
| ----------------- | ----------------------------------------------------------- |
| Base image        | `ubuntu:24.04` (digest-pinned)                              |
| Compilers         | GCC 13 (`gcc-13`, `g++-13`, `gfortran-13`)                  |
| Linear algebra    | **OpenBLAS** (system `libopenblas0-pthread`)                |
| BLAS parallelism  | **non-parallel BLAS inside OMP-parallel loops** (`BLAS_PAR=nopt`) — single-threaded matrix multiplies, threading driven by NOPT itself |
| ECP / GRPP        | [libgrpp](https://github.com/aoleynichenko/libgrpp) @ commit `6e63e88f` (statically linked into the binary) |
| Two-electron ints | [libint](https://github.com/evaleev/libint) @ commit `5fb07b48` (release **v2.13.1**, shared library) |
| libint config     | `LIBINT2_MAX_AM=7`, ERI3/ERI2 engines enabled at 0-th derivative order (energies only) |
| NOPT angular max  | `L_MAX=7`, `RI_L_MAX=7`, `GTO_MAX=20`                       |
| NOPT itself       | cloned at `master` HEAD on every build (not pinned)         |

The runtime stage ships only the shared-library closure of `run_sQM`:
`libstdc++6`, `libgcc-s1`, `libgomp1`, `libgfortran5`, `libopenblas0-pthread`,
plus the libint shared libs and the NOPT basis-set library. Total image size is
a small fraction of the builder.

## Pre-built images

Pull whichever microarchitecture level matches your CPU:

```
docker pull pupeza/nopt:x86-64-v1     # baseline x86_64 (any AMD64 CPU)
docker pull pupeza/nopt:x86-64-v2     # Nehalem+ (SSE4.2) — covers most hardware in service today
docker pull pupeza/nopt:x86-64-v3     # Haswell+ (AVX2)
docker pull pupeza/nopt:x86-64-v4     # AVX-512 capable CPUs (TO BE ADDED)
docker pull pupeza/nopt:latest        # alias for x86-64-v2
```

`latest` deliberately points to the `v2` build so a blind `docker pull` works
on essentially any x86 server from the last decade. For best throughput, pick
the highest level your CPU supports — `cat /proc/cpuinfo | grep -m1 flags`
plus the [psABI level table](https://en.wikipedia.org/wiki/X86-64#Microarchitecture_levels)
will tell you.

## Usage

```
docker run --rm -v "$PWD:/data" pupeza/nopt:x86-64-v3 -i input.inp
```

The entrypoint inspects the ownership of `/data` and drops privileges to that
uid/gid via `gosu` before execing `run_sQM`, so any files written into the
mount inherit the caller's ownership instead of root's.

Override the detected uid/gid with `-e HOST_UID=$(id -u) -e HOST_GID=$(id -g)`
when the mount is root-owned (typical on shared HPC scratch).

## Native build (best performance)

The published images are deliberately portable, so they exclude any
instructions newer than the chosen microarchitecture level. On the host where
you actually run calculations, compiling with `ARCH=native` lets GCC use every
instruction set the CPU supports (AVX-512-VNNI, AMX, etc. where present) and
typically yields measurable speedups in the integral and CI kernels.

```
curl -O https://raw.githubusercontent.com/glebov-io/NOPT/main/dockerbuild/Dockerfile
docker build --build-arg ARCH=native -t nopt:native -f Dockerfile .
docker run --rm -v "$PWD:/data" nopt:native -i input.inp
```

`ARCH=native` is the default, so the `--build-arg` is optional; it is shown
here for clarity. The resulting image is **not portable** between CPU models —
build it on the same machine you'll run on (or one with the same CPU family).

## Build arguments

| Arg              | Default                                       | Notes                                            |
| ---------------- | --------------------------------------------- | ------------------------------------------------ |
| `ARCH`           | `native`                                      | `native`, `x86-64-v1..v4`, `arm64-v8-a/-v8.2-a/-v9-a` (ARM untested) |
| `LIBINT_MAX_AM`  | `7`                                           | libint angular-momentum cap                      |
| `GCC_SERIES`     | `13`                                          | major version of `gcc-N` / `g++-N` / `gfortran-N` |
| `FIXED`          | `true`                                        | `true` = pinned Ubuntu digest + pinned libint/libgrpp commits; `false` = floating Ubuntu tag + upstream HEADs |
| `UBUNTU_BASE`    | `ubuntu:24.04@sha256:c4a8d5503dfb…`           | full base image ref                              |

## License

NOPT is MIT-licensed — see [LICENSE](../LICENSE). Third-party dependencies
retain their own licenses; see
[THIRD_PARTY_NOTICES.md](../THIRD_PARTY_NOTICES.md).

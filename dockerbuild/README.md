# NOPT — Docker

Containerised build of [NOPT](https://github.com/glebov-io/NOPT) with every
dependency pre-compiled. Mount your working directory, run the binary, get
outputs with your own uid/gid.

## Build configuration

| Item              | Value                                                       |
| ----------------- | ----------------------------------------------------------- |
| Base image        | `ubuntu:24.04` (digest-pinned)                              |
| Platforms         | **linux/amd64** and **linux/arm64**                         |
| Compilers         | GCC 13 (`gcc-13`, `g++-13`, `gfortran-13`)                  |
| Linear algebra    | **OpenBLAS** (system `libopenblas0-pthread`)                |
| BLAS parallelism  | **non-parallel BLAS inside OMP-parallel loops** (`BLAS_PAR=nopt`) — single-threaded matrix multiplies, threading driven by NOPT itself |
| ECP / GRPP        | [libgrpp](https://github.com/aoleynichenko/libgrpp) @ commit `6e63e88f` (statically linked into the binary) |
| Two-electron ints | [libint](https://github.com/evaleev/libint) @ commit `5fb07b48` (release **v2.13.1**, shared library; built with `BUILD_TESTING=OFF`) |
| libint config     | `LIBINT2_MAX_AM=7`, ERI3/ERI2 engines enabled at 0-th derivative order (energies only) |
| NOPT angular max  | `L_MAX=7`, `RI_L_MAX=7`, `GTO_MAX=20`                       |
| NOPT itself       | cloned at `master` HEAD on every build (not pinned)         |

The runtime stage ships only the shared-library closure of `run_sQM`:
`libstdc++6`, `libgcc-s1`, `libgomp1`, `libgfortran5`, `libopenblas0-pthread`,
plus the libint shared libs and the NOPT basis-set library. Total image size is
a small fraction of the builder.

## Building from source

The published images are deliberately portable (they exclude any instructions
newer than the chosen microarchitecture level). On the host where you actually
run calculations you can do better.

### Native (best performance)

```
docker build --build-arg ARCH=native -t nopt:native -f Dockerfile .
docker run --rm -v "$PWD:/data" nopt:native -i input.inp
```

`ARCH=native` is the default, so the `--build-arg` is optional. It compiles with
`-march=native -mtune=native`, letting GCC use every instruction the build CPU
supports — the resulting image is **not portable** between CPU models, so build
it on the machine you'll run on (or one with the same CPU family).

## Pre-built images

Just pull `latest` — it automatically gives you the fastest build your CPU can
run, on both x86-64 and ARM:

```
docker pull pupeza/nopt          # picks the optimal build for the host CPU
docker run --rm -v "$PWD:/data" pupeza/nopt -i input.inp
```

`latest` is a **manifest list**, not a single image. Each microarchitecture is
compiled separately (different `-march`/`-mtune`); the list routes your pull to
the right one:

| Host | What you get |
| ---- | ------------ |
| **amd64**, containerd-based runtime (Kubernetes, recent Docker w/ containerd image store, `nerdctl`) | the **slim** image whose amd64 microarch variant (`v2`/`v3`/`v4`) matches the CPU — smallest download, exact fit |
| **amd64**, anything that ignores amd64 variants (older Docker, classic image store) | the **`x86-64-dispatch`** image (bundles every amd64 build, `exec`s the right one at container start) |
| **arm64** (Graviton, Ampere, Apple Silicon under Docker Desktop, …) | the **`arm64-dispatch`** image — runtime selection of the best arm64 build |

Either way you end up running the best binary for the machine.

**Why amd64 has microarch *variants* but arm64 does not:** container runtimes
can match the amd64 `v2/v3/v4` variants against the host CPU, but they do **not**
match arm64 sub-levels — arm64 gets a single `linux/arm64` slot. So the arm64
slot is always the dispatch image, which inspects the CPU (`/proc/cpuinfo`
features, or Apple implementer id) and picks among `arm64-v8-a`, `arm64-v8.2-a`,
`arm64-v9-a`, and `apple-m1`, falling back to baseline `arm64-v8-a` if unsure.
Baseline-only amd64 ("v1") CPUs are likewise served by `x86-64-dispatch` (there
is no `v1` *variant* — baseline x86-64 *is* plain amd64).

> arm64 images are newly added. All four build cleanly and are runtime-tested on Apple Silicon (Apple M4) 
> through the arm64-dispatch auto-selection. Other ARM64 hosts (AWS Graviton, Ampere) are built but **not yet runtime-tested**  — 
> they should work, but if auto-selection misbehaves, pin a specific tag or set NOPT_ARCH (see below).

### Pinning a specific build

All per-level tags are published if you want to pin one explicitly (reproducible
benchmarks, air-gapped mirrors, forcing a level):

```
# amd64
docker pull pupeza/nopt:x86-64-v1        # baseline x86_64 (any AMD64 CPU)
docker pull pupeza/nopt:x86-64-v2        # Nehalem+ (SSE4.2)
docker pull pupeza/nopt:x86-64-v3        # Haswell+ (AVX2)
docker pull pupeza/nopt:x86-64-v4        # AVX-512
docker pull pupeza/nopt:x86-64-dispatch  # heavy all-builds image (runtime select)

# arm64
docker pull pupeza/nopt:arm64-v8-a       # baseline armv8-a (any AArch64 CPU)
docker pull pupeza/nopt:arm64-v8.2-a     # armv8.2-a (FP16/DotProd)
docker pull pupeza/nopt:arm64-v9-a       # armv9-a (SVE)
docker pull pupeza/nopt:apple-m1         # Apple M-series (armv8.4-a, neoverse-n1 tune)
docker pull pupeza/nopt:arm64-dispatch   # heavy all-builds image (runtime select)
```

`cat /proc/cpuinfo` plus the
[psABI level table](https://en.wikipedia.org/wiki/X86-64#Microarchitecture_levels)
(amd64) tells you which level your CPU supports. Inside either dispatch image
you can force a level with `-e NOPT_ARCH=<tag>` (e.g. `x86-64-v3`, `arm64-v9-a`)
and add `-e NOPT_DEBUG=1` to log which build was selected.

## Usage

```
docker run --rm -v "$PWD:/data" pupeza/nopt -i input.inp
```

The entrypoint inspects the ownership of `/data` and drops privileges to that
uid/gid via `gosu` before execing `run_sQM`, so any files written into the
mount inherit the caller's ownership instead of root's.

Override the detected uid/gid with `-e HOST_UID=$(id -u) -e HOST_GID=$(id -g)`
when the mount is root-owned (typical on shared HPC scratch).

### Apple Silicon (arm64)

For Apple M-series use the named `apple-m1` target rather than `native`: it pins
the ARMv8.4-A ISA floor shared by every Apple Silicon chip (M1–M4) with a wide
big-core tune (`neoverse-n1`), so the result is portable across all of them. A
named target is also the correct choice when **cross-building on an x86 host via
emulation** — there `native` would detect the QEMU CPU, not the Apple core.

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

## Build arguments

| Arg              | Default                                       | Notes                                            |
| ---------------- | --------------------------------------------- | ------------------------------------------------ |
| `ARCH`           | `native`                                      | `native`; amd64: `x86-64-v1..v4`; arm64: `arm64-v8-a`, `arm64-v8.2-a`, `arm64-v9-a`, `apple-m1` (alias `apple-silicon`) |
| `LIBINT_MAX_AM`  | `7`                                           | libint angular-momentum cap                      |
| `GCC_SERIES`     | `13`                                          | major version of `gcc-N` / `g++-N` / `gfortran-N` |
| `FIXED`          | `true`                                        | `true` = pinned Ubuntu digest + pinned libint/libgrpp commits; `false` = floating Ubuntu tag + upstream HEADs |
| `UBUNTU_BASE`    | `ubuntu:24.04@sha256:c4a8d5503dfb…`           | full base image ref                              |

Targeting an arm64 microarch builds for that ISA; pick the one matching your
CPU. `BLAS_DIR` is derived from the compiler triplet, so the same Dockerfile
links OpenBLAS correctly on both x86-64 and AArch64.

## License

NOPT is MIT-licensed — see [LICENSE](../LICENSE). Third-party dependencies retain their own licenses; see
[THIRD_PARTY_NOTICES.md](../THIRD_PARTY_NOTICES.md).

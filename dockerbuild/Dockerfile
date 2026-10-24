# syntax=docker/dockerfile:1.7
#
# NOPT — multi-stage build
#
# Build args (override with --build-arg):
#   ARCH            target microarchitecture (see table below); default: native
#   LIBINT_MAX_AM   libint angular-momentum max (and ERI3/ERI2); default: 7
#   GCC_SERIES      major compiler series (gcc-N, g++-N, gfortran-N). Default: 13.
#   DMRG            "true" (default) builds block2 (libblock2.so) and links NOPT
#                   with USE_BLOCK2=yes, enabling cisolver=dmrg. "false" skips it.
#   GTO_MAX         max primitives per contracted shell (libint max_nprim);
#                   default: 27, the longest contraction in NOPT_lib (5Z sets of
#                   Sc-Zn). Each thread's libint buffer grows as GTO_MAX^4
#                   (~150 MB at 20, ~510 MB at 27, ~780 MB at 30 without RI).
#   NOPT_GIT_REF    NOPT commit/branch to build; default: empty = upstream HEAD.
#   OPENBLAS_FROM_SOURCE  auto (default) | yes | no. auto = build OpenBLAS
#                   OPENBLAS_VERSION from source on arm64 only; amd64 keeps the
#                   Ubuntu package (libopenblas0-openmp). See the openblas stage.
#   OPENBLAS_VERSION OpenBLAS release built from source; default: 0.3.34.
#   JOBS            parallel compile jobs; default: auto = bounded by available
#                   RAM (per-job budgets: libint 0.5 GB, NOPT 2.5 GB, block2
#                   3 GB) and nproc. A Docker Desktop VM (often 8 GB) otherwise
#                   OOM-kills cc1plus at -j$(nproc). Set a number to override.
#   FIXED           "true" (default) additionally pins the Ubuntu base image
#                   digest and the libgrpp/libint git commits. "false" lets
#                   those float to upstream HEAD / latest official tag.
#   UBUNTU_BASE     full base image ref. Default: ubuntu:24.04 pinned by
#                   digest (see ARG below). Override with --build-arg
#                   UBUNTU_BASE=ubuntu:24.04 to use the floating tag.
#
# ARCH values:
#   native        -> -march=native -mtune=native        (binary tied to build host)
#   x86-64-v1     -> -march=x86-64       -mtune=generic (baseline x86_64, ~K8/SSE2)
#   x86-64-v2     -> -march=x86-64-v2    -mtune=generic (~Nehalem, SSE4.2)
#   x86-64-v3     -> -march=x86-64-v3    -mtune=generic (~Haswell, AVX2)
#   x86-64-v4     -> -march=x86-64-v4    -mtune=generic (AVX-512)
#   arm64-v8-a    -> -march=armv8-a      -mtune=generic !!!!!!!!!! ARM architectures were not tested at all
#   arm64-v8.2-a  -> -march=armv8.2-a    -mtune=generic
#   arm64-v9-a    -> -march=armv9-a      -mtune=generic
#   apple-m1      -> -march=armv8.4-a    -mtune=neoverse-n1 (Apple M-series, M1+)
#                    (alias: apple-silicon). GCC 13 has no Apple core model, so
#                    this targets the ARMv8.4-A ISA floor shared by every Apple
#                    Silicon chip with a wide big-core tune as a proxy.
#   apple-m4      -> -march=armv8.7-a    -mtune=neoverse-n1 (Apple M4 and newer;
#                    the published arm64 image picks it on Apple CPUs
#                    reporting SME or the ARMv8.7 markers afp+rpres)
#
# Reproducibility model
# ---------------------
# Always pinned (major-version level — survives security updates):
#   * Compilers          -> versioned meta-packages: gcc-${GCC_SERIES} etc.
#                           (default GCC_SERIES=13, the default on 24.04 LTS)
#   * Runtime libraries  -> SOVERSION-named packages (libgfortran5,
#                           libopenblas0-openmp, libstdc++6, ...).
#                           Ubuntu LTS guarantees ABI stability within these.
#
# FIXED=true (default) additionally pins:
#   * Ubuntu base image  -> pinned by sha256 digest in UBUNTU_BASE below.
#   * libgrpp / libint / block2 -> pinned to specific git commits (*_GIT_REF).
#
# NOPT itself is not pinned — cloned at HEAD unless NOPT_GIT_REF is given
# (publish/build-all.sh passes the current HEAD so every variant gets the same commit).
#
# With FIXED=false the RUN steps fall back to upstream branch tips. The FROM
# line cannot read FIXED, so unpinning the base image needs an explicit
# UBUNTU_BASE override (publish/build-all.sh handles this automatically).
#
# Why no exact apt version pins: pinning gcc=4:13.2.0-7ubuntu1 breaks the
# moment Ubuntu publishes a security update; the old build is purged from the
# archive. Pinning via package name (gcc-13, libgfortran5) survives that.
#
# Runtime: mount the working directory and run; the entrypoint drops to the
# host user so generated files inherit the caller's uid/gid:
#   docker run --rm --shm-size=4g -v "$PWD:/data" nopt -i input.inp
# (--shm-size: DMRG scratch lives in /dev/shm, which Docker caps at 64 MB.)
#

# ---- base image (read by both FROM lines) ----------------------------------
# Pinned to ubuntu:24.04 as of 2026-05; override to unpin.
ARG UBUNTU_BASE=ubuntu:24.04@sha256:c4a8d5503dfb2a3eb8ab5f807da5bc69a85730fb49b5cfca2330194ebcc41c7b

############################################################
# Stage 0: libint source generation (NATIVE $BUILDPLATFORM)
############################################################
# libint's code generator and export packaging run here on the *build* platform
# (natively), never under QEMU. Why: when cross-building arm64 on an x86 host,
# emulating the generator is slow AND — at the production LIBINT2_MAX_AM=7, whose
# export tree is ~7600 files — `cmake -E tar` pack/unpack under QEMU silently
# drops files (e.g. include/libint2/config2.h.cmake.in), so the downstream
# configure dies with "File ... does not exist. (configure_file)". It works at
# low MAX_AM (few files) and natively, which is exactly the arm64+scale+emulation
# signature. Generating natively avoids tar-under-emulation entirely; libint's
# emitted C++ is architecture-independent and gets the target -march when the
# builder stage compiles it.
FROM --platform=$BUILDPLATFORM ${UBUNTU_BASE} AS libint-gen

ARG FIXED=true
ARG LIBINT_MAX_AM=7
ARG GCC_SERIES=13
ARG LIBINT_GIT_REF=5fb07b4862f219749c51d59ee08073d20a56d506

ENV DEBIAN_FRONTEND=noninteractive
RUN set -eux; \
    apt-get update; \
    apt-get install -y --no-install-recommends \
        gcc-${GCC_SERIES} g++-${GCC_SERIES} \
        git cmake ninja-build libeigen3-dev libboost-dev libgmp-dev ca-certificates; \
    rm -rf /var/lib/apt/lists/*
ENV CC=gcc-${GCC_SERIES}
ENV CXX=g++-${GCC_SERIES}

WORKDIR /gen
# Full clone (no --depth 1): libint's dynamic_version runs `git describe`.
RUN set -eux; \
    git clone https://github.com/evaleev/libint.git; \
    if [ "$FIXED" = "true" ]; then cd libint && git checkout --quiet "$LIBINT_GIT_REF"; fi

# Configure with GENERATION-time options only (MAX_AM / ERI enables are baked
# into the emitted source; -march / shared-vs-static are compile-time and are
# applied in the builder stage). `--target export` builds the generator, runs
# it, and stages the complete standalone library source under build/libint-<ver>/.
# We copy that directory (not the tarball) so no `cmake -E tar` is ever needed.
RUN set -eux; \
    cd libint && cmake -S . -B build -G Ninja \
        -D CMAKE_CXX_COMPILER="${CXX}" \
        -D CMAKE_BUILD_TYPE=Release \
        -D LIBINT2_ENABLE_ERI3=0 \
        -D LIBINT2_ENABLE_ERI2=0 \
        -D LIBINT2_MAX_AM=${LIBINT_MAX_AM} \
        -D LIBINT2_ERI3_MAX_AM=${LIBINT_MAX_AM} \
        -D LIBINT2_ERI2_MAX_AM=${LIBINT_MAX_AM} \
        -D LIBINT2_ENABLE_GENERIC_CODE=ON && \
    cmake --build build --target export -j "$(nproc)"; \
    d=$(find build -maxdepth 1 -type d -name 'libint-*' | head -1); \
    test -n "$d"; \
    test -f "$d/include/libint2/config2.h.cmake.in"; \
    cp -a "$d" /libint-src

############################################################
# Stage 0b: NOPT source (NATIVE $BUILDPLATFORM)
############################################################
# Unpinned by default — HEAD unless NOPT_GIT_REF. A stage of its own so that
# (a) the fetch happens once for all arches and never under QEMU, and (b) the
# NOPT_GIT_REF arg stays out of the builder stage, where it would invalidate
# the libint/block2 layers on every new commit. Shallow fetch of just the
# wanted commit (GitHub serves any reachable SHA) with retries — full clones of
# this repo have flaked mid-transfer (curl 92 / early EOF).
FROM --platform=$BUILDPLATFORM ${UBUNTU_BASE} AS nopt-src
RUN set -eux; \
    apt-get update; \
    apt-get install -y --no-install-recommends git ca-certificates; \
    rm -rf /var/lib/apt/lists/*
ARG NOPT_GIT_REF=
RUN set -eux; \
    for i in 1 2 3 4 5; do \
        rm -rf /NOPT; \
        git init --quiet /NOPT && \
        git -C /NOPT fetch --quiet --depth 1 https://github.com/glebov-io/NOPT.git "${NOPT_GIT_REF:-HEAD}" && \
        git -C /NOPT checkout --quiet FETCH_HEAD && break; \
        echo "NOPT fetch attempt $i failed; retrying" >&2; sleep 20; \
    done; \
    git -C /NOPT rev-parse --verify HEAD

############################################################
# Stage 0c: OpenBLAS from source (arm64 by default)
############################################################
# Ubuntu 24.04's OpenBLAS 0.3.26 is fine on x86-64, but on arm64 under Docker
# Desktop's VM (Apple Silicon) block2's many small concurrent GEMMs serialize in
# it: DMRG ran 2-28x slower than with 0.3.32/0.3.34 (CDAS-PT2 on a DMRG
# reference: ~900 s -> ~45 s on an M5). On x86-64 (Zen 5) newer releases bring
# nothing (0.92-1.04x), so amd64 keeps the package. The library goes to
# /opt/openblas/lib, first on LD_LIBRARY_PATH at runtime: same soname
# (libopenblas.so.0), so NOPT and block2 link against the package as usual.
# OpenMP threading and DYNAMIC_ARCH (runtime kernel choice), like Ubuntu's build.
FROM ${UBUNTU_BASE} AS openblas
ARG TARGETARCH
ARG GCC_SERIES=13
ARG OPENBLAS_FROM_SOURCE=auto
ARG OPENBLAS_VERSION=0.3.34
RUN set -eux; mkdir -p /opt/openblas/lib; \
    case "$OPENBLAS_FROM_SOURCE" in \
        yes) ;; \
        auto) [ "$TARGETARCH" = arm64 ] || exit 0 ;; \
        *) exit 0 ;; \
    esac; \
    apt-get update; \
    apt-get install -y --no-install-recommends \
        gcc-${GCC_SERIES} gfortran-${GCC_SERIES} make curl ca-certificates; \
    rm -rf /var/lib/apt/lists/*; \
    curl -fsSL --retry 5 -o /tmp/ob.tgz \
        "https://github.com/OpenMathLib/OpenBLAS/releases/download/v${OPENBLAS_VERSION}/OpenBLAS-${OPENBLAS_VERSION}.tar.gz"; \
    tar -C /tmp -xzf /tmp/ob.tgz; cd "/tmp/OpenBLAS-${OPENBLAS_VERSION}"; \
    OB="CC=gcc-${GCC_SERIES} FC=gfortran-${GCC_SERIES} USE_OPENMP=1 DYNAMIC_ARCH=1 NO_STATIC=1 NUM_THREADS=128"; \
    make -j"$(nproc)" $OB > /tmp/openblas-build.log 2>&1 || { tail -40 /tmp/openblas-build.log; exit 1; }; \
    make $OB PREFIX=/opt/openblas install > /dev/null; \
    rm -rf /opt/openblas/include /opt/openblas/lib/cmake /opt/openblas/lib/pkgconfig /tmp/OpenBLAS-* /tmp/ob.tgz; \
    ls -l /opt/openblas/lib

############################################################
# Stage 1: builder
############################################################
FROM ${UBUNTU_BASE} AS builder

ARG ARCH=native
ARG FIXED=true
ARG LIBINT_MAX_AM=7

# Compiler major-version series. 13 is the default on 24.04 LTS.
ARG GCC_SERIES=13

# Pinned upstream refs (honoured only when FIXED=true).
ARG LIBGRPP_GIT_REF=6e63e88f75385b811837efbbb143b2dcd83b00e8
ARG LIBINT_GIT_REF=5fb07b4862f219749c51d59ee08073d20a56d506
# libint commit above corresponds to release tag v2.13.1 (2026-02-03).
ARG BLOCK2_GIT_REF=688dc708f172d91655d4674c92c59d492234dce3
# block2 commit above corresponds to release tag v0.5.4 (2026-09-24).

ENV DEBIAN_FRONTEND=noninteractive

# Translate ARCH into the compiler flag string used by every sub-build.
# Written to a file so subsequent RUN layers pick the same value via $(cat).
RUN set -eux; \
    case "${ARCH}" in \
        native)       MARCH="-march=native -mtune=native" ;; \
        x86-64-v1)    MARCH="-march=x86-64       -mtune=generic" ;; \
        x86-64-v2)    MARCH="-march=x86-64-v2    -mtune=generic" ;; \
        x86-64-v3)    MARCH="-march=x86-64-v3    -mtune=generic" ;; \
        x86-64-v4)    MARCH="-march=x86-64-v4    -mtune=generic" ;; \
        arm64-v8-a)   MARCH="-march=armv8-a      -mtune=generic" ;; \
        arm64-v8.2-a) MARCH="-march=armv8.2-a    -mtune=generic" ;; \
        arm64-v9-a)   MARCH="-march=armv9-a      -mtune=generic" ;; \
        apple-m4)     MARCH="-march=armv8.7-a    -mtune=neoverse-n1" ;; \
        apple-m1|apple-silicon) MARCH="-march=armv8.4-a    -mtune=neoverse-n1" ;; \
        *) echo "Unsupported ARCH='${ARCH}'. See Dockerfile header for valid values." >&2; exit 1 ;; \
    esac; \
    mkdir -p /opt/arch; \
    printf '%s\n' "$MARCH" > /opt/arch/march_flags; \
    echo "ARCH=${ARCH}  ->  $MARCH"

RUN set -eux; \
    apt-get update; \
    apt-get install -y --no-install-recommends \
        gcc-${GCC_SERIES} g++-${GCC_SERIES} gfortran-${GCC_SERIES} \
        libc6-dev make \
        git cmake ninja-build \
        libeigen3-dev libboost-dev \
        libopenblas-openmp-dev liblapacke-dev libgmp-dev \
        python3 ca-certificates; \
    rm -rf /var/lib/apt/lists/*

# `build-jobs <MB per job>`: the -j for a compile whose jobs each peak at about
# that much RAM — min(nproc, MemAvailable / MB), at least 1. JOBS=<n> overrides.
# Measured peaks (x86-64, GCC 13): NOPT's heaviest TUs ~2.3 GB, libint's
# largest generated TU ~0.35 GB, block2 instantiation TUs ~3 GB.
ARG JOBS=auto
COPY --chmod=755 <<'EOF' /usr/local/bin/build-jobs
#!/bin/sh
if [ -n "${JOBS:-}" ] && [ "$JOBS" != auto ]; then echo "$JOBS"; exit 0; fi
n=$(awk -v mb="$1" '/MemAvailable/ {print int($2 / 1024 / mb)}' /proc/meminfo)
c=$(nproc)
[ "$n" -ge 1 ] || n=1
[ "$n" -le "$c" ] || n=$c
echo "$n"
EOF

# Compiler names used by every downstream build. Set via ENV so they survive
# across RUN layers. The versioned binaries (gcc-13 etc.) guarantee that
# GCC_SERIES is actually honoured — the plain `gcc`/`g++`/`gfortran` aliases
# from the meta-packages always point at the distro default and ignore it.
ENV CC=gcc-${GCC_SERIES}
ENV CXX=g++-${GCC_SERIES}
ENV FC=gfortran-${GCC_SERIES}

WORKDIR /deps

# ----- libgrpp (static; gets embedded in run_sQM) -----------------------
RUN set -eux; \
    git clone https://github.com/aoleynichenko/libgrpp.git; \
    if [ "$FIXED" = "true" ]; then cd libgrpp && git checkout --quiet "$LIBGRPP_GIT_REF"; fi
RUN set -eux; MARCH=$(cat /opt/arch/march_flags); \
    mkdir libgrpp/build && cd libgrpp/build && \
    cmake .. \
        -DCMAKE_C_COMPILER="${CC}" \
        -DCMAKE_Fortran_COMPILER="${FC}" \
        -DCMAKE_C_FLAGS="${MARCH}" \
        -DCMAKE_Fortran_FLAGS="${MARCH}" && \
    make -j"$(nproc)"

# ----- libint (shared) — source generated natively in the libint-gen stage,
# compiled here with the target -march. Copying the pre-staged source tree
# (not a tarball) keeps `cmake -E tar` out of the emulated build entirely, which
# is what was corrupting the arm64/QEMU build at MAX_AM=7 (see libint-gen).
COPY --from=libint-gen /libint-src /deps/libint-src
RUN set -eux; MARCH=$(cat /opt/arch/march_flags); \
    cmake -S /deps/libint-src -B /deps/libint-src/build -G Ninja \
        -D CMAKE_CXX_COMPILER="${CXX}" \
        -D CMAKE_INSTALL_PREFIX=/opt/libint \
        -D CMAKE_BUILD_TYPE=Release \
        -D CMAKE_CXX_FLAGS="${MARCH}" \
        -D BUILD_SHARED_LIBS=ON \
        -D BUILD_TESTING=OFF \
        -D LIBINT2_REQUIRE_CXX_API=ON \
        -D LIBINT2_REQUIRE_CXX_API_COMPILED=ON && \
    cmake --build /deps/libint-src/build --target install -j "$(build-jobs 512)"

# ----- block2 (shared; DMRG CI solver) ----------------------------------
# Built as the C++ library (BUILD_CLIB) against OpenBLAS; the exported ABI
# macros (_F77UNDERSCORE, _HAS_OPENBLAS, _LARGE_BOND, ...) are what NOPT's
# Makefile compiles its block2 TUs with — see examples/dmrg/DMRG_README.md.
# configure needs a python3 interpreter; runtime does not. With DMRG=false an
# empty /opt/block2/lib is left behind so the runtime COPY stays unconditional.
# Variant args are declared as late as possible: a changed ARG value misses the
# cache from its first RUN on, so DMRG/GTO_MAX variants all reuse the (slow,
# especially under QEMU) libgrpp + libint layers.
ARG DMRG=true
RUN set -eux; \
    mkdir -p /opt/block2/lib; \
    if [ "$DMRG" = "true" ]; then \
        if [ "$FIXED" = "true" ]; then REF="$BLOCK2_GIT_REF"; else REF=HEAD; fi; \
        for i in 1 2 3 4 5; do \
            rm -rf block2; \
            git init --quiet block2 && \
            git -C block2 fetch --quiet --depth 1 https://github.com/block-hczhai/block2-preview.git "$REF" && \
            git -C block2 checkout --quiet FETCH_HEAD && break; \
            echo "block2 fetch attempt $i failed; retrying" >&2; sleep 20; \
        done; \
        git -C block2 rev-parse --verify HEAD; \
    fi
RUN set -eux; MARCH=$(cat /opt/arch/march_flags); \
    if [ "$DMRG" != "true" ]; then exit 0; fi; \
    cmake -S block2 -B block2/build -G Ninja \
        -D CMAKE_BUILD_TYPE=Release \
        -D CMAKE_CXX_FLAGS="${MARCH}" \
        -D CMAKE_INSTALL_PREFIX=/opt/block2 \
        -D BLA_VENDOR=OpenBLAS \
        -D BUILD_LIB=OFF -D BUILD_CLIB=ON \
        -D USE_MKL=OFF -D USE_OPENBLAS=ON -D F77UNDERSCORE=ON \
        -D LARGE_BOND=ON && \
    cmake --build block2/build --target install -j "$(build-jobs 3072)"; \
    rm -rf block2/build
# block2 is GPL-3.0: its license travels with libblock2.so (and so into every
# image built from it). Separate layer to keep the build cached.
RUN set -eux; \
    if [ "$DMRG" = "true" ]; then cp block2/LICENSE /opt/block2/lib/LICENSE.block2; fi

# ----- NOPT (source from the nopt-src stage) ----------------------------
COPY --from=nopt-src /NOPT /deps/NOPT
WORKDIR /deps/NOPT

# NOPT_ARCH_FLAGS is supplied via the build-time environment so the same
# ARCH translation flows through to NOPT's own compilation units.
COPY <<'EOF' /deps/NOPT/Makefile.vars
CXX:=g++
NOPT_lib_DIR=/deps/NOPT/NOPT_lib
EIGEN_DIR:=/usr/include/eigen3
MANUAL_BOOST_PATH:=no
L_MAX:=7
RI_L_MAX:=7
LIBINT_H_DIR:=/opt/libint/include
LIBINT_L_DIR:=/opt/libint/lib
BLAS_LIB:=openblas
BLAS_DIR:=/usr/lib/$(shell $(CXX) -dumpmachine)
BLAS_PAR:=nopt
USE_GRPP:=yes
GRPP_H_DIR:=/deps/libgrpp/libgrpp
GRPP_L_DIR:=/deps/libgrpp/build/libgrpp
BLOCK2_DIR:=/opt/block2
OPT_LEVEL:=-O3 -g $(NOPT_ARCH_FLAGS) -Wno-sign-compare -Wreturn-type -std=c++17 -ftree-vectorize
EOF

ARG GTO_MAX=27
# `make clean` first: a local tree passed with --build-context nopt-src=<dir> may carry objects of
# a native host build (e.g. macOS Mach-O), which make would otherwise take as up to date.
RUN set -eux; MARCH=$(cat /opt/arch/march_flags); \
    if [ "$DMRG" = "true" ]; then USE_BLOCK2=yes; else USE_BLOCK2=no; fi; \
    make clean CXX="${CXX}" >/dev/null; \
    NOPT_ARCH_FLAGS="${MARCH}" make CXX="${CXX}" USE_BLOCK2="${USE_BLOCK2}" \
        GTO_MAX="${GTO_MAX}" run_sQM -j"$(build-jobs 2560)"

############################################################
# Stage 2: runtime (slim — only shared libs the binary actually loads)
############################################################
FROM ${UBUNTU_BASE} AS runtime

ENV DEBIAN_FRONTEND=noninteractive

# Runtime closure of run_sQM. All packages are SOVERSION-named, so the
# package name itself pins the ABI-relevant major version; security updates
# inside that series are picked up freely.
#   libstdc++6, libgcc-s1       — C++/runtime support
#   libgomp1                    — OpenMP
#   libgfortran5                — Fortran runtime (libgfortran linked in)
#   libopenblas0-openmp         — cblas_* used by NOPT and block2 (matches the build side).
#                                 The OpenMP build, not pthreads: called from inside OpenMP
#                                 regions (block2's batched GEMMs, BLAS_PAR=nopt) it runs on the
#                                 calling thread, where the pthreads build serializes on its
#                                 own pool -- DMRG runs 1.3-2.1x faster (x86-64, 8-16 threads).
#   libint2 shared lib          — copied from builder under /opt/libint
#   libblock2 shared lib        — copied from builder under /opt/block2 (DMRG;
#                                 run_sQM carries an rpath to it)
# libgrpp is statically linked into run_sQM, so nothing to ship for it.
RUN set -eux; \
    apt-get update; \
    apt-get install -y --no-install-recommends \
        libstdc++6 libgcc-s1 libgomp1 libgfortran5 libopenblas0-openmp \
        gosu ca-certificates; \
    rm -rf /var/lib/apt/lists/*

# libint shared libraries
COPY --from=builder /opt/libint/lib /opt/libint/lib

# block2 shared library (empty directory when built with DMRG=false)
COPY --from=builder /opt/block2/lib /opt/block2/lib

# NOPT's license and third-party notices (required with binary redistribution).
COPY --from=builder /deps/NOPT/LICENSE /deps/NOPT/THIRD_PARTY_NOTICES.md /usr/share/doc/nopt/
COPY --from=builder /deps/NOPT/licenses /usr/share/doc/nopt/licenses

# NOPT binary + basis-set library (the path /deps/NOPT/NOPT_lib is baked into
# the binary via -DNOPT_LIB at compile time, so it must live here at runtime).
COPY --from=builder /deps/NOPT/run_sQM       /deps/NOPT/run_sQM
COPY --from=builder /deps/NOPT/NOPT_lib      /deps/NOPT/NOPT_lib

# OpenBLAS built from source (arm64 by default; empty directory otherwise), ahead
# of the package's libopenblas.so.0 -- see the openblas stage.
COPY --from=openblas /opt/openblas/lib /opt/openblas/lib

ENV LD_LIBRARY_PATH=/opt/openblas/lib:/opt/libint/lib:/opt/block2/lib

# ----- entrypoint: drop privileges to the owner of /data ----------------
# gosu preserves the environment (incl. LD_LIBRARY_PATH) when switching uid.
COPY <<'EOF' /usr/local/bin/docker-entrypoint.sh
#!/bin/sh
set -eu

# Prefer explicit overrides (HOST_UID / HOST_GID); otherwise inherit the
# ownership of the mounted /data so generated files match the caller.
TARGET_UID="${HOST_UID:-$(stat -c '%u' /data)}"
TARGET_GID="${HOST_GID:-$(stat -c '%g' /data)}"

# Refuse to silently run as root just because /data is root-owned; that is
# almost always a mount-permissions mistake.
if [ "${TARGET_UID}" = "0" ] && [ -z "${HOST_UID:-}" ]; then
    echo "warning: /data is owned by root; pass -e HOST_UID=\$(id -u) -e HOST_GID=\$(id -g) to override" >&2
fi

# Create matching group/user on the fly (idempotent — re-using an existing
# uid/gid is fine, gosu accepts numeric ids too).
if ! getent group "${TARGET_GID}" >/dev/null 2>&1; then
    groupadd -g "${TARGET_GID}" runner
fi
if ! getent passwd "${TARGET_UID}" >/dev/null 2>&1; then
    useradd -M -N -u "${TARGET_UID}" -g "${TARGET_GID}" -s /usr/sbin/nologin runner
fi

exec gosu "${TARGET_UID}:${TARGET_GID}" /deps/NOPT/run_sQM "$@"
EOF
RUN chmod +x /usr/local/bin/docker-entrypoint.sh

# Variant metadata last, so these ARGs don't invalidate the apt layer above.
ARG ARCH=native
ARG DMRG=true
ARG GTO_MAX=27
ARG NOPT_GIT_REF=
LABEL org.nopt.arch="${ARCH}" org.nopt.dmrg="${DMRG}" org.nopt.gto_max="${GTO_MAX}" \
      org.nopt.ref="${NOPT_GIT_REF}"

WORKDIR /data
ENTRYPOINT ["/usr/local/bin/docker-entrypoint.sh"]

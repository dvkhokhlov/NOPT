# syntax=docker/dockerfile:1.7
#
# NOPT — multi-stage build
#
# Build args (override with --build-arg):
#   ARCH            target microarchitecture (see table below); default: native
#   LIBINT_MAX_AM   libint angular-momentum max (and ERI3/ERI2); default: 7
#   UBUNTU_VERSION  base image tag; default: 24.04
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
#
# Runtime: mount the working directory and run; the entrypoint drops to the
# host user so generated files inherit the caller's uid/gid:
#   docker run --rm -v "$PWD:/data" nopt -i input.inp
#
ARG UBUNTU_VERSION=24.04

############################################################
# Stage 1: builder
############################################################
FROM ubuntu:${UBUNTU_VERSION} AS builder

ARG ARCH=native
ARG LIBINT_MAX_AM=7

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
        *) echo "Unsupported ARCH='${ARCH}'. See Dockerfile header for valid values." >&2; exit 1 ;; \
    esac; \
    mkdir -p /opt/arch; \
    printf '%s\n' "$MARCH" > /opt/arch/march_flags; \
    echo "ARCH=${ARCH}  ->  $MARCH"

RUN apt-get update && apt-get install -y --no-install-recommends \
        gcc g++ gfortran build-essential git cmake ninja-build \
        libeigen3-dev libboost-dev \
        libopenblas-dev liblapacke-dev libgmp-dev \
        ca-certificates \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /deps

# ----- libgrpp (static; gets embedded in run_sQM) -----------------------
RUN git clone --depth 1 https://github.com/aoleynichenko/libgrpp.git
RUN set -eux; MARCH=$(cat /opt/arch/march_flags); \
    mkdir libgrpp/build && cd libgrpp/build && \
    CC=gcc FC=gfortran \
    cmake .. \
        -DCMAKE_C_FLAGS="${MARCH}" \
        -DCMAKE_Fortran_FLAGS="${MARCH}" && \
    make -j"$(nproc)"

# ----- libint (shared; copied into runtime stage) -----------------------
# Full clone (no --depth 1): libint's CMakeLists uses dynamic_version, which
# runs `git describe` against tags. A shallow clone strips tags and configure
# fails with "No names found, cannot describe anything."
RUN git clone https://github.com/evaleev/libint.git
RUN set -eux; MARCH=$(cat /opt/arch/march_flags); \
    cd libint && cmake -S . -B build -G Ninja \
        -D CMAKE_CXX_COMPILER=g++ \
        -D CMAKE_INSTALL_PREFIX=/opt/libint \
        -D CMAKE_BUILD_TYPE=Release \
        -D CMAKE_CXX_FLAGS="${MARCH}" \
        -D LIBINT2_ENABLE_ERI3=0 \
        -D LIBINT2_ENABLE_ERI2=0 \
        -D LIBINT2_MAX_AM=${LIBINT_MAX_AM} \
        -D LIBINT2_ERI3_MAX_AM=${LIBINT_MAX_AM} \
        -D LIBINT2_ERI2_MAX_AM=${LIBINT_MAX_AM} \
        -D LIBINT2_ENABLE_GENERIC_CODE=ON \
        -D BUILD_SHARED_LIBS=ON && \
    cmake --build build --target install -j "$(nproc)"

# ----- NOPT --------------------------------------------------------------
RUN git clone https://github.com/glebov-io/NOPT.git
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
GTO_MAX:=20
LIBINT_H_DIR:=/opt/libint/include
LIBINT_L_DIR:=/opt/libint/lib
BLAS_LIB:=openblas
BLAS_DIR:=/usr/lib/x86_64-linux-gnu
BLAS_PAR:=nopt
USE_GRPP:=yes
GRPP_H_DIR:=/deps/libgrpp/libgrpp
GRPP_L_DIR:=/deps/libgrpp/build/libgrpp
OPT_LEVEL:=-O3 -g $(NOPT_ARCH_FLAGS) -Wno-sign-compare -Wreturn-type -std=c++17 -ftree-vectorize
EOF

RUN set -eux; MARCH=$(cat /opt/arch/march_flags); \
    NOPT_ARCH_FLAGS="${MARCH}" make run_sQM -j"$(nproc)"

############################################################
# Stage 2: runtime (slim — only shared libs the binary actually loads)
############################################################
FROM ubuntu:${UBUNTU_VERSION} AS runtime

ENV DEBIAN_FRONTEND=noninteractive

# Runtime closure of run_sQM:
#   libstdc++6, libgcc-s1       — C++/runtime support
#   libgomp1                    — OpenMP
#   libgfortran5                — Fortran runtime (libgfortran linked in)
#   libopenblas0-pthread        — cblas_* used by NOPT (matches build-side variant)
#   libint2 shared lib          — copied from builder under /opt/libint
# libgrpp is statically linked into run_sQM, so nothing to ship for it.
RUN apt-get update && apt-get install -y --no-install-recommends \
        libstdc++6 libgcc-s1 libgomp1 libgfortran5 libopenblas0-pthread \
        gosu ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# libint shared libraries
COPY --from=builder /opt/libint/lib /opt/libint/lib

# NOPT binary + basis-set library (the path /deps/NOPT/NOPT_lib is baked into
# the binary via -DNOPT_LIB at compile time, so it must live here at runtime).
COPY --from=builder /deps/NOPT/run_sQM       /deps/NOPT/run_sQM
COPY --from=builder /deps/NOPT/NOPT_lib      /deps/NOPT/NOPT_lib

ENV LD_LIBRARY_PATH=/opt/libint/lib

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

WORKDIR /data
ENTRYPOINT ["/usr/local/bin/docker-entrypoint.sh"]

#!/usr/bin/env bash
#
# Build every slim per-arch NOPT image (../Dockerfile), for every variant, into the local image
# store (--load; nothing is pushed — build-multiarch.sh does that with PUSH=1).
#
#   variants : GTO_MAX in {20,27,30} x DMRG {with,without}     (see tags.sh)
#   arches   : x86-64-v1..v4 (amd64), arm64-v8-a .. apple-m4 (arm64, via QEMU)
#   tag      : REPO:gto<N>-dmrg-<arch> | REPO:gto<N>-<arch>
#
# Arch is the outer loop so each arch's libgrpp/libint/block2 layers are built
# once and reused by all its variants (the Dockerfile declares the variant args
# late for exactly this). All variants are built from the same NOPT commit:
# NOPT_GIT_REF defaults to upstream HEAD resolved once at start. A failed image
# is retried (RETRIES, default 2 extra attempts) — cheap, as finished layers come
# from cache — to ride out transient network errors.
#
# Usage:
#   ./build-all.sh                          # everything
#   ARM=0 ./build-all.sh                    # amd64 only
#   GTOS=27 DMRGS=with ./build-all.sh       # a subset (see tags.sh for knobs)
#   SKIP_EXISTING=1 ./build-all.sh          # resume: skip tags already built
#                                           # from the same NOPT commit
#   FIXED=false ./build-all.sh              # float base image + upstream HEADs
# Per-image logs go to build-logs/<tag>.log.
set -euo pipefail

cd "$(dirname "$0")"
. ./tags.sh

AMD=${AMD:-1}
ARM=${ARM:-1}
FIXED=${FIXED:-true}
SKIP_EXISTING=${SKIP_EXISTING:-0}
RETRIES=${RETRIES:-2}
LOG_DIR=${LOG_DIR:-build-logs}
NOPT_GIT_REF=${NOPT_GIT_REF:-$(git ls-remote https://github.com/glebov-io/NOPT.git HEAD | cut -f1)}
[ -n "$NOPT_GIT_REF" ] || { echo "could not resolve NOPT HEAD" >&2; exit 1; }

EXTRA=()
if [ "$FIXED" != "true" ]; then
    EXTRA+=(--build-arg FIXED=false --build-arg UBUNTU_BASE=ubuntu:24.04)
fi

mkdir -p "$LOG_DIR"
echo ">> NOPT commit ${NOPT_GIT_REF}"

build_one() {  # <platform> <arch> <gto> <with|without>
    local platform=$1 arch=$2 gto=$3 dmrg=$4 tag log t0 try
    tag="${REPO}:$(arch_tag "$(variant_name "$gto" "$dmrg")" "$arch")"
    log="${LOG_DIR}/${tag##*:}.log"
    if [ "$SKIP_EXISTING" = "1" ] && [ "$(docker image inspect "$tag" \
            --format '{{index .Config.Labels "org.nopt.ref"}}' 2>/dev/null)" = "$NOPT_GIT_REF" ]; then
        echo ">> skip  ${tag} (already built from this commit)"
        return
    fi
    echo ">> build ${tag}"
    t0=$(date +%s)
    for try in $(seq 0 "$RETRIES"); do
        if docker buildx build -f ../Dockerfile .. --progress=plain --platform "$platform" --load \
                -t "$tag" \
                --build-arg ARCH="$arch" \
                --build-arg GTO_MAX="$gto" \
                --build-arg DMRG="$([ "$dmrg" = with ] && echo true || echo false)" \
                --build-arg NOPT_GIT_REF="$NOPT_GIT_REF" \
                ${EXTRA[@]+"${EXTRA[@]}"} > "$log" 2>&1; then
            echo "   done  ${tag} in $(( $(date +%s) - t0 ))s"
            return
        fi
        echo "!! attempt $((try + 1)) failed for ${tag} — see ${log}" >&2
        tail -5 "$log" >&2
    done
    echo "!! FAILED ${tag}" >&2
    exit 1
}

build_arch() {  # <platform> <arch>
    local g d
    for g in "${GTOS[@]}"; do
        for d in "${DMRGS[@]}"; do build_one "$1" "$2" "$g" "$d"; done
    done
}

if [ "$AMD" = "1" ]; then
    for a in "${AMD_ARCHES[@]}"; do build_arch linux/amd64 "$a"; done
fi
if [ "$ARM" = "1" ]; then
    for a in "${ARM_ARCHES[@]}"; do build_arch linux/arm64 "$a"; done
fi

echo ">> all per-arch images built. Next: PUSH=0 ./build-multiarch.sh (dispatch images),"
echo "   then PUSH=1 ./build-multiarch.sh to publish."

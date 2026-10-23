#!/usr/bin/env bash
#
# Build the dispatch images for every variant and (PUSH=1) publish everything.
# Run after build-all.sh has produced the slim per-arch images locally.
#
# Per variant <v> = gto<N>-dmrg | gto<N> (tag scheme: tags.sh):
#
#   REPO:<v>-<arch>              slim per-arch images (from build-all.sh)
#   REPO:<v>-x86-64-dispatch     heavy image: bundles x86-64-v1..v4, picks at runtime
#   REPO:<v>-arm64-dispatch      heavy image: bundles the five arm64 builds
#   REPO:<v>  manifest list:
#     linux/amd64      -> <v>-x86-64-dispatch   (plain Docker; runtime select)
#     linux/amd64/v2   -> <v>-x86-64-v2         (containerd variant match)
#     linux/amd64/v3   -> <v>-x86-64-v3
#     linux/amd64/v4   -> <v>-x86-64-v4
#     linux/arm64      -> <v>-arm64-dispatch    (only one arm64 slot; runtime select)
#   REPO:latest = same manifest list as REPO:<DEFAULT_VARIANT>
#   REPO:<arch>, REPO:x86-64-dispatch, REPO:arm64-dispatch (legacy bare tags)
#              = the DEFAULT_VARIANT per-arch / dispatch images
#
# Why arm64 has no variant entries: container runtimes match amd64 microarch
# variants (v2/v3/v4) but NOT arm64 sub-levels, so arm64 gets a single slot and
# the dispatch image does the CPU selection at container start.
#
# The arm64 dispatch build runs under QEMU on an x86 host (needs
# `tonistiigi/binfmt`); it only assembles files + apt, so it is quick.
#
# Usage:
#   PUSH=0 ./build-multiarch.sh     # build all dispatch images locally only
#   PUSH=1 ./build-multiarch.sh     # + push per-arch, dispatch, manifests, latest
#   ARM=0 ...                       # amd64 half only
#   GTOS=27 DMRGS=with ...          # a subset (see tags.sh)
set -euo pipefail

cd "$(dirname "$0")"
. ./tags.sh

PUSH=${PUSH:-0}
ARM=${ARM:-1}

AMD_VARIANTS=(v2 v3 v4)    # amd64 manifest variant entries (slim images)

# Resolve the linux/<arch> (no-variant) child digest of a tag. Handles both a
# plain image manifest and an OCI index (buildx-with-attestations).
child_digest() {  # <ref> <arch>
    local ref=$1 arch=$2 raw
    raw=$(docker buildx imagetools inspect --raw "$ref")
    if printf '%s' "$raw" | jq -e 'has("manifests")' >/dev/null 2>&1; then
        printf '%s' "$raw" | jq -r --arg a "$arch" '
            [ .manifests[]
              | select(.platform.os=="linux" and .platform.architecture==$a) ]
            | ( map(select((.platform.variant // "")=="")) + . )[0].digest'
    else
        docker buildx imagetools inspect "$ref" --format '{{.Manifest.Digest}}'
    fi
}

build_dispatch() {  # <variant>
    local v=$1
    echo ">> build ${REPO}:$(arch_tag "$v" x86-64-dispatch)"
    docker build -f Dockerfile.dispatch --provenance=false \
        --build-arg SRC="${REPO}" --build-arg TAG_PREFIX="${v}-" \
        --build-arg TAG_1="${AMD_ARCHES[0]}" --build-arg TAG_2="${AMD_ARCHES[1]}" \
        --build-arg TAG_3="${AMD_ARCHES[2]}" --build-arg TAG_4="${AMD_ARCHES[3]}" \
        -t "${REPO}:$(arch_tag "$v" x86-64-dispatch)" .

    [ "$ARM" = "1" ] || return 0
    echo ">> build ${REPO}:$(arch_tag "$v" arm64-dispatch) (QEMU)"
    docker buildx build -f Dockerfile.dispatch --provenance=false \
        --platform linux/arm64 --load \
        --build-arg SRC="${REPO}" --build-arg TAG_PREFIX="${v}-" \
        --build-arg TAG_1="${ARM_ARCHES[0]}" --build-arg TAG_2="${ARM_ARCHES[1]}" \
        --build-arg TAG_3="${ARM_ARCHES[2]}" --build-arg TAG_4="${ARM_ARCHES[3]}" \
        --build-arg TAG_5="${ARM_ARCHES[4]}" \
        -t "${REPO}:$(arch_tag "$v" arm64-dispatch)" .
}

publish_manifest() {  # <list tag> <variant>
    local list=$1 v=$2 refs k
    declare -A DIG
    DIG[amd_dispatch]=$(child_digest "${REPO}:$(arch_tag "$v" x86-64-dispatch)" amd64)
    for k in "${AMD_VARIANTS[@]}"; do DIG[$k]=$(child_digest "${REPO}:$(arch_tag "$v" "x86-64-${k}")" amd64); done
    [ "$ARM" = "1" ] && DIG[arm_dispatch]=$(child_digest "${REPO}:$(arch_tag "$v" arm64-dispatch)" arm64)
    for k in "${!DIG[@]}"; do printf '   %-13s %s\n' "$k" "${DIG[$k]}"; done

    echo ">> assemble manifest list: ${REPO}:${list}"
    docker manifest rm "${REPO}:${list}" 2>/dev/null || true

    # amd64 dispatch first so a client that just grabs the first amd64 entry
    # (rather than matching the variant) still gets the safe fallback.
    refs="${REPO}@${DIG[amd_dispatch]}"
    for k in "${AMD_VARIANTS[@]}"; do refs="${refs} ${REPO}@${DIG[$k]}"; done
    [ "$ARM" = "1" ] && refs="${refs} ${REPO}@${DIG[arm_dispatch]}"
    # shellcheck disable=SC2086
    docker manifest create "${REPO}:${list}" ${refs}

    docker manifest annotate "${REPO}:${list}" "${REPO}@${DIG[amd_dispatch]}" --os linux --arch amd64
    for k in "${AMD_VARIANTS[@]}"; do
        docker manifest annotate "${REPO}:${list}" "${REPO}@${DIG[$k]}" --os linux --arch amd64 --variant "${k}"
    done
    if [ "$ARM" = "1" ]; then
        docker manifest annotate "${REPO}:${list}" "${REPO}@${DIG[arm_dispatch]}" --os linux --arch arm64
    fi
    docker manifest push "${REPO}:${list}"
}

mapfile -t VARIANTS < <(variants)

for v in "${VARIANTS[@]}"; do build_dispatch "$v"; done

if [ "${PUSH}" != "1" ]; then
    echo ">> PUSH=0: dispatch images built locally; nothing pushed."
    exit 0
fi

for v in "${VARIANTS[@]}"; do
    echo ">> push ${v}: per-arch + dispatch"
    for a in "${AMD_ARCHES[@]}"; do docker push "${REPO}:$(arch_tag "$v" "$a")"; done
    docker push "${REPO}:$(arch_tag "$v" x86-64-dispatch)"
    if [ "$ARM" = "1" ]; then
        for a in "${ARM_ARCHES[@]}"; do docker push "${REPO}:$(arch_tag "$v" "$a")"; done
        docker push "${REPO}:$(arch_tag "$v" arm64-dispatch)"
    fi
    publish_manifest "$v" "$v"
done

# latest + the legacy bare tags (x86-64-v3, apple-m1, x86-64-dispatch, ...)
# follow the default variant, so existing pins keep receiving updates.
if printf '%s\n' "${VARIANTS[@]}" | grep -qx "$DEFAULT_VARIANT"; then
    publish_manifest latest "$DEFAULT_VARIANT"
    LEGACY=("${AMD_ARCHES[@]}" x86-64-dispatch)
    [ "$ARM" = "1" ] && LEGACY+=("${ARM_ARCHES[@]}" arm64-dispatch)
    for t in "${LEGACY[@]}"; do
        echo ">> legacy ${REPO}:${t} -> ${REPO}:$(arch_tag "$DEFAULT_VARIANT" "$t")"
        docker tag "${REPO}:$(arch_tag "$DEFAULT_VARIANT" "$t")" "${REPO}:${t}"
        docker push "${REPO}:${t}"
    done
fi

echo ">> done. verify with:"
echo "   docker buildx imagetools inspect ${REPO}:latest"

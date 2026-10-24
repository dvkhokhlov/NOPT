# Tag scheme shared by build-all.sh and build-multiarch.sh (sourced, not run).
#
#   REPO:<variant>                 manifest list (amd64 + arm64, runtime/variant select)
#   REPO:<variant>-<arch>          slim per-arch image, e.g. gto27-dmrg-x86-64-v3
#   REPO:<variant>-x86-64-dispatch heavy amd64 image (all x86-64 levels bundled)
#   REPO:<variant>-arm64-dispatch  heavy arm64 image (all arm64 levels bundled)
#   REPO:latest                    = REPO:<DEFAULT_VARIANT>
#   REPO:<arch>, REPO:*-dispatch   legacy bare tags = the DEFAULT_VARIANT images
#
# <variant> = gto<GTO_MAX>-dmrg (with block2) | gto<GTO_MAX> (without)

REPO=${REPO:-pupeza/nopt}

GTOS=(${GTOS:-20 27 30})
DMRGS=(${DMRGS:-with without})
DEFAULT_VARIANT=${DEFAULT_VARIANT:-gto27-dmrg}

AMD_ARCHES=(${AMD_ARCHES:-x86-64-v1 x86-64-v2 x86-64-v3 x86-64-v4})
ARM_ARCHES=(${ARM_ARCHES:-arm64-v8-a arm64-v8.2-a arm64-v9-a apple-m1 apple-m4})

variant_name() {  # <gto> <with|without>
    if [ "$2" = with ]; then printf 'gto%s-dmrg' "$1"; else printf 'gto%s' "$1"; fi
}

arch_tag() {  # <variant> <arch|x86-64-dispatch|arm64-dispatch>
    printf '%s-%s' "$1" "$2"
}

variants() {  # every variant name, one per line
    local g d
    for g in "${GTOS[@]}"; do
        for d in "${DMRGS[@]}"; do variant_name "$g" "$d"; echo; done
    done
}

# Publishing the pre-built images (maintainers)

Scripts that build every variant of [`../Dockerfile`](../Dockerfile) for every
CPU level and publish them as multi-arch tags. Users building their own image
don't need anything here — see [`../README.md`](../README.md).

| File | Role |
| ---- | ---- |
| `tags.sh` | tag scheme, sourced by both scripts (`REPO`, variants, arch lists) |
| `build-all.sh` | every variant × every arch → local slim images (no push) |
| `Dockerfile.dispatch` | heavy image bundling all levels of one arch, picks one at runtime |
| `build-multiarch.sh` | dispatch images, then (`PUSH=1`) push + manifest lists + `latest` + legacy tags |

```
./build-all.sh                  # every variant x every arch -> local images (no push)
PUSH=0 ./build-multiarch.sh     # dispatch images for every variant (local)
PUSH=1 ./build-multiarch.sh     # push per-arch + dispatch, manifest lists, latest, legacy tags
```

`build-all.sh` builds arch by arch so each arch's libint/block2 layers are
compiled once and reused by its six variants; all variants get the same NOPT
commit (upstream HEAD resolved once, `NOPT_GIT_REF`). `SKIP_EXISTING=1` resumes
an interrupted run, failed images are retried (`RETRIES`). `GTOS`, `DMRGS`,
`ARM=0` restrict either script to a subset. Per-image logs go to `build-logs/`.
arm64 is built under QEMU on an x86 host (`tonistiigi/binfmt --install arm64`
once); libint's code generation runs natively even then (see the Dockerfile).

## Tag scheme

| Tag | What it is |
| --- | ---------- |
| `<variant>` (e.g. `gto27-dmrg`), `latest` | multi-arch manifest list — picks the best build for the CPU |
| `<variant>-<arch>` (e.g. `gto27-dmrg-x86-64-v3`) | slim single-arch image |
| `<variant>-x86-64-dispatch`, `<variant>-arm64-dispatch` | heavy image bundling every level of one arch, picked at runtime |
| `<arch>`, `x86-64-dispatch`, `arm64-dispatch` (bare, legacy) | same images as `gto27-dmrg-<arch>` |

`<variant>` is `gto20`, `gto27`, `gto30`, each optionally with `-dmrg`;
`latest` = `gto27-dmrg`. `<arch>` is one of `x86-64-v1`…`x86-64-v4`,
`arm64-v8-a`, `arm64-v8.2-a`, `arm64-v9-a`, `apple-m1`, `apple-m4`.

## How a pull picks the build

Each `<variant>` tag is a manifest list:

| Host | What it gets |
| ---- | ------------ |
| **amd64**, containerd-based runtime (Kubernetes, recent Docker w/ containerd image store, `nerdctl`) | the **slim** image whose amd64 microarch variant (`v2`/`v3`/`v4`) matches the CPU |
| **amd64**, anything that ignores amd64 variants (older Docker, classic image store) | `<variant>-x86-64-dispatch` (bundles every amd64 build, `exec`s the right one at container start) |
| **arm64** (Graviton, Ampere, Apple Silicon under Docker Desktop, …) | `<variant>-arm64-dispatch` |

Container runtimes match the amd64 `v2/v3/v4` variants against the host CPU,
but not arm64 sub-levels — arm64 gets a single `linux/arm64` slot. So the arm64
slot is always the dispatch image, which inspects the CPU (`/proc/cpuinfo`
features, or the Apple implementer id) and picks among `arm64-v8-a`,
`arm64-v8.2-a`, `arm64-v9-a`, `apple-m1` and `apple-m4` (Apple M4 and newer,
detected by SME or the ARMv8.7 `afp`/`rpres` features), falling back to
`arm64-v8-a` if unsure. Baseline-only amd64 ("v1") CPUs are likewise served by
`x86-64-dispatch` (there is no `v1` variant — baseline x86-64 *is* plain amd64).

Inside a dispatch image `-e NOPT_ARCH=<arch>` forces a level and
`-e NOPT_DEBUG=1` logs which build was selected.

arm64 images are runtime-tested on Apple Silicon (M1, M4, M5) through the
dispatch auto-selection; other arm64 hosts (Graviton, Ampere) are built but not
yet runtime-tested.

#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."

command -v julia >/dev/null || { echo 'Julia required on build machine' >&2; exit 1; }
command -v R >/dev/null || { echo 'R required on build machine' >&2; exit 1; }

julia --startup-file=no scripts/prepare.jl
julia --startup-file=no --project=julia/PMLRBridge scripts/build.jl
R CMD build Rpkg
echo 'Built PMLCFAEFA_0.1.6.tar.gz'

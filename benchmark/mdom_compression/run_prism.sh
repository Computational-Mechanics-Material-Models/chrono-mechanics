#!/usr/bin/env bash
set -euo pipefail

usage() {
    echo "usage: $0 [worker-domains] [quick|full] [extra ldpm_uc arguments...]" >&2
    exit 2
}

worker_domains="${1:-2}"
mode="${2:-quick}"
[[ "$worker_domains" =~ ^[1-9][0-9]*$ ]] || usage
[[ "$mode" == "quick" || "$mode" == "full" ]] || usage
if (($# >= 2)); then
    shift 2
else
    shift "$#"
fi

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/../.." && pwd)"
chrono_build="${MDOM_CHRONO_BUILD:-${repo_root}/build-mdom}"
case_build="${MDOM_CASE_BUILD:-${repo_root}/build-mdom-compression}"
executable="${case_build}/ldpm_uc"
mpi_exec="${MPIEXEC:-mpiexec}"

[[ -x "$executable" ]] || {
    echo "missing executable: $executable" >&2
    echo "build the example by following ${script_dir}/README.md" >&2
    exit 3
}

if [[ "$mode" == "quick" ]]; then
    t_end="0.0001"
    hist_stride="20"
else
    t_end="0.4"
    hist_stride="200"
fi

export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export LD_LIBRARY_PATH="${chrono_build}/lib:${LD_LIBRARY_PATH:-}"

total_ranks=$((worker_domains + 1))
echo "Launching ${worker_domains} worker domain(s) plus one coordinator (${total_ranks} MPI ranks)."

exec "$mpi_exec" -n "$total_ranks" "$executable" \
    --dt=5e-7 \
    "--t_end=${t_end}" \
    --disp_rate=-5 \
    --poly_a=-1250 \
    --switch_time=0.002 \
    --load_scale=1 \
    "--hist_stride=${hist_stride}" \
    --vtk_stride=0 \
    --slice_axis=x \
    --no_diagonal_lumping \
    --partition_update \
    --no_tag_reaction_selector \
    --fix_top_rotations \
    --fix_boundary \
    --clean_output \
    "$@"

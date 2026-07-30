#!/bin/bash

set -euo pipefail

repo="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$repo"
source "$repo/setup_env.sh"

STARTING_FILES="${DIJET_UNFOLDING_PATH}/results/starting_files"
cp -r "$STARTING_FILES"/* "$DIJET_BUILD_PATH"/.

roounfold_lib="${ROOUNFOLD_PATH:-$repo/external/RooUnfold/build}"
roounfold_src="${ROOUNFOLD_ROOT:-$repo/external/RooUnfold}/src"
root_setup="gSystem->SetBuildDir(\"$repo/build\",kTRUE); gSystem->AddIncludePath(\"-I$roounfold_src\"); gSystem->AddLinkedLibs(\"-L$roounfold_lib -lRooUnfold\"); if (gSystem->Load(\"$roounfold_lib/libRooUnfold\") < 0) gSystem->Exit(1);"

# Keep each path as a separate array element. Command substitution would store
# the newline-delimited output as one scalar and pass every macro to a single
# ROOT `.L` command.
mapfile -t macros < <(
  find "$repo" \
    \( -path "$repo/build" -o -path "$repo/external" \) -prune \
    -o -type f \( -name "*.C" -o -name "*.cxx" \) -print | sort
)

for macro in "${macros[@]}"; do
  echo "ACLiC $macro"
  root -l -b -q -e "$root_setup if (!gSystem->CompileMacro(\"$macro\",\"kf\")) gSystem->Exit(1);"
done
echo "Built ${#macros[@]} macros"

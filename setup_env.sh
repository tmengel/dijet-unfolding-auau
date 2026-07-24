#!/usr/bin/env bash

DIJET_UNFOLDING_PATH="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOUNFOLD_ROOT="${DIJET_UNFOLDING_PATH}/external/RooUnfold"
ROOUNFOLD_PATH="$ROOUNFOLD_ROOT/build"
echo "DIJET_UNFOLDING_PATH=$DIJET_UNFOLDING_PATH"

path_prepend_once() {
  local var_name="$1"
  local new_path="$2"
  local old_value="${!var_name:-}"

  case ":$old_value:" in
    *":$new_path:"*) export "$var_name=$old_value" ;;
    *) export "$var_name=$new_path${old_value:+:$old_value}" ;;
  esac
}

export DIJET_UNFOLDING_PATH
export ROOUNFOLD_ROOT
export ROOUNFOLD_PATH

export DIJET_TNTUPLE_PATH="${DIJET_UNFOLDING_PATH}/rootfiles"
export DIJET_BUILD_PATH="${DIJET_UNFOLDING_PATH}/build"
export DIJET_CONFIG_PATH="${DIJET_UNFOLDING_PATH}/configs"
export DIJET_LOG_PATH="${DIJET_UNFOLDING_PATH}/logs"

mkdir -p "$DIJET_BUILD_PATH"

# Make RooUnfold available to interpreted ROOT macros as well as Makefile builds.
path_prepend_once ROOT_INCLUDE_PATH "${ROOUNFOLD_ROOT}/src"
path_prepend_once LD_LIBRARY_PATH "$ROOUNFOLD_PATH"

echo "ROOUnfold ROOT path: $ROOUNFOLD_ROOT"
echo "LD_LIBRARY_PATH=$LD_LIBRARY_PATH"

mkdir -p "${DIJET_UNFOLDING_PATH}"/{auau_plots,final_plots,jer,jer_plots,logs,njet,response_matrices,rootfiles,systematic_plots,truth_hists,uncertainties,unfolding_hists,unfolding_plots,vertex,centrality,sumeT}

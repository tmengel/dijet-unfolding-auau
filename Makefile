SHELL := /bin/bash

.PHONY: build clean preprocess

build:
	./build.sh

clean:
	rm -rf build/*

preprocess:
	./run_preprocess.sh

# full from-ntuples unfold plots closure:
# 	./scripts/run_workflow.sh $@

# help:
# 	./scripts/run_workflow.sh help

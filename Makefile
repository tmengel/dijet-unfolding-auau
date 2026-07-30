SHELL := /bin/bash

.PHONY: build clean preprocess

build:
	./build.sh

clean:
	rm -f *.o *.so *.d *.pcm *ACLiC* build/*



# full from-ntuples unfold plots closure:
# 	./scripts/run_workflow.sh $@

# help:
# 	./scripts/run_workflow.sh help

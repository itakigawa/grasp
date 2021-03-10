all: grasp utils

.PHONY: grasp utils

grasp:
	$(MAKE) -C grasp 

utils:
	$(MAKE) -C utils 

clean:
	$(MAKE) -C grasp clean
	$(MAKE) -C utils clean
	rm -f ../bin/*

install:
	cp ./grasp/grasp ../bin/
	cp ./utils/finder ../bin/
	cp ./utils/evaluator ../bin/
	cp ./utils/predictor ../bin/



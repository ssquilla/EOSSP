INSTANCES = LKH3/TSPTW/INSTANCES/eosm/
TOURS = LKH3/TSPTW/TOURS/eosm/
PAR = LKH3/TSPTW/PAR/eosm/
INIT = LKH3/TSPTW/INIT/eosm/
LKH = LKH3/
clean:
	cd $(INSTANCES);if(($(ls|wc -l)>0)); then rm *;fi
	cd $(TOURS);if(($(ls|wc -l) >0)); then  rm * ;fi
	cd $(PAR);if(($(ls|wc -l) >0)); then  rm * ;fi
	cd $(INIT);if(($(ls|wc -l) >0)); then  rm * ;fi
all:
	cd $(LKH); make all;

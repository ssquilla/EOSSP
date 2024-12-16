INSTANCES = LKH3/TSPTW/INSTANCES/eosm/
TOURS = LKH3/TSPTW/TOURS/eosm/
LKH = LKH3/
clean:
	cd $(INSTANCES);if(($(ls|wc -l)>0)); then rm *;fi
	cd $(TOURS);if(($(ls|wc -l) >0)); then  rm * ;fi
all:
	cd $(LKH); make all;

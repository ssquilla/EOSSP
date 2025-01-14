EXTERNAL_SOLVER = solvers
LKH_SOURCES = LKH3
LKH_TARGET = $(EXTERNAL_SOLVER)/LKH
INSTANCES = $(LKH_TARGET)/INSTANCES/
TOURS = $(LKH_TARGET)/TOURS/
PAR = $(LKH_TARGET)/PAR/
INIT = $(LKH_TARGET)/INIT/
LKH = LKH3/
clean:
	cd $(INSTANCES);if(($(ls|wc -l)>0)); then rm *;fi
	cd $(TOURS);if(($(ls|wc -l) >0)); then  rm * ;fi
	cd $(PAR);if(($(ls|wc -l) >0)); then  rm * ;fi
	cd $(INIT);if(($(ls|wc -l) >0)); then  rm * ;fi
download_lkh:
	rm -rf LKH-3.0.13.tgz
	rm -r $(LKH_SOURCES)
	wget http://akira.ruc.dk/~keld/research/LKH-3/LKH-3.0.13.tgz
	tar -zxvf  LKH-3.0.13.tgz
	mv LKH-3.0.13 $(LKH_SOURCES)
	rm -r LKH-3.0.13.tgz
all:
	cd $(LKH_SOURCES); make all;
	cp $(LKH_SOURCES)/LKH $(LKH_TARGET)/LKH; 

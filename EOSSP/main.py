import sys, importlib
from pathlib import Path
from resource import *

def importParents(level=1):
    global __package__
    file = Path(__file__).resolve()
    parent, top = file.parent, file.parents[level]
        
    sys.path.append(str(top))
    try:
        sys.path.remove(str(parent))
    except ValueError: # already removed
        pass
    
    __package__ = '.'.join(parent.parts[len(top.parts):])
    importlib.import_module(__package__) # won't be needed after that
        

if __name__ == '__main__' and __package__ is None:
    importParents(level=1)

"""
Created on Tue Jan 10 16:29:22 2023

@author: ssquilla, Squillaci Samuel.
"""

from EOSSP.Utils.config import *
from EOSSP.Utils.Utils import *
global config
config = Config()
instance = config.instance
if config.getOptValue("help"):
    if MPI.COMM_WORLD.Get_rank()==0:
        config.displayHelp()
else:  
    from EOSSP.LNS.LNSSolver import runnableLNS
    from EOSSP.CPModel.CPSolver import runnableCPSolver
    from EOSSP.timeDependentSolver.timeDependentSolver import runnableTimeDependentSolver
    from EOSSP.BPCCAS.BPCCAS import runnableBPCCAS
    from EOSSP.UPCCAS.UPCCAS import runnableUPCCAS
    
    from EOSSP.model.transitionModel import *
    from EOSSP.model.components import *
    from EOSSP.model.componentPlan import *
        
    from EOSSP.model.constellation import Constellation
    
    from mpi4py import MPI
    
    if config.getOptValue("help"):
        if MPI.COMM_WORLD.Get_rank()==0:
            config.displayHelp()
    else:
     
        def enumerateAlgorithms():
            algorithms = []
            algorithms.append(("LNS",runnableLNS))
            algorithms.append(("timeDependentSolver",runnableTimeDependentSolver))
            algorithms.append(("CPSolver",runnableCPSolver))
            algorithms.append(("BPCCAS",runnableBPCCAS))
            algorithms.append(("UPCCAS",runnableUPCCAS))
            return algorithms
        
        def choseFileAndSetupProblem():
            folder = config.getOptValue("data")
            path = '../data/'+folder
            comm = MPI.COMM_WORLD
            rank = comm.Get_rank()
            components = None
            
            if rank==0:
                printColor("\n")
                printColor("Resolution parameters:",c='b')
                printColor("| instance: "+str(instance),c='b')
                printColor("| time limit: "+str(config.getOptValue("time")) + " (s)",c='b')
                printColor("\n")
                file = choseAndBroadcastFile(path,instance)
                comm.Barrier()
                startDate = time()
                constellation = Constellation(file,startDate,True)
            else:
                data = None
                data = comm.bcast(data,root=0)
                file = data['file']
                comm.Barrier()
                startDate = time()
                constellation = Constellation(file,startDate,False)
            components,transitionModel = loadComponentsIfPrecomputed(constellation,file)
            id_mpi = MPI.COMM_WORLD.Get_rank()
            seed_option = config.getOptValue("seed")
            rd.seed(id_mpi+seed_option)
            comm.Barrier()
            return startDate,constellation,components,transitionModel
        
        def loadComponentsIfPrecomputed(constellation,file):
            if "components_time_dep" in os.listdir(Path(file).parent):
                printMaster("Reading time-dependent model...",end ="",c='g')
                start = time()
                folder = Path(file).parent.__str__()+"/components_time_dep"
                components = ComposantesStatiquesPrecalculees(constellation,folder)
                if config.isTimeDependentModeOn():
                    precompute_approx = None
                    if config.getOptValue("initTransitionModel") != "time-dep":
                        precompute_approx = config.getOptValue("initTransitionModel")
                    timeDepModel = modelTransitionTimeDependent(folder,constellation,components,precompute_approx,profile=config.getOptValue("profile"))
                else:
                    timeDepModel = None
                end = time()
                printMaster("OK. "+str(round(end-start))+" s.",c='g')
                return components,timeDepModel
            else:
                return None,None
        
        
        def choseTransitionModel(timeDepModel):
            choix = config.getOptValue("initTransitionModel")
            if timeDepModel is not None:
                start = time()
                printMaster("Creating approximated model...",end="",c='g')
                if choix=="fast":
                    transitionModel = modeleOptimisteExtrait(timeDepModel)
                elif choix=="mean":
                    transitionModel = modeleMoyenExtrait(timeDepModel)
                elif choix=="time-dep":
                    transitionModel = timeDepModel
                else:
                    assert(choix=="slow")
                    transitionModel = modelePessimisteExtrait(timeDepModel)
                end = time()
                printMaster("OK. "+str(round(end-start))+ " s.",c='g')
                return transitionModel
            else:
                availableModels = []
                availableModels.append(("fast",FastModel))
                availableModels.append(("slow",SlowModel))
                availableModels.append(("mean",MeanModel))
                for opt,transitionModel in availableModels:
                    if opt == choix:
                        return transitionModel
            raise ValueError("Unknown model:",choix)
                
        def runChosenSolver():
            availableAlgorithms = enumerateAlgorithms()
            solver_opt = config.getOptValue("solver")
            find = False
            for opt,solver_runnable in availableAlgorithms:
                if solver_opt == opt:
                    solver = solver_runnable
                    find = True
                    break
            if not(find):
                raise NameError('Unknown solver: '+solver_opt)
            return solver()
                
        startDate,constellation,components,timeDepModel = choseFileAndSetupProblem()
        deltaTimeTransition = time() - startDate
        
        if not config.getOptValue("profile"):
            startDate = time() # demarrage de la clock apres la construction du modele de transition
            
            transitionModel = choseTransitionModel(timeDepModel)
            printMaster("Chosen transition model:",transitionModel)
            
            if config.isTimeDependentModeOn():
                printMaster("- model size: ",'{:,d}'.format(timeDepModel.taille_time_dep_vect) +" coefficients.")
                proportionFIFO = timeDepModel.proportionCouplesFIFO
                if proportionFIFO != 1:
                    printMaster("Existing non-FIFO elements. Proportion of non-FIFO: ",proportionFIFO,"%",c='y')
            runnableSolver = runChosenSolver()
            printMaster("Chosen solver:",runnableSolver.getName())
            if runnableSolver.getName()=="Time Dependent Solver":
                # ce solver a besoin de deux modèle de transition : approche et reel (time-dep)
                solutionContainer = runnableSolver.execute(constellation,startDate,transitionModel,timeDepModel,deltaTimeTransition,CCAs=components)
            else:
                solutionContainer = runnableSolver.execute(constellation,startDate,transitionModel,deltaTimeTransition,CCAs=components)
            
            masterSolver = runnableSolver.getMasterSolver()
            
            if masterSolver is not None: # not None seulement dans le cas du Master
                if config.getOptValue("sample") is not None:
                    masterSolver.saveSample(constellation)  
                
                if config.getOptValue("cpu"):
                    masterSolver.plotCPU()
                if config.getOptValue("obj"):
                    masterSolver.plotObjective()
                if config.getOptValue("ccaload"):
                    masterSolver.plotCCAsLoad()
            

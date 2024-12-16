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

@author: ssquilla
"""
from .Utils.config import *
from .Utils.Utils import *
global config
config = Config()
instance = config.instance
if config.getOptValue("help"):
    if MPI.COMM_WORLD.Get_rank()==0:
        config.displayHelp()
else:  
    from .LNS.LNSSolver import runnableLNS
    from .CPModel.CPSolver import runnableCPSolver
    from .timeDependentSolver.timeDependentSolver import runnableTimeDependentSolver
    from .BPCCAS.BPCCAS import runnableBPCCAS
    from .UPCCAS.UPCCAS import runnableUPCCAS
    
    from .model.transitionModel import *
    from .model.components import *
    from .model.componentPlan import *
        
    from .model.constellation import Constellation
    
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
                printColor("Paramètres de résolution :",c='b')
                printColor("| instance : "+str(instance),c='b')
                printColor("| temps de calcul max : "+str(config.getOptValue("time")/60),c='b')
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
                printMaster("Lecture du modèle time-dependent ...",end ="",c='g')
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
                printMaster("Création du modèle approché ...",end="",c='g')
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
            raise ValueError("Choix de modèle inconnu :",choix)
        
        def printRandomTransition(constellation,components,timeDepModel,approximatedModel):
            cca = rd.choice(components.getComposantes())
            acts = rd.choices(components.getActivitesComposante(cca),k=2)
            for a1 in acts:
                for a2 in acts:
                    times = np.linspace(constellation.getActivite(a1).getDebut(),constellation.getActivite(a1).getFin(),10)
                    timeDepTransition = [timeDepModel.getTransition(a1,a2,t) for t in times]
                    approximatedTransition = approximatedModel.getTransition(None,a1,a2,4)
                    print("======================================")
                    print(a1,a2,"Time-dep:",timeDepTransition,"approche",approximatedTransition)
                    print(timeDepModel.getPointsDeControle(a1,a2),times)
                    print("======================================")
            
                
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
                raise NameError('Solver inconnu : '+solver_opt)
            return solver()
                
        startDate,constellation,components,timeDepModel = choseFileAndSetupProblem()
        deltaTimeTransition = time() - startDate
        
        if not config.getOptValue("profile"):
            startDate = time() # demarrage de la clock apres la construction du modele de transition
            
            transitionModel = choseTransitionModel(timeDepModel)
            printMaster("Choix du modèle de transition :",transitionModel)
            
            if config.isTimeDependentModeOn():
                printMaster("- taille du modèle : ",'{:,d}'.format(timeDepModel.taille_time_dep_vect) +" coefficients.")
                proportionFIFO = timeDepModel.proportionCouplesFIFO
                if proportionFIFO != 1:
                    printMaster("Présence de couples non FIFO. Proportion de couples FIFO : ",proportionFIFO,"%",c='y')
            runnableSolver = runChosenSolver()
            printMaster("Solver choisi :",runnableSolver.getName())
            if runnableSolver.getName()=="Time Dependent Solver":
                # ce solver a besoin de deux modèle de transition : approche et reel (time-dep)
                solutionContainer = runnableSolver.execute(constellation,startDate,transitionModel,timeDepModel,deltaTimeTransition,CCAs=components)
            else:
                solutionContainer = runnableSolver.execute(constellation,startDate,transitionModel,deltaTimeTransition,CCAs=components)
            
            masterSolver = runnableSolver.getMasterSolver()
                
            if masterSolver is not None: # not None seulement dans le cas du Master
                if config.getOptValue("cpu"):
                    masterSolver.tracerCPU()
                if config.getOptValue("obj"):
                    masterSolver.tracerObjectif()
                if config.getOptValue("sample") is not None:
                    masterSolver.saveSample(constellation)    
                if config.getOptValue("charge"):
                    masterSolver.tracerChargeCCAs()

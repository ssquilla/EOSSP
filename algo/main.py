import sys, importlib
from pathlib import Path
from resource import *

def import_parents(level=1):
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
    import_parents(level=1)

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
        config.afficherAide()
else:  
    from .LNS.LNSSolver import runnableLNS
    from .CPModel.CPSolver import runnableCPSolver
    from .timeDependentSolver.timeDependentSolver import runnableTimeDependentSolver
    from .BPCCAS.BPCCAS import runnableBPCCAS
    from .UPCCAS.UPCCAS import runnableUPCCAS
    
    from .model.modelTransition import *
    from .model.composantes import *
    from .model.solution_composantes import *
    
    
    
    from .model.constellation import Constellation
    
    from mpi4py import MPI
    
    if config.getOptValue("help"):
        if MPI.COMM_WORLD.Get_rank()==0:
            config.afficherAide()
    else:
     
        def listerPossibilites():
            possibilites = []
            possibilites.append(("LNS",runnableLNS))
            possibilites.append(("timeDependentSolver",runnableTimeDependentSolver))
            possibilites.append(("CPSolver",runnableCPSolver))
            possibilites.append(("BPCCAS",runnableBPCCAS))
            possibilites.append(("UPCCAS",runnableUPCCAS))
            return possibilites
        
        def choseFileAndSetupProblem():
            folder = config.getOptValue("data")
            path = '../data/'+folder
            comm = MPI.COMM_WORLD
            rank = comm.Get_rank()
            composantes = None
            
            if rank==0:
                printColor("\n")
                printColor("Paramètres de résolution :",c='b')
                printColor("| instance : "+str(instance),c='b')
                printColor("| temps de calcul max : "+str(config.getOptValue("time")/60),c='b')
                printColor("\n")
                file = choseAndBroadcastFile(path,instance)
                comm.Barrier()
                start_date = time()
                constellation = Constellation(file,start_date,True)
            else:
                data = None
                data = comm.bcast(data,root=0)
                file = data['file']
                comm.Barrier()
                start_date = time()
                constellation = Constellation(file,start_date,False)
            composantes,modeleDeTransition = chargerCCASiPrecalculees(constellation,file)
            id_mpi = MPI.COMM_WORLD.Get_rank()
            seed_option = config.getOptValue("seed")
            rd.seed(id_mpi+seed_option)
            comm.Barrier()
            return start_date,constellation,composantes,modeleDeTransition
        
        def chargerCCASiPrecalculees(constellation,file):
            if "components_time_dep" in os.listdir(Path(file).parent):
                printMaster("Lecture du modèle time-dependent ...",end ="",c='g')
                start = time()
                folder = Path(file).parent.__str__()+"/components_time_dep"
                composantes = ComposantesStatiquesPrecalculees(constellation,folder)
                if config.isTimeDependentModeOn():
                    precompute_approx = None
                    if config.getOptValue("initTransitionModel") != "time-dep":
                        precompute_approx = config.getOptValue("initTransitionModel")
                    modeleTimeDep = modelTransitionTimeDependent(folder,constellation,composantes,precompute_approx,profile=config.getOptValue("profile"))
                else:
                    modeleTimeDep = None
                end = time()
                printMaster("OK. "+str(round(end-start))+" s.",c='g')
                return composantes,modeleTimeDep
            else:
                return None,None
        
        
        
        def choisirModeleTransition(modeleTimeDep):
            choix = config.getOptValue("initTransitionModel")
            if modeleTimeDep is not None:
                start = time()
                printMaster("Création du modèle approché ...",end="",c='g')
                if choix=="fast":
                    modele = modeleOptimisteExtrait(modeleTimeDep)
                elif choix=="mean":
                    modele = modeleMoyenExtrait(modeleTimeDep)
                elif choix=="time-dep":
                    modele = modeleTimeDep
                else:
                    assert(choix=="slow")
                    modele = modelePessimisteExtrait(modeleTimeDep)
                end = time()
                printMaster("OK. "+str(round(end-start))+ " s.",c='g')
                return modele
            else:
                possibilites = []
                possibilites.append(("fast",modeleRapide))
                possibilites.append(("slow",modeleLent))
                possibilites.append(("mean",modeleMoyen))
                for opt,modele in possibilites:
                    if opt == choix:
                        return modele
            raise ValueError("Choix de modèle inconnu :",choix)
        
        def printRandomTransition(constellation,composantes,modeleTimeDep,modeleApproche):
            cca = rd.choice(composantes.getComposantes())
            acts = rd.choices(composantes.getActivitesComposante(cca),k=2)
            for a1 in acts:
                for a2 in acts:
                    times = np.linspace(constellation.getActivite(a1).getDebut(),constellation.getActivite(a1).getFin(),10)
                    trans_time_dep = [modeleTimeDep.getTransition(a1,a2,t) for t in times]
                    trans_approchee = modeleApproche.getTransition(None,a1,a2,4)
                    print("======================================")
                    print(a1,a2,"Time-dep:",trans_time_dep,"approche",trans_approchee)
                    print(modeleTimeDep.getPointsDeControle(a1,a2),times)
                    print("======================================")
            
                
        def lancerSolverChoisi():
            possibilites = listerPossibilites()
            solver_opt = config.getOptValue("solver")
            find = False
            for opt,solver_runnable in possibilites:
                if solver_opt == opt:
                    solver = solver_runnable
                    find = True
                    break
            if not(find):
                raise NameError('Solver inconnu : '+solver_opt)
            return solver()
                
        start_date,constellation,composantes,modeleTimeDep = choseFileAndSetupProblem()
        dt_construction_transition = time() - start_date
        
        if not config.getOptValue("profile"):
            start_date = time() # demarrage de la clock apres la construction du modele de transition
            
            modeleDeTransition = choisirModeleTransition(modeleTimeDep)
            printMaster("Choix du modèle de transition :",modeleDeTransition)
            
            if config.isTimeDependentModeOn():
                printMaster("- taille du modèle : ",'{:,d}'.format(modeleTimeDep.taille_time_dep_vect) +" coefficients.")
                proportionFIFO = modeleTimeDep.proportionCouplesFIFO
                if proportionFIFO != 1:
                    printMaster("Présence de couples non FIFO. Proportion de couples FIFO : ",proportionFIFO,"%",c='y')
            runnableSolver = lancerSolverChoisi()
            printMaster("Solver choisi :",runnableSolver.getName())
            if runnableSolver.getName()=="Time Dependent Solver":
                # ce solver a besoin de deux modèle de transition : approche et reel (time-dep)
                solution_container = runnableSolver.execute(constellation,start_date,modeleDeTransition,modeleTimeDep,dt_construction_transition,CCAs=composantes)
            else:
                solution_container = runnableSolver.execute(constellation,start_date,modeleDeTransition,dt_construction_transition,CCAs=composantes)
            
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

from solution import *
from constellation import *
from composantes import *
from solution_composantes import *
from mpi4py import MPI
from Utils import *
import math
from config import *
from dispatch import *
from time import time
from time import sleep
import docplex.cp as cp
from docplex.cp.model import *

class CpSolver(Solution):
    def __init__(self,constellation):
        super().__init__(constellation)
        activites = constellation.extraireActivitesRequetes()
        self.grapheDependances = GroupeComposantesActivitesStatiques(constellation,activites)
        printColor(self.grapheDependances,c='b')
        self.objectif = (0,0)
        self.initRequetes(constellation)
        printColor("Durée d'initialisation :",time.time()-self.start_date,c='b')

    def resetModeSolver(self,constellation,failed):
        self.model = CpoModel()
        modes = []
        for r in constellation.getRequetes():
            m = constellation.getRequete(r).getModeCourant(constellation)
            if r in failed:
                m = constellation.getRequete(r).getModeSuivantSansExp(constellation)
            if m is not None:
                modes.append((r,m.getId()))
        return modes
    
    def initVars(self,constellation):
        self.interval_vars = {}
        self.mode_var = {}
        for (r,m) in self.modes:
            self.interval_vars[(r,m)] = {}
            if r!=-1:
                self.mode_var[(r,m)] = self.model.binary_var(name="y_"+str((r,m)))
                self.model.add(self.mode_var[(r,m)])
            for (s,a) in constellation.getRequete(r).getMode(m).getCouples():
                start = int(math.floor(constellation.getSatellite(s).getActivite(a).getDebut()))
                end = int(math.ceil(constellation.getSatellite(s).getActivite(a).getFin()))
                duree = int(math.ceil(constellation.getSatellite(s).getActivite(a).getDuree()))
                assert(start+duree<=end)
                self.interval_vars[(r,m)][a] = self.model.interval_var(start=(start,end),end=(start,end),length=duree,name="Ia_{}_{}".format(s,a))
                self.model.add(self.interval_vars[(r,m)][a])
        obj1 = self.model.sum([self.recompenseBruitee(constellation,r,m)[0]*self.mode_var[(r,m)] for (r,m) in self.mode_var])
        obj2 = self.model.sum([self.model.presence_of(self.interval_vars[x][a]) for x in self.interval_vars for a in self.interval_vars[x]])
        #self.model.maximize_static_lex([obj1,obj2])
        self.model.maximize(obj1)
        
    def initConstraints(self,constellation):
        for (r,m) in self.interval_vars:
            if r!=-1:
                for o in self.interval_vars[(r,m)]:
                    self.interval_vars[(r,m)][o].set_optional()
                    self.model.add(self.model.presence_of(self.interval_vars[(r,m)][o]) == self.mode_var[(r,m)])
                    #print(self.interval_vars[(r,m)][o].is_optional())
            else:
                for d in self.interval_vars[(-1,0)]:
                    self.model.add(self.model.presence_of(self.interval_vars[(-1,0)][d]) == 1)
                    #print(self.interval_vars[(r,m)][d].is_optional())
        self.decompositionActivitesCCA()
        for cca in self.cca_act:
            s = self.grapheDependances.getSatelliteCCA(cca)
            n = len(self.cca_act[cca])
            distance_matrix = np.zeros((n,n),dtype=int)
            vars_list = []
            for i,(x1,a1) in enumerate(self.cca_act[cca]):
                vars_list.append(self.interval_vars[x1][a1])
                for j,(x2,a2) in enumerate(self.cca_act[cca]):
                    if i==j:
                        distance_matrix[i][j] = 0
                    else:
                        transition = constellation.getSatellite(s).getTransition(a1,a2)
                        distance_matrix[i][j] = int(math.ceil(transition))
            mat = self.model.transition_matrix(distance_matrix)
            seq = self.model.sequence_var(vars_list)
            self.model.add(self.model.no_overlap(seq,distance_matrix=mat))
    
    def decompositionActivitesCCA(self):
        mapping_cca_activites = {}
        for x in self.interval_vars:
            for a in self.interval_vars[x]:
                cca = self.grapheDependances.getActiviteCCA(a)
                if cca not in mapping_cca_activites:
                    mapping_cca_activites[cca] = []
                mapping_cca_activites[cca].append((x,a))
        self.cca_act = mapping_cca_activites
    
    def initIteration(self,constellation,failed):
        self.solCCAs = {s : {} for s in constellation.satellites} # liste de solutions CCA
        for cca in self.grapheDependances.getComposantes():
            s = self.grapheDependances.getSatelliteCCA(cca)
            self.solCCAs[s][cca] = SolCCA(cca,s)
        self.modes = self.resetModeSolver(constellation,failed)
        self.initVars(constellation)
        self.initConstraints(constellation)

    """
        =============================================== 
                        RESOLUTION
        =============================================== 
    """

    
    def insererSequences(self,constellation,sequences):
        requetes_retenues = [x[0] for x in self.modes_retenus]
        for cca in sequences:
            s = self.grapheDependances.getSatelliteCCA(cca)
            seq = [a for a in sequences[cca] if constellation.getRequeteActivite(a) in requetes_retenues]
            #print("set sequence",s,cca)
            self.solCCAs[s][cca].setSequence(constellation,seq)
    
    def remplirSolCCAs(self,sol):
        obs_cca = {}
        for (r,m) in self.interval_vars:
            if r==-1 or sol.get_value(self.mode_var[(r,m)])==1:
                for a in self.interval_vars[(r,m)]:
                    cca = self.grapheDependances.getActiviteCCA(a)
                    if cca not in obs_cca:
                        obs_cca[cca] = []
                    val = sol.get_value(self.interval_vars[(r,m)][a])
                    obs_cca[cca].append((val[0],a))
        for cca in obs_cca:
            seq = [x[1] for x in sorted(obs_cca[cca])] 
            s = self.grapheDependances.getSatelliteCCA(cca)
            self.solCCAs[s][cca].setSequence(constellation,seq)
        
    def resoudre(self,constellation,mailbox,afficher=True):
        temps = time.time()
        failed = []
        it = 0
        while(time.time()-self.start_date<config.glob.limite_temps and (it==0 or len(failed)>0)):
            self.initIteration(constellation,failed)
            tlim = min(int(math.floor(config.glob.limite_temps-(time.time()-self.start_date) )),config.cpSolver.iteration_limit)
            if tlim<0:
                break
            if config.cpSolver.threads is not None:
                sol = self.model.solve(TimeLimit=tlim,execfile=config.glob.docplexPath,Workers=config.cpSolver.threads,log_output=None)
            else:
                sol = self.model.solve(TimeLimit=tlim,execfile=config.glob.docplexPath,log_output=None)
            if str(sol.get_solve_status())!="Unknown":
                self.modes_retenus = [(-1,0)]+[(r,m) for (r,m) in self.mode_var if sol.get_value(self.mode_var[(r,m)])==1]
                failed = [r for (r,m) in self.mode_var if sol.get_value(self.mode_var[(r,m)])==0]
                self.objectif = self.calculerObjectif(constellation)
                self.remplirSolCCAs(sol)
                #self.verifierSolution(constellation)
                self.historique.MAJHistorique(time.time()-self.start_date,1,self.objectif,self.solCCAs,self.modes_retenus,self.grapheDependances,constellation)
            self.afficherInfo(time.time(),self.start_date,afficher,color='c',title='ITERATION '+str(it))
            it += 1
        self.afficherInfo(time.time(),self.start_date,afficher,color='y',title='FIN')
        self.objectif,self.solCCAs,self.modes_retenus = self.getBest()
        self.MAJHistorique(min(time.time()-self.start_date,config.glob.limite_temps),2,self.objectif,self.solCCAs,self.modes_retenus,self.grapheDependances,constellation)
        
class Master:
    def __init__(self):
        pass
        
    def resoudre(self,constellation,afficher=True):
        # iterer insertion-réparation
        self.initSolution(constellation)
        self.solution.resoudre(constellation,afficher)
        self.solution.construirePlan(constellation)

    def meilleureSolution(self):
        return self.solution.meilleureSolution()
    
    def releverObjectif(self):
        rel = config.glob.releves
        sol = self.meilleureSolution()
        points = []
        i_rel = 0
        for i,(t,obj,modes) in enumerate(sol):
            if i==len(sol)-1:
                points.append((t//60,t%60,obj,modes))
            elif t>rel[i_rel]:
                points.append((t//60,t%60,obj,modes))
                i_rel += 1
        return points               
    
    def tracerActivite(self,constellation,annoter=False):
        return self.solution.tracerActivite(constellation,annoter)
        
    def tracerHistorique(self,init=True):
        return self.solution.tracerHistorique(init)
    
    def saveSample(self,constellation):
        self.solution.saveSample(constellation)
        
    def initSolution(self,constellation):
        self.solution = CpSolver(constellation)
        
    def getSolution(self):
        return self.solution.getSolution()
    
    def setSolution(self,sol,vid,modes_retenus,modes_candidats,modes_retires,objectif):
        self.solution.setSol(sol,vid,modes_retenus,modes_candidats,modes_retires,objectif)
    
    def verifierSolution(self):
        self.solution.verifierSolution()
        
    def getModesSolution(self):
        return self.solution.getModes()


path = '../data'
global config
config = Config()
instance = config.instance


printColor("\n")
printColor("Paramètres de résolution :",c='b')
printColor("| instance : "+str(instance),c='b')
printColor("| temps de calcul max : "+str(config.glob.limite_temps/60),c='b')
printColor("\n")

files = choseAndBroadcastFile(path,instance)
for file in files:
    constellation = Constellation(file,True)
    process = Master()
    process.resoudre(constellation)
    process.solution.verifierSolution(constellation)
    if config.glob.saveSample:
        process.saveSample(constellation)



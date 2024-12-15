from config import *
global config
config = Config()
instance = config.instance
from mpi4py import MPI
if config.getOptValue("help"):
    if MPI.COMM_WORLD.Get_rank()==0:
        config.afficherAide()
else: 
    from solution import *
    from constellation import *
    from composantes import *
    from solution_composantes import *
    from Utils import *
    import math
    from time import time
    from time import sleep
    import docplex.cp as cp
    from docplex.cp.model import CpoModel,interval_var,binary_var
    
    
    class GlobalSolver(Solver):
        def __init__(self,constellation,start_date,grapheDependances=None):
            super().__init__(constellation,start_date)
            activites = constellation.extraireActivitesRequetes()
            if grapheDependances is None:
                self.grapheDependances = GroupeComposantesActivitesStatiques(constellation,activites)
            else:
                self.grapheDependances = grapheDependances
            printColor(self.grapheDependances,c='b')
            self.objectif = (0,0)
            self.initRequetes(constellation,initFirstMode=False)
            printColor("Durée d'initialisation :",time()-self.start_date,c='b')
    
        def resetModeSolver(self,constellation):
            self.model = CpoModel()
            modes = []
            for r in constellation.getRequetes():
                for m in constellation.getRequete(r).getKMeilleursModes(constellation,config.getOptValue("modes")):
                    modes.append((r,m.getId()))
            return modes
        
        def initVars(self,constellation):
            self.interval_vars = {}
            self.mode_var = {}
            self.requetes_activites = {}
            for (r,m) in self.modes:
                if r!=-1:
                    self.mode_var[(r,m)] = self.model.binary_var(name="y_"+str((r,m)))
                    self.model.add(self.mode_var[(r,m)])
                for (s,a) in constellation.getRequete(r).getMode(m).getCouples():
                    if a not in self.requetes_activites:
                        self.requetes_activites[a] = []
                    self.requetes_activites[a].append((r,m))
                    if a not in self.interval_vars:
                        start = int(math.floor(constellation.getSatellite(s).getActivite(a).getDebut()))
                        end = int(math.ceil(constellation.getSatellite(s).getActivite(a).getFin()))
                        duree = int(math.ceil(constellation.getSatellite(s).getActivite(a).getDuree()))
                        assert(start+duree<=end)
                        self.interval_vars[a] = self.model.interval_var(start=(start,end),end=(start,end),length=duree,name="Ia_{}_{}".format(s,a))
                        self.model.add(self.interval_vars[a])
            obj1 = sum([self.recompenseBruitee(constellation,r,m)[0]*self.mode_var[(r,m)] for (r,m) in self.mode_var])
            self.model.maximize(obj1)
            #self.model.maximize(obj2)
            
        def initConstraints(self,constellation):
            # 1 mode par requete
            for r in constellation.getRequetes():
                modes = [x for x in self.mode_var if x[0]==r]
                if len(modes)>0:
                    self.model.add(sum([self.mode_var[x] for x in modes])<=1)
            # presence des activites = presence d'un mode
            for a in self.interval_vars:
                #print(self.requetes_activites[a],a)
                r = self.requetes_activites[a][0][0]
                if r != -1:    
                    self.interval_vars[a].set_optional()
                    self.model.add(self.model.presence_of(self.interval_vars[a]) == sum([self.mode_var[x] for x in self.requetes_activites[a]]))
                else:
                    self.model.add(self.model.presence_of(self.interval_vars[a]) == 1)
            self.decompositionActivitesCCA()
            # flot sur les cca
            for id_cca in self.cca_act:
                s,cca = id_cca
                n = len(self.cca_act[id_cca])
                distance_matrix = np.zeros((n,n),dtype=int)
                vars_list = []
                for i,a1 in enumerate(self.cca_act[id_cca]):
                    vars_list.append(self.interval_vars[a1])
                    for j,a2 in enumerate(self.cca_act[id_cca]):
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
            for a in self.interval_vars:
                id_cca = self.grapheDependances.getActiviteCCA(a)
                if id_cca not in mapping_cca_activites:
                    mapping_cca_activites[id_cca] = []
                mapping_cca_activites[id_cca].append(a)
            self.cca_act = mapping_cca_activites
        
        def initIteration(self,constellation):
            self.setAllSolCCAs ({s : {} for s in constellation.satellites}) # liste de solutions CCA
            for s,cca in self.grapheDependances.getComposantes():
                self.setSolCCA(s,cca, SolCCA(cca,s))
            self.modes = self.resetModeSolver(constellation)
            self.initVars(constellation)
            self.initConstraints(constellation)
    
        """
            =============================================== 
                            RESOLUTION
            =============================================== 
        """
        
        def insererSequences(self,constellation,sequences):
            requetes_retenues = [x[0] for x in self.getModesRetenus()]
            for s,cca in sequences:
                seq = [a for a in sequences[(s,cca)] if constellation.getRequeteActivite(a) in requetes_retenues]
                #print("set sequence",s,cca)
                self.getSolCCA(s,cca).setSequence(constellation,seq)
        
        def remplirSolCCAs(self,sol):
            obs_cca = {}
            for a in self.interval_vars:
                r = self.requetes_activites[a][0][0]
                if r==-1 or sum([sol.get_value(self.mode_var[x]) for x in self.requetes_activites[a]])==1:
                    id_cca = self.grapheDependances.getActiviteCCA(a)
                    if id_cca not in obs_cca:
                        obs_cca[id_cca] = []
                    val = sol.get_value(self.interval_vars[a])
                    obs_cca[id_cca].append((val[0],a))
            for id_cca in obs_cca:
                s,cca = id_cca
                seq = [x[1] for x in sorted(obs_cca[id_cca])] 
                self.getSolCCA(s,cca).setSequence(constellation,seq)
            
        def resoudre(self,constellation,mailbox,afficher=True):
            temps = time()
            self.initIteration(constellation)
            if config.getOptValue("time")<np.Inf:
                tlim = min(int(math.floor(config.getOptValue("time")-(time()-self.start_date) )),config.getOptValue("time"))
            else:
                tlim = np.Inf
            if tlim>0:
                if config.getOptValue("threads")is not None:
                    sol = self.model.solve(TimeLimit=tlim,execfile=config.glob.docplexPath,Workers=config.getOptValue("threads"),log_output=None)
                else:
                    if config.getOptValue("time")<np.Inf:
                        sol = self.model.solve(TimeLimit=tlim,execfile=config.glob.docplexPath,log_output=None)
                    else:
                        sol = self.model.solve(execfile=config.glob.docplexPath,log_output=None)
                self.solution.historique.gap = sol.solution.get_objective_gap()
                if str(sol.get_solve_status())!="Unknown":
                    self.setModesRetenus([(-1,0)]+[(r,m) for (r,m) in self.mode_var if sol.get_value(self.mode_var[(r,m)])==1])
                    self.objectif = self.calculerObjectif(constellation)
                    self.remplirSolCCAs(sol)
                    #self.verifierSolution(constellation)
            else:
                raise Exception("pas le temps de lancer le solver")
            self.solution.historique.notifierFinExecution(min(time()-self.start_date,config.getOptValue("time")),self.objectif,self.getSolCCAs(),self.getModesRetenus(),self.grapheDependances,constellation)
            self.afficherInfo(time(),self.start_date,constellation,add={'gap':self.solution.historique.gap},color='y',title='FIN')
                
    class Master:
        def __init__(self,start_date):
            self.start_date = start_date
            
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
            self.solution = GlobalSolver(constellation,self.start_date)
            
        def getSolution(self):
            return self.solution.getSolution()
        
        def setSolution(self,sol,vid,modes_retenus,modes_candidats,modes_retires,objectif):
            self.solution.setSol(sol,vid,modes_retenus,modes_candidats,modes_retires,objectif)
        
        def verifierSolution(self):
            self.solution.verifierSolution()
            
        def getModesSolution(self):
            return self.solution.getModes()
    
    
    path = '../data'
    
    printColor("\n")
    printColor("Paramètres de résolution :",c='b')
    printColor("| instance : "+str(instance),c='b')
    printColor("| temps de calcul max : "+str(config.getOptValue("time")/60),c='b')
    printColor("\n")
    
    rd.seed(0)
    
    files = choseAndBroadcastFile(path,instance)
    for file in files:
        start_date = time()
        constellation = Constellation(file,start_date,True)
        process = Master(start_date)
        process.resoudre(constellation)
        process.solution.verifierSolution(constellation)
        if config.getOptValue("sample") is not None:
            process.saveSample(constellation)
    
    

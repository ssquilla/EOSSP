from ..Utils.Utils import *
from ..Utils.config import *

global config
config = Config()
instance = config.instance
from mpi4py import MPI
if config.getOptValue("help"):
    if MPI.COMM_WORLD.Get_rank()==0:
        config.afficherAide()
else: 
    from ..model.solution import *
    from ..model.constellation import *
    from ..model.components import *
    from ..model.componentPlan import *
    
    import math
    from time import time
    from time import sleep
    import docplex.cp as cp
    from docplex.cp.model import CpoModel,interval_var,binary_var
    
    
    class CPSolver(Solver):
        def __init__(self,constellation,start_date,modeleDeTransition,dt_construction_transition,tlim=np.Inf,CCAs=None,solution=None):
            super().__init__(constellation,start_date,modeleDeTransition,dt_construction_transition,tlim=tlim,CCAs=CCAs,solution=solution)
            self.objectif = (0,0)
            self.initRequetes(constellation,initFirstMode=False)
            
        def resetModeSolver(self,constellation):
            printOpen("Génération des modes")
            self.model = CpoModel()
            modes = []
            for r in constellation.getRequetes():
                printOpen("Durée de génération des modes pour la requête ",r,"("+constellation.getRequete(r).getType()+")")
                for m in constellation.getRequete(r).getKMeilleursModes(constellation,config.getOptValue("modes")):
                    modes.append((r,m.getId()))
                printClose()
            printClose()
            #self.solution.historique.tracerImplicationObs(self.grapheDependances,constellation,tri_cca=True)
            return modes
        
        def initVars(self,constellation):
            printOpen("Initialisation des variables")
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
            printClose()
            
        def initConstraints(self,constellation):
            printOpen("Initialisation des contraintes")
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
                            transition = constellation.getSatellite(s).getTransition(a1,a2,self.modeleDeTransition)
                            distance_matrix[i][j] = int(math.ceil(transition))
                mat = self.model.transition_matrix(distance_matrix)
                seq = self.model.sequence_var(vars_list)
                self.model.add(self.model.no_overlap(seq,distance_matrix=mat))
            printClose()
        
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
        
        def remplirSolCCAs(self,constellation,sol):
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
                self.getSolCCA(s,cca).setSequence(constellation,seq,self.modeleDeTransition)
            
        def resoudre(self,constellation,mailbox,afficher=True):
            temps = time()
            self.initIteration(constellation)
            printColor("Résolution ...")
            if config.getOptValue("time")<np.Inf:
                tlim = math.floor(self.getTempsRestant())
            else:
                tlim = np.Inf
            printColor("Temps restant pour la résolution : "+str(tlim),"(s)")
            if tlim>0:
                with open("test_output",'w') as file:
                    if config.getOptValue("threads")is not None:
                        sol = self.model.solve(TimeLimit=tlim,execfile=config.glob.docplexPath,Workers=config.getOptValue("threads"),log_output=file)
                    else:
                        if config.getOptValue("time")<np.Inf:
                            sol = self.model.solve(TimeLimit=tlim,execfile=config.glob.docplexPath,log_output=file)
                        else:
                            sol = self.model.solve(execfile=config.glob.docplexPath,log_output=file)
                    self.solution.historique.gap = sol.solution.get_objective_gap()
                    if str(sol.get_solve_status())!="Unknown":
                        self.setModesRetenus([(-1,0)]+[(r,m) for (r,m) in self.mode_var if sol.get_value(self.mode_var[(r,m)])==1],constellation)
                        self.remplirSolCCAs(constellation,sol)
                        #self.verifierSolution(constellation)
            else:
                raise Exception("pas le temps de lancer le solver")
            self.notifierFinExecution(constellation)
            self.afficherInfo(time(),self.start_date,constellation,add={'gap':self.solution.historique.gap},color='y',title='FIN')
                
    class Master:
        def __init__(self,constellation,start_date,modeleDeTransition,dt_construction_transition,tlim=np.Inf,CCAs=None,solution=None):
            self.start_date = start_date
            self.initSolution(constellation,tlim,modeleDeTransition,dt_construction_transition,CCAs=CCAs,solution=solution)
            
        def resoudre(self,constellation,afficher=True):
            # iterer insertion-réparation
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
            
        def initSolution(self,constellation,tlim,modeleDeTransition,dt_construction_transition,CCAs=None,solution=None):
            self.solution = CPSolver(constellation,self.start_date,modeleDeTransition,dt_construction_transition,tlim=tlim,CCAs=CCAs,solution=solution)
            
        def getSolution(self):
            return self.solution.getSolution()

        def verifierSolution(self):
            self.solution.verifierSolution()
            
        def getModesSolution(self):
            return self.solution.getModesRetenus()
    
    

    class runnableCPSolver:
        def execute(self,constellation,start_date,modeleDeTransition,dt_construction_transition,tlim=np.Inf,CCAs=None,solution=None):
            coeurs = MPI.COMM_WORLD.Get_size()
            config.setOptValue("threads",coeurs)
            if MPI.COMM_WORLD.Get_rank()==0:
                self.process = Master(constellation,start_date,modeleDeTransition,dt_construction_transition,tlim=tlim,CCAs=CCAs,solution=solution)
                self.process.resoudre(constellation)
                for i in range(1,MPI.COMM_WORLD.Get_size()):
                    MPI.COMM_WORLD.send({"sol":self.process.solution.getSolutionContainer()},dest=i)
                if config.verifMode():
                    self.process.solution.verifierSolution(constellation)
                return self.process.solution.getSolutionContainer()
            else:
                activites = constellation.extraireActivitesRequetesActives(None)
                if CCAs is not None:
                    self.grapheDependances = CCAs
                else:
                    self.grapheDependances = GroupeComposantesActivitesStatiques(constellation,activites,modeleDeTransition)

                data = None
                data = MPI.COMM_WORLD.recv(data)
                return data['sol']
            
        def getName(self):
            return "CP Solver"
        
        def getMasterSolver(self):
            if MPI.COMM_WORLD.Get_rank()==0:
                return self.process.solution
            else:
                return None
        
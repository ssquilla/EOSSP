from EOSSP.Utils.Utils import *
from EOSSP.Utils.config import *

global config
config = Config()
instance = config.instance
from mpi4py import MPI
if config.getOptValue("help"):
    if MPI.COMM_WORLD.Get_rank()==0:
        config.afficherAide()
else: 
    from EOSSP.model.solution import *
    from EOSSP.model.constellation import *
    from EOSSP.model.components import *
    from EOSSP.model.componentPlan import *
    
    import math
    from time import time
    from time import sleep
    import docplex.cp as cp
    from docplex.cp.model import CpoModel,interval_var,binary_var
    
    
    class CPSolver(Solver):
        def __init__(self,constellation,startDate,transitionModel,dt_construction_transition,tlim=np.Inf,CCAs=None,solution=None):
            super().__init__(constellation,startDate,transitionModel,dt_construction_transition,tlim=tlim,CCAs=CCAs,solution=solution)
            self.objective = (0,0)
            self.initRequests(constellation,initFirstMode=False)
            
        def resetModeSolver(self,constellation):
            printOpen("Generate modes")
            self.model = CpoModel()
            modes = []
            for r in constellation.getRequests():
                printOpen("Generation duration for request ",r,"("+constellation.getRequest(r).getType()+")")
                for m in constellation.getRequest(r).getKBestModes(constellation,config.getOptValue("modes")):
                    modes.append((r,m.getId()))
                printClose()
            printClose()
            return modes
        
        def initVars(self,constellation):
            printOpen("Variables initialization")
            self.intervalVars = {}
            self.modeVars = {}
            self.requestActivities = {}
            for (r,m) in self.modes:
                if r!=-1:
                    self.modeVars[(r,m)] = self.model.binary_var(name="y_"+str((r,m)))
                    self.model.add(self.modeVars[(r,m)])
                for (s,a) in constellation.getRequest(r).getMode(m).getCouples():
                    if a not in self.requestActivities:
                        self.requestActivities[a] = []
                    self.requestActivities[a].append((r,m))
                    if a not in self.intervalVars:
                        start = int(math.floor(constellation.getSatellite(s).getActivity(a).getStartDate()))
                        end = int(math.ceil(constellation.getSatellite(s).getActivity(a).getEndDate()))
                        duration = int(math.ceil(constellation.getSatellite(s).getActivity(a).getDuration()))
                        assert(start+duration<=end)
                        self.intervalVars[a] = self.model.interval_var(start=(start,end),end=(start,end),length=duration,name="Ia_{}_{}".format(s,a))
                        self.model.add(self.intervalVars[a])
            obj1 = sum([self.noiseUtility(constellation,r,m)[0]*self.modeVars[(r,m)] for (r,m) in self.modeVars])
            self.model.maximize(obj1)
            #self.model.maximize(obj2)
            printClose()
            
        def initConstraints(self,constellation):
            printOpen("Constraints initialization")
            # 1 mode par requete
            for r in constellation.getRequests():
                modes = [x for x in self.modeVars if x[0]==r]
                if len(modes)>0:
                    self.model.add(sum([self.modeVars[x] for x in modes])<=1)
            # presence des activites = presence d'un mode
            for a in self.intervalVars:
                #print(self.requestActivities[a],a)
                r = self.requestActivities[a][0][0]
                if r != -1:    
                    self.intervalVars[a].set_optional()
                    self.model.add(self.model.presence_of(self.intervalVars[a]) == sum([self.modeVars[x] for x in self.requestActivities[a]]))
                else:
                    self.model.add(self.model.presence_of(self.intervalVars[a]) == 1)
            self.getActivitiesByCCA()
            # flot sur les cca
            for idCCA in self.ccaActivities:
                s,cca = idCCA
                n = len(self.ccaActivities[idCCA])
                distance_matrix = np.zeros((n,n),dtype=int)
                vars_list = []
                for i,a1 in enumerate(self.ccaActivities[idCCA]):
                    vars_list.append(self.intervalVars[a1])
                    for j,a2 in enumerate(self.ccaActivities[idCCA]):
                        if i==j:
                            distance_matrix[i][j] = 0
                        else:
                            transition = constellation.getSatellite(s).getTransitionDuration(a1,a2,self.transitionModel)
                            distance_matrix[i][j] = int(math.ceil(transition))
                mat = self.model.transition_matrix(distance_matrix)
                seq = self.model.sequence_var(vars_list)
                self.model.add(self.model.no_overlap(seq,distance_matrix=mat))
            printClose()
        
        def getActivitiesByCCA(self):
            mappingActivitiesByCCA = {}
            for a in self.intervalVars:
                idCCA = self.getGraphComponents().getActivityCCA(a)
                if idCCA not in mappingActivitiesByCCA:
                    mappingActivitiesByCCA[idCCA] = []
                mappingActivitiesByCCA[idCCA].append(a)
            self.ccaActivities = mappingActivitiesByCCA
        
        def initIteration(self,constellation):
            self.setAllSolCCAs ({s : {} for s in constellation.satellites}) # liste de solutions CCA
            for s,cca in self.getGraphComponents().getComponents():
                self.setSolCCA(s,cca, SolCCA(cca,s))
            self.modes = self.resetModeSolver(constellation)
            self.initVars(constellation)
            self.initConstraints(constellation)
    
        """
            =============================================== 
                            RESOLUTION
            =============================================== 
        """
        
        def insertSequences(self,constellation,sequences):
            fulfilledRequests = [x[0] for x in self.getSelectedModes()]
            for s,cca in sequences:
                seq = [a for a in sequences[(s,cca)] if constellation.getRequestActivite(a) in fulfilledRequests]
                #print("set sequence",s,cca)
                self.getSolCCA(s,cca).setSequence(constellation,seq)
        
        def fillSolCCAs(self,constellation,sol):
            obsCCA = {}
            for a in self.intervalVars:
                r = self.requestActivities[a][0][0]
                if r==-1 or sum([sol.get_value(self.modeVars[x]) for x in self.requestActivities[a]])==1:
                    idCCA = self.getGraphComponents().getActivityCCA(a)
                    if idCCA not in obsCCA:
                        obsCCA[idCCA] = []
                    val = sol.get_value(self.intervalVars[a])
                    obsCCA[idCCA].append((val[0],a))
            for idCCA in obsCCA:
                s,cca = idCCA
                seq = [x[1] for x in sorted(obsCCA[idCCA])] 
                self.getSolCCA(s,cca).setSequence(constellation,seq,self.transitionModel)
            
        def resolve(self,constellation,mailbox,afficher=True):
            temps = time()
            self.initIteration(constellation)
            printColor("Solving ...")
            if config.getOptValue("time")<np.Inf:
                tlim = math.floor(self.getTempsRestant())
            else:
                tlim = np.Inf
            printColor("Time left: "+str(tlim),"(s)")
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
                        self.setSelectedModes([(-1,0)]+[(r,m) for (r,m) in self.modeVars if sol.get_value(self.modeVars[(r,m)])==1],constellation)
                        self.fillSolCCAs(constellation,sol)
                        #self.verifySolution(constellation)
            else:
                raise Exception("not enough time to run the solver")
            self.notifyEndExecution(constellation)
            self.displayInformation(time(),self.startDate,constellation,add={'gap':self.solution.historique.gap},color='y',title='END')
                
    class Master:
        def __init__(self,constellation,startDate,transitionModel,dt_construction_transition,tlim=np.Inf,CCAs=None,solution=None):
            self.startDate = startDate
            self.initSolution(constellation,tlim,transitionModel,dt_construction_transition,CCAs=CCAs,solution=solution)
            
        def resolve(self,constellation,afficher=True):
            # iterer insertion-réparation
            self.solution.resolve(constellation,afficher)
            self.solution.buildPlan(constellation)
    
        def getBestSolution(self):
            return self.solution.getBestSolution()
        
        def registerObjective(self):
            rel = config.glob.releves
            sol = self.getBestSolution()
            points = []
            i_rel = 0
            for i,(t,obj,modes) in enumerate(sol):
                if i==len(sol)-1:
                    points.append((t//60,t%60,obj,modes))
                elif t>rel[i_rel]:
                    points.append((t//60,t%60,obj,modes))
                    i_rel += 1
            return points               
        
        def plotActivity(self,constellation,annoter=False):
            return self.solution.plotActivity(constellation,annoter)
            
        def plotHistory(self,init=True):
            return self.solution.plotHistory(init)
        
        def saveSample(self,constellation):
            self.solution.saveSample(constellation)
            
        def initSolution(self,constellation,tlim,transitionModel,dt_construction_transition,CCAs=None,solution=None):
            self.solution = CPSolver(constellation,self.startDate,transitionModel,dt_construction_transition,tlim=tlim,CCAs=CCAs,solution=solution)
            
        def getSolution(self):
            return self.solution.getSolution()

        def verifySolution(self):
            self.solution.verifySolution()
            
        def getModesSolution(self):
            return self.solution.getSelectedModes()
    
    

    class runnableCPSolver:
        def execute(self,constellation,startDate,transitionModel,dt_construction_transition,tlim=np.Inf,CCAs=None,solution=None):
            coeurs = MPI.COMM_WORLD.Get_size()
            config.setOptValue("threads",coeurs)
            if MPI.COMM_WORLD.Get_rank()==0:
                self.process = Master(constellation,startDate,transitionModel,dt_construction_transition,tlim=tlim,CCAs=CCAs,solution=solution)
                self.process.resolve(constellation)
                for i in range(1,MPI.COMM_WORLD.Get_size()):
                    MPI.COMM_WORLD.send({"sol":self.process.solution.getSolutionContainer()},dest=i)
                if config.verifMode():
                    self.process.solution.verifySolution(constellation)
                return self.process.solution.getSolutionContainer()
            else:
                activites = constellation.extraireActivitesRequetesActives(None)
                if CCAs is not None:
                    self.dependenciesGraph = CCAs
                else:
                    self.dependenciesGraph = StaticCCAsGroup(constellation,activites,transitionModel)

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
        
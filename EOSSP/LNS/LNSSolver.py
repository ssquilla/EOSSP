from EOSSP.Utils.config import *
global config
config = Config()
instance = config.instance
from mpi4py import MPI
if config.getOptValue("help"):
    if MPI.COMM_WORLD.Get_rank()==0:
        config.displayHelp()
else:
    
    from EOSSP.LNS import operateursLNS as op
    from EOSSP.LNS.LNSTools import *
    
    from EOSSP.model.solution import *
    from EOSSP.model.constellation import *
    from EOSSP.model.components import *
    from EOSSP.model.componentPlan import *
    #from ..model.oracle import Oracle
    
    from EOSSP.Utils.Utils import *
    from EOSSP.Utils.Utils import printColor,warn,alert,printOpen,printClose,printMaster,choseAndBroadcastFile
    from EOSSP.Utils.Utils import getDisplayDepth,shiftLeftDisplay,shiftRightDisplay
    
    from time import time,sleep
    import math
    import numpy as np
    from copy import deepcopy,copy
    
    class LNS(Solver):
        
        class ActivitiesHeuristic:
            def __init__(self,nom):
                self.nom = nom
            
            def pickActivity(self,constellation,r,activitiesToInsert):
                raise ValueError('abstract class')
        
        class RatioHeuristic(ActivitiesHeuristic):
            def __init__(self,nom,alpha):
                super().__init__(nom)
                self.alpha = alpha
                
            def scoreActivitySet(self,LNSsolver,constellation,r,activitiesToInsert,cheap=False):        
                #printOpen("Scoring des activités",c='b')
                timeSlots = {a : constellation.getRequest(r).getTimeSlot(a) for a in activitiesToInsert}
                listeObs = [constellation.getRequest(r).findActivite(a) for a in activitiesToInsert]
                reward,temporalScore = constellation.getRequest(r).scoreObservationsList(listeObs,constellation,timeSlots)
                duration = sum([constellation.getSatellite(constellation.getSatelliteActivity(a)).getActivity(a).getDuration() for a in activitiesToInsert])
                if not cheap and self.alpha!=0:
                    trans = LNSsolver.requestState[r].diffTransitionCost(LNSsolver,constellation,activitiesToInsert)
                else:
                    trans = 0
                res = (reward,temporalScore,duration+trans)
                return self.aggregateNodes(res)
            
            def aggregateNodes(self,lambda_v):
                return lambda_v[0]/(1+lambda_v[2])**self.alpha
                    
        class RequestState:
            def __init__(self,LNS,constellation,r):
                self.explainations = []
                self.inactives = []
                self.r = r
                idMPI = MPI.COMM_WORLD.Get_rank()
                configuredSeed = config.getOptValue("seed")
                self.randominitInsertTmp = rd.Random(idMPI+configuredSeed)
                self.initModeResearch(LNS,constellation)
            def diffTransitionCost(self,LNS,constellation,activitiesToInsert):
                cout = 0
                for idCCA in self.tmp[activitiesToInsert]:
                    s = idCCA[0]
                    cout += (self.tmp[activitiesToInsert][idCCA].getTotalTransitionCostOfComponent(constellation) - LNS.getSolCCA(s,cca).getTotalTransitionCostOfComponent(constellation))
                return cout
            def getIdRequest(self):
                return self.r            
            def initModeResearch(self,LNS,constellation):
                self.tmp = {} # les solutions temporaires des cca
                r = self.r
                self.inactives = []
                self.feasibleActivities = {}
                # trouver le meilleur mode privé des activités inactives
                mode = constellation.getRequest(r).getBestModeWithoutInactives(constellation,self.inactives)
                if mode is not None:
                    self.currentScore = mode.getUtility()
                else:
                    self.currentScore = 0
                return self.inactives
            def getCurrentScore(self):
                return self.currentScore
            
            def getActivityByCCA(self,activites,LNS):
                activityByCCA = {}
                for (s,a) in activites:
                    idCCA = LNS.getGraphComponents().getActivityCCA(a)
                    if idCCA not in activityByCCA:
                        activityByCCA[idCCA] = []
                    activityByCCA[idCCA].append(a)
                for idCCA in activityByCCA:
                    activityByCCA[idCCA] = tuple(activityByCCA[idCCA])
                return activityByCCA
                
            def tryInsertRequest(self,LNS,constellation,r,forbidSolver,transitionModel):
                #shiftRightDisplay(1)
                self.solution = {}
                #print(r,sorted(self.inactives))
                mode = constellation.getRequest(r).getBestModeWithoutInactives(constellation,self.inactives)
                if mode is None:
                    printColor("Requete "+str(self.r) +" : aucun mode candidat trouvé",c='r')
                    return None,None,True  
                self.currentScore = mode.getUtility()
                return self.tryInsertMode(LNS,constellation,r,mode,forbidSolver,transitionModel)
            
            def saveActivitiesState(self):
                self.inactivesCopy = copy(self.inactives)
                self.feasibleActivitiesCopy = copy(self.feasibleActivities)
            
            def restoreActivitiesState(self):
                self.inactives = self.inactivesCopy
                supp = []
                for a in self.feasibleActivities:
                    if not self.feasibleActivities[a]:
                        if a not in self.feasibleActivitiesCopy or self.feasibleActivitiesCopy[a]:
                            supp.append(a)
                for a in supp:
                    del self.feasibleActivities[a]
                    del self.tmp[a]
            
            def tryInsertMode(self,LNS,constellation,r,mode,forbidSolver,transitionModel):
                #shiftRightDisplay(1)
                self.solution = {}
                if mode is None:
                    printColor("Request "+str(self.r) +": no path found",c='r')
                    return None,None,True  
                self.currentScore = mode.getUtility()
                realisable = True
                #shiftRightDisplay(2)
                printOpen("Testing path feasibility (+ precise score evaluation)")
                activityCCA = self.getActivityByCCA(mode.getPairs(),LNS)
                for idCCA in activityCCA:
                    faisable = self.testFeasibility(constellation,LNS,r,activityCCA[idCCA],idCCA,forbidSolver,transitionModel)
                    if not faisable:
                        realisable = False
                        self.inactives+=list(activityCCA[idCCA])
                        self.feasibleActivities[activityCCA[idCCA]] = False
                        printColor("Request "+str(self.r) +": infeasible mode",c='r')
                        printClose()
                        return None,None,False
                printClose()
                assert(realisable)
                if config.verifMode():
                    for idCCA in activityCCA:
                        s,cca = idCCA
                        assert(LNS.getSolCCA(s,cca).isSequenceFeasible(constellation,transitionModel))
                printColor("Request "+str(self.r) +" feasible mode: ",mode,c='g')
                mode = constellation.getRequest(r).validateCandidateMode(constellation)

                for idCCA in activityCCA:
                    CCAS = list(self.tmp[activityCCA[idCCA]].keys())
                    assert(len(CCAS)==1)
                    if idCCA in self.solution:
                        die("CCA already present",idCCA)
                    self.solution[idCCA] = self.tmp[activityCCA[idCCA]][idCCA]
                return mode,self.solution,False    
    
            def duplicatedCCA(self,noeuds):
                doublons = []
                CCAs = [list(self.tmp[x].keys()) for x in noeuds]
                for CCAList in CCAs:
                    if not(len(CCAList)==1):
                        print(CCAList,CCAs)
                for i,cca1 in enumerate(CCAs):
                    for j,cca2 in enumerate(CCAs):
                        if i!=j and cca1==cca2 and cca1[0] not in doublons:
                            doublons.append(cca1[0])
                return doublons
            
            def testFeasibility(self,constellation,LNS,r,X,idCCA,forbidSolver,transitionModel):
                if X not in self.feasibleActivities:
                        assert(X not in self.tmp)
                        self.feasibleActivities[X] = self.temporaryInsertion(LNS,constellation,r,X,forbidSolver,transitionModel)
                        assert(X in self.tmp)
                        for a in self.tmp[X][idCCA].sequence:
                            if not LNS.getGraphComponents().getActivityCCA(a)==idCCA:
                                print(X,idCCA,self.tmp[X][idCCA].sequence,a,"process",MPI.COMM_WORLD.Get_rank())
                            assert(LNS.getGraphComponents().getActivityCCA(a)==idCCA)
                return self.feasibleActivities[X]
            def checkActivitiesCCA(self,LNS,activitiesToInsert):
                ccaActivities = {}
                for a in activitiesToInsert:
                    idCCA = LNS.getGraphComponents().getActivityCCA(a)
                    if idCCA not in ccaActivities:
                        ccaActivities[idCCA] = []
                    ccaActivities[idCCA].append(a)
                idCCA = list(ccaActivities.keys())[0]
                assert(len(ccaActivities)<=1)
                return ccaActivities,idCCA
            def initTemporaryInsertion(self,LNS,constellation,r,activitiesToInsert):
                self.tmp[activitiesToInsert] = {}
                ccaActivities,idCCA = self.checkActivitiesCCA(LNS,activitiesToInsert)
                # plan earliest,latest
                s,cca = idCCA
                groups = {0 : copy(LNS.getSolCCA(s,cca).getSequence())}
                #if(not LNS.transitionModel.isTimeDependent()):
                try:
                    LNS.computeCriticalPlans(constellation,idCCA)
                except NumericalError:
                    assert(LNS.transitionModel.isTimeDependent())
                self.tmp[activitiesToInsert][idCCA] = LNS.getSolCCA(s,cca).getCopy()
                lengthBeforeInsertion = len(self.tmp[activitiesToInsert][idCCA].sequence)
                cp = ccaActivities[idCCA].copy()
                self.randominitInsertTmp.shuffle(cp)
                ccaActivities[idCCA] = cp
                for a in self.tmp[activitiesToInsert][idCCA].sequence:
                    assert(idCCA==LNS.getGraphComponents().getActivityCCA(a))
                return lengthBeforeInsertion,groups,s,idCCA,ccaActivities
        
            def performTemporaryGreedyInsertion(self,constellation,ccaActivities,activitiesToInsert,idCCA,solverUsedAfter,transitionModel):
                if config.verifMode():
                    assert(self.tmp[activitiesToInsert][idCCA].isSequenceFeasible(constellation,transitionModel))
                feasible = True
                printOpen("insertion of activities",activitiesToInsert," (greedy method) in CCA",idCCA,c='c')
                for i,p in enumerate(ccaActivities[idCCA]):
                    printOpen("greedy insertion")
                    if feasible:
                        assert(self.tmp[activitiesToInsert][idCCA].isSequenceFeasible(constellation,transitionModel))
                    #if not transitionModel.isTimeDependent() and feasible:
                        if self.tmp[activitiesToInsert][idCCA].arePlansUpToDate():
                            insert = self.tmp[activitiesToInsert][idCCA].insertActivityTCriticalPlansMethod(constellation,p,transitionModel)
                        else:
                            assert(transitionModel.isTimeDependent())
                            insert = self.tmp[activitiesToInsert][idCCA].insertActivityTimeConsumingMethod(constellation,p,transitionModel)
                    else:
                        insert = self.tmp[activitiesToInsert][idCCA].insertActivityTimeConsumingMethod(constellation,p,transitionModel)
                    if insert:
                        assert(self.tmp[activitiesToInsert][idCCA].isSequenceFeasible(constellation,transitionModel))
                    
                    feasible = insert and feasible
                    printClose()
                    if not feasible and not solverUsedAfter:
                        break
                    printOpen("Update critical plans")
                    if feasible:
                        try:
                            self.computeCriticalPlans(activitiesToInsert,idCCA,constellation,transitionModel)
                        except NumericalError:
                            assert(transitionModel.isTimeDependent())
                            
                    printClose()
                printClose()
                return feasible        
            def computeCriticalPlans(self,activitiesToInsert,idCCA,constellation,transitionModel):    
                if config.getOptValue("use_solver"):
                    self.tmp[activitiesToInsert][idCCA].updateCriticalPlans(constellation,transitionModel)
                else: # a remplacer par une propagation
                    self.tmp[activitiesToInsert][idCCA].updateCriticalPlans(constellation,transitionModel)
            def resetCCA(self,idCCA,grapheDep):
                del_x = []
                for X in self.tmp:
                    if idCCA in self.tmp[X]:
                        del_x.append(X)
                for x in del_x:
                    del self.tmp[x]
                    if x in self.feasibleActivities:
                        del self.feasibleActivities[x]
                
                activitiesToRelease = grapheDep.getActivitiesOfComponent(idCCA)
                for a in self.inactives:
                    if a in activitiesToRelease:
                        self.inactives.remove(a)
                return len(del_x)
                
            def temporayInsertionSolver(self,constellation,LNS,ccaActivities,feasible,lengthBeforeInsertion,activitiesToInsert,idCCA,groups,forbidSolver,transitionModel):
                try:
                    copy_sequence = self.tmp[activitiesToInsert][idCCA].getSequence()
                    assert(len(activitiesToInsert)>0)
                    printOpen("local search",c='c')
                    groups[1] = ccaActivities[idCCA]
                    assert(len(groups[1])>0)
                    for a in self.tmp[activitiesToInsert][idCCA].sequence:
                        assert(LNS.getGraphComponents().getActivityCCA((a)==idCCA))
                    self.tmp[activitiesToInsert][idCCA].localSearch(constellation,'OPTW',transitionModel,groups=groups,allowNoSolution=True)#'OPTWGroups'?
                    for a in self.tmp[activitiesToInsert][idCCA].sequence:
                        assert(LNS.getGraphComponents().getActivityCCA((a)==idCCA))
                    assumedLength = lengthBeforeInsertion + len(activitiesToInsert)
                    feasible = (assumedLength == len(self.tmp[activitiesToInsert][idCCA].sequence))
                    printClose()
                    return feasible
                except InfeasibleSolutionfromSolverException:
                    if config.verifMode():
                        print("Solver error: infeasible solution.")
                    printClose()
                    return False
        
            # insere les activites dans tmp (copies des sol CCA) et renvoie la liste des cca modifiées
            # sur un noeud : une seule cca est modifiee
            def temporaryInsertion(self,LNS,constellation,r,activitiesToInsert,forbidSolver,transitionModel):
                lengthBeforeInsertion,groups,s,idCCA,ccaActivities = self.initTemporaryInsertion(LNS,constellation,r,activitiesToInsert)
                isSolverUsed = config.getOptValue("use_solver") and not forbidSolver
                feasible = False
                feasible = self.performTemporaryGreedyInsertion(constellation,ccaActivities,activitiesToInsert, idCCA, isSolverUsed,transitionModel)
                assert(feasible is not None)
                if not feasible and isSolverUsed:
                    feasible = self.temporayInsertionSolver(constellation,LNS,ccaActivities, feasible, lengthBeforeInsertion, activitiesToInsert, idCCA, groups,forbidSolver,transitionModel)
                return feasible
        
        class SolutionSaver:
            def __init__(self,LNS,constellation):
                self.solCCAs = {s : {cca : LNS.getSolCCA(s,cca).getCopy() for cca in LNS.getSolCCAs()[s]} for s in LNS.getSolCCAs()}
                self.selectedModes = LNS.getSelectedModes().copy()
                self.requestState = deepcopy(LNS.requestState)
                self.objective = LNS.getObjective(constellation,recompute=True)
            
            def backupSolution(self,LNS):
                LNS.setAllSolCCAs(self.solCCAs)
                LNS.overwriteModes(self.selectedModes)
                LNS.setObjective(self.objective)
                LNS.requestState = self.requestState
        
        def __init__(self,constellation,startDate,transitionModel,dt_construction_transition,tlim=np.inf,CCAs=None,solution=None):
            super().__init__(constellation,startDate,transitionModel,dt_construction_transition,tlim=tlim,CCAs=CCAs,solution=solution)
            shift = getDisplayDepth()-1
            shiftLeftDisplay(shift)
            #self.explication = {}
            # Initialisation des requetes
            self.initRequests(constellation)
            # composantes et solutions de composantes
            self.oracle = {}
            self.ccaPositions = {}
            self.ccaObservations = []
            self.nDegradation = 0
            self.nNoDegradation = 0
            self.tlim = min(self.tlim,config.getOptValue("time"))
            # cca qui contiennent au moins une obs (écarter les cca avec qu'un téléchargement)
            for idCCA in self.getGraphComponents().getComponents():
                s,cca = idCCA
                self.ccaPositions[idCCA] = np.mean(np.array([constellation.getSatellite(s).getActivity(a).getCoordinates() for a in self.getGraphComponents().getActivitiesOfComponent(idCCA)]))
                for a in self.getGraphComponents().getActivitiesOfComponent(idCCA):
                    if constellation.getSatellite(constellation.getSatelliteActivity(a)).isObservation(a):
                        self.ccaObservations.append(idCCA)
                        break
                #self.oracle[idCCA] = Oracle(idCCA)
            self.updatePresentRequests(constellation)
            # si 1 coeur : slave local
            comm = MPI.COMM_WORLD
            comm.Barrier()
            self.notifyPreprocessing(constellation)
            self.initOperators()
            self.initRandomizers()
            shiftRightDisplay(shift)

        def initRandomizers(self):
            idMPI = MPI.COMM_WORLD.Get_rank()
            configuredSeed = config.getOptValue("seed")
            self.requestChoiceRandomizer = rd.Random(idMPI+configuredSeed)
            
        # retourne les CCA qui contiennent au moins une observation
        def getCCAObservations(self):
            return self.ccaObservations
        
        def distanceCCA(self,cca1,cca2):
            return math.dist(self.ccaPositions[cca1]-self.ccaPositions[cca2])
        
        """
            =============================================== 
                            RESOLUTION
            =============================================== 
        """
        def initRequests(self,constellation):
            super().initRequests(constellation,initFirstMode=False)           
    
        def insertSequences(self,constellation,sequences):
            selectedRequests = [x[0] for x in self.getSelectedModes()]
            for idCCA in sequences:
                s,cca = idCCA
                seq = [a for a in sequences[idCCA] if constellation.getRequestActivity(a) in selectedRequests]
                self.getSolCCA(s,cca).setSequence(constellation,seq)
    
        def createBackup(self,constellation):
            return self.SolutionSaver(self,constellation)
                
        def LNSsolve(self,constellation,mailbox,lastRecord,it):
            filtre = ['time','obj','modes','requetes','best']
            operatorCalls = 0
            stableIt = 0
            while stableIt<config.getOptValue("stableIt") and operatorCalls<config.getOptValue("max_operateur") and time()-self.startDate<self.tlim:
                if config.getOptValue("verif"):
                    self.verifySolution(constellation)
                if config.getOptValue("dynamic"):
                    change,requestList = constellation.releaseNewRequests(self.getGraphComponents())
                    if change:
                        self.updateNewRequests(constellation,time()-self.startDate,requestList)
                        self.operator.updateNewRequests(self)
                improvment = self.applyOperator(constellation)
                self.notifyEndOperator(constellation)
                if not improvment:
                    stableIt += 1
                else:
                    stableIt = 0
                operatorCalls += 1
            self.notifyEndIteration(constellation)
            if time()-lastRecord>=config.glob.periode_affichage_info:
                lastRecord = time()
                self.displayInformation(time(),self.startDate,constellation,title="RESULT ITERATION " +str(it)+" CPU "+str(MPI.COMM_WORLD.Get_rank()),filtre=filtre,color='c',add={"number of calls to the operator'":operatorCalls})
            return operatorCalls,lastRecord
        
        def initOperators(self):
            self.operators = []
            # choix de la méthode de perturbation
            if config.getOptValue("version") == "greedy-request":
                perturb = config.getOptValue("perturb_rate")
                destroy = config.getOptValue("destroy")
                acceptWorseSolution=True
                self.perturbator = op.DestroyAndRepairGreedyRequest(self.transitionModel,perturb,acceptWorseSolution)
                acceptWorseSolution=False
                self.operator = op.DestroyAndRepairGreedyRequest(self.transitionModel,destroy,acceptWorseSolution)
            elif config.getOptValue("version") == "greedy-cca":
                perturb = config.getOptValue("perturb_rate")
                acceptWorseSolution=True
                self.perturbator = op.DestroyAndRepairGreedyRequest(self.transitionModel,perturb,acceptWorseSolution)
                acceptWorseSolution=False
                self.operator = op.DestroyAndRepairGreedyCCA(self,transitionModel,acceptWorseSolution,k=config.getOptValue("n_cca"))
            elif config.getOptValue("version") in ["hybrid","coop"]:
                perturb = config.getOptValue("perturb_rate")
                acceptWorseSolution=True
                self.perturbator = op.DestroyAndRepairGreedyRequest(self.transitionModel,perturb,acceptWorseSolution)
                self.operator = op.VoisinageCP(MPI.COMM_WORLD.Get_size(),self.transitionModel,k=config.getOptValue("n_cca"))
            else:
                raise NameError("Unknown perturbation method")
                
        def notifyCCAChangement(self,idCCA):
            ccaTmpSupp = 0
            for r in self.requestState:
                ccaTmpSupp += self.requestState[r].resetCCA(idCCA,self.getGraphComponents())
            return ccaTmpSupp
        
        def checkSequence(self,constellation):
            for s in self.getSolCCAs():
                for cca in self.getSolCCAs():        
                    for a in self.getSolCCA(s,cca).getSequence() :
                        assert(s==constellation.getSatelliteActivity(a))
        
        def overwriteCCASolution(self,constellation,idCCA,sequence):
            s,cca = idCCA
            self.getSolCCA(s,cca).setSequence(constellation,sequence,self.transitionModel)
            return self.notifyCCAChangement(idCCA)
            
        def createInitialSolution(self,constellation):
            repetitions = 1 # determiste désormais donc inutile de chercher plus loin
            self.greedyFill(constellation,requestsLimit=False,forbidSolver=True)
            if config.getOptValue("use_solver"):
                self.resetRequestState(constellation)
            # quand on passe d'un greedy sans solver à greedy avec solver il faut reset les inactives
            return repetitions

        def initRequestState(self,constellation):
            self.requestState = {}
            for r in constellation.getAllRequests():
                self.requestState[r] = self.RequestState(self,constellation,r)
            
        def resetRequestState(self,constellation):
            for r in self.requestState:
                self.requestState[r].initModeResearch(self,constellation)
        
        def saveActivitiesState(self):
            for r in self.requestState:
                self.requestState[r].saveActivitiesState()
                
        def restoreActivitiesState(self):
            for r in self.requestState:
                self.requestState[r].restoreActivitiesState()
                
        def resolve(self,constellation,mailbox,afficher=True):
            comm = MPI.COMM_WORLD
            size = comm.Get_size()
            #filtre = ['time','obj','modes','requetes','best'] # infos à afficher
            it=0
            lastRecord = time()
            operatorsApplicationCounter = 0
            nPerturbations = 0
            self.initRequestState(constellation)
            self.createInitialSolution(constellation)
            self.notifyNoEvent(constellation)
            self.sumDeltaPertubation = 0
            while time()-self.startDate<self.tlim and it<config.getOptValue("max_iteration"):
                it += 1 
                if time()-self.startDate<self.tlim and it<config.getOptValue("max_iteration"):
                    operatorCalls,lastRecord = self.LNSsolve(constellation,mailbox,lastRecord,it)
                    operatorsApplicationCounter += operatorCalls
                    utilityBefore = self.getSolutionContainer().objective[0]
                    if time()-self.startDate<self.tlim:
                        printOpen("Perturb",c='y')
                        self.perturbator.apply(self,constellation)
                        nPerturbations += 1
                        printClose()
                    if config.getOptValue("verif"):
                        self.verifyCCAs(constellation)
                    utilityAfter = self.getSolutionContainer().objective[0]
                    self.sumDeltaPertubation += (utilityAfter-utilityBefore)
                    if utilityBefore>utilityAfter:  
                         self.nDegradation += 1
                    else:
                        self.nNoDegradation += 1  
            try:
                self.meanPertubation = self.sumDeltaPertubation/(self.nDegradation+self.nNoDegradation)           
            except:
                self.meanPertubation = 0
            self.verifySolutionIfVerifyMode(constellation)
            objective,solCCAs,selectedModes = self.getBest()
            self.setAllSolCCAs(solCCAs)
            self.setSelectedModes(selectedModes,constellation)
            self.notifyEndExecution(constellation)
            return operatorsApplicationCounter,it
                
        def restart(self,constellation):
            self.initRequestState(constellation)
            for r in self.requestState:
                self.requestState[r].initModeResearch(self,constellation)
            super().restart(constellation)
    
        def isStable(self,objective):
            return objective[0]<=self.solution.history.getBest()[0][0]
        
        def isActive(self,lambda_v):
            return lambda_v[3]<=config.LNS.Tmax
    
        def repair(self,constellation,destroyedRequests,entryOpt):
            printOpen("repair: "+str(self.operatorReparation),c='b')
            requests = self.operatorReparation.repair(self,constellation,destroyedRequests,entryOpt)
            printColor(str(len(requests))+" rebuilt requests",c='g')
            printClose()
            
        def destroy(self,constellation):
            printOpen("destruction : "+str(self.operatorDestruction),c='c')
            destroyedRequests,sortie_opt = self.operatorDestruction.destroy(self,constellation)
            printClose()
            return destroyedRequests,sortie_opt
        
        def applyOperator(self,constellation):
            printColor('Apply '+str(self.operator),c='c')
            res = self.operator.apply(self,constellation)
            return res
            
        def validateSolution(self,sol):
            ccaTmpSupp = 0
            for idCCA in sol:
                printColor("validate cca",idCCA,c='g')
                s,cca = idCCA
                self.setSolCCA(s,cca,sol[idCCA])
                for r in self.requestState:
                    ccaTmpSupp += self.requestState[r].resetCCA(idCCA,self.getGraphComponents())
            #print(ccaTmpSupp,"cca temporaires supprimées")
    
        # confusion cout de transition et retard
        def aggregateScore(self,score,methode='r1'):
            if methode=='r1':
                return score[0]/(1+score[2])
            elif methode=='r2':
                return score[0]**1.25/(1+score[2])
            elif methode=='r0':
                return score[0]
            else:
                raise ValueError("methode d'aggregation inconnue")
            
        def heuristicReward(self,constellation,s,a):
            return constellation.getSatellite(s).getActivity(a).getScore()/(constellation.getSatellite(s).getActivity(a).getDuration()+1)
        
        def choseRequest(self,candidates,random=False):
            if -1 in candidates:
                return -1
            if len(candidates)>0:
                if not random:
                    return max(candidates,key=lambda r : self.requestState[r].getCurrentScore())
                else:
                    select = self.requestChoiceRandomizer().randInt(len(candidates))
                    return candidates[select]
            else:
                return None
            
        def greedyFill(self,constellation,requestsLimit=True,forbidSolver=False,random=False):
            printOpen("greedy filling")
            candidates = [r for r in self.requestState if constellation.getRequest(r).isActive() and r not in self.fulfilledRequests()]
            nCandidatesRequests = len(candidates)
            succes = 0
            r = self.choseRequest(candidates,random=random)
            nFailures = 0
            counter = 0
            limNotif = 50
            while r is not None and self.getTimeElapsed() < self.tlim and (not requestsLimit or nFailures<config.getOptValue("max_echecs")):
                if time()-self.startDate>self.tlim:
                    printColor("interrupting greedy",c='y')
                    break
                printOpen("Searching new mode for the request",r)
                mode,solutionCCA,vide = self.requestState[r].tryInsertRequest(self,constellation,r,forbidSolver,self.transitionModel)
                if mode is not None:
                    succes += 1
                    self.validateSolution(solutionCCA)
                    self.addSelectedMode((r,mode.getId()),constellation)
                    candidates.remove(r)
                    if config.getOptValue("verif"):
                        self.verifyCCAs(constellation)
                    printClose("Success",c='g')
                else:
                    nFailures += 1
                    if vide:
                        printColor("No more modes",c='m')
                        candidates.remove(r)
                    else:
                        printColor("Infeasible modes",c='m')
                    printClose("Failure",c='r')
                r = self.choseRequest(candidates,random=random)
                if counter==limNotif or r is None:
                    counter = 0
                    self.notifyNoEvent(constellation)
                counter += 1
                #self.notifyNoEvent(constellation)
            printClose(str(succes)+"/"+str(nCandidatesRequests)+" requests inserted. "+str(nFailures)+" failures.",c='g')
            
            
    class Processus:
        def __init__(self,role,constellation,startDate,transitionModel,dt_construction_transition,tlim=np.inf,CCAs=None,solution=None):
            self.role = role
            self.initSolution(constellation,startDate,transitionModel,dt_construction_transition,tlim=tlim,CCAs=CCAs,solution=solution)
            
        def initSolution(self,constellation,startDate,transitionModel,dt_construction_transition,tlim=np.Inf,CCAs=None,solution=None):
            self.solution = LNS(constellation,startDate,transitionModel,dt_construction_transition,tlim=tlim,CCAs=CCAs,solution=solution)
            
        def resolve(self,constellation,mailbox,afficher=True):
            operatorsApplicationCounter,it = self.solution.resolve(constellation,mailbox,afficher)
            if config.getOptValue("version")!="coop":
                MPI.COMM_WORLD.Barrier()
            return operatorsApplicationCounter,it
    
        def bestSolution(self):
            return self.solution.bestSolution()
        
        def recordObjective(self):
            rel = config.glob.releves
            sol = self.bestSolution()
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
            self.solution.saveSample(constellation,add=self.additionnalInfo)
            
        def getSolution(self):
            return self.solution.getSolution()
        
        def setSolution(self,sol,vid,selectedModes,candidatesModes,removedModes,objective):
            self.solution.setSol(sol,vid,selectedModes,candidatesModes,removedModes,objective)
        
        def verifySolution(self):
            self.solution.verifySolution()
            
        def getModesSolution(self):
            return self.solution.getModes()
        
    class Master(Processus):
        def __init__(self,constellation,startDate,transitionModel,dt_construction_transition,tlim=np.Inf,CCAs=None,solution=None):
            super().__init__("master",constellation,startDate,transitionModel,dt_construction_transition,tlim=tlim,CCAs=CCAs,solution=solution)
            
        def resolve(self,constellation,mailbox,afficher=True):
            operatorsApplicationCounter,it = super().resolve(constellation,mailbox,afficher)
            self.collectResults(mailbox)
            #self.terminerProcessus()
            filtre = ['best','obs' ,'requete',"size"]
            self.solution.displayInformation(time(),self.solution.startDate,constellation,color='y',title='END', filtre=filtre)
            self.solution.buildPlan(constellation)
         
        def collectResults(self,mailbox):
            hist = {}
            self.additionnalInfo = {}
            if config.getOptValue("version") in ["hybrid","coop"]:
                gaps,gains = self.solution.operator.getGapsInfos()
            if config.getOptValue("version")!="coop":
                for data in mailbox.readMessages():
                    cpu = data['cpu']
                    hist[cpu] = data['hist']
                    if config.getOptValue("version") in ["hybrid","coop"]:
                        gains += data["gains"]
                        gaps += data["gaps"]
                if config.getOptValue("version") in ["hybrid","coop"] :
                    if len(gains)>0:
                        positiveMeanReward = np.mean([x for x in gains if x>0])
                        successPercentage = len([x for x in gains if x>0])/len(gains)
                    else:
                        positiveMeanReward = 0
                        successPercentage = 0
                    nProblems = len(gains)
                    meanReward = np.mean(gains)
                    self.additionnalInfo["gaps_moyen"] = meanReward
                    self.additionnalInfo["gains"] = positiveMeanReward
                    self.additionnalInfo["fréquence_succès"] = successPercentage
                    self.additionnalInfo["nombre_problèmes"] = nProblems
                self.solution.solution.history.merge(hist)
    
    class Slave(Processus):
        def __init__(self,constellation,startDate,transitionModel,dt_construction_transition,tlim=np.inf,CCAs=None,solution=None):
            comm = MPI.COMM_WORLD
            rank = comm.Get_rank()
            super().__init__("slave "+str(rank),constellation,startDate,transitionModel,dt_construction_transition,tlim=tlim,CCAs=CCAs,solution=solution)
        
        def additionnalInformation(self):
            data = {}
            if config.getOptValue("version") in ["hybrid","coop"]:
                gaps,gains = self.solution.operator.getGapsInfos()
                data["gaps"] = gaps
                data["gains"] = gains
            return data
        
        def formatResult(self):
            add = self.additionnalInformation()
            data = {'hist':self.solution.solution.history,'cpu':MPI.COMM_WORLD.Get_rank()}
            hist = data['hist']
            for key in add:
                data[key] = add[key]
            return data
        
        def resolve(self,constellation,mailbox):
            if config.getOptValue("version")!="coop":
                super().resolve(constellation,mailbox)
                data = self.formatResult()
                mailbox.postMessage(data)
        
        
class runnableLNS:
    def execute(self,constellation,startDate,transitionModel,dt_construction_transition,tlim=np.inf,CCAs=None,solution=None):
        mailbox = BlockingCommunication()
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        if rank == 0:
            self.process = Master(constellation,startDate,transitionModel,dt_construction_transition,tlim=tlim,CCAs=CCAs,solution=solution)
            self.process.resolve(constellation,mailbox)
            
        else:
            self.process = Slave(constellation,startDate,transitionModel,dt_construction_transition,tlim=tlim,CCAs=CCAs,solution=solution)
            self.process.resolve(constellation,mailbox) 
        return self.process.solution.getSolutionContainer()
    
    def getName(self):
        return "Large Neighborhood Search"
    
    def getMasterSolver(self):
        if MPI.COMM_WORLD.Get_rank()==0:
            return self.process.solution
        else:
            return None
        
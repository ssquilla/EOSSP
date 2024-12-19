from EOSSP.Utils.Utils import *
from EOSSP.Utils.config import *
config = Config()

from EOSSP.model.componentPlan import *
from EOSSP.model.history import History
from EOSSP.model.history import PREPROCESSING,NO_EVENT,END_ITERATION,END_RUN,BACKUP,END_OPERATOR

from time import time
from mpi4py import *
from mpi4py import MPI
import pickle

from copy import deepcopy,copy

#import matplotlib
#from matplotlib import patches
#import matplotlib.pyplot as plt
import os

class Solver:
    class Solution:
        def __init__(self,constellation,startDate):
            self.objective = (0,0)
            self.selectedModes = []
            self.solCCAs = {s : {} for s in constellation.satellites} # liste de solutions CCA
            self.history = History(constellation,startDate)
        
        def resetHistory(self,constellation,startDate):
            self.historique = History(constellation,startDate)
            
        def restart(self,constellation):
            self.objective = (0,0)
            self.selectedModes = []
            self.solCCAs = {s : {} for s in constellation.satellites} # liste de solutions CCA
            
        def getSelectedModes(self):
            return self.selectedModes
        
        def getSolCCA(self,s,cca):
            return self.solCCAs[s][cca]
        
        def getSolCCAs(self):
            return self.solCCAs
        
        def updateObjective(self,constellation):
            self.objective = sum([constellation.getRequest(r).getMode(m).getUtility() for (r,m) in self.selectedModes]),sum([constellation.getRequest(r).getMode(m).getTemporalUtility() for (r,m) in self.selectedModes])
            
        def replaceMode(self,r,m,constellation):
            for i,(rr,mm) in enumerate(self.selectedModes):
                if rr==r:
                    self.selectedModes[i] = (r,m)
                    self.updateObjective(constellation)
                    return
            raise ValueError("Requête non présente")
            
        def addSelectedMode(self,mode,constellation):
            r,m = mode
            if config.verifMode():
                assert(r not in [rr for (rr,mm) in self.selectedModes])
            self.selectedModes.append((r,m))
            self.updateObjective(constellation)
            
        def removeRequest(self,r,constellation):
            for i,(rr,mm) in enumerate(self.selectedModes):
                if rr == r:
                    self.selectedModes.pop(i)
                    self.updateObjective(constellation)
                    return
            raise ValueError("requête non présente dans la solution : "+str(r))
    
        def removeSelectedMode(self,mode,constellation):
            self.selectedModes.remove(mode)
            self.updateObjective(constellation)
            
        def setSelectedModes(self,modes,constellation):
            self.selectedModes = modes
            requetes = []
            for (r,m) in modes:
                if(r in requetes):
                    print(modes_a_retirer)
                    raise ValueError("Demande d'ajout de requêtes en double : ",r)
                requetes.append(r)
            
            self.updateObjective(constellation)
        
        def setSolCCAs(self,sol):
            self.solCCAs = sol
        
        def getObjective(self):
            return self.objective
        
        def addInformation(self,cle,valeur):
            self.history.addInformation(cle,valeur)
        
        def __str__(self):
            msg = "---------- Solution ----------\n"
            msg += "| selected requests : "+str(sorted([x[0] for x in self.selectedModes]))+' ' +str(len(self.selectedModes)) + ' modes\n'
            msg += "| objective : " + str(self.objective) +'\n'
            return msg
        
    def __init__(self,constellation,start,transitionModel,dt_construction_transition,tlim=np.Inf,CCAs=None,solution=None):
        self.transitionModel = transitionModel
        self.tlim = min(tlim,config.getOptValue("time"))
        self.shiftDisplay = 1 # paramètre d'affichage
        self.startDate = start
        id_mpi = MPI.COMM_WORLD.Get_rank()
        configuredSeed = config.getOptValue("seed")
        self.randomizerPerturbation = rd.Random(id_mpi+configuredSeed)
        global START_DATE
        START_DATE = self.startDate
        # initialisation des CCAs
        activities = constellation.extractActivitiesActiveRequests(None)
        
        if CCAs is not None:
            self.dependenciesGraph = CCAs
        else:
            self.dependenciesGraph = StaticCCAsGroup(constellation,activities,transitionModel)
            printMaster(self.dependenciesGraph,c='b')
            printMaster("Precomputation duration (CCA etc...) :",time()-self.startDate,"(s)",c='g')
        
        # amorcage solution initiale
        if solution is not None:
            self.solution = solution
            # on indique que les plans earliest latest ne sont peut etre plus a jour
            # => un changement de modele de transition est possible
            for s in self.solution.getSolCCAs():
                for cca in self.solution.getSolCCAs()[s]:
                    self.solution.getSolCCA(s,cca).notifierPlanAJour(False,False)
                    if not transitionModel.isTimeDependent():
                        self.computeCriticalPlans(constellation,(s,cca)) # on remet a jour les plans
            self.verifySolutionIfVerifyMode(constellation,transitionModel)
        else:
            self.solution = self.Solution(constellation,start)
            for (s,cca) in self.dependenciesGraph.getComponents():
                self.setSolCCA(s,cca,SolCCA(cca,s))       
        self.solution.addInformation("transition model building duration: ",dt_construction_transition)
        self.dt_construction_transition = dt_construction_transition
        
        if MPI.COMM_WORLD.Get_rank()==0:
            self.addComment(config.getOptValue("comm"))
        # a voir si on peut se passer de cette liste de candidates
        self.candidateRequests = constellation.getAllRequests()
        # creer les modes courants
        self.currentModes = {}
        if config.getOptValue("dynamic"):
            self.solution.history.registerRequestArrival(0,constellation.getRequestsDepart())
        shiftRightDisplay(self.shiftDisplay)
    
    
    def computeCriticalPlans(self,constellation,id_cca,force=False,print_info=False):
        s,cca = id_cca
        if not self.getSolCCA(s,cca).arePlansUpToDate() or force:
            self.getSolCCA(s,cca).updateCriticalPlans(constellation,self.transitionModel,force=force)
        assert(len(self.getSolCCA(s,cca).getSequence())==len(self.getSolCCA(s,cca).planEst))
        if self.transitionModel.isTimeDependent():
            if not(self.getSolCCA(s,cca).arePlansUpToDate()):
                raise OutdatedCriticalPlans
        
    def updateCriticalPlansSolution(self,constellation,force=False,print_info=False):
        for s in self.solution.solCCAs:
            for cca in self.solution.solCCAs[s]:
                self.computeCriticalPlans(constellation,(s,cca),force=force,print_info=print_info)
            
    def getTimeElapsed(self):
        return time()-self.startDate
    
    def getRemainingTime(self):
        return min(config.getOptValue("time"),self.tlim) - self.getTimeElapsed()
    
    def updateNewRequests(self,constellation,date,requestList):
        self.solution.history.registerRequestArrival(date,requestList)
        self.updatePresentRequests(constellation)
    
    def rejectBestSolution(self):
        self.solution.history.rejectBestSolution()
    
    def getPresentRequests(self,constellation,cca):
        requetes = []
        for a in self.dependenciesGraph.getActivitiesOfComponent(cca):
            r = constellation.getRequestActivity(a)
            if r not in requetes:
                bisect.insort_right(requetes,r)
        return requetes
    
    def updatePresentRequests(self,constellation):
        self.requetes_cca = {}
        for cca in self.dependenciesGraph.getComponents():
            self.requetes_cca[cca] = self.getPresentRequests(constellation,cca)
        self.requetes_en_commun = {}
    
    # les requetes presentes sur la cca meme si les activites ne sont pas dans la solution
    def getPresentRequestsStatic(self,cca):
        return self.requetes_cca[cca]
    
    # chercher les requetes presentes sur les deux cca meme si les activites ne sont pas dans la solution
    def staticCommonRequests(self,cca1,cca2):
        couple = tuple(sorted((cca1,cca2)))
        if couple in self.requetes_en_commun:
            return self.requetes_en_commun[couple]
        else:
            r1 = self.getPresentRequestsStatic(cca1)
            r2 = self.getPresentRequestsStatic(cca2)
            if len(r1)<len(r2):
                res = [r for r in r1 if isInSortedList(r2,r)]
            else:
                res = [r for r in r2 if isInSortedList(r1,r)]
            self.requetes_en_commun[couple] = res
            return res
    
    def getPresentRequestsDynamic(self,constellation,id_cca):
        req = []
        s,cca = id_cca
        for a in self.getSolCCA(s,cca).getSequence():
            r = constellation.getRequestActivity(a)
            if not isInSortedList(req,r):
                bisect.insort_right(req,r)
        return req
                    
    # calcule les requêtes qui peuvent être déplacées d'une cca à l'autre
    # dynaMique(cca1,cca2) = requetes présente dans la séquence de l'une avec des opportunites dans l'autre
    def commonRequestsDynamic(self,constellation,cca1,cca2):
        s1 = cca1[0]
        s2 = cca2[0]
        req_en_commun = []
        req_1 = []
        req_2 = []
        static1 = self.getPresentRequestsStatic(cca1)
        static2 = self.getPresentRequestsStatic(cca2)
        for a in self.getSolCCAs(s1,cca1).getSequence():
            r = constellation.getRequestActivity(a)
            if not isInSortedList(req_1,r):
                bisect.insort_right(req_1,r)
                if not isInSortedList(req_en_commun,r) and r in static2:
                    bisect.insort_right(req_en_commun,r)
        for a in self.getSolCCAs(s2,cca2).getSequence():
            r = constellation.getRequestActivity(a)
            if not isInSortedList(req_2,r):
                bisect.insort_right(req_2,r)
                if not isInSortedList(req_en_commun,r) and r in static1:
                    bisect.insort_right(req_en_commun,r)
        return req_en_commun,req_1,req_2
    
    def getModeIfPresent(self,r):
        for (rr,m) in self.solution.selectedModes:
            if rr==r:
                return m
        return None
    
    def addComment(self,commentaire):
        self.solution.history.addComment(commentaire)
        
    def updateHistory(self,event,constellation):
        self.solution.updateObjective(constellation)
        mesure = time()-self.startDate
        self.solution.history.updateHistory(mesure,event,self.solution.getObjective(),self.solution.getSolCCAs(),self.solution.getSelectedModes(),self.dependenciesGraph,constellation)
     
    def notifyEndOperator(self,constellation):
        if config.verifMode():
            printMaster("end operator",time()-self.startDate,force=True)
        
        self.updateHistory(END_OPERATOR,constellation)        
     
    def notifyPreprocessing(self,constellation):    
        self.updateHistory(PREPROCESSING, constellation)
        
    def notifyEndIteration(self,constellation):
        self.updateHistory(END_ITERATION, constellation)

    def notifyNoEvent(self,constellation):
        self.updateHistory(NO_EVENT,constellation)  

    def notifyEndExecution(self,constellation):
        self.updateHistory(END_RUN, constellation) 

    def notifyBackupSolution(self,constellation):
        self.updateHistory(BACKUP, constellation)              
    
    def deleteLastPoint(self):
        self.solution.history.deleteLastPoint()   
    """
    def gapRecompense(self,constellation):
        up = self.majorantRecompense(constellation)[0]
        return (up - self.objective[0])/up<=config.glob.tolerance_opti
            
    def gapTemps(self,constellation):
        up2 = self.majorantRecompense(constellation)[1]
        return (up2 - self.objective[1])/up2<=config.glob.tolerance_temps
    
    def majorantRecompense(self,constellation):
        return self.solution.history.majorantRecompense(constellation)
    """
    def initRequests(self,constellation,noise=0,initFirstMode=True,conserveOld=True):
        if not conserveOld:
            self.candidateRequests = list([r for r in constellation.getRequests()])
        else:
            self.candidateRequests = list([r for r in constellation.getRequests() if r not in self.fulfilledRequests()])
        for r in constellation.getAllRequests():
            constellation.getRequest(r).init = False
            if initFirstMode:
                mode = constellation.getRequest(r).getCurrentMode(constellation)
                if mode is None:
                    raise ValueError("Mode requête",r,None)
                m = mode.getId()
                self.currentModes[r] = m
            else:
                constellation.getRequest(r).resetModes(constellation,initFirstMode=False,conserveOld=conserveOld)
        # Trier les req candidates
        if initFirstMode:
            rec_courante = lambda r:constellation.getRequest(r).getMode(self.currentModes[r]).getUtility()
            prio = lambda r:constellation.getRequest(r).getPriority()
            cle = lambda r : (prio(r),rec_courante(r))
            self.candidateRequests = sorted(self.candidateRequests,key=cle,reverse=True)
        
        # score candidats modes avant la planif
        modes = []
        for r in self.candidateRequests:
            modes.append((r,0))
        self.resetNoise(constellation,noise)
        self.solution.history.meanScoreBefore = constellation.meanObservationScore(modes)
        self.solution.history.meanObservationScore = self.solution.history.meanScoreBefore
        
        
    def resetNoise(self,constellation,amplitude):
        self.requestNoise = {}
        for r in constellation.getAllRequests():
            self.requestNoise[r] = (self.randomizerPerturbation.random()-0.5)*amplitude
            
    def resetScoreDestruction(self,constellation,fulfilledRequests):
        unfulfilledRequests = []
        i = 0
        for r in constellation.getAllRequests():
            if r not in fulfilledRequests:
                unfulfilledRequests.append((i,r))
                i += 1
        weight = []
        for i,r in unfulfilledRequests:
            weight.append(constellation.getRequest(r).getCurrentMode(constellation).getUtility())
        rmax = len(unfulfilledRequests)
        for i in range(rmax):
            (i,selected_request) = rd.choices(unfulfilledRequests,k=1,weights=weight)[0]
            self.requestNoise[r] = -weight[i] + rmax-i
            weight[i] = 0
    
    def writeSol(self,instance):
        filename = config.donnees.algo_name+'-'+instance
        with open(filename,'w') as file:
            for s in sorted(self.getSolCCAs()):
                for cca in sorted(self.getSolCCAs()[s]):
                    file.write("CCA"+str(cca)+str(sorted(self.getSolCCA(s,cca).sequence))+'\n')
        
    def sortModes(self,modes): # modes = [(r,m,w)]
        cle = lambda rm : (rm[0]==-1,rm[2] + self.requestNoise[rm[0]])
        return sorted(modes,key=cle,reverse=True)
        
    def saveSample(self,constellation,add=None):
        self.solution.history.saveSample(constellation,add=add)
    
    def plotObjective(self):
        self.solution.history.plotObjective()
        
    def plotCPU(self):
        self.solution.history.plotCPU()
        
    def plotCCAsLoad(self):   
        self.solution.history.plotCCAsLoad(self.dependenciesGraph)
        
    def createMappingCCASlaves(self):
        if MPI.COMM_WORLD.Get_size()>1:
            self.ccaSlaves = {cpu : [] for cpu in range(1,MPI.COMM_WORLD.Get_size())}
            sizes = []
            for cca in self.dependenciesGraph.getComponents():
                size =  self.dependenciesGraph.getComponentSize(cca)
                reverse_insort(sizes,(taille,cca))
            cpu = 0
            for (size,cca) in sizes:
                self.ccaSlaves[cpu+1].append(cca)
                cpu = (cpu + 1) % (MPI.COMM_WORLD.Get_size()-1)
            for cpu in self.ccaSlaves:
                MPI.COMM_WORLD.send({"cca":self.ccaSlaves[cpu]},dest=cpu)
        else:
            self.ccaSlaves = {1 : self.dependenciesGraph.getComponents()}
     
    def mappingCCASize(self,sizes):
        # sizes : liste (cca,nb activities a inserer)
        self.ccaSlaves = {cpu : [] for cpu in range(1,MPI.COMM_WORLD.Get_size())}
        ccas = sorted(sizes,key=itemgetter(0),reverse=True)
        cpu = 0
        for (size,cca) in ccas:
            self.ccaSlaves[cpu+1].append(cca)
            cpu = (cpu + 1) % (MPI.COMM_WORLD.Get_size()-1)
    
    def restart(self,constellation):
        self.solution.selectedModes = []
        for (s,cca) in self.dependenciesGraph.getComponents():
                self.setSolCCA(s,cca,SolCCA(cca,s))
        
        
    """
        =============================================== 
                        GETTER
        =============================================== 
    """
    
    def getSelectedModes(self):
        return self.solution.getSelectedModes()
    
    def addSelectedMode(self,mode,constellation):
        self.solution.addSelectedMode(mode,constellation)
        
    def removeSelectedMode(self,mode,constellation):
        self.solution.removeSelectedMode(mode,constellation)
        
    def getBest(self):
        return self.solution.history.getBest()
    
    def restartBestSolution(self,constellation):
        bestObjective,bestSolution,bestSelectedModes = self.solution.history.getBest()
        self.solution.setSolCCAs(bestSolution)
        self.solution.setSelectedModes(bestSelectedModes, constellation)
    
    def getSolutionObservations(self):
        sol = {s : [] for s in self.getSolCCAs()}
        for s in sol:
            ccas = sorted([self.getSolCCA(s,cca) for cca in self.getSolCCAs()[s]],key=lambda solcca : solcca.getStartDate())
            for cca in ccas:
                for a in cca.getSequence():
                    sol[s].append(a)
        return sol   
    
    def extractSequences(self):
        seq = {}
        for s in self.getSolCCAs():
            if s not in seq:
                seq[s] = []
            ccas = sorted([self.getSolCCA(s,cca) for cca in self.getSolCCAs()[s]],key=lambda solcca : solcca.getStartDate())
            for cca in ccas:
                seq[s].append(cca.getSequence())
        return seq
    
    def getSolutionContainer(self):
        return self.solution
    
    def getSolution(self):
        sol = {s : [] for s in self.getSolCCAs()}
        for s in sol:
            ccas = sorted([self.getSolCCA(s,cca) for cca in self.getSolCCAs()[s]],key=lambda solcca : solcca.getStartDate())
            for cca in ccas:
                for a in cca.getSequence():
                    sol[s].append(a)
        return sol            
    
    def getModes(self):
        return self.solution.selectedModes
    
    def getPlan(self):
        try:
            return self.plan.copy()
        except Exception:
            return None
        
    def fulfilledRequests(self):
        return [x[0] for x in self.solution.getSelectedModes()]
    
    def removeRequest(self,r):
        del self.candidateRequests[r]
    
    def removeFulfilledRequest(self,request):
        indice = -1
        for i,(r,m) in enumerate(self.solution.selectedModes):
            if r == request:
                indice = i
                mode = (r,m)
                break
        if indice == -1:
            raise ValueError("Requête non présente")
        else:
            self.solution.selectedModes.pop(indice)
            return mode
            
    def replaceFulfilledRequest(self,request):
        indice = -1
        for i,(r,m) in enumerate(self.solution.selectedModes):
            if r == request:
                indice = i
                break
        if indice == -1:
            self.solution.selectedModes.append((r,m))
        else:
            self.solution.selectedModes[indice] = (r,m)
    
    def setObjective(self,objective):
        self.solution.objective = objective
    
    def setSelectedModes(self,selectedModes,constellation):
        self.solution.setSelectedModes(selectedModes,constellation)
            
    """
        Version qui ne recalcule l'objective. N'utiliser que si necessaire
    """
    def overwriteModes(self,selectedModes):
        self.solution.selectedModes = copy(selectedModes)  
    
    def getObjective(self,constellation=None,recompute=False):
        if recompute:
            assert(constellation is not None)
            return sum([constellation.getRequest(r).getMode(m).getUtility() for (r,m) in self.solution.selectedModes])
        else:
            return self.solution.objective
        
    def getSolCCAs(self):
        return self.solution.solCCAs
    
    def getSolCCA(self,s,cca):
        return self.solution.solCCAs[s][cca]
    
    def setSolCCA(self,s,cca,solCCA):
        self.solution.solCCAs[s][cca] = solCCA
        
    def setAllSolCCAs(self,solCCAs):
        self.solution.solCCAs = solCCAs
        
    """
        =============================================== 
                    CONSTRUCTION DU PLAN
        =============================================== 
    """
    def planEarliestSequence(self,constellation,solution,s,transitionModel=None):
        plan = []
        for i,activity in enumerate(solution):
            if(i==0):
                t = constellation.getSatellite(s).getActivity(activity).getStartDate()
            else:
                prec = solution[i-1]
                duration = constellation.getSatellite(s).getActivity(prec).getDuration()
                start = constellation.getSatellite(s).getActivity(activity).getStartDate()
                if transitionModel is None:
                    transition = self.getTransitionDuration(constellation,s,prec,activity,t+duration,self.transitionModel)
                else:
                    transition = self.getTransitionDuration(constellation,s,prec,activity,t+duration,transitionModel)
                t = max(t + duration + transition,start)
            plan.append((activity,t))
        return plan
    
    def planEarliest(self,constellation,transitionModel=None):
        solution = self.getSolution()
        self.plan = {s : [] for s in solution}
        for s in solution:
            for i,activity in enumerate(solution[s]):
                if(i==0):
                    t = constellation.getSatellite(s).getActivity(activity).getStartDate()
                else:
                    prec = solution[s][i-1]
                    duration = constellation.getSatellite(s).getActivity(prec).getDuration()
                    start = constellation.getSatellite(s).getActivity(activity).getStartDate()
                    if transitionModel is None:
                        transition = self.getTransitionDuration(constellation,s,prec,activity,t+duration,self.transitionModel)
                    else:
                        transition = self.getTransitionDuration(constellation,s,prec,activity,t+duration,transitionModel)
                    t = max(t + duration + transition,start)
                self.plan[s].append((activity,t))
                
    def planLatest(self,constellation):
        assert(not self.transitionModel.isTimeDependent())
        solution = self.getSolution()
        self.plan = {s : [] for s in solution}
        for s in solution:
            for i in range(len(solution[s])-1,-1,-1):
                a = solution[s][i]
                duration = constellation.getSatellite(s).getActivity(a).getDuration()
                if(i==len(solution[s])-1):
                    t = constellation.getSatellite(s).getActivity(activity).getEndDate() - duration
                else:
                    suivante = solution[s][i+1]
                    transition = self.getTransitionDuration(constellation,s,a,suivante,t+duration,self.transitionModel)
                    die("TO DO")
                    t = min(t - duration - transition,constellation.getSatellite(s).getActivity(activity).getEndDate())
                self.plan[s].append((a,t))
                
    # deduit le temps le plus tot pour chaque activite
    def buildPlan(self,constellation,transitionModel=None):
        self.planEarliest(constellation,transitionModel=transitionModel)
        if config.verifMode():
            self.verifyFeasability(constellation,transitionModel=transitionModel)
        
    def totalTardiness(self):
        pi = 0
        solution = self.getSolution()
        for s in solution:
            for i,activity in enumerate(solution[s]):
                if(i==0):
                    t = constellation.getSatellite(s).getActivity(activity).getStartDate()
                else:
                    prec = solution[s][i-1]
                    duration = constellation.getSatellite(s).getActivity(prec).getDuration()
                    start = constellation.getSatellite(s).getActivity(activity).getStartDate()
                    if transitionModel is None:
                        transition = self.getTransitionDuration(constellation,s,prec,activity,t+duration,self.transitionModel)
                    else:
                        transition = self.getTransitionDuration(constellation,s,prec,activity,t+duration,transitionModel)
                    t = max(t + duration + transition,start)
                pi_obs = max(0,t+duration-self.constellation.getSatellite(s).getActivity(a).getEndDate())
                pi += pi_obs
        return pi
    
    # renvoie un dict : satellites -> liste de (retard,flexibilite,activite) triée par retard décroissant
    def cumulatedTardinessLeftToRight(self,constellation,transitionModel=None):
        tardinesses = {}
        solution = self.getSolution()
        for s in solution:
            retard[s] = []
            for i,activity in enumerate(solution[s]):
                start = constellation.getSatellite(s).getActivity(activity).getStartDate()
                duration = constellation.getSatellite(s).getActivity(prec).getDuration()
                if(i==0):
                    t = constellation.getSatellite(s).getActivity(activity).getStartDate()
                    durationPredecessor = 0
                    transition = 0                
                else:
                    prec = solution[s][i-1]
                    if transitionModel is None:
                        transition = self.getTransitionDuration(constellation,s,prec,activity,t+duration,self.transitionModel)
                    else:
                        transition = self.getTransitionDuration(constellation,s,prec,activity,t+duration,transitionModel)
                    t = max(t + duration + transition,start)
                right = max(0,t-start)
                left = max(0,t+duration-constellation.getSatellite(s).getActivity(activity).getEndDate())
                reverse_insort(tardinesses[s],(right,left,activity))
        return tardinesses 
    
    def isTemporallyFeasible(self,constellation):
        self.planEarliest()
        for s in self.plan:
            for i,(p,t) in enumerate(self.plan[s]):
                if(i<len(self.plan[s])-1):
                    if(constellation.getSatellite(s).estObservation(p)):
                        duration = self.constellation.getSatellite(s).getObservation(p).getDuration()
                    else:
                        duration = self.constellation.getSatellite(s).getDownload(p).getDuration()
                    transition = self.getTransitionDuration(constellation,s,p,self.plan[s][i+1][0],t+duree,self.transitionModel)
                    if(t+transition+duration >self.plan[s][i+1][1] + 1e-5):
                        return False
                    if(t+duration>self.constellation.getSatellite(s).getActivity(p).getEndDate()):
                        return False
        return True

    def planEarliestSequence(self,constellation,solution,s,transitionModel):
        plan = []
        for i,activity in enumerate(solution):
            if(i==0):
                t = constellation.getSatellite(s).getActivity(activity).getStartDate()
            else:
                predecessor = solution[i-1]
                duration = constellation.getSatellite(s).getActivity(predecessor).getDuration()
                start = constellation.getSatellite(s).getActivity(activity).getStartDate()
                if transitionModel is None:
                    transition = self.getTransitionDuration(constellation,s,predecessor,activity,t+duration,self.transitionModel)
                else:
                    transition = self.getTransitionDuration(constellation,s,predecessor,activity,t+duration,transitionModel)
                t = max(t + duration + transition,start)
            plan.append((activity,t))
        return plan
    
    def getTransitionDuration(self,constellation,s,a1,a2,start,transitionModel):
        if transitionModel.isTimeDependent():
            transition = constellation.getSatellite(s).getTransitionTimeDependent(a1,a2,start,transitionModel)
        else:
            transition = constellation.getSatellite(s).getTransitionDuration(a1,a2,transitionModel)                    
        return transition
        
    def isTemporallyFeasibleSequence(self,constellation,sequence,s,transitionModel=None):
        for i,activity in enumerate(sequence):
            if(i==0):
                t = constellation.getSatellite(s).getActivity(activity).getStartDate()
            else:
                predecessor = sequence[i-1]
                duration = constellation.getSatellite(s).getActivity(predecessor).getDuration()
                start = constellation.getSatellite(s).getActivity(activity).getStartDate()
                if transitionModel is None:
                    transition = self.getTransitionDuration(constellation,s,predecessor,activity,t+duration,self.transitionModel)
                else:
                    transition = self.getTransitionDuration(constellation,s,predecessor,activity,t+duration,transitionModel)
                t = max(t + duration + transition,start)
                if t+constellation.getSatellite(s).getActivity(activity).getDuration() > constellation.getSatellite(s).getActivity(activity).getEndDate():
                    return False
        return True        
    
    """
        =============================================== 
                        RESOLUTION
        =============================================== 
    """
    def getGraphComponents(self):
        return self.dependenciesGraph
    
    """
        afficher les informations de la solution courante.
        filtres possibles : time,objective,modes,ratio,best,repartition,obs,request,size
        title : indiquer un titre
        temps : date courante
        startDate : date de début d'execution
        constellation : la constellation
        add = dictionnaire de messages optionnels (titre,contenu)
    
    """
    def displayInformation(self,temps,startDate,constellation,core=None,color='c',add={},title='Info',filtre=None):
        if core is None or MPI.COMM_WORLD.Get_rank()==core:
            try:
                self.setObjective(self.computeObjective(constellation))
            except Exception as e:
                pass
            shift = getDisplayDepth()
            shiftLeftDisplay(shift-1)
            self.filtre = filtre
            self.solution.history.filtre = filtre
            title_size_max = 30
            title_string = title[0:max(len(title),title_size_max)]
            title_before = (title_size_max - len(title_string))//2
            title_after = title_size_max - len(title_string) - title_before
            center = 4 + title_size_max 
            width = 120
            left = (width-center)//2
            right = width - left - center
            printColor("\n")
            title_msg = "[" + " "*title_before + title + title_after*" " + " ]"
            left_msg = " " + (left-2)*"-" + " " 
            right_msg = " " + (right-1)*"-"
            printColor( left_msg + title_msg +  right_msg ,c=color)
            for text in add:
                printColor("| "+str(text)+": " + str(add[text]),c=color)
            if filtre is None or 'time' in filtre:
                printColor("| time elapsed: ",(temps-startDate)//60,"minutes",(temps-startDate)%60,c=color)
            if filtre is None or 'obj' in filtre:
                printColor("| objective: ",self.getObjective(),c=color)
            if filtre is None or "modes" in filtre:
                printColor("| selected modes: ",len(self.solution.selectedModes),c=color)
            lines = str(self.solution.history).split('\n')
            for line in lines:
                if len(line)>0:
                    printColor(line,c=color)
            printColor(' ' + (width-1)*"-",c=color)
            shiftRightDisplay(shift-1)
            
    def insertSequences(self,constellation,sequences):
        selectedRequests = [x[0] for x in self.solution.selectedModes]
        for s,cca in sequences:
            seq = [a for a in sequences[cca] if constellation.getRequestActivity(a) in selectedRequests]
            #print("set sequence",s,cca)
            self.getSolCCA(s,cca).setSequence(constellation,seq)
    
    def intersectionActivities(self,constellation,r,composantes):
        newMode = self.currentModes.get(r)
        oldMode = newMode - 1
        act_ancien = [x[1] for x in constellation.getRequest(r).getMode(oldMode).getPairs()]
        if newMode is None:
            return act_ancien,[],[]
        act_newMode = [x[1] for x in constellation.getRequest(r).getMode(newMode).getPairs()]
        removals = []
        prevalidations = {}
        intersectionWorking = {}
        for o in act_newMode:
            cca = composantes.find(o)
            if o in act_ancien and o in self.etat_recherche.activities_done[(r,oldMode)]:
                if cca not in prevalidations:
                    prevalidations[cca] = []
                prevalidations[cca].append(o)
            elif o in self.etat_recherche.activities_working[(r,oldMode)]:
                if cca not in intersectionWorking:
                    intersectionWorking[cca] = []
                intersectionWorking[cca].append(o)
        for o in act_ancien:
            if o not in act_newMode:
                removals.append(o)
        
        return removals,prevalidations,intersectionWorking

    def terminateProcess(self):
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        size = comm.Get_size()
        for i in range(size):
            if i!=0:
                comm.send({'fin':True},dest=i)
        comm.Barrier()
    
    def getModeActivities(self,constellation,r,m):
        activities = []
        for s in self.constellation.getRequest(r).getMode(m).getActivities():
            for o in self.constellation.getRequest(r).getMode(m).getActivities()[s]:
                if o not in activities:
                    activities.append(o)
        return activities
    
    def getAvailableProcesses(self,freeCPUs):
        return freeCPUs!=[]
    
    def choseCPU(self,freeCPUs):
        return freeCPUs.pop()
    
    def noiseUtility(self,constellation,r,m):
        return (constellation.getRequest(r).getMode(m).getUtility() + self.requestNoise[r],constellation.getRequest(r).getMode(m).getTemporalUtility())
    
    def mappingCCACurrentObservations(self,constellation,failed):
        mappingObservations = {}
        for r in self.candidateRequests:
            m = constellation.getRequest(r).getCurrentMode().getId()
            w = constellation.getRequest(r).getCurrentMode().getUtility()
            for (s,o) in constellation.getRequest(r).getMode(m).getPairs():
                cca = self.dependenciesGraph.getActivityCCA(o)
                if cca not in mappingObservations:
                    mappingObservations[cca] = {}
                if (r,m,w) not in mappingObservations[cca]:
                    mappingObservations[cca][(r,m,w)] = []
                mappingObservations[cca][(r,m,w)].append(o)
        
        mappingObservations,sequencesCCA = self.filtrer(constellation,mappingObservations,failed)    
        return mappingObservations,sequencesCCA
        
    # si removals est None on retire toutes les activities du mode. Sinon on retire la liste 'removals'
    def removeModeFromSolution(self,constellation,mode,removals=None):
        r,m = mode
        cca = {}                     
        for (s,o) in constellation.getRequest(r).getMode(m).getPairs():
            if removals is None or o in removals:
                cca_o = self.dependenciesGraph.getActivityCCA(o)
                if (s,cca_o) not in cca:
                    cca[(s,cca_o)] = []
                cca[(s,cca_o)].append(o)
        for (s,cca_a) in cca:
            self.getSolCCA(s,cca_a).retirerListeActivites(constellation,cca[(s,cca_a),self.transitionModel])
                              
    def deleteRequest(self,constellation,r):
        if r in self.currentModes:
            del self.currentModes[r]
        if r in self.candidateRequests:
            self.candidateRequests.remove(r)
            for cca in self.ccaRequetes[r]:
                self.retraitRequeteCCA(cca,r)
            del self.ccaRequetes[r]

        
    """
        =============================================== 
                        GRAPHIQUES
        =============================================== 
    """
    """
    def plotActivity(self,constellation,annotate=False):
        f,ax = plt.subplots(figsize=(15,6))
        colorObs = 'c'
        colorVid = 'g'
        colorTrans = 'r'
        
        for s in self.plan:
            #horizon obs
            for p in constellation.getSatellite(s).getObservations():
                startDate = constellation.getSatellite(s).getObservation(p).getStartDate()
                endDate = constellation.getSatellite(s).getObservation(p).getEndDate()
                ax.plot([startDate,endDate],[s,s],'-'+colorObs,alpha=1,linewidth=2)
            for p in constellation.getSatellite(s).getDownloads():
                startDate = constellation.getSatellite(s).getDownload(p).getStartDate()
                endDate = constellation.getSatellite(s).getDownload(p).getEndDate()
                ax.plot([startDate,endDate],[s,s],'-'+colorVid,alpha=1,linewidth=2)
            
            for i,(p,t) in enumerate(self.plan[s]):
                if(constellation.getSatellite(s).isDownload(p)):
                    color = colorVid
                    duration = constellation.getSatellite(s).getDownload(p).getDuration()
                    lab = 'vidage'
                else:
                    color = colorObs
                    duration = constellation.getSatellite(s).getObservations()[p].getDuration()
                    lab = 'observation'
                ax.add_patch(patches.Rectangle((t,s-1/2),duration,1,color=color,fill=True,label=lab))
                
                if i<len(self.plan[s])-1:
                    transition = self.getTransitionDuration(constellation,s,p,self.plan[s][i+1][0],t+duration,self.transitionModel)
                    ax.add_patch(patches.Rectangle((t+duration,s-1/2),transition,1,color=colorTrans,fill=True,label='transition'))
                if(annotate):
                    ax.annotate('{}'.format(self.sol[s][i][0]),xy=(self.sol[s][i][1], s),xytext=(0, 3),textcoords="offset points",ha='center', va='bottom',color='k')           
        
        plt.yticks(sorted(list(constellation.satellites.keys())))
        if(annotate):
            plt.xlabel("time (s)")
            plt.ylabel("satellite id")
        plt.grid(axis='y')
        ax.xaxis.set_tick_params(labelsize=15)
        ax.yaxis.set_tick_params(labelsize=15)
        a,b,modes = self.solution.history.getBest()
        file = config.getOptValue("instance")["file"]
        folder = config.getOptValue("instance")["folder"]
        algo = config.donnees.algo_name
        #plt.title(config.glob.filename +  " - " +str(config.glob.size) + ' processes - Plan : '+str(len(modes))+" requests")        
        plt.savefig('results/plan/plan_'+folder+"_"+file+"_"+algo+"_"+str(config.glob.size)+".png")
        return f
            
    def plotModesLosses(self):
        rangs = {}
        perte = []
        f,(ax1,ax2) = plt.subplots(1,2)
        for r,m in self.solution.selectedModes:
            tri = sorted([self.constellation.getRequest(r).getMode(mm).getUtility() for mm in self.constellation.getRequest(r).getModes()],reverse=True)
            best = tri[0]
            for i,rec in enumerate(tri):
                if(rec==self.constellation.getRequest(r).getMode(m).getUtility()):
                    if i not in rangs:
                        rangs[i] = 1
                    else:
                        rangs[i] = rangs[i]+1
                    perte.append(100*(best-rec)/best)
                    break
        ax1.set(xlabel='rank of selected modes', ylabel='count')
        ax1.set_title("rank of selected modes - "+str(len(self.solution.selectedModes))+"/"+str(len(self.constellation.getAllRequests())) + " covered")
        ax2.set(xlabel='reward loss compared to the best mode (%)', ylabel='count')
        ax2.set_title("reward loss of selected modes")
           
        n1,bins1,p1 = ax1.hist(rangs,bins=np.arange(0,max(rangs)))
        n2,bins2,p2 = ax2.hist(perte)
        return f
    
    """
    """
        =============================================== 
                        VERIF/DEBUG
        =============================================== 
    """
    def verifyLists(self):
        selectedRequests = [x[0] for x in self.solution.selectedModes]
        assert(len(np.unique(selectedRequests))==len(selectedRequests))
        for request in self.candidateRequests:
            assert(request not in [x[0] for x in self.solution.selectedModes])
        for request in [x[0] for x in self.solution.selectedModes]:
            assert(request not in [x[0] for x in self.candidateModes])
        for (r,m) in self.solution.selectedModes:
            assert((r,m) not in self.candidateModes)
            assert((r,m) not in [x[0] for x in self.removedModes[r]])
        for r in self.candidateRequests:
            for m in self.getRequest(r).candidateModes:
                assert((r,m) not in [x[0] for x in self.removedModes[r]])
        
        for r in self.removedModes:
            for ((rr,m),i) in self.removedModes[r]:
                if(r in self.candidateRequests):
                    for m in self.getRequest(r).candidateModes:
                        print((r,m))
                assert(m not in self.candidateModes[r])
    
    def verifyMode(self,constellation,r,m):
        for s,o in constellation.getRequest(r).getMode(m).getPairs():
            _,cca = self.dependenciesGraph.getActivityCCA(o)
            if not o in self.getSolCCA(s,cca).getSequence():
                print(self.getSolCCA(s,cca),o)
                comm = MPI.COMM_WORLD
                rank = comm.Get_rank()
                coeur = "process: "+str(rank)
                raise ValueError(coeur,"time: ",time()-self.startDate,"missing activity",s,o,'mode',(r,m))
        if r != -1:
            if not constellation.getRequest(r).getMode(m).acceptable(constellation):
                raise ValueError(coeur,"time: ",time()-self.startDate,"content of the mode is unacceptable (inconsistant intialization?)",constellation.getRequest(r).getMode(m))
    
    def verifyDownloadModeSelection(self,constellation):
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        coeur = "process: "+str(rank)        
        for (s,d) in constellation.getRequest(-1).getPairs():
            _,cca = self.dependenciesGraph.getActivityCCA(d)
            if d not in self.getSolCCA(s,cca).getSequence():
                raise ValueError(coeur,"time: ",time()-self.startDate,"Missing download slot",d,self.getSolCCA(s,cca).getSequence())
        
    def verifyModes(self,constellation):
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        coeur = "process: "+str(rank)         
        self.verifyDownloadModeSelection(constellation)
        for (r,m) in self.solution.selectedModes:
            self.verifyMode(constellation,r,m)
        for s in self.getSolCCAs():
            for cca in self.getSolCCAs()[s]:
                for o in self.getSolCCA(s,cca).getSequence():
                    find = False
                    for (r,m) in self.solution.selectedModes:
                        if (s,o) in constellation.getRequest(r).getMode(m).getPairs():
                            find = True
                            break
                    if not find:
                        r = constellation.getRequestActivity(o)
                        typ = constellation.getRequest(r).getType()
                        print("current mode:",constellation.getRequest(r).getCurrentMode(constellation))
                        print("modes of the request")
                        for m in constellation.getRequest(r).modes:
                            print(constellation.getRequest(r).getMode(m))
                        raise ValueError(coeur,"time: ",time()-self.startDate,'activity with no mode',s,o,cca,typ,r,".Fulfilled requests:",r in self.fulfilledRequests())
    
    def verifyCCAs(self,constellation):
        for s in self.getSolCCAs():
            for cca in self.getSolCCAs()[s]:
                faisable = self.getSolCCA(s,cca).isSequenceFeasible(constellation,self.transitionModel)
                if not faisable:
                    die('infeasible CCA',self.getSolCCA(s,cca),"tardiness=",self.getSolCCA(s,cca).retardSequence(constellation,self.getSolCCA(s,cca).sequence,self.transitionModel))
                
    def verifyModeUniqueness(self):
        unique = []
        for (r,m) in self.solution.selectedModes:
            if(r in unique):
                raise ValueError('request covered by several modes')
            else:
                unique.append(r)
    
    def verifyFeasability(self,constellation,transitionModel=None):
        self.verifyCCAs(constellation)
        for s in self.plan:
            transitions = []
            if transitionModel is None:
                transitionModel = self.transitionModel
            for i,(p,t) in enumerate(self.plan[s]):
                if(i<len(self.plan[s])-1):
                    duration = constellation.getSatellite(s).getActivity(p).getDuration()
                    transition = self.getTransitionDuration(constellation,s,p,self.plan[s][i+1][0],t+duration,transitionModel)
                    transitions.append(transition)
                    if(t+transition+duration >self.plan[s][i+1][1] + 1e-5):
                        raise ValueError('transition finishing after next observation start date')
                    if(t+duration>constellation.getSatellite(s).getActivity(p).getEndDate()):
                        print("tardiness:",t+duration-constellation.getSatellite(s).getActivity(p).getEndDate())
                        startDate = constellation.getSatellite(s).getActivity(p).getStartDate()
                        endDate = constellation.getSatellite(s).getActivity(p).getEndDate()
                        printColor("activity", (s,p),"date",t,'duration',duration,"window",startDate,endDate,c='r')
                        print([(p,self.dependenciesGraph.getActivityCCA(p),t) for (p,t) in self.plan[s]])
                        raise ValueError('window tardiness:CCAs solutions are maybe disordered')

    def verifyActiveRequests(self,constellation):
        if config.getOptValue("dynamic"):
            for (r,m) in self.solution.selectedModes:
                if not constellation.getRequest(r).estActif():
                    raise ValueError("Requête inactive")
        
    def verifySolutionIfVerifyMode(self,constellation,transitionModel=None):
        if config.getOptValue("verif"):
            self.verifySolution(constellation,transitionModel=transitionModel)
        
    def verifySolution(self,constellation,transitionModel=None):
        self.buildPlan(constellation,transitionModel=transitionModel)
        self.verifyActiveRequests(constellation)
        self.verifyModes(constellation)
        #self.verifyModeUniqueness()

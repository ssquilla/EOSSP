#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 10:31:30 2022

@author: ssquilla
"""

from EOSSP.model.solution import *
from EOSSP.model.constellation import *
from EOSSP.model.components import *
from EOSSP.model.componentPlan import *
from EOSSP.Utils.Utils import *
from EOSSP.Utils.Utils import printColor,printOpen,printClose,alert,warn
from EOSSP.Utils.Utils import shiftRightDisplay,shiftLeftDisplay
from EOSSP.Utils.config import Config

global config
config = Config()

import random as rd
import math
import docplex.cp as cp
from docplex.cp.model import CpoModel,interval_var,binary_var


class Operator:
    def __init__(self,name,transitionModel,k=2):
        self.k=k
        self.transitionModel = transitionModel
        self.name = name
        self.tabuList = {} # pairs de cca tabu
        self.observationsCCA = None
        self.connectCCAs = None
        mpiId = MPI.COMM_WORLD.Get_rank()
        configuredSeed = config.getOptValue("seed")
        self.randomInstanceCCAChoice = rd.Random(mpiId+configuredSeed)
        
    def acceptable(self,paircca):
        sortedPairs = tuple(sorted(paircca))
        return sortedPairs not in self.tabuList
    def getName(self):
            return self.name
    def tabuUpdate(self):
        pairsToDelete = []
        for pair in self.tabuList:
            self.tabuList[pair] -= 1
            if self.tabuList[pair]==0:
                pairsToDelete.append(pair)
        for pair in pairsToDelete:
            del self.tabuList[pair]
    def __str__(self):
        return 'Operator '+self.name
    def _initializeCCAData(self,LNS,constellation):
        self.observationsCCA = LNS.getCCAObservations()
        self.computeConnectivity(LNS)
        
    def checkEmptyCCA(self,LNS,newActivities,CCAs):
         CCAsObservations = []
         for r in newActivities:
             for m in newActivities[r]:
                 for a in newActivities[r][m]:
                     idCCA = LNS.getGraphComponents().getActivityCCA(a)
                     if idCCA not in CCAsObservations:
                         assert(idCCA in CCAs)
                         CCAsObservations.append(idCCA)
                         if len(CCAs)==len(CCAsObservations):
                             assert(sorted(CCAs)==sorted(CCAsObservations))
                             return False
         return len(CCAs)!=len(CCAsObservations)
    
    def computeConnectivity(self,LNS):
        self.connectCCAs = {}
        for cca1 in LNS.getGraphComponents().getComponents():
            for cca2 in LNS.getGraphComponents().getComponents():
                if cca1<cca2:
                    if cca1 not in self.connectCCAs:
                        self.connectCCAs[cca1] = {}
                    self.connectCCAs[cca1][cca2] = LNS.staticCommonRequests(cca1,cca2)
    
    def updateNewRequests(self,LNS):
        if config.getOptValue("n_cca")>1:
            self.computeConnectivity(LNS)
        self.observationsCCA = LNS.getCCAObservations()
        
    def getConnectivities(self,cca1,cca2):
        c1,c2 = sorted((cca1,cca2))
        return self.connectCCAs[c1][c2]
    
    def _choseFirstCCA(self,LNS,constellation):
        # choix 1ere CCA
        if config.getOptValue("nbh")=="load":
            CCAMeasures = [len(LNS.getPresentRequestsDynamic(constellation,idCCA)) for idCCA in self.observationsCCA]
        elif config.getOptValue("nbh")=="random":
            CCAMeasures = [1 for idCCA in self.observationsCCA]
        else:
            raise ValueError("unknown CP neighborhood CP inconnu : "+config.getOptValue("nbh"))
        return self.randomInstanceCCAChoice.choices([idCCA for idCCA in self.observationsCCA], weights=CCAMeasures,k=1)[0]    
    
    def _initCCAsResearch(self,LNS,constellation):
        if self.connectCCAs is None:
            self._initializeCCAData(LNS, constellation)
    
    def _choseNextCCA(self,LNS,CCAs,candidats,presentRequests):    
        assert(len(CCAs)>0)
        weights = [sum([len(self.getConnectivities(cca2,idCCA)) for cca2 in CCAs]) for idCCA in candidats]
        if sum(weights)==0:
            return False  # alors on ne trouvera personne
        last = CCAs[-1]
        next_cca = self.randomInstanceCCAChoice.choices(candidats,weights=weights,k=1)[0]
        assert(next_cca not in CCAs)
        CCAs.append(next_cca)
        candidats.remove(next_cca)
        for r in LNS.getPresentRequestsStatic(next_cca):
            if not isInSortedList(presentRequests,r):
                bisect.insort_left(presentRequests, r)
        return True
        
    # Long la 1ere fois : initialisation des connées sur les CCA
    def selectCCAs(self,LNS,constellation):  
       t1 = time()
       self._initCCAsResearch(LNS, constellation)
       t2 = time()
       CCAs = [self._choseFirstCCA(LNS,constellation)]
       t3 = time()
       presentRequests = sorted((LNS.getPresentRequestsStatic(CCAs[0])))
       succes = True
       t4 = time()
       candidats = [idCCA for idCCA in self.observationsCCA if idCCA not in CCAs]
       t5 = time()
       while (succes and len(CCAs)<self.k):
           succes = self._choseNextCCA(LNS,CCAs,candidats,presentRequests)
       t6 = time()
       #print("test",t2-t1,t3-t2,t4-t3,t5-t4,t6-t5)
       return CCAs,presentRequests
   
# Operateur qui ne fait rien.
class IdleOperator(Operator):
    def __init__(self):
        super().__init__("Idle",None)
    def apply(self,LNS,constellation):
        pass

class RepairingOperator(Operator):
        def __init__(self,name,transitionModel,k=2):
            super().__init__(name,transitionModel,k=k)
class GlobalGreedyFilling(Operator):
        def __init__(self,transitionModel):
            super().__init__("greedy constructor",transitionModel)
        
        def apply(self,LNS,constellation):
            oldObjective = LNS.getObjective(constellation)
            LNS.greedyFill(constellation,requestsLimit=config.getOptValue("use_solver"),forbidSolver=False)
            return LNS.getObjective(constellation)>oldObjective
class GlobalGreedyFillingAccelerated(Operator):
        def __init__(self,transitionModel):
            super().__init__("greedy constructor",transitionModel)
        
        def apply(self,LNS,constellation):
            oldObjective = LNS.getObjective(constellation)
            LNS.greedyFill(constellation,requestsLimit=config.getOptValue("use_solver"),forbidSolver=False)
            return LNS.getObjective(constellation)>oldObjective
class DestructionOperator(Operator):
        def __init__(self,name,transitionModel,k=2):
            super().__init__(name,transitionModel,k=k)
        def destroy(self,LNS,constellation):
            raise ValueError("not implemented (abstract class).")
        def removeModes(self,LNS,constellation,modesToRemove):
            lengthBefore= len(LNS.getSelectedModes())
            LNS.setSelectedModes([(r,m) for (r,m) in LNS.getSelectedModes() if (r, m) not in modesToRemove],constellation)
            len_apres = len(LNS.getSelectedModes())
            assert(lengthBefore-len_apres==len(modesToRemove))
            activitiesToRemove = [(s,a) for (r,m) in modesToRemove for (s,a) in constellation.getRequest(r).getMode(m).getPairs()]
            activitiesByCCA = {}
            for (s,a) in activitiesToRemove:
                idCCA = LNS.getGraphComponents().getActivityCCA(a)
                if idCCA not in activitiesByCCA:
                    activitiesByCCA[idCCA] = []
                activitiesByCCA[idCCA].append(a)
            for idCCA in activitiesByCCA:
                s,cca = idCCA
                LNS.getSolCCA(s,cca).removeActivityList(constellation,activitiesByCCA[idCCA],self.transitionModel)            

class Perturbator(DestructionOperator):
    def __init__(self,nom,transitionModel):
        super().__init__(nom,transitionModel)
    def apply(self,LNS,constellation):
        raise ValueError("abstract class")
class DestroyAndRepairGreedyRequest(Perturbator):
        def __init__(self,transitionModel,destroyRate,acceptWorseSolution):
            super().__init__("Destroy and repair greedy",None) # None : pas besoin de stoquer un doublon
            self.acceptWorseSolution = acceptWorseSolution
            self.destroy = RequestDestructor(transitionModel,destroyRate)
            self.filling = GlobalGreedyFilling(transitionModel)
            
        def apply(self,LNS,constellation):
            #printOpen("Perturbation :",c='y')
            oldObjective = LNS.getObjective(constellation)
            if not self.acceptWorseSolution:
                printOpen("Creating a backup of the solution")
                self.backup = LNS.createBackup(constellation)
                printClose()
            printOpen("- removing requests")
            self.destroy.apply(LNS,constellation)
            LNS.notifyNoEvent(constellation)
            printClose()
            printOpen("- greedy filling")
            self.filling.apply(LNS,constellation)
            LNS.notifyNoEvent(constellation)
            newObjective = LNS.getObjective(constellation)
            if not self.acceptWorseSolution and newObjective<oldObjective:
                printOpen("New solution is worse. Backup ancient solution ...")
                self.backup.backupSolution(LNS)
                LNS.notifyBackupSolution(constellation)
                printClose()
            printClose()
            if not self.acceptWorseSolution:
                assert(LNS.getObjective()[0]>=oldObjective[0])    
            return oldObjective<newObjective
            #printClose()
  
class RequestDestructor(Perturbator):
        # filling rapide = filling post destruction peu couteuse
        def __init__(self,transitionModel,destroyRate,fastFilling=False):
            super().__init__("request destructor",transitionModel)
            self.fastFilling = fastFilling
            self.destroyRate = destroyRate
            mpiId = MPI.COMM_WORLD.Get_rank()
            configuredSeed = config.getOptValue("seed")
            self.randomInstance = rd.Random(mpiId+configuredSeed)
            
        def destroy(self,LNS,constellation):
            couvertes = LNS.fulfilledRequests()
            destroy = self.randomInstance.sample(couvertes,int(min(len(couvertes)*self.destroyRate,config.LNS.max_destruction)))
            modesToRemove = [(r,m) for (r,m) in LNS.getSelectedModes() if r in destroy and r!=-1]
            self.removeModes(LNS,constellation,modesToRemove)
            # mettre a jour les cca concernées = notifier les changements dans l'objet 'etat requete'
            listOfCCAs = []
            for (r,m) in modesToRemove:
                for (s,a) in constellation.getRequest(r).getPairs():
                    idCCA = LNS.getGraphComponents().getActivityCCA(a)
                    if idCCA not in listOfCCAs:
                        listOfCCAs.append(idCCA)
                        LNS.notifyCCAChangement(idCCA)          
            return [x[0] for x in modesToRemove],None
        # I. retirer des requetes
        # II. potentiellement re-remplir rapidement
        # III. annuler les modifications sur les etats des requetes :
        #   => si le solver est actif hors destruction il ne faut pas indiquer 
        #   les activites comme infaisable car le solver ici est glouton 
        #   (moins performant)
        def apply(self,LNS,constellation):
            oldObjective = LNS.getObjective(constellation)
            self.destroy(LNS,constellation)
            if self.fastFilling:
                if config.getOptValue("use_solver"):
                    LNS.saveActivitiesState()
                LNS.greedyFill(constellation,requestsLimit=config.getOptValue("use_solver"),forbidSolver=True,random=False)
                if config.getOptValue("use_solver"):
                    LNS.restoreActivitiesState()
            return LNS.getObjective(constellation)>oldObjective
class VoisinageCP(Operator):
    def __init__(self,nProcesses,transitionModel,k=2,iterOperator=1):
        super().__init__("CP neighborhood of a CCAs set",transitionModel,k=k)
        self.nProcesses = nProcesses
        #self.positions_cca = None
        self.iterOperator = iterOperator
        # initialisé à la 1ere utilisation
        self.connectCCAs = None
        self.requestsCCA = None
        self.observationsCCA = None
        self.gaps = []
        self.gains = []
        mpiId = MPI.COMM_WORLD.Get_rank()
        configuredSeed = config.getOptValue("seed")
        self.randomFixedModes = rd.Random(mpiId+configuredSeed)
    
    def addInfos(self,gap,gain):
        self.gaps.append(gap)
        self.gains.append(gain)
    
    def getGapsInfos(self):
        return self.gaps,self.gains
    
    def incrementNeighborhoodSize(self):
        self.k += 1
        
    def decrementNeighborhoodSize(self):
        self.k -= 1
        assert(self.k >0)
    
    def _extractModes(self,LNS,constellation,CCAs,presentRequests):
        modes = {}
        externalModes = {}
        externalActivities = {}
        newActivities = {}
        for r in presentRequests:
                newActivities[r] = {}
                m = LNS.getModeIfPresent(r)
                if m is None:
                    act = []
                else:
                    act = [x[1] for x in constellation.getRequest(r).getMode(m).getPairs()]
                currentContent = [x for x in act if LNS.getGraphComponents().getActivityCCA(x) not in CCAs]
                allowExternal = True # permet d'initialiser la sol avec ce mode
                modes[r],externalModes[r],externalActivities[r] = constellation.getRequest(r).generateModesInCCA(LNS.getGraphComponents(),constellation,currentContent,CCAs,allowExternalModes=allowExternal)
                for mode in modes[r]:
                    print(mode)
                    newActivities[r][mode.getId()] = [a for a in [x[1] for x in mode.getPairs()] if a not in currentContent]
                if len(newActivities[r])==0:
                    del newActivities[r]
        # supprimer les requêtes sans modes
        requestsToDelete = []
        for r in modes:
            if len(modes[r])==0:
                requestsToDelete.append(r)
        for r in requestsToDelete:
            del modes[r]
        if config.getOptValue("verif"):
            for r in newActivities:
                for m in newActivities[r]:
                    for a in newActivities[r][m]:
                        assert(LNS.getGraphComponents().getActivityCCA(a) in CCAs)
        return modes,newActivities,externalModes,externalActivities
    
    def checkEmptyCCA(self,LNS,newActivities,CCAs):
         CCAsObservations = []
         for r in newActivities:
             for m in newActivities[r]:
                 for a in newActivities[r][m]:
                     idCCA = LNS.getGraphComponents().getActivityCCA(a)
                     if idCCA not in CCAsObservations:
                         assert(idCCA in CCAs)
                         CCAsObservations.append(idCCA)
                         if len(CCAs)==len(CCAsObservations):
                             assert(sorted(CCAs)==sorted(CCAsObservations))
                             return False
         return len(CCAs)!=len(CCAsObservations)
     
    # gérer les modes qui ont des activités en dehors des CCA considérées
    def updateExternalModes(self,LNS,CCAs,constellation,externalModes,selectedModes):
        for r in externalModes:
            mode = externalModes[r]
            externalActivities = [(s,a) for (s,a) in constellation.getRequest(r).getMode(mode).getPairs() if LNS.getGraphComponents().getActivityCCA(a) not in CCAs]
            if constellation.getRequest(r).acceptable(externalActivities,constellation):
                mode = constellation.getRequest(r).addMode(constellation,externalActivities)
                selectedModes.append((r,mode.getId()))
            else:
                for (s,a) in externalActivities:
                    _,cca = LNS.getGraphComponents().getActivityCCA(a)
                    LNS.getSolCCA(s,cca).removeActivity(constellation,a)
                    
    def evalNewSolution(self,LNS,constellation,sol,externalModes,externalActivities,CCAs,newActivities,fixedModes):
        # retirer les anciens modes
        #print(sorted(list(self.modesVariables.keys())))
        newListOfModes = [(r,m) for (r,m) in LNS.getSelectedModes() if r not in self.modesVariables]
        externalModes = {r:m for (r,m) in LNS.getSelectedModes() if r in self.modesVariables and externalActivities[r]}
        # rajouter les nouveaux modes
        # MAJ la liste de modes satisfaits et eliminer les modes problematiques (hors CCA)
        for r in self.modesVariables:
            for m in self.modesVariables[r]:
                value = sol[self.modesVariables[r][m]]
                if value==1:
                    newListOfModes.append((r,m))
                    if r in externalModes:
                        del externalModes[r]
                        break
        # mettre a jour les modes problématiques (hors CCA)            
        self.updateExternalModes(LNS,CCAs,constellation,externalModes,newListOfModes)
        return self.updateSolutionIfBetter(LNS,constellation,sol,newListOfModes)
        
    def updateSolutionIfBetter(self,LNS,constellation,sol,newListOfModes):
        utilityBefore = LNS.getObjective()[0]
        copyModes = deepcopy(LNS.getSelectedModes())
        LNS.setSelectedModes(newListOfModes,constellation)
        delta_score = round(LNS.getObjective()[0]-utilityBefore,3)
        if not delta_score>=0:
            LNS.setSelectedModes(copyModes,constellation)
            LNS.verifySolutionIfVerifyMode(constellation) 
            return False,delta_score
        
        self.fillSolCCAs(constellation,LNS,sol)
        LNS.verifySolutionIfVerifyMode(constellation) 
        return True,delta_score
    
    def recordError(self,LNS,problem,sol,CCAs):                
        if len(os.listdir("error_logs/solver_errors/"))==0:
            idFile = 0
        else:
            idFile = max([int(file.split('_')[-1]) for file in os.listdir("error_logs/solver_errors")])+1
        path = "error_logs/solver_errors/model_"+problem+"_"+str(CCAs)+"_id_"+str(idFile)
        self.model.export_model(path)
        alert("EchangeCP : pas de solution trouvée pour "+str(CCAs)+" : id="+str(idFile))
        with open(path,"a") as file:
            for (s,cca) in CCAs:
                file.write(str(LNS.getSolCCA(s,cca))+"\n")
            file.write("modes retenus"+str(sorted(LNS.getSelectedModes()))+'\n')
            file.write("sol is None : "+str(sol is None)+'\n')
            if sol is not None:
                file.write("sol status : "+str(sol.get_solve_status())+'\n')
    
    def verifyRequests(self,externalModes,externalActivities,newActivities,CCAs):
        for r in newActivities:
            if externalActivities[r] and externalModes[r]:
                assert(len(newActivities)>1) # au moins le mode externe + des combinaisons sur les cca a explorer
            else:
                assert(len(newActivities)>=1) # au moins un mode
    
    def selectFixedModes(self,LNS,modes,newActivities,externalModes,externalActivities):
        requestsInSolution = [(r,m) for (r,m) in LNS.getSelectedModes() if r in newActivities]
        max_r = config.getOptValue("quota_requests")
        if max_r>=len(requestsInSolution):
            return {}
        else: # on retire -1 car le vidage est traité à part
            fixedRequests = self.randomFixedModes.choices(requestsInSolution,k=max_r)
            res = {}
            for (r,m) in requestsInSolution:
                if (r,m) not in fixedRequests:
                    res[r] = m
                    del modes[r]
                    del externalModes[r]
                    del externalActivities[r]
                    del newActivities[r]
            return res
    
    def choseCCAAndGenerateModes(self,LNS,constellation):
         existsCCAEmpty = True
         shiftRightDisplay(2)
         while existsCCAEmpty:
             printOpen("Sélection des CCA")
             CCAs,presentRequests = self.selectCCAs(LNS,constellation)
             printClose()
             printOpen("Génération des modes")
             modes,newActivities,externalModes,externalActivities = self._extractModes(LNS,constellation,CCAs,presentRequests)       
             fixedModes = self.selectFixedModes(LNS, modes,newActivities,externalModes,externalActivities)
             existsCCAEmpty = self.checkEmptyCCA(LNS,newActivities,CCAs)
             printClose()
             if existsCCAEmpty:
                 printColor("An empty CCA exists. Searching for a new candidate.",c='m')
         shiftLeftDisplay(2)
         return CCAs,presentRequests,fixedModes,modes,newActivities,externalModes,externalActivities
    
    def runSolver(self,LNS,CCAs):
         tlim = min(LNS.getRemainingTime(),config.getOptValue("CPLim"))
         succes = False
         sol = None
         if tlim>0:
            printOpen("CP model resolution",tlim,"(s)")
            t1 = time()
            maxTries = 2
            triesCount = 0
            while not succes and triesCount<maxTries:
                try:
                    if config.glob.docplexPath!="":
                        sol = self.model.solve(TimeLimit = tlim,execfile=config.glob.docplexPath,Workers=self.nProcesses,log_output=None)
                    else:
                        sol = self.model.solve(TimeLimit = tlim,Workers=self.nProcesses,log_output=None)
                    succes = True
                except Exception as solver_error:
                    warn(solver_error)
                    alert("Memory error in the CP Solver",CCAs)
                    self.recordError(LNS,"memory_pb",None,CCAs)
                    triesCount += 1
            t2 = time()
            printClose()
         return sol,succes
    
    def displayGapsInfo(self,LNS,constellation,sol,acceptable,delta_score):
        self.gains.append(delta_score)
        if acceptable:
            gap = sol.solution.get_objective_gap()
            self.gaps.append(gap)
            printColor("Gap : ",gap,c='y')
            if delta_score<=10e-4:
                printColor("Gain : " + str(delta_score),c='r')
            else:
                printColor("Gain : " + str(delta_score),c='g')
        else:
            LNS.verifySolutionIfVerifyMode(constellation) 
            self.gaps.append(-1)
            printColor("No better solution found before reaching the time limit.",c='y')
                        
    def apply(self,LNS,constellation):
        LNS.verifySolutionIfVerifyMode(constellation)
        self.choseCCAAndGenerateModes(LNS,constellation)
        oldObjective = LNS.getObjective(constellation)
        printOpen("Operator initialization.")
        CCAs,presentRequests,fixedModes,modes,newActivities,externalModes,externalActivities = self.choseCCAAndGenerateModes(LNS,constellation)
        printOpen("CP model construction")
        self.buildModel(LNS,constellation,CCAs,newActivities,modes,externalModes,externalActivities,fixedModes)
        printClose()
        printClose()
        sol,succes = self.runSolver(LNS,CCAs)
        if succes:
            if not (sol is not None and sol.get_solve_status()!="Infeasible" and sol.get_solve_status()!="Unknown"):
                self.recordError(LNS,"no-solution",sol,CCAs)
            else:
                LNS.verifySolutionIfVerifyMode(constellation)
                acceptable,delta_score = self.evalNewSolution(LNS,constellation,sol,externalModes,externalActivities,CCAs,newActivities,fixedModes)
                self.displayGapsInfo(LNS,constellation,sol,acceptable,delta_score)
        return LNS.getObjective(constellation)>oldObjective
            
    def getStats(self):
        return self.gaps,self.gains
    
    def buildModel(self,LNS,constellation,CCAs,newActivities,modes,externalModes,externalActivities,fixedModes):
        printOpen("Model construction")
        self.model = CpoModel()
        printClose()
        printOpen("IVariables initialization")
        self.initVars(LNS,constellation,modes,CCAs,newActivities,fixedModes)
        printClose()
        printOpen("Constraints initialization")
        self.initConstraints(LNS,constellation,CCAs,modes,externalModes,externalActivities,fixedModes)
        printClose()
        printOpen("Initial solution creation")
        self.initSol(LNS,constellation,CCAs,modes,fixedModes)
        printClose()
            
    def initVars(self,LNS,constellation,modes,CCAs,newActivities,fixedModes):
        self.modesVariables = {}                
        for r in modes:
            self.modesVariables[r] = {}
            for m in modes[r]:
                self.modesVariables[r][m.getId()] = self.model.binary_var(name="y_"+str((r,m.getId())))

        self.intervalVariables = {}
        for r in newActivities:
            self.intervalVariables[r] = {}
            for m in newActivities[r]:
                for a in newActivities[r][m]:
                    assert(LNS.getGraphComponents().getActivityCCA(a) in CCAs)
                    self.addIntervalVariable(constellation,r,a,True)
                    
        for r in fixedModes:
            self.intervalVariables[r] = {}
            m = fixedModes[r]
            for (s,a) in constellation.getRequest(r).getMode(m).getPairs():
                if LNS.getGraphComponents().getActivityCCA(a) in CCAs:
                    self.addIntervalVariable(constellation,r,a,False)

    def addIntervalVariable(self,constellation,r,a,optional):                    
        if a not in self.intervalVariables[r]:
            s = constellation.getSatelliteActivity(a)
            start = int(math.ceil(config.glob.getScale()*constellation.getSatellite(s).getActivity(a).getStartDate()))
            end = int(math.floor(config.glob.getScale()*constellation.getSatellite(s).getActivity(a).getEndDate()))
            duree = int(math.ceil(config.glob.getScale()*constellation.getSatellite(s).getActivity(a).getDuration()))
            self.intervalVariables[r][a] = self.model.interval_var(start=(start,end),end=(start,end),length=duree,name="Ia_"+str((r,a)))
            if optional:
                self.intervalVariables[r][a].set_optional()
                    
    def initNoOverlap(self,constellation,LNS,CCAs):
        self.sequenceCCA = {}
        for idCCA in CCAs:
            s = idCCA[0]
            printOpen("Searching related activities")
            activites = [(a,r,self.intervalVariables[r][a]) for r in self.intervalVariables for a in self.intervalVariables[r] if LNS.getGraphComponents().getActivityCCA(a) == idCCA]
            if len(activites)==0:
                #print("CCAs vides",self.checkEmptyCCA(LNS,newActivities,CCAs))
                die("No activity in the CCA",idCCA)
            
            printClose()
            n = len(activites)
            distanceMatrix = np.zeros((n,n),dtype=int)
            variableList = []
            printOpen("Transition matrix creation",str(n)+"x"+str(n),"activities")
            for i,(a1,r,var) in enumerate(activites):
                variableList.append(var)
                for j,(a2,r,var2) in enumerate(activites):
                    if i==j:
                        distanceMatrix[i][j] = 0
                    else:
                        transition = config.glob.getScale()*constellation.getSatellite(s).getTransitionDuration(a1,a2,self.transitionModel)
                        distanceMatrix[i][j] = int(math.ceil(transition))
            mat = self.model.transition_matrix(distanceMatrix)
            self.sequenceCCA[idCCA] = self.model.sequence_var(variableList)
            self.model.add(self.model.no_overlap(self.sequenceCCA[idCCA],distance_matrix=mat))
            printClose()
            
    def initConstraints(self,LNS,constellation,CCAs,modes,externalModes,externalActivities,fixedModes):
        printOpen("Creating coupling constraints modes - observations")
        self.initConstraintsPresenceModes(LNS,modes,externalModes,externalActivities,fixedModes)
        printClose()
        printOpen("Creating noOverlap constraints")
        self.initNoOverlap(constellation,LNS,CCAs)
        printClose()
        obj1 = self.model.sum([self.modesVariables[r][m.getId()]*m.getUtility() for r in modes for m in modes[r]])
        self.model.maximize(obj1)

    def initConstraintsPresenceModes(self,LNS,modes,externalModes,externalActivities,fixedModes):
        for r in self.modesVariables:
            if r==-1:
                self.model.add_constraint(self.model.sum([self.modesVariables[r][m] for m in self.modesVariables[r]])==1)
            else:
                self.model.add_constraint(self.model.sum([self.modesVariables[r][m] for m in self.modesVariables[r]])<=1)
        # presence mode = présence activites    
        for r in self.intervalVariables:
            if r not in fixedModes:
                for a in self.intervalVariables[r]:
                    activityModes = [mode.getId() for mode in modes[r] if a in [x[1] for x in mode.getPairs()]]
                    self.model.add_constraint(self.model.presence_of(self.intervalVariables[r][a])==self.model.sum([self.modesVariables[r][m] for m in activityModes]))
    
    def findDatePlan(self,plan,a):
        for activite,t in plan:
            if activite==a:
               return t
        return None
    
    def findMode(self,constellation,LNS,modes):
        initialMode = {}
        for r in modes:
           found = False
           for m in modes[r]:
               pairs = sorted(m.getPairs())
               for (rr,mm) in LNS.getSelectedModes():
                   if r ==rr:
                       if sorted(constellation.getRequest(rr).getMode(mm).getPairs())==pairs:
                           initialMode[r] = m.getId()
                           found = True
                           break
               if found:
                    break
        return initialMode
                   
    def initSol(self,LNS,constellation,CCAs,modes,fixedModes):
        warmstart=self.model.create_empty_solution()
        plan = {}
        for idCCA in CCAs:
            s,cca = idCCA
            plan[idCCA] = LNS.getSolCCA(s,cca).getEarliestPlan(constellation,LNS.getSolCCA(s,cca).getSequence(),self.transitionModel)
        initialMode = self.findMode(constellation,LNS,modes)
        for r in self.modesVariables:
            if r not in fixedModes:
                for m in self.modesVariables[r]:
                    if r in initialMode and initialMode[r]==m:
                        warmstart[self.modesVariables[r][m]]=1
                    else:
                        warmstart[self.modesVariables[r][m]]=0
        
        if config.getOptValue("verif"):
            for r in initialMode:
                    m = initialMode[r]
                    for (s,a) in constellation.getRequest(r).getMode(m).getPairs():
                        _,cca = LNS.getGraphComponents().getActivityCCA(a)
                        assert(a in LNS.getSolCCA(s,cca).sequence)
                        
            for s,cca in CCAs:
                for a in LNS.getSolCCA(s,cca).sequence:
                        r = constellation.getRequestActivity(a)
                        if r in initialMode:
                            m = initialMode[r]
                        else:
                            m = fixedModes[r]
                        obs_mode = [x[1] for x in constellation.getRequest(r).getMode(m).getPairs()]
                        assert(a in obs_mode)
        
        
        # presence des activites
        for r in self.intervalVariables:
            for a in self.intervalVariables[r]:
                idCCA = LNS.getGraphComponents().getActivityCCA(a)
                t = self.findDatePlan(plan[idCCA], a)
                s = constellation.getSatelliteActivity(a)
                if t is not None:
                    start = int(math.ceil(t*config.glob.getScale()))
                    duree = int(math.ceil(constellation.getSatellite(s).getActivity(a).getDuration()*config.glob.getScale()))
                    end = start+duree
                    warmstart.add_interval_var_solution(self.intervalVariables[r][a],presence=True,start=start,end=end)
                else:
                    warmstart.add_interval_var_solution(self.intervalVariables[r][a],presence=False)
        # impossible d'indiquer les modes de départ : les indices des modes présent ne sont pas les mêmes que les nouveaux    
        # indiquer les séquences
        seq = []
        for idCCA in CCAs:
            s,cca = idCCA
            for a in LNS.getSolCCA(s,cca).getSequence():
                assert(idCCA==LNS.getGraphComponents().getActivityCCA(a))
                r = constellation.getRequestActivity(a)
                if a not in self.intervalVariables[r]:
                    print(CCAs,idCCA)
                    print("a in the sequence but no related variable has been created?",r,a,LNS.getGraphComponents().getActivityCCA(a),constellation.getRequest(r).getType())
                    print(r in LNS.fulfilledRequests())
                    m = [mm for (rr,mm) in LNS.getSelectedModes() if rr==r][0]
                    print(constellation.getRequest(r).getMode(m))
                    print(r in fixedModes,r in modes)
                    print("modes of request r :")
                    for m in modes[r]:
                        print(m)
                seq.append(self.intervalVariables[r][a])
            warmstart[self.sequenceCCA[idCCA]] = seq
        self.warmstart = warmstart
        self.model.set_starting_point(warmstart)
        
    def fillSolCCAs(self,constellation,LNS,sol):
        for idCCA in self.sequenceCCA:
            sequence = []
            for it in sol.get_value(self.sequenceCCA[idCCA]):
                res = eval(it.get_name().split("_")[1])
                r,a = res[0],res[1]
                sequence.append(a)
            LNS.overwriteCCASolution(constellation,idCCA,sequence)
            s,cca = idCCA
            if config.verifMode():
                for a in sequence:
                        find = None
                        r = constellation.getRequestActivity(a)
                        for m in self.modesVariables[r]:
                            if sol[self.modesVariables[r][m]]==1:
                                if find is not None:
                                    raise ValueError("several modes for",r)
                                find = m
                        assert(find is not None)
                        assert(sol[self.modesVariables[r][find]]==1)
                        assert((r,find) in LNS.getSelectedModes())
                        s = constellation.getSatelliteActivity(a)
   
                
        if config.getOptValue("verif"):
            for idCCA in self.sequenceCCA:
                s,cca = idCCA
                for a in LNS.getSolCCA(s,cca).getSequence():
                    r = constellation.getRequestActivity(a)
                    find = None
                    for m in self.modesVariables[r]:
                        if sol[self.modesVariables[r][m]]==1:
                            find = m
                            break
                    assert((r,find) in LNS.getSelectedModes())
                    s = constellation.getSatelliteActivity(a)
                    assert((s,a) in constellation.getRequest(r).getMode(find).getPairs())
            LNS.verifySolution(constellation)
            
class DestroyAndRepairGreedyCCA(Operator):
    def __init__(self,LNS,transitionModel,acceptWorseSolution,k=2):
        super().__init__("Greedy windows on CCAs",transitionModel)
        self.acceptWorseSolution = acceptWorseSolution
        self.destroy = CCAsSetDestructor(LNS.positions_cca,transitionModel,k=k)
        self.filling = fillingGreedyWindow(transitionModel)        
        
    def apply(self,LNS,constellation):
        #printOpen("Perturbation :",c='y')
        oldObjective = LNS.getObjective(constellation)
        if not self.acceptWorseSolution:
            printOpen("Creating backup of the solution")
            self.backup = LNS.createBackup(constellation)
            printClose()
        printOpen("- removing requests")
        presentRequests,nouveaux_modes,CCAs = self.destroy.destroy(LNS,constellation)
        LNS.notifyNoEvent(constellation)
        printClose()
        printOpen("- greedy filling")
        self.filling.apply(LNS,constellation,CCAs,presentRequests)
            
        newObjective = LNS.getObjective(constellation)
        if not self.acceptWorseSolution and newObjective<oldObjective:
            printOpen("Solution is worse. Backup the previous solution.")
            self.backup.backupSolution(LNS)
            LNS.notifyBackupSolution(constellation)
            printClose()
        printClose()
        return oldObjective<newObjective

            
class CCAsSetDestructor(DestructionOperator):
    def __init__(self,positions_cca,transitionModel,k=2):
        super().__init__("destructor of a set of CCAs",transitionModel,k=k)
        self.positions_cca = positions_cca
       
    def destroy(self, LNS, constellation):
        CCAs,presentRequests = self.selectCCAs(LNS,constellation)
        satLSTE = [idCCA[0] for idCCA in CCAs]
        newModeContent = {}
        currentModes = {}
        LNS.verifySolution(constellation)
        for s,cca in CCAs:
            for a in LNS.getSolCCA(s,cca).getSequence():
                r = constellation.getRequestActivity(a)
                if r not in newModeContent:
                    selectedModes = [m for (rr,m) in LNS.getSelectedModes() if rr==r]
                    assert(len(selectedModes)==1)
                    currentModes[r] = selectedModes[0]
                    pairs = constellation.getRequest(r).getMode(selectedModes[0]).getPairs()
                    newModeContent[r] = [(sat,act) for (sat,act) in pairs if act!=a]
                else:
                    newModeContent[r] = [(sat,act) for (sat,act) in newModeContent[r] if act!=a]
        for s,cca in CCAs:    
            LNS.getSolCCA(s,cca).reset()
        for r in LNS.requestState:
            for idCCA in CCAs:
                LNS.requestState[r].resetCCA(idCCA)
        for r in newModeContent:
            m = currentModes[r] 
            if r==-1:
                for (s,a) in constellation.getRequest(r).getMode(m).getPairs():
                    idCCA = LNS.getGraphComponents().getActivityCCA(a)
                    if idCCA in CCAs:
                        LNS.solCCAs[s][idCCA[[1]]].insererActivitePlansCritiques(constellation,a)
                LNS.addSelectedMode((-1,0))
            else:
                LNS.retirerModeRetenu((r,m))
                if constellation.getRequest(r).acceptable(newModeContent[r],constellation):
                    mode = constellation.getRequest(r).addMode(constellation,newModeContent[r])
                    LNS.addSelectedMode((r,mode.getId()))
        if -1 in newModeContent:
            del newModeContent[-1]
        return presentRequests,list(newModeContent.keys()),CCAs

        
class fillingGreedyWindow(RepairingOperator):
    def __init__(self,transitionModel,k=2):      
        super("filling glouton time window",transitionModel)
        self.k = k # choix d'une cca1 et d'une cca2 parmi les k plus proches
        # initialisé à la 1ere utilisation
        self.requestsCCA = None
        self.gains = []
        mpiId = MPI.COMM_WORLD.Get_rank()
        configuredSeed = config.getOptValue("seed")
        self.randomFixedModes = rd.Random(mpiId+configuredSeed)
    
    def addInfos(self,gain):
        self.gains.append(gain)

    def _extractModes(self,LNS,constellation,CCAs,presentRequests):
        modes = {}
        externalModes = {}
        externalActivities = {}
        newActivities = {}
        for r in presentRequests:
                newActivities[r] = {}
                m = LNS.getModeIfPresent(r)
                if m is None:
                    act = []
                else:
                    act = [x[1] for x in constellation.getRequest(r).getMode(m).getPairs()]
                currentContent = [x for x in act if LNS.getGraphComponents().getActivityCCA(x) not in CCAs] # activités externes
                allowExternal = True # permet d'initialiser la sol avec ce mode
                modes[r],externalModes[r],externalActivities[r] = constellation.getRequest(r).generateModesInCCA(LNS.getGraphComponents(),constellation,currentContent,CCAs,allowExternalModes=allowExternal)
                for mode in modes[r]:
                    newActivities[r][mode.getId()] = [a for a in [x[1] for x in mode.getPairs()] if a not in currentContent]
                if len(newActivities[r])==0:
                    del newActivities[r]
        # supprimer les requêtes sans modes
        requestsToDelete = []
        for r in modes:
            if len(modes[r])==0:
                requestsToDelete.append(r)
        for r in requestsToDelete:
            del modes[r]
        if config.getOptValue("verif"):
            for r in newActivities:
                for m in newActivities[r]:
                    for a in newActivities[r][m]:
                        assert(LNS.getGraphComponents().getActivityCCA(a) in CCAs)
        return modes,newActivities,externalModes,externalActivities
    
    def evalNewSolution(self,LNS,constellation,sol,externalModes,externalActivities,CCAs,newActivities,fixedModes):
        # retirer les anciens modes
        #print(sorted(list(self.modesVariables.keys())))
        newListOfModes = [(r,m) for (r,m) in LNS.getSelectedModes() if r not in self.modesVariables]
        retraits = {r:m for (r,m) in LNS.getSelectedModes() if r in self.modesVariables and self.hasProblematicExternActivities(externalActivities,externalModes,r)}
        # rajouter les nouveaux modes
        for r in self.modesVariables:
            for m in self.modesVariables[r]:
                value = sol[self.modesVariables[r][m]]
                if value==1:
                    newListOfModes.append((r,m))
                    if r in retraits:
                        del retraits[r]
                        break
        utilityBefore = LNS.getObjective()[0]
        copyModes = deepcopy(LNS.getSelectedModes())
        LNS.setSelectedModes(newListOfModes,constellation)
        delta_score = round(LNS.calculerObjectif(constellation)[0]-utilityBefore,3)
        if not delta_score>=0:
            LNS.setSelectedModes(copyModes,constellation)
            return False,delta_score
        
        for r in retraits:
            mode = retraits[r]
            for (s,a) in constellation.getRequest(r).getMode(mode).getPairs():
                idCCA = LNS.getGraphComponents().getActivityCCA(a)
                LNS.solCCAs[s][idCCA].removeActivity(constellation,a)
        return True,delta_score
    

    def verifyRequests(self,externalModes,externalActivities,newActivities,CCAs):
        for r in newActivities:
            if externalActivities[r] and externalModes[r]:
                assert(len(newActivities)>1) # au moins le mode externe + des combinaisons sur les cca a explorer
            else:
                assert(len(newActivities)>=1) # au moins un mode

    def selectMode(self,LNS,constellation,modes):
        currentUtility = {}
        for (r,m) in LNS.getSelectedModes():
            if r in modes:
                currentUtility[r] = constellation.getRequest(r).getMode(m).getUtility()
        cle = lambda rm:rm[1].getUtility()-currentUtility.get(r,0)
        (req,mode) = max([(r,m) for r in modes for m in modes[r]],key=cle)
        modes[req].remove(mode)
        if len(modes[req])==0:
            del modes[req]
        return req,mode

    def validateMode(self,modes,r):
        if r in modes:
            del modes[r]

    def apply(self,LNS,constellation,CCAs,presentRequests):
         oldObjective = LNS.getObjective(constellation)
         shiftRightDisplay(2)
         printOpen("Modes generation")
         modes,newActivities,externalModes,externalActivities = self._extractModes(LNS,constellation,CCAs,presentRequests)       
         printClose()
         printOpen("Modes insertion")
         if -1 in modes:
             del modes[-1]
         while len(modes.keys())>0 and LNS.getTempsEcoule() < config.getOptValue("time"):
             r,mode_candidat = self.selectMode(LNS,constellation,modes) 
             printOpen("Attempt to insert request",r)
             mode,solutionCCA,vide = LNS.requestState[r].tryInsertMode(LNS,constellation,r,mode_candidat,False)
             #printOpen("Insérer la requete",constellation.getRequest(r).getType(),r,c='g')
             if mode is not None:
                 #printColor("valider",noeuds_a_valider,c='g')
                 self.validateMode(modes,r)
                 LNS.validerSol(solutionCCA)
                 LNS.setSelectedModes([(rr,m) for (rr,m) in LNS.getSelectedModes() if r!=r],constellation)
                 LNS.addSelectedMode((r,mode.getId()))
                 if config.getOptValue("verif"):
                     LNS.verifyCCA(constellation)
                 printClose("Success",c='g')
             else:
                printColor("Infeasible mode",c='m')
                printClose("Failure",c='r')
         printClose()
         shiftLeftDisplay(2)
         newObjective = LNS.getObjective(constellation)
         delta_score = newObjective - oldObjective
         if config.getOptValue("verif"):
             LNS.verifySolution(constellation)
         self.gains.append(delta_score)
         return newObjective>oldObjective

    def findDatePlan(self,plan,a):
        for activite,t in plan:
            if activite==a:
               return t
        return None
    
    def findMode(self,constellation,LNS,modes):
        initialMode = {}
        for r in modes:
           found = False
           for m in modes[r]:
               pairs = sorted(m.getPairs())
               for (rr,mm) in LNS.getSelectedModes():
                   if r ==rr:
                       if sorted(constellation.getRequest(rr).getMode(mm).getPairs())==pairs:
                           initialMode[r] = m.getId()
                           found = True
                           break
               if found:
                    break
        return initialMode
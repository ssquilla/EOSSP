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
    from EOSSP.model.graph import *
    
    from time import time
    from time import sleep
    
    class EtatRecherche:
        def __init__(self):
            self.activitiesDone = {}
            self.activitiesWorking = {}
            self.activitiesToDo = {}
            self.workingModes = []
            self.modeState = {}
            self.failedRequests = []
            self.WORKING_STATE = 0
            self.DONE_STATE = 1
            self.FAILED_STATE = 2
            self.WAITING_RESTART_STATE = 3
            self.WAITING_VALIDATION_STATE = 4
            self.openModes = []
            self.waitingModes = []
            self.CCAsUnderComputation = []
            self.failureResponsibleModes = {} # (r,m) ceux qui le bloquent
            self.toRestart = {}
            self.toDelete = []
            self.blockings = {} # failureResponsibleModes mais dans l'autre sens : (r,m) ceux que (r,m) bloque
            #self.waitingRestart = {}    
                
        def freeCCA(self,idCCA):
            self.CCAsUnderComputation.remove(idCCA)
            
        def notifyRequestInsertFailure(self,r):
            if r not in self.failedRequests:
                self.failedRequests.append(r)
        
        def setModeWaitingForValidation(self,mode):
            self.modeState[mode] = self.WAITING_VALIDATION_STATE
            if mode not in self.waitingModes:
                self.waitingModes.append(mode)
        
        # activitiesList: à indiquer comme deja en travail (pour le transfert de mode)
        def notifyActivitiesPreWorking(self,r,m,activitiesList):
            for a in activitiesList:
                self.activitiesToDo[(r,m)].remove(a)
                self.activitiesWorking[(r,m)].append(a)
                
        def getModesWaitingForValidation(self):
            return self.waitingModes

        def setModeWaitingForRestart(self,mode):
            self.modeState[mode] = self.WAITING_RESTART_STATE
        
        def getWorkingModes(self):
            return self.workingModes

        def getActivitiesToPlan(self,r,m):
            return self.activitiesToDo[(r,m)]!=[]
        
        def getStartedModes(self):
            return self.openModes
        
        def getActivitiesToDelete(self):
            return self.toDelete
        
        def notifyObservationDeletion(self,s,o):
            self.toDelete.remove((s,o))
            
        # suppression du mode pour la version avec anticipation
        def deleteModesForAnticipationVersion(self,constellation,r,m,grapheDep,listToDelete=None):
            self.finishModeScheduling(r,m,False)
            for (s,o) in constellation.getRequest(r).getMode(m).getPairs():
                if listToDelete is None or o in listToDelete:
                    if((s,o) not in self.toDelete):
                        idCCA = grapheDep.getActivityCCA(o)
                        if idCCA in self.CCAsUnderComputation:
                            self.toDelete.append((s,o))
            
        # annule un mode et renvoie la liste de ceux a restart. SUPRESSION DU MODE
        def notifyCertainFailure(self,constellation,r,m,grapheDep,listToDelete=None):
            self.deleteModesForAnticipationVersion(constellation,r,m,grapheDep,listToDelete)
            if not config.getOptValue("anticipation"):
                return []
            # MAJ des requests fautives/activites a restart
            if r in self.failureResponsibleModes:
                for fautif in self.failureResponsibleModes[r]:
                    self.blockings[fautif].remove(r)
                del self.failureResponsibleModes[r]
            if (r,m) in self.toRestart:
                del self.toRestart[(r,m)]
            if r in self.blockings:
                restartList = self.blockings[r].copy()
                for bloc in self.blockings[r]:
                    self.failureResponsibleModes[bloc].remove(r)
                del self.blockings[r]
            else:
                return []
            #clean les modes a restart : eviter un 2e appel au redemarrage
            for r in restartList:
                for fautif in self.failureResponsibleModes[r]:
                    self.blockings[fautif].remove(r)
                if r in self.failureResponsibleModes:
                    del self.failureResponsibleModes[r]
            return restartList 
        
        # indiquer un echec d'activités causé probablement par des requests fautives.
        def notifyUncertainFailure(self,r,m,failureResponsibleModes,activitiesToRestart,waitingForRestart=True):
            if len(failureResponsibleModes)==0:
                if(waitingForRestart):
                    die(r,m,failureResponsibleModes,activitiesToRestart,waitingForRestart)
                printOpen('Mode insertion failure',(r,m),'case 1 (transfer): restart',activitiesToRestart,c='m')
            else:
                assert(waitingForRestart)
                printOpen('Mode insertion failure',(r,m),'case 2 (not top-ranked mode): restart scheduling',activitiesToRestart,c='m')
            #printColor("Reprogrammer",(r,m),activitiesToRestart,c='c')
            assert(config.getOptValue("anticipation"))
            for request in failureResponsibleModes:
                if r not in self.failureResponsibleModes:
                    self.failureResponsibleModes[r] = []
                if request not in self.failureResponsibleModes[r]:
                    self.failureResponsibleModes[r].append(request)
                if request not in self.blockings:
                    self.blockings[request] = []
                if r not in self.blockings[request]:
                    self.blockings[request].append(r)
            for a in activitiesToRestart:
                if (r,m) not in self.toRestart:
                    self.toRestart[(r,m)] = []
                if a not in self.toRestart[(r,m)]:
                    self.toRestart[(r,m)].append(a)
                    self.activitiesWorking[(r,m)].remove(a)
            self.setModeWaitingForRestart((r,m))
            if not waitingForRestart:
                self.restartModeSchedule(r,m)
            printClose()
            
        def finishModeScheduling(self,r,m,sucess):
            if sucess:
                self.modeState[(r,m)] = self.DONE_STATE
            else:
                assert(self.modeState[(r,m)]) != self.FAILED_STATE
                self.modeState[(r,m)] = self.FAILED_STATE
            if (r,m) in self.activitiesDone:
                del self.activitiesDone[(r,m)]
            if (r,m) in self.activitiesToDo:
                del self.activitiesToDo[(r,m)]
            if (r,m) in self.activitiesWorking:
                del self.activitiesWorking[(r,m)]
            if (r,m) in self.openModes:
                self.openModes.remove((r,m))
            if (r,m) in self.waitingModes:
                self.waitingModes.remove((r,m))
            self.workingModes.remove((r,m))
        
        def ccaToDo(self,r,m,idCCA,grapheDep):
            for a in self.activitiesToDo[(r,m)]:
                to_do_cca = grapheDep.getActivityCCA(a)
                if to_do_cca == idCCA:
                    return True
            return False
        """
                features
        """
        def filterModesInScheduling(self,parents):
            return [x for x in parents if x in self.activitiesToDo and self.activitiesToDo[x]!=[] and self.modeState[x] == self.WORKING_STATE]
        
        def getModesInScheduling(self):
            return len(self.workingModes)>0
         
        def getCCAsUnderComputation(self,idCCA):
            return idCCA in self.CCAsUnderComputation
        
        def getCCAAvailableForComputation(self,r,m,grapheDep):
            ccas = [grapheDep.getActivityCCA(a) for a in self.activitiesToDo[(r,m)]]
            freeCCAs = [idCCA for idCCA in ccas if idCCA not in self.CCAsUnderComputation]
            return freeCCAs
        
        def getCCA(self,r,m,grapheDep,idCCA):
            if (r,m) not in self.openModes:
                self.openModes.append((r,m))
            assert(idCCA not in self.CCAsUnderComputation)
            act = [a for a in self.activitiesToDo[(r,m)] if grapheDep.getActivityCCA(a)==idCCA]
            supprimerListeElements(self.activitiesToDo[(r,m)],act)
            self.CCAsUnderComputation.append(idCCA)
            for a in act:
                self.activitiesWorking[(r,m)].append(a)
            return act
            
        # depile les activites et les renvoie
        """
        def getCCASuivante(self,r,m,grapheDep):
            assert(len(self.activitiesToDo[(r,m)])>0)
            if (r,m) not in self.openModes:
                self.openModes.append((r,m))
            ccas = [grapheDep.getActivityCCA(a) for a in self.activitiesToDo[(r,m)]]
            freeCCAs = [idCCA for idCCA in ccas if idCCA not in self.CCAsUnderComputation]
            idCCA = freeCCAs[0]
            self.CCAsUnderComputation.append(idCCA)
            liste_a = []
            act = [a for a in self.activitiesToDo[(r,m)] if grapheDep.getActivityCCA(a)==idCCA]
            supprimerListeElements(self.activitiesToDo[(r,m)],act)
            for a in act:
                assert(a not in self.activitiesToDo[(r,m)])
            for a in act:
                assert (a not in self.activitiesToDo[(r,m)])
            for a in act:
                self.activitiesWorking[(r,m)].append(a)
            return act
        """
        def renameCCA(self,ccaBefore,ccaAfter):
            if ccaBefore in self.CCAsUnderComputation:
                self.CCAsUnderComputation.remove(ccaBefore)
                self.CCAsUnderComputation.append(ccaAfter)
                
        def splitCCA(self,ccaBefore,ccaList):
            if ccaBefore in self.CCAsUnderComputation:
                self.CCAsUnderComputation.remove(ccaBefore)
                for idCCA in ccaList:
                    self.CCAsUnderComputation.append(idCCA)
        
        def isRequestScheduleFailed(self,r):
            return r in self.failedRequests
        
        def getModeState(self,r,m):
            return self.modeState[(r,m)]
        
        # ajoute un mode directement actif
        def addActiveMode(self,r,m,activites):
            self.activitiesDone[(r,m)] = []
            self.activitiesWorking[(r,m)] = []
            self.activitiesToDo[(r,m)] = activites
            self.workingModes.append((r,m))
            self.modeState[(r,m)] = self.WORKING_STATE
        
        def restartModeSchedule(self,r,m):
            assert(self.modeState[(r,m)] in [self.WAITING_RESTART_STATE,self.WORKING_STATE])
            #self.activitiesDone[(r,m)] = []
            #self.activitiesWorking[(r,m)] = []
            self.activitiesToDo[(r,m)] += self.toRestart[(r,m)]
            if (r,m) in self.waitingModes:
                self.waitingModes.remove((r,m))
            self.modeState[(r,m)] = self.WORKING_STATE
            del self.toRestart[(r,m)]
            #printColor('mode redémarré : ',(r,m),'restantes :',self.activitiesToDo[(r,m)],'réussies :',self.activitiesDone[(r,m)],'en cours :',self.activitiesWorking[(r,m)],c='y')
            #self.workingModes.append((r,m))
        
        # si anticiper_modes : renvoie la liste de requests a restart
        def validateMode(self,r,m):
            self.finishModeScheduling(r,m,True)
            if(config.getOptValue("anticipation")):
                requestsToCancel = []
                if r in self.blockings:
                    for requete in self.blockings[r]:
                        self.failureResponsibleModes[requete].remove(r)
                        if self.failureResponsibleModes[requete] == []:
                            requestsToCancel.append(requete)
                    del self.blockings[r]
                return requestsToCancel
            else:
                return []
        
        def prevalidateActivities(self,activites,r,m):
            self.openModes.append((r,m))
            printColor('Prevalidation of activities',activites,"from mode",(r,m),c='m')
            supprimerListeElements(self.activitiesToDo[(r,m)],activites)
            if (r,m) not in self.activitiesDone:
                self.activitiesDone[(r,m)] = []
            self.activitiesDone[(r,m)] += activites
            if len(self.activitiesToDo[(r,m)])==0 and len(self.activitiesWorking[(r,m)])==0:
                return True
            else:
                return False
            
        def validateActivities(self,activites,r,m):
            printColor('Validation of activities',activites,"from mode",(r,m),c='c')
            supprimerListeElements(self.activitiesWorking[(r,m)],activites)
            if (r,m) not in self.activitiesDone:
                self.activitiesDone[(r,m)] = []
            self.activitiesDone[(r,m)] += activites
            if len(self.activitiesToDo[(r,m)])==0 and len(self.activitiesWorking[(r,m)])==0:
                return True
            else:
                return False
    
        def cancelMode(self,r,m):
            self.finishModeScheduling(r,m,False)
    
    
    class UnitParallelCCASearch(Solver):
        def __init__(self,constellation,startDate,transitionModel,dt_construction_transition,tlim,CCAs,solution):
            super().__init__(constellation,startDate,transitionModel,dt_construction_transition,tlim=tlim,CCAs=CCAs,solution=solution)
            # si 1 coeur : slave local
            comm = MPI.COMM_WORLD
            if comm.Get_size()==1:
                self.localSlave = Slave(constellation,transitionModel)
            printColor("Initialization duration:",time()-self.startDate,c='b')
      
        def _addActiveModeToResearchState(self,constellation,r,m,prevalidation=[]):
            activites = []
            for (s,o) in constellation.getRequest(r).getMode(m).getPairs():
                if o not in activites:
                    activites.append(o)
            self.researchState.addActiveMode(r,m,activites) 
            if len(prevalidation)>0:
                self.researchState.prevalidateActivities(prevalidation,r,m)
                if not self.researchState.getActivitiesToPlan(r,m):
                    self.validateMode(constellation,(r,m))
    
        def addActiveMode(self,constellation,r,m,prevalidation=[]):
            self._addActiveModeToResearchState(constellation,r,m,prevalidation)
            
        def isPresentCCA(self,constellation,r,m,idCCA,dependencyGraph):
            for s in constellation.getRequest(r).getMode(m).getActivities():
                for o in constellation.getRequest(r).getMode(m).getSatellite(s).getActivities():
                    if dependencyGraph.getVertex(o).getConnectedComponent()==idCCA:
                        return True
            return False
        
        def getConnectedComponentOfActivity(self,a):
            return self.getGraphComponents().getActivityCCA(a)
        
        def setCCA(self,constellation,s,seq,idCCA,r,m):
            s,cca = idCCA
            self.solution.getSolCCA(s,cca).setSequence(constellation,seq,self.transitionModel)
    
        def validateActivities(self,constellation,activites,mode):
            s = constellation.getSatelliteActivity(activites[0])
            act = {a : None for a in activites}
            idCCA = self.getConnectedComponentOfActivity(activites[0])
            seq = []
            for i,a in enumerate(activites):
                cca_a = self.getConnectedComponentOfActivity(a)
                if cca_a!=idCCA:
                    self.setCCA(constellation,s,seq,idCCA,mode[0],mode[1])
                    idCCA = cca_a
                    seq = []
                seq.append(a)
                if i==len(activites)-1:
                    self.setCCA(constellation,s,seq,cca_a,mode[0],mode[1])
                
            
        """
            =============================================== 
                            RESOLUTION
            =============================================== 
        """
       
        def validateMode(self,constellation,mode):
            printOpen("Validation of mode "+str(mode),c='g')
            startVal = time()
            # retirer la requete de toutes ses cca
            for idCCA in self.ccaRequests[mode[0]]:
                self.removeRequestFromCCA(idCCA,mode[0])
            del self.ccaRequests[mode[0]]
            
            # retirer les obs de la requete non presentes dans le mode
            removals = []
            for s in constellation.getRequest(mode[0]).getActivities():
                for o in constellation.getRequest(mode[0]).getActivitiesBySatellite(s):
                    presentGraphe = o in self.getGraphComponents().getVertices() # explainations retirees...
                    nonPresentMode = o not in constellation.getRequest(mode[0]).getMode(mode[1]).getActivitiesBySatellite(s)
                    if presentGraphe and nonPresentMode:
                        removals.append(o)
            del self.currentModes[mode[0]]
            self.solution.history.registerValidatedRequest(constellation,mode[0])
            requestsToCancel = self.researchState.validateMode(mode[0],mode[1])
            self.solution.addSelectedMode(mode,constellation)

            if len(requestsToCancel)>0:
                printColor("Cancel requests",requestsToCancel,c='m')
            for r in requestsToCancel:
                mode = (r,self.currentModes[r])
                explaination = self.explaination[mode]
                self.scheduleNextMode(constellation,mode,explaination)
                del self.explaination[mode]
            printClose()
        
        def isModeCurrentlyCandidate(self,mode):
            rCurrent = mode[0]
            return rCurrent in self.currentModes and mode[1] == self.currentModes[rCurrent]
        
        def checkPropagateSuccessPrevalidationMode(self,mode,success):
            return config.getOptValue("preval") and not self.isModeCurrentlyCandidate(mode) and success
        
        def testPropagateFailurePrevalidationMode(self,mode,success):
            return not success and not self.isModeCurrentlyCandidate(mode)
            
        def analyzeMessage(self,data,constellation,freeCPUs,mailbox,iteration):
            start = time()
            self.solution.history.registerMessageReceiving(data)
            i = data['source']
            assert(i not in freeCPUs)
            bisect.insort(freeCPUs,i)
            it = data['iteration']
            if iteration==it:
                cpuTime,cca,succes,mode = data['cpu_time'],data['cca'],data['faisable'],data['mode']
                explaination = data['explaination']
                seq = cca.getSequence()
                # si message d'un mode passé : transfert potentiel d'activités validées ou à restart
                if self.checkPropagateSuccessPrevalidationMode(mode,succes):
                    mode = self.transferSuccessActivities(explaination,seq,mode)
                if self.testPropagateFailurePrevalidationMode(mode,succes):
                    removals = self.deleteAncientActivities(seq) # pour version anticipee de UPCCAS
                    mode,activitiesToRestart = self.transferFailedActivities(constellation,explaination,removals,mode)
                else:
                    activitiesToRestart = []
                # indiquer la cca libre et supprimer les activités n'étant plus à programmer
                self.freeCCA(seq,explaination)
                if self.isModeCurrentlyCandidate(mode):
                    # traiter le message d'un mode en cours
                    if succes and not self.researchState.isRequestScheduleFailed(mode[0]) and self.researchState.getModeState(mode[0],mode[1])!=self.researchState.FAILED_STATE:
                        statValidation = time()
                        if(len(seq)==0):
                            print(mode)
                        assert(len(seq)>0)
                        self.validateActivities(constellation,seq,mode)
                        startValidation = time()
                        termine = self.researchState.validateActivities(explaination,mode[0],mode[1])
                        self.solution.history.registerSuccess(i,cpuTime)

                        if termine:
                            if config.getOptValue("anticipation"):
                                self.notifyAnticipationSuccess(constellation,mode)
                            else:
                                self.validateMode(constellation,mode)
                    else:
                        if self.researchState.getModeState(mode[0],mode[1])!=self.researchState.FAILED_STATE:
                            if config.getOptValue("anticipation"):
                                self.notiyAnticipationFailure(constellation,mode,explaination,activitiesToRestart)
                            else:
                                self.scheduleNextMode(constellation,mode,explaination)
                        self.solution.history.registerFailure(i,cpuTime)
                    #die(self.solution)
                        
        def analyzeResults(self,constellation,freeCPUs,mailbox,iteration):
            comm = MPI.COMM_WORLD
            rank = comm.Get_rank()
            startIteration = time()
            if comm.Get_size()==1:
                messages = [self.localSlave.getResult()]
            else:
                messages = mailbox.readMessages()
            for data in messages:
                if data is not None:
                    self.analyzeMessage(data,constellation,freeCPUs,mailbox,iteration)
                    
        def deleteAncientActivities(self,seq):
            supp = []
            for (s,o) in self.researchState.getActivitiesToDelete():
                if o in seq:
                    seq.remove(o)
                    supp.append(o)
                    self.researchState.notifyObservationDeletion(s,o)
            if len(supp)>0:
                printColor('deleting activities in received sequence',supp,"(deleted modes)",c='r')
            return supp
        
        def transferFailedActivities(self,constellation,explaination,removal,mode):
            if mode[0] in self.currentModes:
                newMode = self.currentModes[mode[0]]
                if sum([a in [x[1] for x in constellation.getRequest(mode[0]).getMode(newMode).getPairs()] for a in explaination])==0:
                     return mode,[]
                if newMode!=mode[1]:
                    toRestart = explaination.copy()
                    transfert = []
                    for a in removal:
                        if a in toRestart:
                            toRestart.remove(a)
                    printColor('mode',str(mode)+': activities transfer failed',transfert,'- restart',toRestart,c='r')
                    return (mode[0],newMode),toRestart
                return mode,[]
            return mode,[]
        
        # voir si la planification d'une activité de l'ancien mode peut etre conservée pour le nouveau
        def transferSuccessActivities(self,explaination,seq,mode):
            if mode[0] in self.currentModes:
                newMode = self.currentModes[mode[0]]
                if newMode!=mode[1]:
                    find = False
                    for new in explaination:
                        if new in self.researchState.activitiesWorking[(mode[0],newMode)]:
                            find = True
                            break # suppression deja effectuée au dessus
                    if find:
                        printColor('mode',str(mode)+": activities transfer succeeded",c='g')
                        return (mode[0],newMode)
                    else:
                        return mode
                else:
                    return mode
            return mode
                
        def notifyAnticipationSuccess(self,constellation,mode):
            self.researchState.setModeWaitingForValidation(mode)
            self.releaseWaitingModes(constellation)    
        
        # Traiter un echec : considerer un redemarrage si le mode est prioritaire : il n'y a pas de fautif surclassé
        # activites a restart : car l'echec est acompagnée de la suppression d'autres activités.
        # => on leur redonne alors une chance
        def notiyAnticipationFailure(self,constellation,mode,explaination,activitiesToRestart):
            # I) si on a des activites a restart : on les redemarre et ne cherche pas d'autre mode fautif
            if len(activitiesToRestart): # pas d'attente de redemarrage : on redemarre de suite
                self.researchState.notifyUncertainFailure(mode[0],mode[1],[],activitiesToRestart,False)
                return
            
            # II) si on a pas d'activités à redémarrer de l'ancien mode : on cherche des modes failureResponsibleModes et on redemmarre
            # cca des activites calculees
            ccas = []
            for a in explaination:
                idCCA = self.getGraphComponents().getActivityCCA(a)
                if idCCA not in ccas:
                    ccas.append(idCCA)
            # detection de failureResponsibleModes : en attente s'il y en a
            topRankedOnAllCCAs = True
            for idCCA in ccas:
                if not self.premierCCA(mode[0],idCCA):
                    topRankedOnAllCCAs = False
                    self.explaination[mode] = self.explaination.get(mode,[])
                    for a in explaination:
                        if a not in self.explaination[mode]:
                            self.explaination[mode].append(a)
                    failureResponsibleModes = self.findfailureResponsibleModesCCA(constellation,mode[0],idCCA)
                    self.researchState.notifyUncertainFailure(mode[0],mode[1],failureResponsibleModes,explaination,True)
                    break
                    
            # III) si pas de fautif : annuler car le mode ne passera jamais
            if topRankedOnAllCCAs:
                self.scheduleNextMode(constellation,mode,explaination)
                self.explaination[mode] = []
        
        def findfailureResponsibleModesCCA(self,constellation,r,idCCA):
            m = self.currentModes[r]
            orderRequest = RequestOrdering(r,m,self.noiseUtility(constellation,r,m))
            return [req.getId() for req in self.requestsCCA[idCCA] if req>orderRequest]
        
        def premierCCA(self,r,idCCA):
            return r==max(self.requestsCCA[idCCA]).getId()
        
        # libere la cca et supprimer de la séquence les activités des modes supprimés
        def freeCCA(self,seq,explaination):
            ccas = []
            for a in seq+explaination:
                idCCA = self.getGraphComponents().getActivityCCA(a)
                if idCCA not in ccas:
                    ccas.append(idCCA)
                    self.researchState.freeCCA(idCCA)
                            
        def releaseWaitingModes(self,constellation):
            validation = True
            while validation:
                validation = False
                for waitingMode in self.researchState.getModesWaitingForValidation():
                    if self.isRequestPrioritary(waitingMode[0]):
                        validation = True
                        self.validateMode(constellation,waitingMode)
                        if not(waitingMode not in self.researchState.getModesWaitingForValidation()):
                            die(waitingMode,self.researchState.getModesWaitingForValidation())
                        if not(waitingMode not in self.researchState.getWorkingModes()):
                            die(waitingMode,self.researchState.getWorkingModes())
                        assert(self.researchState.getModeState(waitingMode[0],waitingMode[1])==self.researchState.DONE_STATE)
                        break
             
        def intersectionActivites(self,constellation,r):
            newMode = self.currentModes.get(r)
            previousMode = newMode - 1
            activitiesOldMode = [x[1] for x in constellation.getRequest(r).getMode(previousMode).getPairs()]
            if newMode is None:
                return activitiesOldMode,[],[]
            activitiesNewMode = [x[1] for x in constellation.getRequest(r).getMode(newMode).getPairs()]
            
            removal = []
            activitiesToPrevalidate = []
            intersectionWorkingActivities = []
            for o in activitiesNewMode:
                if o in activitiesOldMode and o in self.researchState.activitiesDone.get((r,previousMode),[]):
                    if config.getOptValue("preval") or config.getOptValue("anticipation"):
                        activitiesToPrevalidate.append(o)
                    else:
                        removal.append(o)
                elif o in self.researchState.activitiesWorking.get((r,previousMode),[]):
                    if config.getOptValue("preval") or config.getOptValue("anticipation"):
                        intersectionWorkingActivities.append(o)
            for o in activitiesOldMode:
                if o not in activitiesNewMode:
                    removal.append(o)
            return removal,activitiesToPrevalidate,intersectionWorkingActivities
        
        def scheduleNextMode(self,constellation,mode,explaination):
            exp_cca = [(o,self.getGraphComponents().getActivityCCA(o)) for o in explaination]
            printOpen('Mode insertion failure',mode,'case 3: cancel. cause : [activity, (sat_cca,n_cca) ... ] = ',exp_cca,c='r')
            self.solution.history.registerFailure(constellation,mode[0])          
            for idCCA in self.ccaRequests[mode[0]]:
                self.removeRequestFromCCA(idCCA,mode[0])
            modesToRestart = self.moveOnRequestState(constellation,mode,explaination,self.getGraphComponents()) # MAJ etat recherche MAJ graphe (pending)
            self.insertRequestInCCA(constellation,mode[0])
            if len(modesToRestart)>0:
                printColor('Restart scheduling for requests',modesToRestart,c='m')
            if config.getOptValue("anticipation"):
                for mode in modesToRestart:
                    self.researchState.restartModeSchedule(mode,self.currentModes[mode])
            printClose()
    
        def insertRequestInCCA(self,constellation,r):
            if r in self.currentModes:
                requete = constellation.getRequest(r)
                m = self.currentModes[r]
                orderRequest = RequestOrdering(r,m,self.noiseUtility(constellation,r,m))
                ccas = []
                for (s,o) in requete.getPairs():
                    try:
                        idCCA = self.getGraphComponents().getActivityCCA(o)
                        if idCCA not in ccas:
                            ccas.append(idCCA)
                    except:
                        pass
                self.ccaRequests[r] = ccas
                #prio = constellation.getRequest(r).getPriorite()
                #rec = self.noiseUtility(constellation,r,self.currentModes[r])
                for idCCA in ccas:
                    if idCCA not in self.requestsCCA:
                        self.requestsCCA[idCCA] = []
                    reverse_insort(self.requestsCCA[idCCA],orderRequest)
        
        def isRequestPrioritaryCCA(self,r,idCCA):
            return r==max(self.requestsCCA[idCCA]).getId()
        
        def isRequestPrioritary(self,r):
            for idCCA in self.ccaRequests[r]:
                if not self.isRequestPrioritaryCCA(r,idCCA):
                    return False
            return True
            
        def extractActivitiesFromMode(self,constellation,r,m):
            act = []
            for (s,o,d) in constellation.getRequest(r).getMode(m).getPairs():
                if o not in act:
                    act.append(o)
                if d not in act:
                    act.append(d)
            return act        
               
        def registerExplaination(self,explaination):
            for ex in explaination:
                if explaination not in self.explainations:
                    self.explainations.append(ex)
                
        def cleanUpSolution(self,constellation):
            for mode in self.researchState.getWorkingModes():
                self.cancelModeSolution(constellation,mode)
           
        def initIteration(self,constellation,noise):
            self.resetNoise(constellation,noise)
            self.researchState = EtatRecherche()
            # ajouter modes actifs
            self.explaination = {}
            self.restart(constellation)
            self.initRequests(constellation,noise,conserveOld=False)
            updateModes = []
            for r in constellation.getRequests():
                m = constellation.getRequest(r).getCurrentMode(constellation).getId()
                updateModes.append((r,m))
            for (r,m) in updateModes:
                self.addActiveMode(constellation,r,m) # etat recherche + CCA sol
            # stoquage etat des requests dans les CCA
            self.updateRequestsAndCCAMappings(constellation) # cca => requests,req => ccas 
            
        def resolve(self,constellation,mailbox):
            timeVar = time()
            comm = MPI.COMM_WORLD
            rank = comm.Get_rank()
            size = comm.Get_size()
            timeVar = time()
            if size==1:
                freeCPUs = [0] # un seul coeur => slave local
            else:
                freeCPUs = list(range(1,size))
            lastRecord = 0
            superIter = 0
            while time()-self.startDate<config.getOptValue("time"):
                sigma = config.getOptValue("noise")*superIter
                self.initIteration(constellation,sigma)
                while time()-self.startDate<config.getOptValue("time") and self.researchState.getModesInScheduling():
                    start_it = time()
                    self.sendTask(constellation,freeCPUs,mailbox,superIter)
                    self.analyzeResults(constellation,freeCPUs,mailbox,superIter)

                    step()
                    if config.verifMode():
                        printColor("working modes:",self.researchState.getWorkingModes(),c='r')
                        printColor("modes waiting for validation:",self.researchState.getModesWaitingForValidation(),c='y')
                    
                    if config.getOptValue("verif"):
                        self.verifyCCAs(constellation)
                    """
                        stockage de la solution courante
                    """
                    if lastRecord<( (time()-self.startDate)//10):
                        lastRecord = (time()-self.startDate)//10
                        self.notifyNoEvent(constellation)
                        self.displayInformation(time(),self.startDate,constellation,title="INFORMATION")
                self.notifyEndIteration(constellation)
                self.displayInformation(time(),self.startDate,constellation,title="END ITERATION "+str(superIter),color='r')
                superIter += 1
                timeVar = time()
                if config.getOptValue("noise") == 0:
                    break
            
            self.cleanUpSolution(constellation)
            self.terminateProcesses()
            self.notifyEndExecution(constellation)
            self.displayInformation(time(),self.startDate,constellation,title='FIN',color='y')
        
        def terminateProcesses(self):
            comm = MPI.COMM_WORLD
            rank = comm.Get_rank()
            size = comm.Get_size()
            for i in range(size):
                if i!=0:
                    comm.send({'fin':True},dest=i)
        
        def updateRequestsAndCCAMappings(self,constellation):
            self.requestsCCA,self.ccaRequests = {},{}
            req,ccas = {},{}
            for r in constellation.getRequests():
                req[r] = []
                m = self.currentModes[r]
                orderRequest = RequestOrdering(r,m,self.noiseUtility(constellation,r,m))
                for (s,o) in constellation.getRequest(r).getPairs():
                    idCCA = self.getGraphComponents().getActivityCCA(o)
                    if idCCA not in ccas:
                        ccas[idCCA] = []
                    if idCCA not in req[r]:
                        req[r].append(idCCA)
                        reverse_insort(ccas[idCCA],orderRequest)
            self.requestsCCA = ccas
            self.ccaRequests = req
        
        def getFreeRequetsAndCCAsPairs(self):
            requests = {} # cca ou il y a du travail
            freeCCAss = []
            for idCCA in self.requestsCCA:
                rec_max = max(self.requestsCCA[idCCA]).getUtility()
                r_max = [r.getId() for r in self.requestsCCA[idCCA] if r.getUtility()==rec_max or r.getId()==-1]
                for r in r_max:
                    m = self.currentModes[r]
                    if r not in requests:
                        requests[r] = self.researchState.getCCAAvailableForComputation(r,m,self.getGraphComponents())
                    if not self.researchState.getCCAsUnderComputation(idCCA):
                        if idCCA not in freeCCAss:
                            freeCCAss.append(idCCA)
            return requests,freeCCAss
        
        # renvoie les (r,cca) tels que r est premiere sur sa CCA et (r,cca) n'est pas encore calculée
        def candidatsTop(self,requests,freeCCAss):
            topCandidates = []
            for idCCA in freeCCAss:
                if len(self.requestsCCA[idCCA])>0:
                    rmax = max(self.requestsCCA[idCCA])
                    assert(self.isRequestPrioritaryCCA(rmax.getId(),idCCA))
                    if idCCA in requests[rmax.getId()]:
                        topCandidates.append((0,rmax,idCCA))
            return topCandidates
        
        # les candidats de levelx inferieurs
        def getCandidatesQueue(self,requests,freeCCAss):
            candidatesQueue = []
            for idCCA in freeCCAss:
                candidats = sorted(self.requestsCCA[idCCA],reverse=True)
                for level,req in enumerate(candidats):
                    r = req.getId()
                    m = req.getIdMode()
                    if idCCA in self.researchState.getCCAAvailableForComputation(r,m,self.getGraphComponents()):
                        candidatesQueue.append((level,req,idCCA))
                        break
            return candidatesQueue
        
        def pickCandidateInQueue(self,candidatesQueue,constellation):
            if len(candidatesQueue)>0:
                res = max(candidatesQueue,key=lambda requestInfo:(-requestInfo[0],requestInfo[1])) # min level,max prio,max rec
                return res
            else:
                return None
            
        def pickTopCandidate(self,topCandidates,constellation):
            if len(topCandidates)>0: # candidats top = liste de (level,(r,m),cca)
                res = max(topCandidates,key=itemgetter(1)) # item[1] => la classe de relation d'ordre
                return res
            else:
                return None
            
        def pickMode(self,constellation):
            start = time()
            # deteminer les ccas libres et les cca libres sur lesquelles les requests ont du travail
            requests,freeCCAs = self.getFreeRequetsAndCCAsPairs()
            #printColor(requests,freeCCAs,c='b')
            start = time()
            topCandidates = self.candidatsTop(requests,freeCCAs)
            #printColor(topCandidates,c='r')
            #printColor(self.researchState.CCAsUnderComputation,c='m')
            res = self.pickTopCandidate(topCandidates,constellation)
               
            if res is not None:
                return res
            # si pas de candidat et mode "anticiper" activé : on cherche les candidats dans la queue
            if not config.getOptValue("anticipation"):
                return None
            else:
                candidatesQueue = self.getCandidatesQueue(requests,freeCCAs)
                return self.pickCandidateInQueue(candidatesQueue,constellation)
        
        def existsAvailableProcesses(self,freeCPUs):
            return freeCPUs!=[]
        
        def pickAvailableProcess(self,freeCPUs):
            return freeCPUs.pop()
        
        def sendTask(self,constellation,freeCPUs,mailbox,iteration):
            while self.existsAvailableProcesses(freeCPUs) and self.researchState.getModesInScheduling():
                res = self.pickMode(constellation)
                if res is not None:
                    level,RequestOrdering,idCCA = res
                    mode = (RequestOrdering.getId(),RequestOrdering.getIdMode())
                    cpu = self.pickAvailableProcess(freeCPUs)
                    activites = self.researchState.getCCA(mode[0],mode[1],self.getGraphComponents(),idCCA)
                    if(len(activites)==0):
                        die(mode[0],mode[1],self.getGraphComponents(),idCCA)
                    assert(len(activites)>0)
                    idCCA = self.getGraphComponents().getActivityCCA(activites[0])
                    cca = self.getSolCCA(idCCA[0],idCCA[1])
                    messageLevel = config.getOptValue("anticipation")*("rank: "+str(level))
                    printColor("picked mode:",str(mode),"on CCA",str(idCCA)+". Activities: "+str(activites),messageLevel,c='b')
                
                    for a in activites:
                        if a in cca.sequence:
                            printColor(a,cca,c='r')
                            assert(False)
                    # demander au slave d'inserer les activites dans la cca
                    if cpu==0:
                        data = {'iteration':iteration,'mode':mode,'cca':deepcopy(cca),'source':0,'activites':activites,'time':time()}
                        self.localSlave.insertActivities(constellation,data,self.transitionModel)
                    else:
                        mailbox.askForSchedulingTask(mode,activites,cca,cpu,iteration)
                else:
                    break
        
        def cancelModeSolution(self,constellation,mode):
            r,m = mode
            cca = []
            for (s,o) in constellation.getRequest(r).getMode(m).getPairs():
                idCCA = self.getGraphComponents().getActivityCCA(o)
                if idCCA not in cca:
                    cca.append(idCCA)
            for (s,cca_a) in cca:
                self.getSolCCA(s,cca_a).cancelMode(constellation,r,m,self.transitionModel)
            
        # si removal est None on retire toutes les activites du mode. Sinon on retire la liste 'removal'
        def removeModeFromSolution(self,constellation,mode,removal=None):
            r,m = mode
            cca = {}                     
            for (s,o) in constellation.getRequest(r).getMode(m).getPairs():
                if removal is None or o in removal:
                    idCCA = self.getGraphComponents().getActivityCCA(o)
                    if idCCA not in cca:
                        cca[idCCA] = []
                    if o in self.getSolCCA(idCCA[0],idCCA[1]).getSequence():
                        cca[idCCA].append(o)
            for (s,cca_a) in cca:
                self.getSolCCA(s,cca_a).removeActivityList(constellation,cca[(s,cca_a)],self.transitionModel)
                                   
        # active le mode pending suivant s'il y en a un et rajoute un pending mode s'il y en a encore
        def moveOnRequestState(self,constellation,mode,explaination,grapheDep):
            r,m = mode
            res = constellation.getRequest(r).getNextMode(explaination,constellation)
            modesToRestart = []
            if res is not None:
                m = res.getId()
                self.currentModes[r] = m
                # removal : les activities presentes uniquement dans l'ancien mode. Les retirer
                # activitiesToPrevalidate : presentes dans les deux modes et deja validées dans l'ancien => activitiesToPrevalidate dans le nouveau mode
                removal,activitiesToPrevalidate,intersectionWorkingActivities = self.intersectionActivites(constellation,r)              
                #printColor(removal,activitiesToPrevalidate,intersectionWorkingActivities,c='m')
                if self.researchState.getModeState(r,m-1)!=self.researchState.FAILED_STATE:
                    modesToRestart = self.researchState.notifyCertainFailure(constellation,r,m-1,grapheDep,removal)
                else:
                    modesToRestart  = []
                self.removeModeFromSolution(constellation,mode,removal) # retirer de la solution
            else:
                modesToRestart = self.researchState.notifyCertainFailure(constellation,r,m,grapheDep)
                self.removeModeFromSolution(constellation,mode) # retirer de la solution
            if res is not None:
                self.addActiveMode(constellation,r,m,activitiesToPrevalidate) # prevalidation de l'intersection deja validee
                self.researchState.notifyActivitiesPreWorking(r,m,intersectionWorkingActivities)
            else:
                printColor('request',r,'('+str(len(constellation.getRequest(r).getPairs())) +' activities) has nor more mode',c='m')
                self.researchState.notifyRequestInsertFailure(r)
                self.deleteRequest(constellation,r)
            return modesToRestart
                
        def removeRequestFromCCA(self,idCCA,r):
            if idCCA in self.requestsCCA:
                for i,x in enumerate(self.requestsCCA[idCCA]):
                    if x.getId()==r:
                        self.requestsCCA[idCCA].pop(i)
                        if self.requestsCCA[idCCA]==[]:
                            del self.requestsCCA[idCCA]
    
    class Process:
        def __init__(self,role):
            self.role = role
        
        def resolve(self,constellation):
            pass
        
    class Master(Process):
        def __init__(self):
            super().__init__("master")

        def resolve(self,constellation,startDate,mailbox,transitionModel,dt_construction_transition,tlim,CCAs,solution):
            # iterer insertion-réparation
            self.initSolution(constellation,startDate,transitionModel,dt_construction_transition,tlim,CCAs,solution)
            self.solution.resolve(constellation,mailbox)
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
            
        def initSolution(self,constellation,startDate,transitionModel,dt_construction_transition,tlim,CCAs,solution):
            self.solution = UnitParallelCCASearch(constellation,startDate,transitionModel,dt_construction_transition,tlim,CCAs,solution)
            
        def getSolution(self):
            return self.solution.getSolution()
        
        def setSolution(self,sol,vid,modes_retenus,modes_candidats,modes_retires,objectif):
            self.solution.setSol(sol,vid,modes_retenus,modes_candidats,modes_retires,objectif)
        
        def verifySolution(self):
            self.solution.verifySolution()
            
        def getModesSolution(self):
            return self.solution.getModes()
    
    class Slave(Process):
        # classe pour planifier les sequences des cca sur les slaves
        class CCAScheduler:
            def __init__(self,solCCA,transitionModel):
                self.solution = solCCA
                self.transitionModel = transitionModel
    
            def __str__(self):
                return self.solution.__str__()
    
            def getSatellite(self):
                return self.solution.getSatellite()
    
            def getIdentifier(self):
                return self.solution.getIdentifier()
        
            def getSolution(self):
                return self.solution
            
            def getSequence(self):
                return self.solution.getSequence()
            
            # pas vraiment une explaination. Plutot la liste d'ajouts ?
            # /!\ en cas d'échec : garde les activités (pour ça qu'on récupère l'ancienne solution)
            # c'est le master qui gère le removal des activités ratées
            def scheduleMode(self,constellation,mode,activites):
                r,m = mode
                explaination = self.solution.insertActivities(constellation,activites,self.transitionModel)
                previousSolution = deepcopy(self.solution)
                sucess = self.solution.computeNewSequence(constellation,self.transitionModel)
                if not sucess:
                    self.solution = previousSolution
                return explaination
                    
            def isSequenceFeasible(self,constellation):
                return self.solution.isSequenceFeasible(constellation,self.transitionModel)
            
        def __init__(self,constellation,transitionModel):
            comm = MPI.COMM_WORLD
            rank = comm.Get_rank()
            super().__init__("slave "+str(rank))
            activites = constellation.extractActivitiesFromRequests()
            self.dependenciesGraph = StaticCCAsGroup(constellation,activites,transitionModel)
            
        def isEndMessage(self,data):
            return 'fin' in data
        
        def extendComputationDuration(self):
            durees = []
            for i,date in enumerate(self.dates_cca):
                if i<len(self.dates_cca)-1:
                    durees.append((date,self.dates_cca[i+1]-date))
            return durees
        
        def packSolution(self,cca,mode,faisable,explaination,iteration):
            return {'iteration':iteration,'cpu_time':self.extendComputationDuration(),'reception':time(),'envoi':(self.sendingDate,self.sendingDuration),'mode':mode,'cca':deepcopy(cca),'faisable':faisable,'source':rank,'explaination':explaination}
    
        def insertActivities(self,constellation,data,transitionModel):
            self.sendingDuration = time() - data['time']
            self.startComputation = time()
            self.sendingDate = data['time']
            iteration = data['iteration']
            start = time()
            mode = data['mode'] # [(r,m)]
            activites = data['activites']
            solCCA = data['cca']
            # ajouter la sequence comme sequence initiale
            self.planif = self.CCAScheduler(solCCA,transitionModel)
            # solution initiale = sequence (activites des modes deja planifies)
            explaination = self.planif.scheduleMode(constellation,mode,activites)
            solutionCCA = self.planif.getSolution()
            feasible = self.planif.isSequenceFeasible(constellation)
            self.dates_cca = [self.startComputation,time()]
            self.result = self.packSolution(solutionCCA,mode,feasible,explaination,iteration)
                    
        def resolve(self,constellation,mailbox,transitionModel):
            comm = MPI.COMM_WORLD
            while True:
                data = comm.recv()
                if self.isEndMessage(data):
                    break
                else:
                    self.insertActivities(constellation,data,transitionModel)
                    mailbox.postMessage(self.result)
        
        def getResult(self):
            if MPI.COMM_WORLD.Get_size()==1:
                self.result['reception'] = (self.result['reception'],time() - self.result['reception'])
            return self.result
        
        
    path = '../data'
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    
    def sendMessage(mode,activites,cca,cpu,iteration):
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        assert(rank==0)
        data = {'iteration':iteration,'mode':mode,'cca':cca,'source':rank,'activites':activites,'time':time()}
        comm.send(data, dest=cpu, tag=rank)
    
    class UniqueMessageCommunication:
        def __init__(self):
            pass
        
        def postMessage(self,data):
            assert(MPI.COMM_WORLD.Get_rank()>0)
            MPI.COMM_WORLD.send(data,dest=0)
            
        def readMessages(self):
            assert(MPI.COMM_WORLD.Get_rank()==0)
            data = MPI.COMM_WORLD.recv()
            print("test")
            data['reception'] = (data['reception'],time()-data['reception'])
            return [data]
        
        def askForSchedulingTask(self,mode,activites,cca,cpu,iteration):
            assert(MPI.COMM_WORLD.Get_rank()==0)
            sendMessage(mode,activites,cca,cpu,iteration)
            
        def __str__(self):
            return "Unique message communication"
        
    class SynchronizedCommunication:
        def __init__(self):
            comm = MPI.COMM_WORLD 
            self.size = comm.Get_size()-1 
            itemsize = MPI.BOOL.Get_size() 
            if comm.Get_rank() == 0: 
                nbytes = self.size * itemsize 
            else: 
                nbytes = 0
            self.win = MPI.Win.Allocate_shared(nbytes, itemsize, comm=comm) 
            # create a numpy array whose data points to the shared mem
            self.cpuMapped, itemsize = self.win.Shared_query(0)
            for i in range(self.size):
                self.cpuMapped[i] = 0
        
        # appelé par le master
        def askForSchedulingTask(self,mode,activites,cca,cpu,iteration):
            assert(MPI.COMM_WORLD.Get_rank()==0)
            self.cpuMapped[cpu-1] = 1
            sendMessage(mode,activites,cca,cpu,iteration)
            
        # appelé par les slaves
        def postMessage(self,data):
            assert(MPI.COMM_WORLD.Get_rank()>0)
            MPI.COMM_WORLD.send(data,dest=0)
            rank = MPI.COMM_WORLD.Get_rank()
            #self.cpuMapped[rank-1] = 1
            
        def readMessages(self):
            assert(MPI.COMM_WORLD.Get_rank()==0)
            nCPUMapped = sum([self.cpuMapped[i] for i in range(self.size)])
            for i in range(nCPUMapped):
                data = MPI.COMM_WORLD.recv()
                data['reception'] = (data['reception'],time()-data['reception'])
                source = data["source"]
                self.cpuMapped[source-1] = 0
                yield data
            else:
                return None
        
        def __str__(self):
            return "Synchronized message"
        
    class SharedCommunication:
        def __init__(self):
            comm = MPI.COMM_WORLD 
            self.size = comm.Get_size()-1 
            itemsize = MPI.BOOL.Get_size() 
            if comm.Get_rank() == 0: 
                nbytes = self.size * itemsize 
            else: 
                nbytes = 0
            self.win = MPI.Win.Allocate_shared(nbytes, itemsize, comm=comm) 
            # create a numpy array whose data points to the shared mem
            self.flagMessages, itemsize = self.win.Shared_query(0)
            for i in range(self.size):
                self.flagMessages[i] = False
        
        def askForSchedulingTask(self,mode,activites,cca,cpu,iteration):
            assert(MPI.COMM_WORLD.Get_rank()==0)
            sendMessage(mode,activites,cca,cpu,iteration)
    
        def postMessage(self,data):
            assert(MPI.COMM_WORLD.Get_rank()>0)
            slot = MPI.COMM_WORLD.Get_rank()-1
            self.flagMessages[slot] = True
            MPI.COMM_WORLD.send(data,dest=0)
        
        def readMessages(self):
            assert(MPI.COMM_WORLD.Get_rank()==0)
            for i,x in enumerate(self.flagMessages):
                if x:
                    process = i+1
                    data = MPI.COMM_WORLD.recv(source=process)
                    data['reception'] = (data['reception'],time()-data['reception'])
                    assert(process==data['source'])
                    self.flagMessages[i] = False
                    yield data
                else:
                    yield None
        
        def __str__(self):
            return "Shared communication: "+ str([self.flagMessages[i] for i in range(self.size)])
    
    # Création de l'objet GLOBAL mailbox
    def createMailbox():
        comm = MPI.COMM_WORLD 
        if config.glob.sync:
            mailbox = SynchronizedCommunication()
        else:
            mailbox = SharedCommunication()
            #mailbox = UniqueMessageCommunication()
        return mailbox   

    
    class runnableUPCCAS:
        def execute(self,constellation,startDate,transitionModel,dt_construction_transition,tlim=np.inf,CCAs=None,solution=None):
            mailbox = SharedCommunication()
            if rank == 0:
                self.process = Master()
                self.process.resolve(constellation,startDate,mailbox,transitionModel,dt_construction_transition,tlim,CCAs,solution)
                
                for i in range(1,MPI.COMM_WORLD.Get_size()):
                    MPI.COMM_WORLD.send({"sol":self.process.solution.getSolutionContainer()},dest=i)
                comm.Barrier()
                return self.process.solution.getSolutionContainer()
            else:
                self.process = Slave(constellation,transitionModel)
                self.process.resolve(constellation,mailbox,transitionModel)
                
                data = None
                data = MPI.COMM_WORLD.recv(data)
                comm.Barrier()
                return data['sol']                
                    
        def getName(self):
            return "Unit Parallel CCA Search"
            
        def getMasterSolver(self):
            if MPI.COMM_WORLD.Get_rank()==0:
                return self.process.solution
            else:
                return None
        
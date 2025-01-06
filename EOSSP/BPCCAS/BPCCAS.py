from EOSSP.Utils.config import *
from EOSSP.Utils.Utils import *
global config
config = Config()
instance = config.instance
from mpi4py import MPI
if config.getOptValue("help"):
    if MPI.COMM_WORLD.Get_rank()==0:
        config.displayHelp()
else:    
    from EOSSP.model.solution import *
    from EOSSP.model.constellation import *
    from EOSSP.model.components import *
    from EOSSP.model.componentPlan import *
    
    from time import time
    from time import sleep
    import random as rd
    
    class BatchParallelCCASearch(Solver):
        def __init__(self,constellation,start,transitionModel,dt_construction_transition,colors={'iteration':'c',"no_event":'w','end':'y'},tlim=np.inf,CCAs=None,solution=None):
            super().__init__(constellation,start,transitionModel,dt_construction_transition,tlim=tlim,CCAs=CCAs,solution=solution)
            self.colors = colors
            self.explainations = {}
            if MPI.COMM_WORLD.Get_size()==1: # structure de données des slaves si processus seul
                self.failuresWithoutExplainations = []
                self.undetermined = []
                self.sequences = {}
            
            comm = MPI.COMM_WORLD
            if comm.Get_size()==1:
                self.localSlave = Slave(constellation,start,transitionModel)
            printColor("Initialization duration:",time()-self.startDate,c='b')
            # mapper les cca
            #self.creerMappingCCASlaves()
        
        def isActiveDestructionMode(self):
            return config.getOptValue("destroy_BPCCAS") is not None and config.getOptValue("destroy_BPCCAS")>0
            
        def isActiveNoisingMode(self):
            return config.getOptValue("noise")>0
    
        def restart(self,constellation,sigma):
            self.solution.restart(constellation)
            for (s,cca) in self.getGraphComponents().getComponents():
                self.setSolCCA(s,cca,SolCCA(cca,s))
        
        def initIteration(self,constellation,sigma,superIter):
            self.initRequests(constellation,sigma,conserveOld=False)
            if superIter>0 or self.isActiveNoisingMode():
                self.restart(constellation,sigma)
            if superIter!=0 and self.isActiveDestructionMode():
                self.destroy(constellation)
            self.notifyNoEvent(constellation)
                
                
        def destroy(self,constellation):
            printOpen("Destruction",c='m')
            if config.getOptValue("verif"):
                self.verifySolution(constellation)        
            
            fulfilled = self.fulfilledRequests()
            destroy = rd.choices(fulfilled,k=int(len(fulfilled)*config.getOptValue("destroy_BPCCAS")))
            if -1 in destroy:
                destroy.remove(-1)
            modesToRemove = [(r,m) for (r,m) in self.getSelectedModes() if r in destroy and r!=-1]
            self.setSelectedModes([(r,m) for (r,m) in self.getSelectedModes() if (r, m) not in modesToRemove],constellation)
            activitiesToRemove = [(s,a) for (r,m) in modesToRemove for (s,a) in constellation.getRequest(r).getMode(m).getPairs()]
            activitiesByCCA = {}
            for (s,a) in activitiesToRemove:
                idCCA = self.getGraphComponents().getActivityCCA(a)
                if idCCA not in activitiesByCCA:
                    activitiesByCCA[idCCA] = []
                activitiesByCCA[idCCA].append(a)
            for (s,cca) in activitiesByCCA:
                self.getSolCCA(s,cca).removeActivityList(constellation,activitiesByCCA[cca])
            if config.getOptValue("verif"):
                self.verifySolution(constellation)
            self.updateDestruction(constellation)
            printColor(str(len(destroy))+" destroyed requests",c='y')
            printClose()
            
        def updateDestruction(self,constellation):
            fulfilled = self.fulfilledRequests()
            for r in constellation.getRequests():
                if r not in fulfilled:
                    constellation.getRequest(r).init = False
                    m = constellation.getRequest(r).getCurrentMode(constellation).getId()
                    self.currentModes[r] = m
                    assert(m is not None)
            # score candidats modes avant la planif
            modes = []
            for r in self.candidateRequests:
                modes.append((r,0))
            self.resetScoreDestruction(constellation,fulfilled)
            self.solution.history.meanScoreBefore = constellation.meanObservationScore(modes)
            self.solution.history.meanObservationScore = self.solution.history.meanScoreBefore
            
            
        def creerMappingCCASlaves(self):
            if MPI.COMM_WORLD.Get_size()>1:
                self.ccaSlaves = {cpu : [] for cpu in range(1,MPI.COMM_WORLD.Get_size())}
                sizes = []
                for idCCA in self.getGraphComponents().getComponents():
                    size =  self.getGraphComponents().getComponentSize(idCCA)
                    reverse_insort(sizes,(size,idCCA))
                cpu = 0
                for (size,idCCA) in sizes:
                    self.ccaSlaves[cpu+1].append(idCCA)
                    cpu = (cpu + 1) % (MPI.COMM_WORLD.Get_size()-1)
                for cpu in self.ccaSlaves:
                    MPI.COMM_WORLD.send({"cca":self.ccaSlaves[cpu]},dest=cpu)
            else:
                self.ccaSlaves = {1 : self.getGraphComponents().getComponents()}
         
        def mappingCCASize(self,sizes):
            # sizes : liste (cca,nb activites a inserer)
            ccas = sorted(sizes,key=itemgetter(0),reverse=True)
            if MPI.COMM_WORLD.Get_size()==1:
                self.ccaSlaves = {0:[]}
                for (taille,cca) in ccas:
                    self.ccaSlaves[0].append(cca)
            else:
                self.ccaSlaves = {cpu : [] for cpu in range(1,MPI.COMM_WORLD.Get_size())}
                cpu = 0
                for (taille,idCCA) in ccas:
                    self.ccaSlaves[cpu+1].append(idCCA)
                    cpu = (cpu + 1) % (MPI.COMM_WORLD.Get_size()-1)
    
        """
            =============================================== 
                            RESOLUTION
            =============================================== 
        """
        def reinsertMode(self,constellation,r):
            pass
    
        def analyzeMessage(self,data,constellation,freeCPUs,mailbox):
            # les sequences solutions ne sont pas transmises ici car les cca sont statiques. 
            # Les slaves ont deja un exemplaire
            start = time()
            self.solution.history.registerMessageReceiving(data)
            i = data['source']
            #del self.jobSlaves[i]['time']
            assert(i not in freeCPUs)
            bisect.insort(freeCPUs,i)
            cpu_time = data['cpu_time']
            explainations = data['explainations']
            undetermined = data['undetermined']
            return explainations,undetermined
        
        def removeAncientActivities(self,constellation,r,m):
            for (s,o) in constellation.getRequest(r).getMode(m-1).getPairs():
                if (s,o) not in constellation.getRequest(r).getMode(m).getPairs():
                    (s,cca) = self.getGraphComponents().getActivityCCA(o)
                    if o in self.getSolCCA(s,cca).getSequence():
                        self.getSolCCA(s,cca).removeActivity(constellation,o)
        
        def analyzeAsynchronousResults(self,constellation,freeCPUs,mailbox):
            failed = list(self.explainations.keys())+self.failuresWithoutExplainations
            succeededModes = []
            for i in self.jobSlaves:
                for cca in self.jobSlaves[i]:
                    for (r,m,score) in self.jobSlaves[i][cca]:
                        if (r,m) not in succeededModes and (r,m) not in failed and (r,m) not in self.undetermined:
                            succeededModes.append((r,m))
            for (r,m) in self.explainations:
                if len(self.explainations[(r,m)])>0:
                    m_new = constellation.getRequest(r).getNextMode(self.explainations[(r,m)],constellation)
                    if m_new is None:
                        self.candidateRequests.remove(r)
                    else:
                        self.removeAncientActivities(constellation,r,m_new.getId())
            
            for rm in succeededModes:
                if rm not in self.getSelectedModes():
                    self.addSelectedMode(rm,constellation)
            self.insertSequences(constellation,self.sequences)
            
            printOpen("Results",c='b')
            printOpen("modes that failed with explaination:",len(list(self.explainations.keys())),c='r')
            printColor((list(self.explainations.keys())),c='r')
            printClose()
            printOpen("modes that failed without explaination:",len(self.failuresWithoutExplainations),c='m')
            printColor(self.failuresWithoutExplainations,c='m')
            printClose()
            printOpen("modes that succeeded:",len(succeededModes),c='g')
            printColor((succeededModes),c='g')
            printClose()
            printOpen("undetermined modes:",len(self.undetermined),c='w')
            printColor((self.undetermined),c='w')
            printClose()
            printClose()
            return failed
        
        def analyzeResults(self,constellation,freeCPUs,mailbox):
            comm = MPI.COMM_WORLD
            rank = comm.Get_rank()
            start_iteration = time()
            if comm.Get_size()==1:
                messages = [self.localSlave.getResults()]
            else:
                messages = mailbox.readMessages()
            
            explainations = {}
            sequences = {}
            undetermined = []
            for data in messages:
                if data is not None:
                    for cca in data['sequences']:
                        sequences[cca] = data['sequences'][cca]
                    expCPU,inde = self.analyzeMessage(data,constellation,freeCPUs,mailbox)
                    for (r,m) in inde:
                        if (r,m) not in undetermined:
                            undetermined.append((r,m))
                    for (r,m) in expCPU:
                        if (r,m) not in explainations:
                            explainations[(r,m)] = []
                        explainations[(r,m)] += expCPU[(r,m)]
                   
            succeededModes = []
            for i in self.jobSlaves:
                for cca in self.jobSlaves[i]:
                    for (r,m,score) in self.jobSlaves[i][cca]:
                        if (r,m) not in succeededModes and (r,m) not in explainations and (r,m) not in undetermined:
                            succeededModes.append((r,m))
            
            for (r,m) in explainations:
                m_new = constellation.getRequest(r).getNextMode(explainations[(r,m)],constellation)
                if m_new is None:
                    self.requestCandidates.remove(r)
            
            for rm in succeededModes:
                if rm not in self.getSelectedModes():
                    self.addSelectedMode(rm,constellation)
            self.insertSequences(constellation,sequences)
            
            printOpen("Result",c='b')
            printColor("modes that failed:",list(explainations.keys()),c='r')
            printColor("modes that succeeded:",succeededModes,c='g')
            printColor("undetermined modes:",undetermined,c='c')
            printClose()
            return list(explainations.keys())
        
        def insertSequences(self,constellation,sequences):
            requests_retenues = [x[0] for x in self.getSelectedModes()]
            for (s,cca) in sequences:
                seq = [a for a in sequences[(s,cca)] if constellation.getRequestActivity(a) in requests_retenues]
                self.getSolCCA(s,cca).setSequence(constellation,seq,self.transitionModel)
    
        def resolve(self,constellation,mailbox):
            comm = MPI.COMM_WORLD
            rank = comm.Get_rank()
            size = comm.Get_size()
            #objectif = self.calculerObjectif(constellation)
            if size==1:
                freeCPUs = [0] # un seul coeur => slave local
            else:
                freeCPUs = list(range(1,size))
            lastRecord = 0
            failed = []
            superIter = 0
            it = 0
            while time()-self.startDate<self.tlim:
                sigma = config.getOptValue("noise")*superIter
                self.initIteration(constellation,sigma,superIter)
                if (self.isActiveNoisingMode() or self.isActiveDestructionMode()):
                    couleur = self.colors['iteration']
                    self.displayInformation(time(),self.startDate,constellation,color=couleur,add={"sigma":sigma},title='START ITERATION '+str(superIter))
                while time()-self.startDate<config.getOptValue("time") and (len(failed)!=0 or it==0):
                    start_it = time()
                    self.sendAsynchronousTask(constellation,freeCPUs,mailbox,failed)
                    failed = self.analyzeAsynchronousResults(constellation,freeCPUs,mailbox)
                    if time()-self.startDate<config.getOptValue("time") and (len(failed)!=0 or it==0):
                        self.notifyNoEvent(constellation)
                        couleur = self.colors['no_event']
                        self.displayInformation(time(),self.startDate,constellation,color=couleur,title="RESULT SUB-ITERATION " +str(it))
                    else:
                        break
                    it += 1
                    temps = time()
                    if config.getOptValue("verif"):
                        self.verifyCCAs(constellation)
                    step()
                if config.getOptValue("verif"):
                    self.verifySolution(constellation)
                self.notifyEndIteration(constellation)
                superIter += 1 
                it = 0
                if (not self.isActiveNoisingMode() and not self.isActiveDestructionMode()):
                    break
            
            self.terminateProcess()
            self.notifyEndExecution(constellation)
            couleur = self.colors['end']
            self.displayInformation(time(),self.startDate,constellation,color=couleur,title='END')
    
        def mappingCCACurrentObservations(self,constellation,failed):
            mappingObs = {}
            
            del_r = []
            for r in self.candidateRequests:
                if r in self.fulfilledRequests():
                    del_r.append(r)
            for r in del_r:
                self.candidateRequests.remove(r)
                
            for r in self.candidateRequests:
                m = constellation.getRequest(r).getCurrentMode(constellation).getId()
                score = self.noiseUtility(constellation,r,m)
                for (s,o) in constellation.getRequest(r).getMode(m).getPairs():
                    idCCA = self.getGraphComponents().getActivityCCA(o)
                    if idCCA not in mappingObs:
                        mappingObs[idCCA] = {}
                    if (r,m,score) not in mappingObs[idCCA]:
                        mappingObs[idCCA][(r,m,score)] = []
                    mappingObs[idCCA][(r,m,score)].append(o)
            mappingObs,ccaSequences = self.filterFailedObservations(constellation,mappingObs,failed)
            return mappingObs,ccaSequences
        
        def newRequest(self,constellation,idCCA,r,m):
            testExplaination = r in [x[0] for x in self.explainations.keys()]
            previousAbsentInCCA = True
            if m>0:
                for (s,o) in constellation.getRequest(r).getMode(m-1).getPairs():
                    if self.getGraphComponents().getActivityCCA(o)==idCCA:
                        previousAbsentInCCA = False
                        break
            return testExplaination and previousAbsentInCCA
        
        def rejectMode(self,keep,idCCA,r,m,score,mappingObs,modesToDelete,ccaSequences):
            rejectedModes = 0
            rejectedActivities = 0
            if (r,m,score) in mappingObs[idCCA]:
                for o in mappingObs[idCCA][(r,m,score)]:
                    if o in ccaSequences[idCCA]: # on est pas sur (modes rates, modes undetermined)
                        ccaSequences[idCCA].remove(o) 
            if (r,m) in self.getSelectedModes() and (r,m) not in modesToDelete:
                modesToDelete.append((r,m))
                rejectedModes += 1
                rejectedActivities += len(mappingObs[idCCA][(r,m,score)])
            return rejectedModes,rejectedActivities
        
        def keepMode(self,r,m,score,mappingObs,ccaSequences,idCCA,keptCCAs):
            keptModes = 0
            keptActivities = 0
            if (r,m) in self.getSelectedModes():
                keptModes += 1
                keptActivities += len(mappingObs[idCCA][(r,m,score)])
                del mappingObs[idCCA][(r,m,score)]
            else:
                for o in mappingObs[idCCA][(r,m,score)]:
                    if o in ccaSequences[idCCA]: # modes undetermined (peut etre partiellement satisfaits)
                        mappingObs[idCCA][(r,m,score)].remove(o)
                        if mappingObs[idCCA][(r,m,score)]==[]:
                            del mappingObs[idCCA][(r,m,score)]
            if mappingObs[idCCA] == {}:
                idCCA = idCCA
                keptCCAs.append(idCCA)
            return keptModes,keptActivities
            
        def filterFailedObservations(self,constellation,mappingObs,failed):
                ccaSequences = {}
                modesToDelete = []
                keptModes = 0
                keptActivities = 0
                keptCCAs = []
                rejectedModes = 0
                rejectedActivities = 0
                
                for idCCA in mappingObs:
                    tri = list(mappingObs[idCCA].keys())
                    modes = sorted(tri,key=lambda x : (x[0]==-1,x[2][0],-x[0]),reverse=True)
                    keep = True
                    candidatesModes = [] # ceux qu'on conserve pas besoin de redemander de les inserer
                    s,cca = idCCA
                    ccaSequences[idCCA] = self.getSolCCA(s,cca).getSequence().copy()
                    for (r,m,score) in modes:
                        if keep:
                            res = self.keepMode(r,m,score,mappingObs,ccaSequences,idCCA,keptCCAs)
                            keptModes += res[0]
                            keptActivities += res[1]
                        else:
                            res = self.rejectMode(keep,idCCA,r,m,score,mappingObs,modesToDelete,ccaSequences)
                            rejectedModes += res[0]
                            rejectedActivities += res[1]
                        # si une nouvelle requete arrive l'ordre change potentiellement
                        if not config.getOptValue("conserve") and ((r,m) in self.explainations or self.newRequest(constellation,idCCA,r,m)):
                            keep = False
                            
                for idCCA in keptCCAs:
                    del mappingObs[idCCA]
                for x in modesToDelete:
                    self.removeSelectedMode(x)
                #printColor("Recycling "+str(keptActivities)+" activities "+str(keptModes)+" partial modes and "+str(len(keptCCAs))+" CCAs",c='y')
                #printColor("Rejecting "+str(rejectedActivities)+" activities "+str(rejectedModes)+" partial modes",c='y')
                return mappingObs,ccaSequences
            
        def dispatchWorkToSlaves(self,mappingObs,ccaSequences):
            sequencesSlave = {}
            sizes = [(sum([len(mappingObs[idCCA][x]) for x in mappingObs[idCCA]]),idCCA) for idCCA in mappingObs]
            self.mappingCCASize(sizes)
            jobSlaves = {}
            for cpu in self.ccaSlaves:
                sequencesSlave[cpu] = {}
                jobSlaves[cpu] = {}
                for idCCA in self.ccaSlaves[cpu]:
                    if idCCA in mappingObs:
                        sequencesSlave[cpu][idCCA] = ccaSequences[idCCA]
                        jobSlaves[cpu][idCCA] = mappingObs[idCCA]
            jobsToDelete = []
            for cpu in jobSlaves:
                if jobSlaves[cpu]=={}:
                    jobsToDelete.append(cpu)
            for cpu in jobsToDelete:
                del jobSlaves[cpu]
            self.jobSlaves = jobSlaves
            return jobSlaves,sequencesSlave
    
        def asynchronousReading(self,constellation,freeCPUs,mailbox,failures,undetermined,sequences):
            count = 0
            for data in mailbox.readMessages():
                if data is not None:
                    count += 1
                    for cca in data['sequences']:
                        self.sequences[cca] = data['sequences'][cca]
                    expCPU,inde = self.analyzeMessage(data,constellation,freeCPUs,mailbox)
                    for (r,m) in inde:
                        if (r,m) not in undetermined:
                            undetermined.append((r,m))
                    for (r,m) in expCPU:
                        if (r,m) not in failures:
                            failures[(r,m)] = []
                        failures[(r,m)] += expCPU[(r,m)]
            return count
        
        def asynchronousSend(self,ccas,mappingObs,ccaSequences,freeCPUs,mailbox,jobSlaves):
            count = 0
            while ccas!=[] and freeCPUs!=[]:
                (taille,idCCA) = ccas.pop(0)
                if idCCA in mappingObs:
                    cpu = freeCPUs.pop()
                    if cpu not in jobSlaves:
                        jobSlaves[cpu] = {}
                    jobSlaves[cpu][idCCA] = mappingObs[idCCA]
                    data = {"taches":{idCCA:jobSlaves[cpu][idCCA]},"sequences":{idCCA:ccaSequences[idCCA]}}
                    mailbox.askForTask(data,cpu)
                    self.solution.history.registerMessageSending(cpu)
                else:
                    count += 1
            return count
        
        def stabilityAnalysis(self,constellation,r,CCAStables,CCAInstables,failures):
            m = constellation.getRequest(r).getCurrentMode(constellation).getId()
            if (r,m) in failures:  # r a échoué
                stable = False
            else:
                stable = True
                for cca in constellation.getRequest(r).getCCAPresentes(self.getGraphComponents()):
                    if cca in CCAInstables:
                        stable = False
                        break
            return stable
        
        def generateRequestExplaination(self,constellation,r,failures,explainations,CCAStables):
            m = constellation.getRequest(r).getCurrentMode(constellation).getId()
            if (r,m) in failures:
                explainations[(r,m)] = []
                for cca in constellation.getRequest(r).getCCAPresentes(self.getGraphComponents()):
                    if cca in CCAStables:
                        ccaFailures = [o for (s,o) in constellation.getRequest(r).getCurrentMode(constellation).getPairs() if self.getGraphComponents().getActivityCCA(o)==cca and o in failures[(r,m)]]
                        explainations[(r,m)] += ccaFailures
                if len(explainations[(r,m)])==0:
                    del explainations[(r,m)]
                    self.failuresWithoutExplainations.append((r,m))
         
        def propagateInstability(self,constellation,r,CCAStables,CCAInstables):
            for cca in constellation.getRequest(r).getCCAPresentes(self.getGraphComponents()):
                if cca in CCAStables:
                    CCAStables.remove(cca)
                    CCAInstables.append(cca)
        
        def generateExplainations(self,constellation,failures):
            if not config.getOptValue("stable_exp"):
                self.failuresWithoutExplainations = []
                self.explainations = failures
            else:
                self.failuresWithoutExplainations = []
                explainations = {}
                requests = [r for r in constellation.getRequests() if constellation.getRequest(r).getCurrentMode(constellation) is not None]
                CCAStables = [idCCA for idCCA in self.getGraphComponents().getComponents()]
                CCAInstables = []
                requests_triees = sorted(requests,key = lambda r : (r==-1,self.noiseUtility(constellation,r,constellation.getRequest(r).getCurrentMode(constellation).getId())[0],-r),reverse=True)
                for r in requests_triees:
                    stable = self.stabilityAnalysis(constellation,r,CCAStables,CCAInstables,failures)
                    self.generateRequestExplaination(constellation,r,failures,explainations,CCAStables)
                    if not stable:
                        self.propagateInstability(constellation,r,CCAStables,CCAInstables)
                self.explainations = explainations
                    
        def simultaneousSendAndRead(self,constellation,mappingObs,ccaSequences,mailbox,freeCPUs):
            t1 = time()
            sizes = [(sum([len(mappingObs[idCCA][x]) for x in mappingObs[idCCA]]),idCCA) for idCCA in mappingObs]
            ccas = sorted(sizes,key=itemgetter(0),reverse=True)
            nCCAs = len(ccas)
            ccas_done = 0
            jobSlaves = {}
            failures = {}
            undetermined = []
            self.sequences = {}
            while ccas_done != nCCAs: # cca done par lecture des messages / cca dones car rien a faire (detection a l'envoi)
                ccas_done += self.asynchronousReading(constellation,freeCPUs,mailbox,failures,undetermined,ccaSequences)
                ccas_done += self.asynchronousSend(ccas,mappingObs,ccaSequences,freeCPUs,mailbox,jobSlaves)
            self.jobSlaves = jobSlaves
            self.generateExplainations(constellation,failures)
            self.undetermined = undetermined
        
        def sendAsynchronousTask(self,constellation,freeCPUs,mailbox,failed):
            mappingObs,ccaSequences = self.mappingCCACurrentObservations(constellation,failed) # cca -> (r,m) -> liste obs sur la cca
            if comm.Get_size()==1:
                freeCPUs.remove(0)
                jobSlaves,sequencesSlave = self.dispatchWorkToSlaves(mappingObs,ccaSequences)
                data = {"time":time(),"taches":jobSlaves[0],"sequences":sequencesSlave[0]}
                self.localSlave.insertActivities(constellation,data)
            else:
                self.simultaneousSendAndRead(constellation,mappingObs,ccaSequences,mailbox,freeCPUs)
            
        def sendTask(self,constellation,freeCPUs,mailbox,failed):
            mappingObs,ccaSequences = self.mappingCCACurrentObservations(constellation,failed) # cca -> (r,m) -> liste obs sur la cca
                # ajouter potentiellement un remapping pour repartir la charge de travail
            jobSlaves,sequencesSlave = self.dispatchWorkToSlaves(mappingObs,ccaSequences)
            #printColor(jobSlaves,c='y')
            comm = MPI.COMM_WORLD
            if comm.Get_size()==1:
                freeCPUs.remove(0)
                data = {"time":time(),"taches":jobSlaves[0],"sequences":sequencesSlave[0]}
                self.localSlave.insertActivities(constellation,data)
            else:
                for cpu in jobSlaves:
                    freeCPUs.remove(cpu)
                    data = {"taches":jobSlaves[cpu],"sequences":sequencesSlave[cpu]}
                    mailbox.askForTask(data,cpu)
                    self.solution.history.registerMessageSending(cpu)
    
    class Processus:
        def __init__(self,role):
            self.role = role
        
        def resolve(self,constellation):
            pass
        
    class Master(Processus):
        def __init__(self,start,tlim):
            super().__init__("master")
            self.start = start
            self.tlim = tlim
            
        def resolve(self,constellation,mailbox,CCAs,solution,transitionModel,dt_construction_transition):
            # iterer insertion-réparation
            self.initSolution(constellation,CCAs,solution,transitionModel,dt_construction_transition)
            self.solution.resolve(constellation,mailbox)
            self.solution.buildPlan(constellation)
    
        def bestSolution(self):
            return self.solution.bestSolution()
        
        def registerObjective(self):
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
        
        def plotActivity(self,constellation,annotate=False):
            return self.solution.plotActivity(constellation,annotate)
            
        def plotHistory(self,init=True):
            return self.solution.plotHistory(init)
        
        def saveSample(self,constellation):
            self.solution.saveSample(constellation)
            
        def initSolution(self,constellation,CCAs,solution,transitionModel,dt_construction_transition):
            self.solution = BatchParallelCCASearch(constellation,self.start,transitionModel,dt_construction_transition,tlim=self.tlim,CCAs=CCAs,solution=solution)
            
        def getSolution(self):
            return self.solution.getSolution()

        def verifySolution(self,transitionModel):
            self.solution.verifySolution(transitionModel)
            
        def getModesSolution(self):
            return self.solution.getModes()
    
    class Slave(Processus):
        # classe pour planifier les sequences des cca sur les slaves
        class CCAPlanner:
            def __init__(self,sol_cca,transitionModel):
                self.solution = sol_cca
                self.transitionModel = transitionModel
    
            def __str__(self):
                return self.solution.__str__()
            
            def setSequence(self,constellation,sequence):
                self.solution.setSequence(constellation,sequence,self.transitionModel)
                
            def getSatellite(self):
                return self.solution.getSatellite()
    
            def getIdentifier(self):
                return self.solution.getIdentifier()
        
            def getSolution(self):
                return self.solution
            
            def getSequence(self):
                return self.solution.getSequence()
            
            def wasSolverUsed(self):
                return self.solverUsed
            
            # VERSION LKH
            def planMode(self,constellation,activites):
                sequenceBefore = self.solution.sequence.copy()
                previousSolution = deepcopy(self.solution)
                self.solution.insertActivities(constellation,activites,self.transitionModel,cheap=False)
                seqInter = self.solution.sequence.copy()
                #printColor(MPI.COMM_WORLD.Get_rank(),"insertion",activites,len(self.solution.sequence),c='g')
                success = self.solution.computeNewSequence(constellation,self.transitionModel,'LKH')
                self.solverUsed = self.solution.wasSolverUsed()
                if not(not success or len(sequenceBefore)+len(activites)==len(self.solution.sequence)):
                    print(MPI.COMM_WORLD.Get_rank(),success,sequenceBefore,activites,self.solution.sequence)
                    assert(False)
                if not success:
                    sequence = previousSolution.getSequence().copy()
                    self.solution = previousSolution
                return success
            
            # VERSION OPTW
            # (r,m) => liste d'activites
            def planListOfModes(self,constellation,modes_tries):
                shiftRightDisplay(4)
                sequenceBefore = self.solution.sequence.copy()
                previousSolution = deepcopy(self.solution)
                self.solution.insertListOfModes(constellation,modes_tries,self.transitionModel)
                seqInter = self.solution.sequence.copy()
                groups = {}
                groups[0] = copy(sequenceBefore)
                for i,mode in enumerate(modes_tries):
                    groups[i+1] = modes_tries[mode]
                    assert(modes_tries[mode]!=[])
                #printColor(MPI.COMM_WORLD.Get_rank(),"insertion",activites,len(self.solution.sequence),c='g')
                res = self.solution.computeNewSequence(constellation,self.transitionModel,solver=config.getOptValue("solver_PCCAS"),groups=groups)
                success,fulfilledRequests = res
                self.solverUsed = self.solution.wasSolverUsed()
                assert(success)
                shiftLeftDisplay(4)
                return fulfilledRequests
                    
            def sequenceFaisable(self,constellation):
                return self.solution.sequenceFaisable(constellation,self.transitionModel)
            
        def __init__(self,constellation,start,transitionModel,CCAs=None):
            comm = MPI.COMM_WORLD
            rank = comm.Get_rank()
            super().__init__("slave "+str(rank))
            activites = constellation.extractActivitiesFromRequests()
            if CCAs is not None:
                self.dependenciesGraph = CCAs
            else:
                self.dependenciesGraph = StaticCCAsGroup(constellation,activites,transitionModel)
            self.planners = {}
            self.start = start
            
            for s in constellation.getSatellites():
                for cca in self.dependenciesGraph.getCCASatellite(s):
                    self.planners[(s,cca)] = self.CCAPlanner(SolCCA(cca,s),transitionModel)
        
        def isEndMessage(self,data):
            return 'fin' in data
        
        def extendComputationTime(self):
            durees = []
            for i,date in enumerate(self.ccaDates):
                if i<len(self.ccaDates)-1:
                    durees.append((date,self.ccaDates[i+1]-date))
            return durees
        
        def packSolution(self):
            return {'undetermined':self.undeterminedStates,'sequences':self.sequences,'explainations':self.explainations,'cpu_time':self.extendComputationTime(),'reception':time(),'envoi':(self.sendingDate,self.sendingDuration),'source':rank,'solver time':self.solverTime}
        
        
        def insertActivitiesOPTW(self,constellation,data):
            self.sendingDuration = time() - data['time']
            self.startComputation = time()
            self.sendingDate = data['time']
            del data['time']
            start = time()
            for idCCA in self.planners: 
                self.planners[idCCA].setSequence(constellation,data['sequences'].get(idCCA,[]))
            self.sequences = {}
            self.explainations = {}
            self.undeterminedStates = []
            self.ccaDates = [time()]
            self.solverTime = 0
            for idCCA in data['taches']:
                modes = sorted(data['taches'][idCCA],key=lambda x : (x[0]==-1,x[2][0],-x[0]),reverse=True)
                activites_modes_tries = {(x[0],x[1]) : data['taches'][idCCA][x] for x in modes}
                startSolver = time()
                fulfilledRequests = self.planners[idCCA].planListOfModes(constellation,activites_modes_tries)
                if self.planners[idCCA].wasSolverUsed():
                        self.solverTime += time()-startSolver
                for (r,m,score) in data['taches'][idCCA]:
                    if r not in fulfilledRequests:
                        activites = data['taches'][idCCA][(r,m,score)]
                        if (r,m) not in self.explainations:
                            self.explainations[(r,m)] = []
                        self.explainations[(r,m)] += activites
                
                self.sequences[idCCA] = self.planners[idCCA].solution.getSequence()
                self.ccaDates.append(time())
                # trier les modes et inserer leurs activites
            self.results = self.packSolution()        
            
        def insertActivities(self,constellation,data):
            self.sendingDuration = time() - data['time']
            self.startComputation = time()
            self.sendingDate = data['time']
            del data['time']
            start = time()
            for cca in self.planners:
                self.planners[cca].setSequence(constellation,data['sequences'].get(cca,[]))
            nActivites = 0
            nModes = 0
            nCCA = 0
            self.sequences = {}
            self.explainations = {}
            self.undeterminedStates = []
            self.ccaDates = [time()]
            self.solverTime = 0
            for cca in data['taches']:
                fails = 0
                n_calls = 0
                stop = False
                if len(list(data['taches'][cca].keys()))>0:
                    nCCA += 1
                for (r,m,score) in sorted(data['taches'][cca].keys(),key=lambda x : (x[0]==-1,x[2][0],-x[0]),reverse=True):
                    if not stop:
                        activites = data['taches'][cca][(r,m,score)]
                        nActivites += len(activites)
                        nModes += 1
                        #printColor(MPI.COMM_WORLD.Get_rank(),cca,(r,m),activites,c='b')
                        startSolver = time()
                        success = self.planners[cca].planMode(constellation,activites)
                        if self.planners[cca].wasSolverUsed():
                            self.solverTime += time()-startSolver
                        if not success:
                            fails += 1
                            n_calls += 1
                            #printColor("CPU ",MPI.COMM_WORLD.Get_rank(),"echec insertion",(r,m),activites,"dans la cca",cca,c='r')
                            if (r,m) not in self.explainations:
                                self.explainations[(r,m)] = []
                            self.explainations[(r,m)] += activites
                            if fails >= config.getOptValue("fails") or n_calls >= config.getOptValue("scalls"):
                                stop = True
                        else:
                            fails = 0
                    else:
                        self.undeterminedStates.append((r,m))
                self.sequences[cca] = self.planners[cca].solution.getSequence()                
                self.ccaDates.append(time())
                # trier les modes et inserer leurs activites
            self.results = self.packSolution()
            #printColor("Processus",MPI.COMM_WORLD.Get_rank(),"terminé (",nActivites,"activités",nModes,"modes",nCCA,"CCAs )",c='b')
        
        def resolve(self,constellation,mailbox):
            comm = MPI.COMM_WORLD
            while True:
                data = comm.recv()
                if self.isEndMessage(data):
                    comm.Barrier()
                    break
                else:
                    if config.getOptValue("solver_PCCAS") == 'LKH':
                        self.insertActivities(constellation,data)
                    else:
                        self.insertActivitiesOPTW(constellation,data)
                    mailbox.postMessage(self.results)
        
        def getResults(self):
            if MPI.COMM_WORLD.Get_size()==1:
                self.results['reception'] = (self.results['reception'],time() - self.results['reception'])
            return self.results
                
    path = '../data'
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    
    
    class Communication:
        def __init__(self):
            # eviter d'envoyer les informations qu'on va recuperer à l'identique ensuite
            self.localData = {}
        
        def sendMessage(self,data,cpu):
            comm = MPI.COMM_WORLD
            rank = comm.Get_rank()
            assert(rank==0)
            data['time']=time()
            comm.send(data, dest=cpu, tag=rank)
        
        def extraireDonneesRedondantes(self,data,cpu):
            # {'mode':mode,'cca':cca,'source':rank,'activites':activites,'time':time()}
            idCCA = data['cca'].getIdentifier()
            self.localData[cpu] = {}
            self.localData[cpu]['sequence'] = data['cca'].getSequence().copy()
            self.localData[cpu]["mode"] = data['mode']
            self.localData[cpu]['activites'] = data['activites']
            self.localData[cpu]['cca'] = idCCA
            del data['mode']
            data['cca'] = idCCA
            
        def reformatData(self,data):
            source = data['source']
            data['cca'] = self.localData[source]['cca']
            data['mode'] = self.localData[source]['mode']
            data['sequence'] = self.localData[source]['sequence']
            data['explication'] = self.localData[source]['activites']
            if data['faisable']:
                data['sequence'] += data['explication']
                
    class UniqueMessageCommunication(Communication):
        def __init__(self):
            pass
        
        def postMessage(self,data):
            assert(MPI.COMM_WORLD.Get_rank()>0)
            MPI.COMM_WORLD.send(data,dest=0)
            
        def readMessages(self):
            assert(MPI.COMM_WORLD.Get_rank()==0)
            data = MPI.COMM_WORLD.recv()
            data['reception'] = (data['reception'],time()-data['reception'])
            return [data]
        
        def askForTask(self,data,cpu):
            assert(MPI.COMM_WORLD.Get_rank()==0)
            sendMessage(data,cpu)
            
        def __str__(self):
            return "Unique message communication"
        
    class SynchronizedCommunication(Communication):
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
            self.cpu_mapped, itemsize = self.win.Shared_query(0)
            for i in range(self.size):
                self.cpu_mapped[i] = 0
        
        # appelé par le master
        def askForTask(self,data,cpu):
            assert(MPI.COMM_WORLD.Get_rank()==0)
            self.cpu_mapped[cpu-1] = 1
            self.sendMessage(data,cpu)
            
        # appelé par les slaves
        def postMessage(self,data):
            assert(MPI.COMM_WORLD.Get_rank()>0)
            MPI.COMM_WORLD.send(data,dest=0)
            rank = MPI.COMM_WORLD.Get_rank()
            #self.cpu_mapped[rank-1] = 1
            
        def readMessages(self):
            assert(MPI.COMM_WORLD.Get_rank()==0)
            nCPUMapped = sum([self.cpu_mapped[i] for i in range(self.size)])
            if nCPUMapped==0:
                return None
            res = []
            for i in range(nCPUMapped):
                data = MPI.COMM_WORLD.recv()
                data['reception'] = (data['reception'],time()-data['reception'])
                source = data["source"]
                self.cpu_mapped[source-1] = 0
                res.append(data)
            return res
        
        def __str__(self):
            return "MSynchronized communication"
            
    class SharedCommunication(Communication):
        def __init__(self):
            super().__init__()
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
            #self.flagMessages = np.ndarray(buffer=self.buf, dtype=bool, shape=(size,))
        
        def askForTask(self,data,cpu):
            assert(MPI.COMM_WORLD.Get_rank()==0)
            self.sendMessage(data,cpu)
    
                
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
    
    class runnableBPCCAS:
        def execute(self,constellation,startDate,transitionModel,dt_construction_transition,tlim=np.inf,CCAs=None,solution=None):
            mailbox = SharedCommunication()
            id_mpi = MPI.COMM_WORLD.Get_rank()
            seed_option = config.getOptValue("seed")
            rd.seed(id_mpi+seed_option)
            if MPI.COMM_WORLD.Get_rank()==0:
                self.process = Master(startDate,tlim)
                self.process.resolve(constellation,mailbox,CCAs,solution,transitionModel,dt_construction_transition)
                
                for i in range(1,MPI.COMM_WORLD.Get_size()):
                    MPI.COMM_WORLD.send({"sol":self.process.solution.getSolutionContainer()},dest=i)
                return self.process.solution.getSolutionContainer()
            else:
                process = Slave(constellation,startDate,transitionModel,CCAs=CCAs)
                process.resolve(constellation,mailbox)
                data = None
                data = MPI.COMM_WORLD.recv(data)
                
                return data['sol']

        def getName(self):
            return "Batch Parallel CCA Search"
        
        def getMasterSolver(self):
            if MPI.COMM_WORLD.Get_rank()==0:
                return self.process.solution
            else:
                return None
        
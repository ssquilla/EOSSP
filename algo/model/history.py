# -*- coding: utf-8 -*-
from datetime import datetime
from time import time
from mpi4py import *
from mpi4py import MPI
import pickle
from copy import deepcopy
import operator
import matplotlib.pyplot as plt
from matplotlib import patches

from .constellation import TYPE_LONG_MONO,TYPE_MONO,TYPE_PERIODIC,TYPE_STEREO,TYPE_SYSTEMATIC,TYPE_DOWNLOAD

from ..Utils.config import *
global config
config = Config()

PREPROCESSING = "PREPROCESSING"
END_ITERATION = "END_ITERATION"
NO_EVENT = "NO_EVENT"
END_RUN = "END_RUN"
BACKUP = "BACKUP"
END_OPERATOR = "END_OPERATOR"

class History:
    def __init__(self,constellation,startDate,dico=None):
        self.startDate = startDate
        self.filtre = None
        # flag d'évenements : indiquer l'état de la recherche
        self.extractRequestsTypesCounts(constellation)
        if dico is None:
            self.historic = [(0,NO_EVENT,(0,0),0,0,MPI.COMM_WORLD.Get_rank(),[],{})] # temps,event,objective,nb modes retenus,mpi rank,ccas,sol
            self.startDatetime = datetime.now().strftime("%m/%d/%Y, %H:%M")
            self.bestObjective = (0,0)
            self.bestSolution = {}
            self.bestSelectedModes = []
            self.resetStats(constellation)
            self.ccas = None
            self.requestsArrival = {}
            instance = constellation.getFilename()
            Ncores = str(MPI.COMM_WORLD.Get_size())
            self.Ncores = Ncores
            self.instance = instance
            self.algoName = config.donnees.algo_name
            self.config = str(config)
            self.meanObservationScore = 0
            self.gap = None
            self.commentaires = []
            self.path = os.getcwd()
            self.additionnalInformation = {}
        else:
            self.requestsArrival = self.extractParameter(dico,'requestsArrival')
            self.ccas = self.extractParameter(dico,'cca')
            self.startDatetime = self.extractParameter(dico,'start')
            self.end_datetime = self.extractParameter(dico,'end')
            self.solution.historic = self.extractParameter(dico,'hist') # temps,event,objective,nb modes retenus
            self.bestObjective = self.extractParameter(dico,'obj')
            self.bestSolution = self.extractParameter(dico,'sol')
            self.bestSelectedModes = self.extractParameter(dico,'modes')
            self.stats = self.extractParameter(dico,'stats')
            self.Ncores = self.extractParameter(dico,'coeurs')
            self.instance = self.extractParameter(dico,'instance')
            self.algoName = self.extractParameter(dico,'algo')
            self.config = self.extractParameter(dico,'config')
            self.meanScoreBefore = self.extractParameter(dico,"score_avant")
            self.meanScoreAfter = self.extractParameter(dico,"score_apres")
            self.meanObservationScore = self.meanScoreAfter
            self.commentaires = self.extractParameter(dico,'commentaires')
            self.path = self.extractParameter(dico,'path')
            self.infosHist = self.extractParameter(dico,'infos hist')
            self.additionnalInformation = self.extractParameter(dico, 'info_add')
            if 'gap' in dico:
                self.gap = dico['gap']
            else:
                self.gap = None
    
    def addInformation(self,cle,valeur):
        self.additionnalInformation[cle] = valeur
    
    def extractRequestsTypesCounts(self,constellation):
        types = {}
        types[TYPE_LONG_MONO] = 0
        types[TYPE_MONO] = 0
        types[TYPE_PERIODIC] = 0
        types[TYPE_STEREO] = 0
        types[TYPE_SYSTEMATIC] = 0
        for r in constellation.getRequests():
            type_r = constellation.getRequest(r).getType()
            if type_r != TYPE_DOWNLOAD:
                if r not in types:
                    types[type_r] = 0
                types[type_r] += 1
        self.types = types
    
    def formatInstance(self):
        res = ""
        for i,t in enumerate(self.types):
            res += str(self.types[t])
            if i<len(list(self.types.keys()))-1:
                res += ','
        return res
    
    def extractParameter(self,dico,cle):
        return dico.get(cle,None)
        
    def merge(self,historics_cpu):
        for cpu in historics_cpu:
            hist = historics_cpu[cpu]
            if hist.bestObjective>self.bestObjective:
                self.bestObjective = hist.bestObjective
                self.bestSolution = hist.bestSolution
                self.modes_retenus = hist.bestSelectedModes
        historics_cpu[MPI.COMM_WORLD.Get_rank()] = self
        new_hist = deepcopy(self.historic)
        for cpu in historics_cpu:
            for x in historics_cpu[cpu].historic:
                if x not in new_hist:
                    new_hist.append(x)
        self.historic = sorted(new_hist,key=operator.itemgetter(0))

                        
    def addComment(self,commentaire):
        self.commentaires.append(commentaire)
        
    def sortObsByRequest(self,constellation):
        self.obsParRequete = {}
        seq = self.extractSequences(self.bestSolution)
        for s in seq:
            for a in seq[s]:
                req = constellation.getRequestActivity(a)
                if req not in self.obsParRequete:
                    self.obsParRequete[req] = []
                self.obsParRequete[req].append(a)
        
    def toDict(self):
        dico = {}
        dico['hist'] = self.historic
        dico['infos hist'] = self.getInfos()
        dico['obj'] = self.bestObjective
        dico['sol'] = self.extractSequences(self.bestSolution)
        dico['obs_par_requetes'] = self.obsParRequete
        dico['modes'] = self.bestSelectedModes
        dico['stats'] = self.stats
        dico['instance'] = self.instance
        dico['coeurs'] = self.Ncores
        dico['algo'] = self.algoName
        dico['config'] = str(self.config)
        dico["score_avant"] = self.meanScoreBefore
        dico["score_apres"] = self.meanScoreAfter
        dico['gap'] = self.gap
        dico['start'] = self.startDatetime
        dico['end'] = self.end_datetime
        dico['commentaires'] = self.commentaires
        dico['path'] = self.path
        dico['cca'] = self.ccas
        dico['arrivee_requetes'] = self.requestsArrival
        dico['info_add'] = self.additionnalInformation
        dico['modes_generes'] = self.generatedModes
        return dico
    
    def registerRequestArrival(self,date,liste_requetes):
        self.requestsArrival[date] = liste_requetes
    
    def getInfos(self):
        mess = " ----------- History of the execution ----------\n"
        lecture = [("chemin",self.path)]
        for i in range(len(self.commentaires)):
            lecture.append(("commentaire "+str(i),self.commentaires[i]))
        lecture += [("algorithme",self.algoName),('instance',self.instance),("nombre de coeurs",self.Ncores)]
        lecture += [("date de démarrage",self.startDatetime),("date de fin",self.end_datetime)]
        #lecture += [("objective",self.bestObjective),("modes retenus",len(self.bestSelectedModes))]
        lecture += [("historique",str(len(self.historic)) + " valeurs enregistrées")]
        for message,contenu in lecture:
            mess += " | "+ message + " : " + str(contenu) +'\n'
        mess += ' '
        mess += self.formatStats() +'\n'
        mess += str(self.config)
        mess += " ------------------------------------------------\n"
        return mess
            
    def lastObjective(self):
        return self.historic[-1][2]
    
    def extractSequences(self,solCCAs):
        seq = {}
        for s in solCCAs:
            if s not in seq:
                seq[s] = []
            ccas = sorted([solCCAs[s][cca] for cca in solCCAs[s]],key=lambda solcca : solcca.getDebut())
            for cca in ccas:
                seq[s]+=cca.getSequence()
        return seq

    def saveSample(self,constellation,add=None):
        if config.getOptValue("sample") is not None:
            self.end_datetime = datetime.now().strftime("%m/%d/%Y, %H:%M")
            folder = '../results/samples/'+config.getOptValue("sample")
            if not os.path.exists(folder):
                os.makedirs(folder)
            found = False
            max_run = 0
            self.generatedModes = {}
            for r in constellation.getRequests():
                self.generatedModes[r] = constellation.getRequest(r).getModesGeneres()
            self.meanScoreAfter = constellation.meanObservationScore(self.bestSelectedModes)
            for root,dirs,files in os.walk(folder+"/"):
                max_run = len(files)
            num_run = str(max_run)
            with open(folder+"/"+config.getOptValue("solver")+"_"+self.instance+'_'+str(self.Ncores)+'_'+str(os.getpid())+'_'+str(num_run)+'.pickle','wb') as handle:
                self.sortObsByRequest(constellation)
                sample = self.toDict()
                if add is not None:
                    for key in add:
                        if key in sample:
                            warning("Clé "+str(key) +" déjà présente. Impossible de sauvegarder dans l'échantillon.")
                        else:
                            sample[key] = add[key]
                pickle.dump(sample, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
    def computeMessageRatio(self):
        try:
            tps_message = sum([x[1]  for cpu in self.stats['envoi'] for x in self.stats['envoi'][cpu]])
            tps_message += sum([x[1] for cpu in self.stats['reception'] for x in self.stats['reception'][cpu]])
            tps_calcul = sum([x[1]  for cpu in self.stats['calculs'] for x in self.stats['calculs'][cpu]])
            ratio = tps_message/tps_calcul
            return ratio
        except:
            return 0
        
    def formatStats(self):
        filtre = self.filtre
        message = ""
        # temps de calcul sur chaque cpu
        dico = {}
        #for stat in ['cpu succes','cpu echec','cpu ignorés','cpu total']:
        if filtre is None or 'ratio' in filtre:
            self.stats['ratio'] = self.computeMessageRatio()
            message += "| ratio tps_comm/tps_calcul : " + str(self.stats['ratio']) + '\n'
        for stat in []:# envoi reception calcul
            dico[stat] = {}
            for cpu in self.stats[stat]:
                dico[stat][cpu] = round(sum([x[1] for x in self.stats[stat][cpu]]),10)
            try:
                if len(mincpu)==0:
                    raise ValueError()
                mincpu = min([dico[stat][cpu] for cpu in dico[stat]])
                meancpu = np.mean([dico[stat][cpu] for cpu in dico[stat]])
                maxcpu = max([dico[stat][cpu] for cpu in dico[stat]])
                stdcpu = np.std([dico[stat][cpu] for cpu in dico[stat]])
            except:
                mincpu,meancpu,maxcpu,stdcpu = np.Inf,None,-np.Inf,0
            message += "| " +str(stat) + " : " + str({'min':mincpu,'mean':meancpu,'max':maxcpu,'std':stdcpu}) + '\n'

        if filtre is None or 'best' in filtre:
            message += "| meilleure solution : " + str(self.bestObjective)  + '\n'

        if filtre is None or 'repartition' in filtre:
            attente,comm,sans_sol,solver = self.getTimeUsage()
            message += "| attente : " + str(round(attente,4)) +" %, communication : " +str(round(comm,4))
            message += " %, calcul hors solver : " + str(round(sans_sol,4)) + " %, calcul solver " + str(round(solver,4)) + " %\n"

        if filtre is None or 'obs' in filtre:
            message += "| score moyenne des observations : " + str(round(self.meanObservationScore,4)) + '\n'
        if filtre is None or 'requete' in filtre:
            message += "| requêtes non planifiées : " + str(self.stats['non planifiées']) + '\n'
            message += "| requêtes planifiées : " + str(self.stats['planifiées']) + '\n'

        if filtre is None or "size" in filtre:
            message += "| #CCA : " + str(self.stats["#CCA"]) + '\n'
            message += "| #A : " + str(self.stats["#A"])
        return message
    
    def getCommunicationDuration(self):
        res = 0
        for cpu in self.stats['reception']:
            #for cca in self.stats['reception'][cpu]:
                for x in self.stats['reception'][cpu]:
                    res += x[1]
                                                           
        for cpu in self.stats['envoi']:
            #for cca in self.stats['envoi'][cpu]:
                for x in self.stats['envoi'][cpu]:
                    res += x[1]
        return res 
    
    def getComputationDuration(self):
        res = 0
        for cpu in self.stats['calculs']:
            #for cca in self.stats['calculs'][cpu]:
                for x in self.stats['calculs'][cpu]:
                    res += x[1]
        return res
                    
    def getTotalDuration(self):
        n = MPI.COMM_WORLD.Get_size()-1
        return n*(time()-self.stats['start date'])
    
    def getTimeUsage(self):
        comm = self.getCommunicationDuration()
        calcul_solver = self.stats['solver time']
        calcul_sans_solver = self.getComputationDuration() - calcul_solver
        attente = self.getTotalDuration() - comm - calcul_solver - calcul_sans_solver
        total = self.getTotalDuration()
        try:
            return attente/total,comm/total,calcul_sans_solver/total,calcul_solver/total
        except Exception:
            return np.Inf,np.Inf,np.Inf,np.Inf
    
    def __str__(self):
        return str(self.formatStats())

    def resetStats(self,constellation):
        self.stats = {'instance':constellation.getFilename(),'messages total':0,'succès':0,'échec':0,'messages ignorés':0}
        self.stats['cpu'] = {}
        self.stats['cpu succes'] = {}
        self.stats['cpu echec'] = {}
        self.stats['cpu ignorés'] = {}
        self.stats['cpu total'] = {}
        self.stats['non planifiées'] = {}
        self.stats['echecs'] = {}
        self.stats['planifiées'] = {}
        self.stats["#A"] = None
        self.stats["#M"] = None
        self.stats["#CCA"] = None
        self.stats['#CCM'] = None
        self.stats['envoi'] = {}
        self.stats['calculs'] = {}
        self.stats['reception'] = {}
        self.stats['ratio'] = 0
        self.stats['solver time'] = 0
        self.stats['start date'] = time()
        for r in constellation.getRequests():
            type_requete = constellation.getRequest(r).getType()
            if type_requete not in self.stats['non planifiées']:
                self.stats['non planifiées'][type_requete] = 0
                self.stats['planifiées'][type_requete] = 0
                self.stats['echecs'][type_requete] = 0
            self.stats['non planifiées'][type_requete] += 1
                
        for clef in ['tailles messages']:
            self.stats[clef] = []
    
    def registerModeFailure(self,constellation,r):
        type_requete = constellation.getRequest(r).getType()
        #self.stats['non planifiées'][type_requete] -= 1
        self.stats['echecs'][type_requete] += 1
        
    def registerMessageReceiving(self,data):
        i = data['source']
        if i not in self.stats["envoi"]:
            self.stats['envoi'][i] = []
            self.stats['reception'][i] = []
            self.stats['calculs'][i] = []
        (t_start,duree) = data['envoi']
        self.stats['envoi'][i].append((t_start-self.startDate,duree))
        (t_start,duree) = data['reception']
        self.stats['reception'][i].append((t_start-self.startDate,duree))
        for (t_start,duree) in data['cpu_time']:
            self.stats['calculs'][i].append((t_start-self.startDate,duree))
        data['cpu_time'] = sum([x[1] for x in data['cpu_time']])
        self.stats["tailles messages"].append(sys.getsizeof(data))
        if 'solver time' in data:
            self.stats['solver time'] += data['solver time']
    
    # total = ignorés + succès + échec
    def registerMessageIgnorance(self,i,cpu_time):
        self.stats['messages ignorés'] += 1
        self.stats['messages total'] += 1
        if i not in self.stats['cpu ignorés']:
            self.stats['cpu ignorés'][i] = 0
        self.stats['cpu ignorés'][i] += cpu_time
        if i not in self.stats['cpu total']:
            self.stats['cpu total'][i] = 0
        self.stats['cpu total'][i] += cpu_time
    
    def setValidatedRequests(self,constellation,modes_retenus):
        for type_requete in self.stats['planifiées']:
            nreq = len([r for r in constellation.getToutesRequetes() if constellation.getRequest(r).getType()==type_requete])
            self.stats['non planifiées'][type_requete] = nreq
            self.stats['planifiées'][type_requete] = 0
        for (r,m) in modes_retenus:
            type_requete = constellation.getRequest(r).getType()
            self.stats['non planifiées'][type_requete] -= 1
            self.stats['planifiées'][type_requete] += 1  
    
    def registerValidatedRequest(self,constellation,r):
        type_requete = constellation.getRequest(r).getType()
        self.stats['non planifiées'][type_requete] -= 1
        self.stats['planifiées'][type_requete] += 1
        
    def registerSuccess(self,i,cpu_time):
        self.stats['succès'] += 1
        self.stats['messages total'] += 1
        if i not in self.stats['cpu succes']:
            self.stats['cpu succes'][i] = 0
        self.stats['cpu succes'][i] += cpu_time
        if i not in self.stats['cpu total']:
            self.stats['cpu total'][i] = 0
        self.stats['cpu total'][i] += cpu_time
    
    def registerFailure(self,i,cpu_time):
        self.stats['échec'] += 1
        self.stats['messages total'] += 1
        if i not in self.stats['cpu echec']:
            self.stats['cpu echec'][i] = 0
        self.stats['cpu echec'][i] += cpu_time
        if i not in self.stats['cpu total']:
            self.stats['cpu total'][i] = 0
        self.stats['cpu total'][i] += cpu_time
    
    def registerMessageSending(self,cpu):
        if cpu not in self.stats['cpu']:
            self.stats['cpu'][cpu] = 1
        else:
            self.stats['cpu'][cpu] += 1
            
    def getBest(self):
        return self.bestObjective,self.bestSolution,self.bestSelectedModes
        
    def plotCPU(self):
        f,ax = plt.subplots(figsize=(15,6))
        plt.xlabel('time (s)')
        plt.ylabel('CPU index')
        tps_max = 0
        tps_min = np.Inf
        for cpu in self.stats['calculs']:
            for (t_start,duree) in self.stats['calculs'][cpu]:
                tps_max = max(tps_max,t_start+duree)
                tps_min = min(tps_min,t_start)
                calcul_handle = ax.add_patch(patches.Rectangle((t_start, cpu),duree,1,edgecolor='black',facecolor = 'green',fill=True,label='calculus') )
        for cpu in self.stats['reception']:
            for (t_start,duree) in self.stats['reception'][cpu]:
                tps_max = max(tps_max,t_start+duree)
                tps_min = min(tps_min,t_start)
                rcv_handle = ax.add_patch(patches.Rectangle((t_start, cpu),duree,1,facecolor = 'red',fill=True,label='answer to master') )
        for cpu in self.stats['envoi']:
            for (t_start,duree) in self.stats['envoi'][cpu]:
                tps_max = max(tps_max,t_start+duree)
                tps_min = min(tps_min,t_start)
                send_handle = ax.add_patch(patches.Rectangle((t_start, cpu),duree,1,facecolor = 'orange',fill=True,label='send to slave') )            
        plt.legend(handles=[send_handle,calcul_handle,rcv_handle])
        plt.title('cpu occupation')
        plt.title('cpu occupation : '+self.formatInstance()+"_"+str(self.Ncores)+" ; obj = "+str(self.bestObjective))
        plt.ylim(0,max(list(self.stats["envoi"].keys()))+1)
        plt.xlim(tps_min,tps_max)
        plt.savefig("../results/cpu/"+config.getOptValue("solver")+"_"+self.formatInstance()+"_"+str(self.Ncores)+".png")
        return f
    
    # init = False => ne pas inclure le premier point (obj = 0)
    def plotObjective(self,init=True):
        if(not init):
            hist = self.historic[1:]
        else:
            hist = self.historic

        X = [x[0] for x in hist]
        Y = [x[2][0] for x in hist]
        col = []
        f,ax = plt.subplots(figsize=(15,6))
        ylim = (min([x[2][0] for x in hist]),max([x[2][0] for x in hist]))
        for event,couleur in zip([PREPROCESSING,BACKUP,END_ITERATION,END_RUN],['y','k','m','r']):
            for x in [xx for xx in hist if xx[1]==event]:
                ax.plot([x[0],x[0]],[ylim[0],ylim[1]],c=couleur,label=event)
        ax.plot(X,Y)
        plt.xlabel('time (s)')
        plt.ylabel('reward')
        algo = str(config.getOptValue("solver"))
        instance = self.formatInstance()
        plt.title(instance + " - " +algo + " - " +str(MPI.COMM_WORLD.Get_size()) + ' processes - Objective : ' + str(self.bestObjective[0]) +" " + str(len(self.bestSelectedModes))+" requests")
        handles, labels = plt.gca().get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        plt.legend(by_label.values(), by_label.keys())
        plt.savefig("../results/objective/obj_"+instance + "_" +algo + "_" +str(MPI.COMM_WORLD.Get_size())+"_proc.png")
        return f
    
    def plotCCAsLoad(self,grapheDep):
        f, ax = plt.subplots()
        ccas = []
        charge = []
        i = 0
        for s in self.bestSolution:
            for cca in self.bestSolution[s]:
                ccas.append(i)
                i += 1
                charge.append(len(self.bestSolution[s][cca].getSequence()))        
        ax.bar(ccas, charge)
        ax.set_ylabel('composante')
        ax.set_title('nombre d\'observations')
        #ax.legend(title='charge des composantes')
        instance = self.formatInstance()
        algo = str(config.getOptValue("solver"))
        plt.savefig("../results/charge_cca/charge_"+instance + "_" +algo + "_" +str(MPI.COMM_WORLD.Get_size())+"_proc.png")
        return f    

    def bestUpdate(self,objective,solution,modes_retenus):
        if(objective>=self.bestObjective):
            self.bestSelectedModes = modes_retenus.copy()
            if solution is not None:
                self.bestSolution = deepcopy(solution)
            self.bestObjective = objective
            #assert(solution is not None)

    def componentsBounds(self,composante):
        return composante.getDebut(),composante.getFin()
    
    def rejectBestSolution(self):
        self.bestObjective = (0,0)
        self.bestSelectedModes = []
    
    # Event : 0 = rien, 1 = fin itération
    def updateHistory(self,time,event,objective,solution,modes_retenus,composantes,constellation):
        if self.ccas is None:
            self.ccas = {scca : self.componentsBounds(composantes.getComposanteConnexe(scca)) for scca in composantes.getComposantes()}
        
        time = min(time,config.getOptValue("time"))
        self.setRequetesValidees(constellation,modes_retenus)
        ccas = composantes.getNombreComposantes()
        self.meanObservationScore = constellation.meanObservationScore(modes_retenus)
        self.stats['#CCA'] = ccas
        self.stats["#A"] = composantes.getNombreElements()
        seq = {s : {cca : solution[s][cca].getSequence() for cca in solution[s]} for s in solution}
        if not config.getOptValue("full_sample"):
            self.historic.append((time,event,objective,len(modes_retenus),ccas,MPI.COMM_WORLD.Get_rank(),None,None))
        else:
            self.historic.append((time,event,objective,len(modes_retenus),ccas,MPI.COMM_WORLD.Get_rank(),modes_retenus,seq))
        self.bestUpdate(objective,solution,modes_retenus)
        self.stats['ratio'] = self.computeMessageRatio()
    
    def deleteLastPoint(self):
        self.historic.pop(-1)

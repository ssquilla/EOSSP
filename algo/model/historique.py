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

from .constellation import TYPE_LONG_MONO,TYPE_MONO,TYPE_PERIODIC,TYPE_STEREO,TYPE_SYSTEMATIC,TYPE_VIDAGE

from ..Utils.config import *
global config
config = Config()

PREPROCESSING = "PREPROCESSING"
END_ITERATION = "END_ITERATION"
NO_EVENT = "NO_EVENT"
END_RUN = "END_RUN"
BACKUP = "BACKUP"
END_OPERATEUR = "END_OPERATEUR"

class Historique:
    def __init__(self,constellation,start_date,dico=None):
        self.start_date = start_date
        self.filtre = None
        # flag d'évenements : indiquer l'état de la recherche
        self.extraireQuotaTypesRequetes(constellation)
        if dico is None:
            self.historique = [(0,NO_EVENT,(0,0),0,0,MPI.COMM_WORLD.Get_rank(),[],{})] # temps,event,objectif,nb modes retenus,mpi rank,ccas,sol
            self.start_datetime = datetime.now().strftime("%m/%d/%Y, %H:%M")
            self.best_objectif = (0,0)
            self.best_solution = {}
            self.best_modes_retenus = []
            self.resetStats(constellation)
            self.ccas = None
            self.arrivee_requetes = {}
            instance = constellation.getFilename()
            Ncoeur = str(MPI.COMM_WORLD.Get_size())
            self.Ncoeur = Ncoeur
            self.instance = instance
            self.algo_name = config.donnees.algo_name
            self.config = str(config)
            self.scoreObsMoyenne = 0
            self.gap = None
            self.commentaires = []
            self.path = os.getcwd()
            self.information_additionnelle = {}
        else:
            self.arrivee_requetes = self.extraireParametre(dico,'arrivee_requetes')
            self.ccas = self.extraireParametre(dico,'cca')
            self.start_datetime = self.extraireParametre(dico,'start')
            self.end_datetime = self.extraireParametre(dico,'end')
            self.solution.historique = self.extraireParametre(dico,'hist') # temps,event,objectif,nb modes retenus
            self.best_objectif = self.extraireParametre(dico,'obj')
            self.best_solution = self.extraireParametre(dico,'sol')
            self.best_modes_retenus = self.extraireParametre(dico,'modes')
            self.stats = self.extraireParametre(dico,'stats')
            self.Ncoeur = self.extraireParametre(dico,'coeurs')
            self.instance = self.extraireParametre(dico,'instance')
            self.algo_name = self.extraireParametre(dico,'algo')
            self.config = self.extraireParametre(dico,'config')
            self.score_moyen_avant = self.extraireParametre(dico,"score_avant")
            self.score_moyen_apres = self.extraireParametre(dico,"score_apres")
            self.scoreObsMoyenne = self.score_moyen_apres
            self.commentaires = self.extraireParametre(dico,'commentaires')
            self.path = self.extraireParametre(dico,'path')
            self.infos_hist = self.extraireParametre(dico,'infos hist')
            self.information_additionnelle = self.extraireParametre(dico, 'info_add')
            if 'gap' in dico:
                self.gap = dico['gap']
            else:
                self.gap = None
    
    def ajouterInformationAdditionnelle(self,cle,valeur):
        self.information_additionnelle[cle] = valeur
    
    def extraireQuotaTypesRequetes(self,constellation):
        types = {}
        types[TYPE_LONG_MONO] = 0
        types[TYPE_MONO] = 0
        types[TYPE_PERIODIC] = 0
        types[TYPE_STEREO] = 0
        types[TYPE_SYSTEMATIC] = 0
        for r in constellation.getRequetes():
            type_r = constellation.getRequete(r).getType()
            if type_r != TYPE_VIDAGE:
                if r not in types:
                    types[type_r] = 0
                types[type_r] += 1
        self.types = types
    
    def formaterInstance(self):
        res = ""
        for i,t in enumerate(self.types):
            res += str(self.types[t])
            if i<len(list(self.types.keys()))-1:
                res += ','
        return res
    
    def extraireParametre(self,dico,cle):
        return dico.get(cle,None)
        
    def fusionner(self,historiques_cpu):
        for cpu in historiques_cpu:
            hist = historiques_cpu[cpu]
            if hist.best_objectif>self.best_objectif:
                self.best_objectif = hist.best_objectif
                self.best_solution = hist.best_solution
                self.modes_retenus = hist.best_modes_retenus
        historiques_cpu[MPI.COMM_WORLD.Get_rank()] = self
        new_hist = deepcopy(self.historique)
        for cpu in historiques_cpu:
            for x in historiques_cpu[cpu].historique:
                if x not in new_hist:
                    new_hist.append(x)
        self.historique = sorted(new_hist,key=operator.itemgetter(0))

                        
    def ajouterCommentaire(self,commentaire):
        self.commentaires.append(commentaire)
        
    def trierObsParRequete(self,constellation):
        self.obsParRequete = {}
        seq = self.extraireSequences(self.best_solution)
        for s in seq:
            for a in seq[s]:
                req = constellation.getRequeteActivite(a)
                if req not in self.obsParRequete:
                    self.obsParRequete[req] = []
                self.obsParRequete[req].append(a)
        
    def toDict(self):
        dico = {}
        dico['hist'] = self.historique
        dico['infos hist'] = self.getInfos()
        dico['obj'] = self.best_objectif
        dico['sol'] = self.extraireSequences(self.best_solution)
        dico['obs_par_requetes'] = self.obsParRequete
        dico['modes'] = self.best_modes_retenus
        dico['stats'] = self.stats
        dico['instance'] = self.instance
        dico['coeurs'] = self.Ncoeur
        dico['algo'] = self.algo_name
        dico['config'] = str(self.config)
        dico["score_avant"] = self.score_moyen_avant
        dico["score_apres"] = self.score_moyen_apres
        dico['gap'] = self.gap
        dico['start'] = self.start_datetime
        dico['end'] = self.end_datetime
        dico['commentaires'] = self.commentaires
        dico['path'] = self.path
        dico['cca'] = self.ccas
        dico['arrivee_requetes'] = self.arrivee_requetes
        dico['info_add'] = self.information_additionnelle
        dico['modes_generes'] = self.modes_generes
        return dico
    
    def notifierArriveeRequetes(self,date,liste_requetes):
        self.arrivee_requetes[date] = liste_requetes
    
    def getInfos(self):
        mess = " ----------- Historique de l'execution ----------\n"
        lecture = [("chemin",self.path)]
        for i in range(len(self.commentaires)):
            lecture.append(("commentaire "+str(i),self.commentaires[i]))
        lecture += [("algorithme",self.algo_name),('instance',self.instance),("nombre de coeurs",self.Ncoeur)]
        lecture += [("date de démarrage",self.start_datetime),("date de fin",self.end_datetime)]
        #lecture += [("objectif",self.best_objectif),("modes retenus",len(self.best_modes_retenus))]
        lecture += [("historique",str(len(self.historique)) + " valeurs enregistrées")]
        for message,contenu in lecture:
            mess += " | "+ message + " : " + str(contenu) +'\n'
        mess += ' '
        mess += self.formatStats() +'\n'
        mess += str(self.config)
        mess += " ------------------------------------------------\n"
        return mess
            
    def dernierObjectif(self):
        return self.historique[-1][2]
    
    def extraireSequences(self,solCCAs):
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
            self.modes_generes = {}
            for r in constellation.getRequetes():
                self.modes_generes[r] = constellation.getRequete(r).getModesGeneres()
            self.score_moyen_apres = constellation.scoreObsMoyenne(self.best_modes_retenus)
            for root,dirs,files in os.walk(folder+"/"):
                max_run = len(files)
            num_run = str(max_run)
            with open(folder+"/"+config.getOptValue("solver")+"_"+self.instance+'_'+str(self.Ncoeur)+'_'+str(os.getpid())+'_'+str(num_run)+'.pickle','wb') as handle:
                self.trierObsParRequete(constellation)
                sample = self.toDict()
                if add is not None:
                    for key in add:
                        if key in sample:
                            warning("Clé "+str(key) +" déjà présente. Impossible de sauvegarder dans l'échantillon.")
                        else:
                            sample[key] = add[key]
                pickle.dump(sample, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
    def calculerRatioMessage(self):
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
            self.stats['ratio'] = self.calculerRatioMessage()
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
        #gap_rec,gap_temps = 0,0
        #try:
        #    gap_rec = round(100*(self.majorant[0] - self.best_objectif[0])/self.majorant[0],2)
        #    gap_temps = round(100*(self.majorant[1] - self.best_objectif[1])/self.majorant[1],2)
        #except:
        #    pass
        if filtre is None or 'best' in filtre:
            message += "| meilleure solution : " + str(self.best_objectif)  + '\n'
        #try:
        #    message += " | majorant : " + str(self.majorant) + '\n'
        #except:
        #    pass
        #message += " | durée solver : " + str(self.stats['solver time']) + "\n"
        if filtre is None or 'repartition' in filtre:
            attente,comm,sans_sol,solver = self.getRepartitionTemps()
            message += "| attente : " + str(round(attente,4)) +" %, communication : " +str(round(comm,4))
            message += " %, calcul hors solver : " + str(round(sans_sol,4)) + " %, calcul solver " + str(round(solver,4)) + " %\n"
        
        #message += " | gap récompense : " + str(gap_rec) + " %" + '\n'
        #message += " | gap temps : " + str(gap_temps) + " %" + '\n'
        if filtre is None or 'obs' in filtre:
            message += "| score moyenne des observations : " + str(round(self.scoreObsMoyenne,4)) + '\n'
        if filtre is None or 'requete' in filtre:
            message += "| requêtes non planifiées : " + str(self.stats['non planifiées']) + '\n'
            message += "| requêtes planifiées : " + str(self.stats['planifiées']) + '\n'
        #message += " | modes echecs : " + str(self.stats["echecs"]) + '\n'
        if filtre is None or "size" in filtre:
            message += "| #CCA : " + str(self.stats["#CCA"]) + '\n'
            message += "| #A : " + str(self.stats["#A"])
        return message
    
    def getDureeCommunication(self):
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
    
    def getTempsCalcul(self):
        res = 0
        for cpu in self.stats['calculs']:
            #for cca in self.stats['calculs'][cpu]:
                for x in self.stats['calculs'][cpu]:
                    res += x[1]
        return res
                    
    def getDureeTotale(self):
        n = MPI.COMM_WORLD.Get_size()-1
        return n*(time()-self.stats['start date'])
    
    def getRepartitionTemps(self):
        comm = self.getDureeCommunication()
        calcul_solver = self.stats['solver time']
        calcul_sans_solver = self.getTempsCalcul() - calcul_solver
        attente = self.getDureeTotale() - comm - calcul_solver - calcul_sans_solver
        total = self.getDureeTotale()
        try:
            return attente/total,comm/total,calcul_sans_solver/total,calcul_solver/total
        except Exception:
            return np.Inf,np.Inf,np.Inf,np.Inf
    
    def __str__(self):
        return str(self.formatStats())

    def majorantRecompense(self,constellation):
        return (sum([constellation.getRequete(r).getMode(0).getRecompense() for r in constellation.getToutesRequetes()]),len(self.best_modes_retenus))

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
        for r in constellation.getRequetes():
            type_requete = constellation.getRequete(r).getType()
            if type_requete not in self.stats['non planifiées']:
                self.stats['non planifiées'][type_requete] = 0
                self.stats['planifiées'][type_requete] = 0
                self.stats['echecs'][type_requete] = 0
            self.stats['non planifiées'][type_requete] += 1
                
        for clef in ['tailles messages']:
            self.stats[clef] = []
    
    def enregistrerEchecMode(self,constellation,r):
        type_requete = constellation.getRequete(r).getType()
        #self.stats['non planifiées'][type_requete] -= 1
        self.stats['echecs'][type_requete] += 1
        
    def enregistrerReceptionMessage(self,data):
        i = data['source']
        if i not in self.stats["envoi"]:
            self.stats['envoi'][i] = []
            self.stats['reception'][i] = []
            self.stats['calculs'][i] = []
        (t_start,duree) = data['envoi']
        self.stats['envoi'][i].append((t_start-self.start_date,duree))
        (t_start,duree) = data['reception']
        self.stats['reception'][i].append((t_start-self.start_date,duree))
        for (t_start,duree) in data['cpu_time']:
            self.stats['calculs'][i].append((t_start-self.start_date,duree))
        data['cpu_time'] = sum([x[1] for x in data['cpu_time']])
        self.stats["tailles messages"].append(sys.getsizeof(data))
        if 'solver time' in data:
            self.stats['solver time'] += data['solver time']
    
    # total = ignorés + succès + échec
    def enregistrerIgnoranceMessage(self,i,cpu_time):
        self.stats['messages ignorés'] += 1
        self.stats['messages total'] += 1
        if i not in self.stats['cpu ignorés']:
            self.stats['cpu ignorés'][i] = 0
        self.stats['cpu ignorés'][i] += cpu_time
        if i not in self.stats['cpu total']:
            self.stats['cpu total'][i] = 0
        self.stats['cpu total'][i] += cpu_time
    
    def setRequetesValidees(self,constellation,modes_retenus):
        for type_requete in self.stats['planifiées']:
            nreq = len([r for r in constellation.getToutesRequetes() if constellation.getRequete(r).getType()==type_requete])
            self.stats['non planifiées'][type_requete] = nreq
            self.stats['planifiées'][type_requete] = 0
        for (r,m) in modes_retenus:
            type_requete = constellation.getRequete(r).getType()
            self.stats['non planifiées'][type_requete] -= 1
            self.stats['planifiées'][type_requete] += 1  
    
    def enregistrerRequeteValidee(self,constellation,r):
        type_requete = constellation.getRequete(r).getType()
        self.stats['non planifiées'][type_requete] -= 1
        self.stats['planifiées'][type_requete] += 1
        
    def enregistrerSucces(self,i,cpu_time):
        self.stats['succès'] += 1
        self.stats['messages total'] += 1
        if i not in self.stats['cpu succes']:
            self.stats['cpu succes'][i] = 0
        self.stats['cpu succes'][i] += cpu_time
        if i not in self.stats['cpu total']:
            self.stats['cpu total'][i] = 0
        self.stats['cpu total'][i] += cpu_time
    
    def enregistrerEchec(self,i,cpu_time):
        self.stats['échec'] += 1
        self.stats['messages total'] += 1
        if i not in self.stats['cpu echec']:
            self.stats['cpu echec'][i] = 0
        self.stats['cpu echec'][i] += cpu_time
        if i not in self.stats['cpu total']:
            self.stats['cpu total'][i] = 0
        self.stats['cpu total'][i] += cpu_time
    
    def enregistrerEnvoi(self,cpu):
        if cpu not in self.stats['cpu']:
            self.stats['cpu'][cpu] = 1
        else:
            self.stats['cpu'][cpu] += 1
            
    def getBest(self):
        return self.best_objectif,self.best_solution,self.best_modes_retenus
        
    def tracerCPU(self):
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
        plt.title('cpu occupation : '+self.formaterInstance()+"_"+str(self.Ncoeur)+" ; obj = "+str(self.best_objectif))
        plt.ylim(0,max(list(self.stats["envoi"].keys()))+1)
        plt.xlim(tps_min,tps_max)
        plt.savefig("../results/cpu/"+config.getOptValue("solver")+"_"+self.formaterInstance()+"_"+str(self.Ncoeur)+".png")
        return f
    
    # init = False => ne pas inclure le premier point (obj = 0)
    def tracerObjectif(self,init=True):
        if(not init):
            hist = self.historique[1:]
        else:
            hist = self.historique

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
        instance = self.formaterInstance()
        plt.title(instance + " - " +algo + " - " +str(MPI.COMM_WORLD.Get_size()) + ' processes - Objective : ' + str(self.best_objectif[0]) +" " + str(len(self.best_modes_retenus))+" requests")
        handles, labels = plt.gca().get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        plt.legend(by_label.values(), by_label.keys())
        plt.savefig("../results/objective/obj_"+instance + "_" +algo + "_" +str(MPI.COMM_WORLD.Get_size())+"_proc.png")
        return f
    
    def tracerChargeCCAs(self,grapheDep):
        f, ax = plt.subplots()
        ccas = []
        charge = []
        i = 0
        for s in self.best_solution:
            for cca in self.best_solution[s]:
                ccas.append(i)
                i += 1
                charge.append(len(self.best_solution[s][cca].getSequence()))        
        ax.bar(ccas, charge)
        ax.set_ylabel('composante')
        ax.set_title('nombre d\'observations')
        #ax.legend(title='charge des composantes')
        instance = self.formaterInstance()
        algo = str(config.getOptValue("solver"))
        plt.savefig("../results/charge_cca/charge_"+instance + "_" +algo + "_" +str(MPI.COMM_WORLD.Get_size())+"_proc.png")
        return f    

    # nombre de fois ou une obs est impliquee dans un mode
    def tracerImplicationObs(self,grapheDep,constellation,tri_cca=False):
        f, ax = plt.subplots()
        charge = {}
        for r in constellation.getRequetes():
            for m in constellation.getRequete(r).getModes():
                for (s,a) in constellation.getRequete(r).getMode(m).getCouples():
                    if tri_cca:
                        cca = grapheDep.getActiviteCCA(a)
                        if cca not in charge:
                            charge[cca] = 0
                        charge[cca] += 1
                    else:
                        if a not in charge:
                            charge[a] = 0
                        charge[a] += 1
                        
        
        ax.bar(range(len(list(charge.keys()))), list(charge.values()))
        ax.set_title('nombre d\'observations présentes dans les modes par composante')
        ax.set_ylabel('nombre d\'observations présentes sur la composante')
        ax.set_xlabel('composante')
        #ax.legend(title='charge des composantes')
        instance = self.formaterInstance()
        algo = str(config.getOptValue("solver"))
        plt.savefig("../results/charge_cca/obs_"+instance + "_" +algo + "_" +str(MPI.COMM_WORLD.Get_size())+"_proc.png")
        return f  

    def bestMAJ(self,objectif,solution,modes_retenus):
        if(objectif>=self.best_objectif):
            self.best_modes_retenus = modes_retenus.copy()
            if solution is not None:
                self.best_solution = deepcopy(solution)
            self.best_objectif = objectif
            #assert(solution is not None)

    def bornesComposante(self,composante):
        return composante.getDebut(),composante.getFin()
    
    def invaliderMeilleureSolution(self):
        self.best_objectif = (0,0)
        self.best_modes_retenus = []
    
    # Event : 0 = rien, 1 = fin itération
    def MAJHistorique(self,time,event,objectif,solution,modes_retenus,composantes,constellation):
        if self.ccas is None:
            self.ccas = {scca : self.bornesComposante(composantes.getComposanteConnexe(scca)) for scca in composantes.getComposantes()}
        
        time = min(time,config.getOptValue("time"))
        self.setRequetesValidees(constellation,modes_retenus)
        ccas = composantes.getNombreComposantes()
        self.scoreObsMoyenne = constellation.scoreObsMoyenne(modes_retenus)
        self.stats['#CCA'] = ccas
        self.stats["#A"] = composantes.getNombreElements()
        seq = {s : {cca : solution[s][cca].getSequence() for cca in solution[s]} for s in solution}
        if not config.getOptValue("full_sample"):
            self.historique.append((time,event,objectif,len(modes_retenus),ccas,MPI.COMM_WORLD.Get_rank(),None,None))
        else:
            self.historique.append((time,event,objectif,len(modes_retenus),ccas,MPI.COMM_WORLD.Get_rank(),modes_retenus,seq))
        self.bestMAJ(objectif,solution,modes_retenus)
        #self.majorant = self.majorantRecompense(constellation)
        self.stats['ratio'] = self.calculerRatioMessage()
    
    def supprimerDernierPoint(self):
        self.historique.pop(-1)

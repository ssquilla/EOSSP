from ..Utils.config import *
from ..Utils.Utils import *
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
    from ..model.composantes import *
    from ..model.solution_composantes import *
    
    from time import time
    from time import sleep
    import random as rd
    
    class BatchParallelCCASearch(Solver):
        def __init__(self,constellation,start,modeleDeTransition,dt_construction_transition,colors={'iteration':'c',"no_event":'w','end':'y'},tlim=np.inf,CCAs=None,solution=None):
            super().__init__(constellation,start,modeleDeTransition,dt_construction_transition,tlim=tlim,CCAs=CCAs,solution=solution)
            self.colors = colors
            self.explications = {}
            if MPI.COMM_WORLD.Get_size()==1: # structure de données des slaves si processus seul
                self.echecs_sans_explications = []
                self.indetermines = []
                self.sequences = {}
            
            comm = MPI.COMM_WORLD
            if comm.Get_size()==1:
                self.localSlave = Slave(constellation,start,modeleDeTransition)
            printColor("Durée d'initialisation :",time()-self.start_date,c='b')
            # mapper les cca
            #self.creerMappingCCASlaves()
        
        def modeDestructionActive(self):
            return config.getOptValue("destroy_BPCCAS") is not None and config.getOptValue("destroy_BPCCAS")>0
            
        def modeBruitageActive(self):
            return config.getOptValue("noise")>0
    
        def redemarrer(self,constellation,sigma):
            self.solution.redemarrer(constellation)
            for (s,cca) in self.grapheDependances.getComposantes():
                self.setSolCCA(s,cca,SolCCA(cca,s))
        
        def initIteration(self,constellation,sigma,super_iter):
            self.initRequetes(constellation,sigma,conserve_old=False)
            if super_iter>0 or self.modeBruitageActive():
                self.redemarrer(constellation,sigma)
            if super_iter!=0 and self.modeDestructionActive():
                self.detruire(constellation)
            self.notifierSansEvenement(constellation)
                
                
        def detruire(self,constellation):
            printOpen("Destruction",c='m')
            if config.getOptValue("verif"):
                self.verifierSolution(constellation)        
            
            couvertes = self.requetesCouvertes()
            detruire = rd.choices(couvertes,k=int(len(couvertes)*config.getOptValue("destroy_BPCCAS")))
            if -1 in detruire:
                detruire.remove(-1)
            modes_a_retirer = [(r,m) for (r,m) in self.getModesRetenus() if r in detruire and r!=-1]
            self.setModesRetenus([(r,m) for (r,m) in self.getModesRetenus() if (r, m) not in modes_a_retirer],constellation)
            activites_a_retirer = [(s,a) for (r,m) in modes_a_retirer for (s,a) in constellation.getRequete(r).getMode(m).getCouples()]
            activites_par_cca = {}
            for (s,a) in activites_a_retirer:
                id_cca = self.grapheDependances.getActiviteCCA(a)
                if id_cca not in activites_par_cca:
                    activites_par_cca[id_cca] = []
                activites_par_cca[id_cca].append(a)
            for (s,cca) in activites_par_cca:
                self.getSolCCA(s,cca).retirerListeActivites(constellation,activites_par_cca[cca])
            if config.getOptValue("verif"):
                self.verifierSolution(constellation)
            self.MAJDestruction(constellation)
            printColor(str(len(detruire))+" requêtes detruites",c='y')
            printClose()
            
        def MAJDestruction(self,constellation):
            couvertes = self.requetesCouvertes()
            for r in constellation.getRequetes():
                if r not in couvertes:
                    constellation.getRequete(r).init = False
                    m = constellation.getRequete(r).getModeCourant(constellation).getId()
                    self.modes_courants[r] = m
                    assert(m is not None)
            # score candidats modes avant la planif
            modes = []
            for r in self.requetes_candidates:
                modes.append((r,0))
            self.resetScoreDestruction(constellation,couvertes)
            self.solution.historique.score_moyen_avant = constellation.scoreObsMoyenne(modes)
            self.solution.historique.scoreObsMoyenne = self.solution.historique.score_moyen_avant
            
            
        def creerMappingCCASlaves(self):
            if MPI.COMM_WORLD.Get_size()>1:
                self.cca_slaves = {cpu : [] for cpu in range(1,MPI.COMM_WORLD.Get_size())}
                tailles = []
                for id_cca in self.grapheDependances.getComposantes():
                    taille =  self.grapheDependances.getTailleComposante(id_cca)
                    reverse_insort(tailles,(taille,id_cca))
                cpu = 0
                for (taille,id_cca) in tailles:
                    self.cca_slaves[cpu+1].append(id_cca)
                    cpu = (cpu + 1) % (MPI.COMM_WORLD.Get_size()-1)
                for cpu in self.cca_slaves:
                    MPI.COMM_WORLD.send({"cca":self.cca_slaves[cpu]},dest=cpu)
            else:
                self.cca_slaves = {1 : self.grapheDependances.getComposantes()}
         
        def mappingCCATaille(self,tailles):
            # tailles : liste (cca,nb activites a inserer)
            ccas = sorted(tailles,key=itemgetter(0),reverse=True)
            if MPI.COMM_WORLD.Get_size()==1:
                self.cca_slaves = {0:[]}
                for (taille,cca) in ccas:
                    self.cca_slaves[0].append(cca)
            else:
                self.cca_slaves = {cpu : [] for cpu in range(1,MPI.COMM_WORLD.Get_size())}
                cpu = 0
                for (taille,id_cca) in ccas:
                    self.cca_slaves[cpu+1].append(id_cca)
                    cpu = (cpu + 1) % (MPI.COMM_WORLD.Get_size()-1)
    
        """
            =============================================== 
                            RESOLUTION
            =============================================== 
        """
        def reinsererMode(self,constellation,r):
            pass
    
        def analyserMessage(self,data,constellation,freeCPUs,mailbox):
            # les sequences solutions ne sont pas transmises ici car les cca sont statiques. 
            # Les slaves ont deja un exemplaire
            start = time()
            self.solution.historique.enregistrerReceptionMessage(data)
            i = data['source']
            #del self.job_slaves[i]['time']
            assert(i not in freeCPUs)
            bisect.insort(freeCPUs,i)
            cpu_time = data['cpu_time']
            explications = data['explications']
            indetermines = data['indetermines']
            return explications,indetermines
        
        def retirerAnciennesActivites(self,constellation,r,m):
            for (s,o) in constellation.getRequete(r).getMode(m-1).getCouples():
                if (s,o) not in constellation.getRequete(r).getMode(m).getCouples():
                    (s,cca) = self.grapheDependances.getActiviteCCA(o)
                    if o in self.getSolCCA(s,cca).getSequence():
                        self.getSolCCA(s,cca).retirerActivite(constellation,o)
        
        def analyserResultatsAsynchrone(self,constellation,freeCPUs,mailbox):
            failed = list(self.explications.keys())+self.echecs_sans_explications
            modes_reussis = []
            for i in self.job_slaves:
                for cca in self.job_slaves[i]:
                    for (r,m,score) in self.job_slaves[i][cca]:
                        if (r,m) not in modes_reussis and (r,m) not in failed and (r,m) not in self.indetermines:
                            modes_reussis.append((r,m))
            for (r,m) in self.explications:
                if len(self.explications[(r,m)])>0:
                    m_new = constellation.getRequete(r).getModeSuivant(self.explications[(r,m)],constellation)
                    if m_new is None:
                        self.requetes_candidates.remove(r)
                    else:
                        self.retirerAnciennesActivites(constellation,r,m_new.getId())
                        
                #else:
                #    self.reinsererMode(constellation,r)
            
            for rm in modes_reussis:
                if rm not in self.getModesRetenus():
                    self.ajouterModeRetenu(rm,constellation)
            self.insererSequences(constellation,self.sequences)
            
            printOpen("Résultats",c='b')
            printOpen("modes ratés avec explications :",len(list(self.explications.keys())),c='r')
            printColor((list(self.explications.keys())),c='r')
            printClose()
            printOpen("modes ratés sans explications :",len(self.echecs_sans_explications),c='m')
            printColor(self.echecs_sans_explications,c='m')
            printClose()
            printOpen("modes réussis :",len(modes_reussis),c='g')
            printColor((modes_reussis),c='g')
            printClose()
            printOpen("modes indeterminés:",len(self.indetermines),c='w')
            printColor((self.indetermines),c='w')
            printClose()
            printClose()
            #die("AnalyserResultats")
            return failed
        
        def analyserResultats(self,constellation,freeCPUs,mailbox):
            comm = MPI.COMM_WORLD
            rank = comm.Get_rank()
            start_iteration = time()
            if comm.Get_size()==1:
                messages = [self.localSlave.getResultat()]
            else:
                messages = mailbox.readMessages()
            
            explications = {}
            sequences = {}
            indetermines = []
            for data in messages:
                if data is not None:
                    for cca in data['sequences']:
                        sequences[cca] = data['sequences'][cca]
                    exp_cpu,inde = self.analyserMessage(data,constellation,freeCPUs,mailbox)
                    for (r,m) in inde:
                        if (r,m) not in indetermines:
                            indetermines.append((r,m))
                    for (r,m) in exp_cpu:
                        if (r,m) not in explications:
                            explications[(r,m)] = []
                        explications[(r,m)] += exp_cpu[(r,m)]
                   
            modes_reussis = []
            for i in self.job_slaves:
                for cca in self.job_slaves[i]:
                    for (r,m,score) in self.job_slaves[i][cca]:
                        if (r,m) not in modes_reussis and (r,m) not in explications and (r,m) not in indetermines:
                            modes_reussis.append((r,m))
            
            for (r,m) in explications:
                m_new = constellation.getRequete(r).getModeSuivant(explications[(r,m)],constellation)
                if m_new is None:
                    self.requets_candidates.remove(r)
            
            for rm in modes_reussis:
                if rm not in self.getModesRetenus():
                    self.ajouterModeRetenu(rm,constellation)
            self.insererSequences(constellation,sequences)
            
            printOpen("Résultats",c='b')
            printColor("modes ratés :",list(explications.keys()),c='r')
            printColor("modes réussis :",modes_reussis,c='g')
            printColor("modes indeterminés:",indetermines,c='c')
            printClose()
            #die("AnalyserResultats")
            return list(explications.keys())
        
        def insererSequences(self,constellation,sequences):
            requetes_retenues = [x[0] for x in self.getModesRetenus()]
            for (s,cca) in sequences:
                seq = [a for a in sequences[(s,cca)] if constellation.getRequeteActivite(a) in requetes_retenues]
                self.getSolCCA(s,cca).setSequence(constellation,seq,self.modeleDeTransition)
    
        def resoudre(self,constellation,mailbox):
            comm = MPI.COMM_WORLD
            rank = comm.Get_rank()
            size = comm.Get_size()
            #objectif = self.calculerObjectif(constellation)
            if size==1:
                freeCPUs = [0] # un seul coeur => slave local
            else:
                freeCPUs = list(range(1,size))
            last_record = 0
            failed = []
            super_iter = 0
            it = 0
            while time()-self.start_date<self.tlim:
                sigma = config.getOptValue("noise")*super_iter
                self.initIteration(constellation,sigma,super_iter)
                if (self.modeBruitageActive() or self.modeDestructionActive()):
                    couleur = self.colors['iteration']
                    self.afficherInfo(time(),self.start_date,constellation,color=couleur,add={"sigma":sigma},title='DEBUT ITERATION '+str(super_iter))
                while time()-self.start_date<config.getOptValue("time") and (len(failed)!=0 or it==0):
                    start_it = time()
                    self.envoyerTacheAsynchrone(constellation,freeCPUs,mailbox,failed)
                    failed = self.analyserResultatsAsynchrone(constellation,freeCPUs,mailbox)
                    if time()-self.start_date<config.getOptValue("time") and (len(failed)!=0 or it==0):
                        self.notifierSansEvenement(constellation)
                        couleur = self.colors['no_event']
                        self.afficherInfo(time(),self.start_date,constellation,color=couleur,title="RESULTAT SOUS-ITERATION " +str(it))
                    else:
                        break
                    it += 1
                    temps = time()
                    if config.getOptValue("verif"):
                        self.verifierCCA(constellation)
                    step()
                if config.getOptValue("verif"):
                    self.verifierSolution(constellation)
                self.notifierFinIteration(constellation)
                #couleur = self.colors['iteration']
                #self.afficherInfo(time(),self.start_date,constellation,color=couleur,add={"sigma":sigma},title='FIN ITERATION '+str(super_iter))
                super_iter += 1 
                it = 0
                up = self.majorantRecompense(constellation)[0]
                
                if (not self.modeBruitageActive() and not self.modeDestructionActive()) :#or  self.gapRecompense(constellation) and self.gapTemps(constellation) :
                    break
            
                
            self.terminerProcessus()
            self.notifierFinExecution(constellation)
            couleur = self.colors['end']
            self.afficherInfo(time(),self.start_date,constellation,color=couleur,title='FIN')
    
        def mappingCCAObservationsCourantes(self,constellation,failed):
            mapping_obs = {}
            
            del_r = []
            for r in self.requetes_candidates:
                if r in self.requetesCouvertes():
                    del_r.append(r)
            for r in del_r:
                self.requetes_candidates.remove(r)
                
            for r in self.requetes_candidates:
                m = constellation.getRequete(r).getModeCourant(constellation).getId()
                score = self.recompenseBruitee(constellation,r,m)
                for (s,o) in constellation.getRequete(r).getMode(m).getCouples():
                    id_cca = self.grapheDependances.getActiviteCCA(o)
                    if id_cca not in mapping_obs:
                        mapping_obs[id_cca] = {}
                    if (r,m,score) not in mapping_obs[id_cca]:
                        mapping_obs[id_cca][(r,m,score)] = []
                    mapping_obs[id_cca][(r,m,score)].append(o)
            mapping_obs,sequences_cca = self.filtrer(constellation,mapping_obs,failed)
            return mapping_obs,sequences_cca
        
        def nouvelleRequete(self,constellation,id_cca,r,m):
            test_exp = r in [x[0] for x in self.explications.keys()]
            precedent_absent_cca = True
            if m>0:
                for (s,o) in constellation.getRequete(r).getMode(m-1).getCouples():
                    if self.grapheDependances.getActiviteCCA(o)==id_cca:
                        precedent_absent_cca = False
                        break
            return test_exp and precedent_absent_cca
        
        def rejeterMode(self,garder,id_cca,r,m,score,mapping_obs,del_modes,sequences_cca):
            mode_rejete = 0
            activites_rejetees = 0
            if (r,m,score) in mapping_obs[id_cca]:
                for o in mapping_obs[id_cca][(r,m,score)]:
                    if o in sequences_cca[id_cca]: # on est pas sur (modes rates, modes indetermines)
                        sequences_cca[id_cca].remove(o) 
            if (r,m) in self.getModesRetenus() and (r,m) not in del_modes:
                del_modes.append((r,m))
                mode_rejete += 1
                activites_rejetees += len(mapping_obs[id_cca][(r,m,score)])
            return mode_rejete,activites_rejetees
        
        def conserverMode(self,r,m,score,mapping_obs,sequences_cca,id_cca,cca_conservees):
            modes_conserves = 0
            activites_conservees = 0
            if (r,m) in self.getModesRetenus():
                modes_conserves += 1
                activites_conservees += len(mapping_obs[id_cca][(r,m,score)])
                del mapping_obs[id_cca][(r,m,score)]
            else:
                for o in mapping_obs[id_cca][(r,m,score)]:
                    if o in sequences_cca[id_cca]: # modes indetermines (peut etre partiellement satisfaits)
                        mapping_obs[id_cca][(r,m,score)].remove(o)
                        if mapping_obs[id_cca][(r,m,score)]==[]:
                            del mapping_obs[id_cca][(r,m,score)]
            if mapping_obs[id_cca] == {}:
                id_cca = id_cca
                cca_conservees.append(id_cca)
            return modes_conserves,activites_conservees
            
        def filtrer(self,constellation,mapping_obs,failed):
                sequences_cca = {}
                del_modes = []
                modes_conserves = 0
                activites_conservees = 0
                cca_conservees = []
                modes_rejetes = 0
                activites_rejetees = 0
                
                            
                for id_cca in mapping_obs:
                    tri = list(mapping_obs[id_cca].keys())
                    modes = sorted(tri,key=lambda x : (x[0]==-1,x[2][0],-x[0]),reverse=True)
                    garder = True
                    modes_candidats = [] # ceux qu'on conserve pas besoin de redemander de les inserer
                    s,cca = id_cca
                    sequences_cca[id_cca] = self.getSolCCA(s,cca).getSequence().copy()
                    for (r,m,score) in modes:
                        if garder:
                            res = self.conserverMode(r,m,score,mapping_obs,sequences_cca,id_cca,cca_conservees)
                            modes_conserves += res[0]
                            activites_conservees += res[1]
                        else:
                            res = self.rejeterMode(garder,id_cca,r,m,score,mapping_obs,del_modes,sequences_cca)
                            modes_rejetes += res[0]
                            activites_rejetees += res[1]
                        # si une nouvelle requete arrive l'ordre change potentiellement
                        if not config.getOptValue("conserve") and ((r,m) in self.explications or self.nouvelleRequete(constellation,id_cca,r,m)):
                            garder = False
                            
                for id_cca in cca_conservees:
                    del mapping_obs[id_cca]
                for x in del_modes:
                    self.retirerModeRetenu(x)
                printColor("Récupération de "+str(activites_conservees)+" activités "+str(modes_conserves)+" modes partiels et "+str(len(cca_conservees))+" CCAs",c='y')
                printColor("Rejet de "+str(activites_rejetees)+" activités "+str(modes_rejetes)+" modes partiels",c='y')
                return mapping_obs,sequences_cca
            
        def repartirTravailSlaves(self,mapping_obs,sequences_cca):
            sequences_slaves = {}
            tailles = [(sum([len(mapping_obs[id_cca][x]) for x in mapping_obs[id_cca]]),id_cca) for id_cca in mapping_obs]
            self.mappingCCATaille(tailles)
            job_slaves = {}
            for cpu in self.cca_slaves:
                sequences_slaves[cpu] = {}
                job_slaves[cpu] = {}
                for id_cca in self.cca_slaves[cpu]:
                    if id_cca in mapping_obs:
                        sequences_slaves[cpu][id_cca] = sequences_cca[id_cca]
                        job_slaves[cpu][id_cca] = mapping_obs[id_cca]
            del_job = []
            for cpu in job_slaves:
                if job_slaves[cpu]=={}:
                    del_job.append(cpu)
            for cpu in del_job:
                del job_slaves[cpu]
            self.job_slaves = job_slaves
            return job_slaves,sequences_slaves
    
        def lectureAsynchrone(self,constellation,freeCPUs,mailbox,echecs,indetermines,sequences):
            count = 0
            for data in mailbox.readMessages():
                if data is not None:
                    count += 1
                    for cca in data['sequences']:
                        self.sequences[cca] = data['sequences'][cca]
                    exp_cpu,inde = self.analyserMessage(data,constellation,freeCPUs,mailbox)
                    for (r,m) in inde:
                        if (r,m) not in indetermines:
                            indetermines.append((r,m))
                    for (r,m) in exp_cpu:
                        if (r,m) not in echecs:
                            echecs[(r,m)] = []
                        echecs[(r,m)] += exp_cpu[(r,m)]
            return count
        
        def envoiAsynchrone(self,ccas,mapping_obs,sequences_cca,freeCPUs,mailbox,job_slaves):
            count = 0
            while ccas!=[] and freeCPUs!=[]:
                (taille,id_cca) = ccas.pop(0)
                if id_cca in mapping_obs:
                    cpu = freeCPUs.pop()
                    if cpu not in job_slaves:
                        job_slaves[cpu] = {}
                    job_slaves[cpu][id_cca] = mapping_obs[id_cca]
                    data = {"taches":{id_cca:job_slaves[cpu][id_cca]},"sequences":{id_cca:sequences_cca[id_cca]}}
                    mailbox.demanderTache(data,cpu)
                    self.solution.historique.enregistrerEnvoi(cpu)
                else:
                    count += 1
            return count
        
        def rechercheStabilite(self,constellation,r,CCAStables,CCAInstables,echecs):
            m = constellation.getRequete(r).getModeCourant(constellation).getId()
            if (r,m) in echecs:  # r a échoué
                stable = False
            else:
                stable = True
                for cca in constellation.getRequete(r).getCCAPresentes(self.grapheDependances):
                    if cca in CCAInstables:
                        stable = False
                        break
            return stable
        
        def genererExplicationRequete(self,constellation,r,echecs,explications,CCAStables):
            m = constellation.getRequete(r).getModeCourant(constellation).getId()
            if (r,m) in echecs:
                explications[(r,m)] = []
                for cca in constellation.getRequete(r).getCCAPresentes(self.grapheDependances):
                    if cca in CCAStables:
                        cca_echecs = [o for (s,o) in constellation.getRequete(r).getModeCourant(constellation).getCouples() if self.grapheDependances.getActiviteCCA(o)==cca and o in echecs[(r,m)]]
                        explications[(r,m)] += cca_echecs
                if len(explications[(r,m)])==0:
                    del explications[(r,m)]
                    self.echecs_sans_explications.append((r,m))
         
        def propagerInstabilite(self,constellation,r,CCAStables,CCAInstables):
            for cca in constellation.getRequete(r).getCCAPresentes(self.grapheDependances):
                if cca in CCAStables:
                    CCAStables.remove(cca)
                    CCAInstables.append(cca)
        
        def genererExplications(self,constellation,echecs):
            if not config.getOptValue("stable_exp"):
                self.echecs_sans_explications = []
                self.explications = echecs
            else:
                self.echecs_sans_explications = []
                explications = {}
                requetes = [r for r in constellation.getRequetes() if constellation.getRequete(r).getModeCourant(constellation) is not None]
                CCAStables = [id_cca for id_cca in self.grapheDependances.getComposantes()]
                CCAInstables = []
                requetes_triees = sorted(requetes,key = lambda r : (r==-1,self.recompenseBruitee(constellation,r,constellation.getRequete(r).getModeCourant(constellation).getId())[0],-r),reverse=True)
                for r in requetes_triees:
                    stable = self.rechercheStabilite(constellation,r,CCAStables,CCAInstables,echecs)
                    self.genererExplicationRequete(constellation,r,echecs,explications,CCAStables)
                    if not stable:
                        self.propagerInstabilite(constellation,r,CCAStables,CCAInstables)
                self.explications = explications
                    
        def lectureEnvoiSimultanes(self,constellation,mapping_obs,sequences_cca,mailbox,freeCPUs):
            t1 = time()
            tailles = [(sum([len(mapping_obs[id_cca][x]) for x in mapping_obs[id_cca]]),id_cca) for id_cca in mapping_obs]
            ccas = sorted(tailles,key=itemgetter(0),reverse=True)
            n_cca = len(ccas)
            ccas_done = 0
            job_slaves = {}
            echecs = {}
            indetermines = []
            self.sequences = {}
            while ccas_done != n_cca: # cca done par lecture des messages / cca dones car rien a faire (detection a l'envoi)
                ccas_done += self.lectureAsynchrone(constellation,freeCPUs,mailbox,echecs,indetermines,sequences_cca)
                ccas_done += self.envoiAsynchrone(ccas,mapping_obs,sequences_cca,freeCPUs,mailbox,job_slaves)
            self.job_slaves = job_slaves
            self.genererExplications(constellation,echecs)
            self.indetermines = indetermines
        
        def envoyerTacheAsynchrone(self,constellation,freeCPUs,mailbox,failed):
            mapping_obs,sequences_cca = self.mappingCCAObservationsCourantes(constellation,failed) # cca -> (r,m) -> liste obs sur la cca
            if comm.Get_size()==1:
                freeCPUs.remove(0)
                job_slaves,sequences_slaves = self.repartirTravailSlaves(mapping_obs,sequences_cca)
                data = {"time":time(),"taches":job_slaves[0],"sequences":sequences_slaves[0]}
                self.localSlave.insererActivites(constellation,data)
            else:
                self.lectureEnvoiSimultanes(constellation,mapping_obs,sequences_cca,mailbox,freeCPUs)
        def envoyerTache(self,constellation,freeCPUs,mailbox,failed):
            mapping_obs,sequences_cca = self.mappingCCAObservationsCourantes(constellation,failed) # cca -> (r,m) -> liste obs sur la cca
                # ajouter potentiellement un remapping pour repartir la charge de travail
            job_slaves,sequences_slaves = self.repartirTravailSlaves(mapping_obs,sequences_cca)
            #printColor(job_slaves,c='y')
            comm = MPI.COMM_WORLD
            #die()
            if comm.Get_size()==1:
                freeCPUs.remove(0)
                data = {"time":time(),"taches":job_slaves[0],"sequences":sequences_slaves[0]}
                self.localSlave.insererActivites(constellation,data)
            else:
                for cpu in job_slaves:
                    freeCPUs.remove(cpu)
                    data = {"taches":job_slaves[cpu],"sequences":sequences_slaves[cpu]}
                    mailbox.demanderTache(data,cpu)
                    self.solution.historique.enregistrerEnvoi(cpu)
    
    class Processus:
        def __init__(self,role):
            self.role = role
        
        def resoudre(self,constellation):
            pass
        
    class Master(Processus):
        def __init__(self,start,tlim):
            super().__init__("master")
            self.start = start
            self.tlim = tlim
            
        def resoudre(self,constellation,mailbox,CCAs,solution,modeleDeTransition,dt_construction_transition):
            # iterer insertion-réparation
            self.initSolution(constellation,CCAs,solution,modeleDeTransition,dt_construction_transition)
            self.solution.resoudre(constellation,mailbox)
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
            
        def initSolution(self,constellation,CCAs,solution,modeleDeTransition,dt_construction_transition):
            self.solution = BatchParallelCCASearch(constellation,self.start,modeleDeTransition,dt_construction_transition,tlim=self.tlim,CCAs=CCAs,solution=solution)
            
        def getSolution(self):
            return self.solution.getSolution()

        def verifierSolution(self,modeleDeTransition):
            self.solution.verifierSolution(modeleDeTransition)
            
        def getModesSolution(self):
            return self.solution.getModes()
    
    class Slave(Processus):
        # classe pour planifier les sequences des cca sur les slaves
        class PlanificateurCCA:
            def __init__(self,sol_cca,modeleDeTransition):
                self.solution = sol_cca
                self.modeleDeTransition = modeleDeTransition
    
            def __str__(self):
                return self.solution.__str__()
            
            def setSequence(self,constellation,sequence):
                self.solution.setSequence(constellation,sequence,self.modeleDeTransition)
                
            def getSatellite(self):
                return self.solution.getSatellite()
    
            def getIdentifiant(self):
                return self.solution.getIdentifiant()
        
            def getSolution(self):
                return self.solution
            
            def getSequence(self):
                return self.solution.getSequence()
            
            def solverUsed(self):
                return self.solver_use
            
            # VERSION LKH
            def planifierMode(self,constellation,activites):
                seq_avant = self.solution.sequence.copy()
                derniereSolution = deepcopy(self.solution)
                self.solution.ajouterActivites(constellation,activites,self.modeleDeTransition,cheap=False)
                seq_inter = self.solution.sequence.copy()
                #printColor(MPI.COMM_WORLD.Get_rank(),"insertion",activites,len(self.solution.sequence),c='g')
                succes = self.solution.calculerSequence(constellation,self.modeleDeTransition,'LKH')
                self.solver_use = self.solution.solverUsed()
                if not(not succes or len(seq_avant)+len(activites)==len(self.solution.sequence)):
                    print(MPI.COMM_WORLD.Get_rank(),succes,seq_avant,activites,self.solution.sequence)
                    assert(False)
                if not succes:
                    sequence = derniereSolution.getSequence().copy()
                    self.solution = derniereSolution
                return succes
            
            # VERSION OPTW
            # (r,m) => liste d'activites
            def planifierListModes(self,constellation,modes_tries):
                shiftRightDisplay(4)
                seq_avant = self.solution.sequence.copy()
                derniereSolution = deepcopy(self.solution)
                self.solution.ajouterListeModes(constellation,modes_tries,self.modeleDeTransition)
                seq_inter = self.solution.sequence.copy()
                groups = {}
                groups[0] = copy(seq_avant)
                for i,mode in enumerate(modes_tries):
                    groups[i+1] = modes_tries[mode]
                    assert(modes_tries[mode]!=[])
                #printColor(MPI.COMM_WORLD.Get_rank(),"insertion",activites,len(self.solution.sequence),c='g')
                res = self.solution.calculerSequence(constellation,self.modeleDeTransition,solver=config.getOptValue("solver_PCCAS"),groups=groups)
                succes,requetes_satisfaites = res
                self.solver_use = self.solution.solverUsed()
                assert(succes)
                """
                if not(not succes or len(seq_avant)+sum([len(modes_tries[m]) for m in modes_tries])==len(self.solution.sequence)):
                    print(MPI.COMM_WORLD.Get_rank(),succes,seq_avant,activites,self.solution.sequence)
                    assert(False)
                if not succes:
                    sequence = derniereSolution.getSequence().copy()
                    self.solution = derniereSolution
                """
                shiftLeftDisplay(4)
                return requetes_satisfaites
                    
            def sequenceFaisable(self,constellation):
                return self.solution.sequenceFaisable(constellation,self.modeleDeTransition)
            
        def __init__(self,constellation,start,modeleDeTransition,CCAs=None):
            comm = MPI.COMM_WORLD
            rank = comm.Get_rank()
            super().__init__("slave "+str(rank))
            activites = constellation.extraireActivitesRequetes()
            if CCAs is not None:
                self.grapheDependances = CCAs
            else:
                self.grapheDependances = GroupeComposantesActivitesStatiques(constellation,activites,modeleDeTransition)
            self.planificateurCCA = {}
            self.start = start
            
            for s in constellation.getSatellites():
                for cca in self.grapheDependances.getCCASatellite(s):
                    self.planificateurCCA[(s,cca)] = self.PlanificateurCCA(SolCCA(cca,s),modeleDeTransition)
        
        def estMessageFin(self,data):
            return 'fin' in data
        
        def etendreDureesCalcul(self):
            durees = []
            for i,date in enumerate(self.dates_cca):
                if i<len(self.dates_cca)-1:
                    durees.append((date,self.dates_cca[i+1]-date))
            return durees
        
        def packSolution(self):
            return {'indetermines':self.etat_inconnus,'sequences':self.sequences,'explications':self.explications,'cpu_time':self.etendreDureesCalcul(),'reception':time(),'envoi':(self.date_envoi,self.duree_envoi),'source':rank,'solver time':self.solver_time}
        
        
        def insererActivitesOPTW(self,constellation,data):
            self.duree_envoi = time() - data['time']
            self.start_calcul = time()
            self.date_envoi = data['time']
            del data['time']
            start = time()
            for id_cca in self.planificateurCCA: 
                self.planificateurCCA[id_cca].setSequence(constellation,data['sequences'].get(id_cca,[]))
            self.sequences = {}
            self.explications = {}
            self.etat_inconnus = []
            self.dates_cca = [time()]
            self.solver_time = 0
            for id_cca in data['taches']:
                modes = sorted(data['taches'][id_cca],key=lambda x : (x[0]==-1,x[2][0],-x[0]),reverse=True)
                activites_modes_tries = {(x[0],x[1]) : data['taches'][id_cca][x] for x in modes}
                start_solver = time()
                requetes_satisfaites = self.planificateurCCA[id_cca].planifierListModes(constellation,activites_modes_tries)
                if self.planificateurCCA[id_cca].solverUsed():
                        self.solver_time += time()-start_solver
                for (r,m,score) in data['taches'][id_cca]:
                    if r not in requetes_satisfaites:
                        activites = data['taches'][id_cca][(r,m,score)]
                        if (r,m) not in self.explications:
                            self.explications[(r,m)] = []
                        self.explications[(r,m)] += activites
                
                self.sequences[id_cca] = self.planificateurCCA[id_cca].solution.getSequence()
                self.dates_cca.append(time())
                # trier les modes et inserer leurs activites
                
            self.resultat = self.packSolution()        
            
        def insererActivites(self,constellation,data):
            self.duree_envoi = time() - data['time']
            self.start_calcul = time()
            self.date_envoi = data['time']
            del data['time']
            start = time()
            for cca in self.planificateurCCA: # le bug est surement ici : pb de cle
                self.planificateurCCA[cca].setSequence(constellation,data['sequences'].get(cca,[]))
            nActivites = 0
            nModes = 0
            nCCA = 0
            self.sequences = {}
            self.explications = {}
            self.etat_inconnus = []
            self.dates_cca = [time()]
            self.solver_time = 0
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
                        start_solver = time()
                        succes = self.planificateurCCA[cca].planifierMode(constellation,activites)
                        if self.planificateurCCA[cca].solverUsed():
                            self.solver_time += time()-start_solver
                        if not succes:
                            fails += 1
                            n_calls += 1
                            #printColor("CPU ",MPI.COMM_WORLD.Get_rank(),"echec insertion",(r,m),activites,"dans la cca",cca,c='r')
                            if (r,m) not in self.explications:
                                self.explications[(r,m)] = []
                            self.explications[(r,m)] += activites
                            if fails >= config.getOptValue("fails") or n_calls >= config.getOptValue("scalls"):
                                stop = True
                        else:
                            fails = 0
                    else:
                        self.etat_inconnus.append((r,m))
                self.sequences[cca] = self.planificateurCCA[cca].solution.getSequence()
                #print("CCA",cca,self.planificateurCCA[cca].solution.planEarliest(constellation,self.sequences[cca]))
                
                self.dates_cca.append(time())
                # trier les modes et inserer leurs activites
            self.resultat = self.packSolution()
            #printColor("Processus",MPI.COMM_WORLD.Get_rank(),"terminé (",nActivites,"activités",nModes,"modes",nCCA,"CCAs )",c='b')
        
        def resoudre(self,constellation,mailbox):
            comm = MPI.COMM_WORLD
            while True:
                data = comm.recv()
                if self.estMessageFin(data):
                    comm.Barrier()
                    #MPI.COMM_WORLD.send({cca : self.planificateurCCA[cca].getSequence() for cca in self.planificateurCCA},dest=0)
                    break
                else:
                    if config.getOptValue("solver_PCCAS") == 'LKH':
                        self.insererActivites(constellation,data)
                    else:
                        self.insererActivitesOPTW(constellation,data)
                    mailbox.posterMessage(self.resultat)
        
        def getResultat(self):
            if MPI.COMM_WORLD.Get_size()==1:
                self.resultat['reception'] = (self.resultat['reception'],time() - self.resultat['reception'])
            return self.resultat
        
        
        
    path = '../data'
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    
    
    class Messagerie:
        def __init__(self):
            # eviter d'envoyer les informations qu'on va recuperer à l'identique ensuite
            self.local_data = {}
        
        def envoyerMessage(self,data,cpu):
            comm = MPI.COMM_WORLD
            rank = comm.Get_rank()
            assert(rank==0)
            data['time']=time()
            #self.extraireDonneesRedondantes(data,cpu)
            comm.send(data, dest=cpu, tag=rank)
        
        def extraireDonneesRedondantes(self,data,cpu):
            # {'mode':mode,'cca':cca,'source':rank,'activites':activites,'time':time()}
            id_cca = data['cca'].getIdentifiant()
            self.local_data[cpu] = {}
            self.local_data[cpu]['sequence'] = data['cca'].getSequence().copy()
            self.local_data[cpu]["mode"] = data['mode']
            self.local_data[cpu]['activites'] = data['activites']
            self.local_data[cpu]['cca'] = id_cca
            del data['mode']
            data['cca'] = id_cca
            
        def reformerDonnees(self,data):
            source = data['source']
            data['cca'] = self.local_data[source]['cca']
            data['mode'] = self.local_data[source]['mode']
            data['sequence'] = self.local_data[source]['sequence']
            data['explication'] = self.local_data[source]['activites']
            if data['faisable']:
                data['sequence'] += data['explication']
                
    class MessagerieMessageUnique(Messagerie):
        def __init__(self):
            pass
        
        def posterMessage(self,data):
            assert(MPI.COMM_WORLD.Get_rank()>0)
            MPI.COMM_WORLD.send(data,dest=0)
            
        def readMessages(self):
            assert(MPI.COMM_WORLD.Get_rank()==0)
            data = MPI.COMM_WORLD.recv()
            data['reception'] = (data['reception'],time()-data['reception'])
            return [data]
        
        def demanderTache(self,data,cpu):
            assert(MPI.COMM_WORLD.Get_rank()==0)
            envoyerMessage(data,cpu)
            
        def __str__(self):
            return "Messagerie à message unique"
        
    class MessagerieSynchronisee(Messagerie):
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
        def demanderTache(self,data,cpu):
            assert(MPI.COMM_WORLD.Get_rank()==0)
            self.cpu_mapped[cpu-1] = 1
            self.envoyerMessage(data,cpu)
            
        # appelé par les slaves
        def posterMessage(self,data):
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
            return "Messagerie synchronisée"
            
    class MessageriePartagee(Messagerie):
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
            self.flag_messages, itemsize = self.win.Shared_query(0)
            for i in range(self.size):
                self.flag_messages[i] = False
            #self.flag_messages = np.ndarray(buffer=self.buf, dtype=bool, shape=(size,))
        
        def demanderTache(self,data,cpu):
            assert(MPI.COMM_WORLD.Get_rank()==0)
            self.envoyerMessage(data,cpu)
    
                
        def posterMessage(self,data):
            assert(MPI.COMM_WORLD.Get_rank()>0)
            slot = MPI.COMM_WORLD.Get_rank()-1
            self.flag_messages[slot] = True
            MPI.COMM_WORLD.send(data,dest=0)
        
        def readMessages(self):
            assert(MPI.COMM_WORLD.Get_rank()==0)
            for i,x in enumerate(self.flag_messages):
                if x:
                    process = i+1
                    data = MPI.COMM_WORLD.recv(source=process)
                    data['reception'] = (data['reception'],time()-data['reception'])
                    #self.reformerDonnees(data)
                    assert(process==data['source'])
                    self.flag_messages[i] = False
                    yield data
                else:
                    yield None
        
        def __str__(self):
            return "Messagerie partagée : "+ str([self.flag_messages[i] for i in range(self.size)])
    
    class runnableBPCCAS:
        def execute(self,constellation,start_date,modeleDeTransition,dt_construction_transition,tlim=np.inf,CCAs=None,solution=None):
            mailbox = MessageriePartagee()
            id_mpi = MPI.COMM_WORLD.Get_rank()
            seed_option = config.getOptValue("seed")
            rd.seed(id_mpi+seed_option)
            if MPI.COMM_WORLD.Get_rank()==0:
                self.process = Master(start_date,tlim)
                self.process.resoudre(constellation,mailbox,CCAs,solution,modeleDeTransition,dt_construction_transition)
                
                for i in range(1,MPI.COMM_WORLD.Get_size()):
                    MPI.COMM_WORLD.send({"sol":self.process.solution.getSolutionContainer()},dest=i)
                
                return self.process.solution.getSolutionContainer()
            else:
                process = Slave(constellation,start_date,modeleDeTransition,CCAs=CCAs)
                process.resoudre(constellation,mailbox)
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
        
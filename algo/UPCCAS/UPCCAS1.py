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
    from ..model.composantes import *
    from ..model.solution_composantes import *
    from ..model.graphes import *
    
    from time import time
    from time import sleep
    
    class EtatRecherche:
        def __init__(self):
            self.activities_done = {}
            self.activities_working = {}
            self.activities_to_do = {}
            self.working_modes = []
            self.etat_mode = {}
            self.requetes_echouees = []
            self.WORKING = 0
            self.DONE = 1
            self.FAILED = 2
            self.WAITING_RESTART = 3
            self.WAITING_VALIDATION = 4
            self.started_modes = []
            self.waiting_modes = []
            self.cca_occupees = []
            self.fautifs = {} # (r,m) ceux qui le bloquent
            self.a_redemarrer = {}
            self.a_supprimer = []
            self.bloquants = {} # fautifs mais dans l'autre sens : (r,m) ceux que (r,m) bloque
            self.waiting_restart = {}    
                
        def libererCCA(self,id_cca):
            self.cca_occupees.remove(id_cca)
            
        def notifierEchecRequete(self,r):
            if r not in self.requetes_echouees:
                self.requetes_echouees.append(r)
        
        def mettreEnAttenteValidation(self,mode):
            self.etat_mode[mode] = self.WAITING_VALIDATION
            if mode not in self.waiting_modes:
                self.waiting_modes.append(mode)
        
        # liste_activites : à indiquer comme deja en travail (pour le transfert de mode)
        def preWorking(self,r,m,liste_activites):
            for a in liste_activites:
                self.activities_to_do[(r,m)].remove(a)
                self.activities_working[(r,m)].append(a)
                
        def getModesEnAttenteValidation(self):
            return self.waiting_modes
        
        #def getModesEnAttenteRedemarrage(self):
        #    return self.waiting_restart
        
        def mettreEnAttenteRedemarrage(self,mode):
            self.etat_mode[mode] = self.WAITING_RESTART
        
        def getWorkingModes(self):
            return self.working_modes
        
        def chercherElement(self,liste,elmt):
            for i,x in enumerate(liste):
                if x==elmt:
                    return i
            return -1
        
        def activitesAFaire(self,r,m):
            return self.activities_to_do[(r,m)]!=[]
        
        def getStartedModes(self):
            return self.started_modes
        
        def getASupprimer(self):
            return self.a_supprimer
        
        def notifierSuppression(self,s,o):
            self.a_supprimer.remove((s,o))
            
        # suppression du mode pour la version avec anticipation
        def supprimerAnticipation(self,constellation,r,m,grapheDep,liste_a_supprimer=None):
            self.terminerMode(r,m,False)
            for (s,o) in constellation.getRequete(r).getMode(m).getCouples():
                if liste_a_supprimer is None or o in liste_a_supprimer:
                    if((s,o) not in self.a_supprimer):
                        id_cca = grapheDep.getActiviteCCA(o)
                        if id_cca in self.cca_occupees:
                            self.a_supprimer.append((s,o))
            
        # annule un mode et renvoie la liste de ceux a redemarrer. SUPRESSION DU MODE
        def echecCertain(self,constellation,r,m,grapheDep,liste_a_supprimer=None):
            self.supprimerAnticipation(constellation,r,m,grapheDep,liste_a_supprimer)
            if not config.getOptValue("anticipation"):
                return []
            # MAJ des requetes fautives/activites a redemarrer
            if r in self.fautifs:
                for fautif in self.fautifs[r]:
                    self.bloquants[fautif].remove(r)
                del self.fautifs[r]
            if (r,m) in self.a_redemarrer:
                del self.a_redemarrer[(r,m)]
            if r in self.bloquants:
                redemarrer = self.bloquants[r].copy()
                for bloc in self.bloquants[r]:
                    self.fautifs[bloc].remove(r)
                del self.bloquants[r]
            else:
                return []
            #clean les modes a redemarrer : eviter un 2e appel au redemarrage
            for r in redemarrer:
                for fautif in self.fautifs[r]:
                    self.bloquants[fautif].remove(r)
                if r in self.fautifs:
                    del self.fautifs[r]
            return redemarrer 
        
        # indiquer un echec d'activités causé probablement par des requetes fautives.
        def echecIncertain(self,r,m,fautifs,activites_a_redemarrer,attenteRedemarrage=True):
            if len(fautifs)==0:
                assert(not attenteRedemarrage)
                printOpen('Echec du mode',(r,m),'cas 1 (transfert) : redémarrer ',activites_a_redemarrer,c='m')
            else:
                assert(attenteRedemarrage)
                printOpen('Echec du mode',(r,m),'cas 2 (mode non prioritaire) : redémarrer',activites_a_redemarrer,c='m')
            #printColor("Reprogrammer",(r,m),activites_a_redemarrer,c='c')
            assert(config.getOptValue("anticipation"))
            for requete in fautifs:
                if r not in self.fautifs:
                    self.fautifs[r] = []
                if requete not in self.fautifs[r]:
                    self.fautifs[r].append(requete)
                if requete not in self.bloquants:
                    self.bloquants[requete] = []
                if r not in self.bloquants[requete]:
                    self.bloquants[requete].append(r)
            for a in activites_a_redemarrer:
                if (r,m) not in self.a_redemarrer:
                    self.a_redemarrer[(r,m)] = []
                if a not in self.a_redemarrer[(r,m)]:
                    self.a_redemarrer[(r,m)].append(a)
                    self.activities_working[(r,m)].remove(a)
            self.mettreEnAttenteRedemarrage((r,m))
            if not attenteRedemarrage:
                self.redemarrerMode(r,m)
            printClose()
            
        def terminerMode(self,r,m,succes):
            if succes:
                self.etat_mode[(r,m)] = self.DONE
            else:
                assert(self.etat_mode[(r,m)]) != self.FAILED
                self.etat_mode[(r,m)] = self.FAILED
            if (r,m) in self.activities_done:
                del self.activities_done[(r,m)]
            if (r,m) in self.activities_to_do:
                del self.activities_to_do[(r,m)]
            if (r,m) in self.activities_working:
                del self.activities_working[(r,m)]
            if (r,m) in self.started_modes:
                self.started_modes.remove((r,m))
            if (r,m) in self.waiting_modes:
                self.waiting_modes.remove((r,m))
            self.working_modes.remove((r,m))
        
        def ccaToDo(self,r,m,id_cca,grapheDep):
            for a in self.activities_to_do[(r,m)]:
                to_do_cca = grapheDep.getActiviteCCA(a)
                if to_do_cca == id_cca:
                    return True
            return False
        """
                features
        """
        def filtrerModesEnCours(self,parents):
            return [x for x in parents if x in self.activities_to_do and self.activities_to_do[x]!=[] and self.etat_mode[x] == self.WORKING]
        
        def modesEnCours(self):
            return len(self.working_modes)>0
         
        def ccaOccupee(self,id_cca):
            return id_cca in self.cca_occupees
        
        def getActivitesRestantes(self,r,m):
            return self.activities_to_do[(r,m)]
        
        def getCCALibres(self,r,m,grapheDep):
            ccas = [grapheDep.getActiviteCCA(a) for a in self.activities_to_do[(r,m)]]
            ccas_libre = [id_cca for id_cca in ccas if id_cca not in self.cca_occupees]
            return ccas_libre
        
        def getCCA(self,r,m,grapheDep,id_cca):
            if (r,m) not in self.started_modes:
                self.started_modes.append((r,m))
            assert(id_cca not in self.cca_occupees)
            act = [a for a in self.activities_to_do[(r,m)] if grapheDep.getActiviteCCA(a)==id_cca]
            supprimerListeElements(self.activities_to_do[(r,m)],act)
            self.cca_occupees.append(id_cca)
            for a in act:
                self.activities_working[(r,m)].append(a)
            return act
            
        # depile les activites et les renvoie
        def getCCASuivante(self,r,m,grapheDep):
            assert(len(self.activities_to_do[(r,m)])>0)
            if (r,m) not in self.started_modes:
                self.started_modes.append((r,m))
            ccas = [grapheDep.getActiviteCCA(a) for a in self.activities_to_do[(r,m)]]
            ccas_libre = [id_cca for id_cca in ccas if id_cca not in self.cca_occupees]
            id_cca = ccas_libre[0]
            self.cca_occupees.append(id_cca)
            liste_a = []
            act = [a for a in self.activities_to_do[(r,m)] if grapheDep.getActiviteCCA(a)==id_cca]
            supprimerListeElements(self.activities_to_do[(r,m)],act)
            for a in act:
                assert(a not in self.activities_to_do[(r,m)])
            for a in act:
                assert (a not in self.activities_to_do[(r,m)])
            for a in act:
                self.activities_working[(r,m)].append(a)
            return act
        
        def changementNomCCA(self,cca_avant,cca_apres):
            if cca_avant in self.cca_occupees:
                self.cca_occupees.remove(cca_avant)
                self.cca_occupees.append(cca_apres)
                
        def divisionCCA(self,cca_avant,liste_cca):
            if cca_avant in self.cca_occupees:
                self.cca_occupees.remove(cca_avant)
                for id_cca in liste_cca:
                    self.cca_occupees.append(id_cca)
        
        def requeteEchouee(self,r):
            return r in self.requetes_echouees
        
        def getEtat(self,r,m):
            return self.etat_mode[(r,m)]
        
        # ajoute un mode directement actif
        def ajouterModeActif(self,r,m,activites):
            self.activities_done[(r,m)] = []
            self.activities_working[(r,m)] = []
            self.activities_to_do[(r,m)] = activites
            self.working_modes.append((r,m))
            self.etat_mode[(r,m)] = self.WORKING
        
        def redemarrerMode(self,r,m):
            assert(self.etat_mode[(r,m)] in [self.WAITING_RESTART,self.WORKING])
            #self.activities_done[(r,m)] = []
            #self.activities_working[(r,m)] = []
            self.activities_to_do[(r,m)] += self.a_redemarrer[(r,m)]
            if (r,m) in self.waiting_modes:
                self.waiting_modes.remove((r,m))
            self.etat_mode[(r,m)] = self.WORKING
            del self.a_redemarrer[(r,m)]
            #printColor('mode redémarré : ',(r,m),'restantes :',self.activities_to_do[(r,m)],'réussies :',self.activities_done[(r,m)],'en cours :',self.activities_working[(r,m)],c='y')
            #self.working_modes.append((r,m))
        
        # si anticiper_modes : renvoie la liste de REQUETES a redemarrer
        def validerMode(self,r,m):
            self.terminerMode(r,m,True)
            if(config.getOptValue("anticipation")):
                requetes_a_annuler = []
                if r in self.bloquants:
                    for requete in self.bloquants[r]:
                        self.fautifs[requete].remove(r)
                        if self.fautifs[requete] == []:
                            requetes_a_annuler.append(requete)
                    del self.bloquants[r]
                return requetes_a_annuler
            else:
                return []
        
        def prevaliderActivites(self,activites,r,m):
            self.started_modes.append((r,m))
            printColor('Prévalidation des activités',activites,"du mode",(r,m),c='m')
            supprimerListeElements(self.activities_to_do[(r,m)],activites)
            if (r,m) not in self.activities_done:
                self.activities_done[(r,m)] = []
            self.activities_done[(r,m)] += activites
            if len(self.activities_to_do[(r,m)])==0 and len(self.activities_working[(r,m)])==0:
                return True
            else:
                return False
            
        def validerActivites(self,activites,r,m):
            printColor('Validation des activités',activites,"du mode",(r,m),c='c')
            supprimerListeElements(self.activities_working[(r,m)],activites)
            if (r,m) not in self.activities_done:
                self.activities_done[(r,m)] = []
            self.activities_done[(r,m)] += activites
            if len(self.activities_to_do[(r,m)])==0 and len(self.activities_working[(r,m)])==0:
                return True
            else:
                return False
    
        def annulerMode(self,r,m):
            self.terminerMode(r,m,False)
    
    
    
    class UnitParallelCCASearch(Solver):
        def __init__(self,constellation,start_date,modeleDeTransition,dt_construction_transition,tlim,CCAs,solution):
            super().__init__(constellation,start_date,modeleDeTransition,dt_construction_transition,tlim=tlim,CCAs=CCAs,solution=solution)
            # si 1 coeur : slave local
            comm = MPI.COMM_WORLD
            if comm.Get_size()==1:
                self.localSlave = Slave(constellation,modeleDeTransition)
            printColor("Durée d'initialisation :",time()-self.start_date,c='b')
      
        def _ajouterModeActifEtatRecherche(self,constellation,r,m,prevalidation=[]):
            activites = []
            for (s,o) in constellation.getRequete(r).getMode(m).getCouples():
                if o not in activites:
                    activites.append(o)
            self.etat_recherche.ajouterModeActif(r,m,activites) 
            if len(prevalidation)>0:
                self.etat_recherche.prevaliderActivites(prevalidation,r,m)
                if not self.etat_recherche.activitesAFaire(r,m):
                    self.validerMode(constellation,(r,m))
    
        def ajouterModeActif(self,constellation,r,m,prevalidation=[]):
            self._ajouterModeActifEtatRecherche(constellation,r,m,prevalidation)
            
        def presentCCA(self,constellation,r,m,id_cca,grapheActivites):
            for s in constellation.getRequete(r).getMode(m).getActivites():
                for o in constellation.getRequete(r).getMode(m).getActivitesSatellite(s):
                    if grapheActivites.getNoeud(o).getComposanteConnexe()==id_cca:
                        return True
            return False
        
        def getComposanteActivite(self,a):
            return self.grapheDependances.getActiviteCCA(a)
        
        def setCCA(self,constellation,s,seq,id_cca,r,m):
            s,cca = id_cca
            self.solution.getSolCCA(s,cca).setSequence(constellation,seq,self.modeleDeTransition)
    
        def validerActivites(self,constellation,activites,mode):
            s = constellation.getSatelliteActivite(activites[0])
            act = {a : None for a in activites}
            id_cca = self.getComposanteActivite(activites[0])
            seq = []
            for i,a in enumerate(activites):
                cca_a = self.getComposanteActivite(a)
                if cca_a!=id_cca:
                    self.setCCA(constellation,s,seq,id_cca,mode[0],mode[1])
                    id_cca = cca_a
                    seq = []
                seq.append(a)
                if i==len(activites)-1:
                    self.setCCA(constellation,s,seq,cca_a,mode[0],mode[1])
            
        def MAJPresenceRequetes(self,ancienne_cca,nouvelles_ccas):
            # MAJ sol CCA
            # MAJ Requete CCAs
            requetes = self.requetesCCA[ancienne_cca]
            del self.requetesCCA[ancienne_cca]
            cca_req = {r : [] for (p,w,r) in requetes}
            rec = lambda r : constellation.getRequete(r).getMode(self.modes_courants[r]).getRecompense() + self.bruit_requete[r]
            prio = lambda r : constellation.getRequete(r).getPriorite()
            for id_cca in nouvelles_ccas:
                req = []
                for (p,w,r) in requetes:
                    for s in constellation.getRequete(r).getActivites():
                        for o in constellation.getRequete(r).getActivites()[s]:
                            if o in self.grapheDependances.getComposanteConnexe(id_cca):
                                reverse_insort(req,(prio(r),rec(r),r))
                                cca_req[r].append(id_cca)
                                break
                self.requetesCCA[id_cca] = req
            # MAJ CCA Requetes
            for r in cca_req:
                self.ccaRequetes[r].remove(ancienne_cca)
                self.ccaRequetes[r] += cca_req[r]
            
            for id_cca in nouvelles_ccas:
                if self.requetesCCA[id_cca]==[]:
                    del self.requetesCCA[id_cca]
                
            
        """
            =============================================== 
                            RESOLUTION
            =============================================== 
        """
       
        def validerMode(self,constellation,mode):
            printOpen("Validation du mode "+str(mode),c='g')
            start_val = time()
            # retirer la requete de toutes ses cca
            for id_cca in self.ccaRequetes[mode[0]]:
                self.retraitRequeteCCA(id_cca,mode[0])
            del self.ccaRequetes[mode[0]]
            
            # retirer les obs de la requete non presentes dans le mode
            retrait = []
            for s in constellation.getRequete(mode[0]).getActivites():
                for o in constellation.getRequete(mode[0]).getActivitesSatellite(s):
                    presentGraphe = o in self.grapheDependances.getNoeuds() # explications retirees...
                    nonPresentMode = o not in constellation.getRequete(mode[0]).getMode(mode[1]).getActivitesSatellite(s)
                    if presentGraphe and nonPresentMode:
                        retrait.append(o)
            del self.modes_courants[mode[0]]
            self.solution.historique.enregistrerRequeteValidee(constellation,mode[0])
            requetes_a_annuler = self.etat_recherche.validerMode(mode[0],mode[1])
            self.solution.ajouterModeRetenu(mode,constellation)

            if len(requetes_a_annuler)>0:
                printColor("Annuler les requêtes",requetes_a_annuler,c='m')
            for r in requetes_a_annuler:
                mode = (r,self.modes_courants[r])
                explication = self.explication[mode]
                self.planifierModeSuivant(constellation,mode,explication)
                del self.explication[mode]
            printClose()
        
        def modeCourant(self,mode):
            r_courant = mode[0]
            return r_courant in self.modes_courants and mode[1] == self.modes_courants[r_courant]
        
        def testSuccesATransferer(self,mode,succes):
            return config.getOptValue("preval") and not self.modeCourant(mode) and succes
        
        def testEchecsATransferer(self,mode,succes):
            return not succes and not self.modeCourant(mode)
            
        def analyserMessage(self,data,constellation,freeCPUs,mailbox,iteration):
            start = time()
            self.solution.historique.enregistrerReceptionMessage(data)
            i = data['source']
            assert(i not in freeCPUs)
            bisect.insort(freeCPUs,i)
            it = data['iteration']
            if iteration==it:
                cpu_time,cca,succes,mode = data['cpu_time'],data['cca'],data['faisable'],data['mode']
                explication = data['explication']
                seq = cca.getSequence()
                # si message d'un mode passé : transfert potentiel d'activités validées ou à redemarrer
                if self.testSuccesATransferer(mode,succes):
                    mode = self.transfertSuccesActivites(explication,seq,mode)
                if self.testEchecsATransferer(mode,succes):
                    retraits = self.supprimerAnciennesActivites(seq) # pour version anticipee de UPCCAS
                    mode,activites_a_redemarrer = self.transfertEchecActivites(constellation,explication,retraits,mode)
                else:
                    activites_a_redemarrer = []
                # indiquer la cca libre et supprimer les activités n'étant plus à programmer
                self.libererCCA(seq,explication)
                if self.modeCourant(mode):
                    # traiter le message d'un mode en cours
                    if succes and not self.etat_recherche.requeteEchouee(mode[0]) and self.etat_recherche.getEtat(mode[0],mode[1])!=self.etat_recherche.FAILED:
                        stat_validation = time()
                        if(len(seq)==0):
                            print(mode)
                        assert(len(seq)>0)
                        self.validerActivites(constellation,seq,mode)
                        start_validation = time()
                        termine = self.etat_recherche.validerActivites(explication,mode[0],mode[1])
                        self.solution.historique.enregistrerSucces(i,cpu_time)

                        if termine:
                            if config.getOptValue("anticipation"):
                                self.validerAnticipation(constellation,mode)
                            else:
                                self.validerMode(constellation,mode)
                    else:
                        if self.etat_recherche.getEtat(mode[0],mode[1])!=self.etat_recherche.FAILED:
                            if config.getOptValue("anticipation"):
                                self.notifierEchecAnticipation(constellation,mode,explication,activites_a_redemarrer)
                            else:
                                self.planifierModeSuivant(constellation,mode,explication)
                        self.solution.historique.enregistrerEchec(i,cpu_time)
                    #die(self.solution)
                        
        def analyserResultats(self,constellation,freeCPUs,mailbox,iteration):
            comm = MPI.COMM_WORLD
            rank = comm.Get_rank()
            start_iteration = time()
            if comm.Get_size()==1:
                messages = [self.localSlave.getResultat()]
            else:
                messages = mailbox.readMessages()
            for data in messages:
                if data is not None:
                    self.analyserMessage(data,constellation,freeCPUs,mailbox,iteration)
                    
        def supprimerAnciennesActivites(self,seq):
            supp = []
            for (s,o) in self.etat_recherche.getASupprimer():
                if o in seq:
                    seq.remove(o)
                    supp.append(o)
                    self.etat_recherche.notifierSuppression(s,o)
            if len(supp)>0:
                printColor('suppression dans la séquence reçue',supp,"(modes supprimés)",c='r')
            return supp
        
        def transfertEchecActivites(self,constellation,explication,retrait,mode):
            if mode[0] in self.modes_courants:
                nouveau_mode = self.modes_courants[mode[0]]
                if sum([a in [x[1] for x in constellation.getRequete(mode[0]).getMode(nouveau_mode).getCouples()] for a in explication])==0:
                     return mode,[]
                if nouveau_mode!=mode[1]:
                    a_redemarrer = explication.copy()
                    transfert = []
                    for a in retrait:
                        if a in a_redemarrer:
                            a_redemarrer.remove(a)
                    printColor('mode',mode,': transfert d\'activités échouées',transfert,'- redémarrer',a_redemarrer,c='r')
                    return (mode[0],nouveau_mode),a_redemarrer
                return mode,[]
            return mode,[]
        
        # voir si la planification d'une activité de l'ancien mode peut etre conservée pour le nouveau
        def transfertSuccesActivites(self,explication,seq,mode):
            if mode[0] in self.modes_courants:
                nouveau_mode = self.modes_courants[mode[0]]
                if nouveau_mode!=mode[1]:
                    find = False
                    for new in explication:
                        if new in self.etat_recherche.activities_working[(mode[0],nouveau_mode)]:
                            find = True
                            break # suppression deja effectuée au dessus
                    if find:
                        printColor('mode',mode,": transfert d\'activités réussies",c='g')
                        return (mode[0],nouveau_mode)
                    else:
                        return mode
                else:
                    return mode
            return mode
                
        def validerAnticipation(self,constellation,mode):
            self.etat_recherche.mettreEnAttenteValidation(mode)
            self.debloquerModesEnAttente(constellation)    
        
        # Traiter un echec : considerer un redemarrage si le mode est prioritaire : il n'y a pas de fautif surclassé
        # activites a redemarrer : car l'echec est acompagnée de la suppression d'autres activités.
        # => on leur redonne alors une chance
        def notifierEchecAnticipation(self,constellation,mode,explication,activites_a_redemarrer):
            # I) si on a des activites a redemarrer : on les redemarre et ne cherche pas d'autre mode fautif
            if len(activites_a_redemarrer): # pas d'attente de redemarrage : on redemarre de suite
                self.etat_recherche.echecIncertain(mode[0],mode[1],[],activites_a_redemarrer,False)
                return
            
            # II) si on a pas d'activités à redémarrer de l'ancien mode : on cherche des modes fautifs et on redemmarre
            # cca des activites calculees
            ccas = []
            for a in explication:
                id_cca = self.grapheDependances.getActiviteCCA(a)
                if id_cca not in ccas:
                    ccas.append(id_cca)
            # detection de fautifs : en attente s'il y en a
            premierPartout = True
            for id_cca in ccas:
                if not self.premierCCA(mode[0],id_cca):
                    premierPartout = False
                    self.explication[mode] = self.explication.get(mode,[])
                    for a in explication:
                        if a not in self.explication[mode]:
                            self.explication[mode].append(a)
                    fautifs = self.findFautifsCCA(mode[0],id_cca)
                    self.etat_recherche.echecIncertain(mode[0],mode[1],fautifs,explication,True)
                    break
                    
            # III) si pas de fautif : annuler car le mode ne passera jamais
            if premierPartout:
                self.planifierModeSuivant(constellation,mode,explication)
                self.explication[mode] = []
        
        def findFautifsCCA(self,r,id_cca):
            i = 0
            while self.requetesCCA[id_cca][i][2]!=r:
                i += 1
            return [x[2] for x in self.requetesCCA[id_cca][:i]]
        
        def premierCCA(self,r,id_cca):
            return r==self.requetesCCA[id_cca][0][2]
        
        # libere la cca et supprimer de la séquence les activités des modes supprimés
        def libererCCA(self,seq,explication):
            ccas = []
            for a in seq+explication:
                i_cca = self.grapheDependances.getActiviteCCA(a)
                if i_cca not in ccas:
                    ccas.append(i_cca)
                    self.etat_recherche.libererCCA(i_cca)
                    
        def getModesBloques(self,r):
            if r in self.requetesCCA:
                bloques = []
                for id_cca in self.requetesCCA[r]:
                    for (p,w,r) in self.requetesCCA:
                        m = self.modes_courants[r]
                        if (r,m) in self.etat_recherche.getModesEnAttente():
                            bloques.append((r,m))
                return bloques
            else:
                return []
        
        def debloquerModesEnAttente(self,constellation):
           
            validation = True
            while validation:
                validation = False
                for mode_attente in self.etat_recherche.getModesEnAttenteValidation():
                    if self.requetePrioritaire(mode_attente[0]):
                        validation = True
                        self.validerMode(constellation,mode_attente)
                        if not(mode_attente not in self.etat_recherche.getModesEnAttenteValidation()):
                            die(mode_attente,self.etat_recherche.getModesEnAttenteValidation())
                        if not(mode_attente not in self.etat_recherche.getWorkingModes()):
                            die(mode_attente,self.etat_recherche.getWorkingModes())
                        assert(self.etat_recherche.getEtat(mode_attente[0],mode_attente[1])==self.etat_recherche.DONE)
                        break
             
        def intersectionActivites(self,constellation,r):
            nouveau_mode = self.modes_courants.get(r)
            ancien_mode = nouveau_mode - 1
            act_ancien = [x[1] for x in constellation.getRequete(r).getMode(ancien_mode).getCouples()]
            if nouveau_mode is None:
                return act_ancien,[],[]
            act_nouveau_mode = [x[1] for x in constellation.getRequete(r).getMode(nouveau_mode).getCouples()]
            retrait = []
            prevalider = []
            intersection_working = []
            for o in act_nouveau_mode:
                if o in act_ancien and o in self.etat_recherche.activities_done.get((r,ancien_mode),[]):
                    if config.getOptValue("preval"):
                        prevalider.append(o)
                    else:
                        retrait.append(o)
                elif o in self.etat_recherche.activities_working.get((r,ancien_mode),[]):
                    if config.getOptValue("preval"):
                        intersection_working.append(o)
            for o in act_ancien:
                if o not in act_nouveau_mode:
                    retrait.append(o)
            return retrait,prevalider,intersection_working
        
        def planifierModeSuivant(self,constellation,mode,explication):
            exp_cca = [(o,self.grapheDependances.getActiviteCCA(o)) for o in explication]
            printOpen('Echec du mode',mode,'cas 3 : annuler. cause : [activité, (sat_cca,n_cca) ... ] = ',exp_cca,c='r')
            self.solution.historique.enregistrerEchecMode(constellation,mode[0])
            for id_cca in self.ccaRequetes[mode[0]]:
                self.retraitRequeteCCA(id_cca,mode[0])
            modes_a_redemarrer = self.defilerEtatRequete(constellation,mode,explication,self.grapheDependances) # MAJ etat recherche MAJ graphe (pending)
            self.insererRequeteCCA(constellation,mode[0])
            if len(modes_a_redemarrer)>0:
                printColor('Redémarrer les requêtes',modes_a_redemarrer,c='m')
            if config.getOptValue("anticipation"):
                for mode in modes_a_redemarrer:
                    self.etat_recherche.redemarrerMode(mode,self.modes_courants[mode])
            printClose()
    
        def insererRequeteCCA(self,constellation,r):
            ccas = []
            if r in self.modes_courants:
                for (s,o) in constellation.getRequete(r).getCouples():
                    try:
                        id_cca = self.grapheDependances.getActiviteCCA(o)
                        if id_cca not in ccas:
                            ccas.append(id_cca)
                    except:
                        pass
                self.ccaRequetes[r] = ccas
                prio = constellation.getRequete(r).getPriorite()
                rec = self.recompenseBruitee(constellation,r,self.modes_courants[r])
                for id_cca in ccas:
                    if id_cca not in self.requetesCCA:
                        self.requetesCCA[id_cca] = []
                    reverse_insort(self.requetesCCA[id_cca],(prio,rec,r))
        
        def requetePrioritaireCCA(self,r,id_cca):
            return r==self.requetesCCA[id_cca][0][2]
        
        def requetePrioritaire(self,r):
            for id_cca in self.ccaRequetes[r]:
                if not self.requetePrioritaireCCA(r,id_cca):
                    return False
            return True
            
        def extraireActivites(self,constellation,r,m):
            act = []
            for (s,o,d) in constellation.getRequete(r).getMode(m).getCouples():
                if o not in act:
                    act.append(o)
                if d not in act:
                    act.append(d)
            return act        
    
        # inutile pour l'instant
        def tenterCorriger(self,constellation,mode):
            return False
                                   
        def ajouterExplication(self,explication):
            for ex in explication:
                if explication not in self.explications:
                    self.explications.append(ex)
                
        def explicationsModes(self,CCMs):
            res = []
            for (activites_solution,activites_mode) in self.explications:
                if self.contient(activites_solution,CCMs):
                    res.append(activites_mode)
            return res
         
        def nettoyerSolution(self,constellation):
            for mode in self.etat_recherche.getWorkingModes():
                self.annulerModeSolution(constellation,mode)
        
        def verifCCA(self,constellation):
            for s in self.solution.getSolCCAs():
                for s,cca in  self.solution.getSolCCAs()[s]:
                    for a in self.solution.getSolCCA(s,cca).getSequence():
                        assert(self.grapheDependances.grapheActivites.getNoeud(a).getComposanteConnexe()==(s,cca))
                        
            for a in self.grapheDependances.grapheActivites.getNoeuds():
                s,cca = self.grapheDependances.grapheActivites.getNoeud(a).getComposanteConnexe()
                for c in self.solution.getSolCCAs()[s]:
                    if cca!=c:
                        assert(a not in self.solution.getSolCCA(s,c).getSequence())
    
        def initIteration(self,constellation,noise):
            self.resetNoise(constellation,noise)
            self.etat_recherche = EtatRecherche()
            # ajouter modes actifs
            self.explication = {}
            self.redemarrer(constellation)
            self.initRequetes(constellation,noise)
            update_modes = []
            for r in constellation.getRequetes():
                m = constellation.getRequete(r).getModeCourant(constellation).getId()
                update_modes.append((r,m))
            for (r,m) in update_modes:
                self.ajouterModeActif(constellation,r,m) # etat recherche + CCA sol
            # stoquage etat des requetes dans les CCA
            self.requetesPresentesCCA(constellation) # cca => requetes,req => ccas 
            

        def resoudre(self,constellation,mailbox):
            temps = time()
            comm = MPI.COMM_WORLD
            rank = comm.Get_rank()
            size = comm.Get_size()
            temps = time()
            if size==1:
                freeCPUs = [0] # un seul coeur => slave local
            else:
                freeCPUs = list(range(1,size))
            last_record = 0
            super_iter = 0
            while time()-self.start_date<config.getOptValue("time"):
                sigma = config.getOptValue("noise")*super_iter
                self.initIteration(constellation,sigma)
                while time()-self.start_date<config.getOptValue("time") and self.etat_recherche.modesEnCours():
                    start_it = time()
                    #printColor(self.requetesCCA,self.ccaRequetes,c='y')
                    self.envoyerTache(constellation,freeCPUs,mailbox,super_iter)
                    self.analyserResultats(constellation,freeCPUs,mailbox,super_iter)

                    step()
                    
                    if config.getOptValue("verif"):
                        self.verifierCCA(constellation)
                    """
                        stockage de la solution courante
                    """
                    if last_record<( (time()-self.start_date)//10):
                        last_record = (time()-self.start_date)//10
                        self.notifierSansEvenement(constellation)
                        self.afficherInfo(time(),self.start_date,constellation,title="INFORMATION")
                self.notifierFinIteration(constellation)
                self.afficherInfo(time(),self.start_date,constellation,title="FIN ITERATION "+str(super_iter),color='r')
                super_iter += 1
                temps = time()
                if config.getOptValue("noise") == 0:# or self.gapRecompense(constellation) and self.gapTemps(constellation):
                    break
            
        
            # attendre ce qu'il reste
            #for x in range(len([i for i in range(1,comm.Get_size()) if i not in freeCPUs])):
                #data = comm.recv()
                #fself.analyserMessage(data,constellation,freeCPUs,mailbox)
            
            
            self.nettoyerSolution(constellation)
            self.verifierSolution(constellation)
            self.terminerProcessus()
            self.notifierFinExecution(constellation)
            self.afficherInfo(time(),self.start_date,constellation,title='FIN',color='y')
        
        def terminerProcessus(self):
            comm = MPI.COMM_WORLD
            rank = comm.Get_rank()
            size = comm.Get_size()
            for i in range(size):
                if i!=0:
                    comm.send({'fin':True},dest=i)
        
        def requetesPresentesCCA(self,constellation):
            self.requetesCCA,self.ccaRequetes = {},{}
            req,ccas = {},{}
            for r in constellation.getRequetes():
                req[r] = []
                m = self.modes_courants[r]
                rec = self.recompenseBruitee(constellation,r,m)
                p = constellation.getRequete(r).getPriorite()
                for (s,o) in constellation.getRequete(r).getCouples():
                    id_cca = self.grapheDependances.getActiviteCCA(o)
                    if id_cca not in ccas:
                        ccas[id_cca] = []
                    if id_cca not in req[r]:
                        req[r].append(id_cca)
                        reverse_insort(ccas[id_cca],(p,rec,r))
            self.requetesCCA = ccas
            self.ccaRequetes = req
        
        def ccaLibres(self):
            requetes = {} # cca ou il y a du travail
            ccas_libres = []
            for id_cca in self.requetesCCA:
                rec_max = max([x[1][0] for x in self.requetesCCA[id_cca]])
                r_max = [x[2] for x in self.requetesCCA[id_cca] if x[1][0]==rec_max or x[2]==-1]
                for r in r_max:
                    m = self.modes_courants[r]
                    if r not in requetes:
                        requetes[r] = self.etat_recherche.getCCALibres(r,m,self.grapheDependances)
                    if not self.etat_recherche.ccaOccupee(id_cca):
                        if id_cca not in ccas_libres:
                            ccas_libres.append(id_cca)
            return requetes,ccas_libres
        
        # renvoie les (r,cca) tels que r est premiere sur sa CCA et (r,cca) n'est pas encore calculée
        def candidatsTop(self,requetes,ccas_libres):
            candidats_top = []
            for id_cca in ccas_libres:
                prio = [x[2] for x in self.requetesCCA[id_cca] if self.requetePrioritaire(x[2])]
                assert(len(prio)<=1)
                #r = max(self.requetesCCA[id_cca],key=lambda r : (r[2]==-1,r[1][0],-r[2]))[2]
                #printColor(self.requetesCCA[id_cca],r,c='m')
                #if id_cca in requetes[r]:
                if len(prio)==1:
                    r = prio[0]
                    if id_cca in requetes[r]:
                        m = self.modes_courants[r]
                        candidats_top.append((0,(r,m),id_cca))
            return candidats_top
        
        # les candidats de niveaux inferieurs
        def candidatsQueue(self,requetes,ccas_libres):
            candidats_queue = []
            for id_cca in ccas_libres:
                niveau = 1
                while niveau<len(self.requetesCCA[id_cca]):
                    r = self.requetesCCA[id_cca][niveau][2]
                    if r not in requetes:
                        m = self.modes_courants[r]
                        requetes[r] = self.etat_recherche.getCCALibres(r,m,self.grapheDependances)
                    if id_cca in requetes.get(r,[]):
                        m = self.modes_courants[r]
                        candidats_queue.append((niveau,(r,m),id_cca))
                        break
                    niveau += 1
            return candidats_queue
        
        def choisirCandidatQueue(self,candidats_queue,constellation):
            if len(candidats_queue)>0:
                prio = lambda rm:rm[1][0]==-1
                rec = lambda rm : self.recompenseBruitee(constellation,rm[1][0],rm[1][1])[0]
                res = max(candidats_queue,key=lambda r:(-r[0],prio(r),rec(r),-r[1][0])) # min niveau,max prio,max rec
                return res
            else:
                return None
            
        def choisirCandidatTop(self,candidats_top,constellation):
            if len(candidats_top)>0: # candidats top = liste de (niveau,(r,m),cca)
                prio = lambda rm : rm[1][0]==-1
                rec = lambda rm : self.recompenseBruitee(constellation,rm[1][0],rm[1][1])[0]
                res = max(candidats_top,key=lambda r :(prio(r),rec(r),-r[1][0]))
                return res
            else:
                return None
            
        def choisirMode(self,constellation):
            #printColor([(cca,self.requetesCCA[cca][:3]) for cca in self.requetesCCA],c='y')
        
            start = time()
            # deteminer les ccas libres et les cca libres sur lesquelles les requetes ont du travail
            requetes,cca_libres = self.ccaLibres()
            #printColor(requetes,cca_libres,c='b')
            start = time()
            candidats_top = self.candidatsTop(requetes,cca_libres)
            #printColor(candidats_top,c='r')
            #printColor(self.etat_recherche.cca_occupees,c='m')
            res = self.choisirCandidatTop(candidats_top,constellation)
            if res is not None:
                return res
            # si pas de candidat et mode "anticiper" activé : on cherche les candidats dans la queue
            if not config.getOptValue("anticipation"):
                #requete_prio = [r for r in self.modes_courants if self.requetePrioritaire(r)][0]
                #r = requete_prio
                #printColor(requete_prio,self.etat_recherche.activities_to_do[(r,self.modes_courants[r])],c='m')
                #printColor(requete_prio,self.etat_recherche.activities_working[(r,self.modes_courants[r])],c='m')
                #printColor(requete_prio,self.etat_recherche.activities_done[(r,self.modes_courants[r])],c='m')
                return None
            else:
                candidats_queue = self.candidatsQueue(requetes,cca_libres)
                return self.choisirCandidatQueue(candidats_queue,constellation)
        
        def getModeActivites(self,constellation,r,m):
            activites = []
            for s in self.constellation.getRequete(r).getMode(m).getActivites():
                for o in self.constellation.getRequete(r).getMode(m).getActivites()[s]:
                    if o not in activites:
                        activites.append(o)
            return activites
        
        def processusDispo(self,freeCPUs):
            return freeCPUs!=[]
        
        def choisirCPU(self,freeCPUs):
            return freeCPUs.pop()
        
        def envoyerTache(self,constellation,freeCPUs,mailbox,iteration):
            while self.processusDispo(freeCPUs) and self.etat_recherche.modesEnCours():
                res = self.choisirMode(constellation)
                if res is not None:
                    niveau,mode,id_cca = res
                    cpu = self.choisirCPU(freeCPUs)
                    activites = self.etat_recherche.getCCA(mode[0],mode[1],self.grapheDependances,id_cca)
                    id_cca = self.grapheDependances.getActiviteCCA(activites[0])
                    cca = self.getSolCCA(id_cca[0],id_cca[1])
                    niveau_msg = config.getOptValue("anticipation")*("rang : "+str(niveau))
                    printColor("mode choisi :",str(mode),"sur la CCA",str(id_cca)+". Activités : "+str(activites),niveau_msg,c='b')
                
                    for a in activites:
                        if a in cca.sequence:
                            printColor(a,cca,c='r')
                            assert(False)
                    # demander au slave d'inserer les activites dans la cca
                    if cpu==0:
                        data = {'iteration':iteration,'mode':mode,'cca':deepcopy(cca),'source':0,'activites':activites,'time':time()}
                        self.localSlave.insererActivites(constellation,data,self.modeleDeTransition)
                    else:
                        mailbox.demanderTache(mode,activites,cca,cpu,iteration)
                    #self.solution.historique.enregistrerEnvoi(cpu)
                else:
                    #printColor(len(freeCPUs),c='m')
                    break
        
        def annulerModeSolution(self,constellation,mode):
            r,m = mode
            cca = []
            for (s,o) in constellation.getRequete(r).getMode(m).getCouples():
                id_cca = self.grapheDependances.getActiviteCCA(o)
                if id_cca not in cca:
                    cca.append(id_cca)
            for (s,cca_a) in cca:
                self.getSolCCA(s,cca_a).annulerMode(constellation,r,m,self.modeleDeTransition)
            
        # si retrait est None on retire toutes les activites du mode. Sinon on retire la liste 'retrait'
        def retirerModeSolution(self,constellation,mode,retrait=None):
            r,m = mode
            cca = {}                     
            for (s,o) in constellation.getRequete(r).getMode(m).getCouples():
                if retrait is None or o in retrait:
                    id_cca = self.grapheDependances.getActiviteCCA(o)
                    if id_cca not in cca:
                        cca[id_cca] = []
                    if o in self.getSolCCA(id_cca[0],id_cca[1]).getSequence():
                        cca[id_cca].append(o)
            for (s,cca_a) in cca:
                self.getSolCCA(s,cca_a).retirerListeActivites(constellation,cca[(s,cca_a)],self.modeleDeTransition)
                                   
        # active le mode pending suivant s'il y en a un et rajoute un pending mode s'il y en a encore
        def defilerEtatRequete(self,constellation,mode,explication,grapheDep):
            r,m = mode
            res = constellation.getRequete(r).getModeSuivant(explication,constellation)
            modes_a_redemarrer = []
            if res is not None:
                m = res.getId()
                self.modes_courants[r] = m
                # retrait : les activities presentes uniquement dans l'ancien mode. Les retirer
                # prevalider : presentes dans les deux modes et deja validées dans l'ancien => Prevalider dans le nouveau mode
                retrait,prevalider,intersection_working = self.intersectionActivites(constellation,r)
                #printColor(retrait,prevalider,intersection_working,c='m')
                if self.etat_recherche.getEtat(r,m-1)!=self.etat_recherche.FAILED:
                    modes_a_redemarrer = self.etat_recherche.echecCertain(constellation,r,m-1,grapheDep,retrait)
                else:
                    modes_a_redemarrer  = []
                self.retirerModeSolution(constellation,mode,retrait) # retirer de la solution
            else:
                modes_a_redemarrer = self.etat_recherche.echecCertain(constellation,r,m,grapheDep)
                self.retirerModeSolution(constellation,mode) # retirer de la solution
            if res is not None:
                self.ajouterModeActif(constellation,r,m,prevalider) # prevalidation de l'intersection deja validee
                self.etat_recherche.preWorking(r,m,intersection_working)
            else:
                printColor('requete',r,'('+str(len(constellation.getRequete(r).getCouples())) +' activités) n\'a plus de modes',c='m')
                self.etat_recherche.notifierEchecRequete(r)
                self.supprimerRequete(constellation,r)
            return modes_a_redemarrer
                
        def retraitRequeteCCA(self,id_cca,r):
            if id_cca in self.requetesCCA:
                for i,x in enumerate(self.requetesCCA[id_cca]):
                    if x[2]==r:
                        self.requetesCCA[id_cca].pop(i)
                        if self.requetesCCA[id_cca]==[]:
                            del self.requetesCCA[id_cca]
    
    class Processus:
        def __init__(self,role):
            self.role = role
        
        def resoudre(self,constellation):
            pass
        
    class Master(Processus):
        def __init__(self):
            super().__init__("master")

        def resoudre(self,constellation,start_date,mailbox,modeleDeTransition,dt_construction_transition,tlim,CCAs,solution):
            # iterer insertion-réparation
            self.initSolution(constellation,start_date,modeleDeTransition,dt_construction_transition,tlim,CCAs,solution)
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
            
        def initSolution(self,constellation,start_date,modeleDeTransition,dt_construction_transition,tlim,CCAs,solution):
            self.solution = UnitParallelCCASearch(constellation,start_date,modeleDeTransition,dt_construction_transition,tlim,CCAs,solution)
            
        def getSolution(self):
            return self.solution.getSolution()
        
        def setSolution(self,sol,vid,modes_retenus,modes_candidats,modes_retires,objectif):
            self.solution.setSol(sol,vid,modes_retenus,modes_candidats,modes_retires,objectif)
        
        def verifierSolution(self):
            self.solution.verifierSolution()
            
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
    
            def getSatellite(self):
                return self.solution.getSatellite()
    
            def getIdentifiant(self):
                return self.solution.getIdentifiant()
        
            def getSolution(self):
                return self.solution
            
            def getSequence(self):
                return self.solution.getSequence()
            
            # pas vraiment une explication. Plutot la liste d'ajouts ?
            # /!\ en cas d'échec : garde les activités (pour ça qu'on récupère l'ancienne solution)
            # c'est le master qui gère le retrait des activités ratées
            def planifierMode(self,constellation,mode,activites):
                r,m = mode
                explication = self.solution.ajouterActivites(constellation,activites,self.modeleDeTransition)
                derniereSolution = deepcopy(self.solution)
                succes = self.solution.calculerSequence(constellation,self.modeleDeTransition)
                if not succes:
                    self.solution = derniereSolution
                return explication
                    
            def sequenceFaisable(self,constellation):
                return self.solution.sequenceFaisable(constellation,self.modeleDeTransition)
            
        def __init__(self,constellation,modeleDeTransition):
            comm = MPI.COMM_WORLD
            rank = comm.Get_rank()
            super().__init__("slave "+str(rank))
            activites = constellation.extraireActivitesRequetes()
            self.grapheDependances = GroupeComposantesActivitesStatiques(constellation,activites,modeleDeTransition)
            
        def estMessageFin(self,data):
            return 'fin' in data
        
        def etendreDureesCalcul(self):
            durees = []
            for i,date in enumerate(self.dates_cca):
                if i<len(self.dates_cca)-1:
                    durees.append((date,self.dates_cca[i+1]-date))
            return durees
        
        def packSolution(self,cca,mode,faisable,explication,iteration):
            return {'iteration':iteration,'cpu_time':self.etendreDureesCalcul(),'reception':time(),'envoi':(self.date_envoi,self.duree_envoi),'mode':mode,'cca':deepcopy(cca),'faisable':faisable,'source':rank,'explication':explication}
    
        def insererActivites(self,constellation,data,modeleDeTransition):
            self.duree_envoi = time() - data['time']
            self.start_calcul = time()
            self.date_envoi = data['time']
            iteration = data['iteration']
            start = time()
            mode = data['mode'] # [(r,m)]
            activites = data['activites']
            sol_cca = data['cca']
            # ajouter la sequence comme sequence initiale
            self.planif = self.PlanificateurCCA(sol_cca,modeleDeTransition)
            # solution initiale = sequence (activites des modes deja planifies)
            explication = self.planif.planifierMode(constellation,mode,activites)
            solution_cca = self.planif.getSolution()
            faisable = self.planif.sequenceFaisable(constellation)
            self.dates_cca = [self.start_calcul,time()]
            self.resultat = self.packSolution(solution_cca,mode,faisable,explication,iteration)
                    
        def resoudre(self,constellation,mailbox,modeleDeTransition):
            comm = MPI.COMM_WORLD
            while True:
                data = comm.recv()
                if self.estMessageFin(data):
                    break
                else:
                    self.insererActivites(constellation,data,modeleDeTransition)
                    mailbox.posterMessage(self.resultat)
        
        def getResultat(self):
            if MPI.COMM_WORLD.Get_size()==1:
                self.resultat['reception'] = (self.resultat['reception'],time() - self.resultat['reception'])
            return self.resultat
        
        
    path = '../data'
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    
    def envoyerMessage(mode,activites,cca,cpu,iteration):
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        assert(rank==0)
        data = {'iteration':iteration,'mode':mode,'cca':cca,'source':rank,'activites':activites,'time':time()}
        comm.send(data, dest=cpu, tag=rank)
    
    class MessagerieMessageUnique:
        def __init__(self):
            pass
        
        def posterMessage(self,data):
            assert(MPI.COMM_WORLD.Get_rank()>0)
            MPI.COMM_WORLD.send(data,dest=0)
            
        def readMessages(self):
            assert(MPI.COMM_WORLD.Get_rank()==0)
            data = MPI.COMM_WORLD.recv()
            print("test")
            data['reception'] = (data['reception'],time()-data['reception'])
            return [data]
        
        def demanderTache(self,mode,activites,cca,cpu,iteration):
            assert(MPI.COMM_WORLD.Get_rank()==0)
            envoyerMessage(mode,activites,cca,cpu,iteration)
            
        def __str__(self):
            return "Messagerie à message unique"
        
    class MessagerieSynchronisee:
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
        def demanderTache(self,mode,activites,cca,cpu,iteration):
            assert(MPI.COMM_WORLD.Get_rank()==0)
            self.cpu_mapped[cpu-1] = 1
            envoyerMessage(mode,activites,cca,cpu,iteration)
            
        # appelé par les slaves
        def posterMessage(self,data):
            assert(MPI.COMM_WORLD.Get_rank()>0)
            MPI.COMM_WORLD.send(data,dest=0)
            rank = MPI.COMM_WORLD.Get_rank()
            #self.cpu_mapped[rank-1] = 1
            
        def readMessages(self):
            assert(MPI.COMM_WORLD.Get_rank()==0)
            nCPUMapped = sum([self.cpu_mapped[i] for i in range(self.size)])
            for i in range(nCPUMapped):
                data = MPI.COMM_WORLD.recv()
                data['reception'] = (data['reception'],time()-data['reception'])
                source = data["source"]
                self.cpu_mapped[source-1] = 0
                yield data
            else:
                return None
        
        def __str__(self):
            return "Messagerie synchronisée"
        
    class MessageriePartagee:
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
            self.flag_messages, itemsize = self.win.Shared_query(0)
            for i in range(self.size):
                self.flag_messages[i] = False
            #self.flag_messages = np.ndarray(buffer=self.buf, dtype=bool, shape=(size,))
        
        def demanderTache(self,mode,activites,cca,cpu,iteration):
            assert(MPI.COMM_WORLD.Get_rank()==0)
            envoyerMessage(mode,activites,cca,cpu,iteration)
    
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
                    assert(process==data['source'])
                    self.flag_messages[i] = False
                    yield data
                else:
                    yield None
        
        def __str__(self):
            return "Messagerie partagée : "+ str([self.flag_messages[i] for i in range(self.size)])
    
    # Création de l'objet GLOBAL mailbox
    def creerMailbox():
        comm = MPI.COMM_WORLD 
        if config.glob.sync:
            mailbox = MessagerieSynchronisee()
        else:
            mailbox = MessageriePartagee()
            #mailbox = MessagerieMessageUnique()
        return mailbox   

    
    class runnableUPCCAS:
        def execute(self,constellation,start_date,modeleDeTransition,dt_construction_transition,tlim=np.inf,CCAs=None,solution=None):
            mailbox = MessageriePartagee()
            if rank == 0:
                self.process = Master()
                self.process.resoudre(constellation,start_date,mailbox,modeleDeTransition,dt_construction_transition,tlim,CCAs,solution)
                
                for i in range(1,MPI.COMM_WORLD.Get_size()):
                    MPI.COMM_WORLD.send({"sol":self.process.solution.getSolutionContainer()},dest=i)
                comm.Barrier()
                return self.process.solution.getSolutionContainer()
            else:
                self.process = Slave(constellation,modeleDeTransition)
                self.process.resoudre(constellation,mailbox,modeleDeTransition)
                
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
        
from solution import *
from constellation import *
from composantes import *
from solution_composantes import *
from mpi4py import MPI
from Utils import *
from config import *
from time import time
from time import sleep
from graphes import *

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
            
    def libererCCA(self,cca):
        self.cca_occupees.remove(cca)
        
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
        #assert(config.glob.anticipation_modes)
        self.a_supprimer.remove((s,o))
        
    # suppression du mode pour la version avec anticipation
    def supprimerAnticipation(self,r,m,grapheDep,liste_a_supprimer=None):
        self.terminerMode(r,m,False)
        for (s,o) in constellation.getRequete(r).getMode(m).getCouples():
            if liste_a_supprimer is None or o in liste_a_supprimer:
                if((s,o) not in self.a_supprimer):
                    cca = grapheDep.getActiviteCCA(o)
                    if cca in self.cca_occupees:
                        self.a_supprimer.append((s,o))
        
    # annule un mode et renvoie la liste de ceux a redemarrer. SUPRESSION DU MODE
    def echecCertain(self,r,m,grapheDep,liste_a_supprimer=None):
        self.supprimerAnticipation(r,m,grapheDep,liste_a_supprimer)
        if not config.glob.anticipation_modes:
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
        assert(config.glob.anticipation_modes)
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
    
    def ccaToDo(self,r,m,cca,grapheDep):
        for a in self.activities_to_do[(r,m)]:
            to_do_cca = grapheDep.getActiviteCCA(a)
            if to_do_cca == cca:
                return True
        return False
    """
            features
    """
    def filtrerModesEnCours(self,parents):
        return [x for x in parents if x in self.activities_to_do and self.activities_to_do[x]!=[] and self.etat_mode[x] == self.WORKING]
    
    def modesEnCours(self):
        return len(self.working_modes)>0
     
    def ccaOccupee(self,cca):
        return cca in self.cca_occupees
    
    def getActivitesRestantes(self,r,m):
        return self.activities_to_do[(r,m)]
    
    def getCCALibres(self,r,m,grapheDep):
        ccas = [grapheDep.getActiviteCCA(a) for a in self.activities_to_do[(r,m)]]
        ccas_libre = [cca for cca in ccas if cca not in self.cca_occupees]
        return ccas_libre
    
    def getCCA(self,r,m,grapheDep,cca):
        if (r,m) not in self.started_modes:
            self.started_modes.append((r,m))
        assert(cca not in self.cca_occupees)
        act = [a for a in self.activities_to_do[(r,m)] if grapheDep.getActiviteCCA(a)==cca]
        supprimerListeElements(self.activities_to_do[(r,m)],act)
        self.cca_occupees.append(cca)
        for a in act:
            self.activities_working[(r,m)].append(a)
        return act
        
    # depile les activites et les renvoie
    def getCCASuivante(self,r,m,grapheDep):
        assert(len(self.activities_to_do[(r,m)])>0)
        if (r,m) not in self.started_modes:
            self.started_modes.append((r,m))
        ccas = [grapheDep.getActiviteCCA(a) for a in self.activities_to_do[(r,m)]]
        ccas_libre = [cca for cca in ccas if cca not in self.cca_occupees]
        cca = ccas_libre[0]
        self.cca_occupees.append(cca)
        liste_a = []
        act = [a for a in self.activities_to_do[(r,m)] if grapheDep.getActiviteCCA(a)==cca]
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
            for cca in liste_cca:
                self.cca_occupees.append(cca)
    
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
        if(config.glob.anticipation_modes):
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
        printColor('Prévalidation des activités',activites,(r,m),c='m')
        supprimerListeElements(self.activities_to_do[(r,m)],activites)
        if (r,m) not in self.activities_done:
            self.activities_done[(r,m)] = []
        self.activities_done[(r,m)] += activites
        if len(self.activities_to_do[(r,m)])==0 and len(self.activities_working[(r,m)])==0:
            return True
        else:
            return False
        
    def validerActivites(self,activites,r,m):
        printColor('Validation des activités',activites,(r,m),c='c')
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



class UnitParallelCCASearch(Solution):
    def __init__(self,constellation):
        super().__init__(constellation)
        # recuperer toutes les observations
        activites = {}
        for r in constellation.getRequetes():
            act = constellation.getRequete(r).getActivites()
            for s in act:
                if s not in activites:
                    activites[s] = []
                for a in act[s]:
                    if a not in activites[s]:
                        activites[s].append(a)
        self.grapheDependances = GroupeComposantesActivitesStatiques(constellation,activites)
        printColor(self.grapheDependances,c='b')
        # si 1 coeur : slave local
        comm = MPI.COMM_WORLD
        if comm.Get_size()==1:
            self.localSlave = Slave()
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
        
    def presentCCA(self,constellation,r,m,cca,grapheActivites):
        for s in constellation.getRequete(r).getMode(m).getActivites():
            for o in constellation.getRequete(r).getMode(m).getActivitesSatellite(s):
                if grapheActivites.getNoeud(o).getComposanteConnexe()==cca:
                    return True
        return False
    
    def getComposanteActivite(self,a):
        return self.grapheDependances.getActiviteCCA(a)
    
    def setCCA(self,constellation,s,seq,cca,r,m):
        self.solCCAs[s][cca].setSequence(constellation,seq)

    def validerActivites(self,constellation,activites,mode):
        s = constellation.getSatelliteActivite(activites[0])
        act = {a : None for a in activites}
        cca = self.getComposanteActivite(activites[0])
        seq = []
        for i,a in enumerate(activites):
            cca_a = self.getComposanteActivite(a)
            if cca_a!=cca:
                self.setCCA(constellation,s,seq,cca,mode[0],mode[1])
                cca = cca_a
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
        for cca in nouvelles_ccas:
            req = []
            for (p,w,r) in requetes:
                for s in constellation.getRequete(r).getActivites():
                    for o in constellation.getRequete(r).getActivites()[s]:
                        if o in self.grapheDependances.getComposanteConnexe(cca):
                            reverse_insort(req,(prio(r),rec(r),r))
                            cca_req[r].append(cca)
                            break
            self.requetesCCA[cca] = req
        # MAJ CCA Requetes
        for r in cca_req:
            self.ccaRequetes[r].remove(ancienne_cca)
            self.ccaRequetes[r] += cca_req[r]
        
        for cca in nouvelles_ccas:
            if self.requetesCCA[cca]==[]:
                del self.requetesCCA[cca]
            
        
    """
        =============================================== 
                        RESOLUTION
        =============================================== 
    """
   
    def validerMode(self,constellation,mode):
        printOpen("Validation du mode "+str(mode),c='g')
        start_val = time()
        # retirer la requete de toutes ses cca
        for cca in self.ccaRequetes[mode[0]]:
            self.retraitRequeteCCA(cca,mode[0])
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
        assert(mode not in self.modes_retenus)
        self.modes_retenus.append(mode)
        self.historique.enregistrerRequeteValidee(constellation,mode[0])
        requetes_a_annuler = self.etat_recherche.validerMode(mode[0],mode[1])
        
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
        return mode[1] == self.modes_courants[r_courant]
        
    def analyserMessage(self,data,constellation,freeCPUs,mailbox,iteration):
        start = time()
        self.historique.enregistrerReceptionMessage(data)
        i = data['source']
        assert(i not in freeCPUs)
        bisect.insort(freeCPUs,i)
        it = data['iteration']
        if iteration==it:
                cpu_time,cca,succes,mode = data['cpu_time'],data['cca'],data['faisable'],data['mode']
                explication = data['explication']
                seq = cca.getSequence()
                # si message d'un mode passé : transfert potentiel d'activités validées ou à redemarrer
                if config.UPCCAS.prevalidation or self.modeCourant(mode):
                    retraits = self.supprimerAnciennesActivites(seq)
                    if succes:
                        mode = self.transfertSuccesActivites(explication,seq,mode)
                    else:
                        mode,activites_a_redemarrer = self.transfertEchecActivites(explication,retraits,mode)
                    # indiquer la cca libre et supprimer les activités n'étant plus à programmer
                self.libererCCA(seq,explication)
                if config.UPCCAS.prevalidation or self.modeCourant(mode):
                    # traiter le message d'un mode en cours
                    if succes and not self.etat_recherche.requeteEchouee(mode[0]) and self.etat_recherche.getEtat(mode[0],mode[1])!=self.etat_recherche.FAILED:
                        stat_validation = time()
                        if(len(seq)==0):
                            print(mode)
                        assert(len(seq)>0)
                        self.validerActivites(constellation,seq,mode)
                        start_validation = time()
                        termine = self.etat_recherche.validerActivites(explication,mode[0],mode[1])
                        self.historique.enregistrerSucces(i,cpu_time)
                        if termine:
                            if config.UPCCAS.anticipation_modes:
                                self.validerAnticipation(constellation,mode)
                            else:
                                self.validerMode(constellation,mode)
                            self.objectif = self.calculerObjectif(constellation)
                            #self.historique.MAJHistorique(time()-self.start_date,0,self.objectif,None,self.modes_retenus,self.grapheDependances,constellation)
                    else:
                        if self.etat_recherche.getEtat(mode[0],mode[1])!=self.etat_recherche.FAILED:
                            if config.UPCCAS.anticipation_modes:
                                self.notifierEchecAnticipation(constellation,mode,explication,activites_a_redemarrer)
                            else:
                                self.planifierModeSuivant(constellation,mode,explication)
                            #self.historique.MAJHistorique(time()-self.start_date,0,self.objectif,None,self.modes_retenus,self.grapheDependances,constellation)
                        self.historique.enregistrerEchec(i,cpu_time)
                    
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
    
    def transfertEchecActivites(self,explication,retrait,mode):
        if mode[0] in self.modes_courants:
            nouveau_mode = self.modes_courants[mode[0]]
            if nouveau_mode!=mode[1]:
                a_redemarrer = explication.copy()
                for a in retrait:
                    if a in a_redemarrer:
                        a_redemarrer.remove(a)
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
            cca = self.grapheDependances.getActiviteCCA(a)
            if cca not in ccas:
                ccas.append(cca)
        # detection de fautifs : en attente s'il y en a
        premierPartout = True
        for cca in ccas:
            if not self.premierCCA(mode[0],cca):
                premierPartout = False
                self.explication[mode] = self.explication.get(mode,[])
                for a in explication:
                    if a not in self.explication[mode]:
                        self.explication[mode].append(a)
                fautifs = self.findFautifsCCA(mode[0],cca)
                self.etat_recherche.echecIncertain(mode[0],mode[1],fautifs,explication,True)
                break
                
        # III) si pas de fautif : annuler car le mode ne passera jamais
        if premierPartout:
            self.planifierModeSuivant(constellation,mode,explication)
            self.explication[mode] = []
    
    def findFautifsCCA(self,r,cca):
        i = 0
        while self.requetesCCA[cca][i][2]!=r:
            i += 1
        return [x[2] for x in self.requetesCCA[cca][:i]]
    
    def premierCCA(self,r,cca):
        return r==self.requetesCCA[cca][0][2]
    
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
            for cca in self.requetesCCA[r]:
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
                #printColor("mode "+str(mode_attente),c='y')
                if self.requetePrioritaire(mode_attente[0]):
                    validation = True
                    self.validerMode(constellation,mode_attente)
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
            if o in act_ancien and o in self.etat_recherche.activities_done[(r,ancien_mode)]:
                if config.UPCCAS.prevalidation:
                    prevalider.append(o)
                else:
                    retrait.append(o)
            elif o in self.etat_recherche.activities_working[(r,ancien_mode)]:
                if config.UPCCAS.prevalidation:
                    intersection_working.append(o)
        for o in act_ancien:
            if o not in act_nouveau_mode:
                retrait.append(o)
        
        return retrait,prevalider,intersection_working
    
    def planifierModeSuivant(self,constellation,mode,explication):
        exp_cca = [(o,self.grapheDependances.getActiviteCCA(o)) for o in explication]
        printOpen('Echec du mode',mode,'cas 3 : annuler. cause :',exp_cca,c='r')
        #printColor("Echec du mode "+str(mode),explication,c='r')
        self.historique.enregistrerEchecMode(constellation,mode[0])
        for cca in self.ccaRequetes[mode[0]]:
            self.retraitRequeteCCA(cca,mode[0])
        modes_a_redemarrer = self.defilerEtatRequete(constellation,mode,explication,self.grapheDependances) # MAJ etat recherche MAJ graphe (pending)
        self.insererRequeteCCA(constellation,mode[0])
        if len(modes_a_redemarrer)>0:
            printColor('Redémarrer les requêtes',modes_a_redemarrer,c='m')
        if config.glob.anticipation_modes:
            for mode in modes_a_redemarrer:
                self.etat_recherche.redemarrerMode(mode,self.modes_courants[mode])
        printClose()

    def insererRequeteCCA(self,constellation,r):
        ccas = []
        if r in self.modes_courants:
            for (s,o) in constellation.getRequete(r).getCouples():
                try:
                    cca = self.grapheDependances.getActiviteCCA(o)
                    if cca not in ccas:
                        ccas.append(cca)
                except:
                    pass
            self.ccaRequetes[r] = ccas
            prio = constellation.getRequete(r).getPriorite()
            rec = self.recompenseBruitee(constellation,r,self.modes_courants[r])
            for cca in ccas:
                if cca not in self.requetesCCA:
                    self.requetesCCA[cca] = []
                reverse_insort(self.requetesCCA[cca],(prio,rec,r))
    
    def requetePrioritaireCCA(self,r,cca):
        return r==self.requetesCCA[cca][0][2]
    
    def requetePrioritaire(self,r):
        for cca in self.ccaRequetes[r]:
            if not self.requetePrioritaireCCA(r,cca):
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
        for s in self.solCCAs:
            for cca in self.solCCAs[s]:
                for a in self.solCCAs[s][cca].sequence:
                    assert(self.grapheDependances.grapheActivites.getNoeud(a).getComposanteConnexe()==cca)
                    
        for a in self.grapheDependances.grapheActivites.getNoeuds():
            cca = self.grapheDependances.grapheActivites.getNoeud(a).getComposanteConnexe()
            s = constellation.mapping_sat[a]
            for c in self.solCCAs[s]:
                if cca!=c:
                    assert(a not in self.solCCAs[s][c].sequence)

    def initIteration(self,constellation,noise):
        self.resetNoise(constellation,noise)
        self.etat_recherche = EtatRecherche()
        # ajouter modes actifs
        self.explication = {}
        self.solCCAs = {s : {} for s in constellation.satellites} # liste de solutions CCA
        for cca in self.grapheDependances.getComposantes():
            s = self.grapheDependances.getSatelliteCCA(cca)
            self.solCCAs[s][cca] = SolCCA(cca,s)
        self.initRequetes(constellation,noise)
        update_modes = []
        for r in constellation.getRequetes():
            m = constellation.getRequete(r).getModeCourant(constellation).getId()
            update_modes.append((r,m))
        for (r,m) in update_modes:
            self.ajouterModeActif(constellation,r,m) # etat recherche + CCA sol
        # stoquage etat des requetes dans les CCA
        self.requetesPresentesCCA(constellation) # cca => requetes,req => ccas 
        self.objectif = (0,0)
        self.modes_retenus = []
            
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
        while time()-self.start_date<config.glob.limite_temps:
            sigma = config.glob.noise*super_iter
            self.initIteration(constellation,sigma)
            while time()-self.start_date<config.glob.limite_temps and self.etat_recherche.modesEnCours():
                start_it = time()
                self.envoyerTache(constellation,freeCPUs,mailbox,super_iter)
                self.analyserResultats(constellation,freeCPUs,mailbox,super_iter)
                if config.glob.verif:
                    self.verifierCCA(constellation)
                """
                    stockage de la solution courante
                """
                self.objectif = self.calculerObjectif(constellation)
                if last_record<( (time()-self.start_date)//10):
                    last_record = (time()-self.start_date)//10
                    self.MAJHistorique(time()-self.start_date,0,self.solCCAs,self.modes_retenus,self.grapheDependances,constellation)
                    self.afficherInfo(time(),self.start_date,constellation,title="INFORMATION")
            self.MAJHistorique(time()-self.start_date,1,self.solCCAs,self.modes_retenus,self.grapheDependances,constellation)
            self.afficherInfo(time(),self.start_date,constellation,title="FIN ITERATION "+str(super_iter),color='r')
            super_iter += 1
            temps = time()
            if config.glob.noise == 0:# or self.gapRecompense(constellation) and self.gapTemps(constellation):
                break
        
        # attendre ce qu'il reste
        #for x in range(len([i for i in range(1,comm.Get_size()) if i not in freeCPUs])):
            #data = comm.recv()
            #fself.analyserMessage(data,constellation,freeCPUs,mailbox)
        
        
        self.nettoyerSolution(constellation)
        self.verifierSolution(constellation)
        self.terminerProcessus()
        self.MAJHistorique(min(time()-self.start_date,config.glob.limite_temps),2,self.solCCAs,self.modes_retenus,self.grapheDependances,constellation)
        #self.objectif,self.solCCAs,self.modes_retenus = self.getBest()
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
                cca = self.grapheDependances.getActiviteCCA(o)
                if cca not in ccas:
                    ccas[cca] = []
                if cca not in req[r]:
                    req[r].append(cca)
                    reverse_insort(ccas[cca],(p,rec,r))
        self.requetesCCA = ccas
        self.ccaRequetes = req
    
    def ccaLibres(self):
        requetes = {} # cca ou il y a du travail
        ccas_libres = []
        #if not config.glob
        for cca in self.requetesCCA:
            rec_max = max([x[1][0] for x in self.requetesCCA[cca]])
            r_max = [x[2] for x in self.requetesCCA[cca] if x[1][0]==rec_max or x[2]==-1]
            for r in r_max:
                m = self.modes_courants[r]
                if r not in requetes:
                    requetes[r] = self.etat_recherche.getCCALibres(r,m,self.grapheDependances)
                if not self.etat_recherche.ccaOccupee(cca):
                    if cca not in ccas_libres:
                        ccas_libres.append(cca)
        return requetes,ccas_libres
    
    # renvoie les (r,cca) tels que r est premiere sur sa CCA et (r,cca) n'est pas encore calculée
    def candidatsTop(self,requetes,ccas_libres):
        candidats_top = []
        for cca in ccas_libres:
            r = max(self.requetesCCA[cca],key=lambda r : (r[2]==-1,r[1][0],-r[2]))[2]
            if cca in requetes[r]:
                m = self.modes_courants[r]
                candidats_top.append((0,(r,m),cca))
        return candidats_top
    
    # les candidats de niveaux inferieurs
    def candidatsQueue(self,requetes,ccas_libres):
        candidats_queue = []
        for cca in ccas_libres:
            niveau = 1
            while niveau<len(self.requetesCCA[cca]):
                r = self.requetesCCA[cca][niveau][2]
                if r not in requetes:
                    m = self.modes_courants[r]
                    requetes[r] = self.etat_recherche.getCCALibres(r,m,self.grapheDependances)
                if cca in requetes.get(r,[]):
                    m = self.modes_courants[r]
                    candidats_queue.append((niveau,(r,m),cca))
                    break
                niveau += 1
        return candidats_queue
    
    def choisirCandidatQueue(self,candidats_queue,constellation):
        if len(candidats_queue)>0:
            prio = lambda rm:rm[1][0]==-1
            rec = lambda rm : self.recompenseBruitee(constellation,rm[1][0],rm[1][0])[0]
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
        """
        req = []
        for cca in self.requetesCCA:
            for niv,rec,r in self.requetesCCA[cca]:
                if r not in req:
                    req.append(r)
        etat = [(r,self.etat_recherche.getEtat(r,self.modes_courants[r])) for r in req]
        printColor(etat,c='m')
        """    
    
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
        if not config.glob.anticipation_modes:
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
                niveau,mode,cca = res
                cpu = self.choisirCPU(freeCPUs)
                activites = self.etat_recherche.getCCA(mode[0],mode[1],self.grapheDependances,cca)
                id_cca = self.grapheDependances.getActiviteCCA(activites[0])
                s = constellation.getSatelliteActivite(activites[0])
                cca = self.solCCAs[s][id_cca]
                niveau_msg = config.glob.anticipation_modes*("rang : "+str(niveau))
                printColor("mode choisi",str(mode),str(id_cca),str(activites),niveau_msg,c='b')
            
                for a in activites:
                    if a in cca.sequence:
                        printColor(a,cca,c='r')
                        assert(False)
                # demander au slave d'inserer les activites dans la cca
                if cpu==0:
                    data = {'iteration':iteration,'mode':mode,'cca':deepcopy(cca),'source':0,'activites':activites,'time':time()}
                    self.localSlave.insererActivites(constellation,data)
                else:
                    mailbox.demanderTache(mode,activites,cca,cpu,iteration)
                #self.historique.enregistrerEnvoi(cpu)
            else:
                #printColor(len(freeCPUs),c='m')
                break
    
    def annulerModeSolution(self,constellation,mode):
        r,m = mode
        cca = []
        for (s,o) in constellation.getRequete(r).getMode(m).getCouples():
            cca_o = self.grapheDependances.getActiviteCCA(o)
            if (s,cca_o) not in cca:
                cca.append((s,cca_o))
        for (s,cca_a) in cca:
            self.solCCAs[s][cca_a].annulerMode(constellation,r,m)
        
    # si retrait est None on retire toutes les activites du mode. Sinon on retire la liste 'retrait'
    def retirerModeSolution(self,constellation,mode,retrait=None):
        r,m = mode
        cca = {}                     
        for (s,o) in constellation.getRequete(r).getMode(m).getCouples():
            if retrait is None or o in retrait:
                cca_o = self.grapheDependances.getActiviteCCA(o)
                if (s,cca_o) not in cca:
                    cca[(s,cca_o)] = []
                if o in self.solCCAs[s][cca_o].getSequence():
                    cca[(s,cca_o)].append(o)
        for (s,cca_a) in cca:
            self.solCCAs[s][cca_a].retirerListeActivites(constellation,cca[(s,cca_a)])
                               
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
            printColor(retrait,prevalider,intersection_working,c='m')
            modes_a_redemarrer = self.etat_recherche.echecCertain(r,m-1,grapheDep,retrait)
            self.retirerModeSolution(constellation,mode,retrait) # retirer de la solution
        else:
            modes_a_redemarrer = self.etat_recherche.echecCertain(r,m,grapheDep)
            self.retirerModeSolution(constellation,mode) # retirer de la solution
        if res is not None:
            self.ajouterModeActif(constellation,r,m,prevalider) # prevalidation de l'intersection deja validee
            self.etat_recherche.preWorking(r,m,intersection_working)
        else:
            printColor('requete',r,'('+str(len(constellation.getRequete(r).getCouples())) +' activités) n\'a plus de modes',c='m')
            self.etat_recherche.notifierEchecRequete(r)
            self.supprimerRequete(constellation,r)
        return modes_a_redemarrer
            
    def retraitRequeteCCA(self,cca,r):
        if cca in self.requetesCCA:
            for i,x in enumerate(self.requetesCCA[cca]):
                if x[2]==r:
                    self.requetesCCA[cca].pop(i)
                    if self.requetesCCA[cca]==[]:
                        del self.requetesCCA[cca]

class Processus:
    def __init__(self,role):
        self.role = role
    
    def resoudre(self,constellation):
        pass
    
class Master(Processus):
    def __init__(self):
        super().__init__("master")
        
    def resoudre(self,constellation,mailbox):
        # iterer insertion-réparation
        self.initSolution(constellation)
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
        
    def initSolution(self,constellation):
        self.solution = UnitParallelCCASearch(constellation)
        
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
        def __init__(self,sol_cca):
            self.solution = sol_cca

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
            explication = self.solution.ajouterActivites(constellation,activites)
            derniereSolution = deepcopy(self.solution)
            succes = self.solution.calculerSequence(constellation)
            if not succes:
                self.solution = derniereSolution
            return explication
                
        def sequenceFaisable(self,constellation):
            return self.solution.sequenceFaisable(constellation)
        
    def __init__(self):
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        super().__init__("slave "+str(rank))
        activites = constellation.extraireActivitesRequetes()
        self.grapheDependances = GroupeComposantesActivitesStatiques(constellation,activites)
        
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

    def insererActivites(self,constellation,data):
        self.duree_envoi = time() - data['time']
        self.start_calcul = time()
        self.date_envoi = data['time']
        iteration = data['iteration']
        start = time()
        mode = data['mode'] # [(r,m)]
        activites = data['activites']
        sol_cca = data['cca']
        # ajouter la sequence comme sequence initiale
        self.planif = self.PlanificateurCCA(sol_cca)
        # solution initiale = sequence (activites des modes deja planifies)
        explication = self.planif.planifierMode(constellation,mode,activites)
        solution_cca = self.planif.getSolution()
        faisable = self.planif.sequenceFaisable(constellation)
        self.dates_cca = [self.start_calcul,time()]
        self.resultat = self.packSolution(solution_cca,mode,faisable,explication,iteration)
                
    def resoudre(self,constellation,mailbox):
        comm = MPI.COMM_WORLD
        while True:
            data = comm.recv()
            if self.estMessageFin(data):
                break
            else:
                self.insererActivites(constellation,data)
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

global config
config = Config()
instance = config.instance


mailbox = MessageriePartagee()

if rank==0 and config.glob.verbose:
    printColor("\n")
    printColor("Paramètres de résolution :",c='b')
    printColor("| instance : "+str(instance),c='b')
    printColor("| temps de calcul max : "+str(config.glob.limite_temps/60),c='b')
    printColor("\n")
comm.Barrier()

rd.seed(0)

if rank == 0:
    files = choseAndBroadcastFile(path,instance)
    for file in files:
        constellation = Constellation(file,False)
        comm.Barrier()
        if config.glob.verbose:
            printColor("Lecture OK.",c='g')
        process = Master()
        process.resoudre(constellation,mailbox)
        #process.solution.historique.tracerCPU()
        process.solution.verifierSolution(constellation)
        if config.glob.saveSample:
            process.saveSample(constellation)
        comm.Barrier()
        process.solution.writeSol(files[0].split('/')[-1])
        #f = process.tracerActivite(constellation)
        #process.tracerHistorique()
else:
    data = None
    data = comm.bcast(data,root=0)
    files = data['files']
    for file in files:
        constellation = Constellation(file,False)
        comm.Barrier()
        process = Slave()
        process.resoudre(constellation,mailbox)
        comm.Barrier()
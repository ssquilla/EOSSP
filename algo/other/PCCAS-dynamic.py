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
        assert(config.glob.anticipation_modes)
        self.a_supprimer.remove((s,o))
        
    # suppression du mode pour la version avec anticipation
    def supprimerAnticipation(self,r,m,grapheDep,liste_a_supprimer=None):
        self.terminerMode(r,m,False)
        for (s,o) in constellation.getRequete(r).getMode(m).getCouples():
            if liste_a_supprimer is None or o in liste_a_supprimer:
                if((s,o) not in self.a_supprimer):
                    cca = grapheDep.grapheActivites.getNoeud(o).getComposanteConnexe()
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
            printOpen('Echec du mode',(r,m),'cas 1 (transfert) : redémarrer ',activites_a_redemarrer,'{',c='m')
        else:
            assert(attenteRedemarrage)
            printOpen('Echec du mode',(r,m),'cas 2 (mode non prioritaire) : redémarrer',activites_a_redemarrer,'{',c='m')
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
        printClose('}',c='m')
        
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
            to_do_cca = grapheDep.grapheActivites.getNoeud(a).getComposanteConnexe()
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
        ccas = [grapheDep.grapheActivites.getNoeud(a).getComposanteConnexe() for a in self.activities_to_do[(r,m)]]
        ccas_libre = [cca for cca in ccas if cca not in self.cca_occupees]
        return ccas_libre
    
    def getCCA(self,r,m,grapheDep,cca):
        if (r,m) not in self.started_modes:
            self.started_modes.append((r,m))
        assert(cca not in self.cca_occupees)
        act = [a for a in self.activities_to_do[(r,m)] if grapheDep.grapheActivites.getNoeud(a).getComposanteConnexe()==cca]
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
        ccas = [grapheDep.grapheActivites.getNoeud(a).getComposanteConnexe() for a in self.activities_to_do[(r,m)]]
        ccas_libre = [cca for cca in ccas if cca not in self.cca_occupees]
        cca = ccas_libre[0]
        self.cca_occupees.append(cca)
        liste_a = []
        act = [a for a in self.activities_to_do[(r,m)] if grapheDep.grapheActivites.getNoeud(a).getComposanteConnexe()==cca]
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



class ParallelCCASearch(Solution):
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
        self.grapheDependances = GrapheDependancesDynamique(constellation,activites)
        printColor(self.grapheDependances,c='b')
        # si 1 coeur : slave local
        comm = MPI.COMM_WORLD
        if comm.Get_size()==1:
            self.localSlave = Slave()
        printColor("Durée d'initialisation :",time()-self.start_date,c='b')
  
    def ajouterModesCCAs(self,constellation,modes):
        activites_cca = {}
        modes_cca = {}
        for (r,m) in modes:
            self.ajouterModeCCASol(constellation,r,m)

    def ajouterActiviteCCASol(self,constellation,s,a,r,m):
        cca = self.grapheDependances.grapheActivites.getNoeud(a).getComposanteConnexe()
        if cca not in self.solCCAs[s]:
            self.solCCAs[s][cca] = SolCCA(cca,s)
        self.solCCAs[s][cca].ajouterActivite(constellation,a)
    
    def ajouterModeCCASol(self,constellation,r,m):
        for s in constellation.getRequete(r).getMode(m).getActivites():
            for o in constellation.getRequete(r).getMode(m).getActivites()[s]:
                self.ajouterActiviteCCASol(constellation,s,o,r,m)
    
    def _ajouterModeActifEtatRecherche(self,constellation,r,m,prevalidation=[]):
        activites = []
        for (s,o) in constellation.getRequete(r).getMode(m).getCouples():
            if o not in activites:
                activites.append(o)
        self.etat_recherche.ajouterModeActif(r,m,activites) 
        if len(prevalidation)>0:
            self.etat_recherche.prevaliderActivites(prevalidation,r,m)

    def ajouterModeActif(self,constellation,r,m,prevalidation=[]):
        self._ajouterModeActifEtatRecherche(constellation,r,m,prevalidation)
        self.ajouterModeCCASol(constellation,r,m)
        
    def presentCCA(self,constellation,r,m,cca,grapheActivites):
        for s in constellation.getRequete(r).getMode(m).getActivites():
            for o in constellation.getRequete(r).getMode(m).getActivitesSatellite(s):
                if grapheActivites.getNoeud(o).getComposanteConnexe()==cca:
                    return True
        return False
    
    def getComposanteActivite(self,a):
        return self.grapheDependances.grapheActivites.getNoeud(a).getComposanteConnexe()
    
    def setCCA(self,constellation,s,seq,cca,r,m):
        new = [a for a in seq if a not in self.solCCAs[s][cca].getSequence()]
        debut = constellation.getSatellite(s).getActivite(seq[0]).getDebut()
        fin = constellation.getSatellite(s).getActivite(seq[-1]).getFin()+config.glob.tau_max
        i_start = None
        i_end = None
        longueur_avant = len(self.solCCAs[s][cca].getSequence())
        for i,a in enumerate(self.solCCAs[s][cca].getSequence()):
            if a in seq and i_start is None:
                i_start = i
            if a not in seq and i_end is None:
                i_end = i-1
        if i_end is None:
            i_end = len(self.solCCAs[s][cca].getSequence())-1
        #printCore("SET CCA "+str(self.solCCAs[s][cca].getSequence())+str(self.solCCAs[s][cca].getSequence()[i_start:i_end])+str(seq))
        self.solCCAs[s][cca].setSequenceIndex(i_start,i_end,seq)

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
                        if o in self.grapheDependances.grapheActivites.getComposanteConnexe(cca).getNoeuds():
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
            
        
    def recalculerListeCCA(self,constellation,recalculs):
        if len(recalculs)>0:
            printColor("Division des CCA "+str(recalculs) +" "+str(self.grapheDependances.grapheActivites),c='y')
            getCCA = lambda a : self.grapheDependances.grapheActivites.getNoeud(a).getComposanteConnexe()
            solution = {s : {} for s in constellation.getSatellites()}
            sat_cca = {}
            # cca => satellite
            for ancienne_cca in recalculs:
                for s in self.solCCAs:
                    if ancienne_cca in self.solCCAs[s]:
                        sat_cca[ancienne_cca] = s
                        break
            for ancienne_cca in recalculs:
                sat_old = sat_cca[ancienne_cca]
                ccas = []
                copy = self.solCCAs[sat_old][ancienne_cca]
                for a in copy.getActivites():
                    c = getCCA(a)
                    if c not in ccas:
                        ccas.append(c)
                # cca => liste des activites    
                activites_cca = {}
                for a in copy.getActivites():
                    cca = getCCA(a)
                    if cca not in activites_cca:
                        activites_cca[cca] = []
                    activites_cca[cca].append(a)
                self.MAJPresenceRequetes(ancienne_cca,ccas)
                #printColor(str(ancienne_cca)+"=>"+str(ccas),c='y')    
                # recreer les ccas
                self.etat_recherche.divisionCCA(ancienne_cca,ccas)
                for c in ccas:
                    activites = activites_cca[c]
                    seq = []
                    for a in copy.getSequence():
                        if getCCA(a)==c:
                            seq.append(a)
                    solution[sat_old][c] = SolCCA(c,sat_old)
                    for a in self.grapheDependances.grapheActivites.getComposanteConnexe(c).getNoeuds():
                        solution[sat_old][c].ajouterActivite(constellation,a)
                    solution[sat_old][c].setSolution(seq)
            
            for cca in recalculs:
                del self.solCCAs[sat_cca[cca]][cca]
            for s in solution:
                for c in solution[s]:
                    self.solCCAs[s][c] = solution[s][c]
        
    def notifierDivisionCCA(self,constellation,ancienne_cca):
        printColor("division CCA",c='r')
        # MAJ sol CCA
        nouvelles_cca = self.s
        # MAJ Requete CCAs
        requetes = self.requetesCCA[ancienne_cca]
        del self.requetesCCA[ancienne_cca]
        cca_req = {r : [] for r in requetes}
        rec = lambda r : constellation.getRequete(r).getMode(self.modes_courants[r]).getRecompense() + self.bruit_requete[r]
        prio = lambda r : constellation.getRequete(r).getPriorite()
        for cca in nouvelles_cca:
            req = []
            for (p,w,r) in requetes:
                for rr,s,o in constellation.getRequete(r).getActivites():
                    if o in self.grapheDependances.grapheActivites.getComposanteConnexe(cca).getNoeuds():
                        reverse_insort(req,(prio(r),rec(r),r))
                        cca_req[r].append(cca)
                        break
            self.requetesCCA[cca] = req
        # MAJ CCA Requetes
        for r in cca_req:
            self.ccaRequetes[r].remove(ancienne_cca)
            self.ccaRequetes[r] += cca_req[r]
            
    def changementNomComposante(self,constellation,ancienne_composante,nouveau_nom=None):
        if nouveau_nom is None:
            new_cca = self.grapheDependances.grapheActivites.getNoeud(noeud).getComposanteConnexe()
        else:
            new_cca = nouveau_nom
        #printColor("changement nom "+str((ancienne_composante,nouveau_nom)),c='y')
        find = False
        for s in self.solCCAs:
            if find:
                break
            for cca in self.solCCAs[s]:
                if cca == ancienne_composante:
                    assert(new_cca not in self.solCCAs[s])
                    noeud = self.solCCAs[s][cca].getActivites()[0]
                    self.solCCAs[s][new_cca] = self.solCCAs[s][cca]
                    self.solCCAs[s][new_cca].renommer(new_cca)
                    del self.solCCAs[s][cca]
                    find = True
                    break
        if ancienne_composante in self.requetesCCA:
            for (p,w,r) in self.requetesCCA[ancienne_composante]:
                self.ccaRequetes[r].remove(ancienne_composante)
                self.ccaRequetes[r].append(new_cca)
            self.requetesCCA[new_cca] = self.requetesCCA[ancienne_composante]
            del self.requetesCCA[ancienne_composante]
        self.etat_recherche.changementNomCCA(ancienne_composante,new_cca)
        

    """
        =============================================== 
                        RESOLUTION
        =============================================== 
    """
   
    def validerMode(self,constellation,mode):
        printOpen("Validation du mode "+str(mode)+' {',c='g')
        start_val = time()
        # retirer la requete de toutes ses cca
        for cca in self.ccaRequetes[mode[0]]:
            self.retraitRequeteCCA(cca,mode[0])
        del self.ccaRequetes[mode[0]]
        
        # retirer les obs de la requete non presentes dans le mode
        retrait = []
        for s in constellation.getRequete(mode[0]).getActivites():
            for o in constellation.getRequete(mode[0]).getActivitesSatellite(s):
                presentGraphe = o in self.grapheDependances.grapheActivites.getNoeuds() # explications retirees...
                nonPresentMode = o not in constellation.getRequete(mode[0]).getMode(mode[1]).getActivitesSatellite(s)
                if presentGraphe and nonPresentMode:
                    retrait.append(o)
        #self.grapheDependances.grapheActivites.retirerListeActivites(retrait,constellation,self)
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
        printClose("}",c='g')
        
    def analyserMessage(self,data,constellation,freeCPUs,mailbox):
        start = time()
        self.historique.enregistrerReceptionMessage(data)
        i = data['source']
        assert(i not in freeCPUs)
        bisect.insort(freeCPUs,i)
        cpu_time,cca,succes,mode = data['cpu_time'],data['cca'],data['faisable'],data['mode']
        explication = data['explication']
        seq = cca.getSequence()
                
        # si message d'un mode passé : transfert potentiel d'activités validées ou à redemarrer
        retraits = self.supprimerAnciennesActivites(seq)
        if succes:
            mode = self.transfertSuccesActivites(explication,seq,mode)
        else:
            mode,activites_a_redemarrer = self.transfertEchecActivites(explication,retraits,mode)
        # indiquer la cca libre et supprimer les activités n'étant plus à programmer
        self.libererCCA(seq,explication)
                
        # traiter le message d'un mode en cours
        if succes and self.etat_recherche.getEtat(mode[0],mode[1])!=self.etat_recherche.FAILED:
            stat_validation = time()
            self.validerActivites(constellation,seq,mode)
            start_validation = time()
            termine = self.etat_recherche.validerActivites(explication,mode[0],mode[1])
            self.historique.enregistrerSucces(i,cpu_time)
            if termine:
                if config.glob.anticipation_modes:
                    self.validerAnticipation(constellation,mode)
                else:
                    self.validerMode(constellation,mode)
                self.objectif = self.calculerObjectif(constellation)
                self.historique.MAJHistorique(time()-self.start_date,0,self.objectif,None,self.modes_retenus,self.grapheDependances,constellation)
        else:
            if self.etat_recherche.getEtat(mode[0],mode[1])!=self.etat_recherche.FAILED:
                if config.glob.anticipation_modes:
                    self.notifierEchecAnticipation(constellation,mode,explication,activites_a_redemarrer)
                else:
                    self.planifierModeSuivant(constellation,mode,explication)
                self.historique.MAJHistorique(time()-self.start_date,0,self.objectif,None,self.modes_retenus,self.grapheDependances,constellation)
            self.historique.enregistrerEchec(i,cpu_time)
                    
    def analyserResultats(self,constellation,freeCPUs,mailbox):
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        start_iteration = time()
        if comm.Get_size()==1:
            messages = [self.localSlave.getResultat()]
        else:
            messages = mailbox.readMessages()
        for data in messages:
            if data is not None:
                self.analyserMessage(data,constellation,freeCPUs,mailbox)
                
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
        nouveau_mode = self.modes_courants[mode[0]]
        if nouveau_mode!=mode[1]:
            a_redemarrer = explication.copy()
            for a in retrait:
                if a in a_redemarrer:
                    a_redemarrer.remove(a)
            return (mode[0],nouveau_mode),a_redemarrer
        return mode,[]
    
    # voir si la planification d'une activité de l'ancien mode peut etre conservée pour le nouveau
    def transfertSuccesActivites(self,explication,seq,mode):
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
            cca = self.grapheDependances.grapheActivites.getNoeud(a).getComposanteConnexe()
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
            i_cca = self.grapheDependances.grapheActivites.getNoeud(a).getComposanteConnexe()
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
                prevalider.append(o)
            elif o in self.etat_recherche.activities_working[(r,ancien_mode)]:
                intersection_working.append(o)
        for o in act_ancien:
            if o not in act_nouveau_mode:
                retrait.append(o)
        
        return retrait,prevalider,intersection_working
    
    def planifierModeSuivant(self,constellation,mode,explication):
        printOpen('Echec du mode',mode,'cas 3 : annuler. cause :',explication,'{',c='r')
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
        printClose('}',c='r')

    def insererRequeteCCA(self,constellation,r):
        ccas = []
        if r in self.modes_courants:
            for (s,o) in constellation.getRequete(r).getCouples():
                try:
                    cca = self.grapheDependances.grapheActivites.getNoeud(o).getComposanteConnexe()
                    if cca not in ccas:
                        ccas.append(cca)
                except:
                    pass
            self.ccaRequetes[r] = ccas
            prio = constellation.getRequete(r).getPriorite()
            rec = constellation.getRequete(r).getMode(self.modes_courants[r]).getRecompense() + self.bruit_requete[r]
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
                               
    def calculerObjectif(self,constellation):
        return sum([constellation.getRequete(x[0]).getMode(x[1]).getRecompense() for x in self.modes_retenus])

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

    def initIteration(self,constellation):
        self.resetNoise(constellation)
        self.etat_recherche = EtatRecherche()
        # ajouter modes actifs
        self.explication = {}
        self.solCCAs = {s : {} for s in constellation.satellites} # liste de solutions CCA
        for cca in self.grapheDependances.getComposantes():
            s = self.grapheDependances.getSatelliteCCA(cca)
            self.solCCAs[s][cca] = SolCCA(cca,s)
            for a in self.grapheDependances.getCCA(cca).getNoeuds():
                self.solCCAs[s][cca].ajouterActivite(constellation,a)
        self.initRequetes(constellation)
        update_modes = []
        for r in constellation.getRequetes():
            m = constellation.getRequete(r).getModeCourant(constellation).getId()
            update_modes.append((r,m))
        for (r,m) in update_modes:
            self.ajouterModeActif(constellation,r,m) # etat recherche + CCA sol
        # stoquage etat des requetes dans les CCA
        self.requetesPresentesCCA(constellation) # cca => requetes,req => ccas 
        self.objectif = 0
            
    def resoudre(self,constellation,mailbox,afficher=True):
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
        while temps-self.start_date<config.glob.limite_temps:
            self.initIteration(constellation)
            while temps-self.start_date<config.glob.limite_temps and self.etat_recherche.modesEnCours():
                start_it = time()
                self.envoyerTache(constellation,freeCPUs,mailbox)
                self.analyserResultats(constellation,freeCPUs,mailbox)
                #self.verifierSolution(constellation)
                old_sol = self.getSolutionObservations()
                """
                    stockage de la solution courante
                """
                self.objectif = self.calculerObjectif(constellation)
                if last_record<( (time()-self.start_date)//10):
                    last_record = (time()-self.start_date)//10
                    # je mets solution à None : inutile de stoquer toutes les solutions intermédiaires 
                    self.historique.MAJHistorique(time()-self.start_date,0,self.objectif,None,self.modes_retenus,self.grapheDependances,constellation)
                    self.afficherInfo(time(),self.start_date,afficher,title="INFOS")
            self.historique.MAJHistorique(time()-self.start_date,1,self.objectif,None,self.modes_retenus,self.grapheDependances,constellation)
            self.afficherInfo(time(),self.start_date,afficher,title="FIN ITERATION "+str(super_iter),color='r')
            super_iter += 1
            temps = time()
            if config.glob.noise == 0 :
                break
        
        # attendre ce qu'il reste
        #for x in range(len([i for i in range(1,comm.Get_size()) if i not in freeCPUs])):
            #data = comm.recv()
            #fself.analyserMessage(data,constellation,freeCPUs,mailbox)
        
        self.nettoyerSolution(constellation)
        self.terminerProcessus()
        self.historique.MAJHistorique(time()-self.start_date,0,self.objectif,self.solCCAs,self.modes_retenus,self.grapheDependances,constellation)
        self.objectif,self.solution,self.modes_retenus = self.getBest()
        printColor(" ----------------------------------------------- [ END ] --------------------------------------------",c='y')
        self.afficherInfo(time(),self.start_date,afficher)
    
    def terminerProcessus(self):
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        size = comm.Get_size()
        for i in range(size):
            if i!=0:
                comm.send({'fin':True},dest=i)
    
    def calculerRequetesPresentes(self,constellation,cca,ccas):
        req = []
        activites = self.grapheDependances.grapheActivites.getComposanteConnexe(cca).getNoeuds()
        for a in activites:
            s = constellation.getSatelliteActivite(a)
            r = constellation.getSatellite(s).getActivite(a).getRequete()
            if r not in ccas:
                ccas[r] = []
            if cca not in ccas[r]:
                ccas[r].append(cca)
            m = self.modes_courants[r]
            rec = constellation.getRequete(r).getMode(m).getRecompense() + self.bruit_requete[r]
            p = constellation.getRequete(r).getPriorite()
            if (p,rec,r) not in req:
                reverse_insort(req,(p,rec,r))
        return req
    
    def requetesPresentesCCA(self,constellation):
        self.requetesCCA,self.ccaRequetes = {},{}
        req,ccas = {},{}
        for r in constellation.getRequetes():
            req[r] = []
            m = self.modes_courants[r]
            rec = constellation.getRequete(r).getMode(m).getRecompense() + self.bruit_requete[r]
            p = constellation.getRequete(r).getPriorite()
            for (s,o) in constellation.getRequete(r).getCouples():
                cca = self.grapheDependances.grapheActivites.getNoeud(o).getComposanteConnexe()
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
            r = self.requetesCCA[cca][0][2]
            m = self.modes_courants[r]
            if r not in requetes:
                requetes[r] = self.etat_recherche.getCCALibres(r,m,self.grapheDependances)
            if not self.etat_recherche.ccaOccupee(cca):
                ccas_libres.append(cca)
        return requetes,ccas_libres
    
    # renvoie les (r,cca) tels que r est premiere sur sa CCA et (r,cca) n'est pas encore calculée
    def candidatsTop(self,requetes,ccas_libres):
        candidats_top = []
        for cca in ccas_libres:
            r = max(self.requetesCCA[cca])[2]
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
            prio = lambda rm:constellation.getRequete(rm[1][0]).getPriorite()
            rec = lambda rm : constellation.getRequete(rm[1][0]).getMode(rm[1][1]).getRecompense()
            res = max(candidats_queue,key=lambda r:(-r[0],prio(r),rec(r))) # min niveau,max prio,max rec
            return res
        else:
            return None
        
    def choisirCandidatTop(self,candidats_top,constellation):
        if len(candidats_top)>0:
            prio = lambda rm : constellation.getRequete(rm[1][0]).getPriorite()
            rec = lambda rm : constellation.getRequete(rm[1][0]).getMode(rm[1][1]).getRecompense()
            res = max(candidats_top,key=lambda r :(prio(r),rec(r)))
            return res
        else:
            return None
        
    def choisirMode(self,constellation):
        #printColor([self.requetesCCA[cca][0] for cca in self.requetesCCA],c='y')
        start = time()
        # deteminer les ccas libres et les cca libres sur lesquelles les requetes ont du travail
        requetes,cca_libres = self.ccaLibres()
        start = time()
        candidats_top = self.candidatsTop(requetes,cca_libres)
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
    
    def envoyerTache(self,constellation,freeCPUs,mailbox):
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
                    data = {'mode':mode,'cca':deepcopy(cca),'source':0,'activites':activites,'time':time()}
                    self.localSlave.insererActivites(constellation,data)
                else:
                    mailbox.demanderTache(mode,activites,cca,cpu)
                #self.historique.enregistrerEnvoi(cpu)
            else:
                #printColor(len(freeCPUs),c='m')
                break
    
    def annulerModeSolution(self,constellation,mode):
        r,m = mode
        cca = []
        for s in constellation.getRequete(r).getMode(m).getActivites():
            for o in constellation.getRequete(r).getMode(m).getActivitesSatellite(s):
                cca_o = self.grapheDependances.grapheActivites.getNoeud(o).getComposanteConnexe()
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
                cca_o = self.grapheDependances.grapheActivites.getNoeud(o).getComposanteConnexe()
                if (s,cca_o) not in cca:
                    cca[(s,cca_o)] = []
                cca[(s,cca_o)].append(o)
        for (s,cca_a) in cca:
            self.solCCAs[s][cca_a].retirerListeActivites(constellation,cca[(s,cca_a)])
                               
    # active le mode pending suivant s'il y en a un et rajoute un pending mode s'il y en a encore
    def defilerEtatRequete(self,constellation,mode,explication,grapheDep):
        r,m = mode
        res = constellation.getRequete(r).getModeSuivant(explication,constellation)
        if res is not None:
            m = res.getId()
            self.modes_courants[r] = m
        # retrait : les activities presentes uniquement dans l'ancien mode. Les retirer
        # prevalider : presentes dans les deux modes et deja validées dans l'ancien => Prevalider dans le nouveau mode
        retrait,prevalider,intersection_working = self.intersectionActivites(constellation,r)
        modes_a_redemarrer = self.etat_recherche.echecCertain(r,m-1,grapheDep,retrait)
        self.retirerModeSolution(constellation,mode,retrait) # retirer de la solution
        self.grapheDependances.grapheActivites.retirerListeActivites(retrait,constellation,self)
        if res is not None:
            self.ajouterModeActif(constellation,r,m,prevalider) # prevalidation de l'intersection deja validee
            self.etat_recherche.preWorking(r,m,intersection_working)
        else:
            printColor('requete',r,'('+str(len(constellation.getRequete(r).getCouples())) +' activités) n\'a plus de modes',c='m')
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
        
    def resoudre(self,constellation,mailbox,afficher=True):
        # iterer insertion-réparation
        self.initSolution(constellation)
        self.solution.resoudre(constellation,mailbox,afficher)
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
        self.solution = ParallelCCASearch(constellation)
        
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
        
    def estMessageFin(self,data):
        return 'fin' in data
    
    def etendreDureesCalcul(self):
        durees = []
        for i,date in enumerate(self.dates_cca):
            if i<len(self.dates_cca)-1:
                durees.append((date,self.dates_cca[i+1]-date))
        return durees
    
    def packSolution(self,cca,mode,faisable,explication):
        return {'cpu_time':self.etendreDureesCalcul(),'reception':time(),'envoi':(self.date_envoi,self.duree_envoi),'mode':mode,'cca':deepcopy(cca),'faisable':faisable,'source':rank,'explication':explication}

    def insererActivites(self,constellation,data):
        self.duree_envoi = time() - data['time']
        self.start_calcul = time()
        self.date_envoi = data['time']
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
        self.resultat = self.packSolution(solution_cca,mode,faisable,explication)
                
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

def envoyerMessage(mode,activites,cca,cpu):
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    assert(rank==0)
    data = {'mode':mode,'cca':cca,'source':rank,'activites':activites,'time':time()}
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
    
    def demanderTache(self,mode,activites,cca,cpu):
        assert(MPI.COMM_WORLD.Get_rank()==0)
        envoyerMessage(mode,activites,cca,cpu)
        
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
    def demanderTache(self,mode,activites,cca,cpu):
        assert(MPI.COMM_WORLD.Get_rank()==0)
        self.cpu_mapped[cpu-1] = 1
        envoyerMessage(mode,activites,cca,cpu)
        
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
    
    def demanderTache(self,mode,activites,cca,cpu):
        assert(MPI.COMM_WORLD.Get_rank()==0)
        envoyerMessage(mode,activites,cca,cpu)

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

if rank == 0:
    files = choseAndBroadcastFile(path,instance)
    for file in files:
        constellation = Constellation(file,False)
        comm.Barrier()
        if config.glob.verbose:
            print("Lecture OK.")
        process = Master()
        process.resoudre(constellation,mailbox)
        process.solution.historique.tracerCPU()
        process.solution.verifierSolution(constellation)
        process.saveSample(constellation)
        comm.Barrier()
        #f = process.tracerActivite(constellation)
        #process.tracerHistorique()
        #plt.show()
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
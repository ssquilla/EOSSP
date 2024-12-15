from ..Utils.Utils import *
from ..Utils.config import *
config = Config()

from .solution_composantes import *
from .historique import Historique
from .historique import PREPROCESSING,NO_EVENT,END_ITERATION,END_RUN,BACKUP,END_OPERATEUR

from time import time
from mpi4py import *
from mpi4py import MPI
import pickle
import matplotlib
from copy import deepcopy,copy
from matplotlib import patches
import matplotlib.pyplot as plt
import os
#from bisect import insort_right
#from bisect import insort_left

class Solver:
    class Solution:
        def __init__(self,constellation,start_date):
            self.objectif = (0,0)
            self.modes_retenus = []
            self.solCCAs = {s : {} for s in constellation.satellites} # liste de solutions CCA
            self.historique = Historique(constellation,start_date)
        
        def resetHistorique(self,constellation,start_date):
            self.historique = Historique(constellation,start_date)
            
        def redemarrer(self,constellation):
            self.objectif = (0,0)
            self.modes_retenus = []
            self.solCCAs = {s : {} for s in constellation.satellites} # liste de solutions CCA
            
        def getModesRetenus(self):
            return self.modes_retenus
        
        def getSolCCA(self,s,cca):
            return self.solCCAs[s][cca]
        
        def getSolCCAs(self):
            return self.solCCAs
        
        def MAJObjectif(self,constellation):
            self.objectif = sum([constellation.getRequete(r).getMode(m).getRecompense() for (r,m) in self.modes_retenus]),sum([constellation.getRequete(r).getMode(m).getScoreTemporel() for (r,m) in self.modes_retenus])
            
        def remplacerMode(self,r,m,constellation):
            for i,(rr,mm) in enumerate(self.modes_retenus):
                if rr==r:
                    self.modes_retenus[i] = (r,m)
                    self.MAJObjectif(constellation)
                    return
            raise ValueError("Requête non présente")
            
        def ajouterModeRetenu(self,mode,constellation):
            r,m = mode
            if config.verifMode():
                assert(r not in [rr for (rr,mm) in self.modes_retenus])
            self.modes_retenus.append((r,m))
            self.MAJObjectif(constellation)
            
        def retirerRequete(self,r,constellation):
            for i,(rr,mm) in enumerate(self.modes_retenus):
                if rr == r:
                    self.modes_retenus.pop(i)
                    self.MAJObjectif(constellation)
                    return
            raise ValueError("requête non présente dans la solution : "+str(r))
    
        def retirerModeRetenu(self,mode,constellation):
            self.modes_retenus.remove(mode)
            self.MAJObjectif(constellation)
            
        def setModesRetenus(self,modes,constellation):
            self.modes_retenus = modes
            requetes = []
            for (r,m) in modes:
                if(r in requetes):
                    print(modes_a_retirer)
                    raise ValueError("Demande d'ajout de requêtes en double : ",r)
                requetes.append(r)
            
            self.MAJObjectif(constellation)
        
        def setSolCCAs(self,sol):
            self.solCCAs = sol
        
        def getObjectif(self):
            return self.objectif
        
        def ajouterInformationAdditionnelle(self,cle,valeur):
            self.historique.ajouterInformationAdditionnelle(cle,valeur)
        
        def __str__(self):
            msg = "---------- Solution ----------\n"
            msg += "| requêtes retenus : "+str(sorted([x[0] for x in self.modes_retenus]))+' ' +str(len(self.modes_retenus)) + ' modes\n'
            msg += "| objectif : " + str(self.objectif) +'\n'
            return msg
        
    def __init__(self,constellation,start,modeleDeTransition,dt_construction_transition,tlim=np.Inf,CCAs=None,solution=None):
        self.modeleDeTransition = modeleDeTransition
        self.tlim = min(tlim,config.getOptValue("time"))
        self.shiftDisplay = 1 # paramètre d'affichage
        self.start_date = start
        id_mpi = MPI.COMM_WORLD.Get_rank()
        seed_option = config.getOptValue("seed")
        self.generateur_aleatoire_perturbation = rd.Random(id_mpi+seed_option)
        global start_date
        start_date = self.start_date
        # initialisation des CCAs
        activites = constellation.extraireActivitesRequetesActives(None)
        
        if CCAs is not None:
            self.grapheDependances = CCAs
        else:
            self.grapheDependances = GroupeComposantesActivitesStatiques(constellation,activites,modeleDeTransition)
            printMaster(self.grapheDependances,c='b')
            printMaster("Durée des précalculs (CCA etc...) :",time()-self.start_date,c='g')
        
        # amorcage solution initiale
        if solution is not None:
            self.solution = solution
            # on indique que les plans earliest latest ne sont peut etre plus a jour
            # => un changement de modele de transition est possible
            for s in self.solution.getSolCCAs():
                for cca in self.solution.getSolCCAs()[s]:
                    self.solution.getSolCCA(s,cca).notifierPlanAJour(False,False)
                    if not modeleDeTransition.estTimeDependent():
                        self.plansCritiques(constellation,(s,cca)) # on remet a jour les plans
            self.verifierSolutionSiVerifMode(constellation,modeleDeTransition)
        else:
            self.solution = self.Solution(constellation,start)
            for (s,cca) in self.grapheDependances.getComposantes():
                self.setSolCCA(s,cca,SolCCA(cca,s))       
        self.solution.ajouterInformationAdditionnelle("duree de construction du modele de transition",dt_construction_transition)
        self.dt_construction_transition = dt_construction_transition
        
        if MPI.COMM_WORLD.Get_rank()==0:
            self.ajouterCommentaire(config.getOptValue("comm"))
        # a voir si on peut se passer de cette liste de candidates
        self.requetes_candidates = constellation.getToutesRequetes()
        # creer les modes courants
        self.modes_courants = {}
        if config.getOptValue("dynamic"):
            self.solution.historique.notifierArriveeRequetes(0,constellation.getRequetesDepart())
        shiftRightDisplay(self.shiftDisplay)
    
    
    def plansCritiques(self,constellation,id_cca,force=False,print_info=False):
        s,cca = id_cca
        if not self.getSolCCA(s,cca).plansAJour() or force:
            self.getSolCCA(s,cca).MAJPlansCritiques(constellation,self.modeleDeTransition,force=force)
        assert(len(self.getSolCCA(s,cca).getSequence())==len(self.getSolCCA(s,cca).planEst))
        if self.modeleDeTransition.estTimeDependent():
            if not(self.getSolCCA(s,cca).plansAJour()):
                raise OutdatedCriticalPlans
        
    def MAJPlansCritiquesSolution(self,constellation,force=False,print_info=False):
        for s in self.solution.solCCAs:
            for cca in self.solution.solCCAs[s]:
                self.plansCritiques(constellation,(s,cca),force=force,print_info=print_info)
            
    def getTempsEcoule(self):
        return time()-self.start_date
    
    def getTempsRestant(self):
        return min(config.getOptValue("time"),self.tlim) - self.getTempsEcoule()
    
    def MAJNouvellesRequetes(self,constellation,date,liste_requetes):
        self.solution.historique.notifierArriveeRequetes(date,liste_requetes)
        self.calculerRequetesPresentes(constellation)
    
    def invaliderMeilleureSolution(self):
        self.solution.historique.invaliderMeilleureSolution()
    
    def requetesPresentes(self,constellation,cca):
        requetes = []
        for a in self.grapheDependances.getActivitesComposante(cca):
            r = constellation.getRequeteActivite(a)
            if r not in requetes:
                bisect.insort_right(requetes,r)
        return requetes
    
    def calculerRequetesPresentes(self,constellation):
        self.requetes_cca = {}
        for cca in self.grapheDependances.getComposantes():
            self.requetes_cca[cca] = self.requetesPresentes(constellation,cca)
        self.requetes_en_commun = {}
    
    # les requetes presentes sur la cca meme si les activites ne sont pas dans la solution
    def getRequetesPresentesStatiques(self,cca):
        return self.requetes_cca[cca]
    
    # chercher les requetes presentes sur les deux cca meme si les activites ne sont pas dans la solution
    def requetesEnCommunStatique(self,cca1,cca2):
        couple = tuple(sorted((cca1,cca2)))
        if couple in self.requetes_en_commun:
            return self.requetes_en_commun[couple]
        else:
            r1 = self.getRequetesPresentesStatiques(cca1)
            r2 = self.getRequetesPresentesStatiques(cca2)
            if len(r1)<len(r2):
                res = [r for r in r1 if isInSortedList(r2,r)]
            else:
                res = [r for r in r2 if isInSortedList(r1,r)]
            self.requetes_en_commun[couple] = res
            return res
    
    def requetesPresentesDynamiques(self,constellation,id_cca):
        req = []
        s,cca = id_cca
        for a in self.getSolCCA(s,cca).getSequence():
            r = constellation.getRequeteActivite(a)
            if not isInSortedList(req,r):
                bisect.insort_right(req,r)
        return req
                    
    # calcule les requêtes qui peuvent être déplacées d'une cca à l'autre
    # dynaMique(cca1,cca2) = requetes présente dans la séquence de l'une avec des opportunites dans l'autre
    def requetesEnCommunDynamique(self,constellation,cca1,cca2):
        s1 = cca1[0]
        s2 = cca2[0]
        req_en_commun = []
        req_1 = []
        req_2 = []
        static1 = self.getRequetesPresentesStatiques(cca1)
        static2 = self.getRequetesPresentesStatiques(cca2)
        for a in self.getSolCCAs(s1,cca1).getSequence():
            r = constellation.getRequeteActivite(a)
            if not isInSortedList(req_1,r):
                bisect.insort_right(req_1,r)
                if not isInSortedList(req_en_commun,r) and r in static2:
                    bisect.insort_right(req_en_commun,r)
        for a in self.getSolCCAs(s2,cca2).getSequence():
            r = constellation.getRequeteActivite(a)
            if not isInSortedList(req_2,r):
                bisect.insort_right(req_2,r)
                if not isInSortedList(req_en_commun,r) and r in static1:
                    bisect.insort_right(req_en_commun,r)
        return req_en_commun,req_1,req_2
    
    def getModeIfPresent(self,r):
        for (rr,m) in self.solution.modes_retenus:
            if rr==r:
                return m
        return None
    
    def ajouterCommentaire(self,commentaire):
        self.solution.historique.ajouterCommentaire(commentaire)
        
    def MAJHistorique(self,event,constellation):
        #self.solution.objectif = self.calculerObjectifListeModes(constellation,modes_retenus)
        #print(self.objectif,modes_retenus)
        self.solution.MAJObjectif(constellation)
        mesure = time()-self.start_date
        self.solution.historique.MAJHistorique(mesure,event,self.solution.getObjectif(),self.solution.getSolCCAs(),self.solution.getModesRetenus(),self.grapheDependances,constellation)
     
    def notifierFinOperateur(self,constellation):
        if config.verifMode():
            printMaster("fin operateur",time()-self.start_date,force=True)
        
        self.MAJHistorique(END_OPERATEUR,constellation)        
     
    def notifierPreprocessing(self,constellation):    
        self.MAJHistorique(PREPROCESSING, constellation)
        
    def notifierFinIteration(self,constellation):
        self.MAJHistorique(END_ITERATION, constellation)

    def notifierSansEvenement(self,constellation):
        self.MAJHistorique(NO_EVENT,constellation)  

    def notifierFinExecution(self,constellation):
        self.MAJHistorique(END_RUN, constellation) 

    def notifierBackupSolution(self,constellation):
        self.MAJHistorique(BACKUP, constellation)              
    
    def supprimerDernierPoint(self):
        self.solution.historique.supprimerDernierPoint()   
    
    def gapRecompense(self,constellation):
        up = self.majorantRecompense(constellation)[0]
        return (up - self.objectif[0])/up<=config.glob.tolerance_opti
            
    def gapTemps(self,constellation):
        up2 = self.majorantRecompense(constellation)[1]
        return (up2 - self.objectif[1])/up2<=config.glob.tolerance_temps
    
    def majorantRecompense(self,constellation):
        return self.solution.historique.majorantRecompense(constellation)
    
    def initRequetes(self,constellation,noise=0,initFirstMode=True,conserve_old=True):
        if not conserve_old:
            self.requetes_candidates = list([r for r in constellation.getRequetes()])# if r not in self.requetesCouvertes()])
        else:
            self.requetes_candidates = list([r for r in constellation.getRequetes() if r not in self.requetesCouvertes()])
        for r in constellation.getToutesRequetes():
            constellation.getRequete(r).init = False
            if initFirstMode:
                mode = constellation.getRequete(r).getModeCourant(constellation)
                if mode is None:
                    raise ValueError("Mode requête",r,None)
                m = mode.getId()
                self.modes_courants[r] = m
            else:
                constellation.getRequete(r).resetModes(constellation,initFirstMode=False,conserve_old=conserve_old)
        # Trier les req candidates
        if initFirstMode:
            rec_courante = lambda r:constellation.getRequete(r).getMode(self.modes_courants[r]).getRecompense()
            prio = lambda r:constellation.getRequete(r).getPriorite()
            cle = lambda r : (prio(r),rec_courante(r))
            self.requetes_candidates = sorted(self.requetes_candidates,key=cle,reverse=True)
        
        # score candidats modes avant la planif
        modes = []
        for r in self.requetes_candidates:
            modes.append((r,0))
        self.resetNoise(constellation,noise)
        self.solution.historique.score_moyen_avant = constellation.scoreObsMoyenne(modes)
        self.solution.historique.scoreObsMoyenne = self.solution.historique.score_moyen_avant
        
        
    def resetNoise(self,constellation,amplitude):
        self.bruit_requete = {}
        for r in constellation.getToutesRequetes():
            #self.bruit_requete[r] = (np.random.rand()-0.5)*amplitude
            self.bruit_requete[r] = (self.generateur_aleatoire_perturbation.random()-0.5)*amplitude
            
    def resetScoreDestruction(self,constellation,req_couvertes):
        req_non_couvertes = []
        i = 0
        for r in constellation.getToutesRequetes():
            if r not in req_couvertes:
                req_non_couvertes.append((i,r))
                i += 1
        poids = []
        for i,r in req_non_couvertes:
            poids.append(constellation.getRequete(r).getModeCourant(constellation).getRecompense())
        rmax = len(req_non_couvertes)
        for i in range(rmax):
            (i,req_choisie) = rd.choices(req_non_couvertes,k=1,weights=poids)[0]
            self.bruit_requete[r] = -poids[i] + rmax-i
            poids[i] = 0
    
    def writeSol(self,instance):
        filename = config.donnees.algo_name+'-'+instance
        with open(filename,'w') as file:
            for s in sorted(self.getSolCCAs()):
                for cca in sorted(self.getSolCCAs()[s]):
                    file.write("CCA"+str(cca)+str(sorted(self.getSolCCA(s,cca).sequence))+'\n')
        
    def sortModes(self,modes): # modes = [(r,m,w)]
        cle = lambda rm : (rm[0]==-1,rm[2] + self.bruit_requete[rm[0]])
        return sorted(modes,key=cle,reverse=True)
        
    def saveSample(self,constellation,add=None):
        self.solution.historique.saveSample(constellation,add=add)
    
    def tracerObjectif(self):
        self.solution.historique.tracerObjectif()
        
    def tracerCPU(self):
        self.solution.historique.tracerCPU()
        
    def tracerChargeCCAs(self):   
        self.solution.historique.tracerChargeCCAs(self.grapheDependances)
        
    def creerMappingCCASlaves(self):
        if MPI.COMM_WORLD.Get_size()>1:
            self.cca_slaves = {cpu : [] for cpu in range(1,MPI.COMM_WORLD.Get_size())}
            tailles = []
            for cca in self.grapheDependances.getComposantes():
                taille =  self.grapheDependances.getTailleComposante(cca)
                reverse_insort(tailles,(taille,cca))
            cpu = 0
            for (taille,cca) in tailles:
                self.cca_slaves[cpu+1].append(cca)
                cpu = (cpu + 1) % (MPI.COMM_WORLD.Get_size()-1)
            for cpu in self.cca_slaves:
                MPI.COMM_WORLD.send({"cca":self.cca_slaves[cpu]},dest=cpu)
        else:
            self.cca_slaves = {1 : self.grapheDependances.getComposantes()}
     
    def mappingCCATaille(self,tailles):
        # tailles : liste (cca,nb activites a inserer)
        self.cca_slaves = {cpu : [] for cpu in range(1,MPI.COMM_WORLD.Get_size())}
        ccas = sorted(tailles,key=itemgetter(0),reverse=True)
        cpu = 0
        for (taille,cca) in ccas:
            self.cca_slaves[cpu+1].append(cca)
            cpu = (cpu + 1) % (MPI.COMM_WORLD.Get_size()-1)
    
    def redemarrer(self,constellation):
        self.solution.modes_retenus = []
        for (s,cca) in self.grapheDependances.getComposantes():
                self.setSolCCA(s,cca,SolCCA(cca,s))
        
        
    """
        =============================================== 
                        GETTER
        =============================================== 
    """
    
    def getModesRetenus(self):
        return self.solution.getModesRetenus()
    
    def ajouterModeRetenu(self,mode,constellation):
        self.solution.ajouterModeRetenu(mode,constellation)
        
    def retirerModeRetenu(self,mode,constellation):
        self.solution.retirerModeRetenu(mode,constellation)
        
    def getBest(self):
        return self.solution.historique.getBest()
    
    def restartBestSolution(self,constellation):
        best_objectif,best_solution,best_modes_retenus = self.solution.historique.getBest()
        self.solution.setSolCCAs(best_solution)
        self.solution.setModesRetenus(best_modes_retenus, constellation)
    
    def getSolutionObservations(self):
        sol = {s : [] for s in self.getSolCCAs()}
        for s in sol:
            ccas = sorted([self.getSolCCA(s,cca) for cca in self.getSolCCAs()[s]],key=lambda solcca : solcca.getDebut())
            for cca in ccas:
                for a in cca.getSequence():
                    sol[s].append(a)
        return sol   
    
    def extraireSequences(self):
        seq = {}
        for s in self.getSolCCAs():
            if s not in seq:
                seq[s] = []
            ccas = sorted([self.getSolCCA(s,cca) for cca in self.getSolCCAs()[s]],key=lambda solcca : solcca.getDebut())
            for cca in ccas:
                seq[s].append(cca.getSequence())
        return seq
    
    def getSolutionContainer(self):
        return self.solution
    
    def getSolution(self):
        sol = {s : [] for s in self.getSolCCAs()}
        for s in sol:
            ccas = sorted([self.getSolCCA(s,cca) for cca in self.getSolCCAs()[s]],key=lambda solcca : solcca.getDebut())
            for cca in ccas:
                for a in cca.getSequence():
                    sol[s].append(a)
        return sol            
    
    def getModes(self):
        return self.solution.modes_retenus
    
    def getPlan(self):
        try:
            return self.plan.copy()
        except Exception:
            return None
        
    def requetesCouvertes(self):
        return [x[0] for x in self.solution.getModesRetenus()]
    
    def retirerRequete(self,r):
        del self.requetes_candidates[r]
    
    def retirerRequeteCouverte(self,requete):
        indice = -1
        for i,(r,m) in enumerate(self.solution.modes_retenus):
            if r == requete:
                indice = i
                mode = (r,m)
                break
        if indice == -1:
            raise ValueError("Requête non présente")
        else:
            self.solution.modes_retenus.pop(indice)
            return mode
            
    def remplacerRequeteCouverte(self,requete):
        indice = -1
        for i,(r,m) in enumerate(self.solution.modes_retenus):
            if r == requete:
                indice = i
                break
        if indice == -1:
            self.solution.modes_retenus.append((r,m))
        else:
            self.solution.modes_retenus[indice] = (r,m)
    
    def setObjectif(self,objectif):
        self.solution.objectif = objectif
    
    def setModesRetenus(self,modes_retenus,constellation):
        self.solution.setModesRetenus(modes_retenus,constellation)
            
    """
        Version qui ne recalcule l'objectif. N'utiliser que si necessaire
    """
    def ecraserModes(self,modes_retenus):
        self.solution.modes_retenus = copy(modes_retenus)  
    
    def getObjectif(self,constellation=None,recompute=False):
        if recompute:
            assert(constellation is not None)
            return sum([constellation.getRequete(r).getMode(m).getRecompense() for (r,m) in self.solution.modes_retenus])
        else:
            return self.solution.objectif
        
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
    def planEarliestSequence(self,constellation,solution,s,modeleDeTransition=None):
        plan = []
        for i,activite in enumerate(solution):
            if(i==0):
                t = constellation.getSatellite(s).getActivite(activite).getDebut()
            else:
                prec = solution[i-1]
                duree = constellation.getSatellite(s).getActivite(prec).getDuree()
                start = constellation.getSatellite(s).getActivite(activite).getDebut()
                if modeleDeTransition is None:
                    transition = self.getTransition(constellation,s,prec,activite,t+duree,self.modeleDeTransition)
                else:
                    transition = self.getTransition(constellation,s,prec,activite,t+duree,modeleDeTransition)
                t = max(t + duree + transition,start)
            plan.append((activite,t))
        return plan
    
    def planEarliest(self,constellation,modeleDeTransition=None):
        solution = self.getSolution()
        self.plan = {s : [] for s in solution}
        for s in solution:
            for i,activite in enumerate(solution[s]):
                if(i==0):
                    t = constellation.satellites[s].getActivite(activite).getDebut()
                else:
                    prec = solution[s][i-1]
                    duree = constellation.getSatellite(s).getActivite(prec).getDuree()
                    start = constellation.satellites[s].getActivite(activite).getDebut()
                    if modeleDeTransition is None:
                        transition = self.getTransition(constellation,s,prec,activite,t+duree,self.modeleDeTransition)
                    else:
                        transition = self.getTransition(constellation,s,prec,activite,t+duree,modeleDeTransition)
                    t = max(t + duree + transition,start)
                self.plan[s].append((activite,t))
                
    def planLatest(self,constellation):
        assert(not self.modeleDeTransition.estTimeDependent())
        solution = self.getSolution()
        self.plan = {s : [] for s in solution}
        for s in solution:
            for i in range(len(solution[s])-1,-1,-1):
                a = solution[s][i]
                duree = constellation.getSatellite(s).getActivite(a).getDuree()
                if(i==len(solution[s])-1):
                    t = constellation.satellites[s].getActivite(activite).getFin() - duree
                else:
                    suivante = solution[s][i+1]
                    transition = self.getTransition(constellation,s,a,suivante,t+duree,self.modeleDeTransition)
                    die("TO DO")
                    t = min(t - duree - transition,constellation.satellites[s].getActivite(activite).getFin())
                self.plan[s].append((a,t))
                
    # deduit le temps le plus tot pour chaque activite
    def construirePlan(self,constellation,modeleDeTransition=None):
        #sequences = self.extraireSequences()
        self.planEarliest(constellation,modeleDeTransition=modeleDeTransition)
        if config.verifMode():
            self.verifierFaisabilite(constellation,modeleDeTransition=modeleDeTransition)
        
    def retardTotal(self):
        pi = 0
        solution = self.getSolution()
        for s in solution:
            for i,activite in enumerate(solution[s]):
                if(i==0):
                    t = constellation.satellites[s].getActivite(activite).getDebut()
                else:
                    prec = solution[s][i-1]
                    duree = constellation.getSatellite(s).getActivite(prec).getDuree()
                    start = constellation.satellites[s].getActivite(activite).getDebut()
                    if modeleDeTransition is None:
                        transition = self.getTransition(constellation,s,prec,activite,t+duree,self.modeleDeTransition)
                    else:
                        transition = self.getTransition(constellation,s,prec,activite,t+duree,modeleDeTransition)
                    t = max(t + duree + transition,start)
                pi_obs = max(0,t+duree-self.constellation.satellites[s].getActivite(a).getFin())
                pi += pi_obs
        return pi
    
    # renvoie un dict : satellites -> liste de (retard,flexibilite,activite) triée par retard décroissant
    def retardGaucheDroitePlanParObs(self,constellation,modeleDeTransition=None):
        retards = {}
        solution = self.getSolution()
        for s in solution:
            retard[s] = []
            for i,activite in enumerate(solution[s]):
                start = constellation.satellites[s].getActivite(activite).getDebut()
                duree = constellation.getSatellite(s).getActivite(prec).getDuree()
                if(i==0):
                    t = constellation.satellites[s].getActivite(activite).getDebut()
                    duree_prec = 0
                    transition = 0                
                else:
                    prec = solution[s][i-1]
                    if modeleDeTransition is None:
                        transition = self.getTransition(constellation,s,prec,activite,t+duree,self.modeleDeTransition)
                    else:
                        transition = self.getTransition(constellation,s,prec,activite,t+duree,modeleDeTransition)
                    t = max(t + duree + transition,start)
                droite = max(0,t-start)
                gauche = max(0,t+duree-constellation.satellites[s].getActivite(activite).getFin())
                reverse_insort(retards[s],(droite,gauche,activite))
        return retards
    

    # renvoie un dict : satellites -> liste de (retard,flexibilite,activite) triée par retard décroissant
    def retardGaucheDroiteSequenceParObs(self,constellation,sequence,s,modeleDeTransition=None):
        retards = []
        for i,activite in enumerate(sequence):
            start = constellation.satellites[s].getActivite(activite).getDebut()
            if(i==0):
                t = constellation.satellites[s].getActivite(activite).getDebut()
                duree_prec = 0
                transition = 0
            else:
                prec = sequence[i-1]
                duree_prec = constellation.getSatellite(s).getActivite(prec).getDuree()
                transition_left = self.getTransition(constellation,s,prec,activite,t+duree_prec,self.modeleDeTransition)
                t = max(t + duree_prec + transition_left,start)
            droite = max(0,t-start)
            gauche = max(0,t+duree_prec-constellation.satellites[s].getActivite(activite).getFin())
            reverse_insort(retards,(droite,gauche,activite))
        return retards    
    
    def faisabiliteTemporelle(self,constellation):
        self.planEarliest()
        for s in self.plan:
            for i,(p,t) in enumerate(self.plan[s]):
                if(i<len(self.plan[s])-1):
                    if(constellation.satellites[s].estObservation(p)):
                        duree = self.constellation.satellites[s].getObservation(p).getDuree()
                    else:
                        duree = self.dureeVidage(s,p) #self.vidages[s][p] 
                    transition = self.getTransition(constellation,s,p,self.plan[s][i+1][0],t+duree,self.modeleDeTransition)
                    if(t+transition+duree >self.plan[s][i+1][1] + 1e-5):
                        return False
                    if(t+duree>self.constellation.satellites[s].getActivite(p).getFin()):
                        return False
        return True

    def planEarliestSequence(self,constellation,solution,s,modeleDeTransition):
        plan = []
        for i,activite in enumerate(solution):
            if(i==0):
                t = constellation.getSatellite(s).getActivite(activite).getDebut()
            else:
                prec = solution[i-1]
                duree = constellation.getSatellite(s).getActivite(prec).getDuree()
                start = constellation.getSatellite(s).getActivite(activite).getDebut()
                if modeleDeTransition is None:
                    transition = self.getTransition(constellation,s,prec,activite,t+duree,self.modeleDeTransition)
                else:
                    transition = self.getTransition(constellation,s,prec,activite,t+duree,modeleDeTransition)
                t = max(t + duree + transition,start)
            plan.append((activite,t))
        return plan
    
    def getTransition(self,constellation,s,a1,a2,start,modeleDeTransition):
        if modeleDeTransition.estTimeDependent():
            transition = constellation.getSatellite(s).getTransitionTimeDependent(a1,a2,start,modeleDeTransition)
        else:
            transition = constellation.getSatellite(s).getTransition(a1,a2,modeleDeTransition)                    
        return transition
        
    def faisabiliteTemporelleSequence(self,constellation,sequence,s,modeleDeTransition=None):
        for i,activite in enumerate(sequence):
            if(i==0):
                t = constellation.getSatellite(s).getActivite(activite).getDebut()
            else:
                prec = sequence[i-1]
                duree = constellation.getSatellite(s).getActivite(prec).getDuree()
                start = constellation.getSatellite(s).getActivite(activite).getDebut()
                if modeleDeTransition is None:
                    transition = self.getTransition(constellation,s,prec,activite,t+duree,self.modeleDeTransition)
                else:
                    transition = self.getTransition(constellation,s,prec,activite,t+duree,modeleDeTransition)
                t = max(t + duree + transition,start)
                if t+constellation.getSatellite(s).getActivite(activite).getDuree() > constellation.getSatellite(s).getActivite(activite).getFin():
                    return False
        return True        
    
    """
        =============================================== 
                        RESOLUTION
        =============================================== 
    """
    def getGrapheComposantes(self):
        return self.grapheDependances
    
    """
        afficher les informations de la solution courante.
        filtres possibles : time,objectif,modes,ratio,best,repartition,obs,requete,size
        title : indiquer un titre
        temps : date courante
        start_date : date de début d'execution
        constellation : la constellation
        add = dictionnaire de messages optionnels (titre,contenu)
    
    """
    def afficherInfo(self,temps,start_date,constellation,core=None,color='c',add={},title='Info',filtre=None):
        if core is None or MPI.COMM_WORLD.Get_rank()==core:
            try:
                self.setObjectif(self.calculerObjectif(constellation))
            except Exception as e:
                pass
            shift = getDisplayDepth()
            shiftLeftDisplay(shift-1)
            self.filtre = filtre
            self.solution.historique.filtre = filtre
            title_size_max = 30
            title_string = title[0:max(len(title),title_size_max)]
            title_before = (title_size_max - len(title_string))//2
            title_after = title_size_max - len(title_string) - title_before
            center = 4 + title_size_max 
            largeur = 120
            left = (largeur-center)//2
            right = largeur - left - center
            printColor("\n")
            title_msg = "[" + " "*title_before + title + title_after*" " + " ]"
            left_msg = " " + (left-2)*"-" + " " 
            right_msg = " " + (right-1)*"-"
            printColor( left_msg + title_msg +  right_msg ,c=color)
            for text in add:
                printColor("| "+str(text)+" : " + str(add[text]),c=color)
            if filtre is None or 'time' in filtre:
                printColor("| temps écoulé : ",(temps-start_date)//60,"minutes",(temps-start_date)%60,c=color)
            if filtre is None or 'obj' in filtre:
                printColor("| objectif : ",self.getObjectif(),c=color)
            if filtre is None or "modes" in filtre:
                printColor("| modes retenus : ",len(self.solution.modes_retenus),c=color)
            lines = str(self.solution.historique).split('\n')
            for line in lines:
                if len(line)>0:
                    printColor(line,c=color)
            printColor(' ' + (largeur-1)*"-",c=color)
            shiftRightDisplay(shift-1)
            
    def insererSequences(self,constellation,sequences):
        requetes_retenues = [x[0] for x in self.solution.modes_retenus]
        for s,cca in sequences:
            seq = [a for a in sequences[cca] if constellation.getRequeteActivite(a) in requetes_retenues]
            #print("set sequence",s,cca)
            self.getSolCCA(s,cca).setSequence(constellation,seq)
    
    def intersectionActivites(self,constellation,r,composantes):
        nouveau_mode = self.modes_courants.get(r)
        ancien_mode = nouveau_mode - 1
        act_ancien = [x[1] for x in constellation.getRequete(r).getMode(ancien_mode).getCouples()]
        if nouveau_mode is None:
            return act_ancien,[],[]
        act_nouveau_mode = [x[1] for x in constellation.getRequete(r).getMode(nouveau_mode).getCouples()]
        retrait = []
        prevalider = {}
        intersection_working = {}
        for o in act_nouveau_mode:
            cca = composantes.find(o)
            if o in act_ancien and o in self.etat_recherche.activities_done[(r,ancien_mode)]:
                if cca not in prevalider:
                    prevalider[cca] = []
                prevalider[cca].append(o)
            elif o in self.etat_recherche.activities_working[(r,ancien_mode)]:
                if cca not in intersection_working:
                    intersection_working[cca] = []
                intersection_working[cca].append(o)
        for o in act_ancien:
            if o not in act_nouveau_mode:
                retrait.append(o)
        
        return retrait,prevalider,intersection_working
    """
    def calculerObjectif(self,constellation):
        return self.calculerObjectifListeModes(constellation,self.solution.modes_retenus)
    
    def calculerObjectifListeModes(self,constellation,modes):
        return sum([constellation.getRequete(r).getMode(m).getRecompense() for (r,m) in modes]),sum([constellation.getRequete(r).getMode(m).getScoreTemporel() for (r,m) in modes])
   """
    def terminerProcessus(self):
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        size = comm.Get_size()
        for i in range(size):
            if i!=0:
                comm.send({'fin':True},dest=i)
        comm.Barrier()
    
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
    
    def recompenseBruitee(self,constellation,r,m):
        return (constellation.getRequete(r).getMode(m).getRecompense() + self.bruit_requete[r],constellation.getRequete(r).getMode(m).getScoreTemporel())
    
    def mappingCCAObservationsCourantes(self,constellation,failed):
        mapping_obs = {}
        for r in self.requetes_candidates:
            m = constellation.getRequete(r).getModeCourant().getId()
            w = constellation.getRequete(r).getModeCourant().getRecompense()
            for (s,o) in constellation.getRequete(r).getMode(m).getCouples():
                cca = self.grapheDependances.getActiviteCCA(o)
                if cca not in mapping_obs:
                    mapping_obs[cca] = {}
                if (r,m,w) not in mapping_obs[cca]:
                    mapping_obs[cca][(r,m,w)] = []
                mapping_obs[cca][(r,m,w)].append(o)
        
        mapping_obs,sequences_cca = self.filtrer(constellation,mapping_obs,failed)    
        return mapping_obs,sequences_cca
        
    # si retrait est None on retire toutes les activites du mode. Sinon on retire la liste 'retrait'
    def retirerModeSolution(self,constellation,mode,retrait=None):
        r,m = mode
        cca = {}                     
        for (s,o) in constellation.getRequete(r).getMode(m).getCouples():
            if retrait is None or o in retrait:
                cca_o = self.grapheDependances.getActiviteCCA(o)
                if (s,cca_o) not in cca:
                    cca[(s,cca_o)] = []
                cca[(s,cca_o)].append(o)
        for (s,cca_a) in cca:
            self.getSolCCA(s,cca_a).retirerListeActivites(constellation,cca[(s,cca_a),self.modeleDeTransition])
                              
    def supprimerRequete(self,constellation,r):
        if r in self.modes_courants:
            del self.modes_courants[r]
        if r in self.requetes_candidates:
            self.requetes_candidates.remove(r)
            for cca in self.ccaRequetes[r]:
                self.retraitRequeteCCA(cca,r)
            del self.ccaRequetes[r]

        
    """
        =============================================== 
                        GRAPHIQUES
        =============================================== 
    """
    def tracerHistorique(self,init=True):
        return self.solution.historique.tracer(init)
    
    def tracerActivite(self,constellation,annoter=False):
        f,ax = plt.subplots(figsize=(15,6))
        couleurObs = 'c'
        couleurVid = 'g'
        couleurTrans = 'r'
        neutre = 'k'
        for s in self.plan:
            
            #horizon obs
            for p in constellation.satellites[s].getObservations():
                debut = constellation.satellites[s].getObservation(p).getDebut()
                fin = constellation.satellites[s].getObservation(p).getFin()
                ax.plot([debut,fin],[s,s],'-'+couleurObs,alpha=1,linewidth=2)
            for p in constellation.satellites[s].getVidages():
                debut = constellation.satellites[s].getVidage(p).getDebut()
                fin = constellation.satellites[s].getVidage(p).getFin()
                ax.plot([debut,fin],[s,s],'-'+couleurVid,alpha=1,linewidth=2)
            
            for i,(p,t) in enumerate(self.plan[s]):
                #ax.scatter([t],[s],c=neutre,marker="+")
                if(constellation.satellites[s].estVidage(p)):
                    couleur = couleurVid
                    duree = self.vidages[s][p]
                    lab = 'vidage'
                else:
                    couleur = couleurObs
                    duree = constellation.satellites[s].getObservations()[p].getDuree()
                    lab = 'observation'
                ax.add_patch(patches.Rectangle((t,s-1/2),duree,1,color=couleur,fill=True,label=lab))
                #ax.plot([t,t+duree],[s,s],color=couleur,label=lab,linewidth=50)
                
                if i<len(self.plan[s])-1:
                    transition = self.getTransition(constellation,s,p,self.plan[s][i+1][0],t+duree,self.modeleDeTransition)
                    ax.add_patch(patches.Rectangle((t+duree,s-1/2),transition,1,color=couleurTrans,fill=True,label='transition'))
                    #ax.plot([t+duree,t+duree+transition],[s,s],'--'+couleurTrans,label='transition',linewidth=50)
                if(annoter):
                    ax.annotate('{}'.format(self.sol[s][i][0]),xy=(self.sol[s][i][1], s),xytext=(0, 3),textcoords="offset points",ha='center', va='bottom',color='k')           
        
        plt.yticks(sorted(list(constellation.satellites.keys())))
        if(annoter):
            plt.xlabel("time (s)")
            plt.ylabel("satellite id")
        plt.grid(axis='y')
        ax.xaxis.set_tick_params(labelsize=15)
        ax.yaxis.set_tick_params(labelsize=15)
        a,b,modes = self.solution.historique.getBest()
        file = config.getOptValue("instance")["file"]
        folder = config.getOptValue("instance")["folder"]
        algo = config.donnees.algo_name
        #plt.title(config.glob.filename +  " - " +str(config.glob.size) + ' processes - Plan : '+str(len(modes))+" requests")
        
        plt.savefig('results/plan/plan_'+folder+"_"+file+"_"+algo+"_"+str(config.glob.size)+".png")
        #legend_without_duplicate_labels(ax)
        return f
            
    def tracerPertesModes(self):
        rangs = {}
        perte = []
        f,(ax1,ax2) = plt.subplots(1,2)
        for r,m in self.solution.modes_retenus:
            tri = sorted([self.constellation.getRequete(r).getMode(mm).getRecompense() for mm in self.constellation.getRequete(r).getModes()],reverse=True)
            best = tri[0]
            for i,rec in enumerate(tri):
                if(rec==self.constellation.getRequete(r).getMode(m).getRecompense()):
                    if i not in rangs:
                        rangs[i] = 1
                    else:
                        rangs[i] = rangs[i]+1
                    perte.append(100*(best-rec)/best)
                    break
        ax1.set(xlabel='rank of selected modes', ylabel='count')
        ax1.set_title("rank of selected modes - "+str(len(self.solution.modes_retenus))+"/"+str(len(self.constellation.getToutesRequetes())) + " covered")
        ax2.set(xlabel='reward loss compared to the best mode (%)', ylabel='count')
        ax2.set_title("reward loss of selected modes")
           
        n1,bins1,p1 = ax1.hist(rangs,bins=np.arange(0,max(rangs)))
        n2,bins2,p2 = ax2.hist(perte)
        
        #print("LATEX")
        #print(histToLatex(n1,bins1,"mode rank","mode count"))
        #print(histToLatex(n2,bins2,"loss compared to the best mode (\%)","mode count"))
        return f
    
    
    """
        =============================================== 
                        VERIF/DEBUG
        =============================================== 
    """
    def verifierListes(self):
        req_retenus = [x[0] for x in self.solution.modes_retenus]
        assert(len(np.unique(req_retenus))==len(req_retenus))
        for requete in self.requetes_candidates:
            assert(requete not in [x[0] for x in self.solution.modes_retenus])
        for requete in [x[0] for x in self.solution.modes_retenus]:
            assert(requete not in [x[0] for x in self.modes_candidats])
        for (r,m) in self.solution.modes_retenus:
            assert((r,m) not in self.modes_candidats)
            assert((r,m) not in [x[0] for x in self.modes_retires[r]])
        for r in self.requetes_candidats:
            for m in self.getRequete(r).modes_candidats:
                assert((r,m) not in [x[0] for x in self.modes_retires[r]])
        
        for r in self.modes_retires:
            for ((rr,m),i) in self.modes_retires[r]:
                if(r in self.requetes_candidates):
                    for m in self.getRequete(r).modes_candidats:
                        print((r,m))
                assert(m not in self.modes_candidats[r])
    
    def verifierMode(self,constellation,r,m):
        for s,o in constellation.getRequete(r).getMode(m).getCouples():
            _,cca = self.grapheDependances.getActiviteCCA(o)
            if not o in self.getSolCCA(s,cca).getSequence():
                print(self.getSolCCA(s,cca),o)
                comm = MPI.COMM_WORLD
                rank = comm.Get_rank()
                coeur = "coeur : "+str(rank)
                raise ValueError(coeur,"time : ",time()-self.start_date,"activité manquante",s,o,'mode',(r,m))
        if r != -1:
            if not constellation.getRequete(r).getMode(m).acceptable(constellation):
                raise ValueError(coeur,"time : ",time()-self.start_date,"contenu du mode non acceptable (construction non consistante ?)",constellation.getRequete(r).getMode(m))
    
    def verifierPresenceModeVidage(self,constellation):
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        coeur = "coeur : "+str(rank)        
        for (s,d) in constellation.getRequete(-1).getCouples():
            _,cca = self.grapheDependances.getActiviteCCA(d)
            if d not in self.getSolCCA(s,cca).getSequence():
                raise ValueError(coeur,"time : ",time()-self.start_date,"Vidage non présente",d,self.getSolCCA(s,cca).getSequence())
        
    def verifierModes(self,constellation):
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        coeur = "coeur : "+str(rank)         
        self.verifierPresenceModeVidage(constellation)
        for (r,m) in self.solution.modes_retenus:
            self.verifierMode(constellation,r,m)
        for s in self.getSolCCAs():
            for cca in self.getSolCCAs()[s]:
                for o in self.getSolCCA(s,cca).getSequence():
                    find = False
                    for (r,m) in self.solution.modes_retenus:
                        if (s,o) in constellation.getRequete(r).getMode(m).getCouples():
                            find = True
                            break
                    if not find:
                        r = constellation.getRequeteActivite(o)
                        typ = constellation.getRequete(r).getType()
                        print("mode courant :",constellation.getRequete(r).getModeCourant(constellation))
                        print("modes de la requete")
                        for m in constellation.getRequete(r).modes:
                            print(constellation.getRequete(r).getMode(m))
                        raise ValueError(coeur,"time : ",time()-self.start_date,'activité sans mode',s,o,cca,typ,r,".Requête couverte :",r in self.requetesCouvertes())
    
    def verifierCCA(self,constellation):
        for s in self.getSolCCAs():
            for cca in self.getSolCCAs()[s]:
                #printColor(self.getSolCCA(s,cca).debut,self.getSolCCA(s,cca).fin,c='m')
                faisable = self.getSolCCA(s,cca).sequenceFaisable(constellation,self.modeleDeTransition)
                if not faisable:
                    die('CCA infaisable',self.getSolCCA(s,cca),"retard=",self.getSolCCA(s,cca).retardSequence(constellation,self.getSolCCA(s,cca).sequence,self.modeleDeTransition))
                
    def verifierUniciteModes(self):
        unique = []
        for (r,m) in self.solution.modes_retenus:
            if(r in unique):
                raise ValueError('requête couverte par plusieurs modes')
            else:
                unique.append(r)
    
    def verifierFaisabilite(self,constellation,modeleDeTransition=None):
        self.verifierCCA(constellation)
        for s in self.plan:
            transitions = []
            if modeleDeTransition is None:
                modeleDeTransition = self.modeleDeTransition
            for i,(p,t) in enumerate(self.plan[s]):
                if(i<len(self.plan[s])-1):
                    duree = constellation.satellites[s].getActivite(p).getDuree()
                    transition = self.getTransition(constellation,s,p,self.plan[s][i+1][0],t+duree,modeleDeTransition)
                    transitions.append(transition)
                    if(t+transition+duree >self.plan[s][i+1][1] + 1e-5):
                        #print(t+transition+duree,self.plan[s][i+1][1],s,i)
                        raise ValueError('transition debordante')
                    if(t+duree>constellation.satellites[s].getActivite(p).getFin()):
                        #printColor("transitions : ",transitions)
                        print("retard :",t+duree-constellation.satellites[s].getActivite(p).getFin())
                        debut = constellation.satellites[s].getActivite(p).getDebut()
                        fin = constellation.satellites[s].getActivite(p).getFin()
                        printColor("activité", (s,p),"date",t,'durée',duree,"fenetre",debut,fin,c='r')
                        print([(p,self.grapheDependances.getActiviteCCA(p),t) for (p,t) in self.plan[s]])
                        raise ValueError('retard fenetre : les CCA sont probablement dans le désordre')

    def verifierRequetesActives(self,constellation):
        if config.getOptValue("dynamic"):
            for (r,m) in self.solution.modes_retenus:
                if not constellation.getRequete(r).estActif():
                    raise ValueError("Requête inactive")
        
    def verifierSolutionSiVerifMode(self,constellation,modeleDeTransition=None):
        if config.getOptValue("verif"):
            self.verifierSolution(constellation,modeleDeTransition=modeleDeTransition)
        
    def verifierSolution(self,constellation,modeleDeTransition=None):
        self.construirePlan(constellation,modeleDeTransition=modeleDeTransition)
        self.verifierRequetesActives(constellation)
        self.verifierModes(constellation)
        #self.verifierUniciteModes()

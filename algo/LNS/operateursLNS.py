#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 10:31:30 2022

@author: ssquilla
"""

from ..model.solution import *
from ..model.constellation import *
from ..model.components import *
from ..model.componentPlan import *
from ..Utils.Utils import *
from ..Utils.Utils import printColor,printOpen,printClose,alert,warn
from ..Utils.Utils import shiftRightDisplay,shiftLeftDisplay
from ..Utils.config import Config

global config
config = Config()

import random as rd
import math
import docplex.cp as cp
from docplex.cp.model import CpoModel,interval_var,binary_var


class Operateur:
    def __init__(self,name,modeleDeTransition,k=2):
        self.k=k
        self.modeleDeTransition = modeleDeTransition
        self.name = name
        self.liste_tabu = {} # pairs de cca tabu
        self.cca_observations = None
        self.connectivites_cca = None
        id_mpi = MPI.COMM_WORLD.Get_rank()
        seed_option = config.getOptValue("seed")
        self.randomInstanceChoixCCA = rd.Random(id_mpi+seed_option)
        
    def acceptable(self,paircca):
        pair_triee = tuple(sorted(paircca))
        return pair_triee not in self.liste_tabu
    def getName(self):
            return self.name
    def MAJTabu(self):
            del_pair = []
            for couple in self.liste_tabu:
                self.liste_tabu[couple] -= 1
                if self.liste_tabu[couple]==0:
                    del_pair.append(couple)
            for couple in del_pair:
                del self.liste_tabu[couple]
    def __str__(self):
        return 'Opérateur '+self.name
    def _initialiserDonneesCCA(self,LNS,constellation):
        self.cca_observations = LNS.getCCAObservations()
        self.calculerConnectivite(LNS)
        
    def checkCCAVide(self,LNS,nouvelles_activites,CCAs):
         CCAs_obs = []
         for r in nouvelles_activites:
             for m in nouvelles_activites[r]:
                 for a in nouvelles_activites[r][m]:
                     id_cca = LNS.grapheDependances.getActiviteCCA(a)
                     if id_cca not in CCAs_obs:
                         assert(id_cca in CCAs)
                         CCAs_obs.append(id_cca)
                         if len(CCAs)==len(CCAs_obs):
                             assert(sorted(CCAs)==sorted(CCAs_obs))
                             return False
         return len(CCAs)!=len(CCAs_obs)
    
    def calculerConnectivite(self,LNS):
        self.connectivites_cca = {}
        for cca1 in LNS.grapheDependances.getComposantes():
            for cca2 in LNS.grapheDependances.getComposantes():
                if cca1<cca2:
                    if cca1 not in self.connectivites_cca:
                        self.connectivites_cca[cca1] = {}
                    self.connectivites_cca[cca1][cca2] = LNS.requetesEnCommunStatique(cca1,cca2)
    
    def MAJNouvellesRequetes(self,LNS):
        if config.getOptValue("n_cca")>1:
            self.calculerConnectivite(LNS)
        self.cca_observations = LNS.getCCAObservations()
        
    def getConnectivites(self,cca1,cca2):
        c1,c2 = sorted((cca1,cca2))
        return self.connectivites_cca[c1][c2]
    
    def _choix1ereCCA(self,LNS,constellation):
        # choix 1ere CCA
        if config.getOptValue("nbh")=="load":
            mesures_cca = [len(LNS.requetesPresentesDynamiques(constellation,id_cca)) for id_cca in self.cca_observations]
        elif config.getOptValue("nbh")=="random":
            mesures_cca = [1 for id_cca in self.cca_observations]
        else:
            raise ValueError("type de voisinage CP inconnu : "+config.getOptValue("nbh"))
        return self.randomInstanceChoixCCA.choices([id_cca for id_cca in self.cca_observations], weights=mesures_cca,k=1)[0]    
    
    def _initRechercheCCAs(self,LNS,constellation):
        if self.connectivites_cca is None:
            self._initialiserDonneesCCA(LNS, constellation)
    
    def _choixCCASuivante(self,LNS,CCAs,candidats,requetes_presentes):    
        assert(len(CCAs)>0)
        weights = [sum([len(self.getConnectivites(cca2,id_cca)) for cca2 in CCAs]) for id_cca in candidats]
        if sum(weights)==0:
            return False  # alors on ne trouvera personne
        last = CCAs[-1]
        next_cca = self.randomInstanceChoixCCA.choices(candidats,weights=weights,k=1)[0]
        assert(next_cca not in CCAs)
        CCAs.append(next_cca)
        candidats.remove(next_cca)
        for r in LNS.getRequetesPresentesStatiques(next_cca):
            if not isInSortedList(requetes_presentes,r):
                bisect.insort_left(requetes_presentes, r)
        return True
        
    # Long la 1ere fois : initialisation des connées sur les CCA
    def selectionnerCCAs(self,LNS,constellation):  
       t1 = time()
       self._initRechercheCCAs(LNS, constellation)
       t2 = time()
       CCAs = [self._choix1ereCCA(LNS,constellation)]
       t3 = time()
       requetes_presentes = sorted((LNS.getRequetesPresentesStatiques(CCAs[0])))
       succes = True
       t4 = time()
       candidats = [id_cca for id_cca in self.cca_observations if id_cca not in CCAs]
       t5 = time()
       while (succes and len(CCAs)<self.k):
           succes = self._choixCCASuivante(LNS,CCAs,candidats,requetes_presentes)
       t6 = time()
       #print("test",t2-t1,t3-t2,t4-t3,t5-t4,t6-t5)
       return CCAs,requetes_presentes
# Operateur qui ne fait rien.
class OperateurPassif(Operateur):
    def __init__(self):
        super().__init__("Passif",None)
    def appliquer(self,LNS,constellation):
        pass

class OperateurReparation(Operateur):
        def __init__(self,name,modeleDeTransition,k=2):
            super().__init__(name,modeleDeTransition,k=k)
class RemplissageGloutonGlobal(Operateur):
        def __init__(self,modeleDeTransition):
            super().__init__("constructeur glouton",modeleDeTransition)
        
        def appliquer(self,LNS,constellation):
            old_objectif = LNS.getObjectif(constellation)
            LNS.greedyFill(constellation,limit_req=config.getOptValue("use_solver"),forbid_solver=False)
            return LNS.getObjectif(constellation)>old_objectif
class RemplissageGloutonGlobalAccelere(Operateur):
        def __init__(self,modeleDeTransition):
            super().__init__("constructeur glouton",modeleDeTransition)
        
        def appliquer(self,LNS,constellation):
            old_objectif = LNS.getObjectif(constellation)
            LNS.greedyFill(constellation,limit_req=config.getOptValue("use_solver"),forbid_solver=False)
            return LNS.getObjectif(constellation)>old_objectif
class OperateurDestruction(Operateur):
        def __init__(self,name,modeleDeTransition,k=2):
            super().__init__(name,modeleDeTransition,k=k)
        def detruire(self,LNS,constellation):
            raise ValueError("non implémenté")
        def retirerModes(self,LNS,constellation,modes_a_retirer):
            len_avant= len(LNS.getModesRetenus())
            LNS.setModesRetenus([(r,m) for (r,m) in LNS.getModesRetenus() if (r, m) not in modes_a_retirer],constellation)
            len_apres = len(LNS.getModesRetenus())
            assert(len_avant-len_apres==len(modes_a_retirer))
            activites_a_retirer = [(s,a) for (r,m) in modes_a_retirer for (s,a) in constellation.getRequete(r).getMode(m).getCouples()]
            activites_par_cca = {}
            for (s,a) in activites_a_retirer:
                id_cca = LNS.grapheDependances.getActiviteCCA(a)
                if id_cca not in activites_par_cca:
                    activites_par_cca[id_cca] = []
                activites_par_cca[id_cca].append(a)
            for id_cca in activites_par_cca:
                s,cca = id_cca
                LNS.getSolCCA(s,cca).retirerListeActivites(constellation,activites_par_cca[id_cca],self.modeleDeTransition)            

class Perturbateur(OperateurDestruction):
    def __init__(self,nom,modeleDeTransition):
        super().__init__(nom,modeleDeTransition)
    def appliquer(self,LNS,constellation):
        raise ValueError("classe abstraite")
class DestroyAndRepairGreedyRequest(Perturbateur):
        def __init__(self,modeleDeTransition,destroy_rate,accept_worse_solution):
            super().__init__("Destroy and repair greedy",None) # None : pas besoin de stoquer un doublon
            self.accept_worse_solution = accept_worse_solution
            self.destroy = DestructeurRequete(modeleDeTransition,destroy_rate)
            self.remplissage = RemplissageGloutonGlobal(modeleDeTransition)
            
        def appliquer(self,LNS,constellation):
            #printOpen("Perturbation :",c='y')
            old_objectif = LNS.getObjectif(constellation)
            if not self.accept_worse_solution:
                printOpen("Création d'un backup de la solution")
                self.backup = LNS.creerBackup(constellation)
                printClose()
            printOpen("- retraits de requetes")
            self.destroy.appliquer(LNS,constellation)
            LNS.notifierSansEvenement(constellation)
            printClose()
            printOpen("- remplissage glouton")
            self.remplissage.appliquer(LNS,constellation)
            LNS.notifierSansEvenement(constellation)
            new_objectif = LNS.getObjectif(constellation)
            if not self.accept_worse_solution and new_objectif<old_objectif:
                printOpen("Solution moins bonne. Récupération de l'ancienne solution.")
                self.backup.backupSolution(LNS)
                LNS.notifierBackupSolution(constellation)
                printClose()
            printClose()
            if not self.accept_worse_solution:
                assert(LNS.getObjectif()[0]>=old_objectif[0])    
            return old_objectif<new_objectif
            #printClose()
"""
class Redemarrage(Perturbateur):
    def __init__(self,modeleDeTransition,use_solver=False):
        self.use_solver = use_solver
    def appliquer(self,LNS,constellation):
        old_objectif = LNS.getObjectif(constellation)
        LNS.redemarrer(constellation)
        LNS.resetEtatsRequetes(constellation)
        LNS.greedyFill(constellation,limit_req=self.use_solver,forbid_solver=(not self.use_solver))
        if not self.use_solver:
            LNS.resetEtatsRequetes(constellation)
        return LNS.getObjectif(constellation)>old_objectif
"""        
class DestructeurRequete(Perturbateur):
        # remplissage rapide = remplissage post destruction peu couteuse
        def __init__(self,modeleDeTransition,destroy_rate,remplissage_rapide=False):
            super().__init__("destructeur de requêtes",modeleDeTransition)
            self.remplissage_rapide = remplissage_rapide
            self.destroy_rate = destroy_rate
            id_mpi = MPI.COMM_WORLD.Get_rank()
            seed_option = config.getOptValue("seed")
            self.randomInstance = rd.Random(id_mpi+seed_option)
            
        def detruire(self,LNS,constellation):
            couvertes = LNS.requetesCouvertes()
            detruire = self.randomInstance.sample(couvertes,int(min(len(couvertes)*self.destroy_rate,config.LNS.max_destruction)))
            modes_a_retirer = [(r,m) for (r,m) in LNS.getModesRetenus() if r in detruire and r!=-1]
            self.retirerModes(LNS,constellation,modes_a_retirer)
            # mettre a jour les cca concernées = notifier les changements dans l'objet 'etat requete'
            liste_cca = []
            for (r,m) in modes_a_retirer:
                for (s,a) in constellation.getRequete(r).getCouples():
                    id_cca = LNS.grapheDependances.getActiviteCCA(a)
                    if id_cca not in liste_cca:
                        liste_cca.append(id_cca)
                        LNS.notifierChangementCCA(id_cca)          
            return [x[0] for x in modes_a_retirer],None
        # I. retirer des requetes
        # II. potentiellement re-remplir rapidement
        # III. annuler les modifications sur les etats des requetes :
        #   => si le solver est actif hors destruction il ne faut pas indiquer 
        #   les activites comme infaisable car le solver ici est glouton 
        #   (moins performant)
        def appliquer(self,LNS,constellation):
            old_objectif = LNS.getObjectif(constellation)
            self.detruire(LNS,constellation)
            if self.remplissage_rapide:
                if config.getOptValue("use_solver"):
                    LNS.sauvegarderEtatActivites()
                LNS.greedyFill(constellation,limit_req=config.getOptValue("use_solver"),forbid_solver=True,random=False)
                if config.getOptValue("use_solver"):
                    LNS.restaurerEtatActivites()
            return LNS.getObjectif(constellation)>old_objectif
class VoisinageCP(Operateur):
    def __init__(self,n_coeurs,modeleDeTransition,k=2,iter_op=1):
        super().__init__("Echangeur CP de pairs de CCA",modeleDeTransition,k=k)
        self.n_coeurs = n_coeurs
        #self.positions_cca = None
        self.iter_op = iter_op
        # initialisé à la 1ere utilisation
        self.connectivites_cca = None
        self.requetes_cca = None
        self.cca_observations = None
        self.gaps = []
        self.gains = []
        id_mpi = MPI.COMM_WORLD.Get_rank()
        seed_option = config.getOptValue("seed")
        self.randomModesFixes = rd.Random(id_mpi+seed_option)
    
    def ajouterInfos(self,gap,gain):
        self.gaps.append(gap)
        self.gains.append(gain)
    
    def getGapsInfos(self):
        return self.gaps,self.gains
    
    def incrementerTaillevoisinage(self):
        self.k += 1
        
    def decrementerTaillevoisinage(self):
        self.k -= 1
        assert(self.k >0)
    
    def _extraireModes(self,LNS,constellation,CCAs,requetes_presentes):
        modes = {}
        external_modes = {}
        external_activites = {}
        nouvelles_activites = {}
        for r in requetes_presentes:
                nouvelles_activites[r] = {}
                m = LNS.getModeIfPresent(r)
                if m is None:
                    act = []
                else:
                    act = [x[1] for x in constellation.getRequete(r).getMode(m).getCouples()]
                contenu_actuel = [x for x in act if LNS.grapheDependances.getActiviteCCA(x) not in CCAs]
                allow_external = True # permet d'initialiser la sol avec ce mode
                modes[r],external_modes[r],external_activites[r] = constellation.getRequete(r).genererModesPresentsCCA(LNS.grapheDependances,constellation,contenu_actuel,CCAs,allow_external_modes=allow_external)
                for mode in modes[r]:
                    nouvelles_activites[r][mode.getId()] = [a for a in [x[1] for x in mode.getCouples()] if a not in contenu_actuel]
                if len(nouvelles_activites[r])==0:
                    del nouvelles_activites[r]
        # supprimer les requêtes sans modes
        del_req = []
        for r in modes:
            if len(modes[r])==0:
                del_req.append(r)
        for r in del_req:
            del modes[r]
        if config.getOptValue("verif"):
            for r in nouvelles_activites:
                for m in nouvelles_activites[r]:
                    for a in nouvelles_activites[r][m]:
                        assert(LNS.grapheDependances.getActiviteCCA(a) in CCAs)
        return modes,nouvelles_activites,external_modes,external_activites
    
    def checkCCAVide(self,LNS,nouvelles_activites,CCAs):
         CCAs_obs = []
         for r in nouvelles_activites:
             for m in nouvelles_activites[r]:
                 for a in nouvelles_activites[r][m]:
                     id_cca = LNS.grapheDependances.getActiviteCCA(a)
                     if id_cca not in CCAs_obs:
                         assert(id_cca in CCAs)
                         CCAs_obs.append(id_cca)
                         if len(CCAs)==len(CCAs_obs):
                             assert(sorted(CCAs)==sorted(CCAs_obs))
                             return False
         return len(CCAs)!=len(CCAs_obs)   
     
    # gérer les modes qui ont des activités en dehors des CCA considérées
    def mettreAJourModesExternes(self,LNS,CCAs,constellation,modes_externes,liste_modes_satisfaits):
        for r in modes_externes:
            mode = modes_externes[r]
            activites_externes = [(s,a) for (s,a) in constellation.getRequete(r).getMode(mode).getCouples() if LNS.grapheDependances.getActiviteCCA(a) not in CCAs]
            if constellation.getRequete(r).acceptable(activites_externes,constellation):
                mode = constellation.getRequete(r).ajouterMode(constellation,activites_externes)
                liste_modes_satisfaits.append((r,mode.getId()))
            else:
                for (s,a) in activites_externes:
                    _,cca = LNS.grapheDependances.getActiviteCCA(a)
                    LNS.getSolCCA(s,cca).retirerActivite(constellation,a)
                    
    def evaluerNouvelleSolution(self,LNS,constellation,sol,external_modes,external_activites,CCAs,nouvelles_activites,modes_fixes):
        # retirer les anciens modes
        #print(sorted(list(self.modes_var.keys())))
        nouvelle_liste_modes = [(r,m) for (r,m) in LNS.getModesRetenus() if r not in self.modes_var]
        modes_externes = {r:m for (r,m) in LNS.getModesRetenus() if r in self.modes_var and external_activites[r]}
        # rajouter les nouveaux modes
        # MAJ la liste de modes satisfaits et eliminer les modes problematiques (hors CCA)
        for r in self.modes_var:
            for m in self.modes_var[r]:
                value = sol[self.modes_var[r][m]]
                if value==1:
                    nouvelle_liste_modes.append((r,m))
                    if r in modes_externes:
                        del modes_externes[r]
                        break
        # mettre a jour les modes problématiques (hors CCA)            
        self.mettreAJourModesExternes(LNS,CCAs,constellation,modes_externes,nouvelle_liste_modes)
        return self.mettreAJourSolutionSiMeilleure(LNS,constellation,sol,nouvelle_liste_modes)
        
    def mettreAJourSolutionSiMeilleure(self,LNS,constellation,sol,nouvelle_liste_modes):
        score_avant = LNS.getObjectif()[0]
        copie_modes = deepcopy(LNS.getModesRetenus())
        LNS.setModesRetenus(nouvelle_liste_modes,constellation)
        delta_score = round(LNS.getObjectif()[0]-score_avant,3)
        if not delta_score>=0:
            LNS.setModesRetenus(copie_modes,constellation)
            LNS.verifierSolutionSiVerifMode(constellation) 
            return False,delta_score
        
        self.remplirSolCCAs(constellation,LNS,sol)
        LNS.verifierSolutionSiVerifMode(constellation) 
        return True,delta_score
    
    def recordError(self,LNS,problem,sol,CCAs):                
        if len(os.listdir("error_logs/solver_errors/"))==0:
            id_fichier = 0
        else:
            id_fichier = max([int(file.split('_')[-1]) for file in os.listdir("error_logs/solver_errors")])+1
        path = "error_logs/solver_errors/model_"+problem+"_"+str(CCAs)+"_id_"+str(id_fichier)
        self.model.export_model(path)
        alert("EchangeCP : pas de solution trouvée pour "+str(CCAs)+" : id="+str(id_fichier))
        with open(path,"a") as file:
            for (s,cca) in CCAs:
                file.write(str(LNS.getSolCCA(s,cca))+"\n")
            file.write("modes retenus"+str(sorted(LNS.getModesRetenus()))+'\n')
            file.write("sol is None : "+str(sol is None)+'\n')
            if sol is not None:
                file.write("sol status : "+str(sol.get_solve_status())+'\n')
    
    def verifRequetes(self,external_modes,external_activites,nouvelles_activites,CCAs):
        for r in nouvelles_activites:
            if external_activites[r] and external_modes[r]:
                assert(len(nouvelles_activites)>1) # au moins le mode externe + des combinaisons sur les cca a explorer
            else:
                assert(len(nouvelles_activites)>=1) # au moins un mode
    
    def selectionnerModesFixes(self,LNS,modes,nouvelles_activites,external_modes,external_activites):
        req_in_solution = [(r,m) for (r,m) in LNS.getModesRetenus() if r in nouvelles_activites]
        max_r = config.getOptValue("quota_requests")
        if max_r>=len(req_in_solution):
            return {}
        else: # on retire -1 car le vidage est traité à part
            retenir = self.randomModesFixes.choices(req_in_solution,k=max_r)
            res = {}
            for (r,m) in req_in_solution:
                if (r,m) not in retenir:
                    res[r] = m
                    del modes[r]
                    del external_modes[r]
                    del external_activites[r]
                    del nouvelles_activites[r]
            return res
    
    def choisirCCAEtGenererModes(self,LNS,constellation):
         exist_cca_empty = True
         shiftRightDisplay(2)
         while exist_cca_empty:
             printOpen("Sélection des CCA")
             CCAs,requetes_presentes = self.selectionnerCCAs(LNS,constellation)
             printClose()
             printOpen("Génération des modes")
             modes,nouvelles_activites,external_modes,external_activites = self._extraireModes(LNS,constellation,CCAs,requetes_presentes)       
             modes_fixes = self.selectionnerModesFixes(LNS, modes,nouvelles_activites,external_modes,external_activites)
             exist_cca_empty = self.checkCCAVide(LNS,nouvelles_activites,CCAs)
             printClose()
             if exist_cca_empty:
                 printColor("Il existe une CCA vide. Recherche d'un nouvel ensemble candidat.",c='m')
         shiftLeftDisplay(2)
         return CCAs,requetes_presentes,modes_fixes,modes,nouvelles_activites,external_modes,external_activites
    
    def executerSolver(self,LNS,CCAs):
         tlim = min(LNS.getTempsRestant(),config.getOptValue("CPLim"))
         succes = False
         sol = None
         if tlim>0:
            printOpen("Résolution du modèle CP",tlim,"(s)")
            t1 = time()
            max_try = 2
            essais = 0
            while not succes and essais<max_try:
                try:
                    sol = self.model.solve(TimeLimit = tlim,execfile=config.glob.docplexPath,Workers=self.n_coeurs,log_output=None)
                    succes = True
                except Exception as solver_error:
                    warn(solver_error)
                    alert("Erreur mémoire CpSolver",CCAs)
                    self.recordError(LNS,"memory_pb",None,CCAs)
                    essais += 1
            t2 = time()
            printClose()
         return sol,succes
    
    def afficherInfoGaps(self,LNS,constellation,sol,acceptable,delta_score):
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
            LNS.verifierSolutionSiVerifMode(constellation) 
            self.gaps.append(-1)
            printColor("Aucune meilleure solution trouvée dans le temps imparti",c='y')
                        
    def appliquer(self,LNS,constellation):
        LNS.verifierSolutionSiVerifMode(constellation)
        self.choisirCCAEtGenererModes(LNS,constellation)
        old_objectif = LNS.getObjectif(constellation)
        printOpen("Initialisation de l'opérateur")
        CCAs,requetes_presentes,modes_fixes,modes,nouvelles_activites,external_modes,external_activites = self.choisirCCAEtGenererModes(LNS,constellation)
        printOpen("Construction du modèle CP")
        self.buildModel(LNS,constellation,CCAs,nouvelles_activites,modes,external_modes,external_activites,modes_fixes)
        printClose()
        printClose()
        sol,succes = self.executerSolver(LNS,CCAs)
        if succes:
            if not (sol is not None and sol.get_solve_status()!="Infeasible" and sol.get_solve_status()!="Unknown"):
                self.recordError(LNS,"no-solution",sol,CCAs)
            else:
                LNS.verifierSolutionSiVerifMode(constellation)
                acceptable,delta_score = self.evaluerNouvelleSolution(LNS,constellation,sol,external_modes,external_activites,CCAs,nouvelles_activites,modes_fixes)
                self.afficherInfoGaps(LNS,constellation,sol,acceptable,delta_score)
        return LNS.getObjectif(constellation)>old_objectif
            
    def getStats(self):
        return self.gaps,self.gains
    
    def buildModel(self,LNS,constellation,CCAs,nouvelles_activites,modes,external_modes,external_activites,modes_fixes):
        printOpen("Construction du modèle")
        self.model = CpoModel()
        printClose()
        printOpen("Initialisation des variables")
        self.initVars(LNS,constellation,modes,CCAs,nouvelles_activites,modes_fixes)
        printClose()
        printOpen("Initialisation des contraintes")
        self.initContraintes(LNS,constellation,CCAs,modes,external_modes,external_activites,modes_fixes)
        printClose()
        printOpen("Création d'une solution initiale")
        self.initSol(LNS,constellation,CCAs,modes,modes_fixes)
        printClose()
            
    def initVars(self,LNS,constellation,modes,CCAs,nouvelles_activites,modes_fixes):
        self.modes_var = {}                
        for r in modes:
            self.modes_var[r] = {}
            for m in modes[r]:
                self.modes_var[r][m.getId()] = self.model.binary_var(name="y_"+str((r,m.getId())))

        self.interval_vars = {}
        for r in nouvelles_activites:
            self.interval_vars[r] = {}
            for m in nouvelles_activites[r]:
                for a in nouvelles_activites[r][m]:
                    assert(LNS.grapheDependances.getActiviteCCA(a) in CCAs)
                    self.ajouterVarInterval(constellation,r,a,True)
                    
        for r in modes_fixes:
            self.interval_vars[r] = {}
            m = modes_fixes[r]
            for (s,a) in constellation.getRequete(r).getMode(m).getCouples():
                if LNS.grapheDependances.getActiviteCCA(a) in CCAs:
                    self.ajouterVarInterval(constellation,r,a,False)

    def ajouterVarInterval(self,constellation,r,a,optional):                    
        if a not in self.interval_vars[r]:
            s = constellation.getSatelliteActivite(a)
            start = int(math.ceil(config.glob.getEchelle()*constellation.getSatellite(s).getActivite(a).getDebut()))
            end = int(math.floor(config.glob.getEchelle()*constellation.getSatellite(s).getActivite(a).getFin()))
            duree = int(math.ceil(config.glob.getEchelle()*constellation.getSatellite(s).getActivite(a).getDuree()))
            self.interval_vars[r][a] = self.model.interval_var(start=(start,end),end=(start,end),length=duree,name="Ia_"+str((r,a)))
            if optional:
                self.interval_vars[r][a].set_optional()
                    
    def initNoOverlap(self,constellation,LNS,CCAs):
        self.seq_cca = {}
        for id_cca in CCAs:
            s = id_cca[0]
            printOpen("Récuperer les activites concernées")
            activites = [(a,r,self.interval_vars[r][a]) for r in self.interval_vars for a in self.interval_vars[r] if LNS.grapheDependances.getActiviteCCA(a) == id_cca]
            if len(activites)==0:
                #print("CCAs vides",self.checkCCAVide(LNS,nouvelles_activites,CCAs))
                die("Aucune activités dans la CCA",id_cca)
            
            printClose()
            n = len(activites)
            distance_matrix = np.zeros((n,n),dtype=int)
            vars_list = []
            printOpen("Création de la matrice de distance et séquence",str(n)+"x"+str(n),"activités")
            for i,(a1,r,var) in enumerate(activites):
                vars_list.append(var)
                for j,(a2,r,var2) in enumerate(activites):
                    if i==j:
                        distance_matrix[i][j] = 0
                    else:
                        transition = config.glob.getEchelle()*constellation.getSatellite(s).getTransition(a1,a2,self.modeleDeTransition)
                        distance_matrix[i][j] = int(math.ceil(transition))
            mat = self.model.transition_matrix(distance_matrix)
            self.seq_cca[id_cca] = self.model.sequence_var(vars_list)
            self.model.add(self.model.no_overlap(self.seq_cca[id_cca],distance_matrix=mat))
            printClose()
            
    def initContraintes(self,LNS,constellation,CCAs,modes,external_modes,external_activites,modes_fixes):
        printOpen("Création des contraintes liantes modes - observations")
        self.initContraintesPresenceModes(LNS,modes,external_modes,external_activites,modes_fixes)
        printClose()
        printOpen("Création des contraintes noOverlap")
        self.initNoOverlap(constellation,LNS,CCAs)
        printClose()
        obj1 = self.model.sum([self.modes_var[r][m.getId()]*m.getRecompense() for r in modes for m in modes[r]])
        self.model.maximize(obj1)

    def initContraintesPresenceModes(self,LNS,modes,external_modes,external_activites,modes_fixes):
        for r in self.modes_var:
            if r==-1:
                self.model.add_constraint(self.model.sum([self.modes_var[r][m] for m in self.modes_var[r]])==1)
            else:
                self.model.add_constraint(self.model.sum([self.modes_var[r][m] for m in self.modes_var[r]])<=1)
        # presence mode = présence activites    
        for r in self.interval_vars:
            if r not in modes_fixes:
                for a in self.interval_vars[r]:
                    modes_activite = [mode.getId() for mode in modes[r] if a in [x[1] for x in mode.getCouples()]]
                    self.model.add_constraint(self.model.presence_of(self.interval_vars[r][a])==self.model.sum([self.modes_var[r][m] for m in modes_activite]))
    
    def findDatePlan(self,plan,a):
        for activite,t in plan:
            if activite==a:
               return t
        return None
    
    def findMode(self,constellation,LNS,modes):
        mode_initial = {}
        for r in modes:
           found = False
           for m in modes[r]:
               couples = sorted(m.getCouples())
               for (rr,mm) in LNS.getModesRetenus():
                   if r ==rr:
                       if sorted(constellation.getRequete(rr).getMode(mm).getCouples())==couples:
                           mode_initial[r] = m.getId()
                           found = True
                           break
               if found:
                    break
        return mode_initial
                   
    def initSol(self,LNS,constellation,CCAs,modes,modes_fixes):
        warmstart=self.model.create_empty_solution()
        plan = {}
        for id_cca in CCAs:
            s,cca = id_cca
            plan[id_cca] = LNS.getSolCCA(s,cca).planEarliest(constellation,LNS.getSolCCA(s,cca).getSequence(),self.modeleDeTransition)
        mode_initial = self.findMode(constellation,LNS,modes)
        for r in self.modes_var:
            if r not in modes_fixes:
                for m in self.modes_var[r]:
                    if r in mode_initial and mode_initial[r]==m:
                        warmstart[self.modes_var[r][m]]=1
                    else:
                        warmstart[self.modes_var[r][m]]=0
        
        if config.getOptValue("verif"):
            for r in mode_initial:
                    m = mode_initial[r]
                    for (s,a) in constellation.getRequete(r).getMode(m).getCouples():
                        _,cca = LNS.grapheDependances.getActiviteCCA(a)
                        assert(a in LNS.getSolCCA(s,cca).sequence)
                        
            for s,cca in CCAs:
                for a in LNS.getSolCCA(s,cca).sequence:
                        r = constellation.getRequeteActivite(a)
                        if r in mode_initial:
                            m = mode_initial[r]
                        else:
                            m = modes_fixes[r]
                        obs_mode = [x[1] for x in constellation.getRequete(r).getMode(m).getCouples()]
                        assert(a in obs_mode)
        
        
        # presence des activites
        for r in self.interval_vars:
            for a in self.interval_vars[r]:
                id_cca = LNS.grapheDependances.getActiviteCCA(a)
                t = self.findDatePlan(plan[id_cca], a)
                s = constellation.getSatelliteActivite(a)
                if t is not None:
                    start = int(math.ceil(t*config.glob.getEchelle()))
                    duree = int(math.ceil(constellation.getSatellite(s).getActivite(a).getDuree()*config.glob.getEchelle()))
                    end = start+duree
                    warmstart.add_interval_var_solution(self.interval_vars[r][a],presence=True,start=start,end=end)
                else:
                    warmstart.add_interval_var_solution(self.interval_vars[r][a],presence=False)
        # impossible d'indiquer les modes de départ : les indices des modes présent ne sont pas les mêmes que les nouveaux    
        # indiquer les séquences
        seq = []
        for id_cca in CCAs:
            s,cca = id_cca
            for a in LNS.getSolCCA(s,cca).getSequence():
                assert(id_cca==LNS.grapheDependances.getActiviteCCA(a))
                r = constellation.getRequeteActivite(a)
                if a not in self.interval_vars[r]:
                    print(CCAs,id_cca)
                    print("a dans la séquence, aucune variable d'intervalle créee ?",r,a,LNS.grapheDependances.getActiviteCCA(a),constellation.getRequete(r).getType())
                    print(r in LNS.requetesCouvertes())
                    m = [mm for (rr,mm) in LNS.getModesRetenus() if rr==r][0]
                    print(constellation.getRequete(r).getMode(m))
                    print(r in modes_fixes,r in modes)
                    print("modes de r :")
                    for m in modes[r]:
                        print(m)
                seq.append(self.interval_vars[r][a])
            warmstart[self.seq_cca[id_cca]] = seq
        self.warmstart = warmstart
        self.model.set_starting_point(warmstart)
        
    def remplirSolCCAs(self,constellation,LNS,sol):
        for id_cca in self.seq_cca:
            sequence = []
            for it in sol.get_value(self.seq_cca[id_cca]):
                res = eval(it.get_name().split("_")[1])
                r,a = res[0],res[1]
                sequence.append(a)
            LNS.ecraserSolutionCCA(constellation,id_cca,sequence)
            s,cca = id_cca
            if config.verifMode():
                for a in sequence:
                        find = None
                        r = constellation.getRequeteActivite(a)
                        for m in self.modes_var[r]:
                            if sol[self.modes_var[r][m]]==1:
                                if find is not None:
                                    raise ValueError("plusieurs modes pour",r)
                                find = m
                        assert(find is not None)
                        assert(sol[self.modes_var[r][find]]==1)
                        assert((r,find) in LNS.getModesRetenus())
                        s = constellation.getSatelliteActivite(a)
   
                
        if config.getOptValue("verif"):
            for id_cca in self.seq_cca:
                s,cca = id_cca
                for a in LNS.getSolCCA(s,cca).getSequence():
                    r = constellation.getRequeteActivite(a)
                    find = None
                    for m in self.modes_var[r]:
                        if sol[self.modes_var[r][m]]==1:
                            find = m
                            break
                    assert((r,find) in LNS.getModesRetenus())
                    s = constellation.getSatelliteActivite(a)
                    assert((s,a) in constellation.getRequete(r).getMode(find).getCouples())
            LNS.verifierSolution(constellation)
            
class DestroyAndRepairGreedyCCA(Operateur):
    def __init__(self,LNS,modeleDeTransition,accept_worse_solution,k=2):
        super().__init__("Greedy windows sur des CCA",modeleDeTransition)
        self.accept_worse_solution = accept_worse_solution
        self.destroy = DestructeurEnsembleCCA(LNS.positions_cca,modeleDeTransition,k=k)
        self.remplissage = RemplissageGloutonFenetre(modeleDeTransition)        
        
    def appliquer(self,LNS,constellation):
        #printOpen("Perturbation :",c='y')
        old_objectif = LNS.getObjectif(constellation)
        if not self.accept_worse_solution:
            printOpen("Création d'un backup de la solution")
            self.backup = LNS.creerBackup(constellation)
            printClose()
        printOpen("- retraits de requetes")
        requetes_presentes,nouveaux_modes,CCAs = self.destroy.detruire(LNS,constellation)
        LNS.notifierSansEvenement(constellation)
        printClose()
        printOpen("- remplissage glouton")
        self.remplissage.appliquer(LNS,constellation,CCAs,requetes_presentes)
            
        new_objectif = LNS.getObjectif(constellation)
        if not self.accept_worse_solution and new_objectif<old_objectif:
            printOpen("Solution moins bonne. Récupération de l'ancienne solution.")
            self.backup.backupSolution(LNS)
            LNS.notifierBackupSolution(constellation)
            printClose()
        printClose()
        return old_objectif<new_objectif

            
class DestructeurEnsembleCCA(OperateurDestruction):
    def __init__(self,positions_cca,modeleDeTransition,k=2):
        super().__init__("destructeur de CCA",modeleDeTransition,k=k)
        self.positions_cca = positions_cca
       
    def detruire(self, LNS, constellation):
        CCAs,requetes_presentes = self.selectionnerCCAs(LNS,constellation)
        sat_lste = [id_cca[0] for id_cca in CCAs]
        nouveau_contenu_mode = {}
        modes_courants = {}
        LNS.verifierSolution(constellation)
        for s,cca in CCAs:
            for a in LNS.getSolCCA(s,cca).getSequence():
                r = constellation.getRequeteActivite(a)
                if r not in nouveau_contenu_mode:
                    modes_couverts = [m for (rr,m) in LNS.getModesRetenus() if rr==r]
                    assert(len(modes_couverts)==1)
                    modes_courants[r] = modes_couverts[0]
                    couples = constellation.getRequete(r).getMode(modes_couverts[0]).getCouples()
                    nouveau_contenu_mode[r] = [(sat,act) for (sat,act) in couples if act!=a]
                else:
                    nouveau_contenu_mode[r] = [(sat,act) for (sat,act) in nouveau_contenu_mode[r] if act!=a]
        for s,cca in CCAs:    
            LNS.getSolCCA(s,cca).reset()
        for r in LNS.etatRequetes:
            for id_cca in CCAs:
                LNS.etatRequetes[r].resetCCA(id_cca)
        for r in nouveau_contenu_mode:
            m = modes_courants[r] 
            if r==-1:
                for (s,a) in constellation.getRequete(r).getMode(m).getCouples():
                    id_cca = LNS.grapheDependances.getActiviteCCA(a)
                    if id_cca in CCAs:
                        LNS.solCCAs[s][id_cca[[1]]].insererActivitePlansCritiques(constellation,a)
                LNS.ajouterModeRetenu((-1,0))
            else:
                LNS.retirerModeRetenu((r,m))
                if constellation.getRequete(r).acceptable(nouveau_contenu_mode[r],constellation):
                    mode = constellation.getRequete(r).ajouterMode(constellation,nouveau_contenu_mode[r])
                    LNS.ajouterModeRetenu((r,mode.getId()))
        if -1 in nouveau_contenu_mode:
            del nouveau_contenu_mode[-1]
        return requetes_presentes,list(nouveau_contenu_mode.keys()),CCAs

        
class RemplissageGloutonFenetre(OperateurReparation):
    def __init__(self,modeleDeTransition,k=2):      
        super("Remplissage glouton fenêtre",modeleDeTransition)
        self.k = k # choix d'une cca1 et d'une cca2 parmi les k plus proches
        # initialisé à la 1ere utilisation
        self.requetes_cca = None
        self.gains = []
        id_mpi = MPI.COMM_WORLD.Get_rank()
        seed_option = config.getOptValue("seed")
        self.randomModesFixes = rd.Random(id_mpi+seed_option)
    
    def ajouterInfos(self,gain):
        self.gains.append(gain)

    def _extraireModes(self,LNS,constellation,CCAs,requetes_presentes):
        modes = {}
        external_modes = {}
        external_activites = {}
        nouvelles_activites = {}
        for r in requetes_presentes:
                nouvelles_activites[r] = {}
                m = LNS.getModeIfPresent(r)
                if m is None:
                    act = []
                else:
                    act = [x[1] for x in constellation.getRequete(r).getMode(m).getCouples()]
                contenu_actuel = [x for x in act if LNS.grapheDependances.getActiviteCCA(x) not in CCAs] # activités externes
                allow_external = True # permet d'initialiser la sol avec ce mode
                modes[r],external_modes[r],external_activites[r] = constellation.getRequete(r).genererModesPresentsCCA(LNS.grapheDependances,constellation,contenu_actuel,CCAs,allow_external_modes=allow_external)
                for mode in modes[r]:
                    nouvelles_activites[r][mode.getId()] = [a for a in [x[1] for x in mode.getCouples()] if a not in contenu_actuel]
                if len(nouvelles_activites[r])==0:
                    del nouvelles_activites[r]
        # supprimer les requêtes sans modes
        del_req = []
        for r in modes:
            if len(modes[r])==0:
                del_req.append(r)
        for r in del_req:
            del modes[r]
        if config.getOptValue("verif"):
            for r in nouvelles_activites:
                for m in nouvelles_activites[r]:
                    for a in nouvelles_activites[r][m]:
                        assert(LNS.grapheDependances.getActiviteCCA(a) in CCAs)
        return modes,nouvelles_activites,external_modes,external_activites
    
    def evaluerNouvelleSolution(self,LNS,constellation,sol,external_modes,external_activites,CCAs,nouvelles_activites,modes_fixes):
        # retirer les anciens modes
        #print(sorted(list(self.modes_var.keys())))
        nouvelle_liste_modes = [(r,m) for (r,m) in LNS.getModesRetenus() if r not in self.modes_var]
        retraits = {r:m for (r,m) in LNS.getModesRetenus() if r in self.modes_var and self.hasProblematicExternActivities(external_activites,external_modes,r)}
        # rajouter les nouveaux modes
        for r in self.modes_var:
            for m in self.modes_var[r]:
                value = sol[self.modes_var[r][m]]
                if value==1:
                    nouvelle_liste_modes.append((r,m))
                    if r in retraits:
                        del retraits[r]
                        break
        score_avant = LNS.getObjectif()[0]
        copie_modes = deepcopy(LNS.getModesRetenus())
        LNS.setModesRetenus(nouvelle_liste_modes,constellation)
        delta_score = round(LNS.calculerObjectif(constellation)[0]-score_avant,3)
        if not delta_score>=0:
            LNS.setModesRetenus(copie_modes,constellation)
            return False,delta_score
        
        for r in retraits:
            mode = retraits[r]
            for (s,a) in constellation.getRequete(r).getMode(mode).getCouples():
                id_cca = LNS.grapheDependances.getActiviteCCA(a)
                LNS.solCCAs[s][id_cca].retirerActivite(constellation,a)
        return True,delta_score
    

    def verifRequetes(self,external_modes,external_activites,nouvelles_activites,CCAs):
        for r in nouvelles_activites:
            if external_activites[r] and external_modes[r]:
                assert(len(nouvelles_activites)>1) # au moins le mode externe + des combinaisons sur les cca a explorer
            else:
                assert(len(nouvelles_activites)>=1) # au moins un mode

    def selectionnerMode(self,LNS,constellation,modes):
        score_courant = {}
        for (r,m) in LNS.getModesRetenus():
            if r in modes:
                score_courant[r] = constellation.getRequete(r).getMode(m).getRecompense()
        cle = lambda rm:rm[1].getRecompense()-score_courant.get(r,0)
        (req,mode) = max([(r,m) for r in modes for m in modes[r]],key=cle)
        modes[req].remove(mode)
        if len(modes[req])==0:
            del modes[req]
        return req,mode

    def validerMode(self,modes,r):
        if r in modes:
            del modes[r]

    def appliquer(self,LNS,constellation,CCAs,requetes_presentes):
         old_objectif = LNS.getObjectif(constellation)
         shiftRightDisplay(2)
         printOpen("Génération des modes")
         modes,nouvelles_activites,external_modes,external_activites = self._extraireModes(LNS,constellation,CCAs,requetes_presentes)       
         printClose()
         printOpen("Insertion des modes")
         if -1 in modes:
             del modes[-1]
         while len(modes.keys())>0 and LNS.getTempsEcoule() < config.getOptValue("time"):
             r,mode_candidat = self.selectionnerMode(LNS,constellation,modes) 
             printOpen("Tentative d'insertion de la requête",r)
             mode,sol_cca,vide = LNS.etatRequetes[r].tenterInsererMode(LNS,constellation,r,mode_candidat,False)
             #printOpen("Insérer la requete",constellation.getRequete(r).getType(),r,c='g')
             if mode is not None:
                 #printColor("valider",noeuds_a_valider,c='g')
                 self.validerMode(modes,r)
                 LNS.validerSol(sol_cca)
                 LNS.setModesRetenus([(rr,m) for (rr,m) in LNS.getModesRetenus() if r!=r],constellation)
                 LNS.ajouterModeRetenu((r,mode.getId()))
                 if config.getOptValue("verif"):
                     LNS.verifierCCA(constellation)
                 printClose("Succès",c='g')
             else:
                printColor("Mode infaisable",c='m')
                printClose("Echec",c='r')
         printClose()
         shiftLeftDisplay(2)
         new_obj = LNS.getObjectif(constellation)
         delta_score = new_obj - old_objectif
         if config.getOptValue("verif"):
             LNS.verifierSolution(constellation)
         self.gains.append(delta_score)
         return new_obj>old_objectif

    def findDatePlan(self,plan,a):
        for activite,t in plan:
            if activite==a:
               return t
        return None
    
    def findMode(self,constellation,LNS,modes):
        mode_initial = {}
        for r in modes:
           found = False
           for m in modes[r]:
               couples = sorted(m.getCouples())
               for (rr,mm) in LNS.getModesRetenus():
                   if r ==rr:
                       if sorted(constellation.getRequete(rr).getMode(mm).getCouples())==couples:
                           mode_initial[r] = m.getId()
                           found = True
                           break
               if found:
                    break
        return mode_initial
from ..Utils.config import *
global config
config = Config()
instance = config.instance
from mpi4py import MPI
if config.getOptValue("help"):
    if MPI.COMM_WORLD.Get_rank()==0:
        config.afficherAide()
else:
    
    from . import operateursLNS as op
    from .LNSTools import *
    
    from ..model.solution import *
    from ..model.constellation import *
    from ..model.composantes import *
    from ..model.solution_composantes import *
    from ..model.oracle import Oracle
    
    from ..Utils.Utils import *
    from ..Utils.Utils import printColor,warn,alert,printOpen,printClose,printMaster,choseAndBroadcastFile
    from ..Utils.Utils import getDisplayDepth,shiftLeftDisplay,shiftRightDisplay
    
    from time import time,sleep
    import math
    import numpy as np
    from copy import deepcopy,copy
    
    class LNS(Solver):
        
        class HeuristiqueActivites:
            def __init__(self,nom):
                self.nom = nom
            
            def choisirActivite(self,constellation,r,activites_a_inserer):
                raise ValueError('classe abstraite')
        
        class HeuristiqueRatio(HeuristiqueActivites):
            def __init__(self,nom,alpha):
                super().__init__(nom)
                self.alpha = alpha
                
            def scoreEnsembleActivites(self,LNSsolver,constellation,r,activites_a_inserer,cheap=False):        
                #printOpen("Scoring des activités",c='b')
                timeSlots = {a : constellation.getRequete(r).getTimeSlot(a) for a in activites_a_inserer}
                listeObs = [constellation.getRequete(r).findActivite(a) for a in activites_a_inserer]
                recompense,score_temp = constellation.getRequete(r).scoreListeObservations(listeObs,constellation,timeSlots)
                duree = sum([constellation.getSatellite(constellation.getSatelliteActivite(a)).getActivite(a).getDuree() for a in activites_a_inserer])
                if not cheap and self.alpha!=0:
                    trans = LNSsolver.etatRequetes[r].diffCoutTransition(LNSsolver,constellation,activites_a_inserer)
                else:
                    trans = 0
                """
                if not cheap:
                    retard = LNSsolver.etatRequetes[r].retard(constellation,activites_a_inserer)
                else:
                    retard = 0
                """
                res = (recompense,score_temp,duree+trans)
                return self.aggregerNoeud(res)
            
            def aggregerNoeud(self,lambda_v):
                return lambda_v[0]/(1+lambda_v[2])**self.alpha
                    
        class EtatRequete:
            def __init__(self,LNS,constellation,r):
                self.explications = []
                self.inactifs = []
                self.r = r
                #self.score_courant = {}
                #self.mode_courant
                id_mpi = MPI.COMM_WORLD.Get_rank()
                seed_option = config.getOptValue("seed")
                self.randominitInsertTmp = rd.Random(id_mpi+seed_option)
                self.initRechercheMode(LNS,constellation)
            def diffCoutTransition(self,LNS,constellation,activites_a_inserer):
                cout = 0
                for id_cca in self.tmp[activites_a_inserer]:
                    s = id_cca[0]
                    cout += (self.tmp[activites_a_inserer][id_cca].coutTransitionComposante(constellation) - LNS.getSolCCA(s,cca).coutTransitionComposante(constellation))
                return cout
            def getIdRequete(self):
                return self.r            
            def initRechercheMode(self,LNS,constellation):
                self.tmp = {} # les solutions temporaires des cca
                r = self.r
                self.inactifs = []
                self.activites_faisables = {}
                # trouver le meilleur mode privé des activités inactives
                mode = constellation.getRequete(r).getBestModeWithoutInactives(constellation,self.inactifs)
                if mode is not None:
                    self.score_courant = mode.getRecompense()
                else:
                    self.score_courant = 0
                return self.inactifs
            def getScoreCourant(self):
                return self.score_courant
            
            def getActiviteByCCA(self,activites,LNS):
                act_by_cca = {}
                for (s,a) in activites:
                    id_cca = LNS.grapheDependances.getActiviteCCA(a)
                    if id_cca not in act_by_cca:
                        act_by_cca[id_cca] = []
                    act_by_cca[id_cca].append(a)
                for id_cca in act_by_cca:
                    act_by_cca[id_cca] = tuple(act_by_cca[id_cca])
                return act_by_cca
                
            def tenterInsererRequete(self,LNS,constellation,r,forbid_solver,modeleDeTransition):
                #shiftRightDisplay(1)
                self.solution = {}
                #print(r,sorted(self.inactifs))
                mode = constellation.getRequete(r).getBestModeWithoutInactives(constellation,self.inactifs)
                if mode is None:
                    printColor("Requete "+str(self.r) +" : aucun mode candidat trouvé",c='r')
                    return None,None,True  
                self.score_courant = mode.getRecompense()
                return self.tenterInsererMode(LNS,constellation,r,mode,forbid_solver,modeleDeTransition)
            
            def sauvegarderEtatActivites(self):
                self.copie_inactifs = copy(self.inactifs)
                self.copie_activites_faisables = copy(self.activites_faisables)
            
            def restaurerEtatActivites(self):
                self.inactifs = self.copie_inactifs
                supp = []
                for a in self.activites_faisables:
                    if not self.activites_faisables[a]:
                        if a not in self.copie_activites_faisables or self.copie_activites_faisables[a]:
                            supp.append(a)
                for a in supp:
                    del self.activites_faisables[a]
                    del self.tmp[a]
            
            def tenterInsererMode(self,LNS,constellation,r,mode,forbid_solver,modeleDeTransition):
                #shiftRightDisplay(1)
                self.solution = {}
                if mode is None:
                    printColor("Requete "+str(self.r) +" : aucun chemin trouvé",c='r')
                    return None,None,True  
                self.score_courant = mode.getRecompense()
                realisable = True
                #shiftRightDisplay(2)
                printOpen("Test de faisabilité du chemin (+ évaluation précise des scores)")
                act_cca = self.getActiviteByCCA(mode.getCouples(),LNS)
                for id_cca in act_cca:
                    faisable = self.testFaisabilite(constellation,LNS,r,act_cca[id_cca],id_cca,forbid_solver,modeleDeTransition)
                    if not faisable:
                        realisable = False
                        self.inactifs+=list(act_cca[id_cca])
                        self.activites_faisables[act_cca[id_cca]] = False
                        printColor("Requete "+str(self.r) +" : mode infaisable",c='r')
                        printClose()
                        return None,None,False
                printClose()
                assert(realisable)
                if config.verifMode():
                    for id_cca in act_cca:
                        s,cca = id_cca
                        assert(LNS.getSolCCA(s,cca).sequenceFaisable(constellation,modeleDeTransition))
                printColor("Requete "+str(self.r) +" mode faisable : ",mode,c='g')
                mode = constellation.getRequete(r).validerModeCandidat(constellation)
                """
                doublons = self.doublonsCCAs(noeuds)
                if len(doublons)>0:
                    warn("Plusieurs noeuds sur la même CCA : ",doublons)
                """
                for id_cca in act_cca:
                    CCAS = list(self.tmp[act_cca[id_cca]].keys())
                    assert(len(CCAS)==1)
                    if id_cca in self.solution:
                        die("CCA déjà présente",id_cca)
                    self.solution[id_cca] = self.tmp[act_cca[id_cca]][id_cca]
                return mode,self.solution,False    
    
            def doublonsCCAs(self,noeuds):
                doublons = []
                CCAs = [list(self.tmp[x].keys()) for x in noeuds]
                for cca_liste in CCAs:
                    if not(len(cca_liste)==1):
                        print(cca_liste,CCAs)
                for i,cca1 in enumerate(CCAs):
                    for j,cca2 in enumerate(CCAs):
                        if i!=j and cca1==cca2 and cca1[0] not in doublons:
                            doublons.append(cca1[0])
                return doublons
            
            def testFaisabilite(self,constellation,LNS,r,X,id_cca,forbid_solver,modeleDeTransition):
                if X not in self.activites_faisables:
                        assert(X not in self.tmp)
                        self.activites_faisables[X] = self.insertionTemporaire(LNS,constellation,r,X,forbid_solver,modeleDeTransition)
                        assert(X in self.tmp)
                        for a in self.tmp[X][id_cca].sequence:
                            if not LNS.grapheDependances.getActiviteCCA(a)==id_cca:
                                print(X,id_cca,self.tmp[X][id_cca].sequence,a,"coeur",MPI.COMM_WORLD.Get_rank())
                            assert(LNS.grapheDependances.getActiviteCCA(a)==id_cca)
                return self.activites_faisables[X]
            def checkActivitesCCA(self,LNS,activites_a_inserer):
                cca_activites = {}
                for a in activites_a_inserer:
                    id_cca = LNS.grapheDependances.getActiviteCCA(a)
                    if id_cca not in cca_activites:
                        cca_activites[id_cca] = []
                    cca_activites[id_cca].append(a)
                id_cca = list(cca_activites.keys())[0]
                assert(len(cca_activites)<=1)
                return cca_activites,id_cca
            def initInsertionTemporaire(self,LNS,constellation,r,activites_a_inserer):
                self.tmp[activites_a_inserer] = {}
                cca_activites,id_cca = self.checkActivitesCCA(LNS,activites_a_inserer)
                # plan earliest,latest
                s,cca = id_cca
                groups = {0 : copy(LNS.getSolCCA(s,cca).getSequence())}
                #if(not LNS.modeleDeTransition.estTimeDependent()):
                try:
                    LNS.plansCritiques(constellation,id_cca)
                except NumericalError:
                    assert(LNS.modeleDeTransition.estTimeDependent())
                self.tmp[activites_a_inserer][id_cca] = LNS.getSolCCA(s,cca).copieTemporaire()
                longueur_avant_insertion = len(self.tmp[activites_a_inserer][id_cca].sequence)
                cp = cca_activites[id_cca].copy()
                self.randominitInsertTmp.shuffle(cp)
                cca_activites[id_cca] = cp
                for a in self.tmp[activites_a_inserer][id_cca].sequence:
                    assert(id_cca==LNS.grapheDependances.getActiviteCCA(a))
                return longueur_avant_insertion,groups,s,id_cca,cca_activites
        
            def insertionTemporaireGloutonne(self,constellation,cca_activites,activites_a_inserer,id_cca,solver_used_after,modeleDeTransition):
                if config.verifMode():
                    assert(self.tmp[activites_a_inserer][id_cca].sequenceFaisable(constellation,modeleDeTransition))
                faisable = True
                printOpen("insertion des activités",activites_a_inserer," (méthode gloutonne) dans la cca",id_cca,c='c')
                for i,p in enumerate(cca_activites[id_cca]):
                    printOpen("insertion gloutonne")
                    if faisable:
                        assert(self.tmp[activites_a_inserer][id_cca].sequenceFaisable(constellation,modeleDeTransition))
                    #if not modeleDeTransition.estTimeDependent() and faisable:
                        if self.tmp[activites_a_inserer][id_cca].plansAJour():
                            inserer = self.tmp[activites_a_inserer][id_cca].insererActivitePlansCritiques(constellation,p,modeleDeTransition)
                        else:
                            assert(modeleDeTransition.estTimeDependent())
                            inserer = self.tmp[activites_a_inserer][id_cca].insererActivite(constellation,p,modeleDeTransition)
                    else:
                        inserer = self.tmp[activites_a_inserer][id_cca].insererActivite(constellation,p,modeleDeTransition)
                    if inserer:
                        assert(self.tmp[activites_a_inserer][id_cca].sequenceFaisable(constellation,modeleDeTransition))
                    
                    faisable = inserer and faisable
                    printClose()
                    if not faisable and not solver_used_after:
                        break
                    printOpen("MAJ des plans critiques")
                    if faisable:
                        try:
                        #if not modeleDeTransition.estTimeDependent() and faisable:
                            self.calculerPlansCritiques(activites_a_inserer,id_cca,constellation,modeleDeTransition)
                        except NumericalError:
                            assert(modeleDeTransition.estTimeDependent())
                            
                    printClose()
                printClose()
                return faisable        
            def calculerPlansCritiques(self,activites_a_inserer,id_cca,constellation,modeleDeTransition):    
                if config.getOptValue("use_solver"):
                    self.tmp[activites_a_inserer][id_cca].MAJPlansCritiques(constellation,modeleDeTransition)
                else: # a remplacer par une propagation
                    self.tmp[activites_a_inserer][id_cca].MAJPlansCritiques(constellation,modeleDeTransition)
            def resetCCA(self,id_cca,grapheDep):
                #printColor("reset cca",cca,"; requête",self.r,c='y')
                del_x = []
                for X in self.tmp:
                    if id_cca in self.tmp[X]:
                        del_x.append(X)
                for x in del_x:
                    del self.tmp[x]
                    if x in self.activites_faisables:
                        del self.activites_faisables[x]
                
                activites_a_liberer = grapheDep.getActivitesComposante(id_cca)
                for a in self.inactifs:
                    if a in activites_a_liberer:
                        self.inactifs.remove(a)
                return len(del_x)
                
            def insertionTemporaireSolver(self,constellation,LNS,cca_activites,faisable,longueur_avant_insertion,activites_a_inserer,id_cca,groups,forbid_solver,modeleDeTransition):
                try:
                    copy_sequence = self.tmp[activites_a_inserer][id_cca].getSequence()
                    assert(len(activites_a_inserer)>0)
                    printOpen("recherche locale",c='c')
                    groups[1] = cca_activites[id_cca]
                    assert(len(groups[1])>0)
                    for a in self.tmp[activites_a_inserer][id_cca].sequence:
                        assert(LNS.grapheDependances.getActiviteCCA((a)==id_cca))
                    self.tmp[activites_a_inserer][id_cca].rechercheLocale(constellation,'OPTWGroups',modeleDeTransition,groups=groups,allow_no_solution=True)
                    for a in self.tmp[activites_a_inserer][id_cca].sequence:
                        assert(LNS.grapheDependances.getActiviteCCA((a)==id_cca))
                    longueur_supposee = longueur_avant_insertion + len(activites_a_inserer)
                    faisable = (longueur_supposee == len(self.tmp[activites_a_inserer][id_cca].sequence))
                    printClose()
                    return faisable
                except InfeasibleSolutionfromSolverException:
                    if config.verifMode():
                        print("Erreur solver : solution infaisable.")
                    printClose()
                    return False
        
            # insere les activites dans tmp (copies des sol CCA) et renvoie la liste des cca modifiées
            # sur un noeud : une seule cca est modifiee
            def insertionTemporaire(self,LNS,constellation,r,activites_a_inserer,forbid_solver,modeleDeTransition):
                longueur_avant_insertion,groups,s,id_cca,cca_activites = self.initInsertionTemporaire(LNS,constellation,r,activites_a_inserer)
                solver_used = config.getOptValue("use_solver") and not forbid_solver
                faisable = False
                #if not solver_used:
                faisable = self.insertionTemporaireGloutonne(constellation,cca_activites,activites_a_inserer, id_cca, solver_used,modeleDeTransition)
                assert(faisable is not None)
                if not faisable and solver_used:
                    faisable = self.insertionTemporaireSolver(constellation,LNS,cca_activites, faisable, longueur_avant_insertion, activites_a_inserer, id_cca, groups,forbid_solver,modeleDeTransition)
                return faisable
        
        class SolutionSaver:
            def __init__(self,LNS,constellation):
                self.solCCAs = {s : {cca : LNS.getSolCCA(s,cca).copieTemporaire() for cca in LNS.getSolCCAs()[s]} for s in LNS.getSolCCAs()}
                self.modes_retenus = LNS.getModesRetenus().copy()
                self.etatRequetes = deepcopy(LNS.etatRequetes)
                self.objectif = LNS.getObjectif(constellation,recompute=True)
            
            def backupSolution(self,LNS):
                LNS.setAllSolCCAs(self.solCCAs)
                LNS.ecraserModes(self.modes_retenus)
                LNS.setObjectif(self.objectif)
                LNS.etatRequetes = self.etatRequetes
        
        def __init__(self,constellation,start_date,modeleDeTransition,dt_construction_transition,tlim=np.inf,CCAs=None,solution=None):
            super().__init__(constellation,start_date,modeleDeTransition,dt_construction_transition,tlim=tlim,CCAs=CCAs,solution=solution)
            shift = getDisplayDepth()-1
            shiftLeftDisplay(shift)
            #self.explication = {}
            # Initialisation des requetes
            self.initRequetes(constellation)
            # composantes et solutions de composantes
            self.oracle = {}
            self.positions_cca = {}
            self.cca_obs = []
            self.nDegradation = 0
            self.nNonDegrade = 0
            self.tlim = min(self.tlim,config.getOptValue("time"))
            # cca qui contiennent au moins une obs (écarter les cca avec qu'un téléchargement)
            for id_cca in self.grapheDependances.getComposantes():
                s,cca = id_cca
                self.positions_cca[id_cca] = np.mean(np.array([constellation.getSatellite(s).getActivite(a).getCoordonnees() for a in self.grapheDependances.getActivitesComposante(id_cca)]))
                for a in self.grapheDependances.getActivitesComposante(id_cca):
                    if constellation.getSatellite(constellation.getSatelliteActivite(a)).estObservation(a):
                        self.cca_obs.append(id_cca)
                        break
                #self.oracle[id_cca] = Oracle(id_cca)
            self.calculerRequetesPresentes(constellation)
            # si 1 coeur : slave local
            comm = MPI.COMM_WORLD
            comm.Barrier()
            self.notifierPreprocessing(constellation)
            self.initOperateurs()
            self.initSourcesAleatoires()
            shiftRightDisplay(shift)

        def initSourcesAleatoires(self):
            id_mpi = MPI.COMM_WORLD.Get_rank()
            seed_option = config.getOptValue("seed")
            self.randomChoixRequete = rd.Random(id_mpi+seed_option)
            
        # retourne les CCA qui contiennent au moins une observation
        def getCCAObservations(self):
            return self.cca_obs
        
        def distanceCCA(self,cca1,cca2):
            return math.dist(self.positions_cca[cca1]-self.positions_cca[cca2])
        
        """
            =============================================== 
                            RESOLUTION
            =============================================== 
        """
        def initRequetes(self,constellation):
            super().initRequetes(constellation,initFirstMode=False)           
    
        def insererSequences(self,constellation,sequences):
            requetes_retenues = [x[0] for x in self.getModesRetenus()]
            for id_cca in sequences:
                s,cca = id_cca
                seq = [a for a in sequences[id_cca] if constellation.getRequeteActivite(a) in requetes_retenues]
                self.getSolCCA(s,cca).setSequence(constellation,seq)
    
        def creerBackup(self,constellation):
            return self.SolutionSaver(self,constellation)
                
        def LNSsolve(self,constellation,mailbox,dernier_record,it):
            filtre = ['time','obj','modes','requetes','best']
            op_call = 0
            stableIt = 0
            while stableIt<config.getOptValue("stableIt") and op_call<config.getOptValue("max_operateur") and time()-self.start_date<self.tlim:
                if config.getOptValue("verif"):
                    self.verifierSolution(constellation)
                if config.getOptValue("dynamic"):
                    change,liste_requetes = constellation.libererNouvellesRequetes(self.grapheDependances)
                    if change:
                        self.MAJNouvellesRequetes(constellation,time()-self.start_date,liste_requetes)
                        self.operateur.MAJNouvellesRequetes(self)
                        #self.perturbateur.MAJNouvellesRequetes(self)
                improvment = self.appliquerOperateur(constellation)
                self.notifierFinOperateur(constellation)
                #self.setObjectif(self.calculerObjectif(constellation))
                if not improvment:
                    stableIt += 1
                else:
                    stableIt = 0
                op_call += 1
            self.notifierFinIteration(constellation)
            if time()-dernier_record>=config.glob.periode_affichage_info:
                dernier_record = time()
                self.afficherInfo(time(),self.start_date,constellation,title="RESULTAT ITERATION " +str(it)+" CPU "+str(MPI.COMM_WORLD.Get_rank()),filtre=filtre,color='c',add={"nombre d'appels à l'opérateur durant l'itération'":op_call})
            return op_call,dernier_record
        
        def initOperateurs(self):
            #self.heuristiqueGreedy = self.HeuristiqueRatio('ratio',1)
            self.operateurs = []
            # choix de la méthode de perturbation
            if config.getOptValue("version") == "greedy-request":
                #self.perturbateur = op.Redemarrage()
                perturb = config.getOptValue("perturb_rate")
                destroy = config.getOptValue("destroy")
                accept_worse_solution=True
                self.perturbateur = op.DestroyAndRepairGreedyRequest(self.modeleDeTransition,perturb,accept_worse_solution)
                accept_worse_solution=False
                self.operateur = op.DestroyAndRepairGreedyRequest(self.modeleDeTransition,destroy,accept_worse_solution)
            elif config.getOptValue("version") == "greedy-cca":
                #self.perturbateur = op.Redemarrage()
                perturb = config.getOptValue("perturb_rate")
                accept_worse_solution=True
                self.perturbateur = op.DestroyAndRepairGreedyRequest(self.modeleDeTransition,perturb,accept_worse_solution)
                accept_worse_solution=False
                self.operateur = op.DestroyAndRepairGreedyCCA(self,modeleDeTransition,accept_worse_solution,k=config.getOptValue("n_cca"))
            elif config.getOptValue("version") in ["hybrid","coop"]:
                perturb = config.getOptValue("perturb_rate")
                accept_worse_solution=True
                self.perturbateur = op.DestroyAndRepairGreedyRequest(self.modeleDeTransition,perturb,accept_worse_solution)
                self.operateur = op.VoisinageCP(MPI.COMM_WORLD.Get_size(),self.modeleDeTransition,k=config.getOptValue("n_cca"))
            else:
                raise NameError("Methode de perturbation inconnue")
                
        def notifierChangementCCA(self,id_cca):
            cca_tmp_supp = 0
            for r in self.etatRequetes:
                cca_tmp_supp += self.etatRequetes[r].resetCCA(id_cca,self.grapheDependances)
            return cca_tmp_supp
        
        def checkSequence(self,constellation):
            for s in self.getSolCCAs():
                for cca in self.getSolCCAs():        
                    for a in self.getSolCCA(s,cca).getSequence() :
                        assert(s==constellation.getSatelliteActivite(a))
        
        def ecraserSolutionCCA(self,constellation,id_cca,sequence):
            s,cca = id_cca
            self.getSolCCA(s,cca).setSequence(constellation,sequence,self.modeleDeTransition)
            return self.notifierChangementCCA(id_cca)
            
        def creerSolutionInitiale(self,constellation):
            repetitions = 1 # determiste désormais donc inutile de chercher plus loin
            self.greedyFill(constellation,limit_req=False,forbid_solver=True)
            if config.getOptValue("use_solver"):
                self.resetEtatsRequetes(constellation)
            # quand on passe d'un greedy sans solver à greedy avec solver il faut reset les inactifs
            return repetitions

        def initEtatRequetes(self,constellation):
            self.etatRequetes = {}
            for r in constellation.getToutesRequetes():
                self.etatRequetes[r] = self.EtatRequete(self,constellation,r)
            
        def resetEtatsRequetes(self,constellation):
            for r in self.etatRequetes:
                self.etatRequetes[r].initRechercheMode(self,constellation)
        
        def sauvegarderEtatActivites(self):
            for r in self.etatRequetes:
                self.etatRequetes[r].sauvegarderEtatActivites()
                
        def restaurerEtatActivites(self):
            for r in self.etatRequetes:
                self.etatRequetes[r].restaurerEtatActivites()
                
        def resoudre(self,constellation,mailbox,afficher=True):
            comm = MPI.COMM_WORLD
            size = comm.Get_size()
            #objectif = self.calculerObjectif(constellation)
            #filtre = ['time','obj','modes','requetes','best'] # infos à afficher
            it=0
            dernier_record = time()
            application_operateur = 0
            n_perturbations = 0
            self.initEtatRequetes(constellation)
            
            self.creerSolutionInitiale(constellation)
            #self.verifierSolutionSiVerifMode(constellation)
            #self.setObjectif(self.calculerObjectif(constellation))
            self.notifierSansEvenement(constellation)
            self.somme_ecart_perturbation = 0
            while time()-self.start_date<self.tlim and it<config.getOptValue("max_iteration"):
                #self.setObjectif(self.calculerObjectif(constellation))
                it += 1 
                if time()-self.start_date<self.tlim and it<config.getOptValue("max_iteration"):
                    op_call,dernier_record = self.LNSsolve(constellation,mailbox,dernier_record,it)
                    application_operateur += op_call
                    score_avant = self.getSolutionContainer().objectif[0]
                    if time()-self.start_date<self.tlim:
                        printOpen("Perturber",c='y')
                        self.perturbateur.appliquer(self,constellation)
                        n_perturbations += 1
                        printClose()
                    if config.getOptValue("verif"):
                        self.verifierCCA(constellation)
                    score_apres = self.getSolutionContainer().objectif[0]
                    self.somme_ecart_perturbation += (score_apres-score_avant)
                    if score_avant>score_apres:  
                         self.nDegradation += 1
                    else:
                        self.nNonDegrade += 1  
            try:
                self.perturbation_moyenne = self.somme_ecart_perturbation/(self.nDegradation+self.nNonDegrade)           
            except:
                self.perturbation_moyenne = 0
            #print("Degradation_moyenne :",self.perturbation_moyenne)
            #print("Nombre de dégradation :",self.nDegradation)
            #print("Nombre de non dégradations :",self.nNonDegrade)
            self.verifierSolutionSiVerifMode(constellation)
            objectif,solCCAs,modes_retenus = self.getBest()
            #self.setObjectif(objectif)
            self.setAllSolCCAs(solCCAs)
            self.setModesRetenus(modes_retenus,constellation)
            
            self.notifierFinExecution(constellation)
            return application_operateur,it
                
        def redemarrer(self,constellation):
            self.initEtatRequetes(constellation)
            for r in self.etatRequetes:
                self.etatRequetes[r].initRechercheMode(self,constellation)
            super().redemarrer(constellation)
    
        def estStable(self,objectif):
            return objectif[0]<=self.solution.historique.getBest()[0][0]
        
        def estActif(self,lambda_v):
            return lambda_v[3]<=config.LNS.Tmax
    
        def reparer(self,constellation,req_detruites,entree_opt):
            printOpen("réparer : "+str(self.operateurReparation),c='b')
            requetes = self.operateurReparation.reparer(self,constellation,req_detruites,entree_opt)
            printColor(str(len(requetes))+" requêtes reconstruites",c='g')
            printClose()
            
        def detruire(self,constellation):
            printOpen("destruction : "+str(self.operateurDestruction),c='c')
            req_detruites,sortie_opt = self.operateurDestruction.detruire(self,constellation)
            printClose()
            return req_detruites,sortie_opt
        
        def appliquerOperateur(self,constellation):
            printColor('Appliquer '+str(self.operateur),c='c')
            res = self.operateur.appliquer(self,constellation)
            return res
            
        def validerSol(self,sol):
            cca_tmp_supp = 0
            for id_cca in sol:
                printColor("valider cca",id_cca,c='g')
                s,cca = id_cca
                self.setSolCCA(s,cca,sol[id_cca])
                for r in self.etatRequetes:
                    cca_tmp_supp += self.etatRequetes[r].resetCCA(id_cca,self.grapheDependances)
            #print(cca_tmp_supp,"cca temporaires supprimées")
    
        # confusion cout de transition et retard
        def aggregerScore(self,score,methode='r1'):
            if methode=='r1':
                return score[0]/(1+score[2])
            elif methode=='r2':
                return score[0]**1.25/(1+score[2])
            elif methode=='r0':
                return score[0]
            else:
                raise ValueError("methode d'aggregation inconnue")
            
        def rewardHeuristique(self,constellation,s,a):
            return constellation.getSatellite(s).getActivite(a).getScore()/(constellation.getSatellite(s).getActivite(a).getDuree()+1)
        
        def choisirRequete(self,candidats,random=False):
            if -1 in candidats:
                return -1
            if len(candidats)>0:
                if not random:
                    return max(candidats,key=lambda r : self.etatRequetes[r].getScoreCourant())
                else:
                    select = self.randomChoixRequete().randInt(len(candidats))
                    return candidats[select]
            else:
                return None
            
        def greedyFill(self,constellation,limit_req=True,forbid_solver=False,random=False):
            printOpen("remplissage glouton")
            candidats = [r for r in self.etatRequetes if constellation.getRequete(r).estActif() and r not in self.requetesCouvertes()]
            NReqCandidates = len(candidats)
            succes = 0
            r = self.choisirRequete(candidats,random=random)
            Nechecs = 0
            compteur = 0
            lim_notif = 50
            while r is not None and self.getTempsEcoule() < self.tlim and (not limit_req or Nechecs<config.getOptValue("max_echecs")):
                if time()-self.start_date>self.tlim:
                    printColor("interruption du glouton",c='y')
                    break
                printOpen("Recherche de chemin pour la requête ",r)
                mode,sol_cca,vide = self.etatRequetes[r].tenterInsererRequete(self,constellation,r,forbid_solver,self.modeleDeTransition)
                #printOpen("Insérer la requete",constellation.getRequete(r).getType(),r,c='g')
                if mode is not None:
                    #printColor("valider",noeuds_a_valider,c='g')
                    succes += 1
                    self.validerSol(sol_cca)
                    self.ajouterModeRetenu((r,mode.getId()),constellation)
                    candidats.remove(r)
                    if config.getOptValue("verif"):
                        self.verifierCCA(constellation)
                    printClose("Succès",c='g')
                else:
                    Nechecs += 1
                    if vide:
                        printColor("Plus de mode",c='m')
                        candidats.remove(r)
                    else:
                        printColor("Mode infaisable",c='m')
                    printClose("Echec",c='r')
                r = self.choisirRequete(candidats,random=random)
                if compteur==lim_notif or r is None:
                    compteur = 0
                    self.notifierSansEvenement(constellation)
                compteur += 1
                #self.notifierSansEvenement(constellation)
            printClose(str(succes)+"/"+str(NReqCandidates)+" requêtes insérées. "+str(Nechecs)+" échecs.",c='g')
            
            
    class Processus:
        def __init__(self,role,constellation,start_date,modeleDeTransition,dt_construction_transition,tlim=np.inf,CCAs=None,solution=None):
            self.role = role
            self.initSolution(constellation,start_date,modeleDeTransition,dt_construction_transition,tlim=tlim,CCAs=CCAs,solution=solution)
            
        def initSolution(self,constellation,start_date,modeleDeTransition,dt_construction_transition,tlim=np.Inf,CCAs=None,solution=None):
            self.solution = LNS(constellation,start_date,modeleDeTransition,dt_construction_transition,tlim=tlim,CCAs=CCAs,solution=solution)
            
        def resoudre(self,constellation,mailbox,afficher=True):
            application_operateur,it = self.solution.resoudre(constellation,mailbox,afficher)
            if config.getOptValue("version")!="coop":
                MPI.COMM_WORLD.Barrier()
            return application_operateur,it
    
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
            self.solution.saveSample(constellation,add=self.info_additionnelle)
            
        def getSolution(self):
            return self.solution.getSolution()
        
        def setSolution(self,sol,vid,modes_retenus,modes_candidats,modes_retires,objectif):
            self.solution.setSol(sol,vid,modes_retenus,modes_candidats,modes_retires,objectif)
        
        def verifierSolution(self):
            self.solution.verifierSolution()
            
        def getModesSolution(self):
            return self.solution.getModes()
        
    class Master(Processus):
        def __init__(self,constellation,start_date,modeleDeTransition,dt_construction_transition,tlim=np.Inf,CCAs=None,solution=None):
            super().__init__("master",constellation,start_date,modeleDeTransition,dt_construction_transition,tlim=tlim,CCAs=CCAs,solution=solution)
            
        def resoudre(self,constellation,mailbox,afficher=True):
            application_operateur,it = super().resoudre(constellation,mailbox,afficher)
            self.recolterResultats(mailbox)
            #self.terminerProcessus()
            self.solution.afficherInfo(time(),self.solution.start_date,constellation,color='y',title='FIN')
            self.solution.construirePlan(constellation)
         
        def recolterResultats(self,mailbox):
            hist = {}
            self.info_additionnelle = {}
            if config.getOptValue("version") in ["hybrid","coop"]:
                gaps,gains = self.solution.operateur.getGapsInfos()
            if config.getOptValue("version")!="coop":
                for data in mailbox.readMessages():
                    cpu = data['cpu']
                    hist[cpu] = data['hist']
                    if config.getOptValue("version") in ["hybrid","coop"]:
                        gains += data["gains"]
                        gaps += data["gaps"]
                if config.getOptValue("version") in ["hybrid","coop"] :
                    if len(gains)>0:
                        gains_positif_moyen = np.mean([x for x in gains if x>0])
                        pourcentage_succes = len([x for x in gains if x>0])/len(gains)
                    else:
                        gains_positif_moyen = 0
                        pourcentage_succes = 0
                    Nproblemes = len(gains)
                    gain_moyen = np.mean(gains)
                    self.info_additionnelle["gaps_moyen"] = gain_moyen
                    self.info_additionnelle["gains"] = gains_positif_moyen
                    self.info_additionnelle["fréquence_succès"] = pourcentage_succes
                    self.info_additionnelle["nombre_problèmes"] = Nproblemes
                self.solution.solution.historique.fusionner(hist)
    
    class Slave(Processus):
        def __init__(self,constellation,start_date,modeleDeTransition,dt_construction_transition,tlim=np.inf,CCAs=None,solution=None):
            comm = MPI.COMM_WORLD
            rank = comm.Get_rank()
            super().__init__("slave "+str(rank),constellation,start_date,modeleDeTransition,dt_construction_transition,tlim=tlim,CCAs=CCAs,solution=solution)
        
        def informationAdditionnelles(self):
            data = {}
            if config.getOptValue("version") in ["hybrid","coop"]:
                gaps,gains = self.solution.operateur.getGapsInfos()
                data["gaps"] = gaps
                data["gains"] = gains
            return data
        
        def formaterResultat(self):
            add = self.informationAdditionnelles()
            data = {'hist':self.solution.solution.historique,'cpu':MPI.COMM_WORLD.Get_rank()}
            hist = data['hist']
            for key in add:
                data[key] = add[key]
            return data
        
        def resoudre(self,constellation,mailbox):
            if config.getOptValue("version")!="coop":
                super().resoudre(constellation,mailbox)
                data = self.formaterResultat()
                mailbox.posterMessage(data)
        
        
class runnableLNS:
    def execute(self,constellation,start_date,modeleDeTransition,dt_construction_transition,tlim=np.inf,CCAs=None,solution=None):
        mailbox = MessagerieBloquante()
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        if rank == 0:
            self.process = Master(constellation,start_date,modeleDeTransition,dt_construction_transition,tlim=tlim,CCAs=CCAs,solution=solution)
            self.process.resoudre(constellation,mailbox)
            
        else:
            self.process = Slave(constellation,start_date,modeleDeTransition,dt_construction_transition,tlim=tlim,CCAs=CCAs,solution=solution)
            self.process.resoudre(constellation,mailbox) 
        return self.process.solution.getSolutionContainer()
    
    def getName(self):
        return "Large Neighborhood Search"
    
    def getMasterSolver(self):
        if MPI.COMM_WORLD.Get_rank()==0:
            return self.process.solution
        else:
            return None
        
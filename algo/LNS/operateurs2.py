#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 27 16:19:29 2022

@author: ssquilla
"""
            
class ReparateurCCA(OperateurReparation):
        def __init__(self):
            super().__init__("Reparateur de cca")
            id_mpi = MPI.COMM_WORLD.Get_rank()
            seed_option = config.getOptValue("seed")
            self.randomInstance = rd.Random(id_mpi+seed_option)
                
        def reparer(self,LNS,constellation,req_detruites,paircca):
            CCAS = list(paircca)
            self.randomInstance.shuffle(CCAS)
            printColor("Réparation des cca",CCAS,c='c')
            requetes = []
            self.modes_avant = {r : constellation.getRequete(r).getModeCourant(constellation).getId() for r in req_detruites}
            self.meilleur_mode_passe = {r : constellation.getRequete(r).getModeCourant(constellation).getId() for r in req_detruites if (r,constellation.getRequete(r).getModeCourant(constellation).getId()) in LNS.modes_retenus}
            for cca in CCAS:
                req = self.traiterCCA(LNS,constellation,req_detruites,cca,CCAS)
                for r in req:
                    if r not in requetes:
                        requetes.append(r)
            return requetes
        
        def pisterRequete(self,LNS,constellation,r):
            print("====================")
            print("requete ",r)
            for m in constellation.getRequete(r).modes:
                print(constellation.getRequete(r).modes[m])
            present = [m for (rr,m) in LNS.modes_retenus if rr==r]
            print("mode present dans la solution",present)
            if len(present)>0:
                for (s,a) in constellation.getRequete(r).getMode(present[0]).getCouples():
                    cca = LNS.grapheDependances.getActiviteCCA(a)
                    print(a,cca,a in LNS.solCCAs[s][cca].sequence)
            else:
                mode = constellation.getRequete(r).getModeCourant(constellation)
                for (s,a) in mode.getCouples():
                    cca = LNS.grapheDependances.getActiviteCCA(a)
                    print(a,cca,a in LNS.solCCAs[s][cca].sequence)
            print("dernier mode passe",self.meilleur_mode_passe.get(r,[]))
            print("====================")
            
        def insertionActivites(self,LNS,constellation,scores_requete,activites,s,cca):
            requetes_triees = sorted(scores_requete.keys(),key=lambda rr:scores_requete[rr],reverse=True)
            groups = {0:[]}
            faisable = True
            printColor("tentative d'insertion de "+str(sum([len(activites[r]) for r in activites]))+' activités')
            i = 1
            for r in requetes_triees:
                for a in activites[r]:
                    if faisable:
                        shiftRightDisplay(5)
                        faisable = LNS.solCCAs[s][cca].insererActivitePlansCritiques(constellation,a)
                        if faisable:
                            LNS.solCCAs[s][cca].MAJPlansCritiques(constellation)
                        shiftLeftDisplay(5)
                    else:
                        LNS.solCCAs[s][cca].insererActivite(constellation,a)
                groups[i] = activites[r]
                i += 1
            if not LNS.solCCAs[s][cca].sequenceFaisable(constellation):
                requetes_satisfaites = LNS.solCCAs[s][cca].rechercheLocale(constellation,solver='OPTW',groups=groups)
            else:
                requetes_satisfaites = requetes_triees
            return requetes_satisfaites
        
        def completerMode(self,constellation,r,activites,ccapair):
            if r in self.meilleur_mode_passe:
                ancien_mode = self.meilleur_mode_passe[r]
                anciennes_activites = constellation.getRequete(r).getMode(ancien_mode).getCouples()
                ancien_score = constellation.getRequete(r).getMode(ancien_mode).getRecompense()
            else:
                ancien_mode = constellation.getRequete(r).getModeCourant(constellation).getId()
                anciennes_activites = [] # mono/stereo : on considere que les activites de la cca
                ancien_score = 0
            anciennes_activites = [x for x in anciennes_activites if x not in activites]
            return anciennes_activites,ancien_score
        
        def explorerNouveauxModes(self,LNS,constellation,requetes_detruites,activites,s,cca,paircca):
            nouvelles_activites = {}
            scores_requete = {}
            delta = {} # difference de score entre le nouveau et l'ancien mode
            for r in requetes_detruites:
                liste_activites = activites.get(r,[])
                nouvelles_activites[r] = [(s,a) for a in liste_activites]
                if len(nouvelles_activites[r])>0:
                    anciennes_activites,score_ancien = self.completerMode(constellation,r,nouvelles_activites[r],paircca)
                    if constellation.getRequete(r).acceptable(anciennes_activites+nouvelles_activites[r],constellation):
                        score_nouveau_mode = constellation.getRequete(r).estimerScore(constellation,anciennes_activites+nouvelles_activites[r])[0]
                        scores_requete[r] = score_nouveau_mode 
                        delta[r] = score_nouveau_mode - score_ancien
            return nouvelles_activites,scores_requete,delta

        def retirerResteMode(self,LNS,constellation,r):
            m = constellation.getRequete(r).getModeCourant(constellation).getId()
            m_old = self.modes_avant[r]
            if m!= m_old:
                activites_a_retirer = [(s,a) for (s,a) in constellation.getRequete(r).getMode(m).getCouples() if (s,a) in constellation.getRequete(r).getMode(m_old).getCouples()]
                activites_par_cca = {}
                for (s,a) in activites_a_retirer:
                    cca = LNS.grapheDependances.getActiviteCCA(a)
                    if cca not in activites_par_cca:
                        activites_par_cca[cca] = []
                    activites_par_cca[cca].append(a)
                for cca in activites_par_cca:
                    s = LNS.grapheDependances.getSatelliteCCA(cca)
                    LNS.solCCAs[s][cca].retirerListeActivites(constellation,activites_par_cca[cca]) 
                
        def getNouvellesActivites(self,LNS,cca,requetes_detruites):
            activites = {}
            for a in LNS.grapheDependances.getComposanteConnexe(cca):
                r = constellation.getRequeteActivite(a)
                if r in requetes_detruites:
                    if r not in activites:
                        activites[r] = []
                    activites[r].append(a)
            return activites
        
        def validerNouveauxModes(self,LNS,constellation,nouvelles_activites,scores_requete,requetes_satisfaites,paircca):
            for r in requetes_satisfaites:
                old_act,score_ancien = self.completerMode(constellation,r,nouvelles_activites[r],paircca)
                nouveau_mode = constellation.getRequete(r).ajouterMode(constellation,old_act+nouvelles_activites[r])
                if r in self.meilleur_mode_passe and (r,self.meilleur_mode_passe[r]) in LNS.modes_retenus:
                    LNS.modes_retenus.remove((r,self.meilleur_mode_passe[r]))
                self.meilleur_mode_passe[r] = nouveau_mode.getId()
                LNS.modes_retenus.append((r,nouveau_mode.getId()))
                
            if config.getOptValue("verif"):
                for r in requetes_satisfaites:
                    for s,a in nouveau_mode.getCouples():
                        cca = LNS.grapheDependances.getActiviteCCA(a)
                        assert(a in LNS.solCCAs[s][cca].getSequence())
                        self.pisterRequete(LNS,constellation,r)
            
        def traiterCCA(self,LNS,constellation,requetes_detruites,cca,paircca):
            activites = self.getNouvellesActivites(LNS,cca,requetes_detruites)
            s = LNS.grapheDependances.getSatelliteCCA(cca)
            nouvelles_activites,scores_requete,delta = self.explorerNouveauxModes(LNS,constellation,requetes_detruites,activites,s,cca,paircca)
            requetes_satisfaites = self.insertionActivites(LNS,constellation,delta,activites,s,cca)
            self.validerNouveauxModes(LNS,constellation,nouvelles_activites,scores_requete,requetes_satisfaites,paircca)
            return requetes_satisfaites
                    
class EchangeurCCA(Operateur):
        def __init__(self,k=10):
            super().__init__("Echangeur de "+ str(k) + " pairs de cca")
            self.k = k
            self.destructeur = None
            self.constructeur = None
        
        def check(self,LNS,constellation,req,req_reconstruites):
            if config.getOptValue("verif"):
                LNS.verifierSolution(constellation)
                for r in req:
                    if r not in req_reconstruites:
                        pass
                        #act = constellation.getRequete(r).getMode(self.modes_avant[r]).getCouples()
                        #assert( r not in LNS.requetesCouvertes() or constellation.getRequete(r).acceptable(act))
                    else:
                        assert(r in LNS.requetesCouvertes())
            
        def appliquer(self,LNS,constellation):
            if self.destructeur is None:
                positions_cca_obs = {cca : LNS.positions_cca[cca] for cca in LNS.cca_obs}
                self.destructeur = LNS.DestructeurPairCCAGeographique(positions_cca_obs)
                self.constructeur = LNS.ReparateurCCA()
                
            for i in range(self.k):
            
                req,pair = self.destructeur.detruire(LNS,constellation)
                self.modes_avant = {r : constellation.getRequete(r).getModeCourant(constellation).getId() for r in req}
                types = {}
                for r in req:
                    t = constellation.getRequete(r).getType()
                    if t not in types:
                        types[t] = 0
                    types[t] += 1
                printColor(str(len(req))+" requêtes detruites : "+str(types),c='c')
                req_reconstruites = self.constructeur.reparer(LNS,constellation,req,pair)
                types = {}
                types_perdues = {}
                # stats des requetes
                for r in req:
                    t = constellation.getRequete(r).getType()
                    if r in req_reconstruites:
                        if t not in types:
                            types[t] = 0
                        types[t] += 1    
                    else:
                        if t not in types_perdues:
                            types_perdues[t] = 0
                        types_perdues[t] += 1        
                printColor(str(len(req_reconstruites))+" requêtes reconstruites : "+str(types),c='c')
                #requetes_perdues = [r for r in req if r not in req_reconstruites]
                printColor("requêtes perdues :",str(types_perdues),c='y')
                self.check(LNS,constellation,req,req_reconstruites)

            

class ExplorateurCP(Operateur):
        def __init__(self,k=2,iter_op=1):
            super().__init__("Echangeur CP de pairs de CCA")
            #self.k = k # choix d'une cca1 et d'une cca2 parmi les k plus proches
            self.positions_cca = None
            self.iter_op = iter_op
        
        def afficherModel(self,sol):
            #print("======= Variables ========")
            for cca in self.interval_vars:
                for r in self.interval_vars[cca]:
                    for a in self.interval_vars[cca][r]:
                        pass
                        #print(self.interval_vars[cca][r][a])
            #print("====== Contraintes =======")
            for e in self.model.get_all_expressions():
                pass
                #print(e)
            if sol.get_solve_status()=="Infeasible":
                self.model.export_model("model_erreur")
        
        def getEtatRequetes(self,LNS,constellation,cca1,cca2,req_en_commun,req_1_dynamique,req_2_dynamique):
            # dyn 1 : présent dans la séquence de 1 avec opportunités dans 2 = déplacable
            s1 = LNS.grapheDependances.getSatelliteCCA(cca1)
            s2 = LNS.grapheDependances.getSatelliteCCA(cca2)
            seq1 = [constellation.getRequeteActivite(a) for a in LNS.solCCAs[s1][cca1].getSequence()]
            seq2 = [constellation.getRequeteActivite(a) for a in LNS.solCCAs[s2][cca2].getSequence()]
            # requetes optionnelles
            requetes_optionnelles_1 = [ x for x in LNS.getRequetesPresentesSatiques(cca1) if x not in seq1 and x not in LNS.requetesCouvertes() and x!=-1]
            requetes_optionnelles_2 = [ x for x in LNS.getRequetesPresentesSatiques(cca2) if x not in seq2 and x not in LNS.requetesCouvertes() and x!=-1]
            rd.shuffle(requetes_optionnelles_1)
            rd.shuffle(requetes_optionnelles_2)
            N1 = min(config.LNS.max_requetes_echange_CP,len(requetes_optionnelles_1))
            N2 = min(config.LNS.max_requetes_echange_CP,len(requetes_optionnelles_2))
            optionnelles = {cca1:requetes_optionnelles_1[:N1],cca2:requetes_optionnelles_2[:N2]}
            # requetes fixes
            fixees = {}
            fixees[cca1] = [x for x in seq1 if x not in LNS.getRequetesPresentesSatiques(cca2)]
            fixees[cca2] = [x for x in seq2 if x not in LNS.getRequetesPresentesSatiques(cca1)]
            if -1 in req_en_commun:
                fixees[cca1].append(-1)
                fidees[cca2].append(-1)
                req_en_commun.remove(-1)   
            # requetes deplacables
            deplacables = req_en_commun
            
            #print("commun",req_en_commun)
            #print("r commun dans la sequence 1",req_1_dynamique)
            #print("r commun dans la sequence 2",req_2_dynamique)
            if config.glob.verif:
                for (s,cca) in zip ([s1,s2],[cca1,cca2]):
                    #print(cca)
                    #print("seq",[constellation.getRequeteActivite(a) for a in LNS.solCCAs[s][cca].getSequence()])
                    #print("activites dans la sequence non deplacables",fixees[cca])
                    #print("activites dans la sequence deplacables",deplacables)
                    for a in LNS.solCCAs[s][cca].getSequence():
                        r = constellation.getRequeteActivite(a)
                        assert(r in fixees[cca] + deplacables)
                        assert(r not in optionnelles[cca])
                    for r in fixees[cca]:
                        assert(r in [constellation.getRequeteActivite(a) for a in LNS.solCCAs[s][cca].getSequence()])
                        assert(r not in optionnelles[cca])
                    for r in deplacables:
                        #assert(r in [constellation.getRequeteActivite(a) for a in LNS.solCCAs[s][cca].getSequence()])
                        assert(r not in optionnelles[cca])
                    for r in optionnelles[cca]:
                        assert(r not in [constellation.getRequeteActivite(a) for a in LNS.solCCAs[s][cca].getSequence()])
            return deplacables,fixees,optionnelles
            
        def appliquer(self,LNS,constellation):
            open_neighbourhood = config.LNS.open_neighbourhood
            if config.glob.verif:
                for s in LNS.solCCAs:
                    for cca in LNS.solCCAs[s]:
                        assert(LNS.solCCAs[s][cca].sequenceFaisable(constellation))
            if self.positions_cca is None:
                self.positions_cca = {cca : LNS.positions_cca[cca] for cca in LNS.cca_obs}
            for i in range(self.iter_op):
                printOpen("Sélectionner une pair de CCAs")    
                (cca1,cca2),req_en_commun,req_1_dynamique,req_2_dynamique = self.selectionnerPairCCA(LNS, constellation)
                deplacables,fixees,optionnelles = self.getEtatRequetes(LNS,constellation,cca1,cca2,req_en_commun,req_1_dynamique,req_2_dynamique)
                printClose()
                printOpen("Initialiser le modèle")
                self.buildModel(LNS,constellation,deplacables,fixees,optionnelles,open_neighbourhood)
                printClose()
                printOpen("Résolution modèle CP")
                tlim = min(LNS.getTempsRestant(),config.LNS.CPtime)
                if tlim>0:
                    self.initSol(LNS,constellation)
                    sol = self.model.solve(TimeLimit = tlim,execfile=config.glob.docplexPath,Workers=1,log_output=None)
                    printClose()
                    if not (sol is not None and sol.get_solve_status()!="Infeasible" and sol.get_solve_status()!="Unknown"):
                        self.model.export_model("exemple_model/model_erreur_"+str((cca1,cca2)))
                        alert("EchangeCP : pas de solution trouvée pour "+str((cca1,cca2)))
                    gap = sol.solution.get_objective_gap()
                    printColor("Gap : ",gap,c='y')
                    score_avant = LNS.objectif[0]
                    self.remplirSolCCAs(constellation,LNS,sol)
                    self.MAJmodes(LNS, constellation, optionnelles, deplacables, fixees, open_neighbourhood)
                    LNS.verifierSolution(constellation)
                    delta_score = round(LNS.calculerObjectif(constellation)[0]-score_avant,3)
                    if delta_score<=10e-4:
                        printColor("Gain : " + str(delta_score),c='r')
                    else:
                        printColor("Gain : " + str(delta_score),c='g')
                    if delta_score<-10e-3:
                        self.model.export_model("exemple_model/model_erreur_"+str((cca1,cca2)))
                        alert("EchangeCP : pas d'incrément de score pour "+str((cca1,cca2)))    
                    if config.glob.verif:
                            LNS.verifierSolution(constellation)
        
        def infoSolutionInitiale(self,LNS,delta_score,sol):
            avant = []
            apres = []
            for cca in self.seq_cca:
                s = LNS.grapheDependances.getSatelliteCCA(cca)
                seq = LNS.solCCAs[s][cca].getSequence()
                apres.append(seq)
                avant.append(self.warmstart[self.seq_cca[cca]])
            score_avant = self.warmstart.get_objective_value()
            score_apres = sol.get_objective_value()
            print("avant",avant,score_avant)
            print("apres",score_apres)
            
        def initSol(self,LNS,constellation):
            warmstart=self.model.create_empty_solution()
            for cca in self.interval_vars:
                # indiquer le choix des modes
                s = LNS.grapheDependances.getSatelliteCCA(cca)
                plan = LNS.solCCAs[s][cca].planEarliest(constellation,LNS.solCCAs[s][cca].getSequence())
                for r in self.interval_vars[cca]:
                    for a in self.interval_vars[cca][r]:
                        if a in LNS.solCCAs[s][cca].getSequence():
                            if r in self.modes_var: # non présent si la requête est fixe
                                warmstart.set_value(self.modes_var[r][cca],1)
                            t = int(math.ceil([x[1] for x in plan if x[0]==a][0]*config.glob.getEchelle()))
                            duree = int(math.ceil(constellation.getSatellite(s).getActivite(a).getDuree()*config.glob.getEchelle()))
                            end = t+duree
                            warmstart.add_interval_var_solution(self.interval_vars[cca][r][a],presence=True,start=t,end=end)
                        else:
                            if r in self.modes_var: # peut ne pas etre presente sur cette cca
                                warmstart.set_value(self.modes_var[r][cca],0)
                            warmstart.add_interval_var_solution(self.interval_vars[cca][r][a],presence=False)
                # indiquer les séquences
                seq = []
                for a in LNS.solCCAs[s][cca].getSequence():
                    r = constellation.getRequeteActivite(a)
                    seq.append(self.interval_vars[cca][r][a])
                warmstart[self.seq_cca[cca]] = seq
            self.warmstart = warmstart
            self.model.set_starting_point(warmstart)
        
        def stoquerActiviteAjoutees(self,LNS,constellation,list_requetes):
            adds = {}
            for cca in self.seq_cca:
                s = LNS.grapheDependances.getSatelliteCCA(cca)
                for a in LNS.solCCAs[s][cca].getSequence():
                    r = constellation.getRequeteActivite(a)
                    if r not in adds:
                        adds[r] = []
                    adds[r].append((s,a))
            return adds
        
        def stoquerActiviteRetirees(self,LNS,constellation,deplacables,fixees,open_neighbourhood):
            if not open_neighbourhood:
                return {}
            removes = {}
            # activites sur la cca
            activites = {}
            req_fixees = []
            for cca in fixees:
                for r in fixees[cca]:
                    if r not in req_fixees:
                        req_fixees.append(r)
            for r in deplacables+req_fixees:
                m = [x[1] for x in LNS.modes_retenus if x[0]==r][0]
                mode_courant = constellation.getRequete(r).getMode(m)
                activites[r] = [x[1] for x in mode_courant.getCouples() if LNS.grapheDependances.getActiviteCCA(x[1]) in self.seq_cca]
                for cca in self.seq_cca:
                    act_cca = [x for x in activites[r] if LNS.grapheDependances.getActiviteCCA(x)==cca]
                    s = LNS.grapheDependances.getSatelliteCCA(cca)
                    manquantes = [x for x in act_cca if x not in LNS.solCCAs[s][cca].getSequence()]
                    if len(manquantes)>0:
                        if r not in removes:
                            removes[r] = []
                        removes[r] += manquantes
            return removes
          
        def creerNouveauxModes(self,LNS,constellation,adds,removes):
            nouveau_contenu_mode = {}
            CCAS = list(self.seq_cca.keys())
            # requetes a mettre a jour
            requetes = sorted(list(adds.keys()))
            for r in removes:
                if r not in requetes:
                    requetes.append(r)
            for r in requetes:
                if r!=-1:
                    # recuperer l'ancien mode
                    old = [m for (rr,m) in LNS.modes_retenus if rr==r]
                    if len(old)>0:
                        old_mode = constellation.getRequete(r).getMode(old[0])
                        LNS.modes_retenus.remove((r,old[0]))
                        assert(len(old)==1)
                    else:
                        old_mode = None
                    # creer le contenu du nouveau mode 
                    nouveau_contenu_mode[r] = adds.get(r,[])
                    if old_mode is not None:
                        couples = old_mode.getCouples()
                        nouveau_contenu_mode[r] += [(sat,act) for (sat,act) in couples if act not in removes.get(r,[]) and LNS.grapheDependances.getActiviteCCA(act) not in CCAS]
                    if constellation.getRequete(r).acceptable(nouveau_contenu_mode[r],constellation):
                        # creer le mode et l'ajouter
                        new_mode = constellation.getRequete(r).ajouterMode(constellation,nouveau_contenu_mode[r])
                        LNS.modes_retenus.append((r,new_mode.getId()))
                    else:
                        assert(len(nouveau_contenu_mode[r])==0) # parce que ca a été calculé pour etre accapetable. Seule possibilité.
        # mettre a jour les modes APRES le remplissage de la solution
        def MAJmodes(self,LNS,constellation,opt,deplacables,fixees,open_neighbourhood):
            requetes_MAJ = deepcopy(deplacables)
            for cca in opt:
                requetes_MAJ += opt[cca]
            adds = self.stoquerActiviteAjoutees(LNS, constellation, requetes_MAJ)
            removes = self.stoquerActiviteRetirees(LNS,constellation,deplacables,fixees,open_neighbourhood)
            self.creerNouveauxModes(LNS,constellation,adds,removes)
        
        def remplirSolCCAs(self,constellation,LNS,sol):
            for cca in self.seq_cca:
                sequence = []
                for it in sol.get_value(self.seq_cca[cca]):
                    sequence.append(eval(it.get_name().split("_")[1])[2])
                LNS.ecraserSolutionCCA(cca,sequence)
            
        def buildModel(self,LNS,constellation,deplacables,fixees,optionnelles,open_neighbourhood):
            self.model = CpoModel()
            activites_fixes,activites_mobiles,activites_optionnelles = self.initVars(LNS,constellation,deplacables,fixees,optionnelles,open_neighbourhood)#fixed_cca,req_en_commun,opt)
            self.initContraintesSolutionBase(constellation,activites_fixes,activites_mobiles,activites_optionnelles,LNS)
            self.initNoOverlap(constellation,LNS)
            self.initRequetesVariables(LNS,constellation,deplacables,fixees,optionnelles,open_neighbourhood)
            
        def initVars(self,LNS,constellation,deplacables,fixees,opt,open_neighbourhood):
            activites_fixes = {}
            activites_mobiles = {}
            activites_optionnelles = {}
            self.interval_vars = {}
            for cca in fixees:
                self.interval_vars[cca] = {}
                for r in fixees[cca] + deplacables + opt[cca]:
                    if r not in self.interval_vars[cca]:
                        self.interval_vars[cca][r] = {}
                s = LNS.grapheDependances.getSatelliteCCA(cca)
                activites_fixes[cca] = {}
                activites_mobiles[cca] = {}
                activites_optionnelles[cca] = {}
                for a in LNS.grapheDependances.getComposanteConnexe(cca):
                    r = constellation.getRequeteActivite(a)
                    if r in fixees[cca] + deplacables + opt[cca]:
                        start = int(math.ceil(config.glob.getEchelle()*constellation.getSatellite(s).getActivite(a).getDebut()))
                        end = int(math.floor(config.glob.getEchelle()*constellation.getSatellite(s).getActivite(a).getFin()))
                        duree = int(math.ceil(config.glob.getEchelle()*constellation.getSatellite(s).getActivite(a).getDuree()))
                        self.interval_vars[cca][r][a] = self.model.interval_var(start=(start,end),end=(start,end),length=duree,name="Ia_"+str((cca,r,a)))
                        assert(start+duree<=end)
                        if (not open_neighbourhood and r in fixees[cca]) or r==-1:
                            if r not in activites_fixes[cca]:
                                activites_fixes[cca][r] = []
                            activites_fixes[cca][r].append(a)
                        elif r in deplacables:
                            if r not in activites_mobiles[cca]:
                                activites_mobiles[cca][r] = []
                            activites_mobiles[cca][r].append(a)
                        else:
                            if r not in activites_optionnelles[cca]:
                                activites_optionnelles[cca][r] = []
                            activites_optionnelles[cca][r].append(a)
            self.modes_var = {}                
            for r in deplacables:
                self.modes_var[r] = {}
                for cca in fixees:
                    self.modes_var[r][cca] = self.model.binary_var(name="y_"+str((r,cca)))
                    assert(r in self.interval_vars[cca])
            if open_neighbourhood:
                for cca in fixees:
                    for r in fixees[cca]:
                        if r not in self.modes_var:
                            self.modes_var[r] = {}
                        if cca not in self.modes_var[r]:
                            self.modes_var[r][cca] = self.model.binary_var(name="y_"+str((r,cca)))
            for cca in opt:
                for r in opt[cca]:
                    if r not in self.modes_var:
                        self.modes_var[r] = {}
                    self.modes_var[r][cca] = self.model.binary_var(name="y_"+str((r,cca)))
                    assert(r in self.interval_vars[cca])
                
            return activites_fixes,activites_mobiles,activites_optionnelles

    def initContraintesSolutionBase(self,constellation,activites_fixes,activites_mobiles,activites_optionnelles,LNS):
            for cca in activites_mobiles:            
                for r in activites_mobiles[cca]:
                    for a in activites_mobiles[cca][r]:
                        self.interval_vars[cca][r][a].set_optional()
            for cca in activites_optionnelles:            
                for r in activites_optionnelles[cca]:
                    for a in activites_optionnelles[cca][r]:
                        self.interval_vars[cca][r][a].set_optional()
            rec_a = lambda a :constellation.getSatellite(constellation.getSatelliteActivite(a)).getActivite(a).getScore()
            obj1 = self.model.sum([self.model.presence_of(self.interval_vars[cca][r][a])*rec_a(a) for cca in self.interval_vars for r in self.interval_vars[cca] for a in self.interval_vars[cca][r]])
            obj_nouveaute = self.CritereNouveaute(LNS)
            self.model.maximize(obj1)
        
    def CritereNouveaute(self,LNS):
            seq = []
            for cca in self.interval_vars:
                s = LNS.grapheDependances.getSatelliteCCA(cca)
                for r in self.interval_vars[cca]:
                    for a in self.interval_vars[cca][r]:
                        if a not in LNS.solCCAs[s][cca].getSequence():
                            seq.append(self.model.presence_of(self.interval_vars[cca][r][a]))
            return self.model.sum(seq)
                            
    def initNoOverlap(self,constellation,LNS):
            self.seq_cca = {}
            for cca in self.interval_vars:
                s = LNS.grapheDependances.getSatelliteCCA(cca)
                activites = [(a,r,self.interval_vars[cca][r][a]) for r in self.interval_vars[cca] for a in self.interval_vars[cca][r]]
                n = len(activites)
                distance_matrix = np.zeros((n,n),dtype=int)
                vars_list = []
                for i,(a1,r,var) in enumerate(activites):
                    vars_list.append(var)
                    for j,(a2,r,var2) in enumerate(activites):
                        if i==j:
                            distance_matrix[i][j] = 0
                        else:
                            transition = config.glob.getEchelle()*constellation.getSatellite(s).getTransition(a1,a2)
                            distance_matrix[i][j] = int(math.ceil(transition))
                mat = self.model.transition_matrix(distance_matrix)
                self.seq_cca[cca] = self.model.sequence_var(vars_list)
                self.model.add(self.model.no_overlap(self.seq_cca[cca],distance_matrix=mat))
        
    def initRequeteOneShot(self,r,ccas,deplacables,optionnelles,open_neighbourhood):
            if r in deplacables and not open_neighbourhood: 
                self.model.add_constraint(self.model.sum([self.modes_var[r][cca] for cca in self.modes_var[r]])==1)
            else:
                assert(open_neighbourhood or ((r in optionnelles[ccas[0]] or r in optionnelles[ccas[1]])and not open_neighbourhood) )
                self.model.add_constraint(self.model.sum([self.modes_var[r][cca] for cca in self.modes_var[r]])<=1)
        
    def initRequeteSystematic(self,constellation,r,ccas,deplacables,optionnelles,open_neighbourhood):
            REQUETE = constellation.getRequete(r) 
            couples_mode_precedent = REQUETE.getModeCourant(constellation).getCouples()
            sat_mode_precedent,Ncouples = couples_mode_precedent[0][0],len(couples_mode_precedent)
            if Ncouples>=2:
                for cca in ccas:
                    s = LNS.grapheDependances.getSatelliteCCA(cca)
                    if s != sat_mode_precedent:
                        self.model.add_constraint(self.modes_var[r][cca]==0)
                if r in deplacables and not open_neighbourhood:
                    self.model.add_constraint(self.model.sum([self.modes_var[r][cca] for cca in self.modes_var[r]])>=1)
                else:
                    assert((r in optionnelles[ccas[0]] or r in optionnelles[ccas[1]])and not open_neighbourhood or open_neighbourhood ) 
     
    def initRequetePeriodique(self,constellation, r,ccas,deplacables,optionnelles,open_neighbourhood): 
            REQUETE = constellation.getRequete(r)                  
            if r in deplacables and not open_neighbourhood:
                self.model.add_constraint(self.model.sum([self.modes_var[r][cca] for cca in self.modes_var[r]])>=1)
            else:
                assert((r in optionnelles[ccas[0]] or r in optionnelles[ccas[1]])and not open_neighbourhood or open_neighbourhood )  
            TS = REQUETE.getTimeSlots()
            A = [a for cca in self.interval_vars for a in self.interval_vars[cca][r]]
            activites = [a for (s,a) in REQUETE.getModeCourant(constellation).getCouples() if a not in A]
            time_slots_deja_presents = [REQUETE.getTimeSlot(a) for a in activites]
            for t in TS:
                slot_act = [x[1] for x in TS[t]]
                var_act = [self.interval_vars[cca][r][a] for cca in self.interval_vars for a in self.interval_vars[cca].get(r,[]) if a in slot_act]
                if t in time_slots_deja_presents:
                    self.model.add_constraint(self.model.sum([self.model.presence_of(x) for x in var_act])==0)
                else:
                    self.model.add_constraint(self.model.sum([self.model.presence_of(x) for x in var_act])<=1)
    
    def initPresenceRequetesCCA(self,constellation,r):
            REQUETE = constellation.getRequete(r)
            for cca in self.modes_var[r]:
                N = len(self.interval_vars[cca][r])
                if REQUETE.getType() in ["ONE_SHOT_MONO","SYSTEMATIC","PERIODIC"]:
                    if N!=1:
                        raise ValueError("plusieurs candidats par cca",list(self.interval_vars[cca][r].keys()),REQUETE.getType())
                elif REQUETE.getType()=="STEREO":
                    if N!=2:
                        raise ValueError("un seul candidat sur la cca",list(self.interval_vars[cca][r].keys()),REQUETE.getType())
                self.model.add_constraint(self.model.sum([self.model.presence_of(self.interval_vars[cca][r][a]) for a in self.interval_vars[cca][r] ])==N*self.modes_var[r][cca])
       
    def initRequetesVariables(self,LNS,constellation,deplacables,fixees,optionnelles,open_neighbourhood):
            ccas = list(fixees.keys())
            for r in self.modes_var:
                REQUETE = constellation.getRequete(r)
                if REQUETE.getType() in ["ONE_SHOT_MONO","ONE_SHOT_STEREO"]:
                    self.initRequeteOneShot(r,ccas,deplacables,optionnelles,open_neighbourhood)
                elif REQUETE.getType()=="SYSTEMATIC":
                    self.initRequeteSystematic(constellation,r,ccas,deplacables,optionnelles,open_neighbourhood)
                else:
                    assert(REQUETE.getType()=="PERIODIC")
                    self.initRequetePeriodique(constellation,r,ccas,deplacables,optionnelles,open_neighbourhood)
                self.initPresenceRequetesCCA(constellation,r)
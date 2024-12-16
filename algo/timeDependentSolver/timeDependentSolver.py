#!/usr/bin/env python3
# -*- coding: utf-8 -*-
    
from mpi4py import MPI

from ..model.constellation import Constellation
from ..model.solution import Solver
from ..model.transitionModel import *

from ..CPModel.CPSolver import runnableCPSolver
from ..LNS.LNSSolver import runnableLNS
from ..BPCCAS.BPCCAS import runnableBPCCAS

import random as rd

from ..Utils.Utils import *
from ..Utils.config import *
config = Config()
instance = config.instance

from copy import copy

if config.getOptValue("help"):
    if MPI.COMM_WORLD.Get_rank()==0:
        config.displayHelp()
else:
    
    """
    Created on Mon Jan  9 11:31:56 2023
    
    @author: ssquilla
    """
    class Messagerie:
        def __init__(self):
            # eviter d'envoyer les informations qu'on va recuperer à l'identique ensuite
            self.local_data = {}
        
        def envoyerMessage(self,data,cpu):
            comm = MPI.COMM_WORLD
            rank = comm.Get_rank()
            assert(rank==0)
            data['time']=time()
            data['cpu'] = rank
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
                
    
    class MessagerieBloquante(Messagerie):
        def __init__(self):
            pass
        
        def posterMessage(self,data):
            assert(MPI.COMM_WORLD.Get_rank()>0)
            MPI.COMM_WORLD.send(data,dest=0)
            
        def readMessages(self):
            assert(MPI.COMM_WORLD.Get_rank()==0)
            for i in range(MPI.COMM_WORLD.Get_size()-1):
                data = MPI.COMM_WORLD.recv()
                #data['reception'] = (data['reception'],time()-data['reception'])
                yield data
        
        def __str__(self):
            return "Messagerie bloquante"
    
    class TimeDependantSolver(Solver):
        def __init__(self,constellation,start_date,modeleDeTransition,dt_construction_transition,solution=None,CCAs=None):
            super().__init__(constellation,start_date,modeleDeTransition,dt_construction_transition,solution=solution,CCAs=CCAs)
            display_depth = getDisplayDepth()
            shiftLeftDisplay(display_depth-1)
            self.initSolverInitPhase(constellation,start_date)
            self.initSolverFillingPhase(constellation,start_date)
            
        def initSolverInitPhase(self,constellation,start_date):  
            tinit = config.getOptValue("initDur")
            solver_init = config.getOptValue("initSolver")
            if solver_init=="CPSolver":
                self.initSolver = runnableCPSolver()
            elif solver_init =="BPCCAS":
                self.initSolver = runnableBPCCAS()
            elif solver_init =="LNS":
                self.initSolver = runnableLNS()
            else:
                raise NameError("Solver inconnu.")
            printMaster("Durée allouée au solver",self.initSolver.getName(), ":",tinit,"(s)",c='y')
        
        def initSolverFillingPhase(self,constellation,start_date):
            if config.getOptValue("chain"):
                solver_filling = config.getOptValue("fillSolver")
            else:
                solver_filling = None
            if solver_filling is not None:
                if solver_filling=="CPSolver":
                    self.fillSolver = runnableCPSolver()
                elif solver_filling =="BPCCAS":
                    self.fillSolver = runnableBPCCAS()
                elif solver_filling =="LNS":
                    self.fillSolver = runnableLNS()
                else:
                    raise NameError("Solver inconnu.")
            else:
                self.fillSolver = None
        
            
        def afficherCouple(self,constellation,modeleTimeDep):
            print(constellation.getSatellite(1).getActivite(123))
            print(constellation.getSatellite(1).getActivite(125))
        
        def verifierTempsTransition(self,composantes,constellation,modeleTimeDep,modeleApproche):  
            self.afficherCouple(constellation,modeleTimeDep)
            for s in constellation.satellites:
                for a1 in constellation.getSatellite(s).getActivites():
                    for a2 in constellation.getSatellite(s).getActivites():
                        print(s,a1,a2,constellation.getSatellite(s).getActivite(a1).getFin(),modeleTimeDep.getPointsDeControle(a1,a2),modeleApproche.getTransition(s,a1,a2,4))
                        if composantes.getActiviteCCA(a1)==composantes.getActiviteCCA(a2):
                            for (t,duree) in modeleTimeDep.getPointsDeControle(a1,a2):
                                if config.getOptValue("initTransitionModel")=="slow":
                                    assert(duree<=modeleApproche.getTransition(s,a1,a2,4))
                                elif config.getOptValue("initTransitionModel")=="fast":
                                    assert(duree>=modeleApproche.getTransition(s,a1,a2,4))
            
        def resoudre(self,constellation,start_date,modeleTimeDependent):
            add = {'algorithme choisi':self.initSolver.getName(),'durée allouée':str(config.getOptValue("initDur"))+" (s)"}
            self.afficherInfo(time(),start_date,constellation,core=0,color='r',add=add,title='PHASE 1 : '+self.initSolver.getName(),filtre=['time'])
            
            
            tinit = config.getOptValue("initDur")
            self.solution = self.initSolver.execute(constellation,start_date,self.modeleDeTransition,self.dt_construction_transition,tlim=tinit,CCAs=self.grapheDependances)

            if config.verifMode():
                pass
                #self.verifierTempsTransition(self.grapheDependances,constellation,modeleTimeDependent,self.modeleDeTransition)

            printMaster("Changement de modèle de transition : ",modeleTimeDependent)
            restant = round(config.getOptValue("time")-(time()-start_date),0)
            add = {'durée restante':str(restant)+" (s)"}
                        
            self.afficherInfo(time(),start_date,constellation,core=0,color='r',add=add,title='PHASE 2 : réparation',filtre=['best','time','modes','requete'])
            self.reparerSolution(constellation,modeleTimeDependent)
            self.verifierSolutionSiVerifMode(constellation,modeleDeTransition=modeleTimeDependent)
            
            self.notifierFinExecution(constellation)
            
            self.verifierSolutionSiVerifMode(constellation,modeleDeTransition=modeleTimeDependent)
            
            if self.fillSolver is not None:
                restant = round(config.getOptValue("time")-(time()-start_date),0)
                add = {'algorithme choisi':self.fillSolver.getName(),'durée restante':str(restant)+" (s)"}
                self.afficherInfo(time(),start_date,constellation,core=0,color='r',add=add,title='PHASE 3 : '+self.fillSolver.getName(),filtre=['best','time','modes','requete'])
                self.solution = self.fillSolver.execute(constellation,start_date,modeleTimeDependent,self.dt_construction_transition,CCAs=self.grapheDependances,solution=self.solution)
            else:
                self.afficherInfo(time(),start_date,constellation,core=0,color='y',add=add,title='FIN DE LA PHASE 2',filtre=['best','time','modes','requete'])
            return self.solution

        def reparerSolution(self,constellation,modeleDeTransition):
            self.restartBestSolution(constellation)
            sequences_reparees,sequences_intactes,activites_retirees = self.reparerSequences(constellation,modeleDeTransition)
            self.invaliderMeilleureSolution()
            constellation.filtrerModesPresents(self.solution,self.grapheDependances,modeleDeTransition)
            printMaster(str(sequences_reparees)+" séquences réparées. "+str(sequences_intactes) + " séquences intactes. " + str(activites_retirees)+" activités retirées.")
            self.solution.ajouterInformationAdditionnelle("sequences réparées",sequences_reparees)
            self.solution.ajouterInformationAdditionnelle("sequences intactes",sequences_intactes)
            self.solution.ajouterInformationAdditionnelle("activités retirées",activites_retirees)
            self.solution.ajouterInformationAdditionnelle("score après réparation",self.solution.getObjectif())
        
        
        def reparerSequences(self,constellation,modeleDeTransition):
            sequences_reparees = 0
            sequences_intactes = 0
            activites_retirees = 0
            for s in self.getSolCCAs():
                for cca in self.getSolCCAs()[s]:
                    faisable = self.getSolCCA(s,cca).sequenceFaisable(constellation,modeleDeTransition)
                    if not faisable and config.getOptValue("initTransitionModel")=="slow":
                        for i,a in enumerate(self.getSolCCA(s,cca).sequence):
                            if i <len(self.getSolCCA(s,cca).sequence)-1:
                                a2 = self.getSolCCA(s,cca).sequence[i+1]
                                print(modeleDeTransition.getTransitionMax(a,a2),modeleDeTransition.getPointsDeControle(a,a2))
                        
                        
                        self.getSolCCA(s,cca).sequenceFaisable(constellation,modeleDeTransition,print_retard=True)
                        print(self.getSolCCA(s,cca).retardSequence(constellation,self.getSolCCA(s,cca).sequence,modeleDeTransition,print_info=True))
                        raise ValueError("model pessimiste. séquence anormalement infaisable.",self.getSolCCA(s,cca).getSequence())
                    if not faisable:
                        sequence = self.getSolCCA(s,cca).getSequence()
                        sequence_reparee,nRemoves = self.reparationRecursive(constellation, sequence, s, modeleDeTransition)
                        self.getSolCCA(s, cca).setSequence(constellation,sequence_reparee,modeleDeTransition)
                        if config.verifMode():
                            assert(self.getSolCCA(s,cca).sequenceFaisable(constellation,modeleDeTransition))
                        sequences_reparees += 1
                        activites_retirees += nRemoves
                    else:
                        sequences_intactes += 1
            return sequences_reparees,sequences_intactes,activites_retirees        
        
        def decompositionPessimiste(self,constellation,sequence,s,modeleDeTransition):
            return [sequence]
        
        def reparationRecursive(self,constellation,sequence,s,modeleDeTransition):
            faisable = self.faisabiliteTemporelleSequence(constellation,sequence,s,modeleDeTransition=modeleDeTransition)
            nRemoves = 0
            if faisable:
                return sequence,0
            else:
                liste_d = [a for a in sequence if constellation.getSatellite(s).estVidage(a)]
                self.retirerActiviteOptimale(constellation,sequence,s,modeleDeTransition)
                nRemoves += 1
                for d in liste_d:
                    assert(d in sequence)
                decomposition = self.decomposerSequenceV2(constellation,sequence,s,modeleDeTransition=modeleDeTransition)
                res = [self.reparationRecursive(constellation,seq,s,modeleDeTransition=modeleDeTransition) for seq in decomposition]
                sequences = self.fusionnerSequence([x[0] for x in res])
                nRemoves = nRemoves + sum([x[1] for x in res])
                return sequences,nRemoves
            
        def fusionnerSequence(self,sequences_ordonnees):
            seq = []
            for sequence in sequences_ordonnees:
                seq += sequence
            return seq
          
        def decomposerSequenceV2(self,constellation,sequence,s,modeleDeTransition):
            decomposition = []
            sequence_courante = []
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
                    if t==start and len(sequence_courante)>0:
                        decomposition.append(sequence_courante.copy())
                        sequence_courante = []
                sequence_courante.append(activite)
            decomposition.append(sequence_courante)
            return decomposition
        
        def decomposerSequence(self,constellation,sequence,s,modeleDeTransition):
            decomposition = []
            sequence_courante = []
            for i,a in enumerate(sequence):
                debut_a = constellation.getActivite(a).getDebut()
                if i>0:
                    prec = sequence[i-1]
                    duree_prec = constellation.getActivite(prec).getDuree()
                    transition_prec = self.getTransition(constellation,s,prec,a,t_prec+duree_prec,modeleDeTransition)
                else:
                    transition_prec = 0
                    duree_prec = 0
                    t_prec = debut_a
                
                if len(sequence_courante)==0:
                    sequence_courante.append(a)
                else:
                    precision = config.glob.getEchelle()
                    # activites liées
                    if t_prec + duree_prec + transition_prec + precision > debut_a:
                        sequence_courante.append(a)
                    # activités indépendantes
                    else:
                        decomposition.append(copy(sequence_courante))
                        sequence_courante = [a]
                # date de début de a. Utilisée comme la date de l'activité précédente à l'itération suivante de la boucle
                t_prec = max(t_prec + duree_prec+transition_prec,debut_a)
            
            decomposition.append(sequence_courante)
            if config.getOptValue("verif"):
                assert(len(sequence)==sum([len(seq) for seq in decomposition]))
            
            return decomposition       
                    
        def retirerActiviteOptimale(self,constellation,sequence,s,modeleDeTransition):
            critere = lambda rla : (rla[0]+rla[1])/(constellation.getSatellite(s).getActivite(rla[2]).getScore()+1)
            critere_recompense = critere = lambda rla : 1/(constellation.getSatellite(s).getActivite(rla[2]).getScore()+1)
            # retard = liste de (retard droite,retard gauche,a) 
            retard = self.retardGaucheDroiteSequenceParObs(constellation,sequence,s,modeleDeTransition=modeleDeTransition)
            # ne pas considérer les vidages
            for i in range(len(retard)-1,-1,-1):
                if constellation.getSatellite(s).estVidage(retard[i][2]):
                    retard.pop(i)            
            rla_opti = max(retard,key=critere_recompense)
            if config.verifMode():
                assert(not constellation.getSatellite(s).estVidage(rla_opti[2]))
            sequence.remove(rla_opti[2])
            
    class runnableTimeDependentSolver:
        def execute(self,constellation,start_date,modeleApproche,modeleTimeDependent,dt_construction_transition,solution=None,CCAs=None):
            assert(modeleTimeDependent is not None)
            # initDur = infini => pas de limite. fin = temps global
            if config.getOptValue("time")<config.getOptValue("initDur") and not np.isinf(config.getOptValue("initDur")):
                if MPI.COMM_WORLD.Get_rank()==0:
                    warn("Indiquez une durée totale (paramètre --time) supérieure au temps d'execution du solver d'initialisation (paramètre --initDur).")
                    warn("Arrêt.")
                self.solver = None
            else:
                self.solver = TimeDependantSolver(constellation,start_date,modeleApproche,dt_construction_transition,solution=solution,CCAs=CCAs)
                self.solver.resoudre(constellation,start_date,modeleTimeDependent)
        
        def getName(self):
            return "Time Dependent Solver"
    
        def getMasterSolver(self):
            if MPI.COMM_WORLD.Get_rank()==0:
                return self.solver
            else:
                return None
        
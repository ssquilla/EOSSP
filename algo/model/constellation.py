from mpi4py import MPI
from math import pi
from math import cos
from math import sin
from math import sqrt
from math import exp
from math import floor
from math import ceil
import numpy as np
from itertools import chain, combinations

from functools import reduce
import itertools
from operator import itemgetter

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.patches as patches
import matplotlib

from copy import copy
from copy import deepcopy
import random as rd
import collections
from time import time
from time import sleep
import getopt
import sys
import os
import subprocess
import heapq

from .requestDealer import RequestDealer


from ..Utils.Utils import *
from ..Utils.config import *
config = Config()

from .graphes import *

global TYPE_PERIODIC
TYPE_PERIODIC = "PERIODIC"
global TYPE_LONG_MONO
TYPE_LONG_MONO = "LONG_MONO"
global TYPE_MONO
TYPE_MONO = "ONE_SHOT_MONO"
global TYPE_STEREO
TYPE_STEREO = "ONE_SHOT_STEREO"
global TYPE_SYSTEMATIC
TYPE_SYSTEMATIC = "SYSTEMATIC"
global TYPE_VIDAGE
TYPE_VIDAGE = "DOWNLOAD"



"""
        =========================================================
                    Attributs de la constellation 
        =========================================================
"""
class Mode: # obs = {s:[o]}
    def __init__(self,idRequete,idMode,recompense,obs,scoreTemporel):
        if not config.glob.score_temporel:
            if not scoreTemporel==0:
                print("requête",idRequete,"score temporel hors domaine.")
                print("mode :",idRequete,idMode,"contenu :",obs)
            assert(scoreTemporel==0)
        self.idRequete = idRequete
        self.idMode = idMode
        self.recompense = recompense
        self.scoreTemporel = scoreTemporel
        assert(recompense>0 or self.idRequete == -1)
        self.observations = obs
        self.couples = []
        for s in self.observations:
            for o in self.observations[s]:
                self.couples.append((s,o))
    
    def acceptable(self,constellation):
        return constellation.getRequete(self.idRequete).acceptable(self.couples,constellation)
        
    def getScoreTemporel(self):
        return self.scoreTemporel
    
    def getCouples(self):
        return self.couples
    
    def getId(self):
        return self.idMode

    def getRequete(self):
        return self.idRequete

    def getRecompense(self):
        return self.recompense
        
    def getActivites(self):
        return self.observations
    
    def getActivitesSatellite(self,s):
        return self.observations.get(s,[])
    
    def __str__(self):
        return "Mode "+str((self.idRequete,self.idMode)) +" : "+str(self.observations) + " recompense = " +str(self.recompense)


class RelationOrdreRequete:
    def __init__(self,r_id,mode_id,score_bruite):
        self.r_id = r_id
        self.mode_id = mode_id
        self.score_bruite = score_bruite

    def getId(self):
        return self.r_id
    
    def getIdMode(self):
        return self.mode_id
    
    def getRecompense(self):
        return self.score_bruite
    
    def __eq__(self,other):
        if self.getId() ==-1:
            if other.getId() == -1:
                assert(self.getRecompense()[0] == 0 and other.getRecompense()[0] ==0)
                return True
            else:
                return False
        else:
            return self.getId() == other.getId() and self.getIdMode() == other.getIdMode()
    
    def __str__(self):
        return "Order("+str(self.r_id)+","+str(self.mode_id)+":"+str(self.score_bruite)+")"
    
    def __ne__(self,other):
        return not self.__eq__(other)
    
    def __lt__(self, other):
        if self.getId()==-1:
            return False
        else:
            return (self.getRecompense() < other.getRecompense()) or (self.getRecompense() == other.getRecompense() and self.getId() > other.getId())

    def __gt__(self,other):
        if self.getId()==-1:
            if other.getId()==-1:
                return False
            else:
                return True
        else:
            if other.getId()==-1:
                return False
            else:
                return (self.getRecompense() > other.getRecompense()) or (self.getRecompense() == other.getRecompense() and self.getId() < other.getId())

    def __ge__(self,other):
        return self.__gt__(other) or self.__eq__(other)
        
    def __le__(self,other):
        return self.__lt__(other) or self.__eq__(other)

class Requete:
        def __init__(self,idRequete,priorite,type_requete,activites):
            self.actif = not config.getOptValue("dynamic")
            self.idRequete = idRequete
            self.priorite = priorite
            self.setType(type_requete)
            seed_option = config.getOptValue("seed")
            self.generateur_aleatoire_observations = rd.Random(seed_option)
            self.generateur_aleatoire_liste_observations = rd.Random(seed_option)
            self.generateur_aleatoire_modes = rd.Random(seed_option)
            self.init = False # modes déjà initialisés = False
            # time-slot -> liste de triplets (reward,satellite,indice activite)
            self.activitesParTimeSlots = activites # {timeSlot:[(s,o,reward)]}
            # satellites -> liste des indices des activités
            self.liste_activites_par_satellite= {}
            self.couples_satellites_activites = []
            self.map_rso = {}
            self.modes = {}
            self.idMode = 0
            for timeSlot in self.activitesParTimeSlots:
                for i,(r,s,o) in enumerate(self.activitesParTimeSlots[timeSlot]):
                    self.map_rso[o] = (timeSlot,i)
                    if(s,o) not in self.couples_satellites_activites:
                        self.couples_satellites_activites.append((s,o))
                    if s not in self.liste_activites_par_satellite:
                        self.liste_activites_par_satellite [s] = []
                    if o not in self.liste_activites_par_satellite [s]:
                        self.liste_activites_par_satellite [s].append(o)
            self.mode_candidat = None
            self.maintient = {}
            self.modes_generes = 0
            assert(len(self.couples_satellites_activites)==sum([len(self.liste_activites_par_satellite [s]) for s in self.liste_activites_par_satellite]))

        def setType(self,type_requete):
            self.type_requete = type_requete
        
        def estActif(self):
            return self.actif
        
        def activer(self):
            self.actif = True
            
        def desactiver(self):
            self.actif = False
            
        def getCCAPresentes(self,grapheDep):
            res = []
            for (s,a) in self.couples_satellites_activites:
                cca = grapheDep.getActiviteCCA(a)
                if cca not in res:
                    res.append(cca)
            return res
               
        def resetCandidats(self):
            self.maintient = {}
            self.idMode = 0
            self.mode_candidat = None
        
        def getModesGeneres(self):
            return self.modes_generes
        
        def findActivite(self,a):
            timeSlot,i = self.map_rso[a]
            return self.activitesParTimeSlots[timeSlot][i]
        
        def getModes(self):
            return self.modes.keys()
        
        def getMode(self,m):
            return self.modes[m]
            
        def getScoreActivite(self,constellation,a):
            return self.findActivite(a)
            
        def getTimeSlot(self,a):
            return self.map_rso[a][0]
         
        def trierParTimeSlot(self,liste_act):
            tri = {}
            for a in liste_act:
                ts = self.getTimeSlot(a)
                if ts not in tri:
                    tri[ts] = []
                tri[ts].append(a)
            return tri
                
        def creerGrapheRequete(self):
            type_requete = self.getType()
            r = self.idRequete
            if type_requete == TYPE_MONO:
                graphe = GrapheMono(self)
            elif type_requete == TYPE_STEREO:
                graphe = GrapheStereo(self)
            elif type_requete == TYPE_PERIODIC:
                graphe = GraphePeriodic(self)
            elif type_requete == TYPE_SYSTEMATIC:
                graphe = GrapheSystematic(self)
            elif type_requete == TYPE_VIDAGE:
                graphe = GrapheVidage(self)
            else:
                raise ValueError("type inconnu",type_requete)
            self.graphe = graphe
            self.contenu_a_valider = None
        
        def getGraphe(self):
            return self.graphe
        
        def getObservationsStructures(self,contenu):
            listeObs = []
            obs = {}
            timeSlots = {}
            for timeSlot in self.activitesParTimeSlots:
                for (r,s,a) in self.activitesParTimeSlots[timeSlot]:
                    if a in contenu:
                        if s not in obs:
                            obs[s] = []
                        obs[s].append(a)
                        listeObs.append((r,s,a))
                        timeSlots[(r,s,a)]=timeSlot
            return listeObs,obs,timeSlots
        
        """
        # Utilisé dans le LNS : on ne sait pas encore si le mode sera effectivement candidat
        def genererMode(self,rewards,inactives):
            recompense_chemin,contenu,noeuds_a_valider = self.graphe.plusLongChemin(rewards,inactives)
            self.contenu_a_valider = contenu.copy()
            return recompense_chemin,contenu,noeuds_a_valider
        """
        def validerModeCandidat(self,constellation):
            mode = self.mode_a_valider
            self.mode_a_valider = None
            self.modes[self.idMode] = mode
            self.idMode += 1
            self.modes_generes += 1
            self.mode_candidat = mode
            return mode
        
        def annulerMode(self):
            self.contenu_a_valider = None
            self.mode_candidat = None
            self.modes[self.idMode] = None
        
        # Utilisé dans l'ancienne version de LNS
        def genererModeEtValider(self,constellation,rewards,inactives):
            recompense_chemin,contenu,noeuds_a_valider = self.graphe.plusLongChemin(rewards,inactives)
            if len(contenu)>0:
                listeObs,obs,timeSlots = self.getObservationsStructures(contenu)
                recompense,scoreTemporel = self.scoreListeObservations(listeObs,constellation,timeSlots)
                mode = Mode(self.getId(),self.idMode,recompense,obs,scoreTemporel)
                self.modes[self.idMode] = mode
                self.idMode += 1
                self.modes_generes += 1
                self.mode_candidat = mode
                return mode,noeuds_a_valider
            else:
                self.mode_candidat = None
                return None,[]
            
        def maintenirObs(self,explication):
            self.maintient = {}
            for (s,o) in self.mode_candidat.getCouples():
                if o not in explication:
                    if s not in self.maintient:
                        self.maintient[s] = []
                    self.maintient[s].append(o)
            #print("maint",self.maintient)
        
        def stats(self,constellation):
            xyz = [0,0,0]
            minCoord = [np.Inf,np.Inf,np.Inf]
            maxCoord = [-np.Inf,-np.Inf,-np.Inf]
            nPoints = 0
            for timeSlot in self.activitesParTimeSlots:
                for (r,s,o) in self.activitesParTimeSlots[timeSlot]:
                    coords = constellation.getSatellite(s).getActivite(o).getCoordonnees()
                    xyz[0] += coords[0]
                    xyz[1] += coords[1]
                    xyz[2] += coords[2]
                    minCoord = [min(minCoord[i],coords[i]) for i in range(3)]
                    maxCoord = [max(maxCoord[i],coords[i]) for i in range(3)]
                    nPoints += 1
            xyz = tuple([elmt/nPoints for elmt in xyz])
            return xyz,minCoord,maxCoord,nPoints
        
        def getMeanTarget(self,constellation):
            xyz = [0,0,0]
            nPoints = 0
            for timeSlot in self.activitesParTimeSlots:
                for (r,s,o) in self.activitesParTimeSlots[timeSlot]:
                        coords = constellation.getSatellite(s).getActivite(o).getCoordonnees()
                        xyz[0] += coords[0]
                        xyz[1] += coords[1]
                        xyz[2] += coords[2]
                        nPoints += 1
            xyz = tuple([elmt/nPoints for elmt in xyz])
            return xyz           
            
        def statsCoordonnees(self,constellation):
            xyz = [0,0,0]
            minCoord = [np.Inf,np.Inf,np.Inf]
            maxCoord = [-np.Inf,-np.Inf,-np.Inf]
            nPoints = 0
            for timeSlot in self.activitesParTimeSlots:
                for (r,s,o) in self.activitesParTimeSlots[timeSlot]:
                        coords = constellation.getSatellite(s).getActivite(o).getCoordonnees()
                        xyz[0] += coords[0]
                        xyz[1] += coords[1]
                        xyz[2] += coords[2]
                        minCoord = [min(minCoord[i],coords[i]) for i in range(3)]
                        maxCoord = [max(maxCoord[i],coords[i]) for i in range(3)]
                        nPoints += 1
            xyz = tuple([elmt/nPoints for elmt in xyz])
            return xyz,minCoord,maxCoord,nPoints
            
        def getCouples(self):
            return self.couples_satellites_activites
        
        def getPriorite(self):
            return self.priorite
        
        def getModeCourant(self,constellation):
            if not self.init:
                self.resetModes(constellation)
            return self.mode_candidat

        def getActivites(self):
            return self.liste_activites_par_satellite
        
        def getActivitesSatellite(self,s):
            return self.liste_activites_par_satellite.get(s,[])

        def getTimeSlots(self):
            return list(self.activitesParTimeSlots.keys())
        
        def getObservationsParTimeSlot(self,timeSlot):
            return self.activitesParTimeSlots[timeSlot]

        def getId(self):
            return self.idRequete

        def getRecompense(self):
            return self.recompense

        def getMode(self,m):
            return self.modes[m]

        def getModes(self):
            return self.modes

        def getType(self):
            return self.type_requete
        
        def getKModesSuivants(self,constellation,k):
            liste = []
            for i in range(k):
                if i<len(self.modes):
                    liste.append(self.modes[i])
                else:
                    res = self.getModeSuivantSansExp(constellation)
                    if res is not None:
                        liste.append(res)
                    else:
                        break
            return liste
        
        def selectionnerObs(self,constellation,listeObs,cleTemporelle,generateur=None):
            assert(len(listeObs)>0)
            # listeObs : liste de (r,s,o) 
            critere = lambda rso : self.scoreObs(rso,constellation,cleTemporelle)
            r_max = max([critere(rso) for rso in listeObs])
            obs_max = [rso for rso in listeObs if critere(rso)==r_max]
            #printColor("choix parmi",obs_max,'activites',c='m')
            if generateur is None:
                gen = self.generateur_aleatoire_observations
            else:
                gen = generateur
            id_obs = gen.randint(0,len(obs_max)-1)
            return obs_max[id_obs]
        
        def selectionnerListeObs(self,constellation,listeObs,cleTemporelle,k,generateur=None):
            assert(len(listeObs)>0)         
            # listeObs : liste de (r,s,o) 
            critere = lambda rso : self.scoreObs(rso,constellation,cleTemporelle)
            #r_max = max([critere(rso) for rso in listeObs])
            #obs_max = [rso for rso in listeObs if critere(rso)==r_max]
            #printColor("choix parmi",obs_max,'activites',c='m')
            if generateur is None:
                gen = self.generateur_aleatoire_liste_observations
            else:
                gen = generateur
            # stoquer les valeurs
            values = {}
            for elmt in listeObs:
                value = critere(elmt)
                if value not in values:
                    values[value] = []
                values[value].append((value,elmt))
            # conserver les meilleurs groupes de valeurs  
            elements_selectionnes = []
            nElmt = 0
            for value in sorted(values,reverse=True):
                if nElmt + len(values[value])<=k:
                    nElmt += len(values[value])
                    elements_selectionnes += values[value]
                else:
                    gen.shuffle(values[value])
                    NelmtsDesires = k - nElmt
                    elements_selectionnes += values[value][:NelmtsDesires]
                    break
            return elements_selectionnes       
        
        def scoreObs(self,rso,constellation,cleTemporelle):
            # low : temps de reference
            t = (constellation.getSatellite(rso[1]).getActivite(rso[2]).getDebut()+constellation.getSatellite(rso[1]).getActivite(rso[2]).getFin())/2
            if config.scoring.method=='relative':
                low = min([constellation.getSatellite(s).getActivite(a).getDebut() for (s,a) in self.couples_satellites_activites])
                up = max([constellation.getSatellite(s).getActivite(a).getFin() for (s,a) in self.couples_satellites_activites])
            elif config.scoring.method=='global':
                low,up = config.donnees.tmin,config.donnees.tmax
            alpha = config.scoring.alpha_weather
            score = (rso[0],cleTemporelle(rso),0)
            assert(score>(0,0,0))
            return score
        
        def nettoyerModesInutiles(self):
            for i in range(self.idMode):
                if i > 0 and i < self.idMode-1 and i in self.modes:
                    del self.modes[i]
        
        # L'interet est d'ajouter un mode en indice le plus haut
        # Il est la copie du mode d'indice m. Utile quand m est le meilleur mode
        # connu pour les algo qui utilise le derniermode  comme point de départ  
        # dans la recherche locale.
        def decalerMode(self,m):
            if self.idRequete==-1:
                return self.idMode
            mode = deepcopy(self.modes[m])
            mode.idMode = self.idMode
            self.mode_candidat = mode
            assert(self.idMode not in self.modes)
            self.modes[self.idMode] = mode
            self.idMode += 1
            self.modes_generes += 1
            #self.nettoyerModesInutiles()
            return mode.getId()
        
        # version générique : très couteuse
        def ajouterMode(self,constellation,liste_activites):
            if config.getOptValue("verif"):
                assert(self.acceptable(liste_activites,constellation))
            couples = len(liste_activites)==0 or type(liste_activites[0])==tuple
            if self.idRequete==-1:
                raise ValueError("impossible de créer de nouveaux modes pour la requête de vidage")
            rso_liste = []
            obs = {}
            TS = {}
            for timeSlot in self.activitesParTimeSlots:
                for r,s,a in self.activitesParTimeSlots[timeSlot]:
                    if (couples and (s,a) in liste_activites) or (not couples and a in liste_activites):
                        TS[a] = timeSlot
                        rso_liste.append((r,s,a))
                        if s not in obs:
                            obs[s] = []
                        obs[s].append(a)
            assert(len(rso_liste)>0)
            lexico = self.scoreListeObservations(rso_liste, constellation, None)       
            mode = Mode(self.idRequete,self.idMode,lexico[0],obs,lexico[1])
            
            self.modes[self.idMode] = mode
            self.mode_candidat = mode
            self.idMode += 1
            self.modes_generes += 1
            return mode
 
        def estimerScore(self,constellation,liste_activites):
            rso_liste = []
            obs = {}
            TS = {}
            for timeSlot in self.activitesParTimeSlots:
                for r,s,a in self.activitesParTimeSlots[timeSlot]:
                    if (s,a) in liste_activites:
                        TS[a] = timeSlot
                        rso_liste.append((r,s,a))
                        if s not in obs:
                            obs[s] = []
                        obs[s].append(a)
            assert(len(rso_liste)>0)
            return self.scoreListeObservations(rso_liste, constellation, None)
    
        def getActivitesPresentesCCA(self,grapheDep,liste_cca):
            activites_cca = []
            for timeSlot in self.activitesParTimeSlots:
                for (r,s,a) in self.activitesParTimeSlots[timeSlot]:
                    if not config.getOptValue("dynamic") or self.estActif():
                        cca = grapheDep.getActiviteCCA(a)
                        if cca in liste_cca:
                            activites_cca.append(a)
            return activites_cca
        
        # contenu actuel : activites du mode retenu privé de celles sur les cca a explorer
        def genererModesPresentsCCA(self,grapheDep,constellation,contenu_actuel,liste_cca_a_explorer,allow_external_modes=True):
            external_modes = False
            external_activites = len(contenu_actuel)>0
            nouveaux_modes = []
            activites_cca = self.getActivitesPresentesCCA(grapheDep, liste_cca_a_explorer)
            set_cca = set(activites_cca)
            max_size_ajout = len(activites_cca)
            combinaisons = []
            for taille in range(max_size_ajout+1):
                combinaisons += list(map(list,itertools.combinations(activites_cca,taille)))
            for combi in combinaisons:
                liste_nouveau_mode = contenu_actuel + combi
                if len(combi)>0 or (allow_external_modes and self.acceptable(liste_nouveau_mode,constellation)):
                    nouveaux_modes.append(self.ajouterMode(constellation, liste_nouveau_mode))
                    if len(combi)==0:
                        external_modes = True
            return nouveaux_modes,external_modes,external_activites
            
        def filtrerModesPresents(self,liste_activites,constellation):
            self.idMode = 0
            if config.verifMode():
                assert(self.acceptable(liste_activites,constellation))
            return self.ajouterMode(constellation,liste_activites)
    
class RequeteMono(Requete):
    def __init__(self,idRequete,priorite,activites):
        super().__init__(idRequete,priorite,TYPE_MONO,activites)

    def genererModesPresentsCCA(self,grapheDep,constellation,contenu_actuel,liste_cca_a_explorer,allow_external_modes=True):
        external_activites = len(contenu_actuel)>0
        if len(contenu_actuel)==1:
            if allow_external_modes:
                return [],False,external_activites
            else:
                return []
        assert(len(contenu_actuel)==0)
        nouveaux_modes = []
        activites_cca = self.getActivitesPresentesCCA(grapheDep, liste_cca_a_explorer)
        for a in activites_cca:
            liste_nouveau_mode = [a]
            nouveaux_modes.append(self.ajouterMode(constellation,liste_nouveau_mode))
        return nouveaux_modes,False,external_activites
    
    def getBestModeWithoutInactives(self,constellation,inactifs):
        if not config.getOptValue("dynamic") or self.estActif():
            activites = [rso for rso in self.activites_candidates if rso[2] not in inactifs]
            if len(activites)>0:
                cleTemporelle = lambda rso : self.scoreTemporel(constellation,rso)
                rso = self.selectionnerObs(constellation,activites,cleTemporelle)
                recompense,scoreTemporel = self.scoreListeObservations([rso],constellation,cleTemporelle)
                obs =  {rso[1]:[rso[2]]}
                self.mode_a_valider = Mode(self.idRequete,self.idMode,recompense,obs,scoreTemporel)
                return self.mode_a_valider
            else:
                self.mode_a_valider = None
                return None
        else:
            return None
    
    def acceptable(self,liste_activites,constellation):
        couples = len(liste_activites)==0 or type(liste_activites[0])==tuple
        if couples:
            act = self.couples_satellites_activites
        else:
            act = [x[1] for x in self.couples_satellites_activites]
        if not sum([x in act for x in liste_activites])==len(liste_activites):
            return False
        return len(liste_activites)==1
    
    def getKMeilleursModes(self,constellation,k):
        assert(self.modes=={})
        liste_modes = []
        if len(self.activites_candidates)==0:
            self.mode_candidat = None
            return []
        else:
            start = time()
            cleTemporelle = lambda rso : self.scoreTemporel(constellation,rso)
            liste_observations = self.selectionnerListeObs(constellation,self.activites_candidates,cleTemporelle,k,generateur=None)
            end = time()

            for (score,rso) in liste_observations:
                start = time()
                recompense,scoreTemporel,_ = score
                obs_mode =  {rso[1]:[rso[2]]}
                self.mode_candidat = Mode(self.idRequete,self.idMode,recompense,obs_mode,scoreTemporel)
                self.modes[self.idMode] = self.mode_candidat
                self.idMode += 1    
                self.modes_generes += 1
                liste_modes.append(self.mode_candidat)
                end = time()
        return liste_modes

    def resetModes(self,constellation,conserve_old=True,initFirstMode=True):
        self.init = True
        self.resetCandidats()
        self.activites_candidates = deepcopy(self.activitesParTimeSlots[0])
        self.first_date = min([constellation.getSatellite(rso_o[1]).getActivite(rso_o[2]).getDebut() for rso_o in self.activites_candidates])
        self.last_date = max([constellation.getSatellite(rso_o[1]).getActivite(rso_o[2]).getFin() for rso_o in self.activites_candidates])
        if not conserve_old:
            self.modes = {}
            self.idMode = 0        
        if initFirstMode:
            self.construireMode(constellation)
        
    def scoreListeObservations(self,listeObs,constellation,timeSlots):
        assert(len(listeObs)==1)
        cleTemporelle = lambda x : self.scoreTemporel(constellation,x)
        lexico = [self.scoreObs(rso,constellation,cleTemporelle) for rso in listeObs][0]
        return (lexico[0],lexico[1])
     
    def defilerActivites(self,explication):
        assert(explication!=[])
        len_avant = len(self.activites_candidates)
        #copie = self.activites_candidates.copy()
        for i,x in enumerate(self.activites_candidates):
            if x[2] in explication:
                self.activites_candidates.pop(i)
                break
        #self.activites_candidates = [x for x in self.activites_candidates if x[2] not in explication]
        len_apres = len(self.activites_candidates)
        assert(len_avant == len_apres + len(explication))
        
    def construireMode(self,constellation):
        if len(self.activites_candidates)==0:
            self.mode_candidat = None
        else:
            cleTemporelle = lambda rso : self.scoreTemporel(constellation,rso)
            rso = self.selectionnerObs(constellation,self.activites_candidates,cleTemporelle)
            recompense,scoreTemporel = self.scoreListeObservations([rso],constellation,cleTemporelle)
            obs =  {rso[1]:[rso[2]]}
            self.mode_candidat = Mode(self.idRequete,self.idMode,recompense,obs,scoreTemporel)
            self.modes[self.idMode] = self.mode_candidat
            self.idMode += 1
            self.modes_generes += 1
     
    def scoreTemporel(self,constellation,rso):
        if config.glob.score_temporel:
            t = constellation.getSatellite(rso[1]).getObservation(rso[2]).getDebut() - self.first_date
            start = 3*60*60
            end = self.last_date - self.first_date
            if t<start:
                return 1
            else:
                return 1 - (t-start)/(end-start)
        else:
            return 0
     
    def getModeSuivant(self,explication,constellation):
        if self.mode_candidat is not None:
            self.defilerActivites(explication)
            self.construireMode(constellation)
        return self.mode_candidat
    
    def getModeSuivantSansExp(self,constellation):
        if self.mode_candidat is None:
            return None
        exp = [self.mode_candidat.getCouples()[0][1]] # 1er couple (s,o) : obs
        return self.getModeSuivant(exp,constellation)

class RequeteLongMono(RequeteMono):
    def __init__(self,idRequete,priorite,activites):
        self.multiplicateur = 4
        for t in activites:
            for i in range(len(activites[t])):
                activites[t][i] = (activites[t][i][0]*self.multiplicateur,activites[t][i][1],activites[t][i][2])
        super().__init__(idRequete, priorite, activites)
        self.setType(TYPE_LONG_MONO) # ecraser le type
    
class RequeteStereo(Requete):
    def __init__(self,idRequete,priorite,activites):
        super().__init__(idRequete,priorite,TYPE_STEREO,activites)
        assert(len(self.couples_satellites_activites)==sum([len(self.liste_activites_par_satellite [s]) for s in self.liste_activites_par_satellite]))

    def genererModesPresentsCCA(self,grapheDep,constellation,contenu_actuel,liste_cca_a_explorer,allow_external_modes=True):
        if not config.getOptValue("dynamic") or self.estActif():
            external_activites = len(contenu_actuel)>0
            if len(contenu_actuel)==2:
                if allow_external_modes:
                    return [],False,external_activites
                else:
                    return []
            assert(len(contenu_actuel)==0)
            nouveaux_modes = []
            activites_cca = self.getActivitesPresentesCCA(grapheDep, liste_cca_a_explorer)
            tri_par_pairs = self.trierParTimeSlot(activites_cca)
            for idPair in tri_par_pairs:
                assert(len(tri_par_pairs[idPair])==2)
                liste_nouveau_mode = tri_par_pairs[idPair]
                #print(liste_nouveau_mode)
                nouveaux_modes.append(self.ajouterMode(constellation,liste_nouveau_mode))
            return nouveaux_modes,False,external_activites
        else:
            return [],False,len(contenu_actuel)>0
            
            
    def acceptable(self,liste_activites,constellation):
        if len(liste_activites)!=2:
            return False
        couples = len(liste_activites)==0 or type(liste_activites[0])==tuple
        if couples:        
            t1 = self.getTimeSlot(liste_activites[0][1])
            t2 = self.getTimeSlot(liste_activites[1][1])
        else:
            t1 = self.getTimeSlot(liste_activites[0])
            t2 = self.getTimeSlot(liste_activites[1])
        return t1 == t2 
    
    def getCouplesStereo(self):
        couples = []
        for ts in self.activitesParTimeSlots:
            couples.append([(x[1],x[2]) for x in self.activitesParTimeSlots[ts].copy()])
        return couples
    
    def selectionnerCoupleActivites(self,listeCouples,constellation):
        r_max = max(self.scoreCouple(x,constellation) for x in listeCouples)
        couples_max = [x for x in listeCouples if self.scoreCouple(x,constellation) == r_max]
        id_couple = self.generateur_aleatoire_observations.randint(0,len(couples_max)-1)
        return couples_max[id_couple]
    
    def getBestModeWithoutInactives(self,constellation,inactifs):
        if not config.getOptValue("dynamic") or self.estActif():
            activites = [rso_couples for rso_couples in self.couples if rso_couples[0][2] not in inactifs and rso_couples[1][2] not in inactifs]
            if len(activites)==0:
                return None
            cleTemporelle = lambda rso : self.scoreTemporel(constellation,rso)
            rso_liste = self.selectionnerCoupleActivites(activites,constellation)
            recompense,scoreTemporel = self.scoreListeObservations(rso_liste,constellation,cleTemporelle)
            assert(len(rso_liste)==2) # 2 obs
            assert(rso_liste[0][1] == rso_liste[1][1]) # même sat
            obs =  {rso_liste[0][1]:[rso_liste[0][2],rso_liste[1][2]]}
            mode =  Mode(self.idRequete,self.idMode,recompense,obs,scoreTemporel)
            if self.acceptable(mode.getCouples(),constellation):
                self.mode_a_valider = mode
                return self.mode_a_valider
            else:
                return None
        else:
            None
        
    def scoreListeObservations(self,couple,constellation,timeSlots):
        assert(len(couple)==2)
        cleTemporelObs1 = lambda rso : self.scoreTemporelCouple(constellation,(rso,(0,0,0)))
        score1 = self.scoreObs(couple[0],constellation,cleTemporelObs1)
        cleTemporelleObs2 = lambda rso : 0
        score2 = self.scoreObs(couple[1],constellation,cleTemporelleObs2)
        return (score1[0] + score2[0],score1[1])
            
    def getIdPaires(self):
        return list(self.activitesParTimeSlots.keys())

    def getPaireObservation(self,idPaire):
        return self.activitesParTimeSlots[idPaire]
    
    def scoreCouple(self,couple,constellation):
        cleTemporelObs1 = lambda rso : self.scoreTemporelCouple(constellation,(rso,(0,0,0)))
        score1 = self.scoreObs(couple[0],constellation,cleTemporelObs1)
        cleTemporelleObs2 = lambda rso : 0
        score2 = self.scoreObs(couple[1],constellation,cleTemporelleObs2)
        return (score1[0] + score2[0],score1[1])
    
    def resetModes(self,constellation,conserve_old=True,initFirstMode=True):
        self.init = True
        self.resetCandidats()
        self.couples = [self.activitesParTimeSlots[timeSlot] for timeSlot in self.activitesParTimeSlots] 
        self.first_date = min([constellation.getSatellite(rso_o[1]).getActivite(rso_o[2]).getDebut() for couple in self.couples for rso_o in couple])
        self.last_date = max([constellation.getSatellite(rso_o[1]).getActivite(rso_o[2]).getFin() for couple in self.couples for rso_o in couple])
        # a = (r,s,o) : observations differents avec le meme sat
        if not conserve_old:
            self.modes = {}
            self.idMode = 0          
        if initFirstMode:
            self.construireMode(constellation)
    
    def selectionnerListeCouples(self,constellation,listeCouples,cleTemporelle,k,generateur=None):
        assert(len(listeCouples)>0)
        # listeObs : liste de (r,s,o) 
        critere = lambda couple : self.scoreCouple(couple,constellation)
        if generateur is None:
            gen = self.generateur_aleatoire_liste_observations
        else:
            gen = generateur
        # stoquer les valeurs
        values = {}
        for elmt in listeCouples:
            value = critere(elmt)
            if value not in values:
                values[value] = []
            values[value].append((value,elmt))
        # conserver les meilleurs groupes de valeurs  
        elements_selectionnes = []
        nElmt = 0
        for value in sorted(values,reverse=True):
            if nElmt + len(values[value])<=k:
                nElmt += len(values[value])
                elements_selectionnes += values[value]
            else:
                gen.shuffle(values[value])
                NelmtsDesires = k - nElmt
                elements_selectionnes += values[value][:NelmtsDesires]
                break
        return elements_selectionnes  

    def getKMeilleursModes(self,constellation,k):
        assert(self.modes=={})
        liste_modes = []
        if len(self.couples)==0:
            self.mode_candidat = None
            return []
        else:
            start = time()
            cleTemporelle = lambda rso : self.scoreTemporel(constellation,rso)
            liste_observations = self.selectionnerListeCouples(constellation,self.couples,cleTemporelle,k,generateur=None)
            end = time()
            for (score,couple) in liste_observations:
                start = time()
                recompense,scoreTemporel = score
                rso1,rso2 = couple
                obs_mode =  {rso1[1]:[rso1[2],rso2[2]]}
                self.mode_candidat = Mode(self.idRequete,self.idMode,recompense,obs_mode,scoreTemporel)
                self.modes[self.idMode] = self.mode_candidat
                self.idMode += 1 
                self.modes_generes += 1
                liste_modes.append(self.mode_candidat)
                end = time()
        return liste_modes
    
    def testExplicationInCouple(self,explication,couple):
        for a in explication:
            if a == couple[0][2] or a == couple[1][2]:
                return True
        return False
        
    def defilerActivites(self,explication):
        for i in range(len(self.couples)-1,-1,-1):
            couple = self.couples[i]
            if self.testExplicationInCouple(explication,couple):
                self.couples.pop(i)
                if config.verifMode():
                    for c in self.couples:
                        if( self.testExplicationInCouple(explication,c)):
                            print(self.idRequete,explication,c)
                            assert(not self.testExplicationInCouple(explication,c))
                break
        
    def construireMode(self,constellation):
        if len(self.couples)==0:
            self.mode_candidat = None
        else:
            couple = self.selectionnerCoupleActivites(self.couples,constellation)
            recompense,scoreTemporel = self.scoreCouple(couple,constellation)
            obs = {couple[0][1]:[couple[0][2],couple[1][2]]}
            self.mode_candidat = Mode(self.idRequete,self.idMode,recompense,obs,scoreTemporel)
            self.modes[self.mode_candidat.getId()] = self.mode_candidat
            self.idMode += 1
            self.modes_generes += 1
    
    def scoreTemporelCouple(self,constellation,couple):
        if config.glob.score_temporel:
            rso = couple[0]
            t = constellation.getSatellite(rso[1]).getObservation(rso[2]).getDebut() - self.first_date
            start = 3*60*60
            end = self.last_date - self.first_date
            if t<start:
                return 1
            else:
                return 1 - (t-start)/(end-start)
            
        else:
            return 0
        
    def getModeSuivant(self,explication,constellation):
        assert(explication!=[])
        if self.mode_candidat is not None:
            self.defilerActivites(explication)
            self.construireMode(constellation)
        return self.mode_candidat

    def getModeSuivantSansExp(self,constellation):
        if self.mode_candidat is None:
            return None
        exp = [self.mode_candidat.getCouples()[0][1],self.mode_candidat.getCouples()[1][1]] # 1er couple (s,o) : obs
        return self.getModeSuivant(exp,constellation)
    
class RequetePeriodic(Requete):
    def __init__(self,idRequete,priorite,activites):
        super().__init__(idRequete,priorite,TYPE_PERIODIC,activites)
        assert(len(self.couples_satellites_activites)==sum([len(self.liste_activites_par_satellite [s]) for s in self.liste_activites_par_satellite]))
        seed_option = config.getOptValue("seed")
        self.generateur_plot = {plot : rd.Random(plot+seed_option) for plot in self.activitesParTimeSlots}

    def genererModesPresentsCCA(self,grapheDep,constellation,contenu_actuel,liste_cca_a_explorer,allow_external_modes=True):
        if not config.getOptValue("dynamic") or self.estActif():
            external_modes = False
            external_activites = len(contenu_actuel)>0
            # I. créer les obs candidates à l'ajout
            time_slots_occupes = [self.getTimeSlot(a) for a in contenu_actuel]
            nouveaux_modes = []
            activites_cca = self.getActivitesPresentesCCA(grapheDep, liste_cca_a_explorer)
            tri_par_slot = self.trierParTimeSlot(activites_cca)
            for ts in time_slots_occupes:
                if ts in tri_par_slot:
                    del tri_par_slot[ts]
            # II. générer les combinaisons possibles d'ajout d'obs
            combinaisons = [[]]
            for ts in tri_par_slot:
                aux = deepcopy(combinaisons)
                for a in tri_par_slot[ts]:
                    aux_2 = deepcopy(aux)
                    for i in range(len(aux_2)):
                        aux_2[i].append(a)
                    combinaisons += aux_2
            # III. combiner les combinaisons au contenu actuel
            for combi in combinaisons:
                liste_nouveau_mode = combi+contenu_actuel
                if len(combi)>0 or (allow_external_modes and self.acceptable(liste_nouveau_mode,constellation)):
                    nouveaux_modes.append(self.ajouterMode(constellation,liste_nouveau_mode))
                    if len(combi)==0:
                        external_modes = True
            return nouveaux_modes,external_modes,external_activites
        else:
            return [],False,len(contenu_actuel)>0
    
    def acceptable(self,liste_activites,constellation):
        
        couples = len(liste_activites)==0 or type(liste_activites[0])==tuple
        
        # vérifier qu'il y a du contenu
        if len(liste_activites)==0:
            return False
        # vérifier que les activités correspondent bien à la requête
        if couples:
            act = self.couples_satellites_activites
        else:
            act = [x[1] for x in self.couples_satellites_activites]
        if not sum([x in act for x in liste_activites])==len(liste_activites):
            return False
        # vérifier l'unicité des time slots
        plots = {}
        for X in liste_activites:
            if couples:
                (s,a) = X
            else:
                a = X
            ts = self.getTimeSlot(a)
            if ts in plots:
                return False # time slot doublement servi
            plots[ts] = [rso for rso in self.activitesParTimeSlots[ts] if rso[2]==a][0]
        # vérifier la positivité du score    
        if self.scoreListeActivites(plots,constellation)[0]<=0:
            return False

        return True
    
    def getTimeSlots(self):
        return {ts : [(x[1],x[2]) for x in self.activitesParTimeSlots[ts]] for ts in self.activitesParTimeSlots}
        
    def defilerActivites(self,explication):
        assert(explication!=[])
        for timeSlot in self.activites_candidates:
            self.activites_candidates[timeSlot] = [x for x in self.activites_candidates[timeSlot] if x[2] not in explication]
        self.maintenirObs(explication)
        
    def scoreListeActivites(self,plots,constellation):
        liste_score = [self.scoreObs(plots[timeSlot],constellation,lambda rso : self.scorePlot(constellation,rso,timeSlot)) for timeSlot in plots]
        #print(plots,[x[0] for x in liste_score ])
        Nplots = len(list(self.activitesParTimeSlots.keys()))
        max_gap = 0
        gap = 0
        for ts in sorted(self.activitesParTimeSlots):
            if len(plots.get(ts,[]))>0:
                if gap>max_gap:
                    max_gap = gap
                gap = 0
            else:
                gap += 1
        alpha = 0.5
        score_rec = sum([x[0] for x in liste_score ]) + alpha*(Nplots-max_gap)
        if len(liste_score)>0:
            if config.glob.score_temporel:
                return (score_rec,np.mean([x[1] for x in liste_score]))
            else:
                return (score_rec,0)
        else:
            return (0,0)
    
    def scoreListeObservations(self,listeObs,constellation,timeSlots):
        plots = {}
        for rso in listeObs:
            timeSlot = self.getTimeSlot(rso[2])
            plots[timeSlot] = rso
        return self.scoreListeActivites(plots,constellation)
    
    def chargerMaintient(self,obs):
        timeSlot_restant = list(self.activites_candidates.keys())
        listeAct = []
        for s in self.maintient:
            obs[s] = []
            for o in self.maintient[s]:
                obs[s].append(o)
                for timeSlot in self.activites_candidates:
                    match_obs = [x for x in self.activites_candidates[timeSlot] if x[2]==o]
                    if len(match_obs)==1:
                        timeSlot_restant.remove(timeSlot)
                        listeAct.append((timeSlot,match_obs[0]))
                        break
        return timeSlot_restant,listeAct

    def getBestModeWithoutInactives(self,constellation,inactifs):
        if not config.getOptValue("dynamic") or self.estActif():
            rso_liste = []
            obs = {}
            liste_obs = []
            for ts in self.activitesParTimeSlots:
                activites = [rso for rso in self.activitesParTimeSlots[ts] if rso[2] not in inactifs]
                cleTemporelle = lambda rso : self.scorePlot(constellation,rso,ts)
                if len(activites)>0:
                    rso = self.selectionnerObs(constellation,activites,cleTemporelle)
                    liste_obs.append(rso[2])
                    rso_liste.append(rso)
                    if rso[1] not in obs:
                        obs[rso[1]] = []
                    obs[rso[1]].append(rso[2])
            recompense,scoreTemporel = self.scoreListeObservations(rso_liste,constellation,cleTemporelle)
            if recompense<=0:
                return None
            mode =  Mode(self.idRequete,self.idMode,recompense,obs,scoreTemporel)
            if self.acceptable(mode.getCouples(),constellation):
                self.mode_a_valider = mode
                return self.mode_a_valider
            else:
                return None
        else:
            return None
                
    def construireMode(self,constellation,autoriser_degradation=True):
        # si degradation non autorisé et un slot n'a plus d'obs : plus de mode
        if not autoriser_degradation:
            for timeSlot in self.activites_candidates:
                if len(self.activites_candidates[timeSlot])==0:
                    self.mode_candidat = None
                    return
        obs = {}
        timeSlot_restant,listeAct = self.chargerMaintient(obs)
        plots = {x[0] : x[1] for x in listeAct}
        for timeSlot in timeSlot_restant:
            if len(self.activites_candidates[timeSlot])>0:
                cleTemporelle = lambda rso : self.scorePlot(constellation,rso,timeSlot)
                rso = self.selectionnerObs(constellation,self.activites_candidates[timeSlot],cleTemporelle,generateur=self.generateur_plot[timeSlot])
                listeAct.append(rso)
                plots[timeSlot] = rso
                s = rso[1]
                if s not in obs:
                    obs[s] = []
                obs[s].append(rso[2])
        # pas d'obs trouvées => mode nul        
        if len(list(obs.keys()))==0:
            self.mode_candidat = None
            return
        
        recompense,scoreTemporel = self.scoreListeActivites(plots,constellation)
        self.plots_candidat = plots
        self.mode_candidat = Mode(self.idRequete,self.idMode,recompense,obs,scoreTemporel)
        self.modes[self.idMode] = self.mode_candidat
        self.idMode += 1
        self.modes_generes += 1
    
    def scorePlot(self,constellation,rso,timeSlot):
        t = constellation.getSatellite(rso[1]).getObservation(rso[2]).getDebut()
        start = 30*60 + self.plots[timeSlot] # 30 minutes
        end = 90*60 + self.plots[timeSlot] # 1h30 : score nul
        if t<start:
            return 1
        else:
            return max(1 - (abs(t-start))/(end-start),0)
        
    def scoreTemporel(self,constellation,plots):
        if config.glob.score_temporel:
            return np.mean([self.scorePlot(constellation,plots[timeSlot],timeSlot) for timeSlot in plots])
        else:
            return 0
        
    def resetModes(self,constellation,conserve_old=True,initFirstMode=True):
        self.init = True
        self.activites_candidates = deepcopy(self.activitesParTimeSlots)
        self.resetCandidats()
        self.plots = {}
        debut = lambda rso: constellation.getSatellite(rso[1]).getObservation(rso[2]).getDebut()
        fin = lambda rso : constellation.getSatellite(rso[1]).getObservation(rso[2]).getFin()
        milieu = lambda rso : (debut(rso)+fin(rso))/2
        for timeSlot in self.activites_candidates:
            self.plots[timeSlot] = np.mean([milieu(rso) for rso in self.activites_candidates[timeSlot]])
        if not conserve_old:
            self.modes = {}
            self.idMode = 0          
        if initFirstMode:
            self.construireMode(constellation)
       
    def getKMeilleursModes(self,constellation,k):
        plots_reduits = {}
        for plot in self.plots:
            if len(self.activitesParTimeSlots[plot])<=k:
                plots_reduits[plot] = self.activitesParTimeSlots[plot]
            else: 
                plots_reduits[plot] = heapq.nlargest(k, self.activitesParTimeSlots[plot], key=itemgetter(0))
        combinaisons = [ {} ]
        liste_modes = []
        for plot in plots_reduits:
            combinaisons = [dict(chain.from_iterable(d.items() for d in (combi, {plot:x}))) for combi in combinaisons for x in plots_reduits[plot] ]
        if len(combinaisons)<=k:
            best_obs = combinaisons
        else:
            best_obs = heapq.nlargest(k,combinaisons,key=lambda x : self.scoreListeActivites(x,constellation))
        for comb in best_obs:
            obs = self.rsoToDict([comb[plot] for plot in comb])
            recompense,scoreTemporel = self.scoreListeActivites(comb,constellation)
            m = Mode(self.idRequete,self.idMode,recompense,obs,scoreTemporel)
            liste_modes.append(m)
            self.modes[self.idMode] = m
            self.idMode += 1
            self.modes_generes += 1
        return liste_modes
                                  
    def rsoToDict(self,rso_list):
        obs = {}
        for (r,s,o) in rso_list:
            if s not in obs:
                obs[s] = []
            obs[s].append(o)
        return obs
    
    def getModeSuivant(self,explication,constellation):
        assert(explication!=[])
        if self.mode_candidat is not None:
            self.defilerActivites(explication)
            self.construireMode(constellation)
        return self.mode_candidat

    def getModeSuivantSansExp(self,constellation):
        if self.mode_candidat is None:
            return None
        tmp = sorted([o for (s,o) in self.mode_candidat.getCouples()])
        id_mode = self.generateur_aleatoire_modes.randint(0,len(tmp)-1)
        exp = [tmp[id_mode]] # 1er couple (s,o) : obs
        return self.getModeSuivant(exp,constellation)
    
class RequeteSystematic(Requete):
    def __init__(self,idRequete,priorite,activites):
        super().__init__(idRequete,priorite,TYPE_SYSTEMATIC,activites)
        assert(len(self.couples_satellites_activites)==sum([len(self.liste_activites_par_satellite [s]) for s in self.liste_activites_par_satellite]))

    def genererModesPresentsCCA(self,grapheDep,constellation,contenu_actuel,liste_cca_a_explorer,allow_external_modes=True):
        if not config.getOptValue("dynamic") or self.estActif():
            external_modes = False
            external_activites = len(contenu_actuel)>0
            if len(contenu_actuel)>0:
                sat = constellation.getSatelliteActivite(contenu_actuel[0])
                if config.getOptValue("verif"):
                    constellation.getSatelliteActivite(contenu_actuel[0]) == sat
                possible_sat = [sat]
            else:
                possible_sat = list(self.liste_activites_par_satellite.keys())
            nouveaux_modes = []
            activites_cca = self.getActivitesPresentesCCA(grapheDep, liste_cca_a_explorer)
            
            for s in possible_sat:
                activites_candidates = [x for x in self.liste_activites_par_satellite [s] if x in activites_cca]
                combinaisons = []
                for taille in range(len(activites_candidates)+1):
                    combinaisons += list(map(list,itertools.combinations(activites_candidates,taille)))
                for combi in combinaisons:
                    liste_nouveau_mode = combi+contenu_actuel
                    if len(combi)>0 or (allow_external_modes and self.acceptable(liste_nouveau_mode,constellation)):
                        nouveaux_modes.append(self.ajouterMode(constellation,liste_nouveau_mode))
                    if len(combi)==0 and self.acceptable(liste_nouveau_mode,constellation):
                        external_modes = True
            return nouveaux_modes,external_modes,external_activites
        else:
            return [],False,len(contenu_actuel)>0
        
    def acceptable(self,liste_activites,constellation):
        if len(liste_activites)==0:
            return False
        couples = len(liste_activites)==0 or type(liste_activites[0])==tuple
        if couples:        
            act = self.couples_satellites_activites
        else:
            act = [x[1] for x in self.couples_satellites_activites]
                
        if not sum([x in act for x in liste_activites])==len(liste_activites):
            return False
        if couples:
            sat = liste_activites[0][0]
        else:
            sat = constellation.getSatelliteActivite(liste_activites[0])
        if couples:
            for s,a in liste_activites:
                if sat!=s:
                    return False
        else:
            for a in liste_activites:
                s = constellation.getSatelliteActivite(a)
                if sat!=s:
                    return False
        return True
            
    def resetModes(self,constellation,conserve_old=True,initFirstMode=True):
        self.init = True
        self.observations_candidates = {}
        for (r,s,o) in self.activitesParTimeSlots[0]:
            if s not in self.observations_candidates:
                self.observations_candidates[s] = {}
            self.observations_candidates[s][o] = r
        self.resetCandidats()
        self.dureeMin = min([constellation.getSatellite(rso[1]).getObservation(rso[2]).getDuree() for rso in self.activitesParTimeSlots[0]]) 
        self.dureeMax = max([constellation.getSatellite(rso[1]).getObservation(rso[2]).getDuree() for rso in self.activitesParTimeSlots[0]])
        if self.dureeMax==self.dureeMin:
            self.delta = 1
        else:
            self.delta = self.dureeMax - self.dureeMin
        if not conserve_old:
            self.modes = {}
            self.idMode = 0          
        if initFirstMode:
            self.construireMode(constellation)
        #print("reset")
        #print(self.idRequete,self.modes[0].getCouples(),[len(self.observations_candidates[t]) for t in self.observations_candidates])
    
    def getKMeilleursModes(self,constellation,k):
        min_longueur = {s:max(1,len(self.liste_activites_par_satellite [s])-k) for s in self.liste_activites_par_satellite}
        comb = {s: [list(l) for n in range(min_longueur[s],len(self.liste_activites_par_satellite [s])+1) for l in combinations([rso for rso in self.activitesParTimeSlots[0] if rso[1]==s],n)] for s in self.liste_activites_par_satellite}
        #print(comb)
        flatten = []
        for s in comb:
            for x in comb[s]:
                elmt = (x,self.scoreListeObservations(x,constellation))
                flatten.append(elmt)
        #print(flatten)
        liste_obs = [x[0] for x in heapq.nlargest(k, flatten, key=itemgetter(1))]
        liste_modes = []
        for l in liste_obs:
            obs = {}
            for rso in l:
                if rso[1] not in obs:
                    obs[rso[1]] = []
                obs[rso[1]].append(rso[2])
            recompense,score_temporel = self.scoreListeObservations(l,constellation)
            m = Mode(self.idRequete,self.idMode,recompense,obs,score_temporel)
            self.modes[self.idMode] = m
            liste_modes.append(m)
            self.idMode += 1
            self.modes_generes += 1
        return liste_modes
                               
    def scoreListeObservations(self,listeObs,constellation,timeSlots=None):
        cleTemporelle = lambda rso : (constellation.getSatellite(rso[1]).getObservation(rso[2]).getDuree() - self.dureeMin)/self.delta
        liste_score = [self.scoreObs(rso,constellation,cleTemporelle) for rso in listeObs]
        return (sum([x[0] for x in liste_score]),np.mean([x[1] for x in liste_score]))
    
    def construireMode(self,constellation):
        if self.observations_candidates=={}:
            s_max = None
        else:
            r_max = max([self.scoreListeObservations([(self.observations_candidates[s][o],s,o) for o in self.observations_candidates[s]],constellation) for s in self.observations_candidates])
            #if self.mode_candidat is not None and r_max==self.mode_candidat.getRecompense():
            #    return # on s'embete pas a changer de satellite si les autres ne sont pas strictement mieux
            liste_sat_max = [s for s in self.observations_candidates if self.scoreListeObservations([(self.observations_candidates[s][o],s,o) for o in self.observations_candidates[s]],constellation) == r_max]
            #id_sat = self.generateur_aleatoire_modes.randint(0,len(liste_sat_max)-1)
            s_max = min(liste_sat_max) # liste_sat_max[id_sat]

        if s_max is not None:
            recompense,score_temporel = r_max
            obs = {s_max:list(self.observations_candidates[s_max].keys())}
            self.mode_candidat = Mode(self.idRequete,self.idMode,recompense,obs,score_temporel)
            self.modes[self.mode_candidat.getId()] = self.mode_candidat
            self.idMode += 1
            self.modes_generes += 1
        else:
            self.mode_candidat = None
            
    def getObsBySatellites(self,inactifs):
        obs_by_sat = {}
        for ts in self.activitesParTimeSlots:
            for (r,s,a) in self.activitesParTimeSlots[ts]:
                if a not in inactifs:
                    if s not in obs_by_sat:
                        obs_by_sat[s] = []
                    obs_by_sat[s].append((r,s,a))
        return obs_by_sat
        
    def getBestModeWithoutInactives(self,constellation,inactifs):
        if not config.getOptValue("dynamic") or self.estActif():
            cleTemporelle = lambda rso : self.scoreTemporel(constellation,rso)
            rso_liste = []
            obs = {}
            liste_obs = []
            rso_by_sat = self.getObsBySatellites(inactifs)
            if rso_by_sat=={}:
                self.mode_a_valider = None
                return None
            best_sat = None
            score_sat = 0
            for sat in rso_by_sat:
                recompense,scoreTemporel = self.scoreListeObservations(rso_by_sat[sat],constellation,cleTemporelle)
                if recompense>score_sat:
                    score_sat = recompense
                    best_sat = sat
            
            obs = {best_sat:[rso[2] for rso in rso_by_sat[best_sat]]}
            mode =  Mode(self.idRequete,self.idMode,recompense,obs,scoreTemporel)            
            if config.getOptValue("verif"):
                assert(self.acceptable(mode.getCouples(),constellation))
            self.mode_a_valider = mode
            return self.mode_a_valider
        else:
            return None
        
    def defilerActivites(self,explication):
        for a in explication:
            for s in self.observations_candidates:
                if a in self.observations_candidates[s]:
                    del self.observations_candidates[s][a]
        
        sat = list(self.observations_candidates.keys())
        for s in sat:
            if self.observations_candidates[s] == {}:
                del self.observations_candidates[s]

    def getModeSuivantSansExp(self,constellation):
        if self.mode_candidat is None:
            return None
        tmp = sorted([o for (s,o) in self.mode_candidat.getCouples()])
        id_mode = self.generateur_aleatoire_modes.randint(0,len(tmp)-1)
        exp = [tmp[id_mode]] # 1er couple (s,o) : obs
        return self.getModeSuivant(exp,constellation)
    
    def toutesExplicationsEnRetard(self,explication):
        return sum([a not in [x[1] for x in self.mode_candidat.getCouples()] for a in explication]) == len(explication)
        
    def getModeSuivant(self,explication,constellation):
        assert(explication!=[])
        if self.mode_candidat is not None:
            en_retard = self.toutesExplicationsEnRetard(explication)
            self.defilerActivites(explication)
            if not en_retard: # si l'explication provient d'un ancien mode inutile de changer encore de satellite
                self.construireMode(constellation)
        return self.mode_candidat
    
class RequeteVidage(Requete):
    def __init__(self,idRequete,priorite,vidages):
        self.liste_activites_par_satellite = vidages
        self.couples_satellites_activites = []
        rso = []
        for s in vidages:
            for v in vidages[s]:
                self.couples_satellites_activites.append((s,v))
                rso.append((0,s,v))
        assert(len(self.couples_satellites_activites)==sum([len(self.liste_activites_par_satellite [s]) for s in self.liste_activites_par_satellite]))
        super().__init__(idRequete,priorite,TYPE_VIDAGE,{0:rso})
        self.type_requete = TYPE_VIDAGE
        self.mode_candidat = Mode(idRequete,0,0,self.liste_activites_par_satellite,0)
        self.modes = {0:self.mode_candidat}
        self.idMode = 0
    
    def acceptable(self,contenu,constellation):
        for (s,a) in self.couples_satellites_activites:
            if (s,a) not in contenu:
                return False
        return True
        
    def genererModesPresentsCCA(self,grapheDep,constellation,contenu_actuel,liste_cca_a_explorer,allow_external_modes=True):
        assert(not config.getOptValue("dynamic") or self.estActif())
        if config.getOptValue("verif"):
            activites_cca = self.getActivitesPresentesCCA(grapheDep, liste_cca_a_explorer)
            liste_nouveau_mode = activites_cca+contenu_actuel
            assert(sorted(liste_nouveau_mode) == sorted([x[1] for x in self.modes[0].getCouples()]))
        return [self.modes[0]],False,len(contenu_actuel)>0
        
    def scoreListeObservations(self,listeObs,constellation,timeSlots):
        return (0,0)
    
    def getBestModeWithoutInactives(self,constellation,inactifs):
        assert(self.estActif())
        if len(inactifs)>0:
            print([constellation.getRequete(a) for a in inactifs])
        assert(len(inactifs)==0)
        self.mode_a_valider = self.mode_candidat
        return self.mode_candidat

    def getModeSuivantSansExp(self,constellation):
        if self.mode_candidat is None:
            return None
        exp = [] # 1er couple (s,o) : obs
        res =  self.getModeSuivant(exp,constellation)
        self.mode_candidat = res
        return res
    
    def getKMeilleursModes(self,constellation,k):
        exp = [] # 1er couple (s,o) : obs
        res =  self.getModeSuivant(exp,constellation)
        self.mode_candidat = res
        return [res]
    
    def resetModes(self,constellation,conserve_old=True,initFirstMode=True):
        self.init = True
        pass
                  
    def getModeSuivant(self,explication,constellation):
        mode = deepcopy(self.mode_candidat)
        self.idMode += 1
        self.modes_generes += 1
        mode.idMode = self.idMode
        retrait = []
        for (s,o) in mode.couples:
            if o in explication:
                retrait.append((s,o))
        for (s,o) in retrait:
            mode.couples.remove((s,o))
            del mode.observations[s][o]
        self.mode_candidat = mode
        self.modes[1] = mode
        return mode
    
class Activite:
        def __init__(self,id_act,satid,debut,fin,longitude,latitude,altitude):
            self.id = id_act
            self.satId = satid
            self.debut = debut
            self.fin = fin
            self.longitude = longitude
            self.latitude = latitude
            self.altitude = altitude

        def getId(self):
            return self.id

        def getDebut(self):
            return self.debut

        def getFin(self):
            return self.fin

        def getSat(self):
            return self.satId

        def getLongitude(self):
            return self.longitude

        def getLatitude(self):
            return self.latitude

        def getAltitude(self):
            return self.altitude

        def getCoordonnees(self):
            return ( self.longitude,self.latitude,self.altitude)

        def estDependant(self,constellation,s,a):
            if s!=self.satId:
                return False
            autre = constellation.satellites[self.satId].getActivite(a)
            left = autre.getDebut() <= self.getDebut() and autre.getFin() + constellation.satellites[s].getTransition(self.id,a) <= self.getDebut()
            right = self.getDebut() <= autre.getDebut() and self.getFin() + constellation.satellites[s].getTransition(a,self.id) <= autre.getDebut()
            return not(left or right)
        
        def distance(self,s,a2,constellation):
            coords2 = constellation.getSatellite(s).getActivite(a2).getCoordonnees()
            return sqrt(sum([(self.getCoordonnes()[i]-coords2[i])**2 for i in range(3)]))
        
        def dependantTauMax(self,s,a,constellation,tau_max):
            autre = constellation.satellites[s].getActivite(a)
            left = autre.getDebut() <= self.getDebut() and autre.getFin() + tau_max <= self.getDebut()
            right = self.getDebut() <= autre.getDebut() and self.getFin() + tau_max <= autre.getDebut()
            return not(left or right)
        
class Observation(Activite):
        def __init__(self,idObs,idSat,debut,fin,dureeObs,x,y,z,score,idRequete):
            super().__init__(idObs,idSat,debut,fin,x,y,z)
            self.duree = dureeObs
            #self.dureeTelechargement = dureeDl
            #self.memoire = memoire
            self.scoringObs(score)
            #printMaster(idObs,score,self.score)
            self.idRequete = idRequete
        
        def scoringObs(self,score):
            if score > 1:
                raise ValueError("Score de l'observation > 1",score)
            #bounds = [x/10 for x in range(5,11)]
            values = [(0.5,0.05),(0.6,0.1),(0.7,0.4),(0.8,0.7),(0.9,0.95),(1,1)]
            values = sorted(values,key=itemgetter(0))
            for (bound,value) in values:
                if score <= bound:
                    self.score = value
                    return
            raise ValueError("Score de l'observation > 1",score)
            
        def getRequete(self):
            return self.idRequete
            
        def getScore(self):
            return self.score
        
        def getScorePlan(self,composantes,solCCAs):
            s = self.getSat()
            cca = composante.getActiviteCCA(self.getIdentifiant())
            charge = solCCAs[s][cca].scoreCharge()
            return (self.score,-charge)
            
        def getDuree(self):
            return self.duree

        def getCoordonnees(self):
            return (self.getLongitude(),self.getLatitude(),self.getAltitude())
        
        def __str__(self):
            return " " + Fore.BLUE + "[Observation, début : " + str(self.debut) + ", id : " + str(self.id) + ", fin : " + str(self.fin) + ", durée : " + str(self.duree) + ", score : " +str(self.score) +"]" +Style.RESET_ALL

class Vidage(Activite):
        # id,satId,windowStart,windowEnd,longitude,latitude,altitude
        def __init__(self,data):
            super().__init__(int(data[0]),int(data[1]),float(data[2]),float(data[3]), float(data[4]),float(data[5]),float(data[6]))
            #super().__init__(int(data[0])+config.donnees.no,int(data[1]),float(data[2]),float(data[3]), float(data[4]),float(data[5]),float(data[6]))
            
        def __str__(self):
            return " " + Fore.GREEN + "[Vidage, début : "+str(self.debut)+", id : "+str(self.id)+", fin : "+str(self.fin) + ", durée : " + str(self.getDuree()) +" ]" + Style.RESET_ALL
        
        def getDuree(self):
            fenetre = self.fin - self.debut
            return min(config.glob.allocation_vidage,0.95*fenetre)
        
        def getRequete(self):
            return -1
        
        def getScore(self):
            return 0
            
class Satellite:
        def __init__(self,identifiant,obs,vid):
            self.obs = obs
            self.id = identifiant
            self.vid = vid

        def supprimerActivite(self,a):
            if self.estVidage(a):
                del self.vid[a]
            if self.estObservation(a):
                del self.obs[p]

        def estVidage(self,p):
            #return p>=config.donnees.no and p<config.donnees.no+config.donnees.nv
            return p in self.vid    
        
        def estObservation(self,p):
            #return p>=0 and p<config.donnees.no
            return p in self.obs

        def getObservations(self):
            return self.obs

        def getObservation(self,p):
            return self.obs[p]

        def getId(self):
            return self.id

        def getVidage(self,p):
            return self.vid[p]

        def getActivites(self):
            r = dict(self.obs)
            r.update(self.vid)
            return r

        def getActivite(self,a):
            if self.estObservation(a):
                return self.obs[a]
            else:
                return self.vid[a]

        def __str__(self):
            s = "IDENTIFIANT : " + str(self.id) + '\n'
            s += "ACTIVITES : \n"
            for o in self.obs :
                s += "\t" + str(self.obs[o]) + "\n"
            return s

        def getVidages(self):
            return self.vid
        
        def getTransition(self,a1,a2,modeleDeTransition):
            assert(not modeleDeTransition.estTimeDependent())
            res = modeleDeTransition.getTransition(self,a1,a2,config.glob.digits)
            return res
        
        def getTransitionTimeDependent(self,a1,a2,start_maneuvering,modeleDeTransition):
            if modeleDeTransition.estTimeDependent():
                return modeleDeTransition.getTransition(a1,a2,start_maneuvering)
            else:
                return modeleDeTransition.getTransition(self,a1,a2)

class Constellation:

    """
        =========================================================
                Initialisation de la constellation
        =========================================================
    """            
    def getRequeteActivite(self,a):
        s = self.mapping_sat[a]
        return self.dependances[s][a]
    
    def initDependancesActivites(self):
        # mapping dep[s][a] => (requete,mode) qui contiennent a
        self.dependances = {s : {} for s in self.satellites}
        for r in self.getToutesRequetes():
            for s in self.getRequete(r).getActivites():
                if s not in self.satellites:
                    die(s,self.satellites.keys(),self.getRequete(r).getActivites().keys())
                for a in  self.getRequete(r).getActivitesSatellite(s):
                    self.dependances[s][a] = r
                
    def initTimeData(self):
        self.tmin = min([self.satellites[s].getActivite(p).getDebut() for s in self.satellites for p in self.satellites[s].getActivites()])
        self.tmax = max([self.satellites[s].getActivite(p).getFin() for s in self.satellites for p in self.satellites[s].getActivites()])
        mins = {s : min([self.satellites[s].getObservation(p).getDebut() for p in self.satellites[s].getObservations()]) for s in self.satellites}
        maxs = {s : max([self.satellites[s].getObservation(p).getFin() for p in self.satellites[s].getObservations()]) for s in self.satellites}
        self.horizonObs = { s : (mins[s],maxs[s]) for s in mins}
        mins = {s : min([self.satellites[s].getVidage(p).getDebut() for p in self.satellites[s].getVidages()]) for s in self.satellites}
        maxs = {s : max([self.satellites[s].getVidage(p).getFin() for p in self.satellites[s].getVidages()]) for s in self.satellites}
        self.horizonVid = { s : (mins[s],maxs[s]) for s in mins} 
    
    def initGlobalData(self):
        config.donnees.ns = self.ns
        config.donnees.no = self.no
        config.donnees.nv = self.nv
        config.donnees.tmin = self.tmin
        config.donnees.tmax = self.tmax    
    
    def readObs(self,line,offset,idRequete):
        if not(len(line)==9+offset):
            print(line,len(line),offset)
            assert(False)
        idObs,idSat = int(line[0+offset]),int(line[1+offset])
        debut,fin = float(line[2+offset]),float(line[3+offset])
        #dureeObs,dureeDl = float(line[4+offset]),float(line[5+offset])
        #memoire = int(line[6+offset])
        dureeObs = float(line[4+offset])
        x,y,z = float(line[5+offset]),float(line[6+offset]),float(line[7+offset])
        score = float(line[8+offset])
        #downloads = [int(x) for x in line[9+offset:]]
        self.mapping_sat[idObs] = idSat
        #assert(score>0)
        return Observation(idObs,idSat,debut,fin,dureeObs,x,y,z,score,idRequete)
    
    def supprimerSatellitesSansObs(self,vid,obs):
        del_vid = [s for s in vid if s not in obs]
        for s in del_vid:
            del vid[s]
    
    def statsCoordonnees(self):
        xyz = [0,0,0]
        minCoord = [np.Inf,np.Inf,np.Inf]
        maxCoord = [-np.Inf,-np.Inf,-np.Inf]
        nPoints = 0
        for r in self.getToutesRequetes():
            xyz_r,minCoord_r,maxCoord_r,n_r = self.getRequete(r).statsCoordonnees(self)
            xyz = [xyz[i]+xyz_r[i] for i in range(3)]
            minCoord = [min(minCoord[i],minCoord_r[i]) for i in range(3)]
            maxCoord = [max(maxCoord[i],maxCoord_r[i]) for i in range(3)]
            nPoints += n_r
        xyz = tuple([xyz[i]/nPoints for i in range(3)])
        return xyz,minCoord,maxCoord,nPoints
    
    def readDownloads(self,fichier,affichage):
        self.nv = int(fichier.readline())
        vid = {}
        for i in range(self.nv):
            line = fichier.readline().split(",")
            v = Vidage(line)
            s = v.getSat()
            if s not in vid:
                vid[s] = {}
            vid[s][v.getId()] = v
            self.mapping_sat[v.getId()] = s
        return vid
    
    def readSatellites(self,obs,vid,affichage):
        vehicules = {}
        for s in vid:
            if s in obs:
                vehicules[s] = Satellite(s,obs[s],vid[s])
        self.ns = len(list(vehicules.keys()))
        self.satellites = vehicules        
                    
    def readRequests(self,fichier):
        line = fichier.readline().split('\n')[0].split(",")
        nRequetes = int(line[0])
        observations_satellites = {}
        no = 0
        self.requetes = {}
        err = 0
        for r in range(nRequetes):
            line = fichier.readline().split('\n')[0].split(",")
            idRequete,nObs,priorite,typeRequete = int(line[0]),int(line[1]),0,line[2]
            no += nObs
            observations = {}
            for i in range(nObs):
                line = fichier.readline().split('\n')[0].split(",")
                if typeRequete == TYPE_PERIODIC or typeRequete == TYPE_STEREO:
                    timeSlot = int(line[0])
                    offset = 1
                else:
                    timeSlot = 0
                    offset = 0
                obs = self.readObs(line,offset,idRequete)
                if timeSlot not in observations:
                    observations[timeSlot] = []
                if (obs.getId() not in [x[2] for x in observations[timeSlot]]):
                    observations[timeSlot].append((obs.getScore(),obs.getSat(),obs.getId()))
                    s = obs.getSat()
                    if s not in observations_satellites:
                        observations_satellites[s] = {}
                    observations_satellites[s][obs.getId()] = obs
                else:
                    err += 1
            if r<config.getOptValue("test_req") and sum([len(observations[s])!=0 for s in observations]):
                if typeRequete==TYPE_MONO:
                    req = RequeteMono(idRequete,priorite,observations)
                elif typeRequete==TYPE_LONG_MONO:
                    req = RequeteLongMono(idRequete,priorite,observations)
                elif typeRequete==TYPE_STEREO:
                    req = RequeteStereo(idRequete,priorite,observations)
                elif typeRequete==TYPE_PERIODIC:
                    req = RequetePeriodic(idRequete,priorite,observations)
                elif typeRequete==TYPE_SYSTEMATIC:
                    if config.getOptValue("include_systematic"):
                        req = RequeteSystematic(idRequete,priorite,observations)
                    else:
                        req = None
                else:
                    print(line)
                    raise ValueError("type de requête inconnue : "+typeRequete)
                if req is not None:
                    self.requetes[idRequete] = req
        #printColor(err,"obs redondantes",c='r')
        self.no = no
        config.donnees.no = no
        return observations_satellites
    
    def filtrerModesPresents(self,solution,grapheDep,modeleDeTransition):
        contenu_modes = {}
        modes = []
        
        for s in solution.getSolCCAs():
            for cca in solution.getSolCCAs()[s]:
                for a in solution.getSolCCAs()[s][cca].getSequence():
                    r = self.getRequeteActivite(a)
                    if r not in contenu_modes:
                        contenu_modes[r] = []
                    contenu_modes[r].append(a)
        for r in self.getRequetes():
            if r!= -1: # vidage
                self.getRequete(r).resetModes(self,conserve_old=False,initFirstMode=False)
                assert(len(self.getRequete(r).getModes())==0)
                if r in contenu_modes:
                    if self.getRequete(r).acceptable(contenu_modes[r],self):
                        mode = self.getRequete(r).filtrerModesPresents(contenu_modes[r],self)
                        assert(mode.getId()==0)
                        modes.append((r,mode.getId()))
                    else:        
                        for a in contenu_modes[r]:
                            (s,cca) = grapheDep.getActiviteCCA(a)
                            solution.getSolCCA(s,cca).retirerActivite(self,a,modeleDeTransition)
            else:
                modes.append((r,0))
        solution.setModesRetenus(modes,self)
        
    def __init__(self,fileName,start_date,affichage=True):
        config.glob.filename = fileName.split("/")[-1].split(".")[0]
        comm = MPI.COMM_WORLD
        size = comm.Get_size()
        config.glob.size = size
        
        self.start_date = start_date
        self.fichier = fileName
        fichier = open(fileName,"r")
        self.mapping_sat = {}
        obs = self.readRequests(fichier)
        vid = self.readDownloads(fichier,affichage)
        self.supprimerSatellitesSansObs(vid,obs)
        self.requetes[-1] = RequeteVidage(-1,1,vid)
        self.readSatellites(obs,vid,affichage)
        self.initTimeData()
        self.initDependancesActivites()
        self.initGlobalData()
        
        
        if config.getOptValue("verif"):
            for r in self.getRequetes():
                activites = self.extraireActivitesRequetes()
                for s in activites:
                    for a in activites[s]:
                        assert(s==self.getSatelliteActivite(a))
        
        if affichage:
            printColor("Lecture OK.\n",c='g')
            xyz,minCoord,maxCoord,nPoints = self.statsCoordonnees()
            shift = getDisplayDepth()-1
            shiftLeftDisplay(shift)
            printColor("====================== [ Constellation ] ======================",c='y')
            printColor("| Requêtes :",str(len(list(self.requetes.keys()))),c='y')
            printColor("| Observations :",str(nPoints),c='y')
            printColor("| Interval latitude : [",minCoord[0],",",maxCoord[0],"]",c='y')
            printColor("| Interval longitude : [",minCoord[1],",",maxCoord[1],"]",c='y')
            printColor("| Point moyen :",(xyz[0],xyz[1]),c='y')
            printColor("===============================================================",c='y')
            shiftRightDisplay(shift)
        if config.getOptValue("dynamic"):
            minute_dernieres_requetes=config.getOptValue("derniere_requete")
            proportion_depart = config.getOptValue("proportion_depart")
            discretisation_minutes = config.getOptValue("dt_requetes")
            self.requestDealer = RequestDealer(self,start_date,proportion_depart,discretisation_minutes*60,minute_dernieres_requetes*60)
            
    
    """
        =========================================================
                        GETTER
        =========================================================    
    """
    def getRequetesDepart(self):
        return self.requestDealer.getRequetesDepart()
    
    def getFilename(self):
        return config.glob.filename
    
    def getSatellite(self,s):
        return self.satellites[s]
    
    def getActivite(self,a,s=None):
        if s is None:
            s = self.getSatelliteActivite(a)
        return self.getSatellite(s).getActivite(a)
    
    def getNombreActivites(self):
        return self.no+self.nv

    def getIdActivitesMax(self):
        return max(list(self.mapping_sat.keys()))
    
    def getSatellites(self):
        return self.satellites
    
    def getToutesRequetes(self):
        return list(self.requetes.keys())
    
    def libererNouvellesRequetes(self,grapheDependances):
        return self.requestDealer.scruterArriveeNouvellesRequetes(self,grapheDependances)
    
    def getRequetes(self,allow_new=False,grapheDependances=None):  
        if config.getOptValue("dynamic"):
            if allow_new:
                self.requestDealer.scruterArriveeNouvellesRequetes(self,grapheDependances)
            else:
                assert(grapheDependances is None)
            return [r for  r in self.requetes if self.requetes[r].estActif()]
        else:
            assert(not allow_new)
            return list(self.requetes.keys())
    
    def getRequete(self,r):
        return self.requetes[r]
    
    def estVidage(self,p):
        return p>=self.no
    
    def estObservation(self,p):
        return p<self.no
    
    def getSatelliteActivite(self,a):
        return self.mapping_sat[a]
    
    def nombresModes(self):
        modes_easy = [len(list(self.requetes[r].getModes().keys())) for r in self.requetes if self.getRequete(r).getType()=="ONE_SHOT_MONO" or self.getRequete(r).getType()=="ONE_SHOT_STEREO"]
        modes_hard = [len(list(self.requetes[r].getModes().keys())) for r in self.requetes if self.getRequete(r).getType()=="PERIODIC" or self.getRequete(r).getType()=="SYSTEMATIC"]
        return np.mean(modes_easy),np.std(modes_easy),np.mean(modes_hard),np.std(modes_hard)
    
    def taillesModes(self):
        modes_easy = [len(self.requetes[r].getMode(m).getCouples()) for r in self.requetes if self.getRequete(r).getType()=="ONE_SHOT_MONO" or self.getRequete(r).getType()=="ONE_SHOT_STEREO"]
        modes_hard = [len(self.requetes[r].getMode(m).getCouples()) for r in self.requetes if self.getRequete(r).getType()=="PERIODIC" or self.getRequete(r).getType()=="SYSTEMATIC"]
        return np.mean(modes_easy),np.std(modes_easy),np.mean(modes_hard),np.std(modes_hard)

    def melangeBatch(self,liste,quantiles):
        longueur = [int(q*len(liste)) for q in quantiles]
        for i,q in enumerate(longueur):
            if(i<len(longueur)-1):
                q2 = longueur[i+1]
                liste[q:q2] = rd.sample(liste[q:q2], len(liste[q:q2]))

    def scoreObsMoyenne(self,liste_modes_candidats):
        score = []
        for (r,m) in liste_modes_candidats:
            for  (s,o) in self.getRequete(r).getCouples():
                if self.getSatellite(s).estObservation(o):
                    score.append(self.getSatellite(s).getObservation(o).getScore())
        if len(score)==0:
            return 0
        else:
            return np.mean(score)
    
    def afficherModes(self):
        for r in self.getRequetes():
            for m in self.getRequete(r).getModes():
                print(self.getRequete(r).getMode(m))
    
    def extraireActivitesMode(self,constellation,r,m):
        act = []
        for (s,o,d) in constellation.getRequete(r).getMode(m).getCouples():
            if o not in act:
                act.append(o)
            if d not in act:
                act.append(d)
        return act  
    
    def extraireActivitesRequetes(self):
        act = {}
        for r in self.getToutesRequetes():
            act_r = self.getRequete(r).getActivites()
            for s in act_r:
                if s not in act:
                    act[s] = []
                for o in act_r[s]:
                    act[s].append(o)
        return act 

    def extraireActivitesRequetesActives(self,grapheDependances):
        act = {}
        for r in self.getRequetes(grapheDependances):
            act_r = self.getRequete(r).getActivites()
            for s in act_r:
                if s not in act:
                    act[s] = []
                for o in act_r[s]:
                    act[s].append(o)
        return act
    """
        =========================================================  
                    VERIF/DEBUG
        ========================================================= 
    """
    def verifierModes(self):
        n_modes = sum([len(list(self.getRequete(r).getModes().keys())) for r in self.getToutesRequetes()])
        modes_faux = 0
        for r in self.getToutesRequetes():
            for m in self.getRequete(r).getModes():
                for s,o,d in self.getRequete(r).getMode(m).getCouples():
                    if(self.satellites[s].getObservation(o).getDuree()>self.satellites[s].getObservation(o).getFin()-self.satellites[s].getObservation(o).getDebut()):
                        modes_faux += 1
                        break
        return modes_faux/n_modes
    
    
    """
        =========================================================  
                        GRAPHIQUES
        =========================================================  
    """
    def tracerPertesModes(self):
        return self.solution.tracerPertesModes()
        #plt.title(str(len(self.solution.modes_retenus))+'/'+str(len(list(self.requetes.keys())))+" fulfilled requests")
    
    def tracerNombresModes(self,affichage='boxplot'):
        if(affichage=='hist'):
            f,ax = plt.subplots()
            valeurs = np.unique([self.getRequete(r).getType() for r in self.getToutesRequetes()])
            for val in valeurs:
                ax.hist([len(self.getRequete(r).getModes()) for r in self.getToutesRequetes() if self.getRequete(r).getType()==val])
            plt.legend(valeurs)
            plt.title('count modes')
            #return f
        elif affichage=='boxplot':       
            df = pd.DataFrame([],columns=['count','type'])
            for r in self.getToutesRequetes():
                df = df.append({'count':len(self.getRequete(r).getModes()),'type':self.getRequete(r).getType()},ignore_index=True)
            sns.boxplot(y='count',x='type',data=df, palette="colorblind").set_title('count modes')
        else:
            raise NameError()
            
    def tracerTaillesModes(self,affichage='boxplot'):
        if(affichage=='hist'):
            f,ax = plt.subplots()
            valeurs = np.unique([self.getRequete(r).getType() for r in self.getToutesRequetes()])
            for val in valeurs:
                ax.hist([len(self.getRequete(r).getMode(m).getCouples()) for r in self.getToutesRequetes() for m in self.getRequete(r).getModes() if self.getRequete(r).getType()==val])
            plt.legend(valeurs)
            plt.title("mode size (number of observations)")
            #return f
        elif affichage=='boxplot':
            df = pd.DataFrame([],columns=['count','type'])
            for r in self.getToutesRequetes():
                for m in self.getRequete(r).getModes():
                    df = df.append({'count':len(self.getRequete(r).getMode(m).getCouples()),'type':self.getRequete(r).getType()},ignore_index=True)
            sns.boxplot(y='count',x='type',data=df, palette="colorblind").set_title('modes size')
        else:
            raise NameError()  
            
    def tracerDemandeConstellation(self):
        f,ax = plt.subplots(len(self.satellites.keys()),1)
        obs = {s : self.satellites[s].getObservations() for s in self.satellites}
        for s in obs:
            X = []
            Y = []
            temps = []
            for o in obs[s]:
                if obs[s][o].getDebut() not in temps:
                    temps.append(obs[s][o].getDebut())
                if obs[s][o].getFin() not in temps:
                    temps.append(obs[s][o].getFin())
            for t in temps:
                count = len([o for o in obs[s] if t>= obs[s][o].getDebut() and t<= obs[s][o].getFin()])
                X.append(t)
                Y.append(count)
            XY = [(X[i],Y[i]) for i in range(len(X))]
            XY = sorted(XY,key = itemgetter(0))
            ax[s].plot([x[0] for x in XY],[x[1] for x in XY],color='black')
        f.suptitle('charge de la constellation')
        return f
    
    def tracerDemande(self,s):
        f,ax = plt.subplots(1,1)
        obs = {s : self.satellites[s].getObservations()}
        X = []
        Y = []
        temps = []
        for o in obs[s]:
            if obs[s][o].getDebut() not in temps:
                temps.append(obs[s][o].getDebut())
            if obs[s][o].getFin() not in temps:
                temps.append(obs[s][o].getFin())
        for t in temps:
            count = len([o for o in obs[s] if t>= obs[s][o].getDebut() and t<= obs[s][o].getFin()])
            X.append(t)
            Y.append(count)
        XY = [(X[i],Y[i]) for i in range(len(X))]
        XY = sorted(XY,key = itemgetter(0))
        ax.plot([x[0] for x in XY],[x[1] for x in XY],color='black')
        f.suptitle('charge du satellite {}'.format(s))
        return f
        
    def tracerMode(self,f,r,m):
        #f,ax = plt.subplots()
        for s,o,d in self.requetes[r].getMode(m).getCouples():
            debut,fin = self.satellites[s].getObservations()[o].getDebut(),self.satellites[s].getObservations()[o].getFin()
            plt.plot([debut,fin],[s,s],'-b',alpha=0.1)
        return f
    
    def tracerModes(self,f):
        for r in self.requetes:
            for m in self.requetes[r].getModes():
                f = self.tracerMode(f,r,m)
        return f
    
    def tracerEchantillonsRecompenses(self):
        rr = rd.choices(list(self.getToutesRequetes().keys()),k=4)
        for r in rr:
            n,bins,patch = plt.hist([constellation.getRequete(r).getMode(m).getRecompense()  for m in constellation.getRequete(r).getModes()])
            print(histToLatex(n,bins,"reward (request "+str(r)+")","mode count"))
            
    def tracerActivite(self,annoter=False):
        f = self.solution.tracerActivite(annoter)
        #plt.title(str(len(self.getModesSolution())) +'/'+ str(len(self.getRequetes().keys())) + ' requêtes satisfaites')
        return f
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

from .graph import *

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
global TYPE_DOWNLOAD
TYPE_DOWNLOAD = "DOWNLOAD"



"""
        =========================================================
                    Attributs de la constellation 
        =========================================================
"""
class Mode: # obs = {s:[o]}
    def __init__(self,idRequest,idMode,utility,obs,temporalUtility):
        if not config.glob.score_temporel:
            if not temporalUtility==0:
                print("requête",idRequest,"score temporel hors domaine.")
                print("mode :",idRequest,idMode,"contenu :",obs)
            assert(temporalUtility==0)
        self.idRequest = idRequest
        self.idMode = idMode
        self.utility = utility
        self.temporalUtility = temporalUtility
        assert(utility>0 or self.idRequest == -1)
        self.observations = obs
        self.pairs = []
        for s in self.observations:
            for o in self.observations[s]:
                self.pairs.append((s,o))
    
    def acceptable(self,constellation):
        return constellation.getRequest(self.idRequest).acceptable(self.pairs,constellation)
        
    def getTemporalUtility(self):
        return self.temporalUtility
    
    def getPairs(self):
        return self.pairs
    
    def getId(self):
        return self.idMode

    def getRequest(self):
        return self.idRequest

    def getUtility(self):
        return self.utility
        
    def getActivities(self):
        return self.observations
    
    def getActivitiesBySatellite(self,s):
        return self.observations.get(s,[])
    
    def __str__(self):
        return "Mode "+str((self.idRequest,self.idMode)) +" : "+str(self.observations) + " utility = " +str(self.utility)


class RequestOrdering:
    def __init__(self,requestId,modeId,noisedUtility):
        self.requestId = requestId
        self.modeId = modeId
        self.noisedUtility = noisedUtility

    def getId(self):
        return self.requestId
    
    def getIdMode(self):
        return self.modeId
    
    def getUtility(self):
        return self.noisedUtility
    
    def __eq__(self,other):
        if self.getId() ==-1:
            if other.getId() == -1:
                assert(self.getUtility()[0] == 0 and other.getUtility()[0] ==0)
                return True
            else:
                return False
        else:
            return self.getId() == other.getId() and self.getIdMode() == other.getIdMode()
    
    def __str__(self):
        return "Order("+str(self.requestId)+","+str(self.modeId)+":"+str(self.noisedUtility)+")"
    
    def __ne__(self,other):
        return not self.__eq__(other)
    
    def __lt__(self, other):
        if self.getId()==-1:
            return False
        else:
            return (self.getUtility() < other.getUtility()) or (self.getUtility() == other.getUtility() and self.getId() > other.getId())

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
                return (self.getUtility() > other.getUtility()) or (self.getUtility() == other.getUtility() and self.getId() < other.getId())

    def __ge__(self,other):
        return self.__gt__(other) or self.__eq__(other)
        
    def __le__(self,other):
        return self.__lt__(other) or self.__eq__(other)

class Request:
        def __init__(self,idRequest,priority,requestType,activities):
            self.actif = not config.getOptValue("dynamic")
            self.idRequest = idRequest
            self.priority = priority
            self.setType(requestType)
            configuredSeed = config.getOptValue("seed")
            self.observationRandomizer = rd.Random(configuredSeed)
            self.observationListRandomizer = rd.Random(configuredSeed)
            self.modesRandomizer = rd.Random(configuredSeed)
            self.init = False # modes déjà initialisés = False
            # time-slot -> liste de triplets (utility,satellite,indice activite)
            self.activitiesByTimeSlots = activities # {timeSlot:[(s,o,utility)]}
            # satellites -> liste des indices des activités
            self.activitiesListBySatellite= {}
            self.pairsSatellitesActivities = []
            self.mapRSO = {}
            self.modes = {}
            self.idMode = 0
            for timeSlot in self.activitiesByTimeSlots:
                for i,(r,s,o) in enumerate(self.activitiesByTimeSlots[timeSlot]):
                    self.mapRSO[o] = (timeSlot,i)
                    if(s,o) not in self.pairsSatellitesActivities:
                        self.pairsSatellitesActivities.append((s,o))
                    if s not in self.activitiesListBySatellite:
                        self.activitiesListBySatellite [s] = []
                    if o not in self.activitiesListBySatellite [s]:
                        self.activitiesListBySatellite [s].append(o)
            self.candidateMode = None
            self.activitiesHeld = {}
            self.generatedModesCounter = 0
            assert(len(self.pairsSatellitesActivities)==sum([len(self.activitiesListBySatellite [s]) for s in self.activitiesListBySatellite]))

        def setType(self,requestType):
            self.requestType = requestType
        
        def isActive(self):
            return self.actif
        
        def activate(self):
            self.actif = True
            
        def deactivate(self):
            self.actif = False
            
        def getCCAs(self,graphDep):
            res = []
            for (s,a) in self.pairsSatellitesActivities:
                cca = graphDep.getActivityCCA(a)
                if cca not in res:
                    res.append(cca)
            return res
               
        def resetCandidateModes(self):
            self.activitiesHeld = {}
            self.idMode = 0
            self.candidateMode = None
        
        def countGeneratedModes(self):
            return self.generatedModesCounter
        
        def findActivity(self,a):
            timeSlot,i = self.mapRSO[a]
            return self.activitiesByTimeSlots[timeSlot][i]
        
        def getModes(self):
            return self.modes.keys()
        
        def getMode(self,m):
            return self.modes[m]
            
        def getActivityUtility(self,constellation,a):
            return self.findActivity(a)
            
        def getTimeSlot(self,a):
            return self.mapRSO[a][0]
         
        def sortByTimeSlot(self,liste_act):
            tri = {}
            for a in liste_act:
                ts = self.getTimeSlot(a)
                if ts not in tri:
                    tri[ts] = []
                tri[ts].append(a)
            return tri
                
        def createRequestGraph(self):
            requestType = self.getType()
            r = self.idRequest
            if requestType == TYPE_MONO:
                graph = GraphMono(self)
            elif requestType == TYPE_STEREO:
                graph = GraphStereo(self)
            elif requestType == TYPE_PERIODIC:
                graph = GraphPeriodic(self)
            elif requestType == TYPE_SYSTEMATIC:
                graph = GraphSystematic(self)
            elif requestType == TYPE_DOWNLOAD:
                graph = GraphDownloads(self)
            else:
                raise ValueError("type inconnu",requestType)
            self.graph = graph
            self.dataToValidate = None
        
        def getGraph(self):
            return self.graph
        
        def getObservationsStructures(self,data):
            listObs = []
            obs = {}
            timeSlots = {}
            for timeSlot in self.activitiesByTimeSlots:
                for (r,s,a) in self.activitiesByTimeSlots[timeSlot]:
                    if a in data:
                        if s not in obs:
                            obs[s] = []
                        obs[s].append(a)
                        listObs.append((r,s,a))
                        timeSlots[(r,s,a)]=timeSlot
            return listObs,obs,timeSlots
        
        """
        # Utilisé dans le LNS : on ne sait pas encore si le mode sera effectivement candidat
        def genererMode(self,utilitys,inactives):
            cumulatedPathReward,contenu,verticesToValidate = self.graph.plusLongChemin(utilitys,inactives)
            self.dataToValidate = contenu.copy()
            return cumulatedPathReward,contenu,verticesToValidate
        """
        def validateCandidateMode(self,constellation):
            mode = self.modeToValidate
            self.modeToValidate = None
            self.modes[self.idMode] = mode
            self.idMode += 1
            self.generatedModesCounter += 1
            self.candidateMode = mode
            return mode
        
        def cancelMode(self):
            self.dataToValidate = None
            self.candidateMode = None
            self.modes[self.idMode] = None
        
        # Utilisé dans l'ancienne version de LNS
        def generateAndValidateMode(self,constellation,utilities,inactives):
            cumulatedPathReward,data,verticesToValidate = self.graph.plusLongChemin(utilities,inactives)
            if len(data)>0:
                listObs,obs,timeSlots = self.getObservationsStructures(data)
                recompense,temporalUtility = self.scoreObservationsList(listObs,constellation,timeSlots)
                mode = Mode(self.getId(),self.idMode,recompense,obs,temporalUtility)
                self.modes[self.idMode] = mode
                self.idMode += 1
                self.generatedModesCounter += 1
                self.candidateMode = mode
                return mode,verticesToValidate
            else:
                self.candidateMode = None
                return None,[]
            
        def holdObservation(self,explaination):
            self.activitiesHeld = {}
            for (s,o) in self.candidateMode.getPairs():
                if o not in explaination:
                    if s not in self.activitiesHeld:
                        self.activitiesHeld[s] = []
                    self.activitiesHeld[s].append(o)
        
        def stats(self,constellation):
            xyz = [0,0,0]
            minCoord = [np.Inf,np.Inf,np.Inf]
            maxCoord = [-np.Inf,-np.Inf,-np.Inf]
            nPoints = 0
            for timeSlot in self.activitiesByTimeSlots:
                for (r,s,o) in self.activitiesByTimeSlots[timeSlot]:
                    coords = constellation.getSatellite(s).getActivity(o).getCoordinates()
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
            for timeSlot in self.activitiesByTimeSlots:
                for (r,s,o) in self.activitiesByTimeSlots[timeSlot]:
                        coords = constellation.getSatellite(s).getActivity(o).getCoordinates()
                        xyz[0] += coords[0]
                        xyz[1] += coords[1]
                        xyz[2] += coords[2]
                        nPoints += 1
            xyz = tuple([elmt/nPoints for elmt in xyz])
            return xyz           
            
        def statsCoordinates(self,constellation):
            xyz = [0,0,0]
            minCoord = [np.Inf,np.Inf,np.Inf]
            maxCoord = [-np.Inf,-np.Inf,-np.Inf]
            nPoints = 0
            for timeSlot in self.activitiesByTimeSlots:
                for (r,s,o) in self.activitiesByTimeSlots[timeSlot]:
                        coords = constellation.getSatellite(s).getActivity(o).getCoordinates()
                        xyz[0] += coords[0]
                        xyz[1] += coords[1]
                        xyz[2] += coords[2]
                        minCoord = [min(minCoord[i],coords[i]) for i in range(3)]
                        maxCoord = [max(maxCoord[i],coords[i]) for i in range(3)]
                        nPoints += 1
            xyz = tuple([elmt/nPoints for elmt in xyz])
            return xyz,minCoord,maxCoord,nPoints
            
        def getPairs(self):
            return self.pairsSatellitesActivities
        
        def getPriority(self):
            return self.priority
        
        def getCurrentMode(self,constellation):
            if not self.init:
                self.resetModes(constellation)
            return self.candidateMode

        def getActivities(self):
            return self.activitiesListBySatellite
        
        def getActivitiesBySatellite(self,s):
            return self.activitiesListBySatellite.get(s,[])

        def getActivitiesByTimeSlots(self):
            return list(self.activitiesByTimeSlots.keys())
        
        def getObservationsByTimeSlot(self,timeSlot):
            return self.activitiesByTimeSlots[timeSlot]

        def getId(self):
            return self.idRequest

        def getUtility(self):
            return self.utility

        def getMode(self,m):
            return self.modes[m]

        def getModes(self):
            return self.modes

        def getType(self):
            return self.requestType
        
        def getKNextModes(self,constellation,k):
            liste = []
            for i in range(k):
                if i<len(self.modes):
                    liste.append(self.modes[i])
                else:
                    res = self.getNextModeWithoutExp(constellation)
                    if res is not None:
                        liste.append(res)
                    else:
                        break
            return liste
        
        def selectObservation(self,constellation,obsList,temporalKey,generator=None):
            assert(len(obsList)>0)
            # obsList : liste de (r,s,o) 
            criteria = lambda rso : self.scoreObs(rso,constellation,temporalKey)
            rMax = max([criteria(rso) for rso in obsList])
            obsMax = [rso for rso in obsList if criteria(rso)==rMax]
            if generator is None:
                gen = self.observationRandomizer
            else:
                gen = generator
            id_obs = gen.randint(0,len(obsMax)-1)
            return obsMax[id_obs]
        
        def selectObservationList(self,constellation,obsList,temporalKey,k,generator=None):
            assert(len(obsList)>0)         
            # obsList : liste de (r,s,o) 
            criteria = lambda rso : self.scoreObs(rso,constellation,temporalKey)
            if generator is None:
                gen = self.observationListRandomizer
            else:
                gen = generator
            # stoquer les valeurs
            values = {}
            for elmt in obsList:
                value = criteria(elmt)
                if value not in values:
                    values[value] = []
                values[value].append((value,elmt))
            # conserver les meilleurs groupes de valeurs  
            selectedElements = []
            nElmt = 0
            for value in sorted(values,reverse=True):
                if nElmt + len(values[value])<=k:
                    nElmt += len(values[value])
                    selectedElements += values[value]
                else:
                    gen.shuffle(values[value])
                    NelmtsDesires = k - nElmt
                    selectedElements += values[value][:NelmtsDesires]
                    break
            return selectedElements       
        
        def scoreObs(self,rso,constellation,temporalKey):
            # low : temps de reference
            t = (constellation.getSatellite(rso[1]).getActivity(rso[2]).getStartDate()+constellation.getSatellite(rso[1]).getActivity(rso[2]).getEndDate())/2
            if config.scoring.method=='relative':
                low = min([constellation.getSatellite(s).getActivity(a).getStartDate() for (s,a) in self.pairsSatellitesActivities])
                up = max([constellation.getSatellite(s).getActivity(a).getEndDate() for (s,a) in self.pairsSatellitesActivities])
            elif config.scoring.method=='global':
                low,up = config.donnees.tmin,config.donnees.tmax
            alpha = config.scoring.alpha_weather
            score = (rso[0],temporalKey(rso),0)
            assert(score>(0,0,0))
            return score
        
        def removeUnusedModes(self):
            for i in range(self.idMode):
                if i > 0 and i < self.idMode-1 and i in self.modes:
                    del self.modes[i]
        
        # L'interet est d'ajouter un mode en indice le plus haut
        # Il est la copie du mode d'indice m. Utile quand m est le meilleur mode
        # connu pour les algo qui utilise le dernier mode  comme point de départ  
        # dans la recherche locale.
        def shiftMode(self,m):
            if self.idRequest==-1:
                return self.idMode
            mode = deepcopy(self.modes[m])
            mode.idMode = self.idMode
            self.candidateMode = mode
            assert(self.idMode not in self.modes)
            self.modes[self.idMode] = mode
            self.idMode += 1
            self.generatedModesCounter += 1
            #self.removeUnusedModes()
            return mode.getId()
        
        # version générique : très couteuse
        def addMode(self,constellation,activitiesList):
            if config.getOptValue("verif"):
                assert(self.acceptable(activitiesList,constellation))
            pairs = len(activitiesList)==0 or type(activitiesList[0])==tuple
            if self.idRequest==-1:
                raise ValueError("can't create new mode for the download request")
            listRSO = []
            obs = {}
            TS = {}
            for timeSlot in self.activitiesByTimeSlots:
                for r,s,a in self.activitiesByTimeSlots[timeSlot]:
                    if (pairs and (s,a) in activitiesList) or (not pairs and a in activitiesList):
                        TS[a] = timeSlot
                        listRSO.append((r,s,a))
                        if s not in obs:
                            obs[s] = []
                        obs[s].append(a)
            assert(len(listRSO)>0)
            lexico = self.scoreObservationsList(listRSO, constellation, None)       
            mode = Mode(self.idRequest,self.idMode,lexico[0],obs,lexico[1])
            
            self.modes[self.idMode] = mode
            self.candidateMode = mode
            self.idMode += 1
            self.generatedModesCounter += 1
            return mode
 
        def estimateScore(self,constellation,activitiesList):
            listRSO = []
            obs = {}
            TS = {}
            for timeSlot in self.activitiesByTimeSlots:
                for r,s,a in self.activitiesByTimeSlots[timeSlot]:
                    if (s,a) in activitiesList:
                        TS[a] = timeSlot
                        listRSO.append((r,s,a))
                        if s not in obs:
                            obs[s] = []
                        obs[s].append(a)
            assert(len(listRSO)>0)
            return self.scoreObservationsList(listRSO, constellation, None)
    
        def getActivitiesInCCA(self,graphDep,liste_cca):
            actitivitiesCCA = []
            for timeSlot in self.activitiesByTimeSlots:
                for (r,s,a) in self.activitiesByTimeSlots[timeSlot]:
                    if not config.getOptValue("dynamic") or self.isActive():
                        cca = graphDep.getActivityCCA(a)
                        if cca in liste_cca:
                            actitivitiesCCA.append(a)
            return actitivitiesCCA
        
        # contenu actuel : activities du mode retenu privé de celles sur les cca a explorer
        def generateModesInCCA(self,graphDep,constellation,currentContent,listCCAToExplore,allowExternalModes=True):
            externalModes = False
            externalActivities = len(currentContent)>0
            newModes = []
            actitivitiesCCA = self.getActivitiesInCCA(graphDep, listCCAToExplore)
            maxAddSize = len(actitivitiesCCA)
            combination = []
            for size in range(maxAddSize+1):
                combination += list(map(list,itertools.combinations(actitivitiesCCA,size)))
            for combi in combination:
                newModes = currentContent + combi
                if len(combi)>0 or (allowExternalModes and self.acceptable(newModes,constellation)):
                    newModes.append(self.addMode(constellation, newModes))
                    if len(combi)==0:
                        externalModes = True
            return newModes,externalModes,externalActivities
            
        def filterPresentModes(self,activitiesList,constellation):
            self.idMode = 0
            if config.verifMode():
                assert(self.acceptable(activitiesList,constellation))
            return self.addMode(constellation,activitiesList)
    
class MonoRequest(Request):
    def __init__(self,idRequest,priority,activities):
        super().__init__(idRequest,priority,TYPE_MONO,activities)

    def generateModesInCCA(self,graphDep,constellation,currentContent,listCCAToExplore,allowExternalModes=True):
        externalActivities = len(currentContent)>0
        if len(currentContent)==1:
            if allowExternalModes:
                return [],False,externalActivities
            else:
                return []
        assert(len(currentContent)==0)
        newModes = []
        actitivitiesCCA = self.getActivitiesInCCA(graphDep, listCCAToExplore)
        for a in actitivitiesCCA:
            newModes = [a]
            newModes.append(self.addMode(constellation,newModes))
        return newModes,False,externalActivities
    
    def getBestModeWithoutInactives(self,constellation,inactives):
        if not config.getOptValue("dynamic") or self.isActive():
            activities = [rso for rso in self.candidateActivities if rso[2] not in inactives]
            if len(activities)>0:
                temporalKey = lambda rso : self.temporalUtility(constellation,rso)
                rso = self.selectObservation(constellation,activities,temporalKey)
                utility,temporalUtility = self.scoreObservationsList([rso],constellation,temporalKey)
                obs =  {rso[1]:[rso[2]]}
                self.modeToValidate = Mode(self.idRequest,self.idMode,utility,obs,temporalUtility)
                return self.modeToValidate
            else:
                self.modeToValidate = None
                return None
        else:
            return None
    
    def acceptable(self,activitiesList,constellation):
        pairs = len(activitiesList)==0 or type(activitiesList[0])==tuple
        if pairs:
            act = self.pairsSatellitesActivities
        else:
            act = [x[1] for x in self.pairsSatellitesActivities]
        if not sum([x in act for x in activitiesList])==len(activitiesList):
            return False
        return len(activitiesList)==1
    
    def getKBestModes(self,constellation,k):
        assert(self.modes=={})
        modesList = []
        if len(self.candidateActivities)==0:
            self.candidateMode = None
            return []
        else:
            start = time()
            temporalKey = lambda rso : self.temporalUtility(constellation,rso)
            observationsList = self.selectObservationList(constellation,self.candidateActivities,temporalKey,k,generator=None)
            end = time()

            for (score,rso) in observationsList:
                start = time()
                utility,temporalUtility,_ = score
                obs_mode =  {rso[1]:[rso[2]]}
                self.candidateMode = Mode(self.idRequest,self.idMode,utility,obs_mode,temporalUtility)
                self.modes[self.idMode] = self.candidateMode
                self.idMode += 1    
                self.generatedModesCounter += 1
                modesList.append(self.candidateMode)
                end = time()
        return modesList

    def resetModes(self,constellation,conserveOld=True,initFirstMode=True):
        self.init = True
        self.resetCandidateModes()
        self.candidateActivities = deepcopy(self.activitiesByTimeSlots[0])
        self.firstDate = min([constellation.getSatellite(rso_o[1]).getActivity(rso_o[2]).getStartDate() for rso_o in self.candidateActivities])
        self.lastDate = max([constellation.getSatellite(rso_o[1]).getActivity(rso_o[2]).getEndDate() for rso_o in self.candidateActivities])
        if not conserveOld:
            self.modes = {}
            self.idMode = 0        
        if initFirstMode:
            self.buildMode(constellation)
        
    def scoreObservationsList(self,obsList,constellation,timeSlots):
        assert(len(obsList)==1)
        temporalKey = lambda x : self.temporalUtility(constellation,x)
        lexico = [self.scoreObs(rso,constellation,temporalKey) for rso in obsList][0]
        return (lexico[0],lexico[1])
     
    def shiftActivities(self,explaination):
        assert(explaination!=[])
        lenBefore = len(self.candidateActivities)
        #copie = self.candidateActivities.copy()
        for i,x in enumerate(self.candidateActivities):
            if x[2] in explaination:
                self.candidateActivities.pop(i)
                break
        #self.candidateActivities = [x for x in self.candidateActivities if x[2] not in explaination]
        lenAfter = len(self.candidateActivities)
        assert(lenBefore == lenAfter + len(explaination))
        
    def buildMode(self,constellation):
        if len(self.candidateActivities)==0:
            self.candidateMode = None
        else:
            temporalKey = lambda rso : self.temporalUtility(constellation,rso)
            rso = self.selectObservation(constellation,self.candidateActivities,temporalKey)
            utility,temporalUtility = self.scoreObservationsList([rso],constellation,temporalKey)
            obs =  {rso[1]:[rso[2]]}
            self.candidateMode = Mode(self.idRequest,self.idMode,utility,obs,temporalUtility)
            self.modes[self.idMode] = self.candidateMode
            self.idMode += 1
            self.generatedModesCounter += 1
     
    def temporalUtility(self,constellation,rso):
        if config.glob.score_temporel:
            t = constellation.getSatellite(rso[1]).getObservation(rso[2]).getStartDate() - self.firstDate
            start = 3*60*60
            end = self.lastDate - self.firstDate
            if t<start:
                return 1
            else:
                return 1 - (t-start)/(end-start)
        else:
            return 0
     
    def getNextMode(self,explaination,constellation):
        if self.candidateMode is not None:
            self.shiftActivities(explaination)
            self.buildMode(constellation)
        return self.candidateMode
    
    def getNextModeWithoutExp(self,constellation):
        if self.candidateMode is None:
            return None
        exp = [self.candidateMode.getPairs()[0][1]] # 1er couple (s,o) : obs
        return self.getNextMode(exp,constellation)

class LongMonoRequest(MonoRequest):
    def __init__(self,idRequest,priority,activities):
        self.multiplicateur = 4
        for t in activities:
            for i in range(len(activities[t])):
                activities[t][i] = (activities[t][i][0]*self.multiplicateur,activities[t][i][1],activities[t][i][2])
        super().__init__(idRequest, priority, activities)
        self.setType(TYPE_LONG_MONO) # ecraser le type
    
class StereoRequest(Request):
    def __init__(self,idRequest,priority,activities):
        super().__init__(idRequest,priority,TYPE_STEREO,activities)
        assert(len(self.pairsSatellitesActivities)==sum([len(self.activitiesListBySatellite [s]) for s in self.activitiesListBySatellite]))

    def generateModesInCCA(self,graphDep,constellation,currentContent,listCCAToExplore,allowExternalModes=True):
        if not config.getOptValue("dynamic") or self.isActive():
            externalActivities = len(currentContent)>0
            if len(currentContent)==2:
                if allowExternalModes:
                    return [],False,externalActivities
                else:
                    return []
            assert(len(currentContent)==0)
            newModes = []
            actitivitiesCCA = self.getActivitiesInCCA(graphDep, listCCAToExplore)
            sortByPairs = self.sortByTimeSlot(actitivitiesCCA)
            for idPair in sortByPairs:
                assert(len(sortByPairs[idPair])==2)
                newModes = sortByPairs[idPair]
                newModes.append(self.addMode(constellation,newModes))
            return newModes,False,externalActivities
        else:
            return [],False,len(currentContent)>0
            
            
    def acceptable(self,activitiesList,constellation):
        if len(activitiesList)!=2:
            return False
        pairs = len(activitiesList)==0 or type(activitiesList[0])==tuple
        if pairs:        
            t1 = self.getTimeSlot(activitiesList[0][1])
            t2 = self.getTimeSlot(activitiesList[1][1])
        else:
            t1 = self.getTimeSlot(activitiesList[0])
            t2 = self.getTimeSlot(activitiesList[1])
        return t1 == t2 
    
    def getStereoPairs(self):
        pairs = []
        for ts in self.activitiesByTimeSlots:
            pairs.append([(x[1],x[2]) for x in self.activitiesByTimeSlots[ts].copy()])
        return pairs
    
    def selectActivityPairs(self,listepairs,constellation):
        rMax = max(self.scorePair(x,constellation) for x in listepairs)
        pairsMax = [x for x in listepairs if self.scorePair(x,constellation) == rMax]
        id_couple = self.observationRandomizer.randint(0,len(pairsMax)-1)
        return pairsMax[id_couple]
    
    def getBestModeWithoutInactives(self,constellation,inactives):
        if not config.getOptValue("dynamic") or self.isActive():
            activities = [rso_pairs for rso_pairs in self.pairs if rso_pairs[0][2] not in inactives and rso_pairs[1][2] not in inactives]
            if len(activities)==0:
                return None
            temporalKey = lambda rso : self.temporalUtility(constellation,rso)
            listRSO = self.selectActivityPairs(activities,constellation)
            utility,temporalUtility = self.scoreObservationsList(listRSO,constellation,temporalKey)
            assert(len(listRSO)==2) # 2 obs
            assert(listRSO[0][1] == listRSO[1][1]) # même sat
            obs =  {listRSO[0][1]:[listRSO[0][2],listRSO[1][2]]}
            mode =  Mode(self.idRequest,self.idMode,utility,obs,temporalUtility)
            if self.acceptable(mode.getPairs(),constellation):
                self.modeToValidate = mode
                return self.modeToValidate
            else:
                return None
        else:
            None
        
    def scoreObservationsList(self,couple,constellation,timeSlots):
        assert(len(couple)==2)
        cleTemporelObs1 = lambda rso : self.temporalUtilityCouple(constellation,(rso,(0,0,0)))
        score1 = self.scoreObs(couple[0],constellation,cleTemporelObs1)
        temporalKeyObs2 = lambda rso : 0
        score2 = self.scoreObs(couple[1],constellation,temporalKeyObs2)
        return (score1[0] + score2[0],score1[1])
            
    def getIdPairs(self):
        return list(self.activitiesByTimeSlots.keys())

    def getObservationPair(self,idPaire):
        return self.activitiesByTimeSlots[idPaire]
    
    def scorePair(self,couple,constellation):
        cleTemporelObs1 = lambda rso : self.temporalUtilityCouple(constellation,(rso,(0,0,0)))
        score1 = self.scoreObs(couple[0],constellation,cleTemporelObs1)
        temporalKeyObs2 = lambda rso : 0
        score2 = self.scoreObs(couple[1],constellation,temporalKeyObs2)
        return (score1[0] + score2[0],score1[1])
    
    def resetModes(self,constellation,conserveOld=True,initFirstMode=True):
        self.init = True
        self.resetCandidateModes()
        self.pairs = [self.activitiesByTimeSlots[timeSlot] for timeSlot in self.activitiesByTimeSlots] 
        self.firstDate = min([constellation.getSatellite(rso_o[1]).getActivity(rso_o[2]).getStartDate() for couple in self.pairs for rso_o in couple])
        self.lastDate = max([constellation.getSatellite(rso_o[1]).getActivity(rso_o[2]).getEndDate() for couple in self.pairs for rso_o in couple])
        # a = (r,s,o) : observations differents avec le meme sat
        if not conserveOld:
            self.modes = {}
            self.idMode = 0          
        if initFirstMode:
            self.buildMode(constellation)
    
    def selectListOfPairs(self,constellation,listeCouples,temporalKey,k,generator=None):
        assert(len(listeCouples)>0)
        # obsList : liste de (r,s,o) 
        criteria = lambda couple : self.scorePair(couple,constellation)
        if generator is None:
            gen = self.observationListRandomizer
        else:
            gen = generator
        # stoquer les valeurs
        values = {}
        for elmt in listeCouples:
            value = criteria(elmt)
            if value not in values:
                values[value] = []
            values[value].append((value,elmt))
        # conserver les meilleurs groupes de valeurs  
        selectedElements = []
        nElmt = 0
        for value in sorted(values,reverse=True):
            if nElmt + len(values[value])<=k:
                nElmt += len(values[value])
                selectedElements += values[value]
            else:
                gen.shuffle(values[value])
                NelmtsDesires = k - nElmt
                selectedElements += values[value][:NelmtsDesires]
                break
        return selectedElements  

    def getKBestModes(self,constellation,k):
        assert(self.modes=={})
        modesList = []
        if len(self.pairs)==0:
            self.candidateMode = None
            return []
        else:
            start = time()
            temporalKey = lambda rso : self.temporalUtility(constellation,rso)
            observationsList = self.selectListOfPairs(constellation,self.pairs,temporalKey,k,generator=None)
            end = time()
            for (score,couple) in observationsList:
                start = time()
                utility,temporalUtility = score
                rso1,rso2 = couple
                obs_mode =  {rso1[1]:[rso1[2],rso2[2]]}
                self.candidateMode = Mode(self.idRequest,self.idMode,utility,obs_mode,temporalUtility)
                self.modes[self.idMode] = self.candidateMode
                self.idMode += 1 
                self.generatedModesCounter += 1
                modesList.append(self.candidateMode)
                end = time()
        return modesList
    
    def testExplainationInPair(self,explaination,couple):
        for a in explaination:
            if a == couple[0][2] or a == couple[1][2]:
                return True
        return False
        
    def shiftActivities(self,explaination):
        for i in range(len(self.pairs)-1,-1,-1):
            couple = self.pairs[i]
            if self.testExplainationInPair(explaination,couple):
                self.pairs.pop(i)
                if config.verifMode():
                    for c in self.pairs:
                        if( self.testExplainationInPair(explaination,c)):
                            print(self.idRequest,explaination,c)
                            assert(not self.testExplainationInPair(explaination,c))
                break
        
    def buildMode(self,constellation):
        if len(self.pairs)==0:
            self.candidateMode = None
        else:
            couple = self.selectActivityPairs(self.pairs,constellation)
            utility,temporalUtility = self.scorePair(couple,constellation)
            obs = {couple[0][1]:[couple[0][2],couple[1][2]]}
            self.candidateMode = Mode(self.idRequest,self.idMode,utility,obs,temporalUtility)
            self.modes[self.candidateMode.getId()] = self.candidateMode
            self.idMode += 1
            self.generatedModesCounter += 1
    
    def temporalUtilityCouple(self,constellation,couple):
        if config.glob.score_temporel:
            rso = couple[0]
            t = constellation.getSatellite(rso[1]).getObservation(rso[2]).getStartDate() - self.firstDate
            start = 3*60*60
            end = self.lastDate - self.firstDate
            if t<start:
                return 1
            else:
                return 1 - (t-start)/(end-start)
            
        else:
            return 0
        
    def getNextMode(self,explaination,constellation):
        assert(explaination!=[])
        if self.candidateMode is not None:
            self.i(explaination)
            self.buildMode(constellation)
        return self.candidateMode

    def getNextModeWithoutExp(self,constellation):
        if self.candidateMode is None:
            return None
        exp = [self.candidateMode.getPairs()[0][1],self.candidateMode.getPairs()[1][1]] # 1er couple (s,o) : obs
        return self.getNextMode(exp,constellation)
    
class PeriodicRequest(Request):
    def __init__(self,idRequest,priority,activities):
        super().__init__(idRequest,priority,TYPE_PERIODIC,activities)
        assert(len(self.pairsSatellitesActivities)==sum([len(self.activitiesListBySatellite [s]) for s in self.activitiesListBySatellite]))
        configuredSeed = config.getOptValue("seed")
        self.generatorPlot = {plot : rd.Random(plot+configuredSeed) for plot in self.activitiesByTimeSlots}

    def generateModesInCCA(self,graphDep,constellation,currentContent,listCCAToExplore,allowExternalModes=True):
        if not config.getOptValue("dynamic") or self.isActive():
            externalModes = False
            externalActivities = len(currentContent)>0
            # I. créer les obs candidates à l'ajout
            time_slots_occupes = [self.getTimeSlot(a) for a in currentContent]
            newModes = []
            actitivitiesCCA = self.getActivitiesInCCA(graphDep, listCCAToExplore)
            sortBySlot = self.sortByTimeSlot(actitivitiesCCA)
            for ts in time_slots_occupes:
                if ts in sortBySlot:
                    del sortBySlot[ts]
            # II. générer les combination possibles d'ajout d'obs
            combination = [[]]
            for ts in sortBySlot:
                aux = deepcopy(combination)
                for a in sortBySlot[ts]:
                    aux_2 = deepcopy(aux)
                    for i in range(len(aux_2)):
                        aux_2[i].append(a)
                    combination += aux_2
            # III. combiner les combination au contenu actuel
            for combi in combination:
                newModes = combi+currentContent
                if len(combi)>0 or (allowExternalModes and self.acceptable(newModes,constellation)):
                    newModes.append(self.addMode(constellation,newModes))
                    if len(combi)==0:
                        externalModes = True
            return newModes,externalModes,externalActivities
        else:
            return [],False,len(currentContent)>0
    
    def acceptable(self,activitiesList,constellation):
        
        pairs = len(activitiesList)==0 or type(activitiesList[0])==tuple
        
        # vérifier qu'il y a du contenu
        if len(activitiesList)==0:
            return False
        # vérifier que les activités correspondent bien à la requête
        if pairs:
            act = self.pairsSatellitesActivities
        else:
            act = [x[1] for x in self.pairsSatellitesActivities]
        if not sum([x in act for x in activitiesList])==len(activitiesList):
            return False
        # vérifier l'unicité des time slots
        plots = {}
        for X in activitiesList:
            if pairs:
                (s,a) = X
            else:
                a = X
            ts = self.getTimeSlot(a)
            if ts in plots:
                return False # time slot doublement servi
            plots[ts] = [rso for rso in self.activitiesByTimeSlots[ts] if rso[2]==a][0]
        # vérifier la positivité du score    
        if self.scoreActivityList(plots,constellation)[0]<=0:
            return False

        return True
    
    def getTimeSlots(self):
        return {ts : [(x[1],x[2]) for x in self.activitiesByTimeSlots[ts]] for ts in self.activitiesByTimeSlots}
        
    def shiftActivities(self,explaination):
        assert(explaination!=[])
        for timeSlot in self.candidateActivities:
            self.candidateActivities[timeSlot] = [x for x in self.candidateActivities[timeSlot] if x[2] not in explaination]
        self.holdObservation(explaination)
        
    def scoreActivityList(self,plots,constellation):
        liste_score = [self.scoreObs(plots[timeSlot],constellation,lambda rso : self.scorePlot(constellation,rso,timeSlot)) for timeSlot in plots]
        #print(plots,[x[0] for x in liste_score ])
        Nplots = len(list(self.activitiesByTimeSlots.keys()))
        maxGap = 0
        gap = 0
        for ts in sorted(self.activitiesByTimeSlots):
            if len(plots.get(ts,[]))>0:
                if gap>maxGap:
                    maxGap = gap
                gap = 0
            else:
                gap += 1
        alpha = 0.5
        utility = sum([x[0] for x in liste_score ]) + alpha*(Nplots-maxGap)
        if len(liste_score)>0:
            if config.glob.score_temporel:
                return (utility,np.mean([x[1] for x in liste_score]))
            else:
                return (utility,0)
        else:
            return (0,0)
    
    def scoreObservationsList(self,obsList,constellation,timeSlots):
        plots = {}
        for rso in obsList:
            timeSlot = self.getTimeSlot(rso[2])
            plots[timeSlot] = rso
        return self.scoreActivityList(plots,constellation)
    
    def loadHeldActivities(self,obs):
        remainingTimeSlot = list(self.candidateActivities.keys())
        listeAct = []
        for s in self.activitiesHeld:
            obs[s] = []
            for o in self.activitiesHeld[s]:
                obs[s].append(o)
                for timeSlot in self.candidateActivities:
                    match_obs = [x for x in self.candidateActivities[timeSlot] if x[2]==o]
                    if len(match_obs)==1:
                        remainingTimeSlot.remove(timeSlot)
                        listeAct.append((timeSlot,match_obs[0]))
                        break
        return remainingTimeSlot,listeAct

    def getBestModeWithoutInactives(self,constellation,inactives):
        if not config.getOptValue("dynamic") or self.isActive():
            listRSO = []
            obs = {}
            observationList = []
            for ts in self.activitiesByTimeSlots:
                activities = [rso for rso in self.activitiesByTimeSlots[ts] if rso[2] not in inactives]
                temporalKey = lambda rso : self.scorePlot(constellation,rso,ts)
                if len(activities)>0:
                    rso = self.selectObservation(constellation,activities,temporalKey)
                    observationList.append(rso[2])
                    listRSO.append(rso)
                    if rso[1] not in obs:
                        obs[rso[1]] = []
                    obs[rso[1]].append(rso[2])
            utility,temporalUtility = self.scoreObservationsList(listRSO,constellation,temporalKey)
            if utility<=0:
                return None
            mode =  Mode(self.idRequest,self.idMode,utility,obs,temporalUtility)
            if self.acceptable(mode.getPairs(),constellation):
                self.modeToValidate = mode
                return self.modeToValidate
            else:
                return None
        else:
            return None
                
    def buildMode(self,constellation,autoriser_degradation=True):
        # si degradation non autorisé et un slot n'a plus d'obs : plus de mode
        if not autoriser_degradation:
            for timeSlot in self.candidateActivities:
                if len(self.candidateActivities[timeSlot])==0:
                    self.candidateMode = None
                    return
        obs = {}
        remainingTimeSlot,listeAct = self.loadHeldActivities(obs)
        plots = {x[0] : x[1] for x in listeAct}
        for timeSlot in remainingTimeSlot:
            if len(self.candidateActivities[timeSlot])>0:
                temporalKey = lambda rso : self.scorePlot(constellation,rso,timeSlot)
                rso = self.selectObservation(constellation,self.candidateActivities[timeSlot],temporalKey,generator=self.generatorPlot[timeSlot])
                listeAct.append(rso)
                plots[timeSlot] = rso
                s = rso[1]
                if s not in obs:
                    obs[s] = []
                obs[s].append(rso[2])
        # pas d'obs trouvées => mode nul        
        if len(list(obs.keys()))==0:
            self.candidateMode = None
            return
        
        utility,temporalUtility = self.scoreActivityList(plots,constellation)
        self.plots_candidat = plots
        self.candidateMode = Mode(self.idRequest,self.idMode,utility,obs,temporalUtility)
        self.modes[self.idMode] = self.candidateMode
        self.idMode += 1
        self.generatedModesCounter += 1
    
    def scorePlot(self,constellation,rso,timeSlot):
        t = constellation.getSatellite(rso[1]).getObservation(rso[2]).getStartDate()
        start = 30*60 + self.plots[timeSlot] # 30 minutes
        end = 90*60 + self.plots[timeSlot] # 1h30 : score nul
        if t<start:
            return 1
        else:
            return max(1 - (abs(t-start))/(end-start),0)
        
    def temporalUtility(self,constellation,plots):
        if config.glob.score_temporel:
            return np.mean([self.scorePlot(constellation,plots[timeSlot],timeSlot) for timeSlot in plots])
        else:
            return 0
        
    def resetModes(self,constellation,conserveOld=True,initFirstMode=True):
        self.init = True
        self.candidateActivities = deepcopy(self.activitiesByTimeSlots)
        self.resetCandidateModes()
        self.plots = {}
        startDate = lambda rso: constellation.getSatellite(rso[1]).getObservation(rso[2]).getStartDate()
        endDate = lambda rso : constellation.getSatellite(rso[1]).getObservation(rso[2]).getEndDate()
        milieu = lambda rso : (startDate(rso)+endDate(rso))/2
        for timeSlot in self.candidateActivities:
            self.plots[timeSlot] = np.mean([milieu(rso) for rso in self.candidateActivities[timeSlot]])
        if not conserveOld:
            self.modes = {}
            self.idMode = 0          
        if initFirstMode:
            self.buildMode(constellation)
       
    def getKBestModes(self,constellation,k):
        plots_reduits = {}
        for plot in self.plots:
            if len(self.activitiesByTimeSlots[plot])<=k:
                plots_reduits[plot] = self.activitiesByTimeSlots[plot]
            else: 
                plots_reduits[plot] = heapq.nlargest(k, self.activitiesByTimeSlots[plot], key=itemgetter(0))
        combination = [ {} ]
        modesList = []
        for plot in plots_reduits:
            combination = [dict(chain.from_iterable(d.items() for d in (combi, {plot:x}))) for combi in combination for x in plots_reduits[plot] ]
        if len(combination)<=k:
            best_obs = combination
        else:
            best_obs = heapq.nlargest(k,combination,key=lambda x : self.scoreActivityList(x,constellation))
        for comb in best_obs:
            obs = self.rsoToDict([comb[plot] for plot in comb])
            utility,temporalUtility = self.scoreActivityList(comb,constellation)
            m = Mode(self.idRequest,self.idMode,utility,obs,temporalUtility)
            modesList.append(m)
            self.modes[self.idMode] = m
            self.idMode += 1
            self.generatedModesCounter += 1
        return modesList
                                  
    def rsoToDict(self,rso_list):
        obs = {}
        for (r,s,o) in rso_list:
            if s not in obs:
                obs[s] = []
            obs[s].append(o)
        return obs
    
    def getNextMode(self,explaination,constellation):
        assert(explaination!=[])
        if self.candidateMode is not None:
            self.shiftActivities(explaination)
            self.buildMode(constellation)
        return self.candidateMode

    def getNextModeWithoutExp(self,constellation):
        if self.candidateMode is None:
            return None
        tmp = sorted([o for (s,o) in self.candidateMode.getPairs()])
        id_mode = self.modesRandomizer.randint(0,len(tmp)-1)
        exp = [tmp[id_mode]] # 1er couple (s,o) : obs
        return self.getNextMode(exp,constellation)
    
class SystematicRequest(Request):
    def __init__(self,idRequest,priority,activities):
        super().__init__(idRequest,priority,TYPE_SYSTEMATIC,activities)
        assert(len(self.pairsSatellitesActivities)==sum([len(self.activitiesListBySatellite [s]) for s in self.activitiesListBySatellite]))

    def generateModesInCCA(self,graphDep,constellation,currentContent,listCCAToExplore,allowExternalModes=True):
        if not config.getOptValue("dynamic") or self.isActive():
            externalModes = False
            externalActivities = len(currentContent)>0
            if len(currentContent)>0:
                sat = constellation.getSatelliteActivity(currentContent[0])
                if config.getOptValue("verif"):
                    constellation.getSatelliteActivity(currentContent[0]) == sat
                possible_sat = [sat]
            else:
                possible_sat = list(self.activitiesListBySatellite.keys())
            newModes = []
            actitivitiesCCA = self.getActivitiesInCCA(graphDep, listCCAToExplore)
            
            for s in possible_sat:
                candidateActivities = [x for x in self.activitiesListBySatellite [s] if x in actitivitiesCCA]
                combination = []
                for size in range(len(candidateActivities)+1):
                    combination += list(map(list,itertools.combinations(candidateActivities,taille)))
                for combi in combination:
                    newModes = combi+currentContent
                    if len(combi)>0 or (allowExternalModes and self.acceptable(newModes,constellation)):
                        newModes.append(self.addMode(constellation,newModes))
                    if len(combi)==0 and self.acceptable(newModes,constellation):
                        externalModes = True
            return newModes,externalModes,externalActivities
        else:
            return [],False,len(currentContent)>0
        
    def acceptable(self,activitiesList,constellation):
        if len(activitiesList)==0:
            return False
        pairs = len(activitiesList)==0 or type(activitiesList[0])==tuple
        if pairs:        
            act = self.pairsSatellitesActivities
        else:
            act = [x[1] for x in self.pairsSatellitesActivities]
                
        if not sum([x in act for x in activitiesList])==len(activitiesList):
            return False
        if pairs:
            sat = activitiesList[0][0]
        else:
            sat = constellation.getSatelliteActivity(activitiesList[0])
        if pairs:
            for s,a in activitiesList:
                if sat!=s:
                    return False
        else:
            for a in activitiesList:
                s = constellation.getSatelliteActivity(a)
                if sat!=s:
                    return False
        return True
            
    def resetModes(self,constellation,conserveOld=True,initFirstMode=True):
        self.init = True
        self.observations_candidates = {}
        for (r,s,o) in self.activitiesByTimeSlots[0]:
            if s not in self.observations_candidates:
                self.observations_candidates[s] = {}
            self.observations_candidates[s][o] = r
        self.resetCandidateModes()
        self.minDuration = min([constellation.getSatellite(rso[1]).getObservation(rso[2]).getDuration() for rso in self.activitiesByTimeSlots[0]]) 
        self.maxDuration = max([constellation.getSatellite(rso[1]).getObservation(rso[2]).getDuration() for rso in self.activitiesByTimeSlots[0]])
        if self.maxDuration==self.minDuration:
            self.delta = 1
        else:
            self.delta = self.maxDuration - self.minDuration
        if not conserveOld:
            self.modes = {}
            self.idMode = 0          
        if initFirstMode:
            self.buildMode(constellation)
    
    def getKBestModes(self,constellation,k):
        min_longueur = {s:max(1,len(self.activitiesListBySatellite [s])-k) for s in self.activitiesListBySatellite}
        comb = {s: [list(l) for n in range(min_longueur[s],len(self.activitiesListBySatellite [s])+1) for l in combinations([rso for rso in self.activitiesByTimeSlots[0] if rso[1]==s],n)] for s in self.activitiesListBySatellite}
        #print(comb)
        flatten = []
        for s in comb:
            for x in comb[s]:
                elmt = (x,self.scoreObservationsList(x,constellation))
                flatten.append(elmt)
        #print(flatten)
        observationList = [x[0] for x in heapq.nlargest(k, flatten, key=itemgetter(1))]
        modesList = []
        for l in observationList:
            obs = {}
            for rso in l:
                if rso[1] not in obs:
                    obs[rso[1]] = []
                obs[rso[1]].append(rso[2])
            utility,score_temporel = self.scoreObservationsList(l,constellation)
            m = Mode(self.idRequest,self.idMode,utility,obs,score_temporel)
            self.modes[self.idMode] = m
            modesList.append(m)
            self.idMode += 1
            self.generatedModesCounter += 1
        return modesList
                               
    def scoreObservationsList(self,obsList,constellation,timeSlots=None):
        temporalKey = lambda rso : (constellation.getSatellite(rso[1]).getObservation(rso[2]).getDuration() - self.minDuration)/self.delta
        liste_score = [self.scoreObs(rso,constellation,temporalKey) for rso in obsList]
        return (sum([x[0] for x in liste_score]),np.mean([x[1] for x in liste_score]))
    
    def buildMode(self,constellation):
        if self.observations_candidates=={}:
            sMax = None
        else:
            rMax = max([self.scoreObservationsList([(self.observations_candidates[s][o],s,o) for o in self.observations_candidates[s]],constellation) for s in self.observations_candidates])
            listSatMax = [s for s in self.observations_candidates if self.scoreObservationsList([(self.observations_candidates[s][o],s,o) for o in self.observations_candidates[s]],constellation) == rMax]
            sMax = min(listSatMax) # listSatMax[id_sat]

        if sMax is not None:
            utility,score_temporel = rMax
            obs = {sMax:list(self.observations_candidates[sMax].keys())}
            self.candidateMode = Mode(self.idRequest,self.idMode,utility,obs,score_temporel)
            self.modes[self.candidateMode.getId()] = self.candidateMode
            self.idMode += 1
            self.generatedModesCounter += 1
        else:
            self.candidateMode = None
            
    def getObsBySatellites(self,inactives):
        observationsBySat = {}
        for ts in self.activitiesByTimeSlots:
            for (r,s,a) in self.activitiesByTimeSlots[ts]:
                if a not in inactives:
                    if s not in observationsBySat:
                        observationsBySat[s] = []
                    observationsBySat[s].append((r,s,a))
        return observationsBySat
        
    def getBestModeWithoutInactives(self,constellation,inactives):
        if not config.getOptValue("dynamic") or self.isActive():
            temporalKey = lambda rso : self.temporalUtility(constellation,rso)
            listRSO = []
            obs = {}
            observationList = []
            rsoBySat = self.getObsBySatellites(inactives)
            if rsoBySat=={}:
                self.modeToValidate = None
                return None
            bestSat = None
            scoreSat = 0
            for sat in rsoBySat:
                utility,temporalUtility = self.scoreObservationsList(rsoBySat[sat],constellation,temporalKey)
                if utility>scoreSat:
                    scoreSat = utility
                    bestSat = sat
            
            obs = {bestSat:[rso[2] for rso in rsoBySat[bestSat]]}
            mode =  Mode(self.idRequest,self.idMode,utility,obs,temporalUtility)            
            if config.getOptValue("verif"):
                assert(self.acceptable(mode.getPairs(),constellation))
            self.modeToValidate = mode
            return self.modeToValidate
        else:
            return None
        
    def shiftActivities(self,explaination):
        for a in explaination:
            for s in self.observations_candidates:
                if a in self.observations_candidates[s]:
                    del self.observations_candidates[s][a]
        
        sat = list(self.observations_candidates.keys())
        for s in sat:
            if self.observations_candidates[s] == {}:
                del self.observations_candidates[s]

    def getNextModeWithoutExp(self,constellation):
        if self.candidateMode is None:
            return None
        tmp = sorted([o for (s,o) in self.candidateMode.getPairs()])
        id_mode = self.modesRandomizer.randint(0,len(tmp)-1)
        exp = [tmp[id_mode]] # 1er couple (s,o) : obs
        return self.getNextMode(exp,constellation)
    
    def testAllExplainationsAreLate(self,explaination):
        return sum([a not in [x[1] for x in self.candidateMode.getPairs()] for a in explaination]) == len(explaination)
        
    def getNextMode(self,explaination,constellation):
        assert(explaination!=[])
        if self.candidateMode is not None:
            isLate = self.testAllExplainationsAreLate(explaination)
            self.shiftActivities(explaination)
            if not isLate: # si l'explication provient d'un ancien mode inutile de changer encore de satellite
                self.buildMode(constellation)
        return self.candidateMode
    
class DownloadRequest(Request):
    def __init__(self,idRequest,priority,downloads):
        self.activitiesListBySatellite = downloads
        self.pairsSatellitesActivities = []
        rso = []
        for s in downloads:
            for v in downloads[s]:
                self.pairsSatellitesActivities.append((s,v))
                rso.append((0,s,v))
        assert(len(self.pairsSatellitesActivities)==sum([len(self.activitiesListBySatellite [s]) for s in self.activitiesListBySatellite]))
        super().__init__(idRequest,priority,TYPE_DOWNLOAD,{0:rso})
        self.requestType = TYPE_DOWNLOAD
        self.candidateMode = Mode(idRequest,0,0,self.activitiesListBySatellite,0)
        self.modes = {0:self.candidateMode}
        self.idMode = 0
    
    def acceptable(self,contenu,constellation):
        for (s,a) in self.pairsSatellitesActivities:
            if (s,a) not in contenu:
                return False
        return True
        
    def generateModesInCCA(self,graphDep,constellation,currentContent,listCCAToExplore,allowExternalModes=True):
        assert(not config.getOptValue("dynamic") or self.isActive())
        if config.getOptValue("verif"):
            actitivitiesCCA = self.getActivitiesInCCA(graphDep, listCCAToExplore)
            newModes = actitivitiesCCA+currentContent
            assert(sorted(newModes) == sorted([x[1] for x in self.modes[0].getPairs()]))
        return [self.modes[0]],False,len(currentContent)>0
        
    def scoreObservationsList(self,obsList,constellation,timeSlots):
        return (0,0)
    
    def getBestModeWithoutInactives(self,constellation,inactives):
        assert(self.isActive())
        if len(inactives)>0:
            print([constellation.getRequest(a) for a in inactives])
        assert(len(inactives)==0)
        self.modeToValidate = self.candidateMode
        return self.candidateMode

    def getNextModeWithoutExp(self,constellation):
        if self.candidateMode is None:
            return None
        exp = [] # 1er couple (s,o) : obs
        res =  self.getNextMode(exp,constellation)
        self.candidateMode = res
        return res
    
    def getKBestModes(self,constellation,k):
        exp = [] # 1er couple (s,o) : obs
        res =  self.getNextMode(exp,constellation)
        self.candidateMode = res
        return [res]
    
    def resetModes(self,constellation,conserveOld=True,initFirstMode=True):
        self.init = True
        pass
                  
    def getNextMode(self,explaination,constellation):
        mode = deepcopy(self.candidateMode)
        self.idMode += 1
        self.generatedModesCounterCounter += 1
        mode.idMode = self.idMode
        retrait = []
        for (s,o) in mode.pairs:
            if o in explaination:
                retrait.append((s,o))
        for (s,o) in retrait:
            mode.pairs.remove((s,o))
            del mode.observations[s][o]
        self.candidateMode = mode
        self.modes[1] = mode
        return mode
    
class Activity:
    def __init__(self,idActivity,satId,startDate,endDate,longitude,latitude,altitude):
        self.id = idActivity
        self.satId = satId
        self.startDate = startDate
        self.endDate = endDate
        self.longitude = longitude
        self.latitude = latitude
        self.altitude = altitude

    def getId(self):
        return self.id

    def getStartDate(self):
        return self.startDate

    def getEndDate(self):
        return self.endDate

    def getSat(self):
        return self.satId

    def getLongitude(self):
        return self.longitude

    def getLatitude(self):
        return self.latitude

    def getAltitude(self):
        return self.altitude

    def getCoordinates(self):
        return ( self.longitude,self.latitude,self.altitude)

    def isDependant(self,constellation,s,a):
        if s!=self.satId:
            return False
        other = constellation.satellites[self.satId].getActivity(a)
        left = other.getStartDate() <= self.getStartDate() and other.getEndDate() + constellation.satellites[s].getTransition(self.id,a) <= self.getStartDate()
        right = self.getStartDate() <= other.getStartDate() and self.getEndDate() + constellation.satellites[s].getTransition(a,self.id) <= other.getStartDate()
        return not(left or right)
    
    def distance(self,s,a2,constellation):
        coords2 = constellation.getSatellite(s).getActivity(a2).getCoordinates()
        return sqrt(sum([(self.getCoordonnes()[i]-coords2[i])**2 for i in range(3)]))
    
    def dependantTauMax(self,s,a,constellation,tau_max):
        other = constellation.satellites[s].getActivity(a)
        left = other.getStartDate() <= self.getStartDate() and other.getEndDate() + tau_max <= self.getStartDate()
        right = self.getStartDate() <= other.getStartDate() and self.getEndDate() + tau_max <= other.getStartDate()
        return not(left or right)
        
class Observation(Activity):
    def __init__(self,idObs,idSat,startDate,endDate,obsDuration,x,y,z,score,idRequest):
        super().__init__(idObs,idSat,startDate,endDate,x,y,z)
        self.duration = obsDuration
        self.scoringObs(score)
        self.idRequest = idRequest
    
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
        
    def getRequest(self):
        return self.idRequest
        
    def getScore(self):
        return self.score
    
    def getScorePlan(self,composantes,solCCAs):
        s = self.getSat()
        cca = composante.getActivityCCA(self.getIdentifiant())
        charge = solCCAs[s][cca].scoreCharge()
        return (self.score,-charge)
        
    def getDuration(self):
        return self.duration

    def getCoordinates(self):
        return (self.getLongitude(),self.getLatitude(),self.getAltitude())
    
    def __str__(self):
        return " " + Fore.BLUE + "[Observation, start : " + str(self.startDate) + ", id : " + str(self.id) + ", end : " + str(self.endDate) + ", duration : " + str(self.duration) + ", score : " +str(self.score) +"]" +Style.RESET_ALL

class Download(Activity):
    def __init__(self,data):
        super().__init__(int(data[0]),int(data[1]),float(data[2]),float(data[3]), float(data[4]),float(data[5]),float(data[6]))
        #super().__init__(int(data[0])+config.donnees.no,int(data[1]),float(data[2]),float(data[3]), float(data[4]),float(data[5]),float(data[6]))
        
    def __str__(self):
        return " " + Fore.GREEN + "[Download, start : "+str(self.startDate)+", id : "+str(self.id)+", end : "+str(self.endDate) + ", duration : " + str(self.getDuration()) +" ]" + Style.RESET_ALL
    
    def getDuration(self):
        fenetre = self.endDate - self.startDate
        return min(config.glob.allocation_vidage,0.95*fenetre)
    
    def getRequest(self):
        return -1
    
    def getScore(self):
        return 0
            
class Satellite:
    def __init__(self,identifiant,obs,vid):
        self.obs = obs
        self.id = identifiant
        self.downloads = vid

    def deleteActivity(self,a):
        if self.isDownload(a):
            del self.downloads[a]
        if self.isObservation(a):
            del self.obs[p]

    def isDownload(self,p):
        #return p>=config.donnees.no and p<config.donnees.no+config.donnees.nv
        return p in self.downloads    
    
    def isObservation(self,p):
        #return p>=0 and p<config.donnees.no
        return p in self.obs

    def getObservations(self):
        return self.obs

    def getObservation(self,p):
        return self.obs[p]

    def getId(self):
        return self.id

    def getDownload(self,p):
        return self.downloads[p]

    def getActivities(self):
        r = dict(self.obs)
        r.update(self.downloads)
        return r

    def getActivity(self,a):
        if self.isObservation(a):
            return self.obs[a]
        else:
            return self.downloads[a]

    def __str__(self):
        s = "IDENTIFIER : " + str(self.id) + '\n'
        s += "ACTIVITIES : \n"
        for o in self.obs :
            s += "\t" + str(self.obs[o]) + "\n"
        return s

    def getDownloads(self):
        return self.downloads
    
    def getTransition(self,a1,a2,modeleDeTransition):
        assert(not modeleDeTransition.estTimeDependent())
        res = modeleDeTransition.getTransition(self,a1,a2,config.glob.digits)
        return res
    
    def getTransitionTimeDependent(self,a1,a2,startManeuveringDate,modeleDeTransition):
        if modeleDeTransition.estTimeDependent():
            return modeleDeTransition.getTransition(a1,a2,startManeuveringDate)
        else:
            return modeleDeTransition.getTransition(self,a1,a2)

class Constellation:

    """
        =========================================================
                Initialisation de la constellation
        =========================================================
    """            
    def getRequestActivity(self,a):
        s = self.mappingSat[a]
        return self.dependencies[s][a]
    
    def initDependancesActivities(self):
        # mapping dep[s][a] => (requete,mode) qui contiennent a
        self.dependencies = {s : {} for s in self.satellites}
        for r in self.getAllRequests():
            for s in self.getRequest(r).getActivities():
                if s not in self.satellites:
                    die(s,self.satellites.keys(),self.getRequest(r).getActivities().keys())
                for a in  self.getRequest(r).getActivitiesBySatellite(s):
                    self.dependencies[s][a] = r
                
    def initTimeData(self):
        self.tmin = min([self.satellites[s].getActivity(p).getStartDate() for s in self.satellites for p in self.satellites[s].getActivities()])
        self.tmax = max([self.satellites[s].getActivity(p).getEndDate() for s in self.satellites for p in self.satellites[s].getActivities()])
        mins = {s : min([self.satellites[s].getObservation(p).getStartDate() for p in self.satellites[s].getObservations()]) for s in self.satellites}
        maxs = {s : max([self.satellites[s].getObservation(p).getEndDate() for p in self.satellites[s].getObservations()]) for s in self.satellites}
        self.horizonObs = { s : (mins[s],maxs[s]) for s in mins}
        mins = {s : min([self.satellites[s].getDownload(p).getStartDate() for p in self.satellites[s].getDownloads()]) for s in self.satellites}
        maxs = {s : max([self.satellites[s].getDownload(p).getEndDate() for p in self.satellites[s].getDownloads()]) for s in self.satellites}
        self.horizonVid = { s : (mins[s],maxs[s]) for s in mins} 
    
    def initGlobalData(self):
        config.donnees.ns = self.ns
        config.donnees.no = self.no
        config.donnees.nv = self.nv
        config.donnees.tmin = self.tmin
        config.donnees.tmax = self.tmax    
    
    def readObs(self,line,offset,idRequest):
        if not(len(line)==9+offset):
            print(line,len(line),offset)
            assert(False)
        idObs,idSat = int(line[0+offset]),int(line[1+offset])
        startDate,endDate = float(line[2+offset]),float(line[3+offset])
        obsDuration = float(line[4+offset])
        x,y,z = float(line[5+offset]),float(line[6+offset]),float(line[7+offset])
        score = float(line[8+offset])
        self.mappingSat[idObs] = idSat
        return Observation(idObs,idSat,startDate,endDate,obsDuration,x,y,z,score,idRequest)
    
    def deleteSatelliteWithoutObservations(self,vid,obs):
        del_vid = [s for s in vid if s not in obs]
        for s in del_vid:
            del vid[s]
    
    def statsCoordinates(self):
        xyz = [0,0,0]
        minCoord = [np.Inf,np.Inf,np.Inf]
        maxCoord = [-np.Inf,-np.Inf,-np.Inf]
        nPoints = 0
        for r in self.getAllRequests():
            xyz_r,minCoord_r,maxCoord_r,n_r = self.getRequest(r).statsCoordinates(self)
            xyz = [xyz[i]+xyz_r[i] for i in range(3)]
            minCoord = [min(minCoord[i],minCoord_r[i]) for i in range(3)]
            maxCoord = [max(maxCoord[i],maxCoord_r[i]) for i in range(3)]
            nPoints += n_r
        xyz = tuple([xyz[i]/nPoints for i in range(3)])
        return xyz,minCoord,maxCoord,nPoints
    
    def readDownloads(self,fichier,display):
        self.nv = int(fichier.readline())
        vid = {}
        for i in range(self.nv):
            line = fichier.readline().split(",")
            v = Download(line)
            s = v.getSat()
            if s not in vid:
                vid[s] = {}
            vid[s][v.getId()] = v
            self.mappingSat[v.getId()] = s
        return vid
    
    def readSatellites(self,obs,vid,display):
        satellites = {}
        for s in vid:
            if s in obs:
                satellites[s] = Satellite(s,obs[s],vid[s])
        self.ns = len(list(satellites.keys()))
        self.satellites = satellites        
                    
    def readRequests(self,fichier):
        line = fichier.readline().split('\n')[0].split(",")
        nRequests = int(line[0])
        satellitesObservations = {}
        no = 0
        self.requests = {}
        err = 0
        for r in range(nRequests):
            line = fichier.readline().split('\n')[0].split(",")
            idRequest,nObs,priority,requestType = int(line[0]),int(line[1]),0,line[2]
            no += nObs
            observations = {}
            for i in range(nObs):
                line = fichier.readline().split('\n')[0].split(",")
                if requestType == TYPE_PERIODIC or requestType == TYPE_STEREO:
                    timeSlot = int(line[0])
                    offset = 1
                else:
                    timeSlot = 0
                    offset = 0
                obs = self.readObs(line,offset,idRequest)
                if timeSlot not in observations:
                    observations[timeSlot] = []
                if (obs.getId() not in [x[2] for x in observations[timeSlot]]):
                    observations[timeSlot].append((obs.getScore(),obs.getSat(),obs.getId()))
                    s = obs.getSat()
                    if s not in satellitesObservations:
                        satellitesObservations[s] = {}
                    satellitesObservations[s][obs.getId()] = obs
                else:
                    err += 1
            if r<config.getOptValue("test_req") and sum([len(observations[s])!=0 for s in observations]):
                if requestType==TYPE_MONO:
                    req = MonoRequest(idRequest,priority,observations)
                elif requestType==TYPE_LONG_MONO:
                    req = LongMonoRequest(idRequest,priority,observations)
                elif requestType==TYPE_STEREO:
                    req = StereoRequest(idRequest,priority,observations)
                elif requestType==TYPE_PERIODIC:
                    req = PeriodicRequest(idRequest,priority,observations)
                elif requestType==TYPE_SYSTEMATIC:
                    if config.getOptValue("include_systematic"):
                        req = SystematicRequest(idRequest,priority,observations)
                    else:
                        req = None
                else:
                    print(line)
                    raise ValueError("type de requête inconnue : "+requestType)
                if req is not None:
                    self.requests[idRequest] = req
        self.no = no
        config.donnees.no = no
        return satellitesObservations
    
    def filterPresentModes(self,solution,graphDep,modeleDeTransition):
        contenu_modes = {}
        modes = []
        
        for s in solution.getSolCCAs():
            for cca in solution.getSolCCAs()[s]:
                for a in solution.getSolCCAs()[s][cca].getSequence():
                    r = self.getRequestActivity(a)
                    if r not in contenu_modes:
                        contenu_modes[r] = []
                    contenu_modes[r].append(a)
        for r in self.getRequests():
            if r!= -1: # Download
                self.getRequest(r).resetModes(self,conserveOld=False,initFirstMode=False)
                assert(len(self.getRequest(r).getModes())==0)
                if r in contenu_modes:
                    if self.getRequest(r).acceptable(contenu_modes[r],self):
                        mode = self.getRequest(r).filterPresentModes(contenu_modes[r],self)
                        assert(mode.getId()==0)
                        modes.append((r,mode.getId()))
                    else:        
                        for a in contenu_modes[r]:
                            (s,cca) = graphDep.getActivityCCA(a)
                            solution.getSolCCA(s,cca).removeActivity(self,a,modeleDeTransition)
            else:
                modes.append((r,0))
        solution.setModesRetenus(modes,self)
        
    def __init__(self,fileName,startDate,display=True):
        config.glob.filename = fileName.split("/")[-1].split(".")[0]
        comm = MPI.COMM_WORLD
        size = comm.Get_size()
        config.glob.size = size
        
        self.startDate = startDate
        self.fichier = fileName
        fichier = open(fileName,"r")
        self.mappingSat = {}
        obs = self.readRequests(fichier)
        vid = self.readDownloads(fichier,display)
        self.deleteSatelliteWithoutObservations(vid,obs)
        self.requests[-1] = DownloadRequest(-1,1,vid)
        self.readSatellites(obs,vid,display)
        self.initTimeData()
        self.initDependancesActivities()
        self.initGlobalData()
        
        
        if config.getOptValue("verif"):
            for r in self.getRequests():
                activities = self.extraireActivitiesRequetes()
                for s in activities:
                    for a in activities[s]:
                        assert(s==self.getSatelliteActivity(a))
        
        if display:
            printColor("Reading OK.\n",c='g')
            xyz,minCoord,maxCoord,nPoints = self.statsCoordinates()
            shift = getDisplayDepth()-1
            shiftLeftDisplay(shift)
            printColor("====================== [ System summary ] ======================",c='y')
            printColor("| Requests:",str(len(list(self.requests.keys()))),c='y')
            printColor("| Observations:",str(nPoints),c='y')
            printColor("| Latitude range: [",minCoord[0],",",maxCoord[0],"]",c='y')
            printColor("| Longitude range: [",minCoord[1],",",maxCoord[1],"]",c='y')
            printColor("| Mean target:",(xyz[0],xyz[1]),c='y')
            printColor("================================================================",c='y')
            shiftRightDisplay(shift)
        if config.getOptValue("dynamic"):
            minuteLastRequests=config.getOptValue("derniere_requete")
            startingProportion = config.getOptValue("proportion_depart")
            dt_minutes = config.getOptValue("dt_requests")
            self.requestDealer = RequestDealer(self,startDate,startingProportion,dt_minutes*60,minuteLastRequests*60)
            
    
    """
        =========================================================
                        GETTER
        =========================================================    
    """
    def getStartingRequests(self):
        return self.requestDealer.getStartingRequests()
    
    def getFilename(self):
        return config.glob.filename
    
    def getSatellite(self,s):
        return self.satellites[s]
    
    def getActivity(self,a,s=None):
        if s is None:
            s = self.getSatelliteActivity(a)
        return self.getSatellite(s).getActivity(a)
    
    def countActivities(self):
        return self.no+self.nv

    def getIdActivitiesMax(self):
        return max(list(self.mappingSat.keys()))
    
    def getSatellites(self):
        return self.satellites
    
    def getAllRequests(self):
        return list(self.requests.keys())
    
    def releaseNewRequests(self,graphDependances):
        return self.requestDealer.scruterArriveeNouvellesRequetes(self,graphDependances)
    
    def getRequests(self,allow_new=False,graphDependances=None):  
        if config.getOptValue("dynamic"):
            if allow_new:
                self.requestDealer.scruterArriveeNouvellesRequetes(self,graphDependances)
            else:
                assert(graphDependances is None)
            return [r for  r in self.requests if self.requests[r].isActive()]
        else:
            assert(not allow_new)
            return list(self.requests.keys())
    
    def getRequest(self,r):
        return self.requests[r]
    
    def isDownload(self,p):
        return p>=self.no
    
    def isObservation(self,p):
        return p<self.no
    
    def getSatelliteActivity(self,a):
        return self.mappingSat[a]
    
    def shuffleBatch(self,liste,quantiles):
        longueur = [int(q*len(liste)) for q in quantiles]
        for i,q in enumerate(longueur):
            if(i<len(longueur)-1):
                q2 = longueur[i+1]
                liste[q:q2] = rd.sample(liste[q:q2], len(liste[q:q2]))

    def meanObservationScore(self,candidateModes):
        score = []
        for (r,m) in candidateModes:
            for  (s,o) in self.getRequest(r).getPairs():
                if self.getSatellite(s).isObservation(o):
                    score.append(self.getSatellite(s).getObservation(o).getScore())
        if len(score)==0:
            return 0
        else:
            return np.mean(score)
    
    def printModes(self):
        for r in self.getRequests():
            for m in self.getRequest(r).getModes():
                print(self.getRequest(r).getMode(m))
    
    def extractActivitiesFromMode(self,constellation,r,m):
        act = []
        for (s,o,d) in constellation.getRequest(r).getMode(m).getPairs():
            if o not in act:
                act.append(o)
            if d not in act:
                act.append(d)
        return act  
    
    def extractActivitiesFromRequests(self):
        act = {}
        for r in self.getAllRequests():
            act_r = self.getRequest(r).getActivities()
            for s in act_r:
                if s not in act:
                    act[s] = []
                for o in act_r[s]:
                    act[s].append(o)
        return act 

    def extractActivitiesActiveRequests(self,graphDependances):
        act = {}
        for r in self.getRequests(graphDependances):
            act_r = self.getRequest(r).getActivities()
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
    def verifyModes(self):
        n_modes = sum([len(list(self.getRequest(r).getModes().keys())) for r in self.getAllRequests()])
        modes_faux = 0
        for r in self.getAllRequests():
            for m in self.getRequest(r).getModes():
                for s,o,d in self.getRequest(r).getMode(m).getPairs():
                    if(self.satellites[s].getObservation(o).getDuration()>self.satellites[s].getObservation(o).getEndDate()-self.satellites[s].getObservation(o).getStartDate()):
                        modes_faux += 1
                        break
        return modes_faux/n_modes
    
    
    """
        =========================================================  
                        GRAPHIQUES
        =========================================================  
    """
    def plotModesLosses(self):
        return self.solution.plotModesLosses()
    
    def plotNumberOfModes(self,display='boxplot'):
        if(display=='hist'):
            f,ax = plt.subplots()
            valeurs = np.unique([self.getRequest(r).getType() for r in self.getAllRequests()])
            for val in valeurs:
                ax.hist([len(self.getRequest(r).getModes()) for r in self.getAllRequests() if self.getRequest(r).getType()==val])
            plt.legend(valeurs)
            plt.title('count modes')
            #return f
        elif display=='boxplot':       
            df = pd.DataFrame([],columns=['count','type'])
            for r in self.getAllRequests():
                df = df.append({'count':len(self.getRequest(r).getModes()),'type':self.getRequest(r).getType()},ignore_index=True)
            sns.boxplot(y='count',x='type',data=df, palette="colorblind").set_title('count modes')
        else:
            raise NameError()
            
    def plotModesSizes(self,display='boxplot'):
        if(display=='hist'):
            f,ax = plt.subplots()
            valeurs = np.unique([self.getRequest(r).getType() for r in self.getAllRequests()])
            for val in valeurs:
                ax.hist([len(self.getRequest(r).getMode(m).getPairs()) for r in self.getAllRequests() for m in self.getRequest(r).getModes() if self.getRequest(r).getType()==val])
            plt.legend(valeurs)
            plt.title("mode size (number of observations)")
            #return f
        elif display=='boxplot':
            df = pd.DataFrame([],columns=['count','type'])
            for r in self.getAllRequests():
                for m in self.getRequest(r).getModes():
                    df = df.append({'count':len(self.getRequest(r).getMode(m).getPairs()),'type':self.getRequest(r).getType()},ignore_index=True)
            sns.boxplot(y='count',x='type',data=df, palette="colorblind").set_title('modes size')
        else:
            raise NameError()  
            
    def plotConstellationLoad(self):
        f,ax = plt.subplots(len(self.satellites.keys()),1)
        obs = {s : self.satellites[s].getObservations() for s in self.satellites}
        for s in obs:
            X = []
            Y = []
            temps = []
            for o in obs[s]:
                if obs[s][o].getStartDate() not in temps:
                    temps.append(obs[s][o].getStartDate())
                if obs[s][o].getEndDate() not in temps:
                    temps.append(obs[s][o].getEndDate())
            for t in temps:
                count = len([o for o in obs[s] if t>= obs[s][o].getStartDate() and t<= obs[s][o].getEndDate()])
                X.append(t)
                Y.append(count)
            XY = [(X[i],Y[i]) for i in range(len(X))]
            XY = sorted(XY,key = itemgetter(0))
            ax[s].plot([x[0] for x in XY],[x[1] for x in XY],color='black')
        f.suptitle('charge de la constellation')
        return f
    
    def plotSatelliteLoad(self,s):
        f,ax = plt.subplots(1,1)
        obs = {s : self.satellites[s].getObservations()}
        X = []
        Y = []
        temps = []
        for o in obs[s]:
            if obs[s][o].getStartDate() not in temps:
                temps.append(obs[s][o].getStartDate())
            if obs[s][o].getEndDate() not in temps:
                temps.append(obs[s][o].getEndDate())
        for t in temps:
            count = len([o for o in obs[s] if t>= obs[s][o].getStartDate() and t<= obs[s][o].getEndDate()])
            X.append(t)
            Y.append(count)
        XY = [(X[i],Y[i]) for i in range(len(X))]
        XY = sorted(XY,key = itemgetter(0))
        ax.plot([x[0] for x in XY],[x[1] for x in XY],color='black')
        f.suptitle('charge du satellite {}'.format(s))
        return f
        
    def plotMode(self,f,r,m):
        #f,ax = plt.subplots()
        for s,o,d in self.requests[r].getMode(m).getPairs():
            startDate,endDate = self.satellites[s].getObservations()[o].getStartDate(),self.satellites[s].getObservations()[o].getEndDate()
            plt.plot([startDate,endDate],[s,s],'-b',alpha=0.1)
        return f
    
    def plotModes(self,f):
        for r in self.requests:
            for m in self.requests[r].getModes():
                f = self.plotMode(f,r,m)
        return f
    
    def plotUtilitiesSamples(self):
        rr = rd.choices(list(self.getAllRequests().keys()),k=4)
        for r in rr:
            n,bins,patch = plt.hist([constellation.getRequest(r).getMode(m).getUtility()  for m in constellation.getRequest(r).getModes()])
            print(histToLatex(n,bins,"utility (request "+str(r)+")","mode count"))
            
    def plotLoad(self,annoter=False):
        f = self.solution.tracerActivite(annoter)
        return f
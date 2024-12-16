#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 13 10:52:39 2022

@author: ssquilla
"""

import random as rd
from time import time
import numpy as np
from math import floor,ceil
from copy import deepcopy

from ..Utils.config import *
global config
config = Config()
from ..Utils.Utils import printMaster,printColor,printOpen,printClose,printNoSpace


class RequestDealer:
    class BatchRequest:
        def __init__(self,arrivingDate,requestsList):
            self.arrivingDate = arrivingDate
            self.requestsList = requestsList
            
        def __str__(self):
            return "Batch("+str(self.requestsList)+", arriving : "+ str(self.arrivingDate) + ")"
        
        def getArrivingDate(self):
            return self.arrivingDate
        
        def getRequests(self):
            return self.requestsList
    # startingProportion : proportion de requêtes considérée à l'instant 0 du run
    # discretization : temps entre deux arrivée de batch de requêtes
    # lastRequestDate : durée écoulée avant l'arrivée de la dernière requête
    def __init__(self,constellation,startDate,startingProportion,discretization,lastRequestDate):
        self.shuffledRequests = rd.Random(config.getOptValue("seed"))
        self.startingProportion = startingProportion
        self.discretization = discretization
        self.lastRequestDate = lastRequestDate
        self.batchs = []
        self.startDate = startDate
        
        if not self.lastRequestDate%self.discretization==0:
            raise ValueError("Arriving date incompatible with discretization. Indicate another discretization step.")
        
        requestPool = constellation.getAllRequests()
        nInitialRequests = len(requestPool)
        self.shuffledRequests.shuffle(requestPool)
        initialNumber = int(self.startingProportion*len(requestPool))
        requestPool.remove(-1)
        requestPool.append(-1) # mise a la fin du vidage => fera parti du batch 0
        # creation du batch initial
        currentDate = 0
        batch_initial = self.createNewRequestBatch(currentDate,requestPool,initialNumber)
        self.batchs.append(batch_initial)
        # information sur les batchs
        remainingRequests = len(requestPool)
        n_batchs = int(self.lastRequestDate/self.discretization)
        batchSize = floor((nInitialRequests-initialNumber)/n_batchs)
        # creation des batchs
        batchContents = [[] for i in range(n_batchs)]
        j = 0
        while len(requestPool)>0:
            batchContents[j%n_batchs].append(requestPool.pop())
            j += 1
        for i in range(n_batchs):
            currentDate = (i+1)*self.discretization
            currentBatch = self.BatchRequest(currentDate,batchContents[i])
            self.batchs.append(currentBatch)
            
        self.remainingBatchs = deepcopy(self.batchs)
        change,requestsList = self.watchRequestArrival(constellation)
        self.firstRequests = requestsList
        
    def getInitialRequests(self):
        return self.firstRequests
    
    def createNewRequestBatch(self,date,requestPool,batchSize):
        liste_batch = []
        assert(batchSize<=len(requestPool))
        for i in range(batchSize):
            liste_batch.append(requestPool.pop())
        return self.BatchRequest(date,liste_batch)
    
    def watchRequestArrival(self,constellation,grapheDependances=None):
        timeElapsed = time()-self.startDate
        change = False
        requestsList = []
        while len(self.remainingBatchs)>0 and timeElapsed>=self.remainingBatchs[0].getArrivingDate():
            change = True
            batch = self.remainingBatchs.pop(0)
            printColor("Activating requests ",batch.getRequests(),depth=2,c='y')
            for r in batch.getRequests():
                constellation.getRequete(r).activer()
                if grapheDependances is not None:
                    grapheDependances.ajouterRequete(constellation,r)
            requestsList += batch.getRequests()
        return change,requestsList

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
        def __init__(self,date_arrivee,liste_requetes):
            self.date_arrivee = date_arrivee
            self.liste_requetes = liste_requetes
            
        def __str__(self):
            return "Batch("+str(self.liste_requetes)+", arrivée : "+ str(self.date_arrivee) + ")"
        
        def getDateArrivee(self):
            return self.date_arrivee
        
        def getRequetes(self):
            return self.liste_requetes
    # proportion_depart : proportion de requêtes considérée à l'instant 0 du run
    # discretisation : temps entre deux arrivée de batch de requêtes
    # date_derniere_requete : durée écoulée avant l'arrivée de la dernière requête
    def __init__(self,constellation,start_date,proportion_depart,discretisation,date_derniere_requete):
        self.melange_ordre_requetes = rd.Random(config.getOptValue("seed"))
        self.proportion_depart = proportion_depart
        self.discretisation = discretisation
        self.date_derniere_requete = date_derniere_requete
        self.batchs = []
        self.start_date = start_date
        
        if not self.date_derniere_requete%self.discretisation==0:
            raise ValueError("Date d'arrivée de la derniere requête incompatible avec la discrétisation. Indiquez un multiple du pas de temps.")
        
        
        pool_requetes = constellation.getToutesRequetes()
        N_requetes_initial = len(pool_requetes)
        self.melange_ordre_requetes.shuffle(pool_requetes)
        nombre_depart = int(self.proportion_depart*len(pool_requetes))
        pool_requetes.remove(-1)
        pool_requetes.append(-1) # mise a la fin du vidage => fera parti du batch 0
        # creation du batch initial
        t_courant = 0
        batch_initial = self.creerBatchNouvellesRequetes(t_courant,pool_requetes,nombre_depart)
        self.batchs.append(batch_initial)
        # information sur les batchs
        requetes_restantes = len(pool_requetes)
        n_batchs = int(self.date_derniere_requete/self.discretisation)
        taille_batch = floor((N_requetes_initial-nombre_depart)/n_batchs)
        # creation des batchs
        contenu_batchs = [[] for i in range(n_batchs)]
        j = 0
        while len(pool_requetes)>0:
            contenu_batchs[j%n_batchs].append(pool_requetes.pop())
            j += 1
        for i in range(n_batchs):
            t_courant = (i+1)*self.discretisation
            batch_courant = self.BatchRequest(t_courant,contenu_batchs[i])
            self.batchs.append(batch_courant)
            
        self.batchs_restants = deepcopy(self.batchs)
        change,liste_requetes = self.scruterArriveeNouvellesRequetes(constellation)
        self.premiere_requetes = liste_requetes
        
    def getRequetesDepart(self):
        return self.premiere_requetes
    
    def creerBatchNouvellesRequetes(self,date,pool_requetes,taille_batch):
        liste_batch = []
        assert(taille_batch<=len(pool_requetes))
        for i in range(taille_batch):
            liste_batch.append(pool_requetes.pop())
        return self.BatchRequest(date,liste_batch)
    
    def scruterArriveeNouvellesRequetes(self,constellation,grapheDependances=None):
        temps_ecoule = time()-self.start_date
        change = False
        liste_requetes = []
        while len(self.batchs_restants)>0 and temps_ecoule>=self.batchs_restants[0].getDateArrivee():
            change = True
            batch = self.batchs_restants.pop(0)
            printColor("Activation des requêtes ",batch.getRequetes(),depth=2,c='y')
            for r in batch.getRequetes():
                constellation.getRequete(r).activer()
                if grapheDependances is not None:
                    grapheDependances.ajouterRequete(constellation,r)
            liste_requetes += batch.getRequetes()
        return change,liste_requetes

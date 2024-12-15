#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 28 14:58:04 2022

@author: ssquilla
"""

"""
class Sequence:
    def __init__(self,liste_sequence):
        sequence_triee = sorted(liste_sequence)
        self.tete_sequence = NoeudSequence(sequence_triee)
        
class NoeudSequence:
    def __init__(self,list_sequence):
        assert(len(list_sequence)>0)
        self.contenu = list_sequence[0]
        if len(list_sequence)>1:
            self.suivant = NoeudSequence(list_sequence[1:])
        else:
            self.suivant = None
"""

class Oracle:
    def __init__(self,cca):
        self.identifiantCCA = cca
        self.sequence_interdites = []
        self.predictions_certaines = 0
        self.signalements = 0
    
    def getCCA(self):
        return self.identifiantCCA
    
    def AinB(self,a,b):
        for x in a:
            if x not in b:
                return False
        return True
        
    def signalerEchec(self,sequence):
        seq = sorted(sequence)
        for ss in self.sequence_interdites:
            if self.AinB(ss,seq):
                return
        self.signalements += 1
        self.sequence_interdites.append(seq)
    
    def predire(self,sequence):
        for seq in self.sequence_interdites:
            if self.AinB(seq,sequence):
                self.predictions_certaines += 1
                return False
        return True # optimiste : si on ne sait pas on suppose que ca passe
        
    def getSignalements(self):
        return self.signalements
     
    def getPredictions(self):
        return self.predictions_certaines
    
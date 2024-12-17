#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 11 10:29:31 2023

@author: ssquilla
"""
from mpi4py import MPI

from EOSSP.Utils.config import *
global config
config = Config()
instance = config.instance


class Communication:
    def __init__(self):
        # eviter d'envoyer les informations qu'on va recuperer à l'identique ensuite
        self.localData = {}
    
    def sendMessages(self,data,cpu):
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        assert(rank==0)
        data['time']=time()
        data['cpu'] = rank
        comm.send(data, dest=cpu, tag=rank)
    
    def extractDuplicatedData(self,data,cpu):
        # {'mode':mode,'cca':cca,'source':rank,'activites':activites,'time':time()}
        idCCA = data['cca'].getIdentifiant()
        self.localData[cpu] = {}
        self.localData[cpu]['sequence'] = data['cca'].getSequence().copy()
        self.localData[cpu]["mode"] = data['mode']
        self.localData[cpu]['activites'] = data['activites']
        self.localData[cpu]['cca'] = idCCA
        del data['mode']
        data['cca'] = idCCA
        
    def reformatData(self,data):
        source = data['source']
        data['cca'] = self.localData[source]['cca']
        data['mode'] = self.localData[source]['mode']
        data['sequence'] = self.localData[source]['sequence']
        data['explication'] = self.localData[source]['activites']
        if data['faisable']:
            data['sequence'] += data['explication']
                
class BlockingCommunication(Communication):
    def __init__(self):
        pass
    
    def postMessage(self,data):
        assert(MPI.COMM_WORLD.Get_rank()>0)
        MPI.COMM_WORLD.send(data,dest=0)
        
    def readMessages(self):
        assert(MPI.COMM_WORLD.Get_rank()==0)
        for i in range(MPI.COMM_WORLD.Get_size()-1):
            data = MPI.COMM_WORLD.recv()
            yield data
    
    def __str__(self):
        return "Blocking communication"
        
        
        
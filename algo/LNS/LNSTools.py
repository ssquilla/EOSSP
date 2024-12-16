#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 11 10:29:31 2023

@author: ssquilla
"""
from mpi4py import MPI

from ..Utils.config import *
global config
config = Config()
instance = config.instance


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
        
        
        
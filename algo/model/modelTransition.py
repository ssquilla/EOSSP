#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 17 14:00:35 2023

@author: ssquilla
"""
import os
import math
import numpy as np
import sys
from mpi4py import MPI
import unittest
import line_profiler


def arrondiPlusProche(x,digits):
    shift = x*(10**digits)
    integer = int(math.floor(shift))
    if (x-integer)<(integer+1-x):
        return integer/10**digits
    else:
        return (integer+1)/10**digits


MAX_TIME_SLOTS = 6

class modelTransitionTimeIndependent:
    def __init__(self,vlat,vlong,tau_max,nom):
        self.vlat = vlat
        self.vlong = vlong
        self.tau_max = tau_max
        self.nom = nom
    
    def utiliseSolverTimeDep(self):
        return False
    
    def getVitesseLatitude(self):
        return self.vlat
    
    def getVitesseLongitude(self):
        return self.vlong
    
    def getMajorantDureeTransition(self):
        return self.tau_max
    
    def __str__(self):
        msg = "Modèle de transition \"" + self.nom +"\" (vlat="+str(self.vlat)
        msg += " ; vlong="+str(self.vlong)
        msg += " ; tau_max="+str(self.tau_max)+")"
        return msg
    
    def estTimeDependent(self):
        return False
    
    def getTransition(self,satellite,p1,p2,digits):
        return self.calculTransition(satellite,p1,p2,digits)
        
    def calculTransition(self,satellite,p1,p2,digits):
        if(p1==p2):
            return 0
        if(satellite.estVidage(p1)):
            coord1 = satellite.getVidages()[p1].getCoordonnees()
        else:
            coord1 = satellite.getObservations()[p1].getCoordonnees()
        if(satellite.estVidage(p2)):
            coord2 = satellite.getVidages()[p2].getCoordonnees()
        else:
            coord2 = satellite.getObservations()[p2].getCoordonnees()
        vlat = self.getVitesseLatitude()
        vlong = self.getVitesseLongitude()
        tau_max = self.getMajorantDureeTransition()
            
        dlat = (coord1[0]-coord2[0])/vlat
        dlong = (coord1[1]-coord2[1])/vlong
        res = round(min(math.sqrt(dlat**2 + dlong**2)  ,tau_max),digits)
        return  res    


class modeleTimeIndependentExtrait:
    def __init__(self,name,componentsFolder,modeleTimeDep):
        self.name = name
        self.tau_max = modeleTimeDep.tau_max
        
        rep = ""
        lst = modeleTimeDep.getComponentsPath().split("/")
        for i,elmt in enumerate(lst):
            if i<len(lst)-1:
                rep += elmt + "/"
        rep += componentsFolder
        self.componentsFolder = rep
        self.mapping_cca = modeleTimeDep.mapping_cca
        self.indice_start_a = modeleTimeDep.indice_start_a
        self.rang = modeleTimeDep.rang
        self.start_cca_approx = modeleTimeDep.start_cca_approx
        self.taille_approx = modeleTimeDep.taille_approx

    def __str__(self):
        return self.name
        
    def estTimeDependent(self):
        return False    
    
    def getMajorantDureeTransition(self):
        return self.tau_max 
    
    def getComponentsPath(self):
        return self.componentsFolder
    
    def utiliseSolverTimeDep(self):
        return True
    
    def getAbsoluteIndexOfPairApprox(self,a1,a2)  :
        return self.start_cca_approx[self.mapping_cca[a1]] + self.getIndexOfPairInCCA(a1,a2)  
    
    def getIndexOfPairInCCA(self,a1,a2):
        cca1 = self.mapping_cca[a1]
        cca2 = self.mapping_cca[a2]
        if cca1!=cca2:
            return -1
        if a1<a2:
            return self.indice_start_a[a1] + self.rang[a2] -1
        elif a1==a2:
            return -1
        else:
            return self.indice_start_a[a1] + self.rang[a2]    
    
class modeleMoyenExtrait(modeleTimeIndependentExtrait):
    def __init__(self,modeleTimeDep):
        super().__init__("Modèle de transition moyen extrait du modèle time-dependent","components_mean",modeleTimeDep)
        self.tabular = modeleTimeDep.approximation_moyenne
        
    def getTransition(self, satellite, a1, a2, digits):
        if a1==a2:
            return 0
        index = self.getAbsoluteIndexOfPairApprox(a1,a2)
        return round(self.tabular[index],digits)
    
class modeleOptimisteExtrait(modeleTimeIndependentExtrait):
    def __init__(self,modeleTimeDep):
        super().__init__("Modèle de transition optimiste extrait du modèle time-dependent","components_optimistic",modeleTimeDep)
        self.tabular = modeleTimeDep.approximation_min
        
    def getTransition(self, satellite, a1, a2, digits):
        if a1==a2:
            return 0
        index = self.getAbsoluteIndexOfPairApprox(a1,a2)
        return round(self.tabular[index],digits)
    
class modelePessimisteExtrait(modeleTimeIndependentExtrait):
    def __init__(self,modeleTimeDep):
        self.comp_suffixe = "components_max"
        super().__init__("Modèle de transition pessismiste extrait du modèle time-dependent","components_pessimistic",modeleTimeDep)
        self.tabular = modeleTimeDep.approximation_max
        
    def getTransition(self, satellite, a1, a2, digits):
        if a1==a2:
            return 0
        index = self.getAbsoluteIndexOfPairApprox(a1,a2)
        return round(self.tabular[index],digits)


profiler = line_profiler.LineProfiler()
     
class modelTransitionTimeDependent:
    
    def __init__(self,foldername,constellation,composantes,precompute_approx=None,profile=False):
        
        self.createTransitionModelFromFiles(foldername,constellation,composantes,precompute_approx=None)
        if profile:
            profiler.print_stats()
            
    def createTransitionModelFromFiles(self,foldername,constellation,composantes,precompute_approx=None):
        assert(precompute_approx is None or precompute_approx in ["slow","mean","fast"])
        self.precompute_approx = precompute_approx
        self.tau_max = 90
        self.tabular = {}
        self.echelle_fichier_composante = 1e-3
        self.foldername = foldername
        digits = -int(math.log10(self.echelle_fichier_composante))
        mapping_files = self.mapCCAToCores(foldername)
        self.initComposantesData(composantes) 
        self.initSharedStructures(constellation)          
        
        size = constellation.getIdActivitesMax()+1
        self.start_window = np.zeros(size)
        self.end_window = np.zeros(size)
        durees = np.zeros(size)
        
        for s in constellation.getSatellites():
            for a in constellation.getSatellite(s).getActivites():
                self.end_window[a] = constellation.getSatellite(s).getActivite(a).getFin()
                self.start_window[a] = constellation.getSatellite(s).getActivite(a).getDebut()
                
        for filename in mapping_files:
            # infos sur la cca : nom du fichier
            split = filename.split(".comp")[0].split("_")
            id_cca = int(split[1]),int(split[2])
            activites = []
            #durees = {}
            # lecture du fichier
            with open(foldername+'/'+filename,'r') as file:
                # information sur chaque fenetre
                nObs = int(file.readline())
                for i in range(nObs):
                    line = file.readline().split('\n')[0].split(" ")
                    id_act,debut = int(line[0]),int(line[2])
                    if id_act != -1:
                        self.end_window[id_act] = constellation.getActivite(id_act).getFin()
                        self.startsMapping[id_act] = int(debut*self.echelle_fichier_composante)
                        durees[id_act] = int(line[1])*self.echelle_fichier_composante
                        self.NtimeSlotsMapping[id_act] = int(line[6])
                        self.dureeTimeSlotsMapping[id_act] = int(line[7])*self.echelle_fichier_composante
                    activites.append((id_act,int(debut*self.echelle_fichier_composante)))
                # transition de la fenetre
                for a1,debut in activites:
                    nValues = self.NtimeSlotsMapping[a1]
                    timeStep = self.dureeTimeSlotsMapping[a1]
                    dt = self.dureeTimeSlotsMapping[a1]
                    #assert(nValues)<=MAX_TIME_SLOTS
                    # information pour chaque voisin
                    for k,(a2,debut) in enumerate(activites):
                        values = file.readline().split("\n")[0].split(" ")
                        # créer le tableau de valeurs pour le couple (a1,a2)
                        if a1 != -1 and a2 != -1 and a1 in self.mapping_cca and a2 in self.mapping_cca:
                            indice_start_a1_a2 = self.getAbsoluteIndexOfPairTimeDep(a1,a2)
                            if a1!=a2:
                                index_approx = self.getAbsoluteIndexOfPairApprox(a1,a2)
                        min_t,max_t,mean_t = np.Inf,0,0
                        for i in range(nValues): # verifier que les activites sont dans le mapping (car certaines obs peuvent etre ignorees. systematiques par exemple)
                            if a1 != -1 and a2 != -1 and a1 in self.mapping_cca and a2 in self.mapping_cca:
                                value = int(values[i])*self.echelle_fichier_composante
                                valeur = arrondiPlusProche((value-durees[a1]),digits)
                                if a1!=a2 and valeur/dt<-1:
                                    self.couplesFIFO[index_approx] = False
                                if not(valeur>=0):
                                    valeur = 0
                                if a1!=a2:
                                    self.vect_time_dep[indice_start_a1_a2+i]=valeur
                                min_t = min(min_t,valeur)
                                mean_t += valeur
                                max_t = max(max_t,valeur)
                        # creation des modeles approches
                        if a1 != -1 and a2 != -1 and a1!=a2 and a1 in self.mapping_cca and a2 in self.mapping_cca:
                            if np.isinf(min_t):
                                min_t = 0
                            mean_t = arrondiPlusProche(mean_t / nValues,3)
                            index = self.getAbsoluteIndexOfPairApprox(a1,a2)
                            self.approximation_min[index_approx] =  min_t
                            self.approximation_max[index_approx] = max_t
                            self.approximation_moyenne[index_approx] = mean_t
        MPI.COMM_WORLD.Barrier()
        self.proportionCouplesFIFO = 100*sum(self.couplesFIFO)/self.taille_approx
    
    
    #@profiler
    def createTransitionModelFromFilesV2(self,foldername,constellation,composantes,precompute_approx=None):
        assert(precompute_approx is None or precompute_approx in ["slow","mean","fast"])
        self.precompute_approx = precompute_approx
        self.tau_max = 90
        self.tabular = {}
        self.echelle_fichier_composante = 1e-3
        self.foldername = foldername
        digits = -int(math.log10(self.echelle_fichier_composante))
        mapping_files = self.mapCCAToCores(foldername)
        self.initComposantesData(composantes) 
        self.initSharedStructures(constellation)          
        
        self.start_window = {}
        self.end_window = {}
        for s in constellation.getSatellites():
            for a in constellation.getSatellite(s).getActivites():
                self.end_window[a] = constellation.getSatellite(s).getActivite(a).getFin()
                self.start_window[a] = constellation.getSatellite(s).getActivite(a).getDebut()
        
        for filename in mapping_files:
            split = filename.split(".comp")[0].split("_")
            id_cca = int(split[1]),int(split[2])
            activites = []
            durees = {}
            
            with open(foldername+'/'+filename,'r') as file:
                nObs = int(file.readline())
                for i in range(nObs):
                    line = file.readline().split('\n')[0].split(" ")
                    id_act,debut = int(line[0]),int(line[2])
                    if id_act != -1:
                        self.end_window[id_act] = constellation.getActivite(id_act).getFin()
                        self.startsMapping[id_act] = int(debut*self.echelle_fichier_composante)
                        durees[id_act] = int(line[1])*self.echelle_fichier_composante
                        self.NtimeSlotsMapping[id_act] = int(line[6])
                        self.dureeTimeSlotsMapping[id_act] = int(line[7])*self.echelle_fichier_composante
                    activites.append((id_act,int(debut*self.echelle_fichier_composante)))
                for a1,debut in activites:
                    nValues = self.NtimeSlotsMapping[a1]
                    timeStep = self.dureeTimeSlotsMapping[a1]
                    assert(nValues)<=MAX_TIME_SLOTS
                    for k,(a2,debut) in enumerate(activites):
                        values = file.readline().split("\n")[0].split(" ")
                        # créer le tableau de valeurs pour le couple (a1,a2)
                        if a1 != -1 and a2 != -1 and a1 in self.mapping_cca and a2 in self.mapping_cca:
                            indice_start_a1_a2 = self.getAbsoluteIndexOfPairTimeDep(a1,a2)
                            if not len(values)==nValues and (a1,a2)!=(-1,-1):
                                print(a1,a2,values,nValues)
                                assert(len(values)==nValues)
                        min_t,max_t,mean_t = np.Inf,0,0
                        for i in range(nValues): # verifier que les activites sont dans le mapping (car certaines obs peuvent etre ignorees. systematiques par exemple)
                            if a1 != -1 and a2 != -1 and a1 in self.mapping_cca and a2 in self.mapping_cca:
                                value = int(values[i])*self.echelle_fichier_composante
                                valeur = round((value-durees[a1]),digits)
                                dt = self.dureeTimeSlotsMapping[a1]
                                if a1!=a2 and valeur/dt<-1:
                                    index = self.getAbsoluteIndexOfPairApprox(a1,a2)
                                    self.couplesFIFO[index] = False
                                    #print(a1,a2,"pas FIFO :",valeur/dt)
                                if not(valeur>=0):
                                    valeur = 0
                                """ test """    
                                if indice_start_a1_a2+i>= self.taille_time_dep_vect:
                                    print(a1,a2,self.mapping_cca[a1])
                                    print(len(sorted(composantes.getActivitesComposante(self.mapping_cca[a1]))),sorted(composantes.getActivitesComposante(self.mapping_cca[a1])))
                                
                                """ """
                                if a1!=a2:
                                    assert(self.vect_time_dep[indice_start_a1_a2+i]==0)
                                    self.vect_time_dep[indice_start_a1_a2+i]=valeur
                                min_t = min(min_t,valeur)
                                mean_t += valeur
                                max_t = max(max_t,valeur)
                        if a1 != -1 and a2 != -1 and a1!=a2 and a1 in self.mapping_cca and a2 in self.mapping_cca:
                            if np.isinf(min_t):
                                min_t = 0
                            mean_t = round(mean_t / nValues,3)
                            assert(min_t>=0)
                            assert(max_t>=0)
                            assert(mean_t>=0)
                            index = self.getAbsoluteIndexOfPairApprox(a1,a2)
                            self.approximation_min[index] =  min_t
                            self.approximation_max[index] = max_t
                            self.approximation_moyenne[index] = mean_t
        MPI.COMM_WORLD.Barrier()
        self.proportionCouplesFIFO = 100*sum(self.couplesFIFO)/self.taille_approx

    def initComposantesData(self,composantes):
        cca_triees = sorted(composantes.getComposantes())
        # connaitre la cca de chaque activite
        self.mapping_cca = {}
        for cca in cca_triees:
            for a in composantes.getActivitesComposante(cca):
                self.mapping_cca[a] = cca
        # rang des activites dans chaque cca
        self.rang = {}
        for cca in cca_triees:
            rg = 0
            for a in sorted(composantes.getActivitesComposante(cca)):
                self.rang[a] = rg
                rg += 1
        # savoir quand demarre l'indice d'une cca dans le vecteur de transition
        self.start_cca_time_dep = {}
        self.start_cca_approx = {}
        start_time_dep = 0
        start_approx = 0
        for cca in cca_triees:
            taille = composantes.getTailleComposante(cca)
            self.start_cca_time_dep[cca] = start_time_dep
            self.start_cca_approx[cca] = start_approx
            # tailles -1 car on ne compte pas les couples redondant (a,a)
            start_time_dep += (taille-1)*taille*MAX_TIME_SLOTS
            start_approx += (taille-1)*taille
        # savoir l'indice de demarrage d'une activite
        self.indice_start_a = {}
        for cca in cca_triees:
            for a in composantes.getActivitesComposante(cca):
                self.indice_start_a[a] = self.rang[a]*(composantes.getTailleComposante(cca)-1)
        self.taille_time_dep_vect = start_time_dep
        self.taille_approx = start_approx
        # activites pas FIFO
        self.couplesFIFO = self.createSharedArrayOfBool(self.taille_approx,init_value=True)
        
   
    def isFIFO(self,a1,a2):
        index = self.getAbsoluteIndexOfPairApprox(a1,a2)
        return self.couplesFIFO[index]
    
    def getAbsoluteIndexOfPairApprox(self,a1,a2):
        if a1==a2:
            return -1
        return self.start_cca_approx[self.mapping_cca[a1]] + self.getIndexOfPairInCCA(a1,a2)          
   
    def getIndexOfPairInCCA(self,a1,a2):
        cca1 = self.mapping_cca[a1]
        cca2 = self.mapping_cca[a2]
        if cca1!=cca2:
            return -1
        if a1<a2:
            return self.indice_start_a[a1] + self.rang[a2] -1
        elif a1==a2:
            return -1
        else:
            return self.indice_start_a[a1] + self.rang[a2]
        
    def getAbsoluteIndexOfPairTimeDep(self,a1,a2):
        if a1==a2:
            return -1
        return self.start_cca_time_dep[self.mapping_cca[a1]] + self.getIndexOfPairInCCA(a1,a2)*MAX_TIME_SLOTS

    def createSharedArrayOfBool(self,size,init_value=False):
        comm = MPI.COMM_WORLD
        if comm.Get_size()>1:
            itemsize = MPI.BOOL.Get_size()
            
            if comm.Get_rank() == 0: 
                nbytes = size * itemsize 
            else: 
                nbytes = 0
            win = MPI.Win.Allocate_shared(nbytes, itemsize, comm=comm) 
            # create a numpy array whose data points to the shared mem
            buf, itemsize = win.Shared_query(0) 
            assert itemsize == MPI.BOOL.Get_size()
            arr = np.ndarray(buffer=buf, dtype='bool', shape=(size,)) 
        else:
            arr = np.ndarray( dtype='bool', shape=(size,)) 
        for i in range(arr.shape[0]):
            arr[i] = init_value
        return arr

    def createSharedArrayOfDouble(self,size):  
        comm = MPI.COMM_WORLD
        if comm.Get_size()>1:
            itemsize = MPI.DOUBLE.Get_size()
            if comm.Get_rank() == 0: 
                nbytes = size * itemsize 
            else: 
                nbytes = 0
            win = MPI.Win.Allocate_shared(nbytes, itemsize, comm=comm) 
            # create a numpy array whose data points to the shared mem
            buf, itemsize = win.Shared_query(0) 
            assert itemsize == MPI.DOUBLE.Get_size()
            arr = np.ndarray(buffer=buf, dtype='d', shape=(size,)) 
        else:
            arr = np.ndarray(dtype='d', shape=(size,))
        for i in range(arr.shape[0]):
            arr[i] = 0
        return arr
    
    def createSharedArrayOfInt(self,size,init_value=0):
        comm = MPI.COMM_WORLD
        if comm.Get_size()>1:
            itemsize = MPI.INTEGER.Get_size()
            comm = MPI.COMM_WORLD
            if comm.Get_rank() == 0: 
                nbytes = size * itemsize 
            else: 
                nbytes = 0
            win = MPI.Win.Allocate_shared(nbytes, itemsize, comm=comm) 
            # create a numpy array whose data points to the shared mem
            buf, itemsize = win.Shared_query(0) 
            assert itemsize == MPI.INTEGER.Get_size()   
            arr = np.ndarray(buffer=buf, dtype='int32', shape=(size,))
        else:
            arr = np.ndarray(dtype='int32', shape=(size,))
        for i in range(arr.shape[0]):
            arr[i] = init_value
        return arr
    
    def createSharedMatrixOfDouble(self,length):  
        comm = MPI.COMM_WORLD
        if comm.Get_size()>1:
            itemsize = MPI.DOUBLE.Get_size()
            if comm.Get_rank() == 0: 
                nbytes = length * length * itemsize 
            else: 
                nbytes = 0
            win = MPI.Win.Allocate_shared(nbytes, itemsize, comm=comm) 
            # create a numpy array whose data points to the shared mem
            buf, itemsize = win.Shared_query(0) 
            arr = np.ndarray(buffer=buf, dtype='d', shape=(length,length))
            for i in range(arr.shape[0]):
                for j in range(arr.shape[1]):
                    arr[i][j] = 0
        else:
            arr = np.zeros((length,length))
        return arr        
        
    def initSharedStructures(self,constellation):
        size = constellation.getIdActivitesMax()+1
        self.NtimeSlotsMapping = self.createSharedArrayOfInt(size)
        self.dureeTimeSlotsMapping = self.createSharedArrayOfDouble(size)
        self.startsMapping = self.createSharedArrayOfInt(size)
        self.vect_time_dep = self.createSharedArrayOfDouble(self.taille_time_dep_vect)
        self.approximation_min = self.createSharedArrayOfDouble(self.taille_approx)
        self.approximation_max = self.createSharedArrayOfDouble(self.taille_approx)
        self.approximation_moyenne = self.createSharedArrayOfDouble(self.taille_approx)
        
    def mapCCAToCores(self,foldername):
        files = sorted(os.listdir(foldername),key = lambda file : os.stat(foldername+'/'+file).st_size,reverse=True)
        mapping = []
        for i in range(MPI.COMM_WORLD.Get_rank(),len(files),MPI.COMM_WORLD.Get_size()):
            mapping.append(files[i])
        return mapping    
    
    def getComponentsPath(self):
        return self.foldername

    def getTransitionMin(self,a1,a2):
        index = self.getAbsoluteIndexOfPairApprox(a1,a2)
        return self.approximation_min[index]                     

    def utiliseSolverTimeDep(self):
        return True
    
    def getTransitionMax(self,a1,a2):
        index = self.getAbsoluteIndexOfPairApprox(a1,a2)
        return self.approximation_max[index] 
        
    def getTransitionMoyenne(self,a1,a2):
        index = self.getAbsoluteIndexOfPairApprox(a1,a2)
        return self.approximation_mean[index] 
    
    #def getComposantes(self):
    #    return list(self.tabular.keys())
    
    #def getActivitesCCA(self,id_cca):
    #    return list(self.tabular[id_cca].keys())
            
    def getPointsDeControle(self,a1,a2,time=True):
        if self.mapping_cca[a1]!=self.mapping_cca[a2]:
            return [(0,0)]
        points = []
        start = self.startsMapping[a1]
        nTimeSlot = self.NtimeSlotsMapping[a1]
        dureeTimeSlot = self.dureeTimeSlotsMapping[a1]
        start_index = self.getAbsoluteIndexOfPairTimeDep(a1,a2)
        if start_index==-1:
            return [(0,0)]
        for t in range(nTimeSlot):
            time = start+t*dureeTimeSlot
            points.append((time,self.vect_time_dep[start_index+t]))
        return points
        
    def trouverTempsProches(self,a1,a2,t):
        start_a1 = self.startsMapping[a1]
        if t<start_a1:
            return None
        NtimeSlot = self.NtimeSlotsMapping[a1]
        dureeTimeSlot = self.dureeTimeSlotsMapping[a1]
        end_a2 = start_a1 + dureeTimeSlot*NtimeSlot
        if t>end_a2:
            return None
        nSlots = (t-start_a1)//dureeTimeSlot
        assert(nSlots==int(nSlots))
        nSlots = int(nSlots)
        t_left = start_a1 + dureeTimeSlot*nSlots
        t_right = start_a1 + dureeTimeSlot*(nSlots+1)
        return nSlots,t_left,t_right
    
    def getTransition(self,a1,a2,t):
        if a1==a2:
            return 0
        start_index = self.getAbsoluteIndexOfPairTimeDep(a1,a2)
        if start_index==-1:
            return 0
        """ cas particulier :
            certaines obs manquent malheureusement d'une valeur dans leur fenetre 
            (erreur dans le generateur)
            sur le dernier segment on retourne la meme valeur que la precedente"""
        t_larger_than_max_value = t>self.start_window[a1]+self.dureeTimeSlotsMapping[a1]*(self.NtimeSlotsMapping[a1]-1)
        if t_larger_than_max_value and t <= self.end_window[a1]+1e-4:
            return self.vect_time_dep[start_index+self.NtimeSlotsMapping[a1]-1]
        
        res = self.trouverTempsProches(a1, a2, t)
        if res is None:
            return 0
        indice_left,t_left,t_right = res
        value_left = self.vect_time_dep[start_index+indice_left]
        if indice_left+1>=self.NtimeSlotsMapping[a1]:
            return value_left
        else:
            value_right = self.vect_time_dep[start_index+indice_left+1]        
        delta_t = (t-t_left)/(t_right-t_left)
        return delta_t*value_right + (1-delta_t)*value_left
    
    def getMajorantCoupleActivite(self,a1,a2):
        if a1==a2:
            return 0
        indice = self.getAbsoluteIndexOfPairApprox(a1,a2)
        return self.approximation_max[indice]
    
    def getMajorantDureeTransition(self):
        return self.tau_max
                
    def estTimeDependent(self):
        return True
    
    def __str__(self):
        return "Modèle de transition time-dependent."
    
global modeleLent
modeleLent = modelTransitionTimeIndependent(0.8,0.8,90,"slow")

global modeleRapide
modeleRapide = modelTransitionTimeIndependent(10,15,90,"fast")

global modeleMoyen
modeleMoyen = modelTransitionTimeIndependent(1,1,90,"mean")
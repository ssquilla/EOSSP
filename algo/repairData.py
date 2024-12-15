#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 15 15:52:52 2023

@author: ssquilla
"""
import sys, importlib
from pathlib import Path
from resource import *
from mpi4py import MPI
from math import *

def import_parents(level=1):
    global __package__
    file = Path(__file__).resolve()
    parent, top = file.parent, file.parents[level]
        
    sys.path.append(str(top))
    try:
        sys.path.remove(str(parent))
    except ValueError: # already removed
        pass
    
    __package__ = '.'.join(parent.parts[len(top.parts):])
    importlib.import_module(__package__) # won't be needed after that
        


if __name__ == '__main__' and __package__ is None:
    import_parents(level=1)
    
import os

from .model.constellation import *

ONE_SHOT_MONO = "ONE_SHOT_MONO"
LONG_MONO = "LONG_MONO"
ONE_SHOT_STEREO = "ONE_SHOT_STEREO"
SYSTEMATIC = "SYSTEMATIC"
PERIODIC = "PERIODIC"

path = "../data/time-dependent/"
def splitAndAddLine(file,lines,add):
    line = file.readline()
    split = line.split("\n")[0].split(",")
    if add:
        lines.append(split) 
    return line,split

def reconstruire(split):
    msg = ""
    for i,x in enumerate(split):
        msg += str(x)
        if i<len(split)-1:
            msg += ","
        else:
            msg += "\n"
    return msg

def reparerNumerotation(dossier):
    for root,dirs,files in os.walk(dossier):
        for filename in files:
            if filename.endswith(".pb"):
                lines = []
                lines_req = {}
                indices_a_changer = []
                print(filename)
                with open(root+'/'+filename,'r') as file:
                    # 1ere ligne : header fichier
                    lines_req[-1] = []
                    _,header = splitAndAddLine(file, lines_req[-1], True)
                    nReq = int(header[0])
                    for r in range(nReq):
                        obs_rencontrees = []
                        # header requete
                        lines_req[r] = []
                        _,header_r = splitAndAddLine(file, lines_req[r], True)
                        nObs,type_r = int(header_r[1]),header_r[2]
                        for o in range(nObs):
                            idObs,i_obs = None,None
                            line,split_o = splitAndAddLine(file,lines_req[r],False)
                            if type_r in [ONE_SHOT_STEREO,PERIODIC]:
                                idObs = int(split_o[1])
                                i_obs = 1
                            else:
                                idObs = int(split_o[0])
                                i_obs = 0
                            if idObs not in obs_rencontrees:
                                lines_req[r].append(split_o)
                                obs_rencontrees.append(idObs)
                        lines_req[r][0][1] = len(obs_rencontrees)
                    lines_req[nReq] = []
                    _,split = splitAndAddLine(file,lines_req[nReq],True)
                    nDownloads = int(split[0])
                    for d in range(nDownloads):
                        _,split = splitAndAddLine(file,lines_req[nReq],True)
                        lines_req[nReq][-1][0] = d
                    
                    i = nDownloads
                    for r in range(0,nReq):
                        header_r = lines_req[r][0]
                        nObs,type_r = int(header_r[1]),header_r[2]
                        for o in range(1,len(lines_req[r])):
                            idObs,i_obs = None,None
                            if type_r in [ONE_SHOT_STEREO,PERIODIC]:
                                idObs = int(split_o[1])
                                i_obs = 1
                            else:
                                idObs = int(split_o[0])
                                i_obs = 0
                            lines_req[r][o][i_obs] = i
                            i += 1
                        
                with open(root+'/'+filename,'w') as file:
                    for r in lines_req:
                        for o in lines_req[r]:
                            tmp = reconstruire(o)
                            file.write(tmp)


def retirerInfosInutilesObs(info_obs,offset):
    assert(len(info_obs)==offset+11)
    info_obs.pop(6+offset)
    info_obs.pop(5+offset)
                
def effacerElements(dossier):
    for root,dirs,files in os.walk(dossier):
        for filename in files:
            if filename.endswith(".pb"):
                lines = []
                lines_req = {}
                indices_a_changer = []
                print(filename)
                with open(root+'/'+filename,'r') as file:
                    # 1ere ligne : header fichier
                    lines_req[-1] = []
                    _,header = splitAndAddLine(file, lines_req[-1], True)
                    if(len(lines_req[-1][0])==2):
                        lines_req[-1][0].pop(1)
                        nReq = int(header[0])
                        for r in range(nReq):
                            obs_rencontrees = []
                            # header requete
                            lines_req[r] = []
                            _,header_r = splitAndAddLine(file, lines_req[r], False)
                            header_r.pop(2)
                            lines_req[r].append(header_r)
                            nObs,type_r = int(header_r[1]),header_r[2]
                            for o in range(nObs):
                                idObs,i_obs = None,None
                                line,split_o = splitAndAddLine(file,lines_req[r],False)
                                if type_r in [ONE_SHOT_STEREO,PERIODIC]:
                                    idObs = int(split_o[1])
                                    i_obs = 1
                                else:
                                    idObs = int(split_o[0])
                                    i_obs = 0
                                if idObs not in obs_rencontrees:
                                    lines_req[r].append(split_o)
                                    obs_rencontrees.append(idObs)
                                    retirerInfosInutilesObs(lines_req[r][-1],i_obs)
                            lines_req[r][0][1] = len(obs_rencontrees)
                        lines_req[nReq] = []
                        _,split = splitAndAddLine(file,lines_req[nReq],True)
                        nDownloads = int(split[0])
                        for d in range(nDownloads):
                            _,split = splitAndAddLine(file,lines_req[nReq],True)
                            lines_req[nReq][-1][0] = d
                        
                        i = nDownloads
                        for r in range(0,nReq):
                            header_r = lines_req[r][0]
                            nObs,type_r = int(header_r[1]),header_r[2]
                            for o in range(1,len(lines_req[r])):
                                idObs,i_obs = None,None
                                if type_r in [ONE_SHOT_STEREO,PERIODIC]:
                                    idObs = int(split_o[1])
                                    i_obs = 1
                                else:
                                    idObs = int(split_o[0])
                                    i_obs = 0
                                lines_req[r][o][i_obs] = i
                                i += 1
                            
                    with open(root+'/'+filename,'w') as file:
                        for r in lines_req:
                            for o in lines_req[r]:
                                tmp = reconstruire(o)
                                file.write(tmp)


def retirerSystematiques(dossier):
    for root,dirs,files in os.walk(dossier):
        for filename in files:
            if filename.endswith(".pb"):
                lines = []
                lines_req = {}
                indices_a_changer = []
                compteur_req = -1
                print(filename,"retrait des systematiques")
                with open(root+'/'+filename,'r') as file:
                    # 1ere ligne : header fichier
                    lines_req[-1] = []
                    _,header = splitAndAddLine(file, lines_req[-1], True)
                    nReq = int(header[0])
                    for r in range(nReq):
                        obs_rencontrees = []
                        # header requete
                        
                        line,header_r = splitAndAddLine(file, lines_req[compteur_req], False)
                        nObs,type_r = int(header_r[1]),header_r[2]
                        if type_r in [ONE_SHOT_MONO,LONG_MONO,ONE_SHOT_STEREO,PERIODIC]:
                            compteur_req += 1
                            lines_req[compteur_req] = []
                            header_r[0] = compteur_req
                            lines_req[compteur_req].append(header_r)
                        
                        for o in range(nObs):
                            idObs,i_obs = None,None
                            line,split_o = splitAndAddLine(file,lines_req[compteur_req],False)
                            if type_r in [ONE_SHOT_MONO,LONG_MONO,ONE_SHOT_STEREO,PERIODIC]:
                                if type_r in [ONE_SHOT_STEREO,PERIODIC]:
                                    idObs = int(split_o[1])
                                    i_obs = 1
                                else:
                                    idObs = int(split_o[0])
                                    i_obs = 0
                                if idObs not in obs_rencontrees:
                                    lines_req[compteur_req].append(split_o)
                                    obs_rencontrees.append(idObs)
                        if type_r in [ONE_SHOT_MONO,LONG_MONO,ONE_SHOT_STEREO,PERIODIC]:
                            lines_req[compteur_req][0][1] = len(obs_rencontrees)
                    lines_req[-1][0] = len(list(lines_req.keys()))-1 # indiquer nb requetes
                        
                    lines_req[nReq] = []
                    _,split = splitAndAddLine(file,lines_req[nReq],True)
                    nDownloads = int(split[0])
                    for d in range(nDownloads):
                        _,split = splitAndAddLine(file,lines_req[nReq],True)
                        lines_req[nReq][-1][0] = d
                
                split_instance = filename.split("_")
                repartition = split_instance.pop(0).split("-")
                repartition.pop(3) # retirer systematique du nom
                tmp = ""
                for i,x in enumerate(repartition):
                    tmp += x
                    if i < len(repartition)-1:
                        tmp += "-"
                repartition = tmp
                print(repartition)
                nom = repartition
                for elmt in split_instance:
                    nom += "_" + elmt
                with open(root+'/'+nom,'w') as file:
                    for r in lines_req:
                        if r == -1:
                            file.write(reconstruire(lines_req[-1]))
                        else:
                            for o in lines_req[r]:
                                tmp = reconstruire(o)
                                file.write(tmp)
                            
                os.unlink(root+'/'+filename) 

path_dep = "../data/time-dependent/"
path_indep = "../data/time-independent/"
home = "../../data_repairing"
def ajouterExtension(dossier):
    for root,dirs,files in os.walk(dossier):
        for filename in files:
            if "comp" in filename:
                if not filename.endswith(".comp"):
                    os.rename(root+"/"+filename,root+"/"+filename+".comp")
            else:
                if not filename.endswith(".pb"):
                    os.rename(root+"/"+filename,root+"/"+filename+".pb")


class modelTransitionTimeDependent:
    def __init__(self,foldername,constellation,precompute_approx=None):
        assert(precompute_approx is None or precompute_approx in ["slow","mean","fast"])
        self.precompute_approx = precompute_approx
        self.tau_max = 90
        self.tabular = {}
        self.echelle_fichier_composante = 1e-3
        self.foldername = foldername
        mapping_files = self.mapCCAToCores(foldername)        
        for filename in mapping_files:
            print(filename)
            split = filename.split(".comp")[0].split("_")
            id_cca = int(split[1]),int(split[2])
            activites = []
            durees = {}
            write_lines = []
            count_values = {}
            with open(foldername+'/'+filename,'r') as file:
                line = file.readline()
                write_lines.append(line)
                nObs = int(line)
                for i in range(nObs):
                    line = file.readline()
                    write_lines.append(line)
                    splitted = line.split('\n')[0].split(" ")
                    id_act = int(splitted[0])
                    count_values[id_act] = int(splitted[6])
                    activites.append((id_act,0))
                for a1,debut in activites:
                    nValues = count_values[a1]
                    for k,(a2,debut) in enumerate(activites):
                        line = ""
                        if a1 != -1 :
                            duree_obs = constellation.getSatellite(id_cca[0]).getActivite(a1).getDuree()
                        else:
                            duree_obs = 0
                        values = file.readline().split("\n")[0].split(" ")
                        # créer le tableau de valeurs pour le couple (a1,a2)
                        for i in range(nValues): # verifier que les activites sont dans le mapping (car certaines obs peuvent etre ignorees. systematiques par exemple)
                                value = int(values[i])
                                if a1!=a2 and value==0:
                                    line += str(int(duree_obs*1000))
                                else:
                                    line += str(int(values[i]))
                                if i< nValues-1:
                                    line += " "
                        line += '\n'
                        write_lines.append(line)
            with open(foldername+'/'+filename,'w') as file:
                for line in write_lines:
                    file.write(line)
        MPI.COMM_WORLD.Barrier()

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
        
   
    
    def getAbsoluteIndexOfPairApprox(self,a1,a2)  :
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
        itemsize = MPI.BOOL.Get_size()
        comm = MPI.COMM_WORLD
        if comm.Get_rank() == 0: 
            nbytes = size * itemsize 
        else: 
            nbytes = 0
        win = MPI.Win.Allocate_shared(nbytes, itemsize, comm=comm) 
        # create a numpy array whose data points to the shared mem
        buf, itemsize = win.Shared_query(0) 
        assert itemsize == MPI.BOOL.Get_size()
        arr = np.ndarray(buffer=buf, dtype='bool', shape=(size,)) 
        for i in range(arr.shape[0]):
            arr[i] = init_value
        return arr

    def createSharedArrayOfDouble(self,size):  
        itemsize = MPI.DOUBLE.Get_size()
        comm = MPI.COMM_WORLD
        if comm.Get_rank() == 0: 
            nbytes = size * itemsize 
        else: 
            nbytes = 0
        win = MPI.Win.Allocate_shared(nbytes, itemsize, comm=comm) 
        # create a numpy array whose data points to the shared mem
        buf, itemsize = win.Shared_query(0) 
        assert itemsize == MPI.DOUBLE.Get_size()
        arr = np.ndarray(buffer=buf, dtype='d', shape=(size,)) 
        for i in range(arr.shape[0]):
            arr[i] = 0
        return arr
    
    def createSharedArrayOfInt(self,size,init_value=0):  
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
        for i in range(arr.shape[0]):
            arr[i] = init_value
        return arr
    
    def createSharedMatrixOfList(self,length):  
        itemsize = sys.getsizeof([])
        comm = MPI.COMM_WORLD
        if comm.Get_rank() == 0: 
            nbytes = length * length * itemsize 
        else: 
            nbytes = 0
        win = MPI.Win.Allocate_shared(nbytes, itemsize, comm=comm) 
        # create a numpy array whose data points to the shared mem
        buf, itemsize = win.Shared_query(0) 
        arr = np.ndarray(buffer=buf, dtype='object', shape=(length,length))
        for i in range(arr.shape[0]):
            for j in range(arr.shape[1]):
                arr[i][j] = None
        return arr
    
    def createSharedMatrixOfDouble(self,length):  
        itemsize = MPI.DOUBLE.Get_size()
        comm = MPI.COMM_WORLD
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
def reparerComposantes(dossier):
    for root,dirs,files in os.walk(dossier):
        for filename in files:
            if filename.endswith(".pb"):
                print(filename+"...")
                constellation = Constellation(root+"/"+filename,0,True)
                modelTransitionTimeDependent(root+"/components_mean",constellation)
                modelTransitionTimeDependent(root+"/components_optimistic",constellation)
                modelTransitionTimeDependent(root+"/components_pessimistic",constellation)
                modelTransitionTimeDependent(root+"/components_time_dep",constellation)
#reparerComposantes(home)
root = "../../data_repairing/instances_CPAIOR/"

#ajouterExtension(path_dep)
retirerSystematiques(root)
reparerNumerotation(root)  
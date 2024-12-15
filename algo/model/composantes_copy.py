import collections
import bisect
from mpi4py import MPI
from .graphes import *

from ..Utils import *
from ..Utils.config import *
config = Config()

from .modelTransition import modeleMoyen

class EvenementActivite: # event = (0:debut|1:fin)
    def __init__(self,sat,date,activite,event):
        self.satellite = sat
        self.activite = activite
        self.event = event
        self.date = date
        
    def getActivite(self):
        return self.activite
    
    def getEvent(self):
        return self.event
    
    def getDate(self):
        return self.date
    
    def getSatellite(self):
        return self.satellite
    
    def __lt__(self, other):
        return (self.satellite,self.date,self.event,self.activite) < (self.satellite,other.date,other.event,other.activite)
    
    def __gt__(self, other):
        return (self.satellite,self.date,self.event,self.activite) > (self.satellite,other.date,other.event,other.activite)
    
    def __eq__(self, other):
        return (self.satellite,self.date,self.event,self.activite) == (self.satellite,other.date,other.event,other.activite)
    
    def __le__(self,other):
        return (self.satellite,self.date,self.event,self.activite) <= (self.satellite,other.date,other.event,other.activite)
    
    def __ge__(self,other):
        return (self.satellite,self.date,self.event,self.activite) >= (self.satellite,other.date,other.event,other.activite)
        
class Tache:
    def __init__(self,cca,liste_modes):
        self.cca = cca
        self.liste_modes = liste_modes
    
    def __str__(self):
        return 'Tache(cca:'+str(self.cca)+',modes:'+str(self.liste_modes)+')'
    
    def getCCA(self):
        return self.cca
    
    def getListeModes(self):
        return self.liste_modes

def safeMergeDict(d1,d2):
    for key in d1:
        if key in d2:
            printColor(d1,d2,key,c='r')
        assert(key not in d2)
    for key in d2:
        assert(key not in d1)
    return {**d1,**d2}
    
class QuickUnion(object):
    def __init__(self, indices):
        self.lst = {a : a for a in indices}

    def find(self, ind):
        while ind != self.lst[ind]:
            ind = self.lst[ind]
        return ind

    def connect(self, p, q):
        return self.find(p) == self.find(q)

    def union(self, p, q):
        pid = self.find(p)
        self.lst[pid] = self.find(q)
    
    def add(self,p):
        self.lst[p] = p
        
    def contain(self,x):
        return x in self.lst
    
    # les composantes de other doivent être indépendante (aucune verification n'est effectuée ici)
    # utile pour calculer des groupes en parallele dont on sait a priori l'independance (satellites)
    def fusionnerQuickUnion(self,other):
        self.lst = safeMergeDict(self.lst,other.lst)
        
    def ecraserQuickUnion(self,other):
        self.lst = other.lst
    
class GroupeComposantesActivitesStatiques:
    class Composante:
        def __init__(self,id_cca,tau_max):
            self.tau_max = tau_max
            self.id = id_cca
            self.activites = []
            self.debut = np.Inf
            self.fin = -np.Inf
            self.satellite = id_cca[0]
        def estVide(self):
            return len(self.activites)==0
        def getDebut(self):
            return self.debut
        def getFin(self):
            return self.fin
        def __str__(self):
            return "Composante connexe\n| taille : "+str(len(self.activites))+"\n| début : "+str(self.debut)+"\n| fin : " + str(self.fin)
        def __contains__(self,elmt):
            if self.satellite!=elmt.getSat():
                return False
            res = intersect((elmt.getDebut()-self.tau_max,elmt.getFin()+self.tau_max),(self.debut,self.fin))
            return res
        def ajouterActivite(self,activite):
            self.activites.append(activite.getId())
            self.debut = min(self.debut,activite.getDebut())
            self.fin = max(self.fin,activite.getFin())
        def ajouterListeActivites(self,activites):
            for a in activites:
                self.ajouterActivite(a)
        def fusionnerComposante(self,other):
            if self.getDebut() < other.getDebut():
                assert(self.getFin() + self.tau_max>=other.getDebut())
                self.fin = other.getFin()
            else:
                assert(other.getFin() + self.tau_max<=self.getDebut())
                self.debut = other.getDebut()
            self.activites += other.activites
        def getElements(self):
            return self.activites
        def getId(self):
            return self.id
            
    def __init__(self,constellation,activites,modeleDeTransition,parallel=True):
        self.tau_max = modeleDeTransition.getMajorantDureeTransition()
        self.activite_cca = {}
        if parallel:
            # on parallelise le calcul des composantes
            satellites_attribues = self.attribuerSatellites(activites)
            self.calculerComposantes({s : activites.get(s,[]) for s in satellites_attribues},constellation)
            self.fusionnerComposantes()
        else:
            self.calculerComposantes(activites,constellation)
   
    def ecraser(self,other):
        self.sat_cca = other.sat_cca
        #self.cca_sat = other.cca_sat
        self.tailles = other.tailles
        self.activite_cca = other.activite_cca
        self.composantes = other.composantes
                
    def fusionnerComposantes(self):
        MPI.COMM_WORLD.Barrier()
        if MPI.COMM_WORLD.Get_rank()>0:
            data = {MPI.COMM_WORLD.Get_rank():self}
            MPI.COMM_WORLD.send(data,dest=0)
            data = MPI.COMM_WORLD.recv(source=0)
        else:
            for i in range(1,MPI.COMM_WORLD.Get_size()):
                data = MPI.COMM_WORLD.recv()
                groupe_a_fusionner = list(data.values())[0]
                #print(groupe_a_fusionner.composantes)
                self.fusionnerGroupe(groupe_a_fusionner)
            data = {'':self}
            for i in range(1,MPI.COMM_WORLD.Get_size()):
                MPI.COMM_WORLD.send(data,dest=i)
        MPI.COMM_WORLD.Barrier()
        if MPI.COMM_WORLD.Get_rank()>0:
            groupe_complet = list(data.values())[0]
            self.ecraser(groupe_complet)
        MPI.COMM_WORLD.Barrier()
    def ajouterRequete(self,constellation,r):
        for (s,a) in constellation.getRequete(r).getCouples():
            activite = constellation.getActivite(a,s=s)
            self.ajouterNouvelleActivite(activite)
    def ajouterNouvelleActivite(self,activite):
        composantes_trouvees = []
        for c in self.composantes:
            if activite in self.composantes[c]:
                if len(composantes_trouvees)==0:
                    self.ajouterActiviteComposante(activite,self.composantes[c])
                composantes_trouvees.append(c)
                
        if len(composantes_trouvees)==0:
            s = activite.getSat()
            indice = max(self.getCCASatellite(s))+1
            nouvelle_composante = self.Composante((s,indice),self.tau_max)
            self.tailles[(s,indice)] = 1
            self.ajouterActiviteComposante(activite,nouvelle_composante)
            self.composantes[(s,indice)] = nouvelle_composante
            self.sat_cca[s].append(indice)
        else:
            self.fusionnerListeComposantes(composantes_trouvees)
            
    def fusionnerListeComposantes(self,liste_composantes):
        assert(len(liste_composantes)>0)
        if len(liste_composantes)==1:
            return
        else:
            head = liste_composantes.pop()
            s_head = head[0]
            for other in liste_composantes:
                self.composantes[head].fusionnerComposante(self.composantes[other])
                s,cca = other
                assert(s==s_head)
                self.sat_cca[s].remove(cca)
                self.tailles[head]+= self.tailles[(s,cca)]
                del self.tailles[(s,cca)]
            for a in self.activite_cca:
                if self.activite_cca[a] in liste_composantes:
                    self.activite_cca[a] = head
            
    def ajouterActiviteComposante(self,activite,composante):
        assert(activite.getId() not in self.activite_cca)
        self.activite_cca[activite.getId()] = composante.getId()
        composante.ajouterActivite(activite)
        if composante.getId() not in self.tailles:
            self.tailles[composante.getId()] = 0
        self.tailles[composante.getId()] += 1
        
        
    def calculerComposantes(self,activites,constellation):
        self.composantes = {}
        self.tailles = {}
        for s in activites:
            compteur_cca = 0
            composanteCourante = self.Composante((s,compteur_cca),self.tau_max)
            compteur_cca += 1
            for i,a in enumerate(sorted(activites[s],key = lambda a : constellation.getActivite(a).getDebut())):
                activite = constellation.getActivite(a)
                if not (composanteCourante.estVide() or activite in composanteCourante):
                    self.composantes[composanteCourante.getId()] = composanteCourante
                    composanteCourante = self.Composante((s,compteur_cca),self.tau_max)
                    compteur_cca += 1
                self.ajouterActiviteComposante(activite,composanteCourante)
            if composanteCourante.getId() not in self.composantes and not composanteCourante.estVide():
                self.composantes[composanteCourante.getId()] = composanteCourante
        #self.cca_sat = {}
        self.sat_cca = {}
        for a in self.activite_cca:
            s,cca = self.activite_cca[a]
            #self.cca_sat[cca] = s
            if s not in self.sat_cca:
                self.sat_cca[s] = []
            self.sat_cca[s].append(cca)
        
    def attribuerSatellites(self,activites):
        size = MPI.COMM_WORLD.Get_size()
        rank = MPI.COMM_WORLD.Get_rank()
        tailles_triees = sorted([s for s in activites],key=lambda s : len(activites[s]),reverse=True)
        satellites_attribues = [i for i in tailles_triees if (i%size)==rank]
        return satellites_attribues
        
    def fusionnerGroupe(self,other):
        self.activite_cca = safeMergeDict(self.activite_cca,other.activite_cca)
        self.sat_cca = safeMergeDict(self.sat_cca,other.sat_cca)
        #self.cca_sat = safeMergeDict(self.cca_sat,other.cca_sat)
        self.tailles = safeMergeDict(self.tailles,other.tailles)
        self.composantes = safeMergeDict(self.composantes,other.composantes)
            
    def fusionnerListeGroupes(self,others):
        for o in others:
            self.fusionnerGroupe(o)
        
    def getCCASatellite(self,s):
        return self.sat_cca[s]
    """    
    def getSatelliteCCA(self,cca):
        return self.cca_sat[cca]
    """
    def getActiviteCCA(self,a):
        return self.activite_cca[a]

    def getTailleComposante(self,cca):
        return self.tailles[cca]

    def getNoeuds(self):
        return self.activite_cca
    def getActivitesComposante(self,cca):
        return self.composantes[cca].getElements()
    def getComposanteConnexe(self,cca):
        return self.composantes[cca]
    
    def getNombreComposantes(self):
        return len(list(self.tailles.keys()))
    
    def getNombreElements(self):
        return len(list(self.activite_cca.keys()))

    def getComposantes(self):
        return list(self.composantes.keys())
    
    def __str__(self):
        min_taille = min([self.getTailleComposante(cca) for cca in self.composantes])
        max_taille = max([self.getTailleComposante(cca) for cca in self.composantes])
        mean_taille = np.mean([self.getTailleComposante(cca) for cca in self.composantes])
        std_taille = np.std([self.getTailleComposante(cca) for cca in self.composantes])
        mess = "Composantes d'activités : processus "+str(MPI.COMM_WORLD.Get_rank())+"\n"
        mess += "| nombre de composantes : "+str(len(list(self.tailles.keys())))+'\n'
        mess += "| nombre d'activités : "+str(len(list(self.activite_cca.keys())))+'\n'
        mess += "| taille de la composante minimale : "+str(min_taille)+'\n'
        mess += "| taille de la composante maximale : " + str(max_taille)+'\n'
        mess += "| taille moyenne des composantes : " + str(mean_taille)+'\n'
        mess += "| variance de la tailles des composantes :" + str(std_taille)+'\n'
        #mess += str(sorted([self.getTailleComposante(cca) for cca in self.composantes])) + '\n'
        return mess

class ComposantesStatiquesPrecalculees(GroupeComposantesActivitesStatiques):
    def __init__(self,constellation,foldername):
        super().__init__(constellation,[],modeleMoyen) # le modele importe peu. On va tout écraser
        mapping = self.mapCCAToCores(foldername)
        self.initSharedStructures(constellation)
        for filename in os.listdir(foldername):
            split = filename.split(".comp")[0].split("_")
            s,num_cca = int(split[1]),int(split[2])
            activites = []
            with open(foldername+'/'+filename,'r') as file:
                nObs = int(file.readline())
                for i in range(nObs):
                    idObs = int(file.readline().split(" ")[0])
                    if idObs>=0: # -1 = obs fictive
                        try: # verifier que l'obs ne fait pas partie des requetes supprimees (systematiques)
                            constellation.getRequeteActivite(idObs)
                            activites.append(constellation.getActivite(idObs))
                        except Exception as e:
                            pass
                        if(idObs in self.activite_cca):
                            raise ValueError("Erreur de numérotation des activités : ",idObs,self.activite_cca[idObs],(s,num_cca))
                        self.activite_cca[idObs] = (s,num_cca)
                self.tailles[(s,num_cca)] = nObs-1
                if s not in self.sat_cca:
                    self.sat_cca[s] = []
                self.sat_cca[s].append(num_cca)
                composante = self.Composante((s,num_cca),0)
                composante.ajouterListeActivites(activites)
                self.composantes[(s,num_cca)] = composante
    
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
        return np.ndarray(buffer=buf, dtype='d', shape=(size,)) 
        
    def initSharedStructures(self,constellation):
        size = constellation.getNombreActivites()
        self.NtimeSlotsMapping = self.createSharedArrayOfDouble(size)
        self.dureeTimeSlotsMapping = self.createSharedArrayOfDouble(size)
        
    
    def mapCCAToCores(self,foldername):
        files = sorted(os.listdir(foldername),key = lambda file : os.stat(foldername+'/'+file).st_size,reverse=True)
        mapping = []
        for i in range(MPI.COMM_WORLD.Get_rank(),len(files),MPI.COMM_WORLD.Get_size()):
            mapping.append(files[i])
        return mapping
        
class GrapheDependancesDynamique:
    def __init__(self,constellation,activites,modeleDeTransition):
        self.tau_max = modeleDeTransition.getMajorantDureeTransition()
        self.DEBUT = 0
        self.FIN = 1
        events = collections.OrderedDict(sorted({s : [] for s in constellation.getSatellites()}.items()))
        self.mode_listing = {s : {} for s in constellation.getSatellites()}
        #self.activites_listing = {s : {} for s in constellation.getSatellites()}
        
        for s in activites:
            for a in activites[s]:
                debut = constellation.getSatellite(s).getActivite(a).getDebut()
                fin = constellation.getSatellite(s).getActivite(a).getFin() + self.tau_max
                
                points = []
                points.append(EvenementActivite(s,debut,a,self.DEBUT))
                points.append(EvenementActivite(s,fin,a,self.FIN))
                
                for pt in points:
                    if pt not in events[s]: # sinon on a plusieurs fois les vidages
                        bisect.insort_left(events[s],pt)
                
        # stoquage => ajout de nouveau point + rapide
        self.pointsInterets = collections.OrderedDict(sorted({s : collections.OrderedDict({}) for s in constellation.getSatellites()}.items()))
        self.grapheActivites = GrapheActivitesDynamique()
        for s in events:
            activites_actives = []
            for event in events[s]:
                t,a,e = event.getDate(),event.getActivite(),event.getEvent()
                if e==self.DEBUT:
                    self._activerActivite(constellation,activites_actives,a)
                else:
                    self._desactiverActivite(constellation,s,activites_actives,a)
                self.pointsInterets[s][t] = activites_actives.copy()
        self.grapheActivites.calculerComposantes()
        printColor(str(self.grapheActivites),c='b')
        
        
        self.cca_sat = {}
        self.sat_cca = {}
        for s in activites:
            for a in activites[s]:
                cca = self.getActiviteCCA(a)
                self.cca_sat[cca] = s
                if s not in self.sat_cca:
                    self.sat_cca[s] = []
                self.sat_cca[s].append(cca)
    
    def getCCA(self,cca):
        return self.grapheActivites.composantes[cca]
    
    def getCCASatellite(self,s):
        return self.sat_cca[s]
        
    def getSatelliteCCA(self,cca):
        return self.cca_sat[cca]

    def getTailleComposante(self,cca):
        return self.grapheActivites.getTailleComposante(cca)

    def getNombreComposantes(self):
        return self.grapheActivites.getNombreComposantes()
    
    def getNombreElements(self):
        return self.grapheActivites.nombreNoeuds()
    
    def getComposantes(self):
        return list(self.grapheActivites.composantes.keys())
    
    def getCCAs(self):
        return self.grapheActivites.composantes
    
    def getActiviteCCA(self,a):
        return self.grapheActivites.getNoeud(a).getComposanteConnexe()
    
    """
        =================================================
                    fonctions internes
        =================================================
    """
    def compterModes(self,s,a):
        return len(self.mode_listing[s][a])
    
    def verifierActivitePresente(self,constellation,s,a):
        assert(a in self.grapheActivites.getNoeuds())
        t_debut = constellation.getSatellite(s).getActivite(a).getDebut()
        t_fin = constellation.getSatellite(s).getActivite(a).getFin() + self.tau_max
        tps_actifs = [ t for t in list(self.pointsInterets[s].keys()) if t>=t_debut and t <t_fin]
        tps_non_actifs = [ t for t in list(self.pointsInterets[s].keys()) if t<t_debut and t >=t_fin]
        for t in tps_actifs:
            if not a in self.pointsInterets[s][t][0]:
                print(a,"not in ",[self.pointsInterets[s][t][0] for t in tps_actifs])
                assert(False)
        for t in tps_non_actifs:
            assert(a not in self.pointsInterets[s][t][0])
    
    def verifierActiviteNonPresente(self,constellation,s,a):
        assert(a not in self.grapheActivites.getNoeuds())
        for t in self.pointsInterets[s]:
            assert(a not in self.pointsInterets[s][t][0])
        
    def verifications(self,constellation):
        # présence des activités
        for (r,m) in self.grapheModes.getNoeuds():
            for (s,o,d) in constellation.getRequete(r).getMode(m).getCouples():
                assert(o in self.grapheActivites.getNoeuds())
                assert(d in self.grapheActivites.getNoeuds())
        
        # présence des modes
        for a in self.grapheActivites.getNoeuds():
            findMode = False
            for (r,m) in self.grapheModes.getNoeuds():
                for (s,o,d) in constellation.getRequete(r).getMode(m).getCouples():
                    if a==o or a==d:
                        findMode = True
                        break
                if findMode:
                    break
            assert(findMode)
                
        # verification composantes modes
        for (r,m) in self.grapheModes.getNoeuds():
            c = self.grapheModes.getNoeud((r,m)).getComposanteConnexe()
            assert((r,m) in self.grapheModes.getComposanteConnexe(c).getNoeuds())
        
        for c in self.grapheModes.getComposantesConnexes():
            for (r,m) in self.grapheModes.getComposanteConnexe(c).getNoeuds():
                assert(self.grapheModes.getNoeud((r,m)).getComposanteConnexe()==c)
        
        # vérification composantes acitvités
        for a in self.grapheActivites.getNoeuds():
            c = self.grapheActivites.getNoeud(a).getComposanteConnexe()
            assert(a in self.grapheActivites.getComposanteConnexe(c).getNoeuds())
            
        for c in self.grapheActivites.getComposantesConnexes():
            for a in self.grapheActivites.getComposanteConnexe(c).getNoeuds():
                assert(self.grapheActivites.getNoeud(a).getComposanteConnexe()==c)
                
        # satellites différents => activités non liés
        for n in self.grapheActivites.getNoeuds():
            for voisin in self.getNoeuds()[n].getVoisins():
                if (constellation.getSatelliteActivite(a1)!=constellation.getSatelliteActivite(a2)):
                    print(constellation.getSatelliteActivite(a1),constellation.getSatelliteActivite(a2))
                    assert(False)
        
    def _addListingModesActivites(self,constellation,s,o,d,r,m):
        if o not in self.mode_listing[s]:
            self.mode_listing[s][o] = []
        assert ((r,m) not in self.mode_listing[s][o]) # dans 1 mode 1 obs n'apparait qu'une fois
        self.mode_listing[s][o].append((r,m))
        if d not in self.mode_listing[s]:
            self.mode_listing[s][d] = []
        assert ((r,m,o) not in self.mode_listing[s][d]) # les vidages apparaissent plusieurs fois dans un mode
        self.mode_listing[s][d].append((r,m,o))
    
    def _removeListingModesActivites(self,s,o,d,r,m):
        self.mode_listing[s][o].remove((r,m))
        self.mode_listing[s][d].remove((r,m,o))

    def _getActifs(self,s,t):
        if t>-np.Inf:
            return (self.pointsInterets[s][t][0].copy(),self.pointsInterets[s][t][1].copy()) 
        else:
            return ([],[])
    
    # un element de la liste qui est dans l'intervalle a,b
    # si liste[i] n'y est pas alors il n'y en a pas
    def _getElementInterval(self,liste,a,b):
        if len(liste)==0:
            return 0
        l = 0
        r = len(liste)-1
        while l<=r:
            i = floor((l+r)/2)
            if liste[i]<=b:
                l = i + 1
            if liste[i]>=a:
                r = i - 1
        return i
    
    # indice des tps pour lequels on doit modifier la liste d'actifs
    def _leftRightActifs(self,temps_actifs,i,a,b):
        if len(temps_actifs)==0:
            return 1,-1
        l = i
        r = i
        while temps_actifs[l]>a and l>0:
            l -= 1
        l += 1
        while temps_actifs[r]<b and l<len(temps_actifs)-1:
            r += 1
        r -= 1
        return l,r
        
    def _activerActivite(self,constellation,activites_actives,a):
        if a not in activites_actives:
            if a not in self.grapheActivites.getNoeuds():
                self.grapheActivites.ajouterActivite(a)
            for a2 in activites_actives:
                self.grapheActivites.lierActivites(a,a2,MAJComposantes=False)
            activites_actives.append(a)

    def _desactiverActivite(self,constellation,s,activites_actives,a):
        # sinon plusieurs desactivation de l'activité qui peut apparaitre dans plusieurs modes
        if a in activites_actives:
            activites_actives.remove(a)     
    
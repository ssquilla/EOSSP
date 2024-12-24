import collections
import bisect
from mpi4py import MPI
from EOSSP.model.graph import *

from EOSSP.Utils import *
from EOSSP.Utils.config import *
config = Config()

from EOSSP.model.transitionModel import MeanModel

class ActivityDateEvent: # event = (0:debut|1:fin)
    def __init__(self,sat,date,activity,event):
        self.satellite = sat
        self.activity = activity
        self.event = event
        self.date = date
        
    def getActivity(self):
        return self.activity
    
    def getEvent(self):
        return self.event
    
    def getDate(self):
        return self.date
    
    def getSatellite(self):
        return self.satellite
    
    def __lt__(self, other):
        return (self.satellite,self.date,self.event,self.activity) < (self.satellite,other.date,other.event,other.activity)
    
    def __gt__(self, other):
        return (self.satellite,self.date,self.event,self.activity) > (self.satellite,other.date,other.event,other.activity)
    
    def __eq__(self, other):
        return (self.satellite,self.date,self.event,self.activity) == (self.satellite,other.date,other.event,other.activity)
    
    def __le__(self,other):
        return (self.satellite,self.date,self.event,self.activity) <= (self.satellite,other.date,other.event,other.activity)
    
    def __ge__(self,other):
        return (self.satellite,self.date,self.event,self.activity) >= (self.satellite,other.date,other.event,other.activity)
        
class PlanningTask:
    def __init__(self,cca,listOfModes):
        self.cca = cca
        self.listOfModes = listOfModes
    
    def __str__(self):
        return 'PlanningTask(cca:'+str(self.cca)+',modes:'+str(self.listOfModes)+')'
    
    def getCCA(self):
        return self.cca
    
    def getListOfModes(self):
        return self.listOfModes

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
    def mergeQuickUnion(self,other):
        self.lst = safeMergeDict(self.lst,other.lst)
        
    def overwriteQuickUnion(self,other):
        self.lst = other.lst
    
class StaticCCAsGroup:
    class Component:
        def __init__(self,ccaIndex,tau_max):
            self.tau_max = tau_max
            self.id = ccaIndex
            self.activities = []
            self.startDate = np.Inf
            self.endDate = -np.Inf
            self.satellite = ccaIndex[0]
        def isEmpty(self):
            return len(self.activities)==0
        def getStartDate(self):
            return self.startDate
        def getEndDate(self):
            return self.endDate
        def __str__(self):
            return "Connected component\n| size : "+str(len(self.activities))+"\n| start date : "+str(self.startDate)+"\n| end date : " + str(self.endDate)
        def __contains__(self,elmt):
            if self.satellite!=elmt.getSat():
                return False
            res = intersect((elmt.getStartDate()-self.tau_max,elmt.getEndDate()+self.tau_max),(self.startDate,self.endDate))
            return res
        def addActivity(self,activity):
            self.activities.append(activity.getId())
            self.startDate = min(self.startDate,activity.getStartDate())
            self.endDate = max(self.endDate,activity.getEndDate())
        def addActivitiesList(self,activities):
            for a in activities:
                self.addActivity(a)
        def mergeComponents(self,other):
            if self.getStartDate() < other.getStartDate():
                assert(self.getEndDate() + self.tau_max>=other.getStartDate())
                self.endDate = other.getEndDate()
            else:
                assert(other.getEndDate() + self.tau_max<=self.getStartDate())
                self.startDate = other.getStartDate()
            self.activities += other.activities
        def getElements(self):
            return self.activities
        def getId(self):
            return self.id
            
    def __init__(self,constellation,activities,transitionModel,parallel=True):
        self.tau_max = transitionModel.getTransititionDurationUpperBound()
        self.ccaActivity = {}
        self.cca_sat = {}
        if parallel:
            # on parallelise le calcul des composantes
            allocatedSatellites = self.allocateSatellites(activities)
            self.computeComponents({s : activities.get(s,[]) for s in allocatedSatellites},constellation)
            self.mergeComponents()
        else:
            self.computeComponents(activities,constellation)
   
    def overwrite(self,other):
        self.satelliteCCA = other.satelliteCCA
        self.sizes = other.sizes
        self.ccaActivity = other.ccaActivity
        self.components = other.components
                
    def mergeComponents(self):
        MPI.COMM_WORLD.Barrier()
        if MPI.COMM_WORLD.Get_rank()>0:
            data = {MPI.COMM_WORLD.Get_rank():self}
            MPI.COMM_WORLD.send(data,dest=0)
            data = MPI.COMM_WORLD.recv(source=0)
        else:
            for i in range(1,MPI.COMM_WORLD.Get_size()):
                data = MPI.COMM_WORLD.recv()
                groupsToMerge = list(data.values())[0]
                self.mergeGroupOfComponents(groupsToMerge)
            data = {'':self}
            for i in range(1,MPI.COMM_WORLD.Get_size()):
                MPI.COMM_WORLD.send(data,dest=i)
        MPI.COMM_WORLD.Barrier()
        if MPI.COMM_WORLD.Get_rank()>0:
            wholeGroup = list(data.values())[0]
            self.overwrite(wholeGroup)
        MPI.COMM_WORLD.Barrier()
    def addRequest(self,constellation,r):
        for (s,a) in constellation.getRequete(r).getCouples():
            activity = constellation.getActivity(a,s=s)
            self.addNewActivity(activity)
    def addNewActivity(self,activity):
        components_trouvees = []
        for c in self.components:
            if activity in self.components[c]:
                if len(components_trouvees)==0:
                    self.addActivityToComponent(activity,self.components[c])
                components_trouvees.append(c)
                
        if len(components_trouvees)==0:
            s = activity.getSat()
            indice = max(self.getCCASatellite(s))+1
            newComponent = self.Component((s,indice),self.tau_max)
            self.sizes[(s,indice)] = 1
            self.addActivityToComponent(activity,newComponent)
            self.components[(s,indice)] = newComponent
            self.satelliteCCA[s].append(indice)
        else:
            self.mergeListOfComponents(components_trouvees)
                    
    def mergeListOfComponents(self,liste_components):
        assert(len(liste_components)>0)
        if len(liste_components)==1:
            return
        else:
            head = liste_components.pop()
            s_head = head[0]
            for other in liste_components:
                self.components[head].mergeComponents(self.components[other])
                s,cca = other
                assert(s==s_head)
                self.satelliteCCA[s].remove(cca)
                self.sizes[head]+= self.sizes[(s,cca)]
                del self.sizes[(s,cca)]
            for a in self.ccaActivity:
                if self.ccaActivity[a] in liste_components:
                    self.ccaActivity[a] = head
            
    def addActivityToComponent(self,activity,component):
        assert(activity.getId() not in self.ccaActivity)
        self.ccaActivity[activity.getId()] = component.getId()
        component.addActivity(activity)
        if component.getId() not in self.sizes:
            self.sizes[component.getId()] = 0
        self.sizes[component.getId()] += 1
        
    def computeComponents(self,activities,constellation):
        self.components = {}
        self.sizes = {}
        for s in activities:
            ccaCounter = 0
            currentComponent = self.Component((s,ccaCounter),self.tau_max)
            ccaCounter += 1
            for i,a in enumerate(sorted(activities[s],key = lambda a : constellation.getActivity(a).getStartDate())):
                activity = constellation.getActivity(a)
                if not (currentComponent.isEmpty() or activity in currentComponent):
                    self.components[currentComponent.getId()] = currentComponent
                    currentComponent = self.Component((s,ccaCounter),self.tau_max)
                    ccaCounter += 1
                self.addActivityToComponent(activity,currentComponent)
            if currentComponent.getId() not in self.components and not currentComponent.isEmpty():
                self.components[currentComponent.getId()] = currentComponent

        self.satelliteCCA = {}
        for a in self.ccaActivity:
            s,cca = self.ccaActivity[a]
            if s not in self.satelliteCCA:
                self.satelliteCCA[s] = []
            self.satelliteCCA[s].append(cca)
            if (s,cca) not in self.cca_sat:
                self.cca_sat[(s,cca)] = s

    def deleteActivity(self,a):
        cca = self.ccaActivity[a]
        self.sizes[cca] -= 1
        del self.ccaActivity[a]
        self.components[cca].activities.remove(a)
     
    def allocateSatellites(self,activities):
        size = MPI.COMM_WORLD.Get_size()
        rank = MPI.COMM_WORLD.Get_rank()
        sizes_triees = sorted([s for s in activities],key=lambda s : len(activities[s]),reverse=True)
        allocatedSatellites = [i for i in sizes_triees if (i%size)==rank]
        return allocatedSatellites
        
    def mergeGroupOfComponents(self,other):
        self.ccaActivity = safeMergeDict(self.ccaActivity,other.ccaActivity)
        self.satelliteCCA = safeMergeDict(self.satelliteCCA,other.satelliteCCA)
        self.sizes = safeMergeDict(self.sizes,other.sizes)
        self.components = safeMergeDict(self.components,other.components)
            
    def mergeGroupList(self,others):
        for o in others:
            self.mergeGroupOfComponents(o)
        
    def getCCASatellite(self,s):
        return self.satelliteCCA[s]

    def getSatelliteCCA(self,cca):
        return self.cca_sat[cca]
    
    def getActivityCCA(self,a):
        return self.ccaActivity[a]

    def getComponentSize(self,cca):
        return self.sizes[cca]

    def getVertices(self):
        return self.ccaActivity
    def getActivitiesOfComponent(self,cca):
        return self.components[cca].getElements()
    def getConnectedComponent(self,cca):
        return self.components[cca]
    
    def getNumberOfComponents(self):
        return len(list(self.sizes.keys()))
    
    def getNumberOfElements(self):
        return len(list(self.ccaActivity.keys()))

    def getComponents(self):
        return list(self.components.keys())
    
    def __str__(self):
        min_taille = min([self.getComponentSize(cca) for cca in self.components])
        max_taille = max([self.getComponentSize(cca) for cca in self.components])
        mean_taille = np.mean([self.getComponentSize(cca) for cca in self.components])
        std_taille = np.std([self.getComponentSize(cca) for cca in self.components])
        mess = "Components of activities: processus "+str(MPI.COMM_WORLD.Get_rank())+"\n"
        mess += "| number of components: "+str(len(list(self.sizes.keys())))+'\n'
        mess += "| number of activities: "+str(len(list(self.ccaActivity.keys())))+'\n'
        mess += "| size of shortest component: "+str(min_taille)+'\n'
        mess += "| size of largest component: " + str(max_taille)+'\n'
        mess += "| size of mean component: " + str(mean_taille)+'\n'
        mess += "| std of component sizes: " + str(std_taille)+'\n'
        return mess

class PrecomputedStaticComponents(StaticCCAsGroup):
    def __init__(self,constellation,foldername):
        super().__init__(constellation,[],meanModel) # le modele importe peu. On va tout écraser
        for filename in os.listdir(foldername):
            split = filename.split(".comp")[0].split("_")
            s,num_cca = int(split[1]),int(split[2])
            activities = []
            with open(foldername+'/'+filename,'r') as file:
                nObs = int(file.readline())
                count = 0
                for i in range(nObs):
                    idObs = int(file.readline().split(" ")[0])
                    if idObs>=0: # -1 = obs fictive
                        if(idObs in self.ccaActivity):
                            raise ValueError("Wrong activities numbering: ",idObs,self.ccaActivity[idObs],(s,num_cca))
                        try: # verifier que l'obs ne fait pas partie des requetes supprimees (systematiques)
                            constellation.getRequestAtivity(idObs)
                            count += 1
                            activities.append(constellation.getActivity(idObs))
                            self.ccaActivity[idObs] = (s,num_cca)
                        except Exception as e: # obs de requetes ignorees (systematiques par exemple)
                            pass
                self.sizes[(s,num_cca)] = count
                if s not in self.satelliteCCA:
                    self.satelliteCCA[s] = []
                self.satelliteCCA[s].append(num_cca)
                component = self.Component((s,num_cca),0)
                component.addActivitiesList(activities)
                self.components[(s,num_cca)] = component
            
class DynamicDependencyGraphOfActivities:
    def __init__(self,constellation,activities,transitionModel):
        self.tau_max = transitionModel.getTransititionDurationUpperBound()
        self.START = 0
        self.END = 1
        events = collections.OrderedDict(sorted({s : [] for s in constellation.getSatellites()}.items()))
        self.listOfModes = {s : {} for s in constellation.getSatellites()}
        
        for s in activities:
            for a in activities[s]:
                startDate = constellation.getSatellite(s).getActivity(a).getStartDate()
                endDate = constellation.getSatellite(s).getActivity(a).getEndDate() + self.tau_max
                
                points = []
                points.append(ActivityDateEvent(s,startDate,a,self.START))
                points.append(ActivityDateEvent(s,endDate,a,self.END))
                
                for pt in points:
                    if pt not in events[s]: # sinon on a plusieurs fois les vidages
                        bisect.insort_left(events[s],pt)
                
        # stoquage => ajout de nouveau point + rapide
        self.pointsOfInterest = collections.OrderedDict(sorted({s : collections.OrderedDict({}) for s in constellation.getSatellites()}.items()))
        self.graphOfActivities = DynamicGraphOfActivities()
        for s in events:
            activeActivities = []
            for event in events[s]:
                t,a,e = event.getDate(),event.getActivity(),event.getEvent()
                if e==self.START:
                    self._activateActivity(constellation,activeActivities,a)
                else:
                    self._deactivateActivity(constellation,s,activeActivities,a)
                self.pointsOfInterest[s][t] = activeActivities.copy()
        self.graphOfActivities.computeComponents()
        printColor(str(self.graphOfActivities),c='b')
        
        
        self.cca_sat = {}
        self.satelliteCCA = {}
        for s in activities:
            for a in activities[s]:
                cca = self.getActivityCCA(a)
                self.cca_sat[cca] = s
                if s not in self.satelliteCCA:
                    self.satelliteCCA[s] = []
                self.satelliteCCA[s].append(cca)
    
    def getConnectedComponent(self, cca):
        self.getCCA(cca)
        
    def getCCA(self,cca):
        return self.graphOfActivities.components[cca]
    
    def getCCASatellite(self,s):
        return self.satelliteCCA[s]
        
    def getSatelliteCCA(self,cca):
        return self.cca_sat[cca]

    def getComponentSize(self,cca):
        return self.graphOfActivities.getComponentSize(cca)

    def getNumberOfComponents(self):
        return self.graphOfActivities.getNumberOfComponents()
    
    def getNumberOfElements(self):
        return self.graphOfActivities.getNumberOfNodes()
    
    def getComponents(self):
        return list(self.graphOfActivities.components.keys())
    
    def getCCAs(self):
        return self.graphOfActivities.components
    
    def getActivityCCA(self,a):
        return self.graphOfActivities.getVertex(a).getConnectedComponent()
    
    def getConnectedComponents(self):
        return self.getComponents()
    
    """
        =================================================
                    fonctions internes
        =================================================
    """
    def countModes(self,s,a):
        return len(self.listOfModes[s][a])
    
    def checkActivityPresence(self,constellation,s,a):
        assert(a in self.graphOfActivities.getVertices())
        startDate = constellation.getSatellite(s).getActivity(a).getStartDate()
        endDate = constellation.getSatellite(s).getActivity(a).getEndDate() + self.tau_max
        activeDates = [ t for t in list(self.pointsOfInterest[s].keys()) if t>=startDate and t <endDate]
        unactiveDates = [ t for t in list(self.pointsOfInterest[s].keys()) if t<startDate and t >=endDate]
        for t in activeDates:
            if not a in self.pointsOfInterest[s][t][0]:
                print(a,"not in ",[self.pointsOfInterest[s][t][0] for t in activeDates])
                assert(False)
        for t in unactiveDates:
            assert(a not in self.pointsOfInterest[s][t][0])
    
    def checkActivityNonPresence(self,constellation,s,a):
        assert(a not in self.graphOfActivities.getVertices())
        for t in self.pointsOfInterest[s]:
            assert(a not in self.pointsOfInterest[s][t][0])
        
    def verifications(self,constellation):
        # présence des activités
        for (r,m) in self.graphOfModes.getVertices():
            for (s,o,d) in constellation.getRequete(r).getMode(m).getCouples():
                assert(o in self.graphOfActivities.getVertices())
                assert(d in self.graphOfActivities.getVertices())
        
        # présence des modes
        for a in self.graphOfActivities.getVertices():
            findMode = False
            for (r,m) in self.graphOfModes.getVertices():
                for (s,o,d) in constellation.getRequete(r).getMode(m).getCouples():
                    if a==o or a==d:
                        findMode = True
                        break
                if findMode:
                    break
            assert(findMode)
                
        # verification composantes modes
        for (r,m) in self.graphOfModes.getVertices():
            c = self.graphOfModes.getVertex((r,m)).getConnectedComponent()
            assert((r,m) in self.graphOfModes.getConnectedComponent(c).getVertices())
        
        for c in self.graphOfModes.getConnectedComponents():
            for (r,m) in self.graphOfModes.getConnectedComponent(c).getVertices():
                assert(self.graphOfModes.getVertex((r,m)).getConnectedComponent()==c)
        
        # vérification composantes acitvités
        for a in self.graphOfActivities.getVertices():
            c = self.graphOfActivities.getVertex(a).getConnectedComponent()
            assert(a in self.graphOfActivities.getConnectedComponent(c).getVertices())
            
        for c in self.graphOfActivities.getConnectedComponents():
            for a in self.graphOfActivities.getConnectedComponent(c).getVertices():
                assert(self.graphOfActivities.getVertex(a).getConnectedComponent()==c)
                
        # satellites différents => activités non liés
        for n in self.graphOfActivities.getVertices():
            for voisin in self.getVertices()[n].getVoisins():
                if (constellation.getSatelliteActivite(a1)!=constellation.getSatelliteActivite(a2)):
                    print(constellation.getSatelliteActivite(a1),constellation.getSatelliteActivite(a2))
                    assert(False)
        
    def _addListingModesActivites(self,constellation,s,o,d,r,m):
        if o not in self.listOfModes[s]:
            self.listOfModes[s][o] = []
        assert ((r,m) not in self.listOfModes[s][o]) # dans 1 mode 1 obs n'apparait qu'une fois
        self.listOfModes[s][o].append((r,m))
        if d not in self.listOfModes[s]:
            self.listOfModes[s][d] = []
        assert ((r,m,o) not in self.listOfModes[s][d]) # les vidages apparaissent plusieurs fois dans un mode
        self.listOfModes[s][d].append((r,m,o))
    
    def _removeListingModesActivites(self,s,o,d,r,m):
        self.listOfModes[s][o].remove((r,m))
        self.listOfModes[s][d].remove((r,m,o))

    def _getActives(self,s,t):
        if t>-np.Inf:
            return (self.pointsOfInterest[s][t][0].copy(),self.pointsOfInterest[s][t][1].copy()) 
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
    def _leftRightActives(self,temps_actifs,i,a,b):
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
        
    def _activateActivity(self,constellation,activeActivities,a):
        if a not in activeActivities:
            if a not in self.graphOfActivities.getVertices():
                self.graphOfActivities.addActivity(a)
            for a2 in activeActivities:
                self.graphOfActivities.linkActivities(a,a2,MAJComponents=False)
            activeActivities.append(a)

    def _deactivateActivity(self,constellation,s,activeActivities,a):
        # sinon plusieurs desactivation de l'activité qui peut apparaitre dans plusieurs modes
        if a in activeActivities:
            activeActivities.remove(a)     
    
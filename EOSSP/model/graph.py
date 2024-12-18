import numpy as np
from EOSSP.Utils.Utils import *

"""
Graph: util tools to create and manipulated graphs for modes and activities.
    One important tool is connected components computation.
"""
class Graph:
    class ConnectedComponent:
        def __init__(self,verticesList):
            self.verticesList = verticesList
            
        def addVertex(self,vertex):
            self.verticesList.append(vertex)
            
        def removeVertex(self,vertex):
            self.verticesList.remove(vertex)
        
        def addVertices(self,verticesList):
            for vertex in verticesList:
                self.addVertex(vertex)
        
        def getVertices(self):
            return self.verticesList
        
        def getNumberOfVertices(self):
            return len(self.verticesList)
        
        def __str__(self):
            return str(self.verticesList)
        
    class NonOrientedVertex:
        def __init__(self,label):
            self.label = label
            self.connectedComponent = label
            self.neighbors = []
            
        def setConnectedComponent(self,comp):
            self.connectedComponent = comp
        
        def getConnectedComponent(self):
            return self.connectedComponent
        
        def getLabel(self):
            return self.label
        
        def getNeighbors(self):
            return self.neighbors
        
        def addNeighbor(self,s):
            self.neighbors.append(s)
            
        def removeNeighbor(self,s):
            self.neighbors.remove(s)
        
        def __str__(self):
            return str(self.label)

    class OrientedNode:
        def __init__(self,label,data=None):
            self.label = label
            self.connectedComponent = label
            self.successors = []
            self.predecessors = []
            self.data = data
           
        def getData(self):
            return self.data
            
        def setConnectedComponent(self,comp):
            self.connectedComponent = comp
        
        def getConnectedComponent(self):
            return self.composanteConnexe
        
        def getLabel(self):
            return self.label
        
        def getPredecessors(self):
            return self.predecessors
        
        def getSuccessors(self):
            return self.successors
        
        def addSuccessor(self,s):
            self.successors.append(s)
            
        def removeSuccessor(self,s):
            self.successors.remove(s)
            
        def addPredecessor(self,s):
            self.predecessors.append(s)
            
        def removePredecessor(self,s):
            self.predecessors.remove(s)
        
        def __str__(self):
            return str(self.label)    
            
    # classe Graph
    def __init__(self,name,oriented):
        self.oriented = oriented
        self.name = name
        self.vertices = {}
        self.components = {}

    def __str__(self):
        res = str(len(list(self.vertices.keys()))) + "vertices " + str(len(list(self.components.keys()))) + " components"
        return res
    
    def getConnectedComponents(self):
        return self.components
    
    def getConnectedComponent(self,c):
        return self.components[c]
    
    def getNumberOfComponents(self):
        return len(list(self.components.keys()))
    
    def checkVertexPresence(self,label):
        return label in self.vertices
    
    def addVertex(self,label,data=None):
        assert(label not in self.vertices)
        if self.oriented:
            self.vertices[label] = self.OrientedNode(label,data)
        else:
            self.vertices[label] = self.NonOrientedVertex(label,data)
        
        self.components[label] = self.ConnectedComponent([label])
        
    def addEdge(self,label1,label2,label=None,MAJComponents=True):
        assert(label1 != label2)
        left = label1
        right = label2
        if self.oriented:
            self.vertices[left].addSuccessor(right)
            self.vertices[right].addPredecessor(left)
        else:
            self.vertices[left].addNeighbor(right)
            self.vertices[right].addNeighbor(left)
        
        if MAJComponents:
            c1 = self.vertices[label1].getConnectedComponent()
            c2 = self.vertices[label2].getConnectedComponent()
            if c1!=c2:
                self._mergeComponents(c1,c2)
                return (c1,c2)
            else:
                return None
        else:
            return None
        
    def getName(self):
        return self.name
    
    def export(self):
        with open(self.name+".dot","w") as file:
            self._writeFraph(file)
    
    def setComponent(self,verticesList,c):
        self.components[c] = self.ConnectedComponent(verticesList)
        for vertex in verticesList:
            self.vertices[vertex].setConnectedComponent(c)
        
    def _mergeComponents(self,c1,c2):
        vertices = self.components[c2].getVertices()
        for vertex in vertices:
            self.vertices[vertex].setConnectedComponent(c1)
        self.components[c1].addVertices(vertices)
        del self.components[c2]
    
    def getVertices(self):
        return self.vertices
    
    def getVertex(self,vertex):
        return self.vertices[vertex]
    
    # retourne : les composantes a recalculer, les composantes a renommer
    def removeListOfVertices(self,listOfLabels):
        # creer la liste des vertexs a mettre a jour
        connectedVertices = []
        removalCountByComponent = {}
        for label in listOfLabels:
            cca = self.vertices[label].getConnectedComponent()
            if cca not in removalCountByComponent:
                removalCountByComponent[cca] = 0
            removalCountByComponent[cca] += 1
            if self.oriented:
                for pred in self.vertices[label].getPredecessors():
                    self.vertices[pred].removeSuccessor(label)
                    if pred not in connectedVertices and pred not in listOfLabels:
                        connectedVertices.append(pred)
                for succ in self.vertices[label].getSuccessors() and succ not in listOfLabels:
                    self.vertices[succ].removePredecessor(label)
                    if succ not in connectedVertices:
                        connectedVertices.append(succ)
            else:       
                for v in self.vertices[label].getNeighbors():
                    self.vertices[v].removeNeighbor(label)
                    if v not in connectedVertices and v not in listOfLabels:
                        connectedVertices.append(v)
        # effacer les anciennes composantes connexes
        oldComponents = {}
        for label in listOfLabels:
            c = self.vertices[label].getConnectedComponent()
            if c not in oldComponents:
                oldComponents[c] = self.components[c].getNumberOfVertices()
                del self.components[c]
            del self.vertices[label]
        verticesDone = []
        componentsToRecompute = [] # les composantes a recalculer (et a renommer)
        componentsToRename = {} # les composantes a renommer (et pas a recalculer)
        # creer les nouvelles composantes
        while len(connectedVertices)>0:
            currentVertex = connectedVertices.pop()
            comp = currentVertex
            oldComponent = self.vertices[currentVertex].getConnectedComponent()
            oldComponentSize = oldComponents[oldComponent]
            component = self.getComponentVertices(currentVertex)
            newComponentSize = len(component)
            if newComponentSize<oldComponentSize-removalCountByComponent[cca]:
                if oldComponent not in componentsToRecompute:
                    componentsToRecompute.append(oldComponent)
            elif oldComponent!=comp:
                componentsToRename[oldComponent] = comp
            # ne pas retraiter les vertexs connectes qui sont dans la meme composante
            verticesToDelete = []
            for vertex in connectedVertices:
                if vertex in component:
                    verticesToDelete.append(vertex)
            for vertex in verticesToDelete:
                connectedVertices.remove(vertex)
            self.setComponent(component,comp)
        return componentsToRecompute,componentsToRename
    
    # retirer le vertex + les aretes + recalcule les components
    def removeVertex(self,label):
        connectedVertices = []
        if self.oriented:
            for pred in self.vertices[label].getPredecessors():
                self.vertices[pred].removeSuccessor(label)
                if pred not in connectedVertices:
                    connectedVertices.append(pred)
            for succ in self.vertices[label].getsuccessors():
                self.vertices[succ].removePredecessor(label)
                if succ not in connectedVertices:
                    connectedVertices.append(succ)
        else:       
            for v in self.vertices[label].getNeighbors():
                self.vertices[v].removeNeighbor(label)
                if v not in connectedVertices:
                    connectedVertices.append(v)
        #assert(label not in connectedVertices)
        oldComponent = self.vertices[label].getConnectedComponent()
        del self.vertices[label]
        del self.components[oldComponent]
        # il faut recalculer les components connexes ici
        verticesDone = []
        i = 0
        while len(connectedVertices)>0:
            currentVertex = connectedVertices.pop()
            comp = currentVertex
            assert(currentVertex!=label)
            component = self.getComponentVertices(currentVertex)
            # ne pas retraiter les vertexs connectes qui sont dans la meme composante
            verticesToDelete = []
            for vertex in connectedVertices:
                if vertex in component:
                    verticesToDelete.append(vertex)
            for vertex in verticesToDelete:
                connectedVertices.remove(vertex)
            self.setComponent(component,comp)
            i += 1
        changesDetected = (i>1 or (i==1 and currentVertex!=oldComponent))
        return oldComponent,changesDetected,i
    
    def getComponentVertices(self,startingVertex):
        visited = {}
        V = list(self.vertices.keys())
        for i in V:
            visited[i] = False
        return self.computeConnectedComponentRecursively([], startingVertex, visited)
   
    def computeConnectedComponentRecursively(self, temp, v, visited):
        # Mark the current vertex as visited
        visited[v] = True
        # Store the vertex to list
        temp.append(v)
        # Repeat for all vertices adjacent
        # to this vertex v
        if self.oriented:
            adj = self.getVertex(v).getSuccessors() + self.getVertex(v).getPredecessors()
        else:
            adj = self.getVertex(v).getNeighbors()
        for i in adj:
            if visited[i] == False:
                # Update the list
                temp = self.computeConnectedComponentRecursively(temp, i, visited)
        return temp
        
    def computeComponents(self):
        self.components = {}
        visited = {}
        V = list(self.vertices.keys())
        for v in V:
            visited[v] = False
        for v in V:
            if visited[v] == False:
                temp = []
                self.setComponent(self.computeConnectedComponentRecursively(temp, v, visited),v)
                
    def getNumberOfVertices(self):
        return len(list(self.vertices.keys()))

class DynamicGraphOfActivities(Graph):
    def __init__(self):
        super().__init__("Graph of activities",oriented=False)
        
    def addActivity(self,activity):
        self.addVertex(activity)
        
    def removeActivity(self,activity,constellation,solution):
        oldComponent,changesDetected,numberOfComponent = self.removeVertex(activity)
        if changesDetected:
            if numberOfComponent==1:
                solution.renameComponent(constellation,oldComponent)
            else:
                solution.notifyCCASplit(constellation,oldComponent)
    
    def removeListOfActivities(self,activities,constellation,solution):
        recomputed,renamed = self.removeListOfVertices(activities)
        for oldCCA in renamemages:
            nouveau_name = renamed[oldCCA]
            solution.renameComponent(constellation,oldCCA,nouveau_name)
        solution.recomputeListOfCCAs(constellation,recomputed)
        
    def linkActivities(self,activity1,activity2,MAJComponents=True):
        if activity1 not in self.vertices[activity2].getNeighbors():
            return self.addEdge(activity1,activity2,MAJComponents)
        
    def getActivities(self):
        return self.getVertices()
    
    def __str__(self):
        res = "Graph of activities' : " + str(self.getNumberOfVertices()) + " activities "
        res += str(self.getNumberOfComponents()) + " components"
        return res

class GraphOfRequest(Graph):
    def __init__(self,request):
        super().__init__("Graph of request "+str(request.getId()),oriented=True)
        self.alpha = "S"
        self.beta = "E"
        self.empty = "X"
        self.activity = ""
        self.addVertices(request)
        self.topologicalOrder()
        #self.partitions = self._computePartitionning(request)
        
    def getPartitions(self):
        return self.partitions
  
    # A recursive function used by topologicalSort
    def topologicalSortUtil(self,v,visited,stack):
 
        # Mark the current node as visited.
        visited[v] = True
 
        # Recur for all the vertices adjacent to this vertex
        for i in self.getVertex(v).getsuccessors():
            if not visited[i]:
                self.topologicalSortUtil(i,visited,stack)
 
        # Push current vertex to stack which stores result
        stack.insert(0,v)
 
    # The function to do Topological Sort. It uses recursive
    # topologicalSortUtil()
    def topologicalOrder(self,inactives=[]):
        # Mark all the vertices as not visited
        visited = {v : False for v in self.getVertices() if v not in inactives}
        stack =[]
 
        # Call the recursive helper function to store Topological
        # Sort starting from all vertices one by one
        for i in self.getVertices():
            if not visited[i] and i not in inactives:
                self.topologicalSortUtil(i,visited,stack)
 
        # Print contents of stack
        self.topological = stack

    def longestPath(self,rewards,inactives):
        rewardsLabels = {self.vertexLabel(v):rewards[v] for v in rewards}
        inactivesLabels = [self.vertexLabel(v) for v in inactives]
        pathCumulatedReward = {self.alpha:(0,[],[])}
        self.topologicalOrder(inactives=inactives)
        for v in self.topological:
            if not v in inactivesLabels:
                for succ in self.getVertex(v).getSuccessors():
                    if succ not in inactivesLabels:
                        currentReward,chem,vertices = pathCumulatedReward.get(succ,(-np.Inf,[],[]))
                        newSuccessor = rewardsLabels.get(succ,0)
                        if currentReward<newSuccessor+pathCumulatedReward[v][0]:
                            vertices = pathCumulatedReward[v][2].copy()
                            elmt = self.getVertex(succ).getData()
                            if len(elmt)!=0:
                                vertices.append(elmt)
                            pathCumulatedReward[succ] = (newSuccessor+pathCumulatedReward[v][0],pathCumulatedReward[v][1]+list(self.getVertex(succ).getData()),vertices)
        try:                    
            return pathCumulatedReward[self.beta]
        except KeyError:
            print("inactive vertices :",inactivesLabels)
            print("graph :",str(self))
            print(self.topological)
            raise KeyError('no path found')

    def testObs(self,v):
        try:
            tuple(v)
            return True
        except Exception:
            return False
    
    def convertNoeud(self,X):
        data = self.getVertex(X).getData()
        return [v for v in data]
    
    def vertexLabel(self,X):
        if type(X)==list or type(X)==tuple:
            return str(X)
        else:
            return self.label(X)
    
    def label(self,o):
        if o==self.alpha or o==self.beta:
            return str(o)
        return o
    
    def __str__(self):
        mess = "Graph of request:"
        for v in self.topological:
            for succ in self.getVertex(v).getSuccessors():
                mess += str(v)+"->"+str(succ)+"\n"
        return mess
        
class GraphMono(GraphOfRequest):
    def __init__(self,request):
        assert(request.getType()== "ONE_SHOT_MONO")
        super().__init__(request)

    def addVertices(self,request):
        self.addVertex(self.alpha,[])
        self.addVertex(self.beta,[])
        for (s,o) in request.getCouples():
            x = o
            self.addVertex(self.vertexLabel((x,)),(x,))
        for v in self.getVertices():
            if v!=self.alpha:
                self.addEdge(self.alpha,v)
            if v!=self.beta:
                self.addEdge(v,self.beta)
    
    def _computePartitionning(self,request):
        partitions = []
        for x in request.getCouples():
            x = o
            partitions.append( [self.vertexLabel((x,))])
        return partitions
            
class GraphStereo(GraphOfRequest):
    def __init__(self,request):
        assert(request.getType()== "ONE_SHOT_STEREO")
        super().__init__(request)
    
    def iterVertices(self,request):
        vertices = []
        for paire in request.getIdPaires():
            paire_obs = request.getPaireObservation(paire)
            (w,s1,o1) = paire_obs[0]
            (w,s2,o2) = paire_obs[1]
            X = (o1,o2)
            vertices.append(X)
        return vertices
    
    def addVertices(self,request):
        self.addVertex(self.alpha,tuple([]))
        self.addVertex(self.beta,tuple([]))
        for X in self.iterVertices(request):
            self.addVertex(self.vertexLabel(X),X)
        for v in self.getVertices():
            if v!=self.alpha:
                self.addEdge(self.alpha,v)
            if v!=self.beta:
                self.addEdge(v,self.beta)

    def _computePartitionning(self,request):
        partitions = []
        for X in self.iterVertices(request):
            partitions.append([self.vertexLabel(X)])
        return partitions
    
class GraphPeriodic(GraphOfRequest):
    def __init__(self,request):
        assert(request.getType()== "PERIODIC")
        super().__init__(request)
    
    def iterVerticesTimeSlot(self,request,timeSlot):
        vertices = []
        activities = request.getObservationsParTimeSlot(timeSlot)
        for w,s,a in activities:
            vertices.append((a,))
        return vertices
        
    def addVertices(self,request):
        self.addVertex(self.alpha,tuple([]))
        self.addVertex(self.beta,tuple([]))
        timeSlots = sorted(list(request.getTimeSlots()))
        # noeuds et aretes entre time slots et extremites
        for i,timeSlot in enumerate(timeSlots):
            self.addVertex("time slot "+str(timeSlot),tuple([]))
            if i==0:
                self.addEdge(self.alpha,"time slot "+str(timeSlots[i])) 
            else:
                self.addEdge("time slot "+str(timeSlots[i-1]),"time slot "+str(timeSlots[i])) 
        for i,timeSlot in enumerate(timeSlots):
            activities = request.getObservationsParTimeSlot(timeSlot)
            for X in self.iterVerticesTimeSlot(request,timeSlot):
                self.addVertex(self.vertexLabel(X),X)
                if i==0:
                    self.addEdge(self.alpha,self.vertexLabel(X))
                else:
                    self.addEdge("time slot "+str(timeSlots[i-1]),self.vertexLabel(X))
                self.addEdge(self.vertexLabel(X),"time slot "+str(timeSlots[i]))
        self.addEdge("time slot "+str(timeSlots[-1]),self.beta)
        
    def _computePartitionning(self,request):
        partitions = []
        timeSlots = sorted(list(request.getTimeSlots()))
        for timeSlot in timeSlots:
            partitions.append([self.vertexLabel(x) for x in self.iterVerticesTimeSlot(request,timeSlot)])
        return partitions
    
class GraphSystematic(GraphOfRequest):
    def __init__(self,request):
        assert(request.getType() == "SYSTEMATIC")
        super().__init__(request)

    def iterVerticesSat(self,request):
        vertices_sat = {}
        for (s,o) in request.getCouples():
            x = o
            if s not in vertices_sat:
                vertices_sat[s] = []
            vertices_sat[s].append((o,))
        return vertices_sat
    
    def addVertices(self,request):
        self.addVertex(self.alpha,())
        self.addVertex(self.beta,())

        iterVertices = self.iterVerticesSat(request)
        for s in iterVertices:
            for X in iterVertices[s]:
                self.addVertex(self.vertexLabel(X),X)
                self.addVertex("not"+self.vertexLabel(X),())
        
        for s in iterVertices:
            prev = [self.alpha]
            for X in iterVertices[s]:
                for p in prev:
                    self.addEdge(p,self.vertexLabel(X))
                    self.addEdge(p,"not"+self.vertexLabel(X))
                prev = [self.vertexLabel(X),"not"+self.vertexLabel(X)]
            for p in prev:
                self.addEdge(p,self.beta)
                
    def _computePartitionning(self,request):
        iterVertices = self.iterVerticesSat(request)
        partitions = []
        for s in iterVertices:
            partitions.append([self.vertexLabel(X) for X in iterVertices[s]])
        return partitions
                
class GraphDownloads(GraphOfRequest):
    def __init__(self,request):
        assert(request.getType() == "DOWNLOAD")
        super().__init__(request)
        
    def iterVertices(self,request):
        vertices = []
        for (s,o) in request.getCouples():
            vertices.append( (o,) )
        return vertices
        
    def addVertices(self,request):
        self.addVertex(self.alpha,[])
        self.addVertex(self.beta,[])
                
        prev = self.alpha
        for X in self.iterVertices(request):
            self.addVertex(self.vertexLabel(X),X)
            self.addEdge(prev,self.vertexLabel(X))
            prev = self.vertexLabel(X)
        self.addEdge(prev,self.beta)

class GraphOfActivities(Graph):
    def __init__(self):
        super().__init__("Graph of activities",oriented=False)
        
    def addActivity(self,activity):
        self.addVertex(activity)
        
    def removeActivity(self,activity,constellation,solution):
        oldComponent,changement,numberOfComponent = self.removeVertex(activity)
        if changement:
            if numberOfComponent==1:
                noeud_compo = solution.renameComponent(constellation,oldComponent)
            else:
                solution.notifyCCASplit(constellation,oldComponent)
    
    def removeListOfActivities(self,activities,constellation,solution):
        recomputed,renamed = self.removeListOfVertices(activities)
        #printColor("renamed "+str(renamed),c='y')
        for oldCCA in renamed:
            newName = renamed[oldCCA]
            solution.renameComponent(constellation,oldCCA,newName)
        solution.recomputeListOfCCAs(constellation,recomputed)
        
    def linkActivities(self,activity1,activity2,MAJComponents=True):
        if activity1 not in self.vertices[activity2].getNeighbors():
            return self.addEdge(activity1,activity2,MAJComponents)
        
    def getActivities(self):
        return self.getVertices()
    
    def __str__(self):
        res = "Graph of activities : " + str(self.getNumberOfVertices()) + " activities "
        res += str(self.getNombrecomponents()) + " components"
        return res
    
class PlanningTask:
    def __init__(self,cca,listOfModes):
        self.cca = cca
        self.listOfModes = listOfModes
    
    def __str__(self):
        return 'PlanningTask(cca:'+str(self.cca)+',modes:'+str(self.listOfModes)+')'
    
    def getCCA(self):
        return self.cca
    
    def getListeModes(self):
        return self.listOfModes
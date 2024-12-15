import numpy as np
from ..Utils.Utils import *

"""
Graphe : utile pour creer des graphes de modes et d'activites,
    l'objectif etant de creer des composantes connexes.
"""
class Graphe:
    class ComposanteConnexe:
        def __init__(self,liste_noeuds):
            self.liste_noeuds = liste_noeuds
            
        def ajouterNoeud(self,noeud):
            #assert(noeud not in self.liste_noeuds)
            self.liste_noeuds.append(noeud)
            
        def retirerNoeud(self,noeud):
            self.liste_noeuds.remove(noeud)
        
        def ajouterListeNoeuds(self,liste_noeuds):
            for noeud in liste_noeuds:
                self.ajouterNoeud(noeud)
        
        def getNoeuds(self):
            return self.liste_noeuds
        
        def getNombreNoeuds(self):
            return len(self.liste_noeuds)
        
        def __str__(self):
            return str(self.liste_noeuds)
        
    class NoeudNonOriente:
        def __init__(self,etiquette):
            self.etiquette = etiquette
            self.composanteConnexe = etiquette
            self.voisins = []
            
        def setComposanteConnexe(self,comp):
            self.composanteConnexe = comp
        
        def getComposanteConnexe(self):
            return self.composanteConnexe
        
        def getEtiquette(self):
            return self.etiquette
        
        def getVoisins(self):
            return self.voisins
        
        def ajouterVoisin(self,s):
            self.voisins.append(s)
            
        def retirerVoisin(self,s):
            self.voisins.remove(s)
        
        def __str__(self):
            return str(self.etiquette)

    class NoeudOriente:
        def __init__(self,etiquette,contenu=None):
            self.etiquette = etiquette
            self.composanteConnexe = etiquette
            self.successeurs = []
            self.predecesseurs = []
            self.contenu = contenu
           
        def getContenu(self):
            return self.contenu
            
        def setComposanteConnexe(self,comp):
            self.composanteConnexe = comp
        
        def getComposanteConnexe(self):
            return self.composanteConnexe
        
        def getEtiquette(self):
            return self.etiquette
        
        def getPredecesseurs(self):
            return self.predecesseurs
        
        def getSuccesseurs(self):
            return self.successeurs
        
        def ajouterSuccesseur(self,s):
            self.successeurs.append(s)
            
        def retirerSuccesseur(self,s):
            self.successeurs.remove(s)
            
        def ajouterPredecesseur(self,s):
            self.predecesseurs.append(s)
            
        def retirerPredecesseur(self,s):
            self.predecesseurs.remove(s)
        
        def __str__(self):
            return str(self.etiquette)    
            
    # classe graphe
    def __init__(self,nom,oriente):
            self.oriente = oriente
            self.nom = nom
            self.noeuds = {}
            #self.aretes = {}
            self.composantes = {}

    #def getAretes(self):
    #    return self.aretes
    
    def __str__(self):
        res = str(len(list(self.noeuds.keys()))) + "noeuds " + str(len(list(self.composantes.keys()))) + " composantes"
        return res
    
    def getComposantesConnexes(self):
        return self.composantes
    
    def getComposanteConnexe(self,c):
        return self.composantes[c]
    
    def getNombreComposantes(self):
        return len(list(self.composantes.keys()))
    
    def noeudExiste(self,etiquette):
        return etiquette in self.noeuds
    
    def ajouterNoeud(self,etiquette,contenu=None):
        assert(etiquette not in self.noeuds)
        if self.oriente:
            self.noeuds[etiquette] = self.NoeudOriente(etiquette,contenu)
        else:
            self.noeuds[etiquette] = self.NoeudNonOriente(etiquette,contenu)
        
        self.composantes[etiquette] = self.ComposanteConnexe([etiquette])
        
    def ajouterArete(self,etiq1,etiq2,label=None,MAJComposantes=True):
        assert(etiq1 != etiq2)
        #assert(etiq1 in self.noeuds)
        #assert(etiq2 in self.noeuds)
        left = etiq1
        right = etiq2
        #self.aretes[(left,right)] = self.Arete(self.noeuds[etiq1],self.noeuds[etiq2],label)
        if self.oriente:
            #assert(right not in self.noeuds[left].getSuccesseurs())
            #assert(left not in self.noeuds[right].getPredecesseurs())
            self.noeuds[left].ajouterSuccesseur(right)
            self.noeuds[right].ajouterPredecesseur(left)
        else:
            #assert(right not in self.noeuds[left].getVoisins())
            #assert(left not in self.noeuds[right].getVoisins())
            self.noeuds[left].ajouterVoisin(right)
            self.noeuds[right].ajouterVoisin(left)
        
        if MAJComposantes:
            c1 = self.noeuds[etiq1].getComposanteConnexe()
            c2 = self.noeuds[etiq2].getComposanteConnexe()
            if c1!=c2:
                self._fusionnerComposantes(c1,c2)
                return (c1,c2)
            else:
                return None
        else:
            return None
        
    def getNom(self):
        return self.nom
    
    def exporter(self):
        with open(self.nom+".dot","w") as file:
            self._ecrireGraphe(file)
    
    def setComposante(self,liste_noeuds,c):
        self.composantes[c] = self.ComposanteConnexe(liste_noeuds)
        for noeud in liste_noeuds:
            self.noeuds[noeud].setComposanteConnexe(c)
        
    def _fusionnerComposantes(self,c1,c2):
        noeuds = self.composantes[c2].getNoeuds()
        for noeud in noeuds:
            self.noeuds[noeud].setComposanteConnexe(c1)
        self.composantes[c1].ajouterListeNoeuds(noeuds)
        del self.composantes[c2]
    
    def getNoeuds(self):
        return self.noeuds
    
    def getNoeud(self,noeud):
        return self.noeuds[noeud]
    
    # retourne : les composantes a recalculer, les composantes a renommer
    def retirerListeNoeuds(self,liste_etiq):
        # creer la liste des noeuds a mettre a jour
        noeuds_connectes = []
        cca_retraits = {}
        for etiq in liste_etiq:
            cca = self.noeuds[etiq].getComposanteConnexe()
            if cca not in cca_retraits:
                cca_retraits[cca] = 0
            cca_retraits[cca] += 1
            if self.oriente:
                for pred in self.noeuds[etiq].getPredecesseurs():
                    self.noeuds[pred].retirerSuccesseur(etiq)
                    if pred not in noeuds_connectes and pred not in liste_etiq:
                        noeuds_connectes.append(pred)
                for succ in self.noeuds[etiq].getSuccesseurs() and succ not in liste_etiq:
                    self.noeuds[succ].retirerPredecesseur(etiq)
                    if succ not in noeuds_connectes:
                        noeuds_connectes.append(succ)
            else:       
                for v in self.noeuds[etiq].getVoisins():
                    self.noeuds[v].retirerVoisin(etiq)
                    if v not in noeuds_connectes and v not in liste_etiq:
                        noeuds_connectes.append(v)
        # effacer les anciennes composantes connexes
        anciennes_composantes = {}
        for etiq in liste_etiq:
            c = self.noeuds[etiq].getComposanteConnexe()
            if c not in anciennes_composantes:
                anciennes_composantes[c] = self.composantes[c].getNombreNoeuds()
                del self.composantes[c]
            del self.noeuds[etiq]
        noeud_traite = []
        changements = [] # les composantes a recalculer (et a renommer)
        rennomages = {} # les composantes a renommer (et pas a recalculer)
        # creer les nouvelles composantes
        while len(noeuds_connectes)>0:
            noeud_courant = noeuds_connectes.pop()
            comp = noeud_courant
            ancienne_composante = self.noeuds[noeud_courant].getComposanteConnexe()
            taille_ancienne_compo = anciennes_composantes[ancienne_composante]
            composante = self.getNoeudsComposante(noeud_courant)
            taille_nouvelle_compo = len(composante)
            if taille_nouvelle_compo<taille_ancienne_compo-cca_retraits[cca]:
                if ancienne_composante not in changements:
                    changements.append(ancienne_composante)
            elif ancienne_composante!=comp:
                rennomages[ancienne_composante] = comp
            # ne pas retraiter les noeuds connectes qui sont dans la meme composante
            del_noeuds = []
            for noeud in noeuds_connectes:
                if noeud in composante:
                    del_noeuds.append(noeud)
            for noeud in del_noeuds:
                noeuds_connectes.remove(noeud)
            self.setComposante(composante,comp)
        return changements,rennomages
    
    # retirer le noeud + les aretes + recalcule les composantes
    def retirerNoeud(self,etiq):
        noeuds_connectes = []
        if self.oriente:
            for pred in self.noeuds[etiq].getPredecesseurs():
                self.noeuds[pred].retirerSuccesseur(etiq)
                if pred not in noeuds_connectes:
                    noeuds_connectes.append(pred)
            for succ in self.noeuds[etiq].getSuccesseurs():
                self.noeuds[succ].retirerPredecesseur(etiq)
                if succ not in noeuds_connectes:
                    noeuds_connectes.append(succ)
        else:       
            for v in self.noeuds[etiq].getVoisins():
                self.noeuds[v].retirerVoisin(etiq)
                if v not in noeuds_connectes:
                    noeuds_connectes.append(v)
        #assert(etiq not in noeuds_connectes)
        ancienne_composante = self.noeuds[etiq].getComposanteConnexe()
        del self.noeuds[etiq]
        del self.composantes[ancienne_composante]
        # il faut recalculer les composantes connexes ici
        noeud_traite = []
        i = 0
        while len(noeuds_connectes)>0:
            noeud_courant = noeuds_connectes.pop()
            comp = noeud_courant#self.noeuds[noeud_courant].getComposanteConnexe()
            assert(noeud_courant!=etiq)
            composante = self.getNoeudsComposante(noeud_courant)
            #assert(etiq not in composante)
            # ne pas retraiter les noeuds connectes qui sont dans la meme composante
            del_noeuds = []
            for noeud in noeuds_connectes:
                if noeud in composante:
                    del_noeuds.append(noeud)
            for noeud in del_noeuds:
                noeuds_connectes.remove(noeud)
            self.setComposante(composante,comp)
            i += 1
        changement = (i>1 or (i==1 and noeud_courant!=ancienne_composante))
        return ancienne_composante,changement,i
    
    def getNoeudsComposante(self,noeud_depart):
        visited = {}
        V = list(self.noeuds.keys())
        for i in V:
            visited[i] = False
        return self.composanteConnexeRec([], noeud_depart, visited)
   
    def composanteConnexeRec(self, temp, v, visited):
        # Mark the current vertex as visited
        visited[v] = True
        # Store the vertex to list
        temp.append(v)
        # Repeat for all vertices adjacent
        # to this vertex v
        if self.oriente:
            adj = self.getNoeud(v).getSuccesseurs() + self.getNoeud(v).getPredecesseurs()
        else:
            adj = self.getNoeud(v).getVoisins()
        for i in adj:
            if visited[i] == False:
                # Update the list
                temp = self.composanteConnexeRec(temp, i, visited)
        return temp
        
    def calculerComposantes(self):
        self.composantes = {}
        visited = {}
        V = list(self.noeuds.keys())
        for v in V:
            visited[v] = False
        for v in V:
            if visited[v] == False:
                temp = []
                self.setComposante(self.composanteConnexeRec(temp, v, visited),v)
                
    def nombreNoeuds(self):
        return len(list(self.noeuds.keys()))

class GrapheActivitesDynamique(Graphe):
    def __init__(self):
        super().__init__("Graphe d'activites",oriente=False)
        
    def ajouterActivite(self,activite):
        self.ajouterNoeud(activite)
        
    def retirerActivite(self,activite,constellation,solution):
        ancienne_composante,changement,nb_comp = self.retirerNoeud(activite)
        if changement:
            if nb_comp==1:
                noeud_compo = solution.changementNomComposante(constellation,ancienne_composante)
            else:
                solution.notifierDivisionCCA(constellation,ancienne_composante)
    
    def retirerListeActivites(self,activites,constellation,solution):
        recalculs,renommages = self.retirerListeNoeuds(activites)
        #printColor("renommages "+str(renommages),c='y')
        for ancienne_cca in renommages:
            nouveau_nom = renommages[ancienne_cca]
            solution.changementNomComposante(constellation,ancienne_cca,nouveau_nom)
        solution.recalculerListeCCA(constellation,recalculs)
        
    def lierActivites(self,activite1,activite2,MAJComposantes=True):
        if activite1 not in self.noeuds[activite2].getVoisins():
            return self.ajouterArete(activite1,activite2,MAJComposantes)
        
    def getActivites(self):
        return self.getNoeuds()
    
    def __str__(self):
        res = "Graphes d'activités' : " + str(self.nombreNoeuds()) + " activités "
        res += str(self.getNombreComposantes()) + " composantes"
        return res

"""    
class EtiqueteRequete:
    def __init__(self):
        pass
    
    def __str__(self):
        return self.contenu

class EtiqueteActivites(EtiqueteRequete):
    def __init__(self,X):
        super().__init__()
        assert(type(X)==list or type(X)==tuple)
        self.contenu = "Activites{"+str(X)+"}"
        
class EtiqueteEvent(EtiqueteRequete):
    def __init_(self):
        super().__init__()
        assert(type(X)==str)
        self.contenu = X
"""

class GrapheRequete(Graphe):
    def __init__(self,requete):
        super().__init__("Graphe de la requête "+str(requete.getId()),oriente=True)
        self.alpha = "S"
        self.beta = "E"
        self.empty = "X"
        self.activite = ""
        self.ajouterNoeuds(requete)
        self.ordreTopologique()
        #self.partitions = self._partitionner(requete)
        
    def getPartitions(self):
        return self.partitions
  
    # A recursive function used by topologicalSort
    def topologicalSortUtil(self,v,visited,stack):
 
        # Mark the current node as visited.
        visited[v] = True
 
        # Recur for all the vertices adjacent to this vertex
        for i in self.getNoeud(v).getSuccesseurs():
            if not visited[i]:
                self.topologicalSortUtil(i,visited,stack)
 
        # Push current vertex to stack which stores result
        stack.insert(0,v)
 
    # The function to do Topological Sort. It uses recursive
    # topologicalSortUtil()
    def ordreTopologique(self,inactifs=[]):
        # Mark all the vertices as not visited
        visited = {v : False for v in self.getNoeuds() if v not in inactifs}
        stack =[]
 
        # Call the recursive helper function to store Topological
        # Sort starting from all vertices one by one
        for i in self.getNoeuds():
            if not visited[i] and i not in inactifs:
                self.topologicalSortUtil(i,visited,stack)
 
        # Print contents of stack
        self.topologique = stack

    def plusLongChemin(self,rewards,inactives):
        rewards_etiquettes = {self.etiquetteNoeud(v):rewards[v] for v in rewards}
        inactives_etiquettes = [self.etiquetteNoeud(v) for v in inactives]
        score_chemin = {self.alpha:(0,[],[])}
        self.ordreTopologique(inactifs=inactives)
        for v in self.topologique:
            if not v in inactives_etiquettes:
                for succ in self.getNoeud(v).getSuccesseurs():
                    if succ not in inactives_etiquettes:
                        score_courant,chem,noeuds = score_chemin.get(succ,(-np.Inf,[],[]))
                        new_succ = rewards_etiquettes.get(succ,0)
                        if score_courant<new_succ+score_chemin[v][0]:
                            noeuds = score_chemin[v][2].copy()
                            elmt = self.getNoeud(succ).getContenu()
                            if len(elmt)!=0:
                                noeuds.append(elmt)
                            score_chemin[succ] = (new_succ+score_chemin[v][0],score_chemin[v][1]+list(self.getNoeud(succ).getContenu()),noeuds)
        try:                    
            return score_chemin[self.beta]
        except KeyError:
            print("noeuds inactifs :",inactives_etiquettes)
            print("graphe :",str(self))
            print(self.topologique)
            raise KeyError('chemin pas trouvé')

    def testObs(self,v):
        try:
            tuple(v)
            return True
        except Exception:
            return False
    
    def convertNoeud(self,X):
        contenu = self.getNoeud(X).getContenu()
        return [v for v in contenu]
    
    def etiquetteNoeud(self,X):
        if type(X)==list or type(X)==tuple:
            return str(X)
        else:
            return self.etiquette(X)
    
    def etiquette(self,o):
        if o==self.alpha or o==self.beta:
            return str(o)
        return o
    
    def __str__(self):
        mess = "Graphe de requête :"
        for v in self.topologique:
            for succ in self.getNoeud(v).getSuccesseurs():
                mess += str(v)+"->"+str(succ)+"\n"
        return mess
        
class GrapheMono(GrapheRequete):
    def __init__(self,requete):
        assert(requete.getType()== "ONE_SHOT_MONO")
        super().__init__(requete)

    def ajouterNoeuds(self,requete):
        self.ajouterNoeud(self.alpha,[])
        self.ajouterNoeud(self.beta,[])
        for (s,o) in requete.getCouples():
            x = o
            self.ajouterNoeud(self.etiquetteNoeud((x,)),(x,))
        for v in self.getNoeuds():
            if v!=self.alpha:
                self.ajouterArete(self.alpha,v)
            if v!=self.beta:
                self.ajouterArete(v,self.beta)
    
    def _partitionner(self,requete):
        partitions = []
        for x in requete.getCouples():
            x = o
            partitions.append( [self.etiquetteNoeud((x,))])
        return partitions
            
class GrapheStereo(GrapheRequete):
    def __init__(self,requete):
        assert(requete.getType()== "ONE_SHOT_STEREO")
        super().__init__(requete)
    
    def iterNoeuds(self,requete):
        noeuds = []
        for paire in requete.getIdPaires():
            paire_obs = requete.getPaireObservation(paire)
            (w,s1,o1) = paire_obs[0]
            (w,s2,o2) = paire_obs[1]
            X = (o1,o2)
            noeuds.append(X)
        return noeuds
    
    def ajouterNoeuds(self,requete):
        self.ajouterNoeud(self.alpha,tuple([]))
        self.ajouterNoeud(self.beta,tuple([]))
        for X in self.iterNoeuds(requete):
            self.ajouterNoeud(self.etiquetteNoeud(X),X)
        for v in self.getNoeuds():
            if v!=self.alpha:
                self.ajouterArete(self.alpha,v)
            if v!=self.beta:
                self.ajouterArete(v,self.beta)
        #print(self.noeuds)
        #print(self.noeuds[self.alpha].getSuccesseurs())

    def _partitionner(self,requete):
        partitions = []
        for X in self.iterNoeuds(requete):
            partitions.append([self.etiquetteNoeud(X)])
        return partitions
    
class GraphePeriodic(GrapheRequete):
    def __init__(self,requete):
        assert(requete.getType()== "PERIODIC")
        super().__init__(requete)
    
    def iterNoeudsTimeSlot(self,requete,timeSlot):
        noeuds = []
        activites = requete.getObservationsParTimeSlot(timeSlot)
        for w,s,a in activites:
            noeuds.append((a,))
        return noeuds
        
    def ajouterNoeuds(self,requete):
        self.ajouterNoeud(self.alpha,tuple([]))
        self.ajouterNoeud(self.beta,tuple([]))
        timeSlots = sorted(list(requete.getTimeSlots()))
        # noeuds et aretes entre time slots et extremites
        for i,timeSlot in enumerate(timeSlots):
            self.ajouterNoeud("time slot "+str(timeSlot),tuple([]))
            if i==0:
                self.ajouterArete(self.alpha,"time slot "+str(timeSlots[i])) 
            else:
                self.ajouterArete("time slot "+str(timeSlots[i-1]),"time slot "+str(timeSlots[i])) 
        for i,timeSlot in enumerate(timeSlots):
            activites = requete.getObservationsParTimeSlot(timeSlot)
            for X in self.iterNoeudsTimeSlot(requete,timeSlot):
                self.ajouterNoeud(self.etiquetteNoeud(X),X)
                if i==0:
                    self.ajouterArete(self.alpha,self.etiquetteNoeud(X))
                else:
                    self.ajouterArete("time slot "+str(timeSlots[i-1]),self.etiquetteNoeud(X))
                self.ajouterArete(self.etiquetteNoeud(X),"time slot "+str(timeSlots[i]))
        self.ajouterArete("time slot "+str(timeSlots[-1]),self.beta)
        
    def _partitionner(self,requete):
        partitions = []
        timeSlots = sorted(list(requete.getTimeSlots()))
        for timeSlot in timeSlots:
            partitions.append([self.etiquetteNoeud(x) for x in self.iterNoeudsTimeSlot(requete,timeSlot)])
        return partitions
    
class GrapheSystematic(GrapheRequete):
    def __init__(self,requete):
        assert(requete.getType() == "SYSTEMATIC")
        super().__init__(requete)

    def iterNoeudsSat(self,requete):
        noeuds_sat = {}
        for (s,o) in requete.getCouples():
            x = o
            if s not in noeuds_sat:
                noeuds_sat[s] = []
            noeuds_sat[s].append((o,))
        return noeuds_sat
    
    def ajouterNoeuds(self,requete):
        self.ajouterNoeud(self.alpha,())
        self.ajouterNoeud(self.beta,())

        iter_noeuds = self.iterNoeudsSat(requete)
        for s in iter_noeuds:
            for X in iter_noeuds[s]:
                self.ajouterNoeud(self.etiquetteNoeud(X),X)
                self.ajouterNoeud("not"+self.etiquetteNoeud(X),())
        
        for s in iter_noeuds:
            prev = [self.alpha]
            for X in iter_noeuds[s]:
                for p in prev:
                    self.ajouterArete(p,self.etiquetteNoeud(X))
                    self.ajouterArete(p,"not"+self.etiquetteNoeud(X))
                prev = [self.etiquetteNoeud(X),"not"+self.etiquetteNoeud(X)]
            for p in prev:
                self.ajouterArete(p,self.beta)
                
    def _partitionner(self,requete):
        iter_noeuds = self.iterNoeudsSat(requete)
        partitions = []
        for s in iter_noeuds:
            partitions.append([self.etiquetteNoeud(X) for X in iter_noeuds[s]])
        return partitions
                
class GrapheVidage(GrapheRequete):
    def __init__(self,requete):
        assert(requete.getType() == "DOWNLOAD")
        super().__init__(requete)
        
    def iterNoeuds(self,requete):
        noeuds = []
        for (s,o) in requete.getCouples():
            noeuds.append( (o,) )
        return noeuds
        
    def ajouterNoeuds(self,requete):
        self.ajouterNoeud(self.alpha,[])
        self.ajouterNoeud(self.beta,[])
                
        prev = self.alpha
        for X in self.iterNoeuds(requete):
            self.ajouterNoeud(self.etiquetteNoeud(X),X)
            self.ajouterArete(prev,self.etiquetteNoeud(X))
            prev = self.etiquetteNoeud(X)
        self.ajouterArete(prev,self.beta)

class GrapheActivites(Graphe):
    def __init__(self):
        super().__init__("Graphe d'activites",oriente=False)
        
    def ajouterActivite(self,activite):
        self.ajouterNoeud(activite)
        
    def retirerActivite(self,activite,constellation,solution):
        ancienne_composante,changement,nb_comp = self.retirerNoeud(activite)
        if changement:
            if nb_comp==1:
                noeud_compo = solution.changementNomComposante(constellation,ancienne_composante)
            else:
                solution.notifierDivisionCCA(constellation,ancienne_composante)
    
    def retirerListeActivites(self,activites,constellation,solution):
        recalculs,renommages = self.retirerListeNoeuds(activites)
        #printColor("renommages "+str(renommages),c='y')
        for ancienne_cca in renommages:
            nouveau_nom = renommages[ancienne_cca]
            solution.changementNomComposante(constellation,ancienne_cca,nouveau_nom)
        solution.recalculerListeCCA(constellation,recalculs)
        
    def lierActivites(self,activite1,activite2,MAJComposantes=True):
        if activite1 not in self.noeuds[activite2].getVoisins():
            return self.ajouterArete(activite1,activite2,MAJComposantes)
        
    def getActivites(self):
        return self.getNoeuds()
    
    def __str__(self):
        res = "Graphes d'activités' : " + str(self.nombreNoeuds()) + " activités "
        res += str(self.getNombreComposantes()) + " composantes"
        return res
    
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
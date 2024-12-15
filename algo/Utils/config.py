import getopt
import sys
import os
import subprocess
from ..model.modelTransition import modeleLent,modeleMoyen,modeleRapide

import numpy as np

class Config:
    def __init__(self):
        pass
class ConfigGlobal(Config):
    def __init__(self):
        self.docplexPath = '/usr/local/CPLEX_Studio1210/cpoptimizer/bin/x86-64_linux/cpoptimizer'
        self.solver = 'OPTW' # ou 'LKH'
        #self.tau_max = 90
        self.allocation_vidage = 3*60
        self.sync = False
        self.tolerance_opti = 0.0
        self.tolerance_temps = 1
        self.periode_affichage_info = 5
        self.digits = 3
        self.score_temporel = False
        # configuration de l'arrivée dynamique des requêtes
    def getEchelle(self):
        return 10**self.digits
    
    def toDict(self):
        return {'digits':self.digits,"durée de vidage":self.allocation_vidage,"solver_cca":self.solver}
    
class ConfigCPSolver(Config):
    def __init__(self):
        self.n_modes = 5
        
    def toDict(self):
        return {"nombre de modes par requête":self.n_modes}
"""    
class ConfigCpSolver(Config):
    def __init__(self):
        self.threads = None
        self.iteration_limit = 0.75*60

    def toDict(self):
        return {'thread':self.thread,'temps par itération':self.iteration_limit}
"""        
class ConfigScoring(Config):
    def __init__(self):
        self.alpha_weather = 0.7
        self.method = 'relative' # relative | global
class ConfigLKH(Config):
    def __init__(self):
        self.max_trials = 3
        self.runs = 3
        self.echelle = 1e3 # 1en = n chiffres apres la virgule
   
    def toDict(self):
        return {'max trials':self.max_trials,'runs':self.runs,'echelle':self.echelle}
class ConfigOPTW(Config):
    def __init__(self):
        self.tmax = 100
        self.restart = 1
        self.echelle = 1e3
        self.tmaxGlobal = 50
        
    def toDict(self):
        return {'tmax':self.tmax,'restart':self.restart,'echelle':self.echelle}
class ConfigUPCCAS(Config):
    def __init__(self):
        pass
    def toDict(self):
        return {}
class ConfigBatchProg(Config):
    def __init__(self):
        pass
    def toDict(self):
        return {}
    
class ConfigLNS(Config):
    def __init__(self):
        self.Tmax = 0
        self.heuristique_reward = 'pure' # ratio ou pure
        self.bruit_modes = 0.0
        self.max_requetes_echange_CP = 20
        #self.voisinage_externe = True
        self.max_destruction = np.Inf
        
    def toDict(self):
        d = {}
        #d["version"]=self.version
        d['bruit des modes']=self.bruit_modes
        #d["utilisation de l'oracle"]=self.use_oracle
        #d['heuristique de reward']=self.heuristique_reward
        #d['méthode de tri des requêtes']=self.ordre_requetes
        d["Tmax"]=self.Tmax
        #d['max stable iteration']=self.maxStableIt
        #d['ratio de destruction']=self.ratio_destruction
        #d['voisinage externe'] = self.voisinage_externe
        return d
    
class Donnees(Config):
    def __init__(self):
        pass
    
    def toDict(self):
        return {'algo':self.algo_name,'instance':self.instance}

class Domain:
    def isDeactivable(self):
        return self.deactivable

class EnumDomain(Domain):
    def __init__(self,values,deactivable=False):
        self.values = values
        self.deactivable = deactivable
        
    def __contains__(self,key):
        return (key is None and self.deactivable) or key in self.values
    
    def evalValue(self,value):
        try:
            _value = str(value)
            if not self.__contains__(_value):
                raise ValueError("Domain error "+str(value)+". Domain : "+self)
            return _value
        except Exception as e:
            print(e)
            raise e
    
class BooleanDomain(Domain):
    def __init__(self,deactivable=False):
        self.deactivable=False
    def __contains__(self, key):
        return (key is None and self.deactivable) or type(key)==bool 

    def evalValue(self,value):
        try:
            if(value is None or value == ""):
                return True # presence du flag sans valeur => True
            _value = bool(value)
            if not self.__contains__(_value):
                raise ValueError("Domain error "+str(value)+". Domain : "+self)
            return _value
        except Exception as e:
            print(e)
            raise e

class IntDomain(Domain):
    def __init__(self,left_bound=-np.Inf,right_bound=np.Inf,deactivable=False):
        self.left_bound = left_bound
        self.right_bound = right_bound
        self.deactivable = deactivable
    
    def tryConvert(self,key):
        try:
            int(key)
            return True
        except:
            return False
        
    def __contains__(self,key):
        if (key is None and self.deactivable):
            return True
        elif self.tryConvert(key) and key>=self.left_bound and key<=self.right_bound:
            return True
        elif np.isinf(key):
            test_left = (np.isinf(self.left_bound) and self.left_bound < 0) or key >= self.left_bound
            test_right = (np.isinf(self.right_bound) and self.right_bound > 0) or key <= self.right_bound
            return test_left and test_right
        else:
            return False

    def evalValue(self,value):
        try:
            _value = int(value)
            if not self.__contains__(_value):
                raise ValueError("Domain error "+str(value)+". Domain : "+self)
            return _value
        except Exception as e:
            print(e)
            raise e

class FloatDomain(Domain):
    def __init__(self,left_bound=-np.Inf,right_bound=np.Inf,deactivable=False):
        self.left_bound = left_bound
        self.right_bound = right_bound
        self.deactivable = deactivable
    
    def tryConvert(self,key):
        try:
            float(key)
            return True
        except:
            return False
    
    def __contains__(self,key):
        if (key is None and self.deactivable):
            return True
        elif self.tryConvert(key) and key>=self.left_bound and key<=self.right_bound:
            return True
        elif np.isinf(key):
            test_left = (np.isinf(self.left_bound) and self.left_bound < 0) or key >= self.left_bound
            test_right = (np.isinf(self.right_bound) and self.right_bound > 0) or key <= self.right_bound
            return  test_left and test_right
        else:
            return False
     
    def evalValue(self,value):
        try:
            _value = float(value)
            if not self.__contains__(_value):
                raise ValueError("Domain error "+str(value)+". Domain : "+self)
            return _value
        except Exception as e:
            print(e)
            raise e
    
class StrDomain(Domain):
    def __init__(self,values,deactivable=False):
        self.values = values
        self.deactivable = deactivable
    
    def getValues(self):
        return self.values        
    def __contains__(self,key):
        return self.values is None or (self.deactivable and key is None) or type(key)==str  and key in self.values

    def evalValue(self,value):
        try:
            _value = str(value)
            if not self.__contains__(_value):
                raise ValueError("Domain error "+str(value)+". Domain : "+self)
            return _value
        except Exception as e:
            print(e)
            raise e
            
class PathDomain(Domain):
    def __init__(self,deactivable=False):
        self.deactivable = deactivable
        
    def __contains__(self,key):
        return (self.deactivable and key is None) or type(key)==str 

    def evalValue(self,value):
        try:
            _value = str(value)
            if not self.__contains__(_value):
                raise ValueError("Domain error "+str(value)+". Domain : "+self)
            return _value
        except Exception as e:
            print(e)
            raise e
            
class Option:
    def __init__(self,short,name,algos,defaut_value=None,domain=None,description=None):
        self.short = short
        self.name = name
        self.algos = algos
        self.value = defaut_value
        self.description = description
        if domain is None:
            self.inferType(defaut_value)
        else:
            self.domain = domain
        assert(isinstance(self.domain,Domain))
        if not self.domain.isDeactivable():
            assert(defaut_value is not None)
        if self.value is not None:
            assert(self.checkValueInDomain())
    
    def getDomain(self):
        return self.domain
    
    def getDescription(self):
        return self.description
    
    def setValue(self,key):
        try:
            self.value = self.domain.evalValue(key)
        except Exception as e:
            print(self.name + " : " + e)
            raise e

    def getValue(self):
        return self.value
    
    def isFlag(self):
        return isinstance(self.domain,BooleanDomain)
        
    def isActive(self):
        return self.value is not None
     
    def checkValueInDomain(self):
        return self.value in self.domain
        
    def inferType(self,defaut_value):
        if type(defaut_value)==int:
            self.domain = IntDomain()
        elif type(defaut_value)==bool:
            self.domain = BooleanDomain()
        elif type(defaut_value)==float:
            self.domain = FloatDomain()
        else:
            raise ValueError(self.name,"Precise option domain.")
    
    def getValue(self):
        return self.value
    
    def getName(self):
        return self.name
    
    def getShortName(self):
        return self.short
    
    def getAlgos(self):
        return self.algos.copy()
    
    def estFlag(self):
        return self.flag
    
    def estValeur(self):
        return not self.flag

class Config:
    def __init__(self):
        self.glob = ConfigGlobal()
        self.CPSolver = ConfigCPSolver()
        self.UPCCAS = ConfigUPCCAS()
        #self.cpSolver = ConfigCpSolver()
        self.lkh = ConfigLKH()
        self.donnees = Donnees()
        self.batchProg = ConfigBatchProg()
        self.scoring = ConfigScoring()
        self.LNS = ConfigLNS()
        self.OPTW = ConfigOPTW()
        self.algos = ['LNS','BPCCAS','UPCCAS','CPSolver','cpSolver','dataVisu','timeDependentSolver']
        self.commentaire = None
        try:
            self.instance = self.parseParametres()
        except Exception as e:
            print(e)
            
    def restrict(self,o,restrict_option):
        if o in restrict_option:
            testAlgo = self.getOptValue("solver") in restrict_option[o]
            testSolverTimeDep =  self.getOptValue("solver")=="timeDependentSolver" and "LNS" in [self.getOptValue("initSolver"),self.getOptValue("fillSolver")]
            if not testAlgo and not testSolverTimeDep:
                raise ValueError("Option inconnue (pour cet algorithme) : "+str(o))
    
    def getInstance(self):
        return self.instance
    
    def verifMode(self):
        return self.getOptValue("verif")
        
    def restreindreOptions(self):
        for opt in self.options:
            short_name,long = self.options[opt].getShortName(),self.options[opt].getName()
            algos = self.options[opt].getAlgos()
            self.restrict_options["--"+long]=algos
            if short_name is not None:
                self.restrict_options["-"+short_name]=algos
        
    def ajouterOption(self,short_name,nom,algos,defaut_value,domain=None,description=None):
        self.options[nom] = Option(short_name,nom,algos,defaut_value,domain=domain,description=description)
    
    def getOption(self,name):
        return self.options[name]
    
    def initOptions(self):
        self.options = {}
        # name,algos,defaut_value=None,domain=None
        self.ajouterOption(None,"step",["timeDependentSolver","BPCCAS","UPCCAS","LNS"],False,description="Activer le mode pas à pas (quand l'algo comprend des 'steps').")
        self.ajouterOption(None,"cpu",["timeDependentSolver","BPCCAS","UPCCAS","LNS"],False,description="Activer pour créer le graphique de charge de coeurs.")
        self.ajouterOption(None,"charge",["timeDependentSolver","BPCCAS","UPCCAS","LNS","CPSolver"],False,description="Graphe de charge des composantes.")
        self.ajouterOption(None,"obj",["timeDependentSolver","BPCCAS","UPCCAS","LNS","CPSolver"],False,description="Activer pour créer la courbe d'évolution de l'objectif.")
        self.ajouterOption(None,"sample",["timeDependentSolver","BPCCAS","UPCCAS","LNS","CPSolver"],None,PathDomain(deactivable=True),description="Indiquer un chemin pour sauvegarder les résultats.")
        self.ajouterOption("v","verbose",["timeDependentSolver","BPCCAS","UPCCAS","LNS","CPSolver","dataVisu"],0,IntDomain(left_bound=0),description="Indiquer le niveau d'affichage.")
        self.ajouterOption("o","verif",["timeDependentSolver","BPCCAS","UPCCAS","LNS","CPSolver"],False,description="Activer pour effectuer des vérifications en cours d'execution (plus lent, plus sûr).")
        self.ajouterOption(None,"full_sample",["timeDependentSolver","BPCCAS","UPCCAS","LNS","CPSolver"],False,description="Activer pour alléger la quantité d'information sauvegardée.")
        self.ajouterOption("t","time",["timeDependentSolver","BPCCAS","UPCCAS","LNS","CPSolver"],5,FloatDomain(left_bound=0),description="Indiquer le temps d'execution (en secondes).")
        # solvers CP
        self.ajouterOption("m","modes",["CPSolver"],5,IntDomain(left_bound=1),description="Indiquer le nombre de modes considérés par requêtes.")
        self.ajouterOption("w","threads",["CPSolver","cpSolver"],1,IntDomain(left_bound=1),description="Indiquer le nombre de threads travailleurs.")
        # BPCCAS
        self.ajouterOption(None,"conserve",["BPCCAS"],False,description="Mode conservation du BPCCAS (lent) : conserve les modes validées même s'ils perdent en rang (score).")
        self.ajouterOption(None,"stable_exp",["BPCCAS"],False,description="Mode stable explication (lent) : rejette un mode avec explication seulement si aucune autre requête mieux classé et non explorée existe.")
        self.ajouterOption(None,"fails",["BPCCAS"],5,IntDomain(left_bound=1),description="Nombre d'éches consecutifs maximum pour l'appel au solver sur un coeur.")
        self.ajouterOption(None,"scalls",["BPCCAS"],50,IntDomain(left_bound=1),description="Nombre d'appels au solver maximum sur un coeur.")
        self.ajouterOption(None,"solver_PCCAS",["BPCCAS","UPCCAS"],"OPTWGroups",StrDomain(["LKH","OPTW","OPTWGroups"]),description="Solver pour BPCCAS et UPCCAS.")
        #LNS
        self.ajouterOption("k","n_cca",["LNS"],2,IntDomain(left_bound=1),description="Nombre de CCA considérées dans le voisinage (version LNS avec CP).")
        self.ajouterOption("u","use_solver",["LNS"],False,description="Activer pour activer le solver à chaque insertion temporaire (mode glouton).")
        self.ajouterOption(None,"quota_requests",["LNS"],np.Inf,IntDomain(left_bound=1,right_bound=np.Inf),description="Nombre de requêtes à détruire au maximum (LNS hybride).")
        self.ajouterOption(None,"nbh", ["LNS"], "load",StrDomain(["random","load"]),description="Mode de séléction des CCAs durant la recherche locale.")
        self.ajouterOption(None,"version",["LNS"],"greedy-request",StrDomain(["greedy-request","greedy-cca","hybrid","coop"]),description="Choix de la version du LNS.")
        self.ajouterOption(None,"time_insert",["LNS"],np.Inf,FloatDomain(left_bound=0,deactivable=True),description="Temps de temps maximum passé à l'insertion d'un mode (s). (Mode glouton)")
        self.ajouterOption(None,"tries_insert",["LNS"],np.Inf,IntDomain(left_bound=0,deactivable=True),description="Nombre d'insertions maximum considérés à l'étape d'insertion. (Mode glouton)")
        self.ajouterOption("i","stableIt",["LNS"],15,IntDomain(left_bound=1),description="Nombre d'itération stable avant perturbation durant la recherche locale.")
        self.ajouterOption(None,"max_operateur",["LNS"],np.Inf,IntDomain(left_bound=1),description="Nombre maximum d'appels de l'opérateur.")
        self.ajouterOption(None,"max_iteration",["LNS"],np.Inf,IntDomain(left_bound=1),description="Nombre d'itérations maximum.")
        self.ajouterOption(None,"CPLim",["LNS"],3.0,FloatDomain(left_bound=0),description="Temps alloué au solver CP (version LNS avec CP).")
        self.ajouterOption(None,"max_echecs",["LNS"],np.Inf,IntDomain(left_bound=0),description="Nombre d'echecs maximum du glouton.")
        # Time dependent solver
        self.ajouterOption(None,"initSolver",["timeDependentSolver"],"LNS",StrDomain(["LNS","CPSolver","BPCCAS"]),description="Solver d'initialisation.")
        self.ajouterOption(None,"initDur",["timeDependentSolver"],np.Inf,FloatDomain(left_bound=0),description="Temps d'initialisation.")
        self.ajouterOption(None,"fillSolver",["timeDependentSolver"],"LNS",StrDomain(["LNS","CPSolver","BPCCAS"]),description="Solver d'initialisation.")
        self.ajouterOption(None,"initTransitionModel",["timeDependentSolver","BPCCAS","UPCCAS","LNS","CPSolver"],"mean",StrDomain(["slow","fast","mean","time-dep"]),description="Modèle de transition : (slow|fast|mean).")
        self.ajouterOption(None,"chain",["timeDependentSolver"],False,BooleanDomain(),description="Flag pour activer l'enchaînement des modèles approchés et exactes.")
        
        #self.ajouterOption("o","ordering",["LNS"],True)
        #self.ajouterOption("H","heuristique",["LNS"],True)
        # UPCCAS
        self.ajouterOption(None,"preval",["UPCCAS"],False,description="Mode prévalidation : transfère les activités déjà validées d'un nouveau mode (celles du mode précédent). Courcircuite potentiellement le rang de scores.")
        self.ajouterOption("a","anticipation",["UPCCAS"],False,description="Anticipe la planification d'un mode moins bien classé si ceux mieux classés sont déjà planifiés sur la CCA en question.")
        # Autres
        self.ajouterOption(None,"data",["dataVisu","timeDependentSolver","UPCCAS","BPCCAS","LNS","CPSolver"],"time-independent",StrDomain(["time-independent","time-dependent"]),description="Choix du dossier des données d'entrées.")
        #self.ajouterOption(None,"flat",["LNS","BPCCAS","UPCCAS","CPSolver"],False,BooleanDomain(),description="Indiquer que l'arborescence des données est aplatie (un 1 seul dossier parent).")
        self.ajouterOption(None,"solver",["timeDependentSolver","UPCCAS","BPCCAS","LNS","CPSolver"],"LNS",StrDomain(["timeDependentSolver","UPCCAS","BPCCAS","LNS","CPSolver"]),description="Solver à utiliser.") # option de bruitage
        self.ajouterOption(None,"derniere_requete",["UPCCAS","BPCCAS","LNS"],4.5,FloatDomain(left_bound=0),description="Date de la dernière requête en minutes.") # option de bruitage
        self.ajouterOption(None,"proportion_depart",["UPCCAS","BPCCAS","LNS"],0.6,FloatDomain(left_bound=0,right_bound=1),description="Proportion de requêtes au départ (cas dynamique).") # option de bruitage
        self.ajouterOption(None,"dt_requetes",["UPCCAS","BPCCAS","LNS"],0.5,FloatDomain(left_bound=0),description="Durée entre l'arrivée des batchs de requêtes.") # option de bruitage       
        self.ajouterOption(None,"noise",["UPCCAS","BPCCAS"],0,FloatDomain(left_bound=0),description="Amplitude du bruit (itération stochastique).") # option de bruitage
        self.ajouterOption("d","destroy",["LNS"],0.2,FloatDomain(left_bound=0,right_bound=1),description="Amplitude de destruction (itération stochastique).") # ratio de destruction
        self.ajouterOption(None,"perturb_rate",["LNS"],0.2,FloatDomain(left_bound=0,right_bound=1),description="Amplitude de destruction (perturbation).") # ratio de destruction
        self.ajouterOption("d","destroy_BPCCAS",["BPCCAS"],0,FloatDomain(left_bound=0,right_bound=1),description="Amplitude de destruction (itération stochastique).") # ratio de destruction
        self.ajouterOption("h","help",["timeDependentSolver","BPCCAS","UPCCAS","LNS","CPSolver","dataVisu"],False,description="Afficher l'aide.")# help
        self.ajouterOption("c","comm",["BPCCAS","UPCCAS","LNS","CPSolver"],"",StrDomain(None,deactivable=True),"Ajouter un commentaire à l'échantillon de résultats.")# help
        self.ajouterOption("f","folder",["timeDependentSolver","BPCCAS","UPCCAS","LNS","CPSolver","dataVisu"],None,IntDomain(left_bound=0,deactivable=True),"Indiquer le dossier de l'instance.")# help
        self.ajouterOption("n","file",["timeDependentSolver","BPCCAS","UPCCAS","LNS","CPSolver","dataVisu"],None,IntDomain(left_bound=0,deactivable=True),"Indiquer le fichier de l'instance.")# help
        self.ajouterOption("s","seed",["BPCCAS","UPCCAS","LNS","CPSolver","dataVisu"],0,IntDomain(),description="Graine aléatoire.")
        self.ajouterOption(None,"include_systematic",["BPCCAS","UPCCAS","LNS","CPSolver","dataVisu"],False,BooleanDomain(),description="Inclure les requêtes systématiques.")
        self.ajouterOption(None,"restart",["LNS","BPCCAS","UPCCAS"],1,IntDomain(left_bound=0),description="Nombre de restarts du solver OPTW (algos avec solver OPTW).")
        self.ajouterOption(None,"dynamic",["LNS","BPCCAS","UPCCAS","CPSolver"],False,BooleanDomain(),description="Activer l'arrivée dynamique des requêtes.")
        self.ajouterOption(None,"test_req",["LNS","BPCCAS","UPCCAS","CPSolver"],np.Inf,IntDomain(left_bound=1),description="Limiter le nombre de requêtes (pour le debug).")
        self.ajouterOption(None,"profile",["LNS","BPCCAS","UPCCAS","CPSolver"],False,BooleanDomain(),description="Activer le profiling du modèle de transition.")

        #self.ajouterOption(None,"vlat",["LNS","BPCCAS","UPCCAS","CPSolver"],1,FloatDomain(left_bound=0),description="Vitesse angulaire latitude.")
        #self.ajouterOption(None,"vlong",["LNS","BPCCAS","UPCCAS","CPSolver"],1,FloatDomain(left_bound=0),description="Vitesse angulaire longitude.")
        
        # limiter les options aux algos indiqués
        self.restreindreOptions()
    
    def isTimeDependentModeOn(self):
        return self.getOptValue("initTransitionModel")=="time-dep" or self.getOptValue("solver")=="timeDependentSolver"
    
    def isDataTimeDependent(self):
        return self.getOptValue("data")=="time-dependent"
    
    def afficherAide(self):
        for opt in self.options:
            if self.getOptValue("solver") in self.options[opt].getAlgos():
                if self.options[opt].getShortName() is not None:
                    print("-"+self.options[opt].getShortName(),end=" || ")
                print("--"+ self.options[opt].getName()+" : "+self.options[opt].getDescription())
        
    def monterFlagOption(self):
        short_listing = ""
        long_listing = []
        for opt in self.options:
            short = self.options[opt].getShortName()
            if self.options[opt].isFlag():
                if short is not None:
                    short_listing += short
                long_listing.append(opt)
            else:
                if short is not None:
                    short_listing += short+":"
                long_listing.append(opt+"=")
        return short_listing,long_listing
    
    def setOptValue(self,o,value):
        if o in self.options:
            self.options[o].setValue(value)
        else:
            opt = self.findOpt(o)
            self.options[opt].setValue(value)
    
    def findOpt(self,o):
        try_long = o.split("--")
        if len(try_long)>1 and try_long[1] in self.options:
            return try_long[1]
        else:
            try_short = o.split("-")
            assert(len(try_short)>0)
            oo = try_short[1]
            for opt in self.options:
                if self.options[opt].getShortName()==oo:
                    return opt
            raise ValueError("Options not found : "+str(o))
    
    def getOptValue(self,key):
        return self.options[key].getValue()
    
    # temps de calcul : 5mn par defaut
    # instance : num instance si spécifié, -1 pour toutes les instances, None sinon
    def parseParametres(self):
        instance = None
        
        self.restrict_options = {}
        self.initOptions()
        short_listing,long_listing = self.monterFlagOption()
        opts, args = getopt.getopt(sys.argv[1:], short_listing, long_listing)
        
        self.donnees.algo_name = None
        for o, a in opts:
            if o=="--solver":
                self.donnees.algo_name = a
                break
        if self.donnees.algo_name is None:
            self.donnees.algo_name = self.getOptValue("solver")
            
        #print(opts,args)
        op = [x[0] for x in opts]
        if '--destroy' in op and '--noise' in op:
            raise Exception('--destroy et --noise ne peuvent être utilisés simultanémment')
        for o, a in opts:
            self.restrict(o,self.restrict_options)
            opt = self.findOpt(o)
            self.options[opt].setValue(a)
        #self.options["time"].setValue(self.options["time"].getValue()*60)               
        # traiter les incohérences d'options avec le solver utilisé
        if not(self.getOptValue("solver") in self.algos):
            raise NameError("Algorithme inconnu",self.donnees.algo_name)
        
        instance = {"folder":self.getOptValue('folder'),"file":self.getOptValue('file')}
        if instance["folder"] is None:
            instance = None
        self.donnees.instance = instance
        if instance is not None and 'file' in instance and 'folder' not in instance:
            raise ValueError("Spécifiez le dossier de l'instance")
        
        if self.getOptValue("solver")=="timeDependentSolver":
            if not self.getOptValue("data")=="time-dependent":
                raise ValueError("Utilisez les instances time-dependent avec le solver timeDependentSolver. --data=time-dependent")
        
        if self.getOptValue("chain"):
            if not self.getOptValue("initTransitionModel") in ["slow","mean","fast"]:
                raise ValueError("Définissez un modèle approché avec --initTransitionModel=<model>. Valeurs possibles : "+str(["slow","mean","fast"]))
            if self.getOptValue("initDur")>=self.getOptValue("time"):
                raise ValueError("Temps alloué au solver approché dépassant la durée totale. Indiquez --initDur=<duree>.")
        return instance
    
    def __str__(self):
        mess = " ------------ configuration globale ---------------------\n"
        for opt in self.options:
            if self.donnees.algo_name in self.options[opt].getAlgos():
                mess += " | " + opt + ' : ' + str(self.options[opt].getValue()) + '\n'
        
        mess += " ------------ configuration globale ---------------------\n"
        globDict = self.glob.toDict()
        for key,value in globDict.items():
            mess +=" | " + key + ' : ' + str(value) + '\n'
        mess += " ------------ infos générale ---------------------------\n"
        for key,value in self.donnees.toDict().items():
            mess += " | " + key + " : " + str(value) + '\n'
        if self.donnees.algo_name == "CPSolver":
            mess += " ------------ configuration solver global ---------------\n"
            dictGlobalSol = self.CPSolver.toDict()
            for key,value in dictGlobalSol.items():
                mess +=" | " + key + ' : ' + str(value) + '\n'
        """
        if self.donnees.algo_name == "cpSolver":
            dictCp = self.cpSolver.toDict()
            mess += " ------------ configuration solver iteratif -------------\n"
            for key,value in dictCp.items():
                mess +=" | " + key + ' : ' + str(value) + '\n'
        """
        if self.donnees.algo_name == "LNS":
            dictLNS = self.LNS.toDict()
            mess += " ------------ configuration LNS -------------------------\n"
            for key,value in dictLNS.items():
                mess +=" | " + key + ' : ' + str(value) + '\n'
        if self.donnees.algo_name == "BPCCAS":
            dictBatch = self.batchProg.toDict()
            mess += " ------------ configuration BPCCAS ----------------------\n"
            for key,value in dictBatch.items():
                mess +=" | " + key + ' : ' + str(value) + '\n'
        if self.donnees.algo_name == 'UPCCAS':
            mess += " ------------ configuration UPCCAS -----------------------\n"
            for key,value in self.UPCCAS.toDict().items():
                mess +=" | " + key + ' : ' + str(value) + '\n'            
        if self.donnees.algo_name in ['BPCCAS','UPPCAS'] or self.donnees.algo_name == 'LNS':
            if self.glob.solver == 'OPTW':
                dictOPTW = self.OPTW.toDict()
                mess += " ------------ configuration OPTW --------------------------\n"
                for key,value in dictOPTW.items():
                    mess +=" | " + key + ' : ' + str(value) + '\n'
            else:
                dictLKH = self.lkh.toDict()
                mess += " ------------ configuration LKH --------------------------\n"
                for key,value in dictLKH.items():
                    mess +=" | " + key + ' : ' + str(value) + '\n'
        mess += " ------------ modèle de transition --------------------------\n"
        if self.getOptValue("initTransitionModel")=="fast":
            transition_model = modeleRapide
        elif self.getOptValue("initTransitionModel")=="mean":
            transition_model = modeleMoyen
        elif self.getOptValue("initTransitionModel")=="slow":
            transition_model = modeleLent
        else:
            transition_model = "modèle time-dependent"
        mess += " | modèle initial : " + str(transition_model) +"\n"
        mess += " ------------ OPTIONS -----------------------\n"
        for opt in self.options:
            mess += " | " + str(opt) + " : " + str(self.getOptValue(opt)) + "\n"
        return mess
            
        
        
        
        
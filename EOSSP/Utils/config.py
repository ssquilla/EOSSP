import getopt
import sys
import os
import subprocess
from EOSSP.model.transitionModel import SlowModel,MeanModel,FastModel

import numpy as np

class Config:
    def __init__(self):
        pass
class ConfigGlobal(Config):
    def __init__(self):
        #self.docplexPath = '/usr/local/CPLEX_Studio1210/cpoptimizer/bin/x86-64_linux/cpoptimizer'
        self.docplexPath = ''
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
    def getScale(self):
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
        d['bruit des modes']=self.bruit_modes
        d["Tmax"]=self.Tmax
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
                raise ValueError("Domain error "+str(value)+". Domain: "+self)
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
                raise ValueError("Domain error "+str(value)+". Domain: "+self)
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
                raise ValueError("Domain error "+str(value)+". Domain: "+self)
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
                raise ValueError("Domain error "+str(value)+". Domain: "+self)
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
                raise ValueError("Domain error "+str(value)+". Domain: "+self)
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
                raise ValueError("Domain error "+str(value)+". Domain: "+self)
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
        
    def isValeur(self):
        return not self.isFlag()

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
                raise ValueError("Unknwown option for the selected algorithm: "+str(o))
    
    def getInstance(self):
        return self.instance
    
    def verifMode(self):
        return self.getOptValue("verif")
        
    def applyOptionsParsingRules(self):
        for opt in self.options:
            short_name,long = self.options[opt].getShortName(),self.options[opt].getName()
            algos = self.options[opt].getAlgos()
            self.restrict_options["--"+long]=algos
            if short_name is not None:
                self.restrict_options["-"+short_name]=algos
        
    def declareOption(self,short_name,nom,algos,defaut_value,domain=None,description=None):
        self.options[nom] = Option(short_name,nom,algos,defaut_value,domain=domain,description=description)
    
    def getOption(self,name):
        return self.options[name]
    
    def initOptions(self):
        self.options = {}
        # name,algos,defaut_value=None,domain=None
        self.declareOption(None,"step",["timeDependentSolver","BPCCAS","UPCCAS","LNS"],False,description="Activate step by step mode (when compatible)")
        self.declareOption("v","verbose",["timeDependentSolver","BPCCAS","UPCCAS","LNS","CPSolver","dataVisu"],0,IntDomain(left_bound=0),description="Indicates the display depth.")
        self.declareOption("o","verif",["timeDependentSolver","BPCCAS","UPCCAS","LNS","CPSolver"],False,description="Activate to proceed the implemented checks during execution (slower, safer).")
        self.declareOption(None,"full_sample",["timeDependentSolver","BPCCAS","UPCCAS","LNS","CPSolver"],False,description="Activate to add supplementary information in samples.")
        self.declareOption("t","time",["timeDependentSolver","BPCCAS","UPCCAS","LNS","CPSolver"],5,FloatDomain(left_bound=0),description="Set the time limit (in seconds).")
        self.declareOption(None,"sample",["timeDependentSolver","BPCCAS","UPCCAS","LNS","CPSolver"],None,PathDomain(deactivable=True),description="Indicates a path to save the results.")
        self.declareOption(None,"obj",["timeDependentSolver","BPCCAS","UPCCAS","LNS","CPSolver"],False,description="Export the objective curve.")
        self.declareOption(None,"cpu",["timeDependentSolver","BPCCAS","UPCCAS","LNS"],False,description="Export the load of each process.")
        self.declareOption(None,"ccaload",["timeDependentSolver","BPCCAS","UPCCAS","LNS","CPSolver"],False,description="Export the load of connected components.")
        # solvers CP
        self.declareOption("m","modes",["CPSolver"],5,IntDomain(left_bound=1),description="Indicates the number of modes considered by request.")
        self.declareOption("w","threads",["CPSolver","cpSolver"],1,IntDomain(left_bound=1),description="Indicates the number of workers")
        self.declareOption(None,"cplex",["CPSolver","LNS"],self.glob.docplexPath,PathDomain(),description="Indicates the Cplex path.")
        #self.declareOption(None,"cplex",["CPSolver","LNS"],"",PathDomain(),description="Indicates the Cplex path.")
        # BPCCAS
        self.declareOption(None,"conserve",["BPCCAS"],False,description="Conservation mode for BPCCAS (slower): conserve all validated modes even if they violated the greedy order.")
        self.declareOption(None,"stable_exp",["BPCCAS"],False,description="Stable explication mode (slower): reject modes with explaination only if another unexplored and better ranked requests exists.")
        self.declareOption(None,"fails",["BPCCAS"],5,IntDomain(left_bound=1),description="Maximum number of insertion failures within a solver execution (for each process).")
        self.declareOption(None,"scalls",["BPCCAS"],50,IntDomain(left_bound=1),description="Maximum number of call to the solver (for each process).")
        self.declareOption(None,"solver_PCCAS",["BPCCAS","UPCCAS"],"OPTW",StrDomain(["LKH","OPTW","OPTWGroups"]),description="External TSPTW choice for BPCCAS and UPCCAS.")
        #LNS
        self.declareOption("k","n_cca",["LNS"],2,IntDomain(left_bound=1),description="Number of CCAs considered to build a neighborhood (LNS solver hybridized with CP).")
        self.declareOption("u","use_solver",["LNS"],False,description="Activate to exploit the external TSPTW solver during activities insertions.")
        self.declareOption(None,"quota_requests",["LNS"],np.Inf,IntDomain(left_bound=1,right_bound=np.Inf),description="Maximum number of requests to destroy (LNS solver hybridized with CP).")
        self.declareOption(None,"nbh", ["LNS"], "load",StrDomain(["random","load"]),description="Heuristic to select CCAs during local search (random|load: load profile measure of the CCA).")
        self.declareOption(None,"version",["LNS"],"greedy-request",StrDomain(["greedy-request","greedy-cca","hybrid","coop"]),description="Select the LNS version.")
        self.declareOption(None,"time_insert",["LNS"],np.Inf,FloatDomain(left_bound=0,deactivable=True),description="Maximum time spent to insert a given mode (s). (Greedy mode)")
        self.declareOption(None,"tries_insert",["LNS"],np.Inf,IntDomain(left_bound=0,deactivable=True),description="Maximum number of insertions considered at each iteration. (Greedy mode)")
        self.declareOption("i","stableIt",["LNS"],15,IntDomain(left_bound=1),description="Maximum number of stable iterations before performing a perturbation during local search.")
        self.declareOption(None,"max_operateur",["LNS"],np.Inf,IntDomain(left_bound=1),description="Maximum number of calls to the operators (Meant for debug).")
        self.declareOption(None,"max_iteration",["LNS"],np.Inf,IntDomain(left_bound=1),description="Maximum number of iterations (Meant for debug).")
        self.declareOption(None,"CPLim",["LNS"],3.0,FloatDomain(left_bound=0),description="Time limit allocated to the CP solver (LNS hybridized with CP).")
        self.declareOption(None,"max_echecs",["LNS"],np.Inf,IntDomain(left_bound=0),description="Maximum number failed insertions during each greedy procedure.")
        # Time dependent solver
        self.declareOption(None,"initSolver",["timeDependentSolver"],"LNS",StrDomain(["LNS","CPSolver","BPCCAS"]),description="Initial solver.")
        self.declareOption(None,"initDur",["timeDependentSolver"],np.Inf,FloatDomain(left_bound=0),description="Initial duration.")
        self.declareOption(None,"fillSolver",["timeDependentSolver"],"LNS",StrDomain(["LNS","CPSolver","BPCCAS"]),description="Insertion solver.")
        self.declareOption(None,"initTransitionModel",["timeDependentSolver","BPCCAS","UPCCAS","LNS","CPSolver"],"mean",StrDomain(["slow","fast","mean","time-dep"]),description="Transition model: (slow|fast|mean).")
        self.declareOption(None,"chain",["timeDependentSolver"],False,BooleanDomain(),description="Activate to chain approximated and exact model (Hybrid smodel strategy).")
        
        # UPCCAS
        self.declareOption(None,"preval",["UPCCAS"],False,description="Prevalidation mode: transfer activities that are already planned of a failed mode to the new generated mode. Potentially violates the greedy order.")
        self.declareOption("a","anticipation",["UPCCAS"],False,description="Lock the planning of mode that are not top-ranked w.r.t the greedy mode (considering the non-generated modes also). Prevent greedy order violation.")
        # Autres
        self.declareOption(None,"data",["dataVisu","timeDependentSolver","UPCCAS","BPCCAS","LNS","CPSolver"],"time-independent",StrDomain(["time-independent","time-dependent"]),description="Set the input data folder.")
        self.declareOption(None,"solver",["timeDependentSolver","UPCCAS","BPCCAS","LNS","CPSolver"],"LNS",StrDomain(["timeDependentSolver","UPCCAS","BPCCAS","LNS","CPSolver"]),description="Solver to use.") # option de bruitage
        self.declareOption(None,"derniere_requete",["UPCCAS","BPCCAS","LNS"],4.5,FloatDomain(left_bound=0),description="Delay before the first requests arrival (in minutes).") # option de bruitage
        self.declareOption(None,"proportion_depart",["UPCCAS","BPCCAS","LNS"],0.6,FloatDomain(left_bound=0,right_bound=1),description="Initial proportion of available requests (dynamic case).") # option de bruitage
        self.declareOption(None,"dt_requetes",["UPCCAS","BPCCAS","LNS"],0.5,FloatDomain(left_bound=0),description="Duration between the arrival of each batch of requests.") # option de bruitage       
        self.declareOption(None,"noise",["UPCCAS","BPCCAS"],0,FloatDomain(left_bound=0),description="Noise amplitude (stochastic iterations).") # option de bruitage
        self.declareOption("d","destroy",["LNS"],0.2,FloatDomain(left_bound=0,right_bound=1),description="Destruction amplitude (stochastic iterations).") # ratio de destruction
        self.declareOption(None,"perturb_rate",["LNS"],0.2,FloatDomain(left_bound=0,right_bound=1),description="Destruction amplitude (perturbation).") # ratio de destruction
        self.declareOption("d","destroy_BPCCAS",["BPCCAS"],0,FloatDomain(left_bound=0,right_bound=1),description="Destruction amplitude (stochastic iteration).") # ratio de destruction
        self.declareOption("h","help",["timeDependentSolver","BPCCAS","UPCCAS","LNS","CPSolver","dataVisu"],False,description="Display help menu.")# help
        self.declareOption("c","comm",["BPCCAS","UPCCAS","LNS","CPSolver"],"",StrDomain(None,deactivable=True),"Add a comment to the results sample.")# help
        self.declareOption("f","folder",["timeDependentSolver","BPCCAS","UPCCAS","LNS","CPSolver","dataVisu"],None,IntDomain(left_bound=0,deactivable=True),"Indicates the data folder number.")# help
        self.declareOption("n","file",["timeDependentSolver","BPCCAS","UPCCAS","LNS","CPSolver","dataVisu"],None,IntDomain(left_bound=0,deactivable=True),"Indicates the data instance number.")# help
        self.declareOption("s","seed",["BPCCAS","UPCCAS","LNS","CPSolver","dataVisu"],0,IntDomain(),description="Random seed.")
        self.declareOption(None,"include_systematic",["BPCCAS","UPCCAS","LNS","CPSolver","dataVisu"],False,BooleanDomain(),description="Activate to include systematic requests.")
        self.declareOption(None,"restart",["LNS","BPCCAS","UPCCAS"],1,IntDomain(left_bound=0),description="Number of restart of the external TSPTW solvers (when compatible).")
        self.declareOption(None,"dynamic",["LNS","BPCCAS","UPCCAS","CPSolver"],False,BooleanDomain(),description="Activate to set the dynamic arrival of the requests.")
        self.declareOption(None,"test_req",["LNS","BPCCAS","UPCCAS","CPSolver"],np.Inf,IntDomain(left_bound=1),description="Limit the number of requests considered (meant for debug).")
        self.declareOption(None,"profile",["LNS","BPCCAS","UPCCAS","CPSolver"],False,BooleanDomain(),description="Activate the profiliing of the transition model (meant for debug).")
        
        # limiter les options aux algos indiqués
        self.applyOptionsParsingRules()
    
    def isTimeDependentModeOn(self):
        return self.getOptValue("initTransitionModel")=="time-dep" or self.getOptValue("solver")=="timeDependentSolver"
    
    def isDataTimeDependent(self):
        return self.getOptValue("data")=="time-dependent"
    
    def displayHelp(self):
        for opt in self.options:
            if self.getOptValue("solver") in self.options[opt].getAlgos():
                if self.options[opt].getShortName() is not None:
                    print("-"+self.options[opt].getShortName(),end=" || ")
                print("--"+ self.options[opt].getName()+": "+self.options[opt].getDescription())
        
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
            if o=="--cplex" and a is not None:
                self.glob.docplexPath = a
                
        if self.donnees.algo_name is None:
            self.donnees.algo_name = self.getOptValue("solver")
            
        #print(opts,args)
        op = [x[0] for x in opts]
        if '--destroy' in op and '--noise' in op:
            raise Exception('--destroy et --noise options can\'t be used simultaneously')
        for o, a in opts:
            self.restrict(o,self.restrict_options)
            opt = self.findOpt(o)
            self.options[opt].setValue(a)
        #self.options["time"].setValue(self.options["time"].getValue()*60)               
        # traiter les incohérences d'options avec le solver utilisé
        if not(self.getOptValue("solver") in self.algos):
            raise NameError("Unknown algorithm",self.donnees.algo_name)
        
        instance = {"folder":self.getOptValue('folder'),"file":self.getOptValue('file')}
        if instance["folder"] is None:
            instance = None
        self.donnees.instance = instance
        if instance is not None and 'file' in instance and 'folder' not in instance:
            raise ValueError("Specify the data folder")
        
        if self.getOptValue("solver")=="timeDependentSolver":
            if not self.getOptValue("data")=="time-dependent":
                raise ValueError("Use time-dependent instances with the solver 'timeDependentSolver'. --data=time-dependent")
        
        if self.getOptValue("chain"):
            if not self.getOptValue("initTransitionModel") in ["slow","mean","fast"]:
                raise ValueError("Define a new approximated model using --initTransitionModel=<model>. Possible values : "+str(["slow","mean","fast"]))
            if self.getOptValue("initDur")>=self.getOptValue("time"):
                raise ValueError("Time allocated to approximated solve exceeding the overall time limit. Specify --initDur=<duree>.")
        return instance
    
    def __str__(self):
        mess = " ------------ configuration globale ---------------------\n"
        for opt in self.options:
            if self.donnees.algo_name in self.options[opt].getAlgos():
                mess += " | " + opt + ' : ' + str(self.options[opt].getValue()) + '\n'
        globDict = self.glob.toDict()
        for key,value in globDict.items():
            mess +=" | " + key + ' : ' + str(value) + '\n'
       
        for key,value in self.donnees.toDict().items():
            mess += " | " + key + " : " + str(value) + '\n'
        if self.donnees.algo_name == "CPSolver":
            mess += " ------------ CP configuration ---------------\n"
            dictGlobalSol = self.CPSolver.toDict()
            for key,value in dictGlobalSol.items():
                mess +=" | " + key + ' : ' + str(value) + '\n'
        if self.donnees.algo_name == "LNS":
            dictLNS = self.LNS.toDict()
            mess += " ------------ LNS  configuration -------------------------\n"
            for key,value in dictLNS.items():
                mess +=" | " + key + ' : ' + str(value) + '\n'
        if self.donnees.algo_name == "BPCCAS":
            dictBatch = self.batchProg.toDict()
            mess += " ------------ BPCCAS configuration ----------------------\n"
            for key,value in dictBatch.items():
                mess +=" | " + key + ' : ' + str(value) + '\n'
        if self.donnees.algo_name == 'UPCCAS':
            mess += " ------------ UPCCAS configuration -----------------------\n"
            for key,value in self.UPCCAS.toDict().items():
                mess +=" | " + key + ' : ' + str(value) + '\n'            
        if self.donnees.algo_name in ['BPCCAS','UPPCAS'] or self.donnees.algo_name == 'LNS':
            if self.glob.solver == 'OPTW':
                dictOPTW = self.OPTW.toDict()
                mess += " ------------ OPTW configuration --------------------------\n"
                for key,value in dictOPTW.items():
                    mess +=" | " + key + ' : ' + str(value) + '\n'
            else:
                dictLKH = self.lkh.toDict()
                mess += " ------------ LKH configuration --------------------------\n"
                for key,value in dictLKH.items():
                    mess +=" | " + key + ' : ' + str(value) + '\n'
        mess += " ------------ transition model --------------------------\n"
        if self.getOptValue("initTransitionModel")=="fast":
            transitionModel = FastModel
        elif self.getOptValue("initTransitionModel")=="mean":
            transitionModel = MeanModel
        elif self.getOptValue("initTransitionModel")=="slow":
            transitionModel = SlowModel
        else:
            transitionModel = "time-dependent model"
        mess += " | initial model : " + str(transitionModel) +"\n"
        mess += " ------------ OPTIONS -----------------------\n"
        for opt in self.options:
            mess += " | " + str(opt) + " : " + str(self.getOptValue(opt)) + "\n"
        return mess
            
        
        
        
        
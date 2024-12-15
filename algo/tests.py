from os import walk,getpid
from datetime import datetime
import subprocess
import time
import pickle


sample_path = "tests/process_"+str(getpid())+"_"+str(datetime.today()).replace(" ","_")
BPCCAS = "BPCCAS.py"
UPCCAS = "UPCCAS.py"
LNS = "LNS.py"

preval = " --preval "
conserve = " --conserve "
stable = " --stable_exp "
verif = " --verif "
no_option = ""
greedy_ordering = " --ordering=greedy "
uniqueLNS = " --max_iteration=1 "

print("Tests : "+sample_path)

class Test:
    def __init__(self,instance,instance_name,algo,options,nom_test,modes_attendus,obs_attendues,reward_attendu,ncoeurs):
        self.instance = instance
        self.instance_name = instance_name
        self.algo = algo
        self.nom_test = nom_test
        self.obs_attendues = obs_attendues
        self.reward_attendu = reward_attendu
        self.options = options
        self.modes_attendus = modes_attendus
        self.mpiexec_cmd = [LNS,BPCCAS,UPCCAS]
        self.python_cmd = []
        self.ncoeurs = ncoeurs
        
    def buildCmdLine(self):
        cmd_line = ""
        if self.algo in self.mpiexec_cmd:
            cmd_line += "mpiexec -n " + str(self.ncoeurs) + " python3 "
        else:
            assert(self.algo in self.python_cmd)
            cmd_line += "python3 "
        cmd_line += self.algo + " "
        cmd_line += self.options
        cmd_line += " -n " + str(self.instance) + " "
        if self.algo in self.python_cmd:
            cmd_line += " -w " + str(self.ncoeurs)
        cmd_line += " --sample=" + sample_path
        cmd_line += " -m \"" + self.nom_test +"\""
        return cmd_line
    
    def runAlgo(self):
        cmd_line = self.buildCmdLine()
        #print(cmd_line)
        p = subprocess.Popen(cmd_line, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        #for line in p.stdout.readlines():
        #    print(line)
        retval = p.wait()

    def compareList(self,liste1,liste2):
        l1_notl2 = []
        l2_notl1 = []
        for x in liste1:
            if x not in liste2:
                l1_notl2.append(x)
        for x in liste2:
            if x not in liste1:
                l2_notl1.append(x)
        return l1_notl2,l2_notl1
    
    
    def extraireBool(self,dico,key):
        lines = dico['infos hist'].split('\n')
        line = [x for x in lines if "| "+key+" : " in x][0]
        line = line.split('|')[1].split(":")[1]
        return eval(line)
        
    def testStabilite(self,dico):
        return (sum([x in stable for x in self.options.split()]))==self.extraireBool(dico,'stabilités des explications')
                
    def testPreval(self,dico):
        return (sum([x in preval for x in self.options.split()]))==self.extraireBool(dico,'prevalidation')
    
    def testConserve(self,dico):
        return (sum([x in conserve for x in self.options.split()])) == self.extraireBool(dico,'conservation des modes')
                
    def compareFile(self,dico):
        algo = self.algo.split(".py")[0]
        cmp = dico['instance']==self.instance_name and dico['algo']==algo and self.ncoeurs==int(dico['coeurs'])
        if cmp:
            if algo=='UPCCAS':
                if not self.testPreval(dico):
                    return False
            elif algo=='BPCCAS':
                if not self.testConserve(dico):
                    return False
                if not self.testStabilite(dico):
                    return False
            return True
        else:
            return False
    
    def searchFile(self):
        found = False
        for root,dirs,files in walk("../results/samples/"+sample_path):
            for file in files:
                with open("../results/samples/"+sample_path+"/"+file,'rb') as handle:
                    dico_hist  = pickle.load(handle)
                    if self.compareFile(dico_hist):
                        return dico_hist
        raise Exception('fichier non trouvé')
    
    def extraireObs(self,dico):
        obs = []
        for s in dico['sol']:
            obs += dico['sol'][s]
        return obs
        
    def verifierResultats(self):
        print("======== TEST " + self.nom_test+" ============")
        dico_hist  = self.searchFile()
        if round(self.reward_attendu,4) != round(dico_hist['obj'][0],4):
            print("- Reward obtenu :",dico_hist['obj'][0],"/ reward attendu :",self.reward_attendu)
            #print(dico_hist['infos hist'])
        else:
            print("- Reward OK")
        if self.modes_attendus != []:
            x1,x2 = self.compareList(self.modes_attendus,dico_hist['modes'])
            if len(x1)>0 or len(x2) > 0:
                print("- Modes attendus :",self.modes_attendus)
                print("  Modes obtenus :",dico_hist['modes'])
                print("  Modes manquants :",x1)
                print("  Modes en trop :",x2)
            else:
                print("- Modes OK")
            x1,x2 = self.compareList(self.obs_attendues,self.extraireObs(dico_hist))
            if len(x1)>0 or len(x2) > 0:
                print("- Activités attendus :",self.obs_attendues)
                print("  Activités obtenus :",self.extraireObs(dico_hist))
                print("  Activités manquantes :",x1)
                print("  Activités en trop :",x2)
            else:
                print("- Activités OK")               

tests = []




"""
Scoring des obs (classes d'équivalences)
    0 < x <= 0.5 : 0.05
    0.5 < x <= 0.6 : 0.1
    0.6 < x <= 0.7 : 0.4
    0.7 < x <= 0.8 : 0.7
    0.9 < x <= 0.9 : 0.95
    0.95 < x <= 1 : 1
"""

"""
instance 101 : durée obs = 9
    r0 mono
        w0 = 1 [0,10] sat 0
    r1 systematic
        w1 = 0.05 [0,10] sat 0 
        w2 = 0.05 [200,210] sat 0
        w3 = 0.05 [200,210] sat 1
    r2 mono
        w4 = 0.1 [200,210] sat 0
        w5 = 0.05 [0,10] sat 1 

"""
ncoeurs = 10

instance = "test_rejet_premature"
# instance,instance_name,algo,options,nom_test,modes_attendus,obs_attendues,reward_attendu,ncoeurs
"""
    BPCCAS classique : génération optimiste
"""
modes = [(-1,0),(0,0),(1,1),(2,1)]
obs = [0,2,5,6,7] # 6,7 = vidages
tests.append(Test(101,instance,BPCCAS,no_option,"test_rejet_premature BPCCAS",modes,obs,1.1,ncoeurs))


"""
    BPCCAS stable : génération pessimiste (CCR)
"""
modes = [(-1,0),(0,0),(2,0),(1,2)]
obs = [0,3,4,6,7] # 6,7 = vidages
tests.append(Test(101,instance,BPCCAS,stable,"test_rejet_premature BPCCAS stable",modes,obs,1.15,ncoeurs))


"""
    BPCCAS conserve : version qui conserve l'intégralité des modes (ordre inchangé sur les cca)
"""
modes = [(-1,0),(0,0),(1,1),(2,1)]
obs = [0,2,5,6,7] # 6,7 = vidages
tests.append(Test(101,instance,BPCCAS,conserve,"test_rejet_premature BPCCAS conserve",modes,obs,1.1,ncoeurs))


"""
    UPCCAS classique : pas de prévalidation
"""
modes = [(-1,0),(0,0),(1,2),(2,0)]
obs = [0,3,4,6,7] # 6,7 = vidages
tests.append(Test(101,instance,UPCCAS,no_option,"test_rejet_premature UPCCAS",modes,obs,1.15,ncoeurs))

"""
    UPCCAS avec prévalidation :
    Avec la prévalidation on obtient différents résultats selon le nombre de coeurs : 1.15 avec 10 coeurs, 1.1 avec 2 coeurs

"""
modes = [(-1,0),(0,0),(1,1),(2,1)]
obs = [0,2,5,6,7] # 6,7 = vidages
tests.append(Test(101,"test_rejet_premature",UPCCAS,preval,"test_rejet_premature UPCCAS preval 10 coeurs",modes,obs,1.1,10))

modes = [(-1,0),(0,0),(1,2),(2,0)]
obs = [0,3,4,6,7] # 6,7 = vidages
tests.append(Test(101,"test_rejet_premature",UPCCAS,preval,"test_rejet_premature UPCCAS preval 2 coeurs",modes,obs,1.15,2))

"""
    LNS : la méthode de génération n'étant pas la même on a pas de déterminisme sur le chemin
    pris en cas d'égalité lors de la génération de mode. 
    Pour la requête 2 on pourrait soit :
        choisir o2 sur le satellite 0
        choisir o3 sur le satellite 1
    ici l'algo choisit o3 (alors que la génération de BPCCAS choisit o2)
    la requête 3 est alors impactée : choix de o4 (au lieu de o5)
"""
modes = [(-1,0),(0,0),(1,0),(2,0)]
obs = [0,4,3,6,7] # 6,7 = vidages
tests.append(Test(101,instance,LNS,greedy_ordering+uniqueLNS,"test LNS 1ere iteration",modes,obs,1.15,1))


"""
Scoring des obs (classes d'équivalences)
    0 < x <= 0.5 : 0.05
    0.5 < x <= 0.6 : 0.1
    0.6 < x <= 0.7 : 0.4
    0.7 < x <= 0.8 : 0.7
    0.9 < x <= 0.9 : 0.95
    0.95 < x <= 1 : 1
"""

"""
instance 102 : durée obs = 9
    r0 periodique
        w0 = 0.05 [0,10] sat 0 TS0
        w1 = 0.05 [200,210] sat 0 TS1
        w2 = 0.05 [400,410] sat 0 TS2
        w3 = 0.4 [600,610] sat 0 TS3
        w4 = 0.05 [610,620] sat 0 TS3
    r2 mono
        w5 = 1 [600,610] sat 0 
    r3 periodique
        w6 = 0.05 [-210,-200] sat 0
        w7 = 0.05 [0,10] sat 0
        w8 = 0.1 [200,210] sat 0
        w9 = 0.05 [400,410] sat 0 
    d0
"""

instance = 'test_prevalidation'
obs = [0,1,2,4,5,10]
modes = [(0,1),(2,0),(-1,0)]
tests.append(Test(102,instance,BPCCAS,no_option,"test prévalidation BPCCAS",modes,obs,1.2,ncoeurs))

obs = [6,7,8,9,5,10]
modes = [(3,0),(2,0),(-1,0)]
tests.append(Test(102,instance,BPCCAS,stable,"test prévalidation BPCCAS stable",modes,obs,1.25,ncoeurs))

obs = [6,7,8,9,5,10]
modes = [(3,0),(2,0),(-1,0)]
tests.append(Test(102,instance,UPCCAS,no_option,"test prévalidation UPCCAS",modes,obs,1.25,ncoeurs))

obs = [0,1,2,4,5,10]
modes = [(0,1),(2,0),(-1,0)]
tests.append(Test(102,instance,UPCCAS,preval,"test prévalidation UPCCAS preval",modes,obs,1.2,ncoeurs))

# mode 3 dégradé : 1 seule période ...
obs = [0,1,2,4,5,6,10]
modes = [(0,0),(2,0),(3,0),(-1,0)]
tests.append(Test(102,instance,LNS,greedy_ordering+uniqueLNS,"test LNS 1ere iteration",modes,obs,1.25,1))


"""
Scoring des obs (classes d'équivalences)
    0 < x <= 0.5 : 0.05
    0.5 < x <= 0.6 : 0.1
    0.6 < x <= 0.7 : 0.4
    0.7 < x <= 0.8 : 0.7
    0.9 < x <= 0.9 : 0.95
    0.95 < x <= 1 : 1
"""

"""
instance 103 : durée obs = 9
    r0 mono
        w0 = 1 [0,10] sat 0
    r1 mono
        w1 = 0.95 [0,10] sat 0     
        w2 = 0.7 [200,210] sat 0 
    r2 mono
        w3 = 0.05 [200,210] sat 0
    d0

"""
instance = "test_conservation_batch"
obs = [0,2,4]
modes = [(-1,0),(1,1),(0,0)]
tests.append(Test(103,instance,BPCCAS,no_option,"test conservation batch BPCCAS",modes,obs,1.7,ncoeurs))

obs = [0,3,4]
modes = [(-1,0),(2,0),(0,0)]
tests.append(Test(103,instance,BPCCAS,conserve,"test conservation batch BPCCAS conserve",modes,obs,1.05,ncoeurs))

obs = [0,2,4]
modes = [(-1,0),(1,1),(0,0)]
tests.append(Test(103,instance,UPCCAS,no_option,"test conservation batch UPCCAS",modes,obs,1.7,ncoeurs))

obs = [0,2,4]
modes = [(-1,0),(1,1),(0,0)]
tests.append(Test(103,instance,UPCCAS,preval,"test conservation batch UPCCAS preval",modes,obs,1.7,ncoeurs))

obs = [0,2,4]
modes = [(-1,0),(1,0),(0,0)]
tests.append(Test(103,instance,LNS,greedy_ordering+uniqueLNS,"test LNS 1ere iteration",modes,obs,1.7,1))

for test in tests:
    test.runAlgo()

for test in tests:
    test.verifierResultats()

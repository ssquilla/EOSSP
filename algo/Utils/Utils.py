"""
    Quelques fonctions standards utiles.

"""
from .config import *
global config
config = Config()

import colorama
from colorama import Fore
from colorama import Style
from mpi4py import MPI
import random as rd
from time import time
import bisect
from copy import deepcopy 
import math


def genericMax(a,b):
    if a >= b:
        return a
    else:
        return b
    
def genericMaxList(list_generic):
    if len(list_generic)==0:
        raise Exception('list vide : pas de max')
    max_val = list_generic[0]
    for elmt in list_generic:
        if elmt>max_val:
        #if max_val<elmt:
            max_val=elmt
    return max_val
    

def step():
    if config.getOptValue("step"):
        printColor("Appuyez sur n'importe quelle touche pour continuer.",c='y',align_left=True,force=True)
        input()

def arrondiPlusProche(x,digits):
    shift = x*(10**digits)
    integer = int(math.floor(shift))
    if (x-integer)<(integer+1-x):
        return integer/10**digits
    else:
        return (integer+1)/10**digits

def arrondi(x,digits,up=True):
    if up:
        return math.ceil(x*(10**digits))/10**digits
    else:
        return math.floor(x*(10**digits))/10**digits
        
def safe_list_get (l, idx, default):
    try:
        return l[idx]
    except IndexError:
        return default

def assertLessThan(a,b,msg=""):
    if a>b:
        die(a,'<=',b,'failed.',msg)
    
def shuffleSort(liste,key,reverse):
    cle = lambda x : (key(x),rd.random())
    return sorted(liste,key=cle,reverse=reverse)

def reverse_insort(a, x, lo=0, hi=None):
    """Insert item x in list a, and keep it reverse-sorted assuming a
    is reverse-sorted.

    If x is already in a, insert it to the right of the rightmost x.

    Optional args lo (default 0) and hi (default len(a)) bound the
    slice of a to be searched.
    """
    if lo < 0:
        raise ValueError('lo must be non-negative')
    if hi is None:
        hi = len(a)
    while lo < hi:
        mid = (lo+hi)//2
        if x > a[mid]: hi = mid
        else: lo = mid+1
    a.insert(lo, x)

def indexSortedList(a, x):
    'Locate the leftmost value exactly equal to x'
    i = bisect.bisect_left(a, x)
    if i != len(a) and a[i] == x:
        return i
    raise ValueError

def isInSortedList(a,x):
    try:
        indexSortedList(a,x)
        return True
    except ValueError:
        return False
    
def delElement(liste,test):
    suppr = []
    for i,x in enumerate(liste):
        if test(x):
            suppr.append(i)
    suppr.reverse()
    for i in suppr:
        liste.pop(i)

def legend_without_duplicate_labels(ax):
    handles, labels = ax.get_legend_handles_labels()
    unique = [(h, l) for i, (h, l) in enumerate(zip(handles, labels)) if l not in labels[:i]]
    ax.legend(*zip(*unique))

def intersect(I,J):
    a=max(I[0],J[0])
    b=min(I[1],J[1])
    return a<=b
    
def stepIntegral(valeurs,bornes=None):
    integral = 0
    for i,(t,y) in enumerate(valeurs):
        if(i<len(valeurs)-1):
            suivant=valeurs[i+1]
            if bornes is None or t>=bornes[0]:
                gauche = t
            else:
                gauche = 0
                droite = 0
            if bornes is None or t<=suivant[0]:
                droite = suivant[0]
            else:
                droite = bornes[1]
            aire = y*(droite-gauche)
            integral += aire
    return integral

def alert(*args,**kwargs):
    kwargs['align_left'] = True
    kwargs['force'] = True
    printColor(*args,**kwargs,c='r')

def warn(*args,**kwargs):
    prof = getDisplayDepth()
    shiftLeftDisplay(prof)
    kwargs['align_left'] = True
    kwargs['force'] = True
    printColor(*args,**kwargs,c='y')
    shiftRightDisplay(prof)
    
    
def histToLatex(n,bins,xlabel,ylabel):
    chaine = "\\begin{tikzpicture}\n"
    chaine += "\\begin{axis}[area style,xlabel={"+xlabel+"},ylabel={"+ylabel+"}]\n"
    chaine += "\\addplot+[ybar interval,mark=no] plot coordinates { \n"
    for i in range(len(n)):
        chaine += "\t\t\t("+str(bins[i])+str(",")+str(n[i])+")\n"

    chaine += "};\n"
    chaine += "\\end{axis}\n"
    chaine += "\\end{tikzpicture}\n"
    return chaine

def insererListeTriee(liste,e,f=lambda e:e):
    if(len(liste)==0):
        return [e]
    elif(f(e) <= f(liste[0])):
        return [e] + liste
    
    for i,x in enumerate(liste):
        if(i<len(liste)-1 and f(x) <= f(e) and f(e) <= liste[i+1]):
            return liste[:i+1] + [e] + liste[i+1:]
    else:
        return liste + [e]
    
global profondeur_print
profondeur_print = 1
global timers
timers = []
global colors
colors = []

def printMaster(*args,**kwargs):
    if MPI.COMM_WORLD.Get_rank()==0:
        printColor(*args,**kwargs)

def printCore(*args,**kwargs):
    core = kwargs.get('core',0)
    del kwargs['core']
    if MPI.COMM_WORLD.Get_rank()==core:
        printColor(*args,**kwargs)


def getDisplayDepth():
    global profondeur_print
    return profondeur_print

def castColor(text,**kwargs):
    if 'c' in kwargs:
            c = kwargs['c']
            del kwargs['c']
    else:
            c = 'w'
    if c=='w':
            couleur = Style.RESET_ALL
    elif c=='r':
            couleur = Fore.RED
    elif c=='g':
            couleur = Fore.GREEN
    elif c=='y':
            couleur = Fore.YELLOW
    elif c=='b':
            couleur = Fore.BLUE
    elif c=='m':
            couleur = Fore.MAGENTA
    elif c=='c':
            couleur = Fore.CYAN
    else:
            raise ValueError("couleur inconnue")
    mess = couleur + text + Style.RESET_ALL
    return mess


def printColor(*args,**kwargs):
    core = kwargs.get('core',None)
    if core is None or MPI.COMM_WORLD.Get_rank()==core:
        if "align_left" in kwargs:
            align_left = kwargs['align_left']
        else:
            align_left = False
        prof = None
        if 'depth' in kwargs:
            prof = kwargs['depth']
            assert(type(prof)==int)
            assert(prof>=0)
        else:
            prof = profondeur_print
        if not align_left:
            if config.getOptValue("verbose")>=prof or kwargs.get('force',False):
                if 'align_left' in kwargs:
                    del kwargs['align_left']
                print((prof-1)*'\t',end='')
                printNoSpace(*args,**kwargs)
        else:
            printNoSpace(*args,**kwargs)

def remove(param_dict,supp):
    for s in supp:
        if s in param_dict:
            del param_dict[s]
        
def printNoSpace(*args,**kwargs):
    core = kwargs.get('core',None)
    if core is None or MPI.COMM_WORLD.Get_rank()==core:
        if 'depth' in kwargs:
            prof = kwargs['depth']
            assert(type(prof)==int)
            assert(prof>=0)
            del kwargs['depth']
        else:
            prof = profondeur_print
        if kwargs.get('force',False) or config.getOptValue("verbose")>=prof:
            mess = ""
            for i,arg in enumerate(args):
                    replace = str(arg)
                    replace = replace.replace("\n","\n"+(prof-1)*'\t')
                    mess += " " + str(replace)
            
            if 'time' in kwargs:
                mess += " time : " +str(kwargs['time']) +"s"
                del kwargs['time']
            mess = castColor(mess,**kwargs)
            mess += Style.RESET_ALL
            remove(kwargs,['c','force','align_left','core'])
            print(mess,**kwargs)

def printOpen(*args,**kwargs):
    printColor(*args,'{',**kwargs)
    global profondeur_print 
    profondeur_print += 1
    timers.append(time())
    global colors
    if 'c' in kwargs:
        colors.append(kwargs['c'])
    else:
        colors.append('w')

def shiftRightDisplay(n_shift):
    global profondeur_print
    profondeur_print += n_shift
        
def shiftLeftDisplay(n_shift):
    global profondeur_print
    profondeur_print -= n_shift
    
def printClose(*args,**kwargs):
    global profondeur_print
    profondeur_print -= 1
    delay = round(time()-timers.pop(-1),7)
    global colors
    # message de cloture
    mess = ""
    for i,arg in enumerate(args):
        replace = str(arg)
        replace = replace.replace("\n","\n"+(profondeur_print-1)*'\t')
        mess += " " + str(replace)
    mess = castColor(mess,**kwargs)
    
    couleur = colors.pop(-1)
    kwargs['time'] = delay
    kwargs['c'] = couleur
    
    printColor(mess,end='')
    printNoSpace('}',**kwargs)    
      
def die(*args,**kwargs):
    global profondeur_print
    shift = profondeur_print
    shiftLeftDisplay(shift)
    printColor("======== [DEBUG KILL] ==========",c='r')
    printColor('|',c='r',end='')
    printColor(*args,**kwargs,c='r')
    printColor("================================",c='r')
    shiftRightDisplay(shift)
    raise ValueError()
    
    
def estVidage(a):
    return a>=config.donnees.no

def stringSequenceActivite(sequence):
    res = "["
    for a in sequence:
        if estVidage(a):
            res += Fore.GREEN + str(a) + Style.RESET_ALL + " "
        else:
            res += Fore.BLUE + str(a) + Style.RESET_ALL + " "
    res += "]"
    return res

def printSequenceActivite(sequence):
    print(stringSequenceActivite(sequence))

def choseFolder(root,instance):
    dir_list = sorted(os.listdir(root))#[x for x in os.listdir(root) if os.path.isdir(x)])
    while True:
        try:
            profondeur = len(root.split('/'))
            if instance is None:
                printColor(2*max(0,profondeur-1)*' '+root,c='y')
            for i,name in enumerate(dir_list):
                if os.path.isfile(root+'/'+name):
                    couleur = 'w'
                else:
                    couleur = 'y'
                if instance is None:
                    printColor(str(i)+")",(2*profondeur-len(str(i)+")"))*' '+name,c=couleur)
            if instance is None:
                printColor("votre choix : ",end='',c='y')
                res = input()
                x = int(res)
            else:
                x = int(instance['folder'])
            return dir_list[x]
        except Exception as e:
            print(e)
            print("Sélection invalide.")

def cleDossier(nom_dossier):
    try:
        # format legerement different dans le cas time-dependent
        if config.isDataTimeDependent():
            suffixe = nom_dossier.split("/")[-1]
            info = suffixe.split("_")
            req,meteo = [int(x) for x in info[0].split("-")],info[2]
            seed = int(info[1].split("-")[1])
            sat = np.prod(np.array([int(x) for x in info[3].split("x")]))
            POI = eval(info[4].split("-")[1])
            if POI is None:
                POI = np.Inf
            nombre_requetes = sum(req)
            repartition = tuple(req)
            return (meteo,sat,POI,nombre_requetes,repartition,seed)
        else:
            suffixe = nom_dossier.split("/")[-1]
            info = suffixe.split("_")
            req,meteo = [int(x) for x in info[0].split("-")],info[1]
            sat = np.prod(np.array([int(x) for x in info[2].split("x")]))
            POI = eval(info[3].split("-")[1])
            if POI is None:
                POI = np.Inf
            nombre_requetes = sum(req)
            repartition = tuple(req)
            return (meteo,sat,POI,nombre_requetes,repartition)
    except Exception as e:
        print(e)
        if config.isDataTimeDependent():
            return ("",0,np.Inf,0,(0,0,0,0,0),0)
        else:
            return ("",0,np.Inf,0,(0,0,0,0,0))

def getSeed(file_name):
    return int(file_name.split("seed-")[1].split("_")[0])
    
def recursiveChoice(root_choice,instance):
    fileList = []
    i = 0
    for root, dirs, files in os.walk(root_choice):
        profondeur = len(root.split('/'))
        if instance is None or 'file' not in instance:
            printColor(2*max(0,profondeur-1)*' '+root,c='y')
        for name in sorted(files,key=getSeed):
            if instance is None or 'file' not in instance:
                printColor(str(i)+")",(2*profondeur-len(str(i)+")"))*' '+name)
            fileList.append(root+'/'+name)
            i += 1

    if instance is None or 'file' not in instance:
        printColor("votre choix : ",end='',c='y')
        x = int(input())
        
    else:
        x = instance['file']
    return fileList[x]

def sortedRecursiveChoice(root_choice,instance):
    fileList = []
    i = 0
    dir_enum = {}
    for root, dirs, files in os.walk(root_choice):
        if len(files)>0:
            dir_enum[root] = files
    for root in sorted(dir_enum,key=cleDossier):
        profondeur = len(root.split('/'))
        printColor(2*max(0,profondeur-1)*' '+root,c='y')
        for name in sorted(dir_enum[root],key=getSeed):
            if instance is None or 'file' not in instance:
                printColor(str(i)+")",(2*profondeur-len(str(i)+")"))*' '+name)
            fileList.append(root+'/'+name)
            i += 1

    if instance is None or 'file' not in instance:
        printColor("votre choix : ",end='',c='y')
        x = int(input())
        
    else:
        x = instance['file']
    return fileList[x]

def organizedChoice(root_choice,instance):
    def extractParents(profStart,file_name):
        enumerating_parents = []
        listing = file_name.split("/")
        start_root = ""
        for elemt in listing[:profStart+1]:
            start_root += elemt + "/"
        enumerating_parents.append(deepcopy(start_root))
        for elmt in listing[profStart+1:]:
            start_root += elmt + "/"
            enumerating_parents.append(deepcopy(start_root))
        return enumerating_parents
    def getProfondeur(name):
        return len(name.split('/'))    
    def checkFormat(name):
        if name[0]=='.':
            return False
        if not os.path.isfile(name) and "components" not in name.split("/")[-1]:
            return True
        if len(name.split(".pb"))>1:
            return True
        return False
    fileList = []
    i = 0
    dir_enum = {}
    for root, dirs, files in os.walk(root_choice):
        if "components" not in root.split("/")[-1] and sum(["components" in f for f in files])==0:
            if len(files)>0:
                dir_enum[root] = files
    parents_explored = []
    for root in sorted(dir_enum,key=cleDossier):
        names = [name for name in dir_enum[root] if checkFormat(name)]
        for name in sorted(names,key=getSeed):
            # possible potential new folder
            listing_parents = extractParents(2,root)
            for parent in listing_parents:
                if parent not in parents_explored:
                    parents_explored.append(parent)
                    profondeur_parent = getProfondeur(parent)
                    if instance is None or 'file' not in instance:
                        dossier_instance = "components_time_dep" in os.listdir(parent)
                        # l'affichage en mode "time-dep" est shifte d'un cran a gauche.
                        # => sinon affichage d'un nom de dossier supplementaire et inutile
                        # (du au format des donnees different du format time-independent)
                        if not config.isDataTimeDependent() or not dossier_instance:
                            printColor(2*max(0,profondeur_parent-1)*' '+parent,c='y')
            # print the file name and append to the list
            if instance is None or 'file' not in instance:
                profondeur = getProfondeur(root)+1-(config.getOptValue("data")=="time-dependent")
                printColor(str(i)+")",(2*profondeur-len(str(i)+")"))*' '+name,end="")
            if "components_time_dep" in os.listdir(root):
                printColor(" (composantes précalculées)",c='g')
            else:
                printColor()
            fileList.append(root+'/'+name)
            i += 1
    if instance is None or 'file' not in instance:
        printColor("votre choix : ",end='',c='y')
        x = int(input())
    else:
        x = instance['file']
    file = fileList[x]
    return file

          
def choseAndBroadcastFile(path,instance):
    shift = getDisplayDepth()-1
    shiftLeftDisplay(shift)
    comm = MPI.COMM_WORLD
    x = None
    while True:
        try:
            # Sélection à la racine du 1er dossier (potentiel fichier tout de même)
            root_choice = path+'/'+choseFolder(path,instance)
            # choix récursif de n'importe quel fichier dans l'arborescence
            if os.path.isfile(root_choice):
                file = root_choice
                if instance is not None and 'file' in instance:
                    raise ValueError(str(file) +' n\'est pas un dossier : mauvaise utilisation de -f et -n')
            else:
                    file = organizedChoice(root_choice,instance)
                    #file = recursiveChoice(root_choice,instance)
                    printColor(file,c='w')
            data = {'file':file}
            comm.bcast(data,root=0)
            break
        except ValueError as e:
            print(e)
            printColor("Erreur input",c='r')
    shiftRightDisplay(shift)
    return data['file']

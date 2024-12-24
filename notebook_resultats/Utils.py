"""
    Quelques fonctions standards utiles.

"""
import colorama
from colorama import Fore
from colorama import Style
from mpi4py import MPI
import random as rd
from time import time
import bisect
from copy import deepcopy 

def safe_list_get (l, idx, default):
    try:
        return l[idx]
    except IndexError:
        return default
        
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
    if "align_left" in kwargs:
        align_left = kwargs['align_left']
    else:
        align_left = False
    if not align_left:
        if True:
            if 'force' in kwargs:
                del kwargs['force']
            if 'align_left' in kwargs:
                del kwargs['align_left']
            print((profondeur_print-1)*'\t',end='')
            printNoSpace(*args,**kwargs)
    else:
        printNoSpace(*args,**kwargs)
        
def printNoSpace(*args,**kwargs):
    if kwargs.get('force',False) or True:
        mess = ""
        for i,arg in enumerate(args):
                replace = str(arg)
                replace = replace.replace("\n","\n"+(profondeur_print-1)*'\t')
                mess += " " + str(replace)
        
        if 'time' in kwargs:
            mess += " time : " +str(kwargs['time']) +"s"
            del kwargs['time']
        mess = castColor(mess,**kwargs)
        mess += Style.RESET_ALL
        if 'c' in kwargs:
            del kwargs['c']
        if 'force' in kwargs:
            del kwargs['force']
        if 'align_left' in kwargs:
            del kwargs['align_left']
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


def cleDossier(nom_dossier):
    try:
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
        return ("",0,np.Inf,0,(0,0,0,0,0))
    
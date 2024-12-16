from .components import *

from ..Utils.config import *
config = Config()

from mpi4py import MPI
import socket
from time import *
import math
from os.path import exists

class InfeasibleSolutionfromSolverException(Exception):
    pass

def supprimerListeElements(liste,elmts):
    del_idx = []
    for i,a in enumerate(liste):
        if a in elmts:
            del_idx.append(i)
    for j in range(len(del_idx)-1,-1,-1):
        i = del_idx[j]
        liste.pop(i)
 
class NumericalError(Exception):
    pass

""" fonctions pour trouver la maneuvre qui finit a une certaine date (time-dep) """           
# trouve la sate sur le segment telle que la transition arrive sur la date tnext
def inverseLinearManeuver(tnext,points_left,points_right,duree):
    tleft,L = points_left
    tright,R = points_right
    amplitude = tright-tleft
    duration_delta = R-L
    t = ((tnext)*amplitude+tleft*R-tright*L)/(duration_delta+amplitude)
    #t = ((tnext-duree)*amplitude+tleft*R-tright*L)/(duration_delta+amplitude)
    return t
        
# left = right si points_controle[left] vaut exactement la date recherchee
# left + 1 = right sinon. Alors la bascule est dans ce segment
def searchChangeSegment(points_controle,date_next_a,duree):
    left = 0
    right = len(points_controle)-1
        
    start_maneuver,maneuver = points_controle[left]
    #while(left<len(points_controle)-1 and start_maneuver+maneuver+duree<date_next_a):
    while(left<len(points_controle)-1 and start_maneuver+maneuver<date_next_a):
        left += 1
        start_maneuver,maneuver = points_controle[left]
    left -= 1    
        
    start_maneuver,maneuver = points_controle[right]
    #while(right > 0 and start_maneuver+maneuver+duree>date_next_a):
    while(right > 0 and start_maneuver+maneuver>date_next_a):
        right -= 1        
        start_maneuver,maneuver = points_controle[right]
    right += 1
    return left,right    

class InfeasibleSolutionException(Exception):
    def __init__(self,msg=""):
        super().__init__(msg)

class OutdatedCriticalPlans(Exception):
    pass

EPSILON = 1e-4

def findManeuver(constellation,s,p,next_a,date_next_a,modeleDeTransition,ensure_feasible):
        duree = constellation.getSatellite(s).getActivite(p).getDuree()
        points_controle = modeleDeTransition.getPointsDeControle(p,next_a)
        fin_fenetre = constellation.getSatellite(s).getActivite(p).getFin()
        

        
        """ ce cas n'aurait pas du arriver, mais il y a un soucis dans
            les donnees. Il peut arriver qu'il manque un points ..."""
        if points_controle[-1][0]<fin_fenetre:
            points_controle.append((fin_fenetre,points_controle[-1][1]))
            
        if p==558 and next_a==1282:
            print(points_controle)
        # cas ou la transition est toujours faisable
        if points_controle == [(0,0)]:
            return fin_fenetre,0      
        # cas ou la maneuver au plus tard arrive avant date_next_a
        end_window_maneuver = constellation.getSatellite(s).getActivite(p).getFin()
        latest_maneuver = constellation.getSatellite(s).getTransitionTimeDependent(p,next_a,end_window_maneuver,modeleDeTransition)
        if latest_maneuver+end_window_maneuver<=date_next_a:
            assert(latest_maneuver<modeleDeTransition.getMajorantDureeTransition())
            return end_window_maneuver,latest_maneuver
        end_window_maneuver = points_controle[-1][0]
        # cas ou partir au debut de la fenetre ne suffit pas : elager le if suivant
        if points_controle[0][0]+points_controle[0][1]>=date_next_a:
            msg = "transition "+str(p)+"->"+str(next_a)+" tau="+str(date_next_a)+"\n"
            msg += "Start_p + tau(p) >= date_next_a"
            raise InfeasibleSolutionException(msg)
        # cas ou partir au debut de la fenetre + duree (au plus tot) ne suffit pas
        start_earliest_maneuver = points_controle[0][0]+duree
        earliest_maneuver = constellation.getSatellite(s).getTransitionTimeDependent(p,next_a,start_earliest_maneuver,modeleDeTransition)
        if earliest_maneuver+start_earliest_maneuver>date_next_a:
            assert(earliest_maneuver<modeleDeTransition.getMajorantDureeTransition())
            if not ensure_feasible:
                msg = "!earliest_maneuver+start_earliest_maneuver>=date_next_a\n"
                msg += "transition "+str(p)+"->"+str(next_a)+"\n"
                msg += str("earliest_maneuver="+str(earliest_maneuver)+"\n")
                msg += str("start_earliest_maneuver="+str(start_earliest_maneuver)+"\n")
                msg += str("date_next_a="+str(date_next_a)+"\n")
                msg += str("earliest_maneuver+start_earliest_maneuver="+str(earliest_maneuver+start_earliest_maneuver)+"\n")
                raise InfeasibleSolutionException(msg)
            else:
                raise NumericalError()
                #return earliest_maneuver - ecart,start_earliest_maneuver
        # ici on cherche le segment sur lequel on va basculer dans le cas infaisable
        left,right = searchChangeSegment(points_controle,date_next_a,duree)
        if left==right:
            return points_controle[left]
        else:
            points_left = points_controle[left]
            points_right = points_controle[right]
            t = inverseLinearManeuver(date_next_a,points_left,points_right,duree)
            maneuver = constellation.getSatellite(s).getTransitionTimeDependent(p,next_a,t,modeleDeTransition)
            if not(t+maneuver<=date_next_a):
                msg = "Manoeuvre "+str(p)+"->"+str(next_a)
                msg += ". left="+str(points_left)+",right="+str(points_right)+ "points ="+str(points_controle)
                msg += " start="+str(t)+" duree="+str(maneuver)
                assertLessThan(t+maneuver,date_next_a+EPSILON,msg)
                raise NumericalError()
            return t,maneuver

def calculerTransitionToNextTimeDep(constellation,s,p,next_a,date_next_a,modeleDeTransition,ensure_feasible=False):
    assert(modeleDeTransition.estTimeDependent())
    start_maneuver,duration = findManeuver(constellation,s,p,next_a,date_next_a,modeleDeTransition,ensure_feasible)
    computed_maneuver = constellation.getSatellite(s).getTransitionTimeDependent(p,next_a,start_maneuver,modeleDeTransition)
    if not(start_maneuver+computed_maneuver<=date_next_a):
        raise NumericalError()
    if not(start_maneuver+duration<=date_next_a):
        assert(start_maneuver+duration<=date_next_a+EPSILON)
        raise NumericalError()
    assert(duration==computed_maneuver)
    
    return  (start_maneuver,duration)
                    

class SolCCA:  
    def __init__(self,identifiant,satellite):
        self.identifiant = identifiant
        self.satellite = satellite
        self.debut = np.Inf
        self.fin = -np.Inf
        #self.activites = []
        self.sequence = []
        self.temps_LKH = []
        self.planEst = []
        self.planLts = []
        self.notifierPlanAJour(True,True)
        self.couts_transition = 0
        self.couts_transition_a_jour = True
        self.position_dernier_insere = None
        self.retard = 0

    def plansAJour(self):
        return self.planEarliestAJour and self.planLatestAJour
   
    def renommer(self,new):
        self.identifiant = new
        
    def copieTemporaire(self):
        copie = SolCCA(self.identifiant,self.satellite)
        copie.debut = self.debut
        copie.fin = self.fin
        copie.planEst = self.planEst.copy()
        copie.planLts = self.planLts.copy()
        copie.notifierPlanAJour(self.planEarliestAJour,self.planLatestAJour)
        copie.sequence = self.sequence.copy()
        copie.retard = self.retard
        copie.couts_transition = self.couts_transition
        copie.couts_transition_a_jour = self.couts_transition_a_jour
        copie.position_dernier_insere = self.position_dernier_insere
        return copie
    
    def formaterInfo(self):
        msg = "CCA "+str(self.identifiant)
        msg += "(sat. "+str(self.satellite) +")"
        msg += " € ["+ str(self.debut)+"," +str(self.fin) +"]"
        return msg
        
    def reset(self):
        self.debut = np.Inf
        self.fin = -np.Inf
        #self.activites = []
        self.sequence = []
        self.planEst = []
        self.planLts = []
        self.notifierPlanAJour(True,True)
        self.couts_transition = 0
        self.couts_transition_a_jour = True
        self.position_dernier_insere = None
        
    def scoreCharge(self,methode='densite'):
        if methode == 'densite':
            return len(self.sequence)/(self.fin-self.debut)
        elif methode == 'cardinal':
            return len(self.sequence)
        else:
            raise ValueError('Méthode de scoring de CCA inconnu')
    
    def scoreChargeAjoutActivites(self,nActivites,methode='densite'):
        if methode=='densite':
            return (len(self.getSequence())+nActivites)/(self.fin-self.debut)
        elif methode == 'cardinal':
            return len(self.sequence)+nActivites
        else:
            raise ValueError('Méthode de scoring de CCA inconnu')
        
    def renommer(self,new):
        self.identifiant = new
        
    def __str__(self):
        return "SolCCA(" + str(self.identifiant) +":"+str(self.sequence)+")"
        
    def notifierPlanAJour(self,early,late):
        self.planEarliestAJour = early
        self.planLatestAJour = late
        
    def retirerListeActivites(self,constellation,activites,modeleDeTransition):
        self.sequence = [a for a in self.sequence if a not in activites]
        self.planEst = [at for at in self.planEst if at[0] not in activites]
        self.planLts = [at for at in self.planLts if at[0] not in activites]
        self.notifierPlanAJour(False,False)
        self.couts_transition_a_jour = False
        self.MAJDebutFin(constellation,modeleDeTransition)
            
    def retirerActivite(self,constellation,a,modeleDeTransition):
        self.sequence.remove(a)
        self.notifierPlanAJour(False,False)
        self.couts_transition_a_jour = False
        self.MAJDebutFin(constellation,modeleDeTransition)
        
    def estVide(self):
        return len(self.activites)==0
    
    def annulerMode(self,constellation,r,m,modeleDeTransition):
        if self.sequence is not None:
            for (s,o) in constellation.getRequete(r).getMode(m).getCouples():
                if o in self.sequence:
                    self.sequence.remove(o)
        self.notifierPlanAJour(False,False)
        self.couts_transition_a_jour = False
        self.MAJDebutFin(constellation,modeleDeTransition)
    
    def MAJDebutFin(self,constellation,modeleDeTransition):
        tau_max = modeleDeTransition.getMajorantDureeTransition()
        try:
            self.debut = min([constellation.getSatellite(self.satellite).getActivite(aa).getDebut() for aa in self.sequence])
            self.fin = max([constellation.getSatellite(self.satellite).getActivite(aa).getFin()+tau_max for aa in self.sequence])
        except ValueError:
            self.debut = np.Inf
            self.fin = - np.Inf
            
    def setSequence(self,constellation,seq,modeleDeTransition):
        self.sequence = deepcopy(seq)
        self.MAJDebutFin(constellation,modeleDeTransition)
        self.notifierPlanAJour(False,False)
        self.couts_transition_a_jour = False
        
    def setSequenceIndex(self,i_start,i_end,seq,modeleDeTransition):
        #for a in seq:
        #    assert(a in self.activites)
        if i_end<len(self.sequence):
            self.sequence[i_start:] = seq + self.sequence[i_end+1:]
        else:
            self.sequence[i_start:] = seq
        self.MAJDebutFin(constellation,modeleDeTransition)
        self.notifierPlanAJour(False,False)
        self.couts_transition_a_jour = False
        
    """
        ===============================================
                        GETTERS
        ===============================================
    """
    def contient(self,a):
        return a in self.sequence
    
    def getIdentifiant(self):
        return self.identifiant
    
    def getSatellite(self):
        return self.satellite
    
    def getSequence(self):
        return self.sequence
    
    #def getActivites(self):
    #    return self.activites
    
    def getDebut(self):
        return self.debut
    
    def getFin(self):
        return self.fin
    
    # return iLeft,iRight
    # planEarly[:iLeft] = activites a gauche
    # si iLeft<len(planEarly) : planEarly[iLeft:iRight] = activites en conflit
    # si iRight<len(planLatest) : planLatest[iRight:] : activites a droite
    def conflitPlan(self,constellation,activite,modeleDeTransition):
        planEarly,planLatest = self.planEst,self.planLts
        if config.getOptValue("verif"):
            for i,m in enumerate(planEarly):
                if not planEarly[i][1]<=planLatest[i][1]:
                    print(self.identifiant,activite,planEarly,planLatest)
                    print(planEarly[i],planLatest[i])
                    assert(False)
        
        # cas vide
        if len(planEarly)==0:
            return 0,0
        
        debut = activite.getDebut()
        fin = activite.getFin()
        iLeft = len(planLatest) # position acceptable minimum
        iRight = 0 # position acceptable maximum
        conflitIndice = lambda i : intersect( (planEarly[i][1],planLatest[i][1]),(debut,fin))
        isLeft = lambda i : planLatest[i][1] < debut
        isRight = lambda i : fin < planEarly[i][1]
        for i in range(len(planEarly)+1):
            if i == len(planEarly):
                test_right = True
            else:
                test_right = not isLeft(i)
            if iLeft == 0:
                test_left = True
            else:
                test_left = not isRight(i-1)
            if test_left and test_right:
                iLeft = min(iLeft,i)
                iRight = max(iRight,i)
        return iLeft,iRight        

    # TODO : propagation
    def MAJPlansCritiques(self,constellation,modeleDeTransition,force=False):
        #assert(not modeleDeTransition.estTimeDependent())
        if not self.plansAJour() or force:
            self._calculerPlanEarliest(constellation,modeleDeTransition)
            self._calculerPlanLatest(constellation,modeleDeTransition)
            if config.verifMode():
                for i in range(len(self.planEst)):
                    if not(self.planEst[i][1]<=self.planLts[i][1]):
                        die(self.planEst[i][1],self.planLts[i][1])
            
            
        return self.retard  <= 10**(-config.glob.digits)     
            
    def _calculerPlanEarliest(self,constellation,modeleDeTransition):
        self.planEst,self.couts_transition,self.retard = self.planEarliest(constellation,self.sequence,modeleDeTransition,add_transition_tardiness=True)
        self.planEarliestAJour = True
        self.couts_transition_a_jour = True
        assert(len(self.sequence)==len(self.planEst))
        
    def _calculerPlanLatest(self,constellation,modeleDeTransition):
        #assert(not modeleDeTransition.estTimeDependent())
        self.planLts = self.planLatest(constellation,self.sequence,modeleDeTransition)
        self.planLatestAJour = True
  
    def getEarliestTransition(self,constellation,s,a1,a2,start,modeleDeTransition):
        if modeleDeTransition.estTimeDependent():
            transition = constellation.getSatellite(s).getTransitionTimeDependent(a1,a2,start,modeleDeTransition)
        else:
            transition = constellation.getSatellite(s).getTransition(a1,a2,modeleDeTransition)                    
        return transition             

    def planEarliest(self,constellation,sequence,modeleDeTransition,add_transition_tardiness=False,add_transition=False):
        try:
            transition_totale = 0
            plan = []
            retard_total = 0
            s = self.satellite
            for i,activite in enumerate(sequence):
                if(i==0):
                    t = constellation.satellites[s].getActivite(activite).getDebut()
                    
                else:
                    prec = sequence[i-1]
                    duree = constellation.getSatellite(s).getActivite(prec).getDuree()
                    start = constellation.satellites[s].getActivite(activite).getDebut()
                    end = constellation.satellites[s].getActivite(activite).getFin()
                    transition = self.getEarliestTransition(constellation,s,prec,activite,t+duree,modeleDeTransition)
                    transition_totale += transition
                    retard_total +=  round(max(t + duree + transition - end,0),config.glob.digits)
                    t = round(max(t + duree + transition,start),config.glob.digits)
                if add_transition and i>0:
                    plan[-1] = (plan[-1][0],plan[-1][1],transition)
                plan.append((activite,t))
            if add_transition_tardiness:
                return plan,transition_totale,retard_total
            else:
                return plan
        except KeyError:
            print(s,sequence)
            raise KeyError


    def planLatest(self,constellation,sequence,modeleDeTransition):
        #assert(not modeleDeTransition.estTimeDependent())
        plan = []
        s = self.satellite
        for i in range(len(sequence)-1,-1,-1):
            a = sequence[i]
            duree = constellation.getSatellite(s).getActivite(a).getDuree()
            if(i==len(sequence)-1):
                t = constellation.getSatellite(s).getActivite(a).getFin() - duree
                transition = 0
            else:
                suivante = plan[-1][0]
                if modeleDeTransition.estTimeDependent():
                    start,transition = calculerTransitionToNextTimeDep(constellation,s,a,suivante,plan[-1][1]-10**(-6),modeleDeTransition,ensure_feasible=True)
                    if a==558 and suivante==1282:
                        print("=========================",a,suivante,start,transition,plan)
                    
                    transition = transition
                    activity_start = start-duree
                    t = min(activity_start,constellation.getSatellite(s).getActivite(a).getFin()-duree)
                    if t<self.planEst[i][1]:
                        self.notifierPlanAJour(False,False)
                        raise NumericalError()
                else:
                    transition = constellation.getSatellite(s).getTransition(a,suivante,modeleDeTransition)
                    tbefore = t
                    t = min(t - duree - transition,constellation.getSatellite(s).getActivite(a).getFin()-duree)
                # correction arrondi
                
                """
                if modeleDeTransition.estTimeDependent():
                    trans = constellation.getSatellite(s).getTransitionTimeDependent(a,suivante,t+duree,modeleDeTransition)
                    while t + trans + duree >  plan[-1][1]:
                        ecart = (t + trans + duree - plan[-1][1])
                        if ecart>0:
                            print(ecart)
                        t = arrondi(t-(ecart),4,up=False)
                        if not(ecart<1):
                            trans = constellation.getSatellite(s).getTransitionTimeDependent(a,suivante,t,modeleDeTransition)
                            #print(duree)
                            #print(transition)
                            transition = constellation.getSatellite(s).getTransitionTimeDependent(a,suivante,activity_start,modeleDeTransition)
                            print(transition)
                            die(t,activity_start,t+trans+duree,activity_start+transition+duree,plan[-1][1])
                            die(ecart,"le couple n'est peut etre pas FIFO",modeleDeTransition.isFIFO(a,suivante)) 
                        trans = constellation.getSatellite(s).getTransitionTimeDependent(a,suivante,t,modeleDeTransition)
                """
                if modeleDeTransition.estTimeDependent():
                    trans = constellation.getSatellite(s).getTransitionTimeDependent(a,suivante,t+duree,modeleDeTransition)
                    assert(t+trans<=start+transition)
            if config.verifMode():
                if len(plan)>0:
                    assert(t+constellation.getSatellite(s).getTransitionTimeDependent(a,suivante,t+duree,modeleDeTransition)<=plan[-1][1])
            plan.append((a,t))
  
        plan.reverse()
        if config.verifMode():
            for i,(a,t) in enumerate(plan):
                duree = constellation.getSatellite(s).getActivite(a).getDuree()
                assert(t+duree<=constellation.getSatellite(s).getActivite(a).getFin())
                if i<len(plan)-1:
                    suivante = plan[i+1][0]
                    assert(t+duree+constellation.getSatellite(s).getTransitionTimeDependent(a,suivante,t+duree,modeleDeTransition)<=plan[i+1][1])
        return plan
    

    
    """
        =============================================== 
                        RECHERCHE LOCALE
        ===============================================        
    """
    def rechercheLocale(self,constellation,solver,modeleDeTransition,groups=None,allow_no_solution=False):
        if groups is None:
            assert(solver == 'LKH')
        else:
            assert(solver in ['OPTW','OPTWGroups'])
        if solver == 'LKH' and self.sequence == []:
            return
        #printColor('len ',MPI.COMM_WORLD.Get_rank(),len(self.sequence),c='y')
        depart = min([constellation.getSatellite(self.satellite).getActivite(a).getDebut() for a in self.sequence])
        #temps = time()
        if solver=='LKH':
            printOpen("Recherche locale cca",self.identifiant,c='y')
            printOpen("Write LKH cca",self.identifiant,c='c')
            self.writeLKH(constellation,self.sequence,depart,modeleDeTransition)
            printClose()
            printOpen("Exec LKH cca",self.identifiant,c='c')
            self.execLKH()
            printClose()
            printOpen("Read LKH cca",self.identifiant,c='c')
            res = self.readLKH()
            printClose()
            printClose()
            if res is not None:
                self.sequence = res
            self.cleanFolderLKH()
        elif solver=='OPTWGroups':
            t0 = time()
            if modeleDeTransition.utiliseSolverTimeDep():
                self.writeOPTWGroupsTimeDep(constellation,self.sequence,depart,groups,modeleDeTransition)
            else:
                self.writeOPTWGroups(constellation,self.sequence,depart,groups,modeleDeTransition)
            t1 = time()
            if modeleDeTransition.utiliseSolverTimeDep():
                retval = self.execOPTWGroupsTimeDep(modeleDeTransition)
                if retval!=0:
                    die(self.last_call_optw)
            else:
                retval = self.execOPTWGroups(modeleDeTransition)
            t2 = time()
            if modeleDeTransition.utiliseSolverTimeDep():
                res,r_satisfaites,earliest,latest = self.readOPTWGroupsTimeDep(constellation,groups,allow_no_solution=allow_no_solution)
            else:
                res,r_satisfaites = self.readOPTWGroups(constellation,groups,allow_no_solution=allow_no_solution)
            t3 = time()
            if res is not None:
                for i,a in enumerate(res):
                    find = False
                    for g in groups:
                        if a in groups[g]:
                            find = True
                            break
                    if not find:
                        if i<len(res)-1:
                            print("Activités inconnues",a)
                            print("sequence",self.sequence)
                            print("groupes",groups)
                            print(res)
                            name = self.last_call_optw
                            print("fichier",name)
                            assert(find)
                        else:
                            warn("erreur solver : activité fictive",a)
            if retval == 0:
                self.cleanFolderOPTWGroups()
            else:
                print(self.last_call_optw+" failed")
            copy = deepcopy(self.sequence)
            self.sequence = res
            # patch du bug du solver : backup de la solution
            if not self.sequenceFaisable(constellation,modeleDeTransition):
                raise InfeasibleSolutionfromSolverException
            return list(r_satisfaites.keys())
        elif solver=='OPTW':
            t0 = time()
            self.writeOPTW(constellation,self.sequence,depart,groups,modeleDeTransition)
            t1 = time()
            self.execOPTW()
            t2 = time()
            res,r_satisfaites = self.readOPTW(constellation,groups,allow_no_solution=allow_no_solution)
            t3 = time()
            self.cleanFolderOPTW()
            if res is not None:
                for a in res:
                    find = False
                    for g in groups:
                        if a in groups[g]:
                            find = True
                            break
                    if not find:
                        print("Activités inconnues",a)
                        print("sequence",self.sequence)
                        print("groupes",groups)
                        print(res)
                        name = self.last_call_optw
                        print("fichier",name)
                    assert(find)
                self.sequence = res
            t4 = time()
            #print(round(t1-t0,3),round(t2-t1,3),round(t3-t2,3),round(t4-t3,3))
            return list(r_satisfaites.keys())

        #self.temps_LKH.append(time()-temps)
        #printColor('len apres',MPI.COMM_WORLD.Get_rank(),len(self.sequence),c='y')
    
    """ OPTW with groups time-dependent"""
        
    def execOPTWGroupsTimeDep(self,modeleDeTransition):
        solver_path = "./SingleSatSolver"
        output_path = "OPTW/TOURS/eosm/"+self.last_call_optw+".tour"
        goDir = "cd ../OPTW;"
        filename = "comp_"+str(self.satellite)+"_"+str(self.identifiant)+".comp"
        cmd = goDir + solver_path
        cmd += " -multiTimeSteps "
        cmd += " -fRequests OPTW/INSTANCES/eosm/" + self.last_call_optw + ".optw"
        cmd += " -fComponent " + modeleDeTransition.getComponentsPath()+"/"+filename
        cmd += " -tmax "+str(config.OPTW.tmax) +" -restart " +str(config.getOptValue("restart"))
        cmd += " -tmaxGlobal " +str(config.OPTW.tmaxGlobal)
        cmd += " -o " +str(output_path)
        cmd += " -dates"
        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        for line in p.stdout.readlines():
            if config.verifMode():
                print(line)
        retval = p.wait()
        return retval
    
    def writeGroupsOPTWGroupsTimeDep(self,groups,file):
        file.write(str(len(groups.keys())-1)+'\n')
        for g in groups:
            if g>0:
                file.write(str(len(groups[g])))
                for a in groups[g]:
                    file.write(" "+str(a))
                file.write("\n")
        
    # groups : id groupe => liste d'activites
    def writeOPTWGroupsTimeDep(self,constellation,sequence,depart,groups,modeleDeTransition):
        for g in groups:
            for a in groups[g]:
                assert(a in sequence)
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        ECHELLE = config.OPTW.echelle
        #transition = lambda a1,a2 : constellation.getSatellite(self.satellite).getTransition(a1,a2,modeleDeTransition)
        s = self.satellite
        arrivee = int(ECHELLE*max([constellation.getSatellite(self.satellite).getActivite(a).getFin()-constellation.getSatellite(s).getActivite(a).getDuree() for a in self.sequence]))
        #self.mapping_id_optw = {}
        name = config.donnees.algo_name + "_" + str(rank) + "_" + socket.gethostname() + '_sat_' + str(self.satellite) + '_cca_' + str(self.identifiant)
        folder = "../OPTW/OPTW/INSTANCES/eosm/"
        tau_max = modeleDeTransition.getMajorantDureeTransition()
        
        with open(folder + name + '.optw','w') as file:
            self.writeSeqInitOPTWGroups(groups[0],file,modeleDeTransition)
            self.writeGroupsOPTWGroupsTimeDep(groups,file)
                   
        self.last_call_optw = name
        
    def readOPTWGroupsTimeDep(self,constellation,groups,allow_no_solution=False):
        name = self.last_call_optw
        file = "../OPTW/OPTW/TOURS/eosm/"+name+".tour"
        with open(file,'r') as f:
            lines = f.read().splitlines()
            if allow_no_solution and len(lines)==0:
                activites = []
            else:
                assert(len(lines)==3)
                line = lines[0]
                activites = [int(a) for a in line.split(" ")]
            r_satisfaites = {}
            for a in activites:
                r = constellation.getRequeteActivite(a)
                if r not in r_satisfaites:
                    r_satisfaites[r] = 1
                else:
                    r_satisfaites[r] += 1
            if config.getOptValue("verif") and 0 in groups:
                for a in groups[0]:
                    assert(constellation.getRequeteActivite(a) in r_satisfaites)
            earliest = [int(a) for a in lines[1].split(" ")]
            latest = [int(a) for a in lines[2].split(" ")]
        return activites,r_satisfaites,earliest,latest
        
    """ OPTW with GROUPS """
        
    def execOPTWGroups(self,modeleDeTransition):
        solver_path = "./OrienteeringWithGroupsSolver"
        output_path = "OPTW/TOURS/eosm/"+self.last_call_optw+".tour"
        goDir = "cd ../OPTW;"
        cmd = goDir + solver_path
        cmd += " -format ccati "
        cmd += "-f OPTW/INSTANCES/eosm/" + self.last_call_optw + ".optw"
        cmd += " -tmax "+str(config.OPTW.tmax) +" -restart " +str(config.getOptValue("restart"))
        cmd += " -tmaxGlobal " +str(config.OPTW.tmaxGlobal)
        cmd += " -o " +str(output_path)
        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        for line in p.stdout.readlines():
            if config.verifMode():
                print(line)
        retval = p.wait()
        return retval

    def cleanFolderOPTWGroups(self):
        try:
            folder_instance = "../OPTW/OPTW/INSTANCES/eosm/"
            os.unlink(folder_instance+self.last_call_optw+".optw")
        except FileNotFoundError:
            pass
        try:
            #for tour in [filename for filename in os.listdir("../LKH3/TSPTW/TOURS/eosm/") if filename.startswith(self.last_call_lkh)]:
                folder_tour = "../OPTW/OPTW/TOURS/eosm/"
                os.unlink(folder_tour+self.last_call_optw+'.tour')
        except FileNotFoundError:
            pass
        try:
            folder_tour = "../OPTW/OPTW/TSPTW/PAR/"
            os.unlink(folder_tour+self.last_call_optw+'.par')
        except FileNotFoundError:
            pass
        try:
            folder_tour = "../OPTW/OPTW/TSPTW/INIT/"
            os.unlink(folder_tour+self.last_call_optw+'.init')
        except FileNotFoundError:
            pass
    
    def readOPTWGroups(self,constellation,groups,allow_no_solution=False):
        name = self.last_call_optw
        file = "../OPTW/OPTW/TOURS/eosm/"+name+".tour"
        with open(file,'r') as f:
            lines = f.read().splitlines()
            if allow_no_solution and len(lines)==0:
                activites = []
            else:
                assert(len(lines)==1)
                line = lines[0]
                activites = [int(a) for a in line.split(" ")]
            r_satisfaites = {}
            for a in activites:
                r = constellation.getRequeteActivite(a)
                if r not in r_satisfaites:
                    r_satisfaites[r] = 1
                else:
                    r_satisfaites[r] += 1
            if config.getOptValue("verif") and 0 in groups:
                for a in groups[0]:
                    assert(constellation.getRequeteActivite(a) in r_satisfaites)
        return activites,r_satisfaites
    
    def writeActiviteOPTWGroups(self,constellation,file,a,g,modeleDeTransition):
        ECHELLE = config.OPTW.echelle
        file.write(str(a)+" ")
        #file.write(str(g)+" ")
        coord = constellation.getSatellite(self.satellite).getActivite(a).getCoordonnees()
        vlat = modeleDeTransition.getVitesseLatitude()
        vlong = modeleDeTransition.getVitesseLongitude()
        lat,long = coord[0],coord[1]
        lat_effective,long_effective = lat/vlat,long/vlong
        file.write(str(long_effective)+" ")
        file.write(str(lat_effective)+" ")
        duree = constellation.getSatellite(self.satellite).getActivite(a).getDuree()
        file.write(str(int(ECHELLE*math.ceil(duree)))+" ")
        debut = constellation.getSatellite(self.satellite).getActivite(a).getDebut()
        file.write(str(int(ECHELLE*math.ceil(debut)))+" ")
        fin = constellation.getSatellite(self.satellite).getActivite(a).getFin()-constellation.getSatellite(self.satellite).getActivite(a).getDuree()
        file.write(str(int(ECHELLE*math.floor(fin)))+'\n')
        
        
    def writeSeqInitOPTWGroups(self,seq_init,file,modeleDeTransition):
        file.write(str(len(seq_init)))
        for a in seq_init:
            file.write(" "+str(a))
        file.write("\n")
    
    def writeGroupsOPTWGroups(self,groups,file):
        file.write(str(len(groups.keys())-1)+'\n')
        for g in groups:
            if g>0:
                file.write(str(len(groups[g])))
                for a in groups[g]:
                    file.write(" "+str(a))
                file.write("\n")
        
    # groups : id groupe => liste d'activites
    def writeOPTWGroups(self,constellation,sequence,depart,groups,modeleDeTransition):
        for g in groups:
            for a in groups[g]:
                assert(a in sequence)
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        ECHELLE = config.OPTW.echelle
        #transition = lambda a1,a2 : constellation.getSatellite(self.satellite).getTransition(a1,a2,modeleDeTransition)
        s = self.satellite
        arrivee = int(ECHELLE*max([constellation.getSatellite(self.satellite).getActivite(a).getFin()-constellation.getSatellite(s).getActivite(a).getDuree() for a in self.sequence]))
        #self.mapping_id_optw = {}
        name = config.donnees.algo_name + "_" + str(rank) + "_" + socket.gethostname() + '_sat_' + str(self.satellite) + '_cca_' + str(self.identifiant)
        folder = "../OPTW/OPTW/INSTANCES/eosm/"
        tau_max = modeleDeTransition.getMajorantDureeTransition()

        with open(folder + name + '.optw','w') as file:
            if not modeleDeTransition.utiliseSolverTimeDep():
                Nnoeuds = len(sequence)
                file.write(str(Nnoeuds)+'\n')
                file.write(str(1/ECHELLE)+" "+str(tau_max*ECHELLE)+'\n')
                for g in groups:
                    for a in groups[g]:
                        self.writeActiviteOPTWGroups(constellation,file,a,g,modeleDeTransition)
                        #self.mapping_id_optw[a] = a
            self.writeSeqInitOPTWGroups(groups[0],file,modeleDeTransition)
            self.writeGroupsOPTWGroups(groups,file)
                    
        self.last_call_optw = name
       
        
    """       OPTW classique """
        
    def execOPTW(self):
        tmaxGlobal = config.OPTW.tmaxGlobal
        solver_path = "./OrienteeringSolver"
        output_path = "OPTW/TOURS/eosm/"+self.last_call_optw+".tour"
        goDir = "cd ../OPTW;"
        cmd = goDir + solver_path
        cmd += " -format ccati "
        cmd += "-f OPTW/INSTANCES/eosm/" + self.last_call_optw + ".optw"
        cmd += " -tmax "+str(config.OPTW.tmax) +" -restart " +str(config.OPTW.restart)
        cmd += " -tmaxGlobal " +str(tmaxGlobal)
        cmd += " -o " +str(output_path)
        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        for line in p.stdout.readlines():
            print(line)
        retval = p.wait()
    
    
    def writeActiviteOPTW(self,constellation,file,a,g):
        ECHELLE = config.OPTW.echelle
        file.write(str(a)+" ")
        #file.write(str(g)+" ")
        coord = constellation.getSatellite(self.satellite).getActivite(a).getCoordonnees()
        long,lat = coord[0],coord[1]
        file.write(str(long)+" ")
        file.write(str(lat)+" ")
        duree = constellation.getSatellite(self.satellite).getActivite(a).getDuree()
        file.write(str(int(math.ceil(duree*ECHELLE)))+" ")
        debut = constellation.getSatellite(self.satellite).getActivite(a).getDebut()
        file.write(str(int(math.ceil(debut*ECHELLE)))+" ")
        fin = constellation.getSatellite(self.satellite).getActivite(a).getFin()-constellation.getSatellite(self.satellite).getActivite(a).getDuree()
        file.write(str(int(math.floor(fin*ECHELLE)))+'\n')
    
    def writeSeqInitOPTW(self,seq_init,file):
        file.write(str(len(seq_init)))
        for a in seq_init:
            file.write(" "+str(a))
        file.write("\n")
    
    def writeGroupsOPTW(self,groups,file):
        file.write(str(len(groups.keys())-1)+'\n')
        for g in groups:
            if g>0:
                file.write(str(len(groups[g])))
                for a in groups[g]:
                    file.write(" "+str(a))
                file.write("\n")
        
    # groups : id groupe => liste d'activites
    def writeOPTW(self,constellation,sequence,depart,groups,modeleDeTransition):
        for g in groups:
            for a in groups[g]:
                assert(a in sequence)
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        ECHELLE = config.OPTW.echelle
        
        transition = lambda a1,a2 : constellation.getSatellite(self.satellite).getTransition(a1,a2)
        s = self.satellite
        arrivee = int(ECHELLE*max([constellation.getSatellite(self.satellite).getActivite(a).getFin()-constellation.getSatellite(s).getActivite(a).getDuree() for a in self.sequence]))
        self.mapping_id_optw = {}
        name = config.donnees.algo_name + "_" + str(rank) + "_" + socket.gethostname() + '_sat_' + str(self.satellite) + '_cca_' + str(self.identifiant)
        folder = "../OPTW/OPTW/INSTANCES/eosm/"
        tau_max = modeleDeTransition.getMajorantDureeTransition()
        with open(folder + name + '.optw','w') as file:
            Nnoeuds = len(sequence)
            file.write(str(Nnoeuds)+'\n')
            file.write(str(1/ECHELLE)+" "+str(tau_max*ECHELLE)+'\n')

            for g in groups:
                for a in groups[g]:
                    self.writeActiviteOPTW(constellation,file,a,g)
                    self.mapping_id_optw[a] = a
            self.writeSeqInitOPTW(groups[0],file)
            self.writeGroupsOPTW(groups,file)
                    
        self.last_call_optw = name
        #self.writeParFile()
        #self.writeInitialTourFile(seq)
    
    
    def cleanFolderOPTW(self):
        try:
            folder_instance = "../OPTW/OPTW/INSTANCES/eosm/"
            os.unlink(folder_instance+self.last_call_optw+".optw")
        except FileNotFoundError:
            pass
        try:
            #for tour in [filename for filename in os.listdir("../LKH3/TSPTW/TOURS/eosm/") if filename.startswith(self.last_call_lkh)]:
                folder_tour = "../OPTW/OPTW/TOURS/eosm/"
                os.unlink(folder_tour+self.last_call_optw+'.tour')
        except FileNotFoundError:
            pass
        try:
            folder_tour = "../OPTW/OPTW/TSPTW/PAR/"
            os.unlink(folder_tour+self.last_call_optw+'.par')
        except FileNotFoundError:
            pass
        try:
            folder_tour = "../OPTW/OPTW/TSPTW/INIT/"
            os.unlink(folder_tour+self.last_call_optw+'.init')
        except FileNotFoundError:
            pass
    
    def readOPTW(self,constellation,groups,allow_no_solution=False):
        name = self.last_call_optw
        file = "../OPTW/OPTW/TOURS/eosm/"+name+".tour"
        with open(file,'r') as f:
            lines = f.read().splitlines()
            if allow_no_solution and len(lines)==0:
                activites = []
            else:
                assert(len(lines)==1)
                line = lines[0]
                activites = [int(a) for a in line.split(" ")]
            r_satisfaites = {}
            for a in activites:
                r = constellation.getRequeteActivite(a)
                if r not in r_satisfaites:
                    r_satisfaites[r] = 1
                else:
                    r_satisfaites[r] += 1
            if config.getOptValue("verif") and 0 in groups:
                for a in groups[0]:
                    assert(constellation.getRequeteActivite(a) in r_satisfaites)
        return activites,r_satisfaites

    
   
    """      LKH """
       
    # Ecrire le fichier de parametre pour LKH
    def writeParFile(self):
        name = self.last_call_lkh+".par"
        folder = "../LKH3/TSPTW/PAR/"
        with open(folder+name,'w') as file:
            file.write("SPECIAL\n")
            file.write("PROBLEM_FILE = INSTANCES/eosm/"+self.last_call_lkh+".tsptw\n")
            file.write("MAX_TRIALS = "  + str(config.lkh.max_trials)+'\n')
            file.write("RUNS = " + str(config.lkh.runs)+'\n')
            file.write("TRACE_LEVEL = 0\n")
            #file.write("STOP_AT_OPTIMUM = YES\n")
            file.write("OUTPUT_TOUR_FILE = TOURS/eosm/"+self.last_call_lkh+".tour\n")
            file.write("INITIAL_TOUR_FILE = INIT/"+self.last_call_lkh+".init\n")
            
    # Ecrire le fichier de solution initiale pour LKH
    def writeInitialTourFile(self,sequence):
        folder = "../LKH3/TSPTW/INIT/"
        #print(folder+self.last_call_lkh+".init\n")
        with open(folder+self.last_call_lkh+".init",'w') as file:
            file.write('NAME :'+self.last_call_lkh+".init\n")
            file.write("TYPE : TOUR\n")
            file.write("DIMENSION : " + str(len(sequence)+1) +"\n")
            file.write("TOUR_SECTION\n")
            file.write("1\n")
            for a in sequence:
                for id_a,act_a in self.mapping_id_lkh.items():
                    if act_a==a:
                        file.write(str(id_a)+"\n")
                        break
            file.write("-1\n")
            file.write("EOF\n")
    
    # Ecrire e fichier d'entrée pour LKH
    def writeLKH(self,constellation,seq,depart,modeleDeTransition):
        ECHELLE = config.lkh.echelle
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        depart_echelle = ECHELLE*depart
        transition = lambda a1,a2 : constellation.getSatellite(self.satellite).getTransition(a1,a2,modeleDeTransition)
        s = self.satellite
        arrivee = int(ECHELLE*max([constellation.getSatellite(self.satellite).getActivite(a).getFin()-constellation.getSatellite(s).getActivite(a).getDuree() for a in self.sequence]))
        mapping_id = {}
        sequence = [-1] + seq
        name = config.getOptValue("solver") + "_" + str(rank) + "_" + socket.gethostname() + '_sat_' + str(self.satellite) + '_cca_' + str(self.identifiant)
        folder = "../LKH3/TSPTW/INSTANCES/eosm/"
        
        with open(folder + name + '.tsptw','w') as file:
            Nnoeuds = len(sequence)
            #printColor('ecriture',MPI.COMM_WORLD.Get_rank(),Nnoeuds,c='y')
            file.write('NAME : '+name+'.tsptw\n')
            file.write('TYPE : TSPTW\n')
            file.write('DIMENSION : ' + str(Nnoeuds) + '\n')
            file.write('EDGE_WEIGHT_TYPE : EXPLICIT\n')
            file.write('EDGE_WEIGHT_FORMAT : FULL_MATRIX\n')
            file.write('EDGE_WEIGHT_SECTION\n')

            # distance entre chaque noeuds : transition + duree
            for a1 in sequence:
                for a2 in sequence:
                    if a1==a2 or a1==-1 or a2==-1:
                        file.write('0 ')
                    else:
                        duree = constellation.getSatellite(self.satellite).getActivite(a1).getDuree()
                        file.write(str(int(ECHELLE*(transition(a1,a2)+duree)))+' ')
                file.write('\n')
            file.write('TIME_WINDOW_SECTION\n')
            for i,a in enumerate(sequence):
                if a==-1:             # fenetre = (debut,fin-duree)
                    file.write('1 ' + str(depart_echelle) + ' ' + str(arrivee) + '\n')
                else:
                    debut = constellation.getSatellite(self.satellite).getActivite(a).getDebut()
                    fin = constellation.getSatellite(self.satellite).getActivite(a).getFin()-constellation.getSatellite(s).getActivite(a).getDuree()
                    file.write(str(i+1) + ' ' + str(int(ECHELLE*debut)) + ' ' + str(int(ECHELLE*fin)) + '\n')
                    #mapping_id[a] = a
                    mapping_id[i+1] = a
            file.write('DEPOT_SECTION\n')
            file.write('1\n')
            file.write('-1\n')
            file.write('EOF')
        self.last_call_lkh = name
        self.mapping_id_lkh = mapping_id
        self.writeParFile()
        self.writeInitialTourFile(seq)
        
    
    def execLKH(self):
        goDir = "cd ../LKH3/TSPTW ; "
        cmd = goDir + "../LKH PAR/" + self.last_call_lkh + ".par"
        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        for line in p.stdout.readlines():
            print(line)
        retval = p.wait()
 
    def cleanFolderLKH(self):
        try:
            folder_instance = "../LKH3/TSPTW/INSTANCES/eosm/"
            os.unlink(folder_instance+self.last_call_lkh+".tsptw")
        except FileNotFoundError:
            pass
        try:
            #for tour in [filename for filename in os.listdir("../LKH3/TSPTW/TOURS/eosm/") if filename.startswith(self.last_call_lkh)]:
                folder_tour = "../LKH3/TSPTW/TOURS/eosm/"
                os.unlink(folder_tour+self.last_call_lkh+'.tour')
        except FileNotFoundError:
            pass
        try:
            folder_tour = "../LKH3/TSPTW/PAR/"
            os.unlink(folder_tour+self.last_call_lkh+'.par')
        except FileNotFoundError:
            pass
        try:
            folder_tour = "../LKH3/TSPTW/INIT/"
            os.unlink(folder_tour+self.last_call_lkh+'.init')
        except FileNotFoundError:
            pass
    
        
    def readLKH(self):
        try:
            sequence = []
            name = self.last_call_lkh
            mapping = self.mapping_id_lkh
            file = "../LKH3/TSPTW/TOURS/eosm/"+name+".tour"
            if not exists(file):
                return None
            with open(file,'r') as file:
                lines = file.read().splitlines()
                i = 0
                while i < len(lines):
                    try:
                        entier = int(lines[i])
                        if entier==-1:
                            break
                        elif entier!=1:
                            sequence.append(self.mapping_id_lkh[entier])
                    except:
                        pass
                    i+=1
            if len(sequence)<len(self.sequence):
                return None
            return sequence
        except Exception as exc:
            print("Last call",self.last_call_lkh)
            print("Mapping id",self.mapping_id_lkh)
            print(exc)
            sequence = []
            name = self.last_call_lkh
            mapping = self.mapping_id_lkh
            prefixed = [filename for filename in os.listdir("../LKH3/TSPTW/TOURS/eosm/") if filename.startswith(name)]
            if len(prefixed)==0:
                return None
            best = None
            obj = np.Inf
            for filename in prefixed:
                tour = int(filename.split('.')[1])
                if tour<obj:
                    obj=tour
                    best = filename
            print("Echec de la lecture du fichier résultat",best)
            with open("../LKH3/TSPTW/TOURS/eosm/"+best,'r') as file:
                lines = file.read().splitlines()
                i = 0
                while int(lines[i])>0: #-1 -> fin du fichier
                    print(lines[i])
                    i+=1
            raise exc               
        
    def solverUsed(self):
        return self.use_solver
    
    def calculerSequence(self,constellation,modeleDeTransition,solver='LKH',groups=None):
        assert(solver in ['LKH','OPTW','OPTWGroups'])
        shiftRightDisplay(5)
        t1 = time()
        if not self.sequenceFaisable(constellation,modeleDeTransition):
            copy_sequence = deepcopy(self.sequence)
            try:
                self.use_solver = True
                res = self.rechercheLocale(constellation,solver,modeleDeTransition,groups)
                shiftLeftDisplay(5)
                if config.getOptValue("verif") and solver =='OPTW':
                    if not self.sequenceFaisable(constellation,modeleDeTransition):
                        die("CCA : "+str(self.identifiant)+ " " +str(self.sequence)+" non faisable. Retard total : "+ str(self.retardComposante(constellation,modeleDeTransition,print_info=True)))
                if solver in ['OPTW','OPTWGroups']:
                    return True,res # sequence forcement faisable, renvoyer les requetes satisfaites
                else:
                    return self.sequenceFaisable(constellation,modeleDeTransition)
            except InfeasibleSolutionfromSolverException:
                if config.verifMode():
                    print("Erreur solver : solution infaisable.")
                self.sequence = copy_sequence
                return False,None
        else:
            shiftLeftDisplay(5)
            self.use_solver = False
            if solver=='LKH':
                return True
            else:
                req = []
                for g in groups:
                    if len(groups[g])>0:
                        req.append(constellation.getRequeteActivite(groups[g][0]))
                return True,req
            
    def ajouterMode(self,constellation,r,m,modeleDeTransition,cheap=False):
        ajouts = []
        # ajout des activites            
        for s in constellation.getRequete(r).getMode(m).getActivites():
            for o in constellation.getRequete(r).getMode(m).getActivites()[s]:
                if o in self.activites and o not in self.sequence:
                    if cheap:
                        self.insererActiviteCheap(constellation,o,modeleDeTransition,test_faisab=False) # ajouter a un endroit intelligent
                    else:
                        self.insererActivite(constellation,o,modeleDeTransition,test_faisab=False) # ajouter a un endroit intelligent
                ajouts.append(o)
        assert(len(self.sequence) > 0)
        return ajouts
    
    def ajouterListeModes(self,constellation,modes_tries,modeleDeTransition):
        ajouts = []
        #if not config.isTimeDependentModeOn():
        if not self.planEarliestAJour:
                assert(not self.planLatestAJour)
                faisable = self.MAJPlansCritiques(constellation,modeleDeTransition)
        assert(self.planLatestAJour and self.planEarliestAJour)
        
        for (r,m) in modes_tries:
            for a in modes_tries[(r,m)]:
                if config.isTimeDependentModeOn():
                    self.insererActivite(constellation,a,modeleDeTransition,test_faisab=False) # ajouter a un endroit intelligent
                else:
                    if faisable:
                        self.insererActivitePlansCritiques(constellation,a,modeleDeTransition)
                        faisable = self.MAJPlansCritiques(constellation,modeleDeTransition)
                    else:
                        self.insererActiviteCheap(constellation,a,modeleDeTransition,test_faisab=False) # ajouter ou la transition est minimale
                ajouts.append(a)
        assert(len(self.sequence) > 0)
        return ajouts        
        
    def ajouterActivites(self,constellation,activites,modeleDeTransition,cheap=False):
        ajouts = []
        # ajout des activites            
        for o in sorted(activites):
            if o not in self.sequence:
                if cheap:
                    self.insererActiviteCheap(constellation,o,modeleDeTransition,test_faisab=False) # ajouter a un endroit intelligent
                else:
                    self.insererActivite(constellation,o,modeleDeTransition,test_faisab=False) # ajouter a un endroit intelligent
                ajouts.append(o)
        assert(len(self.sequence) >0)
        return ajouts
    
    def retardSequence(self,constellation,sequence,modeleDeTransition,depart=None,print_info=False):       
        s = self.satellite
        retard = 0
        if depart is None and len(sequence)>0:
            depart = constellation.getSatellite(s).getActivite(sequence[0]).getDebut()
            t = depart
            for i,a in enumerate(sequence):
                t = round(max(t,constellation.getSatellite(s).getActivite(a).getDebut()),config.glob.digits)
                duree = constellation.getSatellite(s).getActivite(a).getDuree()
                retard += max(0,-constellation.getSatellite(s).getActivite(a).getFin()+t+duree)# retard tronqué
                if print_info:
                    print(t,constellation.getSatellite(s).getActivite(a),retard)                
                if(i<len(sequence)-1):
                    if modeleDeTransition.estTimeDependent():
                        t += round(duree + modeleDeTransition.getTransition(a,sequence[i+1], t+duree),config.glob.digits)
                    else:
                        t += round(duree + constellation.satellites[s].getTransition(a,sequence[i+1],modeleDeTransition),config.glob.digits)
        return retard
 
    def retardAlgebriqueSequence(self,constellation,sequence,modeleDeTransition,depart):
        s = self.satellite
        retard = 0
        t = depart
        for i,a in enumerate(sequence):
            t += constellation.getSatellite(s).getActivite(a).getDuree()
            retard += t-constellation.satellites[s].getActivite(a).getFin()
            if(i<len(sequence)-1):
                if modeleDeTransition.estTimeDependent():
                    t += modeleDeTransition.getTransition(a,sequence[i+1], t)
                else:
                    t += constellation.satellites[s].getTransition(a,sequence[i+1],modeleDeTransition)
                t = round(max(t,constellation.satellites[s].getActivite(a).getDebut()),config.glob.digits)
        return retard

    # ( retard , retard algébrique )
    def duoRetardSequence(self,constellation,sequence,modeleDeTransition,depart=None):
        s = self.satellite
        retard = 0
        retard_algebrique = 0
        if depart is None and len(sequence)>0:
            depart = constellation.getSatellite(s).getActivite(sequence[0]).getDebut()
        t = depart
        for i,a in enumerate(sequence):
            t += constellation.getSatellite(s).getActivite(a).getDuree()
            retard += max(0,t-constellation.satellites[s].getActivite(a).getFin())# retard tronqué
            retard_algebrique += max(0,t-constellation.satellites[s].getActivite(a).getFin())# retard tronqué
            if(i<len(sequence)-1):
                if modeleDeTransition.estTimeDependent():
                    t += modeleDeTransition.getTransition(a,sequence[i+1], t)
                else:
                    t += constellation.satellites[s].getTransition(a,sequence[i+1],modeleDeTransition)
                t = max(t,constellation.satellites[s].getActivite(a).getDebut())
        return retard,retard_algebrique
    
    """
        ===============================================
                    INSERTION DES ACTIVITES
        ===============================================
    """
    # solution = proposition d'ordre pour la solution, pas la solution courante
    # info supplementaire pour otpimiser le calcul du retard : Retards a chaque position + dates early
    def partitionConflitRetard(self,constellation,p,modeleDeTransition):
        assert(not modeleDeTransition.estTimeDependent())
        s = self.satellite
        a,b = constellation.getSatellite(s).getActivite(p).getDebut(),constellation.getSatellite(s).getActivite(p).getFin()
        L,C,R = [],[],[]
        Tearly = []
        Retards = []
        for i,pp in enumerate(self.sequence):
            # calcul de la duree
            duree = constellation.getSatellite(s).getActivite(p).getDuree()
            # determination des 3 listes
            if i==0:
                Tearly.append(constellation.getSatellite(s).getActivite(pp).getDebut())
                Retards.append((0,0))
            else:
                prev = self.sequence[i-1]
                duree = constellation.getSatellite(s).getActivite(pp).getDuree()
                tearly = Tearly[-1]+constellation.getSatellite(s).getTransition(prev,pp,modeleDeTransition)
                tearly += duree
                tearly = max(tearly,constellation.getSatellite(s).getActivite(pp).getDebut())
                retard_alg = constellation.getSatellite(s).getActivite(pp).getFin()-duree-tearly
                retard = max(retard_alg,0)
                Tearly.append(tearly)
                Retards.append((retard+Retards[-1][0],retard_alg+Retards[-1][1]))
            if(constellation.getSatellite(s).getActivite(pp).getFin()+constellation.getSatellite(s).getTransition(pp,p,modeleDeTransition)<a):
                L.append(pp)
            elif b + constellation.getSatellite(s).getTransition(p,pp,modeleDeTransition) < constellation.getSatellite(s).getActivite(pp).getDebut():
                R.append(pp)
            else:
                C.append(pp)
        return L,C,R,Tearly,Retards

    def conflit(self,a1,a2,constellation,modeleDeTransition):
        s = self.satellite
        debut1,fin1 = constellation.getSatellite(s).getActivite(a1).getDebut(),constellation.getSatellite(s).getActivite(a2).getFin()
        debut2,fin2 = constellation.getSatellite(s).getActivite(a2).getDebut(),constellation.getSatellite(s).getActivite(a2).getFin()
        test1 = fin1+modeleDeTransition.getMajorantCoupleActivite(a1,a2)<debut2
        test2 = fin2+modeleDeTransition.getMajorantCoupleActivite(a2,a1)<debut1
        return not test1 and not test2
    
    # solution = proposition d'ordre pour la solution, pas la solution courante
    def partitionConflit(self,constellation,p,modeleDeTransition):
        s = self.satellite
        a,b = constellation.getSatellite(s).getActivite(p).getDebut(),constellation.getSatellite(s).getActivite(p).getFin()
        L,C,R = [],[],[]
        for pp in self.sequence:
            # calcul de la duree
            duree = constellation.getSatellite(s).getActivite(p).getDuree()
            if modeleDeTransition.estTimeDependent():
                transition_left = modeleDeTransition.getMajorantCoupleActivite(pp,p)
                transition_right = modeleDeTransition.getMajorantCoupleActivite(p,pp)
            else:
                transition_right = constellation.getSatellite(s).getTransition(pp,p,modeleDeTransition)
                transition_left = transition_right
            # determination des 3 listes
            if(constellation.getSatellite(s).getActivite(pp).getFin()+transition_left<a):
                L.append(pp)
            elif b + transition_right < constellation.getSatellite(s).getActivite(pp).getDebut():
                R.append(pp)
            else:
                C.append(pp)
        return L,C,R
    
    # not self.modeleDeTransition.estTimeDependent()
    
    def retardComposante(self,constellation,modeleDeTransition,force=False,print_info=False):
        if not self.planEarliestAJour or force:
            if print_info:
                printColor("MAJ des plans critiques.")
            self.MAJPlansCritiques(constellation,modeleDeTransition,force=force)
        else:
            if print_info:
                printColor("plans critiques à jour.")
        return self.retard

    def coutTransitionComposante(self,constellation,modeleDeTransition):
        if not self.couts_transition_a_jour:
            self.MAJPlansCritiques(constellation,modeleDeTransition)
        return self.couts_transition
    
    def sequenceFaisable(self,constellation,modeleDeTransition,print_retard=False):
        if self.sequence is None:
            return False
        if(len(self.sequence))==0:
            return True
        s = self.satellite
        t = constellation.getSatellite(s).getActivite(self.sequence[0]).getDebut()
        for i,p in enumerate(self.sequence):
            if print_retard:
                print(t,constellation.getSatellite(s).getActivite(p))
            duree = constellation.getSatellite(s).getActivite(p).getDuree()
            if t+duree-constellation.getSatellite(s).getActivite(p).getFin()>10**(-config.glob.digits):
                if print_retard:
                    printColor(t,constellation.getSatellite(s).getActivite(p),c='r')
                return False
            if i<len(self.sequence)-1:
                if modeleDeTransition.estTimeDependent():
                    transition = constellation.getSatellite(s).getTransitionTimeDependent(p,self.sequence[i+1],t+duree,modeleDeTransition)
                else:
                    transition = constellation.getSatellite(s).getTransition(p,self.sequence[i+1],modeleDeTransition)
                if print_retard:
                    print(p,self.sequence[i+1],t+duree,transition,t+duree>constellation.getSatellite(s).getActivite(p).getFin()>10**(-config.glob.digits))
                t = round(max(t + transition + duree,constellation.getSatellite(s).getActivite(self.sequence[i+1]).getDebut()),config.glob.digits)
        return True
    
    def extensionConflit(self,L,C,R,constellation,modeleDeTransition):
        if len(L)>0:
            depart_L = self.planEarliest(constellation,L,modeleDeTransition)[-1][1]
            depart = lambda c : depart_L
        else: # c_extend toujours != []
            depart = lambda c : constellation.getSatellite(self.satellite).getActivite(c[0]).getDebut()
        if L==[]:
            C_left = lambda c : c
        else:
            C_left = lambda c : [L[-1]] + c
        if R==[]:
            C_right = lambda c : c
        else:
            C_right = lambda c : c + [R[0]]
        C_extend = lambda c : C_left(C_right(c))
        return depart,C_extend
   
    # ( retard , retard algébrique )
    def duoRetardInsertionTemporaire(self,constellation,sequence,modeleDeTransition,depart=None):
        s = self.satellite
        retard = 0
        retard_algebrique = 0
        if depart is None and len(sequence)>0:
            depart = constellation.getSatellite(s).getActivite(sequence[0]).getDebut()
        t = depart
        for i,a in enumerate(sequence):
            t += constellation.getSatellite(s).getActivite(a).getDuree()
            retard += round(max(0,t-constellation.satellites[s].getActivite(a).getFin()),config.glob.digits)# retard tronqué
            retard_algebrique +=round( max(0,t-constellation.satellites[s].getActivite(a).getFin()),config.glob.digits)# retard tronqué
            if(i<len(sequence)-1):
                t += constellation.satellites[s].getTransition(a,sequence[i+1],modeleDeTransition)
                t = round(max(t,constellation.satellites[s].getActivite(a).getDebut()),config.glob.digits)
        return retard,retard_algebrique
    
    def critereRetardOptimise(self,constellation,p,i,L,C,R,Tearly,Retards):
        delta_retard = (0,0)
        rang_sequence = len(L)+i
        duree = constellation.getSatellite(self.satellite).getActivite(p).getDuree()
        debut = constellation.getSatellite(self.satellite).getActivite(p).getDebut()
        fin = constellation.getSatellite(self.satellite).getActivite(p).getFin()
        # date early de p
        if rang_sequence==0:
            date_early = constellation.getSatellite(self.satellite).getActivite(p).getDebut()
        else:
            prev = self.sequence[rang_sequence-1]
            date_early = Tearly[rang_sequence-1]
            date_early += constellation.getSatellite(self.satellite).getActivite(prev).getDuree()
            date_early = max(date_early,debut) 
        # retard de p
        retard_alg = fin - duree - date_early
        retard = (round(max(0,retard_alg),config.glob.digits),round(retard_alg,config.glob.digits))
        # date early du successeur
        if rang_sequence==len(self.sequence)-1:
            return retard
        else:
            delta = retard
            rang_suivant = rang_sequence+1
            t_k = date_early
            for k in range(rang_sequence,len(self.sequence)-1):
                ak = constellation.getSatellite(self.satellite).getActivite(self.sequence[k])
                if k==0:
                    t_k = round(max(date_early + duree,ak.getDebut(),0),config.glob.digits)
                else:
                    prev = constellation.getSatellite(self.satellite).getActivite(self.sequence[k-1])
                    t_k = round(max(t_k + prev.getDuree(),ak.getDebut()),config.glob.digits)
                retard_k_alg = ak.getFin()-t_k-ak.getDuree()
                retard_k = (config.glob.digits(max(retard_alg,0),config.glob.digits),round(retard_alg,config.glob.digits))
                delta = (delta[0]+retard[0],delta[1]+retard[1])
                if delta==Retards[k]:
                    return delta
            return delta
    
    def insererActivitePlansCritiques(self,constellation,p,modeleDeTransition):
        #assert(not modeleDeTransition.estTimeDependent())
        def calculerTearly(constellation,p,i,modeleDeTransition):
            s = self.satellite
            if i>0:
                prev = self.planEst[i-1]
                Tearly = prev[1]
                duree = constellation.getSatellite(s).getActivite(prev[0]).getDuree()
                Tearly += duree
                start_previous_maneuver = Tearly
                if not modeleDeTransition.estTimeDependent():
                    Tearly += constellation.getSatellite(s).getTransition(prev[0],p,modeleDeTransition)
                else:
                    Tearly += constellation.getSatellite(s).getTransitionTimeDependent(prev[0],p,prev[1]+duree,modeleDeTransition)
                Tearly = max(Tearly,activite.getDebut())
            else:
                Tearly = activite.getDebut()
                start_previous_maneuver = None
            return Tearly,start_previous_maneuver
        
        def calculerLatest(p,i,modeleDeTransition):
            if i<len(self.planLts):
                next_a = self.planLts[i]
                latest_start_next_activity = next_a[1]
                if modeleDeTransition.estTimeDependent():
                    latest_maneuver_start,trans_to_next_duration = calculerTransitionToNextTimeDep(constellation,s,p,next_a[0],next_a[1]-10e-4,modeleDeTransition)
                    self.verifierCalculTransition(modeleDeTransition,p,next_a[0],latest_maneuver_start,trans_to_next_duration,next_a[1])
                
                else:
                    trans_to_next_duration =  float(constellation.getSatellite(s).getTransition(p,next_a[0],modeleDeTransition))
                    latest_maneuver_start = latest_start_next_activity - trans_to_next_duration
            else:
                latest_maneuver_start = np.Inf
                trans_to_next_duration = 0
                latest_start_next_activity = np.Inf
            return latest_maneuver_start,latest_start_next_activity,trans_to_next_duration
        # a appeler avant l'insertion
                
        def calculerDeltaTransitionNouvellesInsertion(sequence,constellation,p,position_insertion,modeleDeTransition,start_previous_maneuver):
            if not modeleDeTransition.estTimeDependent():
                # ajout de previous(p) -> p
                ajout_transition_gauche = 0
                if position_insertion>0:
                    ajout_transition_gauche += constellation.getSatellite(s).getTransition(sequence[position_insertion-1],p,modeleDeTransition)
                # ajout de p -> next(p)
                ajout_transition_droite = 0
                if position_insertion<len(sequence):
                    ajout_transition_droite += constellation.getSatellite(s).getTransition(p,sequence[position_insertion],modeleDeTransition)
                # perte de la transition previous(p) -> next(p) 
                perte_transition = 0
                if position_insertion > 0 and position_insertion<len(sequence):
                    perte_transition += constellation.getSatellite(s).getTransition(sequence[position_insertion-1],sequence[position_insertion],modeleDeTransition)
                return ajout_transition_gauche + ajout_transition_droite - perte_transition
            else:
                # ajout de previous(p) -> p
                ajout_transition_gauche = 0
                if position_insertion>0:
                    ajout_transition_gauche += constellation.getSatellite(s).getTransitionTimeDependent(sequence[position_insertion-1],p,start_previous_maneuver,modeleDeTransition)
                # ajout de p -> next(p)
                ajout_transition_droite = 0
                if position_insertion<len(sequence):
                    if start_previous_maneuver is not None:
                        start = max(start_previous_maneuver+ajout_transition_gauche,constellation.getSatellite(s).getActivite(p).getDebut())
                        start_maneuver_to_next = start + constellation.getSatellite(s).getActivite(p).getDuree()
                        ajout_transition_droite += constellation.getSatellite(s).getTransitionTimeDependent(p,sequence[position_insertion],start_maneuver_to_next,modeleDeTransition)
                # perte de la transition previous(p) -> next(p) 
                perte_transition = 0
                if position_insertion > 0 and position_insertion<len(sequence):
                    perte_transition += constellation.getSatellite(s).getTransitionTimeDependent(sequence[position_insertion-1],sequence[position_insertion],start_previous_maneuver,modeleDeTransition)
                return ajout_transition_gauche + ajout_transition_droite - perte_transition
                
        def calculerPositionInsertionOptimale(constellation,p,iLeft,iRight,modeleDeTransition):
            assert(self.plansAJour())
            
            position_insertion_optimale,critere_opti,time_opti = -1,(0,-np.Inf),np.Inf  # faisable,flexibilite
            for i in range(iLeft,iRight+1):
                try:
                    t12 = time()
                    assert(self.plansAJour())
                    assert(len(self.sequence)==len(self.planLts))
                    Tearly,start_previous_maneuver = calculerTearly(constellation,p,i,modeleDeTransition)
                    t13 = time()
                    latest_maneuver_start,latest_start_next_activity,trans_to_next_duration = calculerLatest(p,i,modeleDeTransition)
                    t14 = time()
                    faisable = Tearly + activite.getDuree() <= activite.getFin() and Tearly + activite.getDuree() <= latest_maneuver_start

                    if config.verifMode():
                        if i<len(self.sequence)-1:
                            assert(constellation.getSatellite(s).getTransitionTimeDependent(p,self.sequence[i],latest_maneuver_start,modeleDeTransition) <= self.planLts[i][1])
                        
                    if config.verifMode() and faisable:
                        if position_insertion_optimale<len(self.sequence):
                            self.sequence = self.sequence[:i]+[p]+self.sequence[i:]
                        else:
                            self.sequence.append(p)
                        if not(self.sequenceFaisable(constellation,modeleDeTransition)):
                            for a in self.sequence:
                                print(constellation.getSatellite(s).getActivite(a))
                            print("----------------------------------")
                            
                            print("insertion de ",p," a la positon ",i)
                            print("start_tau_latest",latest_maneuver_start)
                            print("latest_start_next",latest_start_next_activity)
                            print("duration_time",trans_to_next_duration)
                            print("early",self.planEst)
                            print("late",self.planLts)
                            print("avant",self.planEarliest(constellation,self.sequence,modeleDeTransition,add_transition=True))
                            self.sequence.remove(p)
                            print("apres",self.planEarliest(constellation,self.sequence,modeleDeTransition,add_transition=True))
                            die(self.sequence,self.retardSequence(constellation,self.sequence,modeleDeTransition),self.sequenceFaisable(constellation,modeleDeTransition))
                        self.sequence.remove(p)
                    
                    #Tearly+activite.getDuree()<=activite.getFin() and Tearly+activite.getDuree()+trans_to_next<=latest_next
                    flexibilite = latest_maneuver_start-(Tearly+activite.getDuree())
                    critere = (faisable,flexibilite)
                    if critere>critere_opti:
                        position_insertion_optimale,critere_opti,time_opti = i,critere,start_previous_maneuver
                    t15 = time()
                    printColor("calcul de Tearly",round(t13-t12,7),"calcul de Tlatest",round(t14-t13,7),"calcul de la flexibilité",round(t15-t14,7),c='m')
                except (InfeasibleSolutionException,NumericalError) as error:
                    pass
            return position_insertion_optimale,critere_opti,time_opti
      
        if config.verifMode() :
            if not(self.sequenceFaisable(constellation,modeleDeTransition)):
                print(self.sequence,position_insertion_optimale,critere_opti,start_previous_maneuver)
                assert(False)        

        printOpen("insertion de l'activité",p,c='m')
        s = self.satellite
        t0 = time() 
        activite = constellation.getSatellite(s).getActivite(p)
        if self.sequence == []:
            self.sequence = [p]
            faisable = True
            self.retard = 0
            self.couts_transition = 0
            self.planEarliestAJour = True
            self.planLatestAJour = True
            self.position_dernier_insere = 0
        else:
            self.MAJPlansCritiques(constellation,modeleDeTransition)
            if not self.plansAJour():
                raise OutdatedCriticalPlans                   
            
            iLeft,iRight = self.conflitPlan(constellation,activite,modeleDeTransition)
            if config.getOptValue("verif"):
                self.checkCoutTransition(constellation,modeleDeTransition)
                       
            t1 = time()
            position_insertion_optimale,critere_opti,start_previous_maneuver = calculerPositionInsertionOptimale(constellation,p,iLeft,iRight,modeleDeTransition)
            t2 = time()
            self.couts_transition += round(calculerDeltaTransitionNouvellesInsertion(self.sequence,constellation,p,position_insertion_optimale,modeleDeTransition,start_previous_maneuver),config.glob.digits)
            if position_insertion_optimale<len(self.sequence):
                self.sequence = self.sequence[:position_insertion_optimale]+[p]+self.sequence[position_insertion_optimale:]
            else:
                self.sequence.append(p)    
            faisable = critere_opti[0]           
            t3 = time()
            printColor("Bornes d'insertions:"+str(round(t1-t0,7))+"/pos. optimale:"+str(round(t2-t1,7))+"/MAJ séquence:"+str(round(t3-t2,7)),c='m') 
            
            if config.getOptValue("verif"):
                self.checkCoutTransition(constellation,modeleDeTransition)
            self.position_dernier_insere = position_insertion_optimale
        self.notifierPlanAJour(False,False)
        self.couts_transition_a_jour = True
        self.debut = min(self.debut,activite.getDebut())
        self.fin = max(self.fin,activite.getFin()+modeleDeTransition.getMajorantDureeTransition())
        printClose()
        
        if config.verifMode() and faisable:
            if not(self.sequenceFaisable(constellation,modeleDeTransition)):
                self.sequenceFaisable(constellation,modeleDeTransition)
                print(constellation.getSatellite(self.satellite).getActivite(p))
                print(self.planEst)
                print(self.planLts)
                for i,a in enumerate(self.sequence):
                    print(constellation.getSatellite(s).getActivite(a))
                    #print(constellaiton.getTransitionTimeDependent(self.sequence[i+1],a,start_maneuvering,modeleDeTransition))
                die(p,self.sequence,position_insertion_optimale,critere_opti,start_previous_maneuver)
                
                assert(False)
        
        if faisable:
            try:
                self.MAJPlansCritiques(constellation,modeleDeTransition,force=True)
            except NumericalError:
                assert(modeleDeTransition.estTimeDependent())
        
        return faisable
    
    def verifierCalculTransition(self,modeleDeTransition,a1,a2,date_depart,longueur_supposee,latest_a2):
        if config.verifMode():
            computed_maneuver = modeleDeTransition.getTransition(a1,a2,date_depart)
            if not computed_maneuver==longueur_supposee:
                die(computed_maneuver,longueur_supposee)
            if date_depart+computed_maneuver>latest_a2:
                die(date_depart+computed_maneuver,latest_a2)
                    
    def checkCoutTransition(self,constellation,modeleDeTransition):
        transition_apres_insertion = self.couts_transition
        retard_apres = self.retard
        self.MAJPlansCritiques(constellation,modeleDeTransition)
        if not abs(self.couts_transition-transition_apres_insertion)<config.glob.getEchelle():
            print(self.couts_transition-transition_apres_insertion)
            print(self.sequence,p,transition_apres_insertion,self.couts_transition)
        assert(abs(self.couts_transition-transition_apres_insertion)<config.glob.getEchelle())        
    
    def insererActivite(self,constellation,p,modeleDeTransition,test_faisab=True):
        len_avant = len(self.sequence)
        s = self.satellite
        L,C,R = self.partitionConflit(constellation,p,modeleDeTransition)
        # c_etendue : c + extremites = sous-sequence a optimiser
        if len(L)>0:
            depart = lambda c : self.planEarliest(constellation,L,modeleDeTransition)[-1][1]
        else: # c_extend toujours != []
            depart = lambda c : constellation.getSatellite(s).getActivite(c[0]).getDebut()
        if L==[]:
            C_left = lambda c : c
        else:
            C_left = lambda c : [L[-1]] + c
        if R==[]:
            C_right = lambda c : c
        else:
            C_right = lambda c : c + [R[0]]
        C_extend = lambda c : C_left(C_right(c))
        # critere 1 = retard. Minimiser = mesurer la faisabilite (faisable si critere = 0)
        # critere 2 : retard algebrique. Minimiser = augmenter la flexibilite
        critere_opti,sequence_opti = None,None
        calculer_critere = lambda sequence : (self.retardSequence(constellation,sequence,modeleDeTransition,depart=depart(sequence)),self.retardAlgebriqueSequence(constellation,sequence,modeleDeTransition,depart(sequence)))
        for i in range(len(C)+1):
            ordre = C[0:i]+[p]+C[i:]
            if(critere_opti is None):
                critere_opti,sequence_opti = calculer_critere(C_extend(ordre)),ordre
            else:
                critere = calculer_critere(C_extend(ordre))
                if(critere<critere_opti):
                    critere_opti,sequence_opti = critere,ordre
        self.sequence = L+sequence_opti+R
        self.notifierPlanAJour(False,False)
        assert(len(self.sequence)==len_avant+1) 
        if test_faisab:
            return self.sequenceFaisable(constellation,modeleDeTransition)
        #self.planEarliestAJour = False
        
    def insererActiviteCheap(self,constellation,p,modeleDeTransition,test_faisab=True):
        assert(not modeleDeTransition.estTimeDependent())
        def evalPos(constellation,i,a):
            s = self.satellite
            if i == 0:
                trans_left = 0
            else:
                trans_left = constellation.getSatellite(s).getTransition(self.sequence[i-1],p,modeleDeTransition)
            if i == len(self.sequence):
                trans_right = 0
            else:
                trans_right = constellation.getSatellite(s).getTransition(p,self.sequence[i],modeleDeTransition)
            return trans_left + trans_right
            
        len_avant = len(self.sequence)
        if len_avant == 0:
            self.sequence = [p]
            if test_faisab:
                return True
        else:
            best_pos = None
            best_critere = None # transition
            for i in range(len(self.sequence)+1):
                critere = evalPos(constellation,i,p)
                if best_critere is None or critere > best_critere:
                    best_critere = critere
                    best_pos = i
                    
            if best_pos ==0:
                self.sequence = [p] + self.sequence[best_pos:]
            elif best_pos != len(self.sequence):
                self.sequence = self.sequence[0:best_pos] + [p] + self.sequence[best_pos:]
            else:
                self.sequence = self.sequence + [p]        
            self.notifierPlanAJour(False,False)
            assert(len(self.sequence)==len_avant+1) 
            if test_faisab:
                return self.sequenceFaisable(constellation,modeleDeTransition)
        #self.planEarliestAJour = False
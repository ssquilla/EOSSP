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
def inverseLinearManeuver(tnext,pointsLeft,pointsRight,duration):
    tleft,L = pointsLeft
    tright,R = pointsRight
    amplitude = tright-tleft
    duration_delta = R-L
    t = ((tnext)*amplitude+tleft*R-tright*L)/(duration_delta+amplitude)
    return t
        
# left = right si controlPoints[left] vaut exactement la date recherchee
# left + 1 = right sinon. Alors la bascule est dans ce segment
def searchChangeSegment(controlPoints,date_next_a,duration):
    left = 0
    right = len(controlPoints)-1
        
    startManeuver,maneuver = controlPoints[left]
    while(left<len(controlPoints)-1 and startManeuver+maneuver<date_next_a):
        left += 1
        startManeuver,maneuver = controlPoints[left]
    left -= 1    
        
    startManeuver,maneuver = controlPoints[right]
    while(right > 0 and startManeuver+maneuver>date_next_a):
        right -= 1        
        startManeuver,maneuver = controlPoints[right]
    right += 1
    return left,right    

class InfeasibleSolutionException(Exception):
    def __init__(self,msg=""):
        super().__init__(msg)

class OutdatedCriticalPlans(Exception):
    pass

EPSILON = 1e-4

def findManeuver(constellation,s,p,next_a,date_next_a,transitionModel,ensureFeasible):
    duration = constellation.getSatellite(s).getActivity(p).getDuration()
    controlPoints = transitionModel.getPointsDeControle(p,next_a)
    endDate = constellation.getSatellite(s).getActivity(p).getEndDate()
            
    """ ce cas n'aurait pas du arriver, mais il y a un soucis dans
        les donnees. Il peut arriver qu'il manque un points ..."""
    if controlPoints[-1][0]<endDate:
        controlPoints.append((endDate,controlPoints[-1][1]))
        
    # cas ou la transition est toujours feasible
    if controlPoints == [(0,0)]:
        return endDate,0   
    
    # cas ou la maneuver au plus tard arrive avant date_next_a
    endWindowManeuver = constellation.getSatellite(s).getActivity(p).getEndDate()
    latestManeuver = constellation.getSatellite(s).getTransitionDurationTimeDependent(p,next_a,endWindowManeuver,transitionModel)
    if latestManeuver+endWindowManeuver<=date_next_a:
        assert(latestManeuver<transitionModel.getTransititionDurationUpperBound())
        return endWindowManeuver,latestManeuver
    endWindowManeuver = controlPoints[-1][0]
    
    # cas ou partir au debut de la fenetre ne suffit pas : elager le if suivant
    if controlPoints[0][0]+controlPoints[0][1]>=date_next_a:
        msg = "transition "+str(p)+"->"+str(next_a)+" tau="+str(date_next_a)+"\n"
        msg += "Start_p + tau(p) >= date_next_a"
        raise InfeasibleSolutionException(msg)
        
    # cas ou partir au debut de la fenetre + duration (au plus tot) ne suffit pas
    startEarliestManeuver = controlPoints[0][0]+duration
    earliestManeuver = constellation.getSatellite(s).getTransitionDurationTimeDependent(p,next_a,startEarliestManeuver,transitionModel)
    if earliestManeuver+startEarliestManeuver>date_next_a:
        assert(earliestManeuver<transitionModel.getTransititionDurationUpperBound())
        if not ensureFeasible:
            msg = "!earliestManeuver+startEarliestManeuver>=date_next_a\n"
            msg += "transition "+str(p)+"->"+str(next_a)+"\n"
            msg += str("earliestManeuver="+str(earliestManeuver)+"\n")
            msg += str("startEarliestManeuver="+str(startEarliestManeuver)+"\n")
            msg += str("date_next_a="+str(date_next_a)+"\n")
            msg += str("earliestManeuver+startEarliestManeuver="+str(earliestManeuver+startEarliestManeuver)+"\n")
            raise InfeasibleSolutionException(msg)
        else:
            raise NumericalError()
            
    # ici on cherche le segment sur lequel on va basculer dans le cas infeasible
    left,right = searchChangeSegment(controlPoints,date_next_a,duration)
    if left==right:
        return controlPoints[left]
    else:
        pointsLeft = controlPoints[left]
        pointsRight = controlPoints[right]
        t = inverseLinearManeuver(date_next_a,pointsLeft,pointsRight,duration)
        maneuver = constellation.getSatellite(s).getTransitionDurationTimeDependent(p,next_a,t,transitionModel)
        if not(t+maneuver<=date_next_a):
            msg = "Maneuver "+str(p)+"->"+str(next_a)
            msg += ". left="+str(pointsLeft)+",right="+str(pointsRight)+ "points ="+str(controlPoints)
            msg += " start="+str(t)+" duration="+str(maneuver)
            assertLessThan(t+maneuver,date_next_a+EPSILON,msg)
            raise NumericalError()
        return t,maneuver

def computeTransitionToNextTimeDepPoint(constellation,s,p,next_a,date_next_a,transitionModel,ensureFeasible=False):
    assert(transitionModel.isTimeDependent())
    startManeuver,duration = findManeuver(constellation,s,p,next_a,date_next_a,transitionModel,ensureFeasible)
    computerManeuver = constellation.getSatellite(s).getTransitionDurationTimeDependent(p,next_a,startManeuver,transitionModel)
    if not(startManeuver+computerManeuver<=date_next_a):
        raise NumericalError()
    if not(startManeuver+duration<=date_next_a):
        assert(startManeuver+duration<=date_next_a+EPSILON)
        raise NumericalError()
    assert(duration==computerManeuver)
    
    return  (startManeuver,duration)
                    

class SolCCA:  
    def __init__(self,identifier,satellite):
        self.identifier = identifier
        self.satellite = satellite
        self.startDate = np.Inf
        self.endDate = -np.Inf
        #self.activities = []
        self.sequence = []
        self.LKHTime = []
        self.planEst = []
        self.planLts = []
        self.notifyPlansUpToDate(True,True)
        self.transitionCosts = 0
        self.transitionCostsUpToDate = True
        self.positionLastInserted = None
        self.tardiness = 0

    def arePlansUpToDate(self):
        return self.isEarliestPlanUpToDate and self.isLatestPlanUpToDate
   
    def rename(self,new):
        self.identifier = new
        
    def getCopy(self):
        copy = SolCCA(self.identifier,self.satellite)
        copy.startDate = self.startDate
        copy.endDate = self.endDate
        copy.planEst = self.planEst.copy()
        copy.planLts = self.planLts.copy()
        copy.notifyPlansUpToDate(self.isEarliestPlanUpToDate,self.isLatestPlanUpToDate)
        copy.sequence = self.sequence.copy()
        copy.tardiness = self.tardiness
        copy.transitionCosts = self.transitionCosts
        copy.transitionCostsUpToDate = self.transitionCostsUpToDate
        copy.positionLastInserted = self.positionLastInserted
        return copy
    
    def formatInformation(self):
        msg = "CCA "+str(self.identifier)
        msg += "(sat. "+str(self.satellite) +")"
        msg += " € ["+ str(self.startDate)+"," +str(self.endDate) +"]"
        return msg
        
    def reset(self):
        self.startDate = np.Inf
        self.endDate = -np.Inf
        self.sequence = []
        self.planEst = []
        self.planLts = []
        self.notifyPlansUpToDate(True,True)
        self.transitionCosts = 0
        self.transitionCostsUpToDate = True
        self.positionLastInserted = None
        
    def loadScore(self,method='density'):
        if method == 'density':
            return len(self.sequence)/(self.endDate-self.startDate)
        elif method == 'cardinal':
            return len(self.sequence)
        else:
            raise ValueError('Unknown scoring method.')
    
    def loadScorePlusActivityAddingCost(self,nActivities,method='density'):
        if method=='density':
            return (len(self.getSequence())+nActivities)/(self.endDate-self.startDate)
        elif method == 'cardinal':
            return len(self.sequence)+nActivities
        else:
            raise ValueError('Unknown scoring method.')
        
    def rename(self,new):
        self.identifier = new
        
    def __str__(self):
        return "SolCCA(" + str(self.identifier) +":"+str(self.sequence)+")"
        
    def notifyPlansUpToDate(self,early,late):
        self.isEarliestPlanUpToDate = early
        self.isLatestPlanUpToDate = late
        
    def removeActivityList(self,constellation,activities,transitionModel):
        self.sequence = [a for a in self.sequence if a not in activities]
        self.planEst = [at for at in self.planEst if at[0] not in activities]
        self.planLts = [at for at in self.planLts if at[0] not in activities]
        self.notifyPlansUpToDate(False,False)
        self.transitionCostsUpToDate = False
        self.updateStartEndDates(constellation,transitionModel)
            
    def removeActivity(self,constellation,a,transitionModel):
        self.sequence.remove(a)
        self.notifyPlansUpToDate(False,False)
        self.transitionCostsUpToDate = False
        self.updateStartEndDates(constellation,transitionModel)
        
    def isEmpty(self):
        return len(self.activities)==0
    
    def cancelMode(self,constellation,r,m,transitionModel):
        if self.sequence is not None:
            for (s,o) in constellation.getRequest(r).getMode(m).getPairs():
                if o in self.sequence:
                    self.sequence.remove(o)
        self.notifyPlansUpToDate(False,False)
        self.transitionCostsUpToDate = False
        self.updateStartEndDates(constellation,transitionModel)
    
    def updateStartEndDates(self,constellation,transitionModel):
        tau_max = transitionModel.getTransititionDurationUpperBound()
        try:
            self.startDate = min([constellation.getSatellite(self.satellite).getActivity(aa).getStartDate() for aa in self.sequence])
            self.endDate = max([constellation.getSatellite(self.satellite).getActivity(aa).getEndDate()+tau_max for aa in self.sequence])
        except ValueError:
            self.startDate = np.Inf
            self.endDate = - np.Inf
            
    def setSequence(self,constellation,seq,transitionModel):
        self.sequence = deepcopy(seq)
        self.updateStartEndDates(constellation,transitionModel)
        self.notifyPlansUpToDate(False,False)
        self.transitionCostsUpToDate = False
        
    def setSequenceIndex(self,i_start,i_end,seq,transitionModel):
        #for a in seq:
        #    assert(a in self.activities)
        if i_end<len(self.sequence):
            self.sequence[i_start:] = seq + self.sequence[i_end+1:]
        else:
            self.sequence[i_start:] = seq
        self.updateStartEndDates(constellation,transitionModel)
        self.notifyPlansUpToDate(False,False)
        self.transitionCostsUpToDate = False
        
    """
        ===============================================
                        GETTERS
        ===============================================
    """
    def contains(self,a):
        return a in self.sequence
    
    def getIdentifier(self):
        return self.identifier
    
    def getSatellite(self):
        return self.satellite
    
    def getSequence(self):
        return self.sequence
        
    def getStartDate(self):
        return self.startDate
    
    def getEndDate(self):
        return self.endDate
    
    # return iLeft,iRight
    # planEarly[:iLeft] = activities a gauche
    # si iLeft<len(planEarly) : planEarly[iLeft:iRight] = activities en conflit
    # si iRight<len(getLatestPlan) : getLatestPlan[iRight:] : activities a droite
    def getConflictingActivitiesIndices(self,constellation,activity,transitionModel):
        planEarly,getLatestPlan = self.planEst,self.planLts
        if config.getOptValue("verif"):
            for i,m in enumerate(planEarly):
                if not planEarly[i][1]<=getLatestPlan[i][1]:
                    print(self.identifier,activity,planEarly,getLatestPlan)
                    print(planEarly[i],getLatestPlan[i])
                    assert(False)
        
        # cas vide
        if len(planEarly)==0:
            return 0,0
        
        startDate = activity.getStartDate()
        endDate = activity.getEndDate()
        iLeft = len(getLatestPlan) # position acceptable minimum
        iRight = 0 # position acceptable maximum
        conflitIndice = lambda i : intersect( (planEarly[i][1],getLatestPlan[i][1]),(startDate,endDate))
        isLeft = lambda i : getLatestPlan[i][1] < startDate
        isRight = lambda i : endDate < planEarly[i][1]
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
    def updateCriticalPlans(self,constellation,transitionModel,force=False):
        #assert(not transitionModel.isTimeDependent())
        if not self.arePlansUpToDate() or force:
            self._computeEarliestPlan(constellation,transitionModel)
            self._computeLatestPlan(constellation,transitionModel)
            if config.verifMode():
                for i in range(len(self.planEst)):
                    if not(self.planEst[i][1]<=self.planLts[i][1]):
                        die(self.planEst[i][1],self.planLts[i][1])
            
            
        return self.tardiness  <= 10**(-config.glob.digits)     
            
    def _computeEarliestPlan(self,constellation,transitionModel):
        self.planEst,self.transitionCosts,self.tardiness = self.getEarliestPlan(constellation,self.sequence,transitionModel,addTransitionTardiness=True)
        self.isEarliestPlanUpToDate = True
        self.transitionCostsUpToDate = True
        assert(len(self.sequence)==len(self.planEst))
        
    def _computeLatestPlan(self,constellation,transitionModel):
        #assert(not transitionModel.isTimeDependent())
        self.planLts = self.getLatestPlan(constellation,self.sequence,transitionModel)
        self.isLatestPlanUpToDate = True
  
    def getEarliestTransition(self,constellation,s,a1,a2,start,transitionModel):
        if transitionModel.isTimeDependent():
            transition = constellation.getSatellite(s).getTransitionDurationTimeDependent(a1,a2,start,transitionModel)
        else:
            transition = constellation.getSatellite(s).getTransitionDuration(a1,a2,transitionModel)                    
        return transition             

    def getEarliestPlan(self,constellation,sequence,transitionModel,addTransitionTardiness=False,addTransition=False):
        try:
            totalTransition = 0
            plan = []
            totalTardiness = 0
            s = self.satellite
            for i,activity in enumerate(sequence):
                if(i==0):
                    t = constellation.satellites[s].getActivity(activity).getStartDate()
                    
                else:
                    prec = sequence[i-1]
                    duration = constellation.getSatellite(s).getActivity(prec).getDuration()
                    start = constellation.satellites[s].getActivity(activity).getStartDate()
                    end = constellation.satellites[s].getActivity(activity).getEndDate()
                    transition = self.getEarliestTransition(constellation,s,prec,activity,t+duration,transitionModel)
                    totalTransition += transition
                    totalTardiness +=  round(max(t + duration + transition - end,0),config.glob.digits)
                    t = round(max(t + duration + transition,start),config.glob.digits)
                if addTransition and i>0:
                    plan[-1] = (plan[-1][0],plan[-1][1],transition)
                plan.append((activity,t))
            if addTransitionTardiness:
                return plan,totalTransition,totalTardiness
            else:
                return plan
        except KeyError:
            print(s,sequence)
            raise KeyError


    def getLatestPlan(self,constellation,sequence,transitionModel):
        plan = []
        s = self.satellite
        for i in range(len(sequence)-1,-1,-1):
            a = sequence[i]
            duration = constellation.getSatellite(s).getActivity(a).getDuration()
            if(i==len(sequence)-1):
                t = constellation.getSatellite(s).getActivity(a).getEndDate() - duration
                transition = 0
            else:
                successor = plan[-1][0]
                if transitionModel.isTimeDependent():
                    start,transition = computeTransitionToNextTimeDepPoint(constellation,s,a,successor,plan[-1][1]-10**(-6),transitionModel,ensureFeasible=True)
                    startActivity = start-duration
                    t = min(startActivity,constellation.getSatellite(s).getActivity(a).getEndDate()-duration)
                    if t<self.planEst[i][1]:
                        self.notifyPlansUpToDate(False,False)
                        raise NumericalError()
                else:
                    transition = constellation.getSatellite(s).getTransitionDuration(a,successor,transitionModel)
                    tbefore = t
                    t = min(t - duration - transition,constellation.getSatellite(s).getActivity(a).getEndDate()-duration)
                # correction arrondi
                if transitionModel.isTimeDependent():
                    trans = constellation.getSatellite(s).getTransitionDurationTimeDependent(a,successor,t+duration,transitionModel)
                    assert(t+trans<=start+transition)
            if config.verifMode():
                if len(plan)>0:
                    assert(t+constellation.getSatellite(s).getTransitionDurationTimeDependent(a,successor,t+duration,transitionModel)<=plan[-1][1])
            plan.append((a,t))
  
        plan.reverse()
        if config.verifMode():
            for i,(a,t) in enumerate(plan):
                duration = constellation.getSatellite(s).getActivity(a).getDuration()
                assert(t+duration<=constellation.getSatellite(s).getActivity(a).getEndDate())
                if i<len(plan)-1:
                    successor = plan[i+1][0]
                    assert(t+duration+constellation.getSatellite(s).getTransitionDurationTimeDependent(a,successor,t+duration,transitionModel)<=plan[i+1][1])
        return plan
    

    
    """
        =============================================== 
                        RECHERCHE LOCALE
        ===============================================        
    """
    def localSearch(self,constellation,solver,transitionModel,groups=None,allowNoSolution=False):
        if groups is None:
            assert(solver == 'LKH')
        else:
            assert(solver in ['OPTW','OPTWGroups'])
        if solver == 'LKH' and self.sequence == []:
            return
        #printColor('len ',MPI.COMM_WORLD.Get_rank(),len(self.sequence),c='y')
        start = min([constellation.getSatellite(self.satellite).getActivity(a).getStartDate() for a in self.sequence])
        #temps = time()
        if solver=='LKH':
            printOpen("Recherche locale cca",self.identifier,c='y')
            printOpen("Write LKH cca",self.identifier,c='c')
            self.writeLKH(constellation,self.sequence,start,transitionModel)
            printClose()
            printOpen("Exec LKH cca",self.identifier,c='c')
            self.execLKH()
            printClose()
            printOpen("Read LKH cca",self.identifier,c='c')
            res = self.readLKH()
            printClose()
            printClose()
            if res is not None:
                self.sequence = res
            self.cleanFolderLKH()
        elif solver=='OPTWGroups':
            t0 = time()
            if transitionModel.usesTimeDependentSolver():
                self.writeOPTWGroupsTimeDep(constellation,self.sequence,start,groups,transitionModel)
            else:
                self.writeOPTWGroups(constellation,self.sequence,start,groups,transitionModel)
            t1 = time()
            if transitionModel.usesTimeDependentSolver():
                retval = self.execOPTWGroupsTimeDep(transitionModel)
                if retval!=0:
                    die(self.lastCallOPTW)
            else:
                retval = self.execOPTWGroups(transitionModel)
            t2 = time()
            if transitionModel.usesTimeDependentSolver():
                res,fulfilledRequests,earliest,latest = self.readOPTWGroupsTimeDep(constellation,groups,allowNoSolution=allowNoSolution)
            else:
                res,fulfilledRequests = self.readOPTWGroups(constellation,groups,allowNoSolution=allowNoSolution)
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
                            print("unknwown activity",a)
                            print("sequence",self.sequence)
                            print("groups",groups)
                            print(res)
                            name = self.lastCallOPTW
                            print("fichier",name)
                            assert(find)
                        else:
                            warn("solver error: unknwon activity",a)
            if retval == 0:
                self.cleanFolderOPTWGroups()
            else:
                print(self.lastCallOPTW+" failed")
            copy = deepcopy(self.sequence)
            self.sequence = res
            # patch du bug du solver : backup de la solution
            if not self.isSequenceFeasible(constellation,transitionModel):
                raise InfeasibleSolutionfromSolverException
            return list(fulfilledRequests.keys())
        elif solver=='OPTW':
            t0 = time()
            self.writeOPTW(constellation,self.sequence,start,groups,transitionModel)
            t1 = time()
            self.execOPTW()
            t2 = time()
            res,fulfilledRequests = self.readOPTW(constellation,groups,allowNoSolution=allowNoSolution)
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
                        print("unknown activity",a)
                        print("sequence",self.sequence)
                        print("groups",groups)
                        print(res)
                        name = self.lastCallOPTW
                        print("file",name)
                    assert(find)
                self.sequence = res
            t4 = time()
            return list(fulfilledRequests.keys())

    
    """ OPTW with groups time-dependent"""
    def execOPTWGroupsTimeDep(self,transitionModel):
        solverPath = "./SingleSatSolver"
        ouputPath = "OPTW/TOURS/eosm/"+self.lastCallOPTW+".tour"
        goDir = "cd ../OPTW;"
        filename = "comp_"+str(self.satellite)+"_"+str(self.identifier)+".comp"
        cmd = goDir + solverPath
        cmd += " -multiTimeSteps "
        cmd += " -fRequests OPTW/INSTANCES/eosm/" + self.lastCallOPTW + ".optw"
        cmd += " -fComponent " + transitionModel.getComponentsPath()+"/"+filename
        cmd += " -tmax "+str(config.OPTW.tmax) +" -restart " +str(config.getOptValue("restart"))
        cmd += " -tmaxGlobal " +str(config.OPTW.tmaxGlobal)
        cmd += " -o " +str(ouputPath)
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
        
    # groups : id groupe => liste d'activities
    def writeOPTWGroupsTimeDep(self,constellation,sequence,start,groups,transitionModel):
        for g in groups:
            for a in groups[g]:
                assert(a in sequence)
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        SCALE = config.OPTW.echelle
        #transition = lambda a1,a2 : constellation.getSatellite(self.satellite).getTransitionDuration(a1,a2,transitionModel)
        s = self.satellite
        endDate = int(SCALE*max([constellation.getSatellite(self.satellite).getActivity(a).getEndDate()-constellation.getSatellite(s).getActivity(a).getDuration() for a in self.sequence]))
        #self.mapping_id_optw = {}
        name = config.donnees.algo_name + "_" + str(rank) + "_" + socket.gethostname() + '_sat_' + str(self.satellite) + '_cca_' + str(self.identifier)
        folder = "../OPTW/OPTW/INSTANCES/eosm/"
        tau_max = transitionModel.getTransititionDurationUpperBound()
        
        with open(folder + name + '.optw','w') as file:
            self.writeSeqInitOPTWGroups(groups[0],file,transitionModel)
            self.writeGroupsOPTWGroupsTimeDep(groups,file)
                   
        self.lastCallOPTW = name
        
    def readOPTWGroupsTimeDep(self,constellation,groups,allowNoSolution=False):
        name = self.lastCallOPTW
        file = "../OPTW/OPTW/TOURS/eosm/"+name+".tour"
        with open(file,'r') as f:
            lines = f.read().splitlines()
            if allowNoSolution and len(lines)==0:
                activities = []
            else:
                assert(len(lines)==3)
                line = lines[0]
                activities = [int(a) for a in line.split(" ")]
            fulfilledRequests = {}
            for a in activities:
                r = constellation.getRequestActivity(a)
                if r not in fulfilledRequests:
                    fulfilledRequests[r] = 1
                else:
                    fulfilledRequests[r] += 1
            if config.getOptValue("verif") and 0 in groups:
                for a in groups[0]:
                    assert(constellation.getRequestActivity(a) in fulfilledRequests)
            earliest = [int(a) for a in lines[1].split(" ")]
            latest = [int(a) for a in lines[2].split(" ")]
        return activities,fulfilledRequests,earliest,latest
        
    """ OPTW with GROUPS """
        
    def execOPTWGroups(self,transitionModel):
        solverPath = "./OrienteeringWithGroupsSolver"
        ouputPath = "OPTW/TOURS/eosm/"+self.lastCallOPTW+".tour"
        goDir = "cd ../OPTW;"
        cmd = goDir + solverPath
        cmd += " -format ccati "
        cmd += "-f OPTW/INSTANCES/eosm/" + self.lastCallOPTW + ".optw"
        cmd += " -tmax "+str(config.OPTW.tmax) +" -restart " +str(config.getOptValue("restart"))
        cmd += " -tmaxGlobal " +str(config.OPTW.tmaxGlobal)
        cmd += " -o " +str(ouputPath)
        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        for line in p.stdout.readlines():
            if config.verifMode():
                print(line)
        retval = p.wait()
        return retval

    def cleanFolderOPTWGroups(self):
        try:
            folder_instance = "../OPTW/OPTW/INSTANCES/eosm/"
            os.unlink(folder_instance+self.lastCallOPTW+".optw")
        except FileNotFoundError:
            pass
        try:
            #for tour in [filename for filename in os.listdir("../LKH3/TSPTW/TOURS/eosm/") if filename.startswith(self.last_call_lkh)]:
                folder_tour = "../OPTW/OPTW/TOURS/eosm/"
                os.unlink(folder_tour+self.lastCallOPTW+'.tour')
        except FileNotFoundError:
            pass
        try:
            folder_tour = "../OPTW/OPTW/TSPTW/PAR/"
            os.unlink(folder_tour+self.lastCallOPTW+'.par')
        except FileNotFoundError:
            pass
        try:
            folder_tour = "../OPTW/OPTW/TSPTW/INIT/"
            os.unlink(folder_tour+self.lastCallOPTW+'.init')
        except FileNotFoundError:
            pass
    
    def readOPTWGroups(self,constellation,groups,allowNoSolution=False):
        name = self.lastCallOPTW
        file = "../OPTW/OPTW/TOURS/eosm/"+name+".tour"
        with open(file,'r') as f:
            lines = f.read().splitlines()
            if allowNoSolution and len(lines)==0:
                activities = []
            else:
                assert(len(lines)==1)
                line = lines[0]
                activities = [int(a) for a in line.split(" ")]
            fulfilledRequests = {}
            for a in activities:
                r = constellation.getRequestActivity(a)
                if r not in fulfilledRequests:
                    fulfilledRequests[r] = 1
                else:
                    fulfilledRequests[r] += 1
            if config.getOptValue("verif") and 0 in groups:
                for a in groups[0]:
                    assert(constellation.getRequestActivity(a) in fulfilledRequests)
        return activities,fulfilledRequests
    
    def writeActivityOPTWGroups(self,constellation,file,a,g,transitionModel):
        SCALE = config.OPTW.echelle
        file.write(str(a)+" ")
        #file.write(str(g)+" ")
        coord = constellation.getSatellite(self.satellite).getActivity(a).getCoordinates()
        vlat = transitionModel.getVitesseLatitude()
        vlong = transitionModel.getVitesseLongitude()
        lat,long = coord[0],coord[1]
        lat_effective,long_effective = lat/vlat,long/vlong
        file.write(str(long_effective)+" ")
        file.write(str(lat_effective)+" ")
        duration = constellation.getSatellite(self.satellite).getActivity(a).getDuration()
        file.write(str(int(SCALE*math.ceil(duration)))+" ")
        startDate = constellation.getSatellite(self.satellite).getActivity(a).getStartDate()
        file.write(str(int(SCALE*math.ceil(startDate)))+" ")
        endDate = constellation.getSatellite(self.satellite).getActivity(a).getEndDate()-constellation.getSatellite(self.satellite).getActivity(a).getDuration()
        file.write(str(int(SCALE*math.floor(endDate)))+'\n')
        
    def writeSeqInitOPTWGroups(self,seq_init,file,transitionModel):
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
        
    # groups : id groupe => liste d'activities
    def writeOPTWGroups(self,constellation,sequence,start,groups,transitionModel):
        for g in groups:
            for a in groups[g]:
                assert(a in sequence)
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        SCALE = config.OPTW.echelle
        #transition = lambda a1,a2 : constellation.getSatellite(self.satellite).getTransitionDuration(a1,a2,transitionModel)
        s = self.satellite
        endDate = int(SCALE*max([constellation.getSatellite(self.satellite).getActivity(a).getEndDate()-constellation.getSatellite(s).getActivity(a).getDuration() for a in self.sequence]))
        #self.mapping_id_optw = {}
        name = config.donnees.algo_name + "_" + str(rank) + "_" + socket.gethostname() + '_sat_' + str(self.satellite) + '_cca_' + str(self.identifier)
        folder = "../OPTW/OPTW/INSTANCES/eosm/"
        tau_max = transitionModel.getTransititionDurationUpperBound()

        with open(folder + name + '.optw','w') as file:
            if not transitionModel.usesTimeDependentSolver():
                sequenceLength = len(sequence)
                file.write(str(sequenceLength)+'\n')
                file.write(str(1/SCALE)+" "+str(tau_max*SCALE)+'\n')
                for g in groups:
                    for a in groups[g]:
                        self.writeActivityOPTWGroups(constellation,file,a,g,transitionModel)
            self.writeSeqInitOPTWGroups(groups[0],file,transitionModel)
            self.writeGroupsOPTWGroups(groups,file)
                    
        self.lastCallOPTW = name
       
        
    """       OPTW classique """
        
    def execOPTW(self):
        tmaxGlobal = config.OPTW.tmaxGlobal
        solverPath = "./OrienteeringSolver"
        ouputPath = "OPTW/TOURS/eosm/"+self.lastCallOPTW+".tour"
        goDir = "cd ../OPTW;"
        cmd = goDir + solverPath
        cmd += " -format ccati "
        cmd += "-f OPTW/INSTANCES/eosm/" + self.lastCallOPTW + ".optw"
        cmd += " -tmax "+str(config.OPTW.tmax) +" -restart " +str(config.OPTW.restart)
        cmd += " -tmaxGlobal " +str(tmaxGlobal)
        cmd += " -o " +str(ouputPath)
        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        for line in p.stdout.readlines():
            print(line)
        retval = p.wait()
    
    def writeActiviteOPTW(self,constellation,file,a,g):
        SCALE = config.OPTW.echelle
        file.write(str(a)+" ")
        #file.write(str(g)+" ")
        coord = constellation.getSatellite(self.satellite).getActivity(a).getCoordinates()
        long,lat = coord[0],coord[1]
        file.write(str(long)+" ")
        file.write(str(lat)+" ")
        duration = constellation.getSatellite(self.satellite).getActivity(a).getDuration()
        file.write(str(int(math.ceil(duration*SCALE)))+" ")
        startDate = constellation.getSatellite(self.satellite).getActivity(a).getStartDate()
        file.write(str(int(math.ceil(startDate*SCALE)))+" ")
        endDate = constellation.getSatellite(self.satellite).getActivity(a).getEndDate()-constellation.getSatellite(self.satellite).getActivity(a).getDuration()
        file.write(str(int(math.floor(endDate*SCALE)))+'\n')
    
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
        
    # groups : id groupe => liste d'activities
    def writeOPTW(self,constellation,sequence,start,groups,transitionModel):
        for g in groups:
            for a in groups[g]:
                assert(a in sequence)
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        SCALE = config.OPTW.echelle
        
        transition = lambda a1,a2 : constellation.getSatellite(self.satellite).getTransitionDuration(a1,a2)
        s = self.satellite
        endDate = int(SCALE*max([constellation.getSatellite(self.satellite).getActivity(a).getEndDate()-constellation.getSatellite(s).getActivity(a).getDuration() for a in self.sequence]))
        self.mapping_id_optw = {}
        name = config.donnees.algo_name + "_" + str(rank) + "_" + socket.gethostname() + '_sat_' + str(self.satellite) + '_cca_' + str(self.identifier)
        folder = "../OPTW/OPTW/INSTANCES/eosm/"
        tau_max = transitionModel.getTransititionDurationUpperBound()
        with open(folder + name + '.optw','w') as file:
            sequenceLength = len(sequence)
            file.write(str(sequenceLength)+'\n')
            file.write(str(1/SCALE)+" "+str(tau_max*SCALE)+'\n')

            for g in groups:
                for a in groups[g]:
                    self.writeActiviteOPTW(constellation,file,a,g)
                    self.mapping_id_optw[a] = a
            self.writeSeqInitOPTW(groups[0],file)
            self.writeGroupsOPTW(groups,file)
                    
        self.lastCallOPTW = name    
    
    def cleanFolderOPTW(self):
        try:
            folder_instance = "../OPTW/OPTW/INSTANCES/eosm/"
            os.unlink(folder_instance+self.lastCallOPTW+".optw")
        except FileNotFoundError:
            pass
        try:
            #for tour in [filename for filename in os.listdir("../LKH3/TSPTW/TOURS/eosm/") if filename.startswith(self.last_call_lkh)]:
                folder_tour = "../OPTW/OPTW/TOURS/eosm/"
                os.unlink(folder_tour+self.lastCallOPTW+'.tour')
        except FileNotFoundError:
            pass
        try:
            folder_tour = "../OPTW/OPTW/TSPTW/PAR/"
            os.unlink(folder_tour+self.lastCallOPTW+'.par')
        except FileNotFoundError:
            pass
        try:
            folder_tour = "../OPTW/OPTW/TSPTW/INIT/"
            os.unlink(folder_tour+self.lastCallOPTW+'.init')
        except FileNotFoundError:
            pass
    
    def readOPTW(self,constellation,groups,allowNoSolution=False):
        name = self.lastCallOPTW
        file = "../OPTW/OPTW/TOURS/eosm/"+name+".tour"
        with open(file,'r') as f:
            lines = f.read().splitlines()
            if allowNoSolution and len(lines)==0:
                activities = []
            else:
                assert(len(lines)==1)
                line = lines[0]
                activities = [int(a) for a in line.split(" ")]
            fulfilledRequests = {}
            for a in activities:
                r = constellation.getRequestActivity(a)
                if r not in fulfilledRequests:
                    fulfilledRequests[r] = 1
                else:
                    fulfilledRequests[r] += 1
            if config.getOptValue("verif") and 0 in groups:
                for a in groups[0]:
                    assert(constellation.getRequestActivity(a) in fulfilledRequests)
        return activities,fulfilledRequests

    
   
    """      LKH """
       
    # Ecrire le fichier de parametre pour LKH
    def writeParameterFileLKH(self):
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
    def writeLKH(self,constellation,seq,start,transitionModel):
        SCALE = config.lkh.echelle
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        startScale = SCALE*start
        transition = lambda a1,a2 : constellation.getSatellite(self.satellite).getTransitionDuration(a1,a2,transitionModel)
        s = self.satellite
        endDate = int(SCALE*max([constellation.getSatellite(self.satellite).getActivity(a).getEndDate()-constellation.getSatellite(s).getActivity(a).getDuration() for a in self.sequence]))
        mapping_id = {}
        sequence = [-1] + seq
        name = config.getOptValue("solver") + "_" + str(rank) + "_" + socket.gethostname() + '_sat_' + str(self.satellite) + '_cca_' + str(self.identifier)
        folder = "../LKH3/TSPTW/INSTANCES/eosm/"
        
        with open(folder + name + '.tsptw','w') as file:
            sequenceLength = len(sequence)
            #printColor('ecriture',MPI.COMM_WORLD.Get_rank(),sequenceLength,c='y')
            file.write('NAME : '+name+'.tsptw\n')
            file.write('TYPE : TSPTW\n')
            file.write('DIMENSION : ' + str(sequenceLength) + '\n')
            file.write('EDGE_WEIGHT_TYPE : EXPLICIT\n')
            file.write('EDGE_WEIGHT_FORMAT : FULL_MATRIX\n')
            file.write('EDGE_WEIGHT_SECTION\n')

            # distance entre chaque noeuds : transition + duration
            for a1 in sequence:
                for a2 in sequence:
                    if a1==a2 or a1==-1 or a2==-1:
                        file.write('0 ')
                    else:
                        duration = constellation.getSatellite(self.satellite).getActivity(a1).getDuration()
                        file.write(str(int(SCALE*(transition(a1,a2)+duration)))+' ')
                file.write('\n')
            file.write('TIME_WINDOW_SECTION\n')
            for i,a in enumerate(sequence):
                if a==-1: # fenetre = (startDate,endDate-duration)
                    file.write('1 ' + str(startScale) + ' ' + str(endDate) + '\n')
                else:
                    startDate = constellation.getSatellite(self.satellite).getActivity(a).getStartDate()
                    endDate = constellation.getSatellite(self.satellite).getActivity(a).getEndDate()-constellation.getSatellite(s).getActivity(a).getDuration()
                    file.write(str(i+1) + ' ' + str(int(SCALE*startDate)) + ' ' + str(int(SCALE*endDate)) + '\n')
                    #mapping_id[a] = a
                    mapping_id[i+1] = a
            file.write('DEPOT_SECTION\n')
            file.write('1\n')
            file.write('-1\n')
            file.write('EOF')
        self.last_call_lkh = name
        self.mapping_id_lkh = mapping_id
        self.writeParameterFileLKH()
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
            print("Failure while reading the result",best)
            with open("../LKH3/TSPTW/TOURS/eosm/"+best,'r') as file:
                lines = file.read().splitlines()
                i = 0
                while int(lines[i])>0: #-1 -> fin du fichier
                    print(lines[i])
                    i+=1
            raise exc               
        
    def wasSolverUsed(self):
        return self.usedSolver
    
    def computeNewSequence(self,constellation,transitionModel,solver='LKH',groups=None):
        assert(solver in ['LKH','OPTW','OPTWGroups'])
        shiftRightDisplay(5)
        t1 = time()
        if not self.isSequenceFeasible(constellation,transitionModel):
            sequenceCopy = deepcopy(self.sequence)
            try:
                self.usedSolver = True
                res = self.localSearch(constellation,solver,transitionModel,groups)
                shiftLeftDisplay(5)
                if config.getOptValue("verif") and solver =='OPTW':
                    if not self.isSequenceFeasible(constellation,transitionModel):
                        die("CCA : "+str(self.identifier)+ " " +str(self.sequence)+" unfeasible. Total tardiness: "+ str(self.getComponentTardiness(constellation,transitionModel,printInfo=True)))
                if solver in ['OPTW','OPTWGroups']:
                    return True,res # sequence forcement feasible, renvoyer les requetes satisfaites
                else:
                    return self.isSequenceFeasible(constellation,transitionModel)
            except InfeasibleSolutionfromSolverException:
                if config.verifMode():
                    print("Solver error: infeasible solution.")
                self.sequence = sequenceCopy
                return False,None
        else:
            shiftLeftDisplay(5)
            self.usedSolver = False
            if solver=='LKH':
                return True
            else:
                req = []
                for g in groups:
                    if len(groups[g])>0:
                        req.append(constellation.getRequestActivity(groups[g][0]))
                return True,req
            
    def insertMode(self,constellation,r,m,transitionModel,cheap=False):
        adds = []
        # ajout des activities            
        for s in constellation.getRequest(r).getMode(m).getActivities():
            for o in constellation.getRequest(r).getMode(m).getActivities()[s]:
                if o in self.activities and o not in self.sequence:
                    if cheap:
                        self.insertActivityCheapMethod(constellation,o,transitionModel,feasibilityCheck=False) # ajouter a un endroit intelligent
                    else:
                        self.insertActivityTimeConsumingMethod(constellation,o,transitionModel,feasibilityCheck=False) # ajouter a un endroit intelligent
                adds.append(o)
        assert(len(self.sequence) > 0)
        return adds
    
    def insertListOfModes(self,constellation,sortedModesList,transitionModel):
        adds = []
        #if not config.isTimeDependentModeOn():
        if not self.isEarliestPlanUpToDate:
            assert(not self.isLatestPlanUpToDate)
            feasible = self.updateCriticalPlans(constellation,transitionModel)
        assert(self.isLatestPlanUpToDate and self.isEarliestPlanUpToDate)
        
        for (r,m) in sortedModesList:
            for a in sortedModesList[(r,m)]:
                if config.isTimeDependentModeOn():
                    self.insertActivities(constellation,a,transitionModel,feasibilityCheck=False) # ajouter a un endroit intelligent
                else:
                    if feasible:
                        self.insertActivityTCriticalPlansMethod(constellation,a,transitionModel)
                        feasible = self.updateCriticalPlans(constellation,transitionModel)
                    else:
                        self.insertActivityCheapMethod(constellation,a,transitionModel,feasibilityCheck=False) # ajouter ou la transition est minimale
                adds.append(a)
        assert(len(self.sequence) > 0)
        return adds        
        
    def insertActivities(self,constellation,activities,transitionModel,cheap=False):
        adds = []
        # ajout des activities            
        for o in sorted(activities):
            if o not in self.sequence:
                if cheap:
                    self.insertActivityCheapMethod(constellation,o,transitionModel,feasibilityCheck=False) # ajouter a un endroit intelligent
                else:
                    self.insertActivityTimeConsumingMethod(constellation,o,transitionModel,feasibilityCheck=False) # ajouter a un endroit intelligent
                adds.append(o)
        assert(len(self.sequence) >0)
        return adds
    
    def getSequenceTardiness(self,constellation,sequence,transitionModel,start=None,printInfo=False):       
        s = self.satellite
        tardiness = 0
        if start is None and len(sequence)>0:
            start = constellation.getSatellite(s).getActivity(sequence[0]).getStartDate()
            t = start
            for i,a in enumerate(sequence):
                t = round(max(t,constellation.getSatellite(s).getActivity(a).getStartDate()),config.glob.digits)
                duration = constellation.getSatellite(s).getActivity(a).getDuration()
                tardiness += max(0,-constellation.getSatellite(s).getActivity(a).getEndDate()+t+duration)# tardiness tronqué
                if printInfo:
                    print(t,constellation.getSatellite(s).getActivity(a),tardiness)                
                if(i<len(sequence)-1):
                    if transitionModel.isTimeDependent():
                        t += round(duration + transitionModel.getTransitionDuration(a,sequence[i+1], t+duration),config.glob.digits)
                    else:
                        t += round(duration + constellation.satellites[s].getTransitionDuration(a,sequence[i+1],transitionModel),config.glob.digits)
        return tardiness
 
    def getSignedTardinessOfSequence(self,constellation,sequence,transitionModel,start):
        s = self.satellite
        tardiness = 0
        t = start
        for i,a in enumerate(sequence):
            t += constellation.getSatellite(s).getActivity(a).getDuration()
            tardiness += t-constellation.satellites[s].getActivity(a).getEndDate()
            if(i<len(sequence)-1):
                if transitionModel.isTimeDependent():
                    t += transitionModel.getTransitionDuration(a,sequence[i+1], t)
                else:
                    t += constellation.satellites[s].getTransitionDuration(a,sequence[i+1],transitionModel)
                t = round(max(t,constellation.satellites[s].getActivity(a).getStartDate()),config.glob.digits)
        return tardiness

    # ( unsigned tardiness , signed tardiness )
    def getTardinessOfSequence(self,constellation,sequence,transitionModel,start=None):
        s = self.satellite
        tardiness = 0
        signedTardiness = 0
        if start is None and len(sequence)>0:
            start = constellation.getSatellite(s).getActivity(sequence[0]).getStartDate()
        t = start
        for i,a in enumerate(sequence):
            t += constellation.getSatellite(s).getActivity(a).getDuration()
            tardiness += max(0,t-constellation.satellites[s].getActivity(a).getEndDate())# tardiness tronqué
            signedTardiness += max(0,t-constellation.satellites[s].getActivity(a).getEndDate())# tardiness tronqué
            if(i<len(sequence)-1):
                if transitionModel.isTimeDependent():
                    t += transitionModel.getTransitionDuration(a,sequence[i+1], t)
                else:
                    t += constellation.satellites[s].getTransitionDuration(a,sequence[i+1],transitionModel)
                t = max(t,constellation.satellites[s].getActivity(a).getStartDate())
        return tardiness,signedTardiness
    
    """
        ===============================================
                    INSERTION DES ACTIVITES
        ===============================================
    """
    # solution = proposition d'ordre pour la solution, pas la solution courante
    # info supplementaire pour otpimiser le calcul du tardiness : tardinesss a chaque position + dates early
    def partitionConflictingNewActivityMinimizingTardiness(self,constellation,p,transitionModel):
        assert(not transitionModel.isTimeDependent())
        s = self.satellite
        a,b = constellation.getSatellite(s).getActivity(p).getStartDate(),constellation.getSatellite(s).getActivity(p).getEndDate()
        L,C,R = [],[],[]
        Tearly = []
        tardinessList = []
        for i,pp in enumerate(self.sequence):
            # calcul de la duration
            duration = constellation.getSatellite(s).getActivity(p).getDuration()
            # determination des 3 listes
            if i==0:
                Tearly.append(constellation.getSatellite(s).getActivity(pp).getStartDate())
                tardinessList.append((0,0))
            else:
                prev = self.sequence[i-1]
                duration = constellation.getSatellite(s).getActivity(pp).getDuration()
                tearly = Tearly[-1]+constellation.getSatellite(s).getTransitionDuration(prev,pp,transitionModel)
                tearly += duration
                tearly = max(tearly,constellation.getSatellite(s).getActivity(pp).getStartDate())
                tardiness_alg = constellation.getSatellite(s).getActivity(pp).getEndDate()-duration-tearly
                tardiness = max(tardiness_alg,0)
                Tearly.append(tearly)
                tardinessList.append((tardiness+tardinessList[-1][0],tardiness_alg+tardinessList[-1][1]))
            if(constellation.getSatellite(s).getActivity(pp).getEndDate()+constellation.getSatellite(s).getTransitionDuration(pp,p,transitionModel)<a):
                L.append(pp)
            elif b + constellation.getSatellite(s).getTransitionDuration(p,pp,transitionModel) < constellation.getSatellite(s).getActivity(pp).getStartDate():
                R.append(pp)
            else:
                C.append(pp)
        return L,C,R,Tearly,tardinessList

    def areActivitiesConflicting(self,a1,a2,constellation,transitionModel):
        s = self.satellite
        startDate1,endDate1 = constellation.getSatellite(s).getActivity(a1).getStartDate(),constellation.getSatellite(s).getActivity(a2).getEndDate()
        startDate2,endDate2 = constellation.getSatellite(s).getActivity(a2).getStartDate(),constellation.getSatellite(s).getActivity(a2).getEndDate()
        test1 = endDate1+transitionModel.getMajorantCoupleActivite(a1,a2)<startDate2
        test2 = endDate2+transitionModel.getMajorantCoupleActivite(a2,a1)<startDate1
        return not test1 and not test2
    
    # solution = proposition d'ordre pour la solution, pas la solution courante
    def partitionConflictWithActivity(self,constellation,p,transitionModel):
        s = self.satellite
        a,b = constellation.getSatellite(s).getActivity(p).getStartDate(),constellation.getSatellite(s).getActivity(p).getEndDate()
        L,C,R = [],[],[]
        for pp in self.sequence:
            # calcul de la duration
            duration = constellation.getSatellite(s).getActivity(p).getDuration()
            if transitionModel.isTimeDependent():
                leftTransition = transitionModel.getMajorantCoupleActivite(pp,p)
                rightTransition = transitionModel.getMajorantCoupleActivite(p,pp)
            else:
                rightTransition = constellation.getSatellite(s).getTransitionDuration(pp,p,transitionModel)
                leftTransition = rightTransition
            # determination des 3 listes
            if(constellation.getSatellite(s).getActivity(pp).getEndDate()+leftTransition<a):
                L.append(pp)
            elif b + rightTransition < constellation.getSatellite(s).getActivity(pp).getStartDate():
                R.append(pp)
            else:
                C.append(pp)
        return L,C,R
    
    # not self.transitionModel.isTimeDependent()
    
    def getComponentTardiness(self,constellation,transitionModel,force=False,printInfo=False):
        if not self.isEarliestPlanUpToDate or force:
            if printInfo:
                printColor("Updating critical plans.")
            self.updateCriticalPlans(constellation,transitionModel,force=force)
        else:
            if printInfo:
                printColor("Critical plans up to date")
        return self.tardiness

    def getTotalTransitionCostOfComponent(self,constellation,transitionModel):
        if not self.transitionCostsUpToDate:
            self.updateCriticalPlans(constellation,transitionModel)
        return self.transitionCosts
    
    def isSequenceFeasible(self,constellation,transitionModel,printTardiness=False):
        if self.sequence is None:
            return False
        if(len(self.sequence))==0:
            return True
        s = self.satellite
        t = constellation.getSatellite(s).getActivity(self.sequence[0]).getStartDate()
        for i,p in enumerate(self.sequence):
            if printTardiness:
                print(t,constellation.getSatellite(s).getActivity(p))
            duration = constellation.getSatellite(s).getActivity(p).getDuration()
            if t+duration-constellation.getSatellite(s).getActivity(p).getEndDate()>10**(-config.glob.digits):
                if printTardiness:
                    printColor(t,constellation.getSatellite(s).getActivity(p),c='r')
                return False
            if i<len(self.sequence)-1:
                if transitionModel.isTimeDependent():
                    transition = constellation.getSatellite(s).getTransitionDurationTimeDependent(p,self.sequence[i+1],t+duration,transitionModel)
                else:
                    transition = constellation.getSatellite(s).getTransitionDuration(p,self.sequence[i+1],transitionModel)
                if printTardiness:
                    print(p,self.sequence[i+1],t+duration,transition,t+duration>constellation.getSatellite(s).getActivity(p).getEndDate()>10**(-config.glob.digits))
                t = round(max(t + transition + duration,constellation.getSatellite(s).getActivity(self.sequence[i+1]).getStartDate()),config.glob.digits)
        return True
    
    def extensionConflict(self,L,C,R,constellation,transitionModel):
        if len(L)>0:
            start_L = self.getEarliestPlan(constellation,L,transitionModel)[-1][1]
            start = lambda c : start_L
        else: # c_extend toujours != []
            start = lambda c : constellation.getSatellite(self.satellite).getActivity(c[0]).getStartDate()
        if L==[]:
            C_left = lambda c : c
        else:
            C_left = lambda c : [L[-1]] + c
        if R==[]:
            C_right = lambda c : c
        else:
            C_right = lambda c : c + [R[0]]
        C_extend = lambda c : C_left(C_right(c))
        return start,C_extend
   
    # ( tardiness , tardiness algébrique )
    def getSignedUnsignedTardiness(self,constellation,sequence,transitionModel,start=None):
        s = self.satellite
        tardiness = 0
        signedTardiness = 0
        if start is None and len(sequence)>0:
            start = constellation.getSatellite(s).getActivity(sequence[0]).getStartDate()
        t = start
        for i,a in enumerate(sequence):
            t += constellation.getSatellite(s).getActivity(a).getDuration()
            tardiness += round(max(0,t-constellation.satellites[s].getActivity(a).getEndDate()),config.glob.digits)# tardiness tronqué
            signedTardiness +=round( max(0,t-constellation.satellites[s].getActivity(a).getEndDate()),config.glob.digits)# tardiness tronqué
            if(i<len(sequence)-1):
                t += constellation.satellites[s].getTransitionDuration(a,sequence[i+1],transitionModel)
                t = round(max(t,constellation.satellites[s].getActivity(a).getStartDate()),config.glob.digits)
        return tardiness,signedTardiness
    
    def getOptimizedTardinessCriteria(self,constellation,p,i,L,C,R,Tearly,tardinessList):
        delta_tardiness = (0,0)
        sequenceRank = len(L)+i
        duration = constellation.getSatellite(self.satellite).getActivity(p).getDuration()
        startDate = constellation.getSatellite(self.satellite).getActivity(p).getStartDate()
        endDate = constellation.getSatellite(self.satellite).getActivity(p).getEndDate()
        # date early de p
        if sequenceRank==0:
            earlyDate = constellation.getSatellite(self.satellite).getActivity(p).getStartDate()
        else:
            prev = self.sequence[sequenceRank-1]
            earlyDate = Tearly[sequenceRank-1]
            earlyDate += constellation.getSatellite(self.satellite).getActivity(prev).getDuration()
            earlyDate = max(earlyDate,startDate) 
        # tardiness de p
        tardiness_alg = endDate - duration - earlyDate
        tardiness = (round(max(0,tardiness_alg),config.glob.digits),round(tardiness_alg,config.glob.digits))
        # date early du successeur
        if sequenceRank==len(self.sequence)-1:
            return tardiness
        else:
            delta = tardiness
            t_k = earlyDate
            for k in range(sequenceRank,len(self.sequence)-1):
                ak = constellation.getSatellite(self.satellite).getActivity(self.sequence[k])
                if k==0:
                    t_k = round(max(earlyDate + duration,ak.getStartDate(),0),config.glob.digits)
                else:
                    prev = constellation.getSatellite(self.satellite).getActivity(self.sequence[k-1])
                    t_k = round(max(t_k + prev.getDuration(),ak.getStartDate()),config.glob.digits)
                tardiness_k_alg = ak.getEndDate()-t_k-ak.getDuration()
                tardiness_k = (config.glob.digits(max(tardiness_alg,0),config.glob.digits),round(tardiness_alg,config.glob.digits))
                delta = (delta[0]+tardiness[0],delta[1]+tardiness[1])
                if delta==tardinessList[k]:
                    return delta
            return delta
    
    def insertActivityTCriticalPlansMethod(self,constellation,p,transitionModel):
        #assert(not transitionModel.isTimeDependent())
        def computesTearly(constellation,p,i,transitionModel):
            s = self.satellite
            if i>0:
                prev = self.planEst[i-1]
                Tearly = prev[1]
                duration = constellation.getSatellite(s).getActivity(prev[0]).getDuration()
                Tearly += duration
                startPreviousManeuver = Tearly
                if not transitionModel.isTimeDependent():
                    Tearly += constellation.getSatellite(s).getTransitionDuration(prev[0],p,transitionModel)
                else:
                    Tearly += constellation.getSatellite(s).getTransitionDurationTimeDependent(prev[0],p,prev[1]+duration,transitionModel)
                Tearly = max(Tearly,constellation.getSatellite(s).getActivity(p).getStartDate())
            else:
                activity = constellation.getSatellite(s).getActivity(p)
                Tearly = activity.getStartDate()
                startPreviousManeuver = None
            return Tearly,startPreviousManeuver
        
        def computeLatestDate(p,i,transitionModel):
            if i<len(self.planLts):
                next_a = self.planLts[i]
                latestStartNextActivity = next_a[1]
                if transitionModel.isTimeDependent():
                    latestStartManeuver,transToNextDuration = computeTransitionToNextTimeDepPoint(constellation,s,p,next_a[0],next_a[1]-10e-4,transitionModel)
                    self.verifyTransitions(transitionModel,p,next_a[0],latestStartManeuver,transToNextDuration,next_a[1])
                
                else:
                    transToNextDuration =  float(constellation.getSatellite(s).getTransitionDuration(p,next_a[0],transitionModel))
                    latestStartManeuver = latestStartNextActivity - transToNextDuration
            else:
                latestStartManeuver = np.Inf
                transToNextDuration = 0
                latestStartNextActivity = np.Inf
            return latestStartManeuver,latestStartNextActivity,transToNextDuration
        # a appeler avant l'insertion
                
        def computeDeltaTransitionForNewInsertion(sequence,constellation,p,insertionPosition,transitionModel,startPreviousManeuver):
            if not transitionModel.isTimeDependent():
                # ajout de previous(p) -> p
                addTransitionLeft = 0
                if insertionPosition>0:
                    addTransitionLeft += constellation.getSatellite(s).getTransitionDuration(sequence[insertionPosition-1],p,transitionModel)
                # ajout de p -> next(p)
                addTransitionRight = 0
                if insertionPosition<len(sequence):
                    addTransitionRight += constellation.getSatellite(s).getTransitionDuration(p,sequence[insertionPosition],transitionModel)
                # perte de la transition previous(p) -> next(p) 
                transitionLoss = 0
                if insertionPosition > 0 and insertionPosition<len(sequence):
                    transitionLoss += constellation.getSatellite(s).getTransitionDuration(sequence[insertionPosition-1],sequence[insertionPosition],transitionModel)
                return addTransitionLeft + addTransitionRight - transitionLoss
            else:
                # ajout de previous(p) -> p
                addTransitionLeft = 0
                if insertionPosition>0:
                    addTransitionLeft += constellation.getSatellite(s).getTransitionDurationTimeDependent(sequence[insertionPosition-1],p,startPreviousManeuver,transitionModel)
                # ajout de p -> next(p)
                addTransitionRight = 0
                if insertionPosition<len(sequence):
                    if startPreviousManeuver is not None:
                        start = max(startPreviousManeuver+addTransitionLeft,constellation.getSatellite(s).getActivity(p).getStartDate())
                        startManeuverToNext = start + constellation.getSatellite(s).getActivity(p).getDuration()
                        addTransitionRight += constellation.getSatellite(s).getTransitionDurationTimeDependent(p,sequence[insertionPosition],startManeuverToNext,transitionModel)
                # perte de la transition previous(p) -> next(p) 
                transitionLoss = 0
                if insertionPosition > 0 and insertionPosition<len(sequence):
                    transitionLoss += constellation.getSatellite(s).getTransitionDurationTimeDependent(sequence[insertionPosition-1],sequence[insertionPosition],startPreviousManeuver,transitionModel)
                return addTransitionLeft + addTransitionRight - transitionLoss
                
        def computeOptimalInsertionPosition(constellation,p,iLeft,iRight,transitionModel):
            assert(self.arePlansUpToDate())
            activity = constellation.getActivity(p)
            optimalInsertionPosition,optiCriteria,timeOpti = -1,(0,-np.Inf),np.Inf  # feasible,flexibilite
            for i in range(iLeft,iRight+1):
                try:
                    t12 = time()
                    assert(self.arePlansUpToDate())
                    assert(len(self.sequence)==len(self.planLts))
                    Tearly,startPreviousManeuver = computesTearly(constellation,p,i,transitionModel)
                    t13 = time()
                    latestStartManeuver,latestStartNextActivity,transToNextDuration = computeLatestDate(p,i,transitionModel)
                    t14 = time()
                    feasible = Tearly + activity.getDuration() <= activity.getEndDate() and Tearly + activity.getDuration() <= latestStartManeuver

                    if config.verifMode():
                        if i<len(self.sequence)-1:
                            assert(constellation.getSatellite(s).getTransitionDurationTimeDependent(p,self.sequence[i],latestStartManeuver,transitionModel) <= self.planLts[i][1])
                        
                    if config.verifMode() and feasible:
                        if optimalInsertionPosition<len(self.sequence):
                            self.sequence = self.sequence[:i]+[p]+self.sequence[i:]
                        else:
                            self.sequence.append(p)
                        if not(self.isSequenceFeasible(constellation,transitionModel)):
                            for a in self.sequence:
                                print(constellation.getSatellite(s).getActivity(a))
                            print("----------------------------------")
                            
                            print("insertion of ",p," at position ",i)
                            print("start_tau_latest",latestStartManeuver)
                            print("latest_start_next",latestStartNextActivity)
                            print("duration_time",transToNextDuration)
                            print("early",self.planEst)
                            print("late",self.planLts)
                            print("before",self.getEarliestPlan(constellation,self.sequence,transitionModel,addTransition=True))
                            self.sequence.remove(p)
                            print("after",self.getEarliestPlan(constellation,self.sequence,transitionModel,addTransition=True))
                            die(self.sequence,self.getSequenceTardiness(constellation,self.sequence,transitionModel),self.isSequenceFeasible(constellation,transitionModel))
                        self.sequence.remove(p)
                    
                    #Tearly+activity.getDuration()<=activity.getEndDate() and Tearly+activity.getDuration()+trans_to_next<=latest_next
                    flexibilite = latestStartManeuver-(Tearly+activity.getDuration())
                    critere = (feasible,flexibilite)
                    if critere>optiCriteria:
                        optimalInsertionPosition,optiCriteria,timeOpti = i,critere,startPreviousManeuver
                    t15 = time()
                    printColor("Tearly computation",round(t13-t12,7),"latest computation",round(t14-t13,7),"flexibility computation",round(t15-t14,7),c='m')
                except (InfeasibleSolutionException,NumericalError) as error:
                    pass
            return optimalInsertionPosition,optiCriteria,timeOpti
      
        if config.verifMode() :
            if not(self.isSequenceFeasible(constellation,transitionModel)):
                print(self.sequence,optimalInsertionPosition,optiCriteria,startPreviousManeuver)
                assert(False)        

        printOpen("insert activity",p,c='m')
        s = self.satellite
        t0 = time() 
        activity = constellation.getSatellite(s).getActivity(p)
        if self.sequence == []:
            self.sequence = [p]
            feasible = True
            self.tardiness = 0
            self.transitionCosts = 0
            self.isEarliestPlanUpToDate = True
            self.isLatestPlanUpToDate = True
            self.positionLastInserted = 0
        else:
            self.updateCriticalPlans(constellation,transitionModel)
            if not self.arePlansUpToDate():
                raise OutdatedCriticalPlans                   
            
            iLeft,iRight = self.getConflictingActivitiesIndices(constellation,activity,transitionModel)
            if config.getOptValue("verif"):
                self.checkTransitionCost(constellation,transitionModel)
                       
            t1 = time()
            optimalInsertionPosition,optiCriteria,startPreviousManeuver = computeOptimalInsertionPosition(constellation,p,iLeft,iRight,transitionModel)
            t2 = time()
            self.transitionCosts += round(computeDeltaTransitionForNewInsertion(self.sequence,constellation,p,optimalInsertionPosition,transitionModel,startPreviousManeuver),config.glob.digits)
            if optimalInsertionPosition<len(self.sequence):
                self.sequence = self.sequence[:optimalInsertionPosition]+[p]+self.sequence[optimalInsertionPosition:]
            else:
                self.sequence.append(p)    
            feasible = optiCriteria[0]           
            t3 = time()
            printColor("Insertion bounds:"+str(round(t1-t0,7))+"/pos. optimal:"+str(round(t2-t1,7))+"/update  sequence:"+str(round(t3-t2,7)),c='m') 
            
            if config.getOptValue("verif"):
                self.checkTransitionCost(constellation,transitionModel)
            self.positionLastInserted = optimalInsertionPosition
        self.notifyPlansUpToDate(False,False)
        self.transitionCostsUpToDate = True
        self.startDate = min(self.startDate,activity.getStartDate())
        self.endDate = max(self.endDate,activity.getEndDate()+transitionModel.getTransititionDurationUpperBound())
        printClose()
        
        if config.verifMode() and feasible:
            if not(self.isSequenceFeasible(constellation,transitionModel)):
                self.isSequenceFeasible(constellation,transitionModel)
                print(constellation.getSatellite(self.satellite).getActivity(p))
                print(self.planEst)
                print(self.planLts)
                for i,a in enumerate(self.sequence):
                    print(constellation.getSatellite(s).getActivity(a))
                die(p,self.sequence,optimalInsertionPosition,optiCriteria,startPreviousManeuver)
                        
        if feasible:
            try:
                self.updateCriticalPlans(constellation,transitionModel,force=True)
            except NumericalError:
                assert(transitionModel.isTimeDependent())
        
        return feasible
    
    def verifyTransitions(self,transitionModel,a1,a2,startDate,assumedDuration,latest_a2):
        if config.verifMode():
            computerManeuver = transitionModel.getTransitionDuration(a1,a2,startDate)
            if not computerManeuver==assumedDuration:
                die(computerManeuver,assumedDuration)
            if startDate+computerManeuver>latest_a2:
                die(startDate+computerManeuver,latest_a2)
                    
    def checkTransitionCost(self,constellation,transitionModel):
        transitionAfterInsertion = self.transitionCosts
        self.updateCriticalPlans(constellation,transitionModel)
        if not abs(self.transitionCosts-transitionAfterInsertion)<config.glob.getScale():
            print(self.transitionCosts-transitionAfterInsertion)
            print(self.sequence,p,transitionAfterInsertion,self.transitionCosts)
        assert(abs(self.transitionCosts-transitionAfterInsertion)<config.glob.getScale())        
    
    def insertActivityTimeConsumingMethod(self,constellation,p,transitionModel,feasibilityCheck=True):
        lengthBefore = len(self.sequence)
        s = self.satellite
        L,C,R = self.partitionConflictWithActivity(constellation,p,transitionModel)
        # c_etendue : c + extremites = sous-sequence a optimiser
        if len(L)>0:
            start = lambda c : self.getEarliestPlan(constellation,L,transitionModel)[-1][1]
        else: # c_extend toujours != []
            start = lambda c : constellation.getSatellite(s).getActivity(c[0]).getStartDate()
        if L==[]:
            C_left = lambda c : c
        else:
            C_left = lambda c : [L[-1]] + c
        if R==[]:
            C_right = lambda c : c
        else:
            C_right = lambda c : c + [R[0]]
        C_extend = lambda c : C_left(C_right(c))
        # critere 1 = tardiness. Minimiser = mesurer la faisabilite (feasible si critere = 0)
        # critere 2 : tardiness algebrique. Minimiser = augmenter la flexibilite
        optiCriteria,optimizedSequence = None,None
        criteriaFunction = lambda sequence : (self.getSequenceTardiness(constellation,sequence,transitionModel,start=start(sequence)),self.getSignedTardinessOfSequence(constellation,sequence,transitionModel,start(sequence)))
        for i in range(len(C)+1):
            ordre = C[0:i]+[p]+C[i:]
            if(optiCriteria is None):
                optiCriteria,optimizedSequence = criteriaFunction(C_extend(ordre)),ordre
            else:
                critere = criteriaFunction(C_extend(ordre))
                if(critere<optiCriteria):
                    optiCriteria,optimizedSequence = critere,ordre
        self.sequence = L+optimizedSequence+R
        self.notifyPlansUpToDate(False,False)
        assert(len(self.sequence)==lengthBefore+1) 
        if feasibilityCheck:
            return self.isSequenceFeasible(constellation,transitionModel)
        
    def insertActivityCheapMethod(self,constellation,p,transitionModel,feasibilityCheck=True):
        assert(not transitionModel.isTimeDependent())
        def evalPos(constellation,i,a):
            s = self.satellite
            if i == 0:
                transitionLeft = 0
            else:
                transitionLeft = constellation.getSatellite(s).getTransitionDuration(self.sequence[i-1],p,transitionModel)
            if i == len(self.sequence):
                transitionRight = 0
            else:
                transitionRight = constellation.getSatellite(s).getTransitionDuration(p,self.sequence[i],transitionModel)
            return transitionLeft + transitionRight
            
        lengthBefore = len(self.sequence)
        if lengthBefore == 0:
            self.sequence = [p]
            if feasibilityCheck:
                return True
        else:
            bestPosition = None
            bestCriteria = None # transition
            for i in range(len(self.sequence)+1):
                critere = evalPos(constellation,i,p)
                if bestCriteria is None or critere > bestCriteria:
                    bestCriteria = critere
                    bestPosition = i
                    
            if bestPosition ==0:
                self.sequence = [p] + self.sequence[bestPosition:]
            elif bestPosition != len(self.sequence):
                self.sequence = self.sequence[0:bestPosition] + [p] + self.sequence[bestPosition:]
            else:
                self.sequence = self.sequence + [p]        
            self.notifyPlansUpToDate(False,False)
            assert(len(self.sequence)==lengthBefore+1) 
            if feasibilityCheck:
                return self.isSequenceFeasible(constellation,transitionModel)
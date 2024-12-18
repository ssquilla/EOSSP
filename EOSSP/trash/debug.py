from solution import *
from constellation import *
from composantes import *
from solution_composantes import *
from mpi4py import MPI
from Utils import *
from config import *
from time import time
from time import sleep




global config
config = Config()
instance = config.instance



rd.seed(MPI.COMM_WORLD.Get_rank())

path = '../data'
file = choseAndBroadcastFile(path,instance)[0]
constellation = Constellation(file,True)

activites = constellation.extraireActivitesRequetes()
grapheDependances = GroupeComposantesActivitesStatiques(constellation,activites)
print(grapheDependances)
print(str(grapheDependances).replace('\n','\n\t\t\t\t'),type(str(grapheDependances))==str)
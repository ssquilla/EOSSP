"""
from config import *
global config
config = Config()
instance = config.instance
from mpi4py import MPI
if config.getOptValue("help"):
    if MPI.COMM_WORLD.Get_rank()==0:
        config.afficherAide()
else:    
    from solution import *
    from constellation import *
    from composantes import *
    from solution_composantes import *
    from Utils import *
    from time import time
    from time import sleep
    import random as rd       
        
    path = '../data'
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    
    output = path+'/recompose_periodiques'
    
    seed_option = config.getOptValue("seed")
    rd.seed(seed_option)
    
    
    
    file = choseAndBroadcastFile(path,instance)[0]
    constellation = Constellation(file,True)
""" 
path = '../data/recompose_periodic'
id_req = 0
id_obs = 0
with open(path,'r') as file:
    lines = file.readlines()
    with open('../data/output_recompose','w') as output_file:
        for line in lines:
            if 'PERIODIC' in line:
                elmt = line.split('\n')[0].split(',')
                wrong_id_req,a,b,c = elmt[0],elmt[1],elmt[2],elmt[3]
                output_file.write(str(id_req)+','+a+','+b+','+c+'\n')
                id_req += 1
            else: 
                line_obs = line.split('\n')[0].split(",")
                if len(line_obs)==11:
                    line_obs[1] = str(id_obs)
                    id_obs += 1
                    for i,x in enumerate(line_obs):
                        output_file.write(x)
                        if i<len(line_obs)-1:
                            print(',')
                        else:
                            print('\n')
                else:
                    output_file.write(line)
    

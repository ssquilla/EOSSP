import sys, importlib
from pathlib import Path
import random as rd

def import_parents(level=1):
    global __package__
    file = Path(__file__).resolve()
    parent, top = file.parent, file.parents[level]
        
    sys.path.append(str(top))
    try:
        sys.path.remove(str(parent))
    except ValueError: # already removed
        pass
    
    __package__ = '.'.join(parent.parts[len(top.parts):])
    importlib.import_module(__package__) # won't be needed after that
        


if __name__ == '__main__' and __package__ is None:
    import_parents(level=1)

    #!/usr/bin/env python3
    # -*- coding: utf-8 -*-
"""
    Created on Thu Nov  3 09:55:49 2022
    
    @author: ssquilla
"""
from EOSSP.Utils.config import *
global config
config = Config()
instance = config.instance
from EOSSP.model.modelTransition import modelTransitionTimeDependent
from EOSSP.model.composantes import ComposantesStatiquesPrecalculees
from mpi4py import MPI
if config.getOptValue("help"):
    if MPI.COMM_WORLD.Get_rank()==0:
        config.afficherAide()
else: 
    from EOSSP.Utils.Utils import *
    from EOSSP.model.constellation import Constellation
    from EOSSP.model.composantes import *
    from EOSSP.model.solution_composantes import *
    import pandas as pd
    import geopandas as gpd
    import matplotlib.pyplot as plt
    from matplotlib.widgets import Slider
    from math import ceil,floor
    import fcntl
    from shapely.geometry import Polygon,Point,LineString
    
    folder = config.getOptValue("data")
    path = '../data/'+folder
    
    class Locker:
        def __enter__ (self):
            self.fp = open("./lockfile.lck",'w')
            fcntl.flock(self.fp.fileno(), fcntl.LOCK_EX)
    
        def __exit__ (self, _type, value, tb):
            fcntl.flock(self.fp.fileno(), fcntl.LOCK_UN)
            self.fp.close()
    
    def drawScoreMap(file,constellation):
        def getDates(constellation):
            times = []
            for r in constellation.getRequetes():
                for (s,o) in constellation.getRequete(r).getCouples():
                    debut = floor(round(constellation.getSatellite(s).getActivite(o).getDebut(),2))
                    fin = ceil(round(constellation.getSatellite(s).getActivite(o).getFin(),2))
                    if debut not in times:
                        bisect.insort_left(times, debut)
                    if fin not in times:
                        bisect.insort_left(times, fin)
            return times            
            
        def getObsByDates(constellation,times):
            obs_by_dates = {}
            for i,t in enumerate(times):
                if i<len(times)-1:
                    time_interval = (t,times[i+1])
                obs_liste = []
                for r in constellation.getRequetes():
                    for (s,a) in constellation.getRequete(r).getCouples():
                        interval_act = (constellation.getSatellite(s).getActivite(a).getDebut(),constellation.getSatellite(s).getActivite(a).getFin())
                        if intersect(interval_act, time_interval):
                            coord = constellation.getSatellite(s).getActivite(a).getCoordonnees()
                            score = constellation.getSatellite(s).getActivite(a).getScore()
                            obs_liste.append((coord,score))
                obs_by_dates[t] = obs_liste
            return obs_by_dates
        
        def couleur(score):
            left = np.array([0,0,1])
            right = np.array([1,0,0])
            assert(score<=1)
            assert(score>=0)
            return score*right+(1-score)*left
        
        def drawMapSubsetObs(constellation,obs_at_date):
            df = pd.DataFrame(columns=["Couleur","Latitude","Longitude"])
            for (coord,score) in obs_at_date:
                df.append({"Couleur":couleur(score),"Latitude":coord[0],"Longitude":coord[1]},ignore_index=True)
            gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.Longitude, df.Latitude))
            print(gdf.head())
            world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
            ax = world[world.continent == 'Europe'].plot(color='green', edgecolor='black',alpha=0.5)
            ax.set_facecolor('#A8C5DD')
            gdf.plot(ax=ax,color='red')
            plt.show()
        if MPI.COMM_WORLD.Get_rank()==0:
            print("Démarrage ...")
            times = getDates(constellation)
            print("Points temporels crées")
            obs_by_date = getObsByDates(constellation,times)
            print("Observations réparties par dates")
            Tstart = times[0]
            print([len(obs_by_date[t]) for t in obs_by_date])
            drawMapSubsetObs(constellation,obs_by_date[Tstart])
            print("Dessin de la carte des scores")
    
    def drawTransition(file,constellation):
        print("Construction du modèle de transition ...")
        splitted = file.split("/")
        folder = ""
        for i,x in enumerate(splitted):
            if i<len(splitted)-1:
                folder += x + "/"
        folder +="components"
        modele = modelTransitionTimeDependent(folder)
        print("Lecture des composantes ...")
        CCAs = ComposantesStatiquesPrecalculees(constellation,folder)
        print("OK.")
        nCouples = 10
        ccas_candidates = [c for c in CCAs.composantes.keys() if len(CCAs.getActivitesComposante(c))>nCouples]
        cca = rd.choice(list(ccas_candidates))
        f,axs = plt.subplots(nCouples,nCouples)
        activites = rd.choices(CCAs.getActivitesComposante(cca),k=nCouples)
        for i,a1 in enumerate(activites):
            for j,a2 in enumerate(activites):
                points = modele.getPointsDeControle(a1,a2)
                axs[i][j].plot([p[0] for p in points],[p[1] for p in points])
                axs[i][j].set_title("activités "+str((a1,a2)))
        plt.show()
        
    def drawTargetsMap(file,constellation,details=True):
        def couleur(type_req):
            col_type = {}
            col_type["ONE_SHOT_MONO"] = 'k'
            col_type['ONE_SHOT_STEREO'] = 'b'
            col_type['DOWNLOAD'] = 'r'
            col_type['PERIODIC'] ='Orange'
            col_type['SYSTEMATIC'] = 'm'
            col_type['LONG_MONO'] = 'c'
            return col_type[type_req]
        
        if MPI.COMM_WORLD.Get_rank()==0:
            df = pd.DataFrame(columns=["Requête","Latitude","Longitude"])
            for r in constellation.getRequetes():
                xyz = constellation.getRequete(r).getMeanTarget(constellation)
                df = df.append({"Requête":int(r),'Latitude':xyz[0],"Longitude":xyz[1]},True)
            df['Requête'] = df['Requête'].astype(int)
            
            xyz,minP,maxP,N = constellation.statsCoordonnees()
            #print("min coord",minP,"max coord",maxP)
            gdf = gpd.GeoDataFrame(
                df, geometry=gpd.points_from_xy(df.Longitude, df.Latitude))
            #gdf['colour'] = [couleur(constellation,r) for r in gdf['Requête']]
            gdf['type'] = [constellation.getRequete(r).getType() for r in gdf['Requête']]
            world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
            ax = world[world.continent == 'Europe'].plot(
                color='#b9f9b2', edgecolor='black',alpha=0.5) 
            ax.set_facecolor('#A8C5DD')
            ax.set_xlim((minP[1]-10,maxP[1]+10))
            ax.set_ylim((minP[0]-10,maxP[0]+10))
            
            if details:
                for type_req, data in gdf.groupby('type'):
                    data.plot(ax=ax,label=type_req,color=couleur(type_req),categorical=True,legend=True)
                plt.legend()
            else:
                gdf.plot(ax=ax,color='black',marker='.',label='target')
                plt.legend(fontsize=30,loc='upper right')
            printColor('Info points :',c='r')
            printColor("| "+str(len(gdf))+" points")
            printColor("| "+str(len(pd.unique(gdf["geometry"])))+" valeurs différentes")
            plt.show()
    
    def drawCCAs(filename,constellation):
        def loopCCAJob(constellation,liste_cca,cca_working,grapheDependances):
            done = False
            info_cca = {}
            while not done:
                with Locker():
                    found = False
                    for i in range(len(cca_working)):
                        if not cca_working[i]:
                            cca_working[i] = 1
                            found = True
                            cca = liste_cca[i]
                            traiterCCA(constellation,info_cca,cca,grapheDependances)
                            break
                if not found:
                    done = True
            return info_cca
                    
        def fusionnerInfoCCA(info_cca):
            if MPI.COMM_WORLD.Get_rank()==0:
                for i in range(1,MPI.COMM_WORLD.Get_size()):
                    data = MPI.COMM_WORLD.recv()
                    for cca in data:
                        assert(cca not in info_cca)
                        info_cca[cca] = data[cca]
            else:
                MPI.COMM_WORLD.send(info_cca,dest=0)
            MPI.COMM_WORLD.Barrier()
            
        def traiterCCA(constellation,info_cca,cca,grapheDependances):
            debut = np.Inf
            fin = -np.Inf
            minLat,minLong = np.Inf,np.Inf
            maxLat,maxLong = -np.Inf,-np.Inf
            liste_lat = []
            liste_long = []
            s = grapheDependances.getSatelliteCCA(cca)
            for a in grapheDependances.getComposanteConnexe(cca):
                activite = constellation.getActivite(a)
                debut = min(debut,activite.getDebut())
                fin = max(activite.getFin(),fin)
                minLat = min(minLat,activite.getCoordonnees()[0])
                minLong = min(minLong,activite.getCoordonnees()[1])
                maxLat = max(maxLat,activite.getCoordonnees()[0])
                maxLong = max(maxLong,activite.getCoordonnees()[1])
                liste_lat.append(activite.getCoordonnees()[0])
                liste_long.append(activite.getCoordonnees()[1])
            N = len(liste_lat)
            if N==1:
                info_cca[cca] = {"sat":s,"cca":cca,"time":(debut,fin),"shape":Point(liste_long[0],liste_lat[0])}
            elif N==2:
                info_cca[cca] = {"sat":s,"cca":cca,"time":(debut,fin),"shape":LineString(zip(liste_long,liste_lat))}
            else:
                liste_coord = [Point(x,y) for (x,y) in zip(liste_long,liste_lat)]
                polygon_geom = Polygon(liste_coord).convex_hull
                info_cca[cca] = {"sat":s,"cca":cca,"time":(debut,fin),"shape":polygon_geom}
            
        def initListeCCAATraite(grapheDependances):
            liste_cca = grapheDependances.getComposantes()
            """ init une liste partagee de booleen sur les cca restantes a traiter """
            size = len(liste_cca) 
            itemsize = MPI.BOOL.Get_size() 
            if MPI.COMM_WORLD.Get_rank() == 0: 
                nbytes = size * itemsize 
            else: 
                nbytes = 0
            win = MPI.Win.Allocate_shared(nbytes, itemsize, comm=MPI.COMM_WORLD) 
            # create a numpy array whose data points to the shared mem
            cca_working, itemsize = win.Shared_query(0)
            if MPI.COMM_WORLD.Get_rank()==0:
                for i in range(size):
                    cca_working[i] = 0
            return liste_cca,cca_working
            
        def buildInfoCCA(constellation):
            activites = constellation.extraireActivitesRequetes()
            printMaster("Création des composantes connexes ...")
            grapheDependances = GroupeComposantesActivitesStatiques(constellation,activites)
            printMaster(grapheDependances,c='b')
            liste_cca,cca_working = initListeCCAATraite(grapheDependances)
            info_cca = loopCCAJob(constellation,liste_cca,cca_working,grapheDependances)
            fusionnerInfoCCA(info_cca)
            return info_cca
        
        
        def couleurTime(debut,fin,startHorizon,endHorizon):
            middle = ((debut+fin)/2-startHorizon)/(endHorizon-startHorizon)
            assert(middle>=0)
            assert(middle <= 1)
            gauche = np.array([0,0,1])
            droite = np.array([1,0,0])
            return middle*droite + (1-middle)*gauche
        
        def couleurSat(sat,Nsat):
            t = sat/Nsat
            assert(t>=0)
            assert(t <= 1)
            gauche = np.array([0.0,0.0,1.0])
            droite = np.array([1.0,0.0,0.0])
            return t*droite + (1-t)*gauche
        def extractDates(times):
            dates = []
            times = sorted([round(x[0],2) for x in times]+[round(x[1],2) for x in times])
            for t in times:
                if len(dates)==0 or t>dates[-1]:
                    dates.append(t)
            return dates
        def plotDate(file,constellation,T,fig,world,gdf,ax,Xlim,Ylim,Ncca):
            ax.cla()
            world[world.continent == 'Europe'].plot(ax=ax,color='green', edgecolor='black',alpha=0.7)
            ax.set_facecolor('#A8C5DD')
            testInt = lambda inter : inter[0]<= T and inter[1]>=T
            gdf_plot = gdf[gdf["time"].apply(testInt)]
            if len(gdf_plot)>0:
                gdf_plot.plot(alpha=0.4,ax=ax ,edgecolor='k',color=gdf_plot['couleur'])
                plt.legend(["satellite " + str(sat) for sat in gdf_plot['sat']],loc='upper right')
            ax.set_xlim(Xlim)
            ax.set_ylim(Ylim)
            ax.set_title(file.split('/')[-1] + " - time : "+str(round(T)) + " (s) - " + str(len(gdf_plot)) +"/" + str(Ncca) +" CCA" )
                    
        def initSlider(fig,dates):
            fig.subplots_adjust(bottom=0.25)
            fig.set_figheight(5)
            fig.set_figwidth(5)
            # Make a horizontal slider to control the frequency.
            axtime = fig.add_axes([0.25,0.1,0.65,0.03])
            slider = Slider(
                ax=axtime,
                label='time',
                valmin=dates[0],
                valmax=dates[-1],
                valinit=dates[0],
            )
            return slider
    
        def get_cmap(n, name='hsv'):
            '''Returns a function that maps each index in 0, 1, ..., n-1 to a distinct 
            RGB color; the keyword argument name must be a standard mpl colormap name.'''
            return plt.cm.get_cmap(name, n)
    
        info_cca = buildInfoCCA(constellation)
        if MPI.COMM_WORLD.Get_rank()==0:
            df = pd.DataFrame(columns=["cca","time","shape"])
            for cca in info_cca:
                df = df.append(info_cca[cca],ignore_index=True)
            df['cca'] = df['cca'].astype(int)
            df['sat'] = df['sat'].astype(int)
            cmap = get_cmap(len(df['sat'].value_counts()))
            Ncca = len(df['cca'])
            gdf = gpd.GeoDataFrame(df, geometry=df["shape"])
            Nsat = len(list(constellation.getSatellites().keys()))
            gdf['couleur'] = gdf["sat"].apply(lambda s : cmap(s))
            world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
    
            start = min(x[0] for x in df['time'])
            end = max(x[1] for x in df['time'])
            dates = extractDates(df["time"])
    
            fig,ax = plt.subplots()
            xyz,minP,maxP,N = constellation.statsCoordonnees()
            Xlim = (minP[1]-10,maxP[1]+10)
            Ylim = (minP[0]-10,maxP[0]+10)
            
            slider = initSlider(fig,dates)
            plotDate(filename,constellation,dates[0],fig,world,gdf,ax,Xlim,Ylim,Ncca)
            slider.on_changed(lambda var : plotDate(filename,constellation,var,fig,world,gdf,ax,Xlim,Ylim,Ncca))
            plt.show()
            
    def drawObsDurations(file,constellation):
        if MPI.COMM_WORLD.Get_rank()==0:
            durees = {}
            for r in constellation.getRequetes():
                type_req = constellation.getRequete(r).getType()
                if type_req not in durees:
                    durees[type_req] = []
                for (s,a) in constellation.getRequete(r).getCouples():
                    duree = constellation.getSatellite(s).getActivite(a).getDuree()
                    durees[type_req].append(duree)
            nAxis = len(list(durees.keys()))
            f,axs = plt.subplots(nrows=nAxis)
            for i,type_req in enumerate(durees):
                n, bins, patches = axs[i].hist(durees[type_req], density=False, facecolor='g', alpha=0.75)
                axs[i].set_title(type_req)
            plt.show()
    
    def choisirEtExecuterOption(file,constellation,choix):
        if MPI.COMM_WORLD.Get_rank()==0:
            retour = False
            while not retour:
                choix_retenu = None
                while(choix_retenu is None):
                    try:
                        liste_choix = []
                        for i,text in enumerate(choix):
                            print(str(i)+") "+text)
                            liste_choix.append(text)
                        x = int(input())
                        if x < len(liste_choix) and x >= 0:
                            choix_retenu = liste_choix[x]
                        else:
                            raise ValueError('Erreur input')
                    except Exception as e:
                        printColor(e,c='r')
                if choix[choix_retenu] is None:
                    retour = True
                    data = {"stop":True}
                    MPI.COMM_WORLD.bcast(data,root=0)
                else:
                    data = {"choix":choix_retenu}
                    MPI.COMM_WORLD.bcast(data,root=0)
                    choix[choix_retenu](file,constellation)
        else:
            while True:
                data = None
                data = MPI.COMM_WORLD.bcast(data,root=0)
                if 'stop' in data:
                    break
                else:
                    choix[data['choix']](file,constellation)
                    
        
    def initVisualisation():    
        choix = {}
        choix["afficher la carte des cibles au sol"] = lambda f,c : drawTargetsMap(f,c,details=False)
        choix["afficher la carte des cibles au sol (détails en couleur)"] = lambda f,c : drawTargetsMap(f,c,details=True)
        choix["afficher la carte des composantes connexes (/!\ lent)"] = drawCCAs
        choix["afficher les scores des observations (/!\ lent)"] = drawScoreMap
        choix["afficher les durées des activités"] = drawObsDurations
        if config.isDataTimeDependent():
            choix["afficher les durées de transition"] = drawTransition
        choix["retour"] = None

        
        rank = MPI.COMM_WORLD.Get_rank()
        while(True):
            if rank==0:
                file = choseAndBroadcastFile(path,instance)
                constellation = Constellation(file,True)
                MPI.COMM_WORLD.Barrier()
            else:
                data = None
                data = MPI.COMM_WORLD.bcast(data,root=0)
                file = data['files'][0]
                constellation = Constellation(file,False)
                MPI.COMM_WORLD.Barrier()
            choisirEtExecuterOption(file,constellation,choix)
    
    initVisualisation()
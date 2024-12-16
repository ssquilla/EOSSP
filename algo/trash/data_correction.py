#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 13:59:35 2023

@author: ssquilla
"""
import os


path = "../data/time-independent"

def corrigerNumerotationActivites(filename):
    print("correcting "+filename+"...")
    new_content = []
    i = 0
    indices_lignes_a_changer = {}
    with open(filename,'r') as file:
        first = file.readline()
        i += 1
        new_content.append(first)
        Nrequetes = int(first.split(",")[0])
        # skip requests
        for r in range(Nrequetes):
            content = file.readline()
            i += 1
            new_content.append(content)
            line = content.split("\n")[0].split(",")
            idRequete,nObs,type_requete = int(line[0]),int(line[1]),line[3]
            for o in range(nObs):
                content = file.readline()
                if type_requete in ["ONE_SHOT_MONO","LONG_MONO","SYSTEMATIC"]:
                    indices_lignes_a_changer[i] = 0
                else:
                    assert(type_requete in ["ONE_SHOT_STEREO","PERIODIC"])
                    indices_lignes_a_changer[i] = 1
                i += 1
                new_content.append(content)
        download = file.readline()
        new_content.append(download)
        nDownloads = int(download.split("\n")[0])
        for l in range(nDownloads):
            new_content.append(file.readline())
  
    with open(filename,'w') as file:
        for i in range(len(new_content)):
            if i in indices_lignes_a_changer:
                line = new_content[i].split("\n")[0].split(",")
                new_line = ""
                indice = indices_lignes_a_changer[i]
                line[indice] = str(int(line[indice])+nDownloads)
                new_content[i] = ""
                for j,elmt in enumerate(line):
                    new_line += elmt
                    if j < len(line)-1:
                        new_line += ","
                    else:
                        new_line += "\n"
                file.write(new_line)
            else:
                file.write(new_content[i])
            i += 1
        
for root, dirs, files in os.walk(path):
   for name in files:
      corrigerNumerotationActivites(root+"/"+name)


            
# petite manip pour les imports relatifs
import sys, importlib
from pathlib import Path
from resource import *
from functools import partial

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
Created on Mon Jul 10 18:10:38 2023

@author: ssquilla
"""

import tkinter as tk
from tkinter import *

import subprocess
import os

from .Utils.config import *

config = Config()

def callback():
    print("Let's go !")

def updateMenuText(menu,variables,opt,value):
    variables[opt] = value
    menu.title = opt+" : "+str(value)
    print("mise a jour d'option ",opt,":",value)
    

def printInfo(variables,opt):
    varInt = variables[opt]
    print(opt + ":" + str(varInt.get()))
    
def MAJcheckbox(variables,opt):
    variables[opt] = varInt.get()
    print(opt + ":" + str(varInt.get()))

class SchedUIler(tk.Tk):
    def __init__(self,width,height,marginright=30):
        tk.Tk.__init__(self)
        
        self.title("SchedUIler")
        self.width = width
        self.height= height
        self.marginright = marginright
        self.geometry(str(self.width)+"x"+str(self.height))
        root_paned = PanedWindow(self, orient=VERTICAL)
        root_paned.pack(side=TOP, expand=Y, fill=BOTH, pady=2, padx=2)
        self.createTopPanel(root_paned)
        self.createBottomPanel(root_paned)
        #r.pack()

    def callback(self,x):
        print(x)
        
    def buildOptionsList(self,top_paned):
        liste = Frame(top_paned,width=self.marginright, background='gray')
        variables = {}
        for opt in config.options:
            value = config.getOptValue(opt)
            domain = config.options[opt].getDomain()
            if isinstance(domain, BooleanDomain):
                variables[opt] = IntVar()
                cb = Checkbutton(liste, text=opt,variable=variables[opt], command=partial(printInfo,variables,opt),bg='gray', anchor="w",padx=0)
                cb.pack()
                
            elif isinstance(domain,StrDomain):
                if domain.getValues() is not None:
                    menuFormat = Menubutton(liste, text=opt+":"+value, width='20', borderwidth=2, bg='gray', activebackground='darkorange',relief = RAISED)
                    menu = Menu(menuFormat,tearoff=0)
                    updateMenuText(menu,variables,opt,value)

                    for val in domain.getValues():
                        menu.add_command(label=val,command=partial(updateMenuText,menu,variables,opt,val))
                    menuFormat.configure(menu=menu)
                    menuFormat.pack()
                else:
                    print("domain None",opt)
                #rb = Radiobutton(liste,text=opt,value=0)
                #rb.pack()
            #liste.insert(i,opt)
        liste.pack(side=RIGHT,fill=BOTH,expand=False)

    def buildTerminal(self,top_paned):
        self.Terminal = Frame(top_paned,width=100-self.marginright)
        self.Terminal.pack(side=LEFT,fill = BOTH, expand = True)
        winfo = self.Terminal.winfo_id()
        os.system("xterm -into %d -maximized &" % winfo) # -sb elimine la scrollbar, %d current dir -maximized
        
    def createTopPanel(self,root_paned):
        top_paned = PanedWindow(root_paned, orient=HORIZONTAL)#,height=w_height-bottom_margin)
        top_paned.pack(side=TOP, expand=Y, fill=BOTH, pady=2, padx=2)
        self.buildTerminal(top_paned)
        self.buildOptionsList(top_paned)

        return top_paned

    def createBottomPanel(self,root_paned):
        bottom_paned = PanedWindow(root_paned, height="1cm", orient=HORIZONTAL)
        bottom_paned.pack(side=BOTTOM, expand=False,fill=BOTH) #, pady=2, padx=2
        start_button = Button(bottom_paned,text='Démarrer', command=callback)
        bottom_paned.add(start_button)
        return bottom_paned



scheduiler = SchedUIler(1200,500)
scheduiler.mainloop()
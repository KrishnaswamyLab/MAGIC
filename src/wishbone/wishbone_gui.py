#!/usr/local/bin/python3

import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle
from matplotlib.path import Path
from functools import reduce
import wishbone
import os
import sys
import platform
import pandas as pd
import tkinter as tk
import numpy as np
from tkinter import filedialog
import pickle
import random

class wishbone_gui(tk.Tk):
    def __init__(self,parent):
        tk.Tk.__init__(self,parent)
        self.parent = parent
        self.initialize()

    def initialize(self):
        self.grid()
        self.vals = None
        self.currentPlot = None
        self.data = {}
        self.scseq = False

        #set up menu bar
        self.menubar = tk.Menu(self)
        self.fileMenu = tk.Menu(self.menubar, tearoff=0)
        self.menubar.add_cascade(label="File", menu=self.fileMenu)
        self.fileMenu.add_command(label="Load data", command=self.loadData)
        self.fileMenu.add_command(label="Save data", state='disabled', command=self.saveData)
        self.fileMenu.add_command(label="Exit Wishbone", command=self.quitWB)

        self.analysisMenu = tk.Menu(self.menubar, tearoff=0)
        self.menubar.add_cascade(label="Analysis", menu=self.analysisMenu)
        self.analysisMenu.add_command(label="Principal component analysis", state='disabled', command=self.runPCA)
        self.analysisMenu.add_command(label="tSNE", state='disabled', command=self.runTSNE)
        self.analysisMenu.add_command(label="Diffusion map", state='disabled', command=self.runDM)
        self.analysisMenu.add_command(label="GSEA", state='disabled', command=self.runGSEA)
        self.analysisMenu.add_command(label="Wishbone", state='disabled', command=self.runWishbone)
        self.analysisMenu.add_command(label="Magic", state='disabled', command=self.runMagic)

        self.visMenu = tk.Menu(self.menubar, tearoff=0)
        self.menubar.add_cascade(label="Visualization", menu=self.visMenu)
        self.visMenu.add_command(label="Principal component analysis", state='disabled', command=self.plotPCA)
        self.visMenu.add_command(label="tSNE", state='disabled', command=self.plotTSNE)
        self.visMenu.add_command(label="Diffusion map", state='disabled', command=self.plotDM)
        self.visMenu.add_command(label="GSEA Results", state='disabled', command=self.showGSEAResults)
        self.wishboneMenu = tk.Menu(self)
        self.visMenu.add_cascade(label="Wishbone", menu=self.wishboneMenu)
        self.wishboneMenu.add_command(label="On tSNE", state='disabled', command=self.plotWBOnTsne)
        self.wishboneMenu.add_command(label="Marker trajectory", state='disabled', command=self.plotWBMarkerTrajectory)
        self.wishboneMenu.add_command(label="Heat map", state='disabled', command=self.plotWBHeatMap)
        self.geneExpMenu = tk.Menu(self)
        self.visMenu.add_cascade(label="Gene expression", menu=self.geneExpMenu)
        self.geneExpMenu.add_command(label="Scatter plot", state='disabled', command=self.scatterGeneExp)
        self.geneExpMenu.add_command(label="On tSNE", state='disabled', command=self.plotGeneExpOntSNE)
        self.geneExpMenu.add_command(label="Compare data sets", state='disabled', command=self.scatterGeneExpAgainstOther)
        self.visMenu.add_command(label="Set gate", state='disabled', command=self.setGate)
        
        self.config(menu=self.menubar)

        #intro screen
        tk.Label(self, text=u"Wishbone", font=('Helvetica', 48), fg="black", bg="white", padx=100, pady=50).grid(row=0)
        tk.Label(self, text=u"To get started, select a data file by clicking File > Load Data", fg="black", bg="white", padx=100, pady=25).grid(row=1)

        #update
        self.protocol('WM_DELETE_WINDOW', self.quitWB)
        self.grid_columnconfigure(0,weight=1)
        self.resizable(True,True)
        self.update()
        self.geometry(self.geometry())       
        self.focus_force()

    def loadData(self):
        self.dataFileName = filedialog.askopenfilename(title='Load data file', initialdir='~/.wishbone/data')
        if(self.dataFileName != ""):
            #pop up data options menu
            self.fileInfo = tk.Toplevel()
            self.fileInfo.title("Data options")
            tk.Label(self.fileInfo, text=u"File name: ").grid(column=0, row=0)
            tk.Label(self.fileInfo, text=self.dataFileName.split('/')[-1]).grid(column=1, row=0)

            tk.Label(self.fileInfo,text=u"Name:" ,fg="black",bg="white").grid(column=0, row=1)
            self.fileNameEntryVar = tk.StringVar()
            tk.Entry(self.fileInfo, textvariable=self.fileNameEntryVar).grid(column=1,row=1)

            self.dataFileType = self.dataFileName.split('.', 1)[1]
            if  self.dataFileType == 'fcs':
                tk.Label(self.fileInfo,text=u"Cofactor:" ,fg="black",bg="white").grid(column=0, row=2)
                self.cofactorVar = tk.IntVar()
                self.cofactorVar.set(5)
                tk.Entry(self.fileInfo, textvariable=self.cofactorVar).grid(column=1,row=2)
            elif self.dataFileType == 'csv':
                self.normalizeVar = tk.BooleanVar()
                tk.Checkbutton(self.fileInfo, text=u"Normalize", variable=self.normalizeVar).grid(column=0, row=2, columnspan=2)
                tk.Label(self.fileInfo, text=u"The normalize parameter is used for correcting for library size among cells.").grid(column=0, row=3, columnspan=2)
            elif self.dataFileType == 'mtx' or self.dataFileType == 'mtx.gz':
                #get gene names file
                tk.Button(self.fileInfo, text="Select gene names file", command=self.getGeneNameFile).grid(column=0, row=2)

                #filter parameters
                self.filterVar = tk.BooleanVar()
                tk.Checkbutton(self.fileInfo, text=u"Filter cells/genes", variable=self.filterVar).grid(column=0, row=3, columnspan=4)
                self.filterCellMinVar = tk.IntVar()
                self.filterCellMinVar.set(0)
                tk.Label(self.fileInfo,text=u"Filter by reads per cell. Min:" ,fg="black",bg="white").grid(column=0, row=4)
                tk.Entry(self.fileInfo, textvariable=self.filterCellMinVar).grid(column=1,row=4)
                self.filterCellMaxVar = tk.IntVar()
                self.filterCellMaxVar.set(0)
                tk.Label(self.fileInfo,text=u" Max:" ,fg="black",bg="white").grid(column=2, row=4)
                tk.Entry(self.fileInfo, textvariable=self.filterCellMaxVar).grid(column=3,row=4)
                self.filterGeneNonzeroVar = tk.IntVar()
                self.filterGeneNonzeroVar.set(0)
                tk.Label(self.fileInfo,text=u"Filter by nonzero cells per gene. Min:" ,fg="black",bg="white").grid(column=0, row=5)
                tk.Entry(self.fileInfo, textvariable=self.filterGeneNonzeroVar).grid(column=1,row=5)
                self.filterGeneMolsVar = tk.IntVar()
                self.filterGeneMolsVar.set(0)
                tk.Label(self.fileInfo,text=u"Filter by reads per gene. Min:" ,fg="black",bg="white").grid(column=0, row=6)
                tk.Entry(self.fileInfo, textvariable=self.filterGeneMolsVar).grid(column=1,row=6)

                #normalize
                self.normalizeVar = tk.BooleanVar()
                tk.Checkbutton(self.fileInfo, text=u"Normalize", variable=self.normalizeVar).grid(column=0, row=7, columnspan=4)
                tk.Label(self.fileInfo, text=u"The normalize parameter is used for correcting for library size among cells.").grid(column=0, row=8, columnspan=4)

            tk.Button(self.fileInfo, text="Cancel", command=self.fileInfo.destroy).grid(column=1, row=9)
            tk.Button(self.fileInfo, text="Load", command=self.processData).grid(column=2, row=9)

            self.wait_window(self.fileInfo)

    def getGeneNameFile(self):
        self.geneNameFile = filedialog.askopenfilename(title='Select gene name file', initialdir='~/.wishbone/data')
        tk.Label(self.fileInfo,text=self.geneNameFile.split('/')[-1] ,fg="black",bg="white").grid(column=1, row=2)

    def processData(self):

        if len(self.data) == 0:
            #clear intro screen
            for item in self.grid_slaves():
                item.grid_forget()

            tk.Label(self,text=u"Data sets:" , fg="black",bg="white").grid(column=0,row=0)

            #set up canvas for plots
            self.fig, self.ax = wishbone.wb.get_fig()
            self.canvas = FigureCanvasTkAgg(self.fig, self)
            self.canvas.show()
            self.canvas.get_tk_widget().grid(column=1, row=1, rowspan=10, columnspan=4,sticky='NSEW')

        #load data based on input type
        if self.dataFileType == 'fcs':    # mass cytometry data
            scdata = wishbone.wb.SCData.from_fcs(os.path.expanduser(self.dataFileName), 
                                                      cofactor=self.cofactorVar.get())
            wb = None
        elif self.dataFileType == 'csv':    # sc-seq data
            scdata = wishbone.wb.SCData.from_csv(os.path.expanduser(self.dataFileName), data_type='sc-seq', 
                                                      normalize=self.normalizeVar.get())
            wb = None
        elif self.dataFileType == 'mtx' or self.dataFileType == 'mtx.gz':   # sparse matrix
            if self.filterVar.get() == True:
                scdata = wishbone.wb.SCData.from_mtx(os.path.expanduser(self.dataFileName), os.path.expanduser(self.geneNameFile),
                                                     filter_cell_min=self.filterCellMinVar.get(), filter_cell_max=self.filterCellMaxVar.get(), 
                                                     filter_gene_nonzero=self.filterGeneNonzeroVar.get(),
                                                     filter_gene_mols=self.filterGeneMolsVar.get(), normalize=self.normalizeVar.get())
            else:
                scdata = wishbone.wb.SCData.from_mtx(os.path.expanduser(self.dataFileName), os.path.expanduser(self.geneNameFile),
                                                          normalize=self.normalizeVar.get())   
            wb = None
        else:   # pickled Wishbone object
            wb = wishbone.wb.Wishbone.load(self.dataFileName)
            scdata = wb.scdata

        self.data[self.fileNameEntryVar.get()] = {'scdata' : scdata, 'wb' : wb, 'state' : tk.BooleanVar(),
                                                  'genes' : scdata.data.columns.values, 'gates' : {}}
        for i in range(len(self.data)):
            tk.Checkbutton(self, text=list(self.data.keys())[i], variable=list(self.data.values())[i]['state']).grid(column=0, row=i+1)
        
        if self.scseq == False and scdata.data_type == 'sc-seq':
            self.scseq = True

        #set up buttons based on data type
        if self.scseq == True:
            self.saveButton = tk.Button(self, text=u"Save plot", state='disabled', command=self.savePlot)
            self.saveButton.grid(column = 3, row=0)
            self.diff_component = tk.StringVar()
            self.diff_component.set('Component 1')
            self.component_menu = tk.OptionMenu(self, self.diff_component,
                                                'Component 1', 'Component 2', 'Component 3',
                                                'Component 4', 'Component 5', 'Component 6', 
                                                'Component 7', 'Component 8', 'Component 9')
            self.component_menu.config(state='disabled')
            self.component_menu.grid(row=0, column=1)
            self.updateButton = tk.Button(self, text=u"Update component", command=self.updateComponent, state='disabled')
            self.updateButton.grid(column=2, row=0)
            self.analysisMenu.entryconfig(5, state='normal')

            #enable buttons based on current state of scdata object
            if scdata.pca:
                self.analysisMenu.entryconfig(1, state='normal')
                self.visMenu.entryconfig(0, state='normal')
            if isinstance(scdata.tsne, pd.DataFrame):
                self.analysisMenu.entryconfig(2, state='normal')
                self.visMenu.entryconfig(1, state='normal')
                self.geneExpMenu.entryconfig(1, state='normal')
            if isinstance(scdata.diffusion_eigenvectors, pd.DataFrame):
                self.analysisMenu.entryconfig(3, state='normal')
                self.analysisMenu.entryconfig(4, state='normal')
                self.visMenu.entryconfig(2, state='normal')
        else:
            self.saveButton = tk.Button(self, text=u"Save plot", state='disabled', command=self.savePlot)
            self.saveButton.grid(column = 4, row=0)

            self.analysisMenu.delete(0)
            self.analysisMenu.delete(2)
            self.visMenu.delete(0)
            self.visMenu.delete(2)
            self.analysisMenu.entryconfig(1, state='normal')
            self.analysisMenu.entryconfig(3, state='normal')

            #enable buttons based on current state of scdata object
            if isinstance(scdata.tsne, pd.DataFrame):
                self.visMenu.entryconfig(0, state='normal')
                self.geneExpMenu.entryconfig(1, state='normal')
            if isinstance(scdata.diffusion_eigenvectors, pd.DataFrame):
                self.analysisMenu.entryconfig(2, state='normal')
                self.visMenu.entryconfig(1, state='normal')

        #enable buttons
        self.analysisMenu.entryconfig(0, state='normal')
        self.fileMenu.entryconfig(1, state='normal')
        self.geneExpMenu.entryconfig(0, state='normal')
        self.concatButton = tk.Button(self, text=u"Concatenate selected datasets", state='disabled', wraplength=80, command=self.concatenateData)
        self.concatButton.grid(column=0, row=10)
        if len(self.data) > 1:
            self.concatButton.config(state='normal')
            self.geneExpMenu.entryconfig(2, state='normal')

        if wb:
                if isinstance(wb.trajectory, pd.Series):
                    self.wishboneMenu.entryconfig(0, state='normal')
                    self.wishboneMenu.entryconfig(1, state='normal')
                    self.wishboneMenu.entryconfig(2, state='normal')

        self.geometry('800x550')
        #destroy pop up menu
        self.fileInfo.destroy()

    def saveData(self):
        for key in self.data:
            if self.data[key]['state'].get() == True:
                pickleFileName = filedialog.asksaveasfilename(title=key + ': save wishbone data', defaultextension='.p', initialfile=key)
                if pickleFileName != None:
                    if self.data[key]['wb'] != None:
                        self.data[key]['wb'].save(pickleFileName)
                    else:
                        self.data[key]['scdata'].save_as_wishbone(pickleFileName)

    def concatenateData(self):
        self.concatOptions = tk.Toplevel()
        self.concatOptions.title("Concatenate data sets")

        tk.Label(self.concatOptions, text=u"New data set name:", fg="black", bg="white").grid(column=0, row=0)
        self.nameVar = tk.StringVar()
        tk.Entry(self.concatOptions, textvariable=self.nameVar).grid(column=1, row=0)

        self.joinVar = tk.BooleanVar()
        self.joinVar.set(True)
        tk.Checkbutton(self.concatOptions, text=u"Outer join", variable=self.joinVar).grid(column=0, row=1, columnspan=2)

        tk.Button(self.concatOptions, text="Concatenate", command=self._concatenateData).grid(column=1, row=2)
        tk.Button(self.concatOptions, text="Cancel", command=self.concatOptions.destroy).grid(column=0, row=2)
        self.wait_window(self.concatOptions)

    def _concatenateData(self):
        to_concat = []
        for key in self.data:
            if self.data[key]['state'].get() == True:
                to_concat.append(self.data[key]['scdata'])

        scdata = to_concat[0].concatenate_data(to_concat[1:], 
                                               join='outer' if self.joinVar.get() == True else 'inner')

        self.data[self.nameVar.get()] = {'scdata' : scdata, 'wb' : None, 'state' : tk.BooleanVar(),
                                         'genes' : scdata.data.columns.values, 'gates' : {}}
        
        #update data sets
        for i in range(len(self.data)):
            tk.Checkbutton(self, text=list(self.data.keys())[i], variable=list(self.data.values())[i]['state']).grid(column=0, row=i+1)
        
        self.concatOptions.destroy()

    def runPCA(self):
        for key in self.data:
            if self.data[key]['state'].get() == True:
                self.data[key]['scdata'].run_pca()

        #enable buttons
        self.analysisMenu.entryconfig(1, state='normal')
        self.visMenu.entryconfig(0, state='normal')

    def runTSNE(self):
        for key in self.data:
            if self.data[key]['state'].get() == True:
                #pop up for # components
                self.tsneOptions = tk.Toplevel()
                self.tsneOptions.title(key + ": tSNE options")
                self.curKey = key
                if self.data[key]['scdata'].data_type == 'sc-seq':
                    tk.Label(self.tsneOptions,text=u"Number of components:" ,fg="black",bg="white").grid(column=0, row=0)
                    self.nCompVar = tk.IntVar()
                    self.nCompVar.set(15)
                    tk.Entry(self.tsneOptions, textvariable=self.nCompVar).grid(column=1,row=0)
                tk.Label(self.tsneOptions,text=u"Perplexity:" ,fg="black",bg="white").grid(column=0, row=1)
                self.perplexityVar = tk.IntVar()
                self.perplexityVar.set(30)
                tk.Entry(self.tsneOptions, textvariable=self.perplexityVar).grid(column=1,row=1)
                tk.Button(self.tsneOptions, text="Run", command=self._runTSNE).grid(column=1, row=2)
                tk.Button(self.tsneOptions, text="Cancel", command=self.tsneOptions.destroy).grid(column=0, row=2)
                self.wait_window(self.tsneOptions)

    def _runTSNE(self):
        print(self.curKey)
        if self.data[self.curKey]['scdata'].data_type == 'sc-seq':
            self.data[self.curKey]['scdata'].run_tsne(n_components=self.nCompVar.get(), perplexity=self.perplexityVar.get())
        else:
            self.data[self.curKey]['scdata'].run_tsne(n_components=None, perplexity=self.perplexityVar.get())
        self.data[self.curKey]['gates'] = {}

        #enable buttons
        if self.data[self.curKey]['scdata'].data_type == 'sc-seq':
            self.analysisMenu.entryconfig(2, state='normal')
            self.visMenu.entryconfig(1, state='normal')
        else:
            self.visMenu.entryconfig(0, state='normal')
        self.geneExpMenu.entryconfig(1, state='normal')
        self.tsneOptions.destroy()

    def runDM(self):
        for key in self.data:
            if self.data[key]['state'].get() == True:
                self.data[key]['scdata'].run_diffusion_map()

            #enable buttons
            if self.data[key]['scdata'].data_type == 'sc-seq':
                self.analysisMenu.entryconfig(3, state='normal')
                self.analysisMenu.entryconfig(4, state='normal')
                self.visMenu.entryconfig(2, state='normal')
            else:
                self.analysisMenu.entryconfig(2, state='normal')
                self.visMenu.entryconfig(1, state='normal')

    def runGSEA(self):

        GSEAFileName = filedialog.askopenfilename(title='Select gmt File', initialdir='~/.wishbone/tools')

        if GSEAFileName != "":
            if 'mouse' in GSEAFileName:
                gmt_file_type = 'mouse'
            else:
                gmt_file_type = 'human'

            for key in self.data:
                if self.data[key]['state'].get() == True: 
                    self.data[key]['scdata'].run_diffusion_map_correlations()
                    self.data[key]['scdata'].data.columns = self.data[key]['scdata'].data.columns.str.upper()
                    outputPrefix = filedialog.asksaveasfilename(title=key + ': input file prefix for saving output', initialdir='~/.wishbone/gsea')
                    
                    self.data[key]['gsea_reports'] = self.data[key]['scdata'].run_gsea(output_stem= os.path.expanduser(outputPrefix), 
                                                                                       gmt_file=(gmt_file_type, GSEAFileName.split('/')[-1]))

            #enable buttons
            self.visMenu.entryconfig(3, state='normal')

    def runWishbone(self):
        for key in self.data:
            if self.data[key]['state'].get() == True: 
                self.curKey = key

                #popup menu for wishbone options
                self.wbOptions = tk.Toplevel()
                self.wbOptions.title(key + ": Wishbone options")

                #s
                tk.Label(self.wbOptions,text=u"Start cell:",fg="black",bg="white").grid(column=0,row=0)
                self.start = tk.StringVar()
                tk.Entry(self.wbOptions, textvariable=self.start).grid(column=1,row=0)
                if(len(self[key]['gates']) > 0):
                    self.cell_gate = tk.StringVar()
                    self.cell_gate.set('Use cell gate')
                    self.gate_menu = tk.OptionMenu(self.wbOptions, self.cell_gate,
                                                   *list(self.data[key]['gates'].keys()))
                    self.gate_menu.grid(row=0, column=2)

                #k
                tk.Label(self.wbOptions,text=u"k:",fg="black",bg="white").grid(column=0,row=1)
                self.k = tk.IntVar()
                tk.Entry(self.wbOptions, textvariable=self.k).grid(column=1,row=1)
                self.k.set(15)
                
                #components list
                tk.Label(self.wbOptions, text=u"Components list:", fg='black', bg='white').grid(column=0, row=2)
                self.compList = tk.StringVar()
                tk.Entry(self.wbOptions, textvariable=self.compList).grid(column=1, row=2)
                self.compList.set("1, 2, 3")

                #num waypoints
                tk.Label(self.wbOptions, text=u"Number of waypoints:", fg='black', bg='white').grid(column=0, row=3)
                self.numWaypoints = tk.IntVar()
                tk.Entry(self.wbOptions, textvariable=self.numWaypoints).grid(column=1, row=3)
                self.numWaypoints.set(250)

                #branch
                self.branch = tk.BooleanVar()
                self.branch.set(True)
                tk.Checkbutton(self.wbOptions, text=u"Branch", variable=self.branch).grid(column=0, row=4, columnspan=2)

                tk.Button(self.wbOptions, text="Run", command=self._runWishbone).grid(column=1, row=5)
                tk.Button(self.wbOptions, text="Cancel", command=self.wbOptions.destroy).grid(column=0, row=5)
                self.wait_window(self.wbOptions)

    def _runWishbone(self):
        self.data[self.curKey]['wb'] = wishbone.wb.Wishbone(self.data[self.curKey]['scdata'])

        if self.cell_gate.get() == 'Use cell gate':
            self.data[self.curKey]['wb'].run_wishbone(start_cell=self.start.get(), k=self.k.get(), 
                                                      components_list=[int(comp) for comp in self.compList.get().split(',')], 
                                                      num_waypoints=self.numWaypoints.get())
        else:
            #randomly select start cell in gate
            print('Using cell gate:')
            print(self.cell_gate.get())
            start_cell = random.sample(list(self.data[self.curKey]['gates'][self.cell_gate.get()]), 1)[0]
            print('Start cell: ' + start_cell)
            self.data[self.curKey]['wb'].run_wishbone(start_cell=start_cell, k=self.k.get(), 
                                                      components_list=[int(comp) for comp in self.compList.get().split(',')], 
                                                      num_waypoints=self.numWaypoints.get())
        
        #enable buttons
        self.wishboneMenu.entryconfig(0, state='normal')
        self.wishboneMenu.entryconfig(1, state='normal')
        self.wishboneMenu.entryconfig(2, state='normal')
        self.wbOptions.destroy()

    def runMagic(self):
        for key in list(self.data):
            if self.data[key]['state'].get() == True:
                #pop up for parameters
                self.magicOptions = tk.Toplevel()
                self.magicOptions.title(key + ": MAGIC options")
                self.curKey = key

                tk.Label(self.magicOptions,text=u"Magic-ed data name:" ,fg="black",bg="white").grid(column=0, row=0)
                self.nameVar = tk.StringVar()
                self.nameVar.set(key + ' Magic')
                tk.Entry(self.magicOptions, textvariable=self.nameVar).grid(column=1,row=0)

                tk.Label(self.magicOptions,text=u"# of PCA components:" ,fg="black",bg="white").grid(column=0, row=1)
                self.nCompVar = tk.IntVar()
                self.nCompVar.set(0)
                tk.Entry(self.magicOptions, textvariable=self.nCompVar).grid(column=1,row=1)

                tk.Label(self.magicOptions,text=u"t:" ,fg="black",bg="white").grid(column=0, row=2)
                self.tVar = tk.IntVar()
                self.tVar.set(8)
                tk.Entry(self.magicOptions, textvariable=self.tVar).grid(column=1,row=2)

                tk.Label(self.magicOptions,text=u"kNN:" ,fg="black",bg="white").grid(column=0, row=3)
                self.kNNVar = tk.IntVar()
                self.kNNVar.set(20)
                tk.Entry(self.magicOptions, textvariable=self.kNNVar).grid(column=1,row=3)

                tk.Label(self.magicOptions,text=u"Epsilon:" ,fg="black",bg="white").grid(column=0, row=4)
                self.epsilonVar = tk.IntVar()
                self.epsilonVar.set(0)
                tk.Entry(self.magicOptions, textvariable=self.epsilonVar).grid(column=1,row=4)

                self.rescaleVar = tk.BooleanVar()
                tk.Checkbutton(self.magicOptions, text=u"Rescale", variable=self.rescaleVar).grid(column=0, row=5, columnspan=2)

                tk.Button(self.magicOptions, text="Cancel", command=self.magicOptions.destroy).grid(column=0, row=6)
                tk.Button(self.magicOptions, text="Run", command=self._runMagic).grid(column=1, row=6)
                self.wait_window(self.magicOptions)

    def _runMagic(self):
        scdata = self.data[self.curKey]['scdata'].run_magic(n_pca_components=self.nCompVar.get() if self.nCompVar.get() > 0 else None,
                                                            t=self.tVar.get(), knn=self.kNNVar.get(), 
                                                            epsilon=self.epsilonVar.get(), rescale=self.rescaleVar.get())
        
        self.data[self.nameVar.get()] = {'scdata' : scdata, 'wb' : None, 'state' : tk.BooleanVar(),
                                         'genes' : scdata.data.columns.values, 'gates' : {}}
        
        #update data sets
        for i in range(len(self.data)):
            tk.Checkbutton(self, text=list(self.data.keys())[i], variable=list(self.data.values())[i]['state']).grid(column=0, row=i+1)
        
        self.magicOptions.destroy()

    def plotPCA(self):
        self.saveButton.config(state='normal')
        if self.scseq == True:
            self.component_menu.config(state='disabled')
            self.updateButton.config(state='disabled')
            self.visMenu.entryconfig(6, state='disabled')
        else:
            self.visMenu.entryconfig(4, state='disabled')

        #pop up for # components
        self.PCAOptions = tk.Toplevel()
        self.PCAOptions.title("PCA Plot Options")
        tk.Label(self.PCAOptions,text=u"Max variance explained (ylim):",fg="black",bg="white").grid(column=0, row=0)
        self.yLimVar = tk.DoubleVar()
        self.yLimVar.set(0.10)
        tk.Entry(self.PCAOptions, textvariable=self.yLimVar).grid(column=1,row=0)
        tk.Label(self.PCAOptions, text=u"Number of components:", fg='black', bg='white').grid(column=0, row=1)
        self.compVar = tk.IntVar()
        self.compVar.set(15)
        tk.Entry(self.PCAOptions, textvariable=self.compVar).grid(column=1, row=1)
        tk.Button(self.PCAOptions, text="Plot", command=self._plotPCA).grid(column=1, row=2)
        tk.Button(self.PCAOptions, text="Cancel", command=self.PCAOptions.destroy).grid(column=0, row=2)
        self.wait_window(self.PCAOptions)

    def _plotPCA(self):
        self.resetCanvas()
        keys = [key for key in self.data if self.data[key]['state'].get() == True]
        if len(keys) > 1:
            self.fig = plt.figure(figsize=[8, 4 * int(np.ceil(len(keys)/2))])
            gs = gridspec.GridSpec(int(np.ceil(len(keys)/2)), 2)
            for i in range(len(keys)):
                self.ax = self.fig.add_subplot(gs[int(i/2), i%2])
                self.data[keys[i]]['scdata'].plot_pca_variance_explained(fig=self.fig, ax=self.ax, ylim=(0, self.yLimVar.get()), n_components=self.compVar.get())
                self.ax.set_title(keys[i])
            gs.tight_layout(self.fig)

        else:
            self.fig, self.ax = self.data[keys[0]]['scdata'].plot_pca_variance_explained(ylim=(0, self.yLimVar.get()), n_components=self.compVar.get())
        
        self.canvas = FigureCanvasTkAgg(self.fig, self)
        self.canvas.show()
        self.canvas.get_tk_widget().grid(column=1, row=1, rowspan=10, columnspan=4) 
        self.currentPlot = 'pca'

        #enable buttons
        self.saveButton.config(state='normal')
        self.PCAOsudoptions.destroy()

    def plotTSNE(self):
        #pop up for color by
        self.TSNEOptions = tk.Toplevel()
        self.TSNEOptions.title("tSNE Plot Options")

        self.colorVar = tk.BooleanVar()
        tk.Checkbutton(self.TSNEOptions, text=u"Color by gene", variable=self.colorVar).grid(column=0, row=0, columnspan=2)
        
        self.data_set = tk.StringVar()
        self.data_set.set('Select data set')
        self.data_menu = tk.OptionMenu(self.TSNEOptions, self.data_set,
                                      *list(self.data.keys()))
        self.data_menu.grid(row=1, column=0)
        tk.Button(self.TSNEOptions, text="Use data set", command=self._updateGeneData).grid(column=1, row=1)

        tk.Label(self.TSNEOptions,text=u"Gene:",fg="black",bg="white").grid(column=0, row=2)
        tk.Entry(self.TSNEOptions).grid(column=1, row=2)

        tk.Button(self.TSNEOptions, text="Plot", command=self._plotTSNE).grid(column=1, row=5)
        tk.Button(self.TSNEOptions, text="Cancel", command=self.TSNEOptions.destroy).grid(column=0, row=5)
        self.wait_window(self.TSNEOptions)

    def _updateGeneData(self):
        self.geneInput = wishbone.autocomplete_entry.AutocompleteEntry(self.data[self.data_set.get()]['genes'], self.TSNEOptions, listboxLength=3)
        self.geneInput.grid(column=1, row=2)

    def _plotTSNE(self):
        keys = [key for key in self.data if self.data[key]['state'].get() == True]

        self.saveButton.config(state='normal')
        if self.scseq == True:
            self.component_menu.config(state='disabled')
            self.updateButton.config(state='disabled')
            if len(keys) == 1:
                self.visMenu.entryconfig(6, state='normal')
        else:
            if len(keys) == 1:
                self.visMenu.entryconfig(4, state='normal')

        if self.colorVar.get() == True:
            color = self.data[self.data_set.get()]['scdata'].data[self.geneInput.var.get()]
        else:
            color = None

        self.resetCanvas()
        if len(keys) > 1:
            self.fig = plt.figure(figsize=[8, 4 * int(np.ceil(len(keys)/2))])
            gs = gridspec.GridSpec(int(np.ceil(len(keys)/2)), 2)
            for i in range(len(keys)):
                self.ax = self.fig.add_subplot(gs[int(i/2), i%2])
                self.data[keys[i]]['scdata'].plot_tsne(fig=self.fig, ax=self.ax, color=color)
                self.ax.set_title(keys[i])
            gs.tight_layout(self.fig)

        else:
            self.fig, self.ax = self.data[keys[0]]['scdata'].plot_tsne(color=color)

        self.canvas = FigureCanvasTkAgg(self.fig, self)
        self.canvas.show()
        self.canvas.get_tk_widget().grid(column=1, row=1, rowspan=10, columnspan=4) 
        self.currentPlot = 'tsne'
        self.TSNEOptions.destroy()

    def plotDM(self):

        keys = [key for key in self.data if self.data[key]['state'].get() == True]
        if len(keys) > 2:
            self.DMError = tk.Toplevel()
            self.DMError.title("Error -- too many datasets selected")
            tk.Label(self.DMError,text=u"Please select up to two datasets to plot diffusion map components",fg="black",bg="white").grid(column=0, row=0)
            tk.Button(self.DMError, text="Ok", command=self.DMError.destroy).grid(column=0, row=1)
            self.wait_window(self.DMError)

        else:
            self.saveButton.config(state='normal')
            if self.scseq == True:
                self.component_menu.config(state='disabled')
                self.updateButton.config(state='disabled')
                self.visMenu.entryconfig(6, state='disabled')
            else:
                self.visMenu.entryconfig(4, state='disabled')

            self.geometry('950x550')
            self.resetCanvas()

            self.fig, self.ax = self.data[keys[0]]['scdata'].plot_diffusion_components(other_data=self.data[keys[1]]['scdata'] if len(keys) > 1 else None)
            for i in range(len(keys)):
                plt.figtext(0.05, 0.75 - (0.5*i), keys[i], rotation='vertical')

            self.canvas = FigureCanvasTkAgg(self.fig, self)
            self.canvas.show()
            self.canvas.get_tk_widget().grid(column=1, row=1, rowspan=10, columnspan=4,sticky='W') 
            self.currentPlot = 'dm_components'

    def showGSEAResults(self):
        keys = [key for key in self.data if self.data[key]['state'].get() == True]
        if len(keys) != 1:
            self.GSEAError = tk.Toplevel()
            self.GSEAError.title("Error -- too many datasets selected")
            tk.Label(self.GSEAError,text=u"Please select exactly one dataset to view GSEA results.",fg="black",bg="white").grid(column=0, row=0)
            tk.Button(self.GSEAError, text="Ok", command=self.GSEAError.destroy).grid(column=0, row=1)
            self.wait_window(self.GSEAError)
        else:
            self.curKey = keys[0]
            self.saveButton.config(state='disabled')
            self.component_menu.config(state='normal')
            self.updateButton.config(state='normal')
            self.visMenu.entryconfig(6, state='disabled')

            self.resetCanvas()
            self.canvas = tk.Canvas(self, width=600, height=300)
            self.canvas.grid(column=1, row=1, rowspan=17, columnspan=4)
            self.outputText(1)
            self.currentPlot = 'GSEA_result_'+self.diff_component.get()

    def updateComponent(self):
        self.resetCanvas()
        self.canvas = tk.Canvas(self, width=600, height=300)
        self.canvas.grid(column=1, row=1, rowspan=17, columnspan=4,sticky='NSEW')
        self.outputText(int(self.diff_component.get().split(' ')[-1]))
        self.currentPlot = 'GSEA_result_'+self.diff_component.get()

    def outputText(self, diff_component):
        pos_text = str(self.data[self.curKey]['gsea_reports'][diff_component]['pos']).split('\n')
        pos_text = pos_text[1:len(pos_text)-1]
        pos_text = '\n'.join(pos_text)
        neg_text = str(self.data[self.curKey]['gsea_reports'][diff_component]['neg']).split('\n')
        neg_text = neg_text[1:len(neg_text)-1]
        neg_text = '\n'.join(neg_text)
        self.canvas.create_text(5, 5, anchor='nw', text='Positive correlations:\n\n', font=('Helvetica', 16, 'bold'))
        self.canvas.create_text(5, 50, anchor='nw', text=pos_text)
        self.canvas.create_text(5, 150, anchor='nw', text='Negative correlations:\n\n', font=('Helvetica', 16, 'bold'))
        self.canvas.create_text(5, 200, anchor='nw', text=neg_text)

    def plotWBOnTsne(self):
        keys = [key for key in self.data if self.data[key]['state'].get() == True]

        if len(keys) > 2:
            self.WBError = tk.Toplevel()
            self.WBError.title("Error -- too many datasets selected")
            tk.Label(self.WBError,text=u"Please select up to two datasets to plot their wishbone trajectories",fg="black",bg="white").grid(column=0, row=0)
            tk.Button(self.WBError, text="Ok", command=self.WBError.destroy).grid(column=0, row=1)
            self.wait_window(self.WBError)

        else:
            self.saveButton.config(state='normal')
            if self.scseq == True:
                self.component_menu.config(state='disabled')
                self.updateButton.config(state='disabled')
                self.visMenu.entryconfig(6, state='disabled')
            else:
                self.visMenu.entryconfig(4, state='disabled')

            self.resetCanvas()

            self.fig, self.ax = self.data[keys[0]]['wb'].plot_wishbone_on_tsne(other_data=self.data[keys[1]]['wb'] if len(keys) > 1 else None)
            for i in range(len(keys)):
                    plt.figtext(0.05, 0.75 - (0.5*i), keys[i], rotation='vertical')
            self.canvas = FigureCanvasTkAgg(self.fig, self)
            self.canvas.show()
            self.canvas.get_tk_widget().grid(column=1, row=1, rowspan=10, columnspan=4)
            self.currentPlot = 'wishbone_on_tsne'

    def plotWBMarkerTrajectory(self):
        self.getGeneSelection()
        if len(self.selectedGenes) < 1:
            print('Error: must select at least one gene')

        else:
            self.saveButton.config(state='normal')
            if self.scseq == True:
                self.component_menu.config(state='disabled')
                self.updateButton.config(state='disabled')
                self.visMenu.entryconfig(6, state='disabled')
            else:
                self.visMenu.entryconfig(4, state='disabled')

            self.resetCanvas()
            keys = [key for key in self.data if self.data[key]['state'].get() == True]

            self.fig = plt.figure(figsize=[14, 4 * len(keys)])
            gs = gridspec.GridSpec(len(keys), 1)

            for i in range(len(keys)):
                self.ax = self.fig.add_subplot(gs[i, 0])
                self.data[keys[i]]['wb'].plot_marker_trajectory(self.selectedGenes, fig=self.fig, ax=self.ax)
                self.ax.set_title(keys[i])

            self.fig.set_size_inches(10, 6, forward=True)
            self.fig.tight_layout()
            self.fig.subplots_adjust(right=0.8)
                
            self.canvas = FigureCanvasTkAgg(self.fig, self)
            self.canvas.show()
            self.canvas.get_tk_widget().grid(column=1, row=1, rowspan=10, columnspan=5, sticky='W')
            self.currentPlot = 'wishbone_marker_trajectory'
            self.geometry('1050x550')

            #enable buttons
            self.wishboneMenu.entryconfig(2, state='normal')

    def plotWBHeatMap(self):
        keys = [key for key in self.data if self.data[key]['state'].get() == True]
        if len(keys) > 2:
            self.WBError = tk.Toplevel()
            self.WBError.title("Error -- too many datasets selected")
            tk.Label(self.WBError,text=u"Please select up to two datasets to plot expression heatmaps",fg="black",bg="white").grid(column=0, row=0)
            tk.Button(self.WBError, text="Ok", command=self.WBError.destroy).grid(column=0, row=1)
            self.wait_window(self.WBError)

        else:
            self.getGeneSelection()
            if len(self.selectedGenes) < 1:
                print('Error: must select at least one gene')
            else:
                self.saveButton.config(state='normal')
                if self.scseq == True:
                    self.component_menu.config(state='disabled')
                    self.updateButton.config(state='disabled')
                    self.visMenu.entryconfig(6, state='disabled')
                else:
                    self.visMenu.entryconfig(4, state='disabled')

                self.resetCanvas()

                vals1, tmp1, tmp2 = self.data[keys[0]]['wb'].plot_marker_trajectory(self.selectedGenes)
                if len(keys) == 2:
                    vals2, tmp1, tmp2 = self.data[keys[1]]['wb'].plot_marker_trajectory(self.selectedGenes)
                self.fig, self.ax = self.data[keys[0]]['wb'].plot_marker_heatmap(vals1, other_data=[self.data[keys[1]]['wb'], vals2])
                self.fig.set_size_inches(10, 4, forward=True)
                self.fig.tight_layout()
                for i in range(len(keys)):
                    plt.figtext(0.01, 0.75 - (0.5*i), keys[i], rotation='vertical')

                self.canvas = FigureCanvasTkAgg(self.fig, self)
                self.canvas.show()
                self.canvas.get_tk_widget().grid(column=1, row=1, rowspan=10, columnspan=5)
                self.currentPlot = 'wishbone_marker_heatmap'

    def scatterGeneExp(self):
        self.getGeneSelection()
        if len(self.selectedGenes) < 1:
            print('Error: must select at least one gene')
        elif len(self.selectedGenes) > 3:
            print('Error: too many genes selected. Must select either 2 or 3 genes to scatter')
        else:
            self.saveButton.config(state='normal')
            if self.scseq == True:
                self.component_menu.config(state='disabled')
                self.updateButton.config(state='disabled')
                self.visMenu.entryconfig(6, state='disabled')
            else:
                self.visMenu.entryconfig(4, state='disabled')

            self.resetCanvas()
            keys = [key for key in self.data if self.data[key]['state'].get() == True]
            if len(keys) > 1:
                self.fig = plt.figure(figsize=[8, 4 * int(np.ceil(len(keys)/2))])
                gs = gridspec.GridSpec(int(np.ceil(len(keys)/2)), 2)
                for i in range(len(keys)):
                    if len(self.selectedGenes) == 3:
                        self.ax = self.fig.add_subplot(gs[int(i/2), i%2], projection='3d')
                    else:
                        self.ax = self.fig.add_subplot(gs[int(i/2), i%2])
                    self.data[keys[i]]['scdata'].scatter_gene_expression(self.selectedGenes, fig=self.fig, ax=self.ax)
                    self.ax.set_title(keys[i])
                gs.tight_layout(self.fig)

            else:
                self.fig, self.ax = self.scdata.scatter_gene_expression(self.selectedGenes)

            self.canvas = FigureCanvasTkAgg(self.fig, self)
            self.canvas.show()
            self.canvas.get_tk_widget().grid(column=1, row=1, rowspan=10, columnspan=4)
            self.currentPlot = '_'.join(self.selectedGenes) + '_scatter'

    def plotGeneExpOntSNE(self):
        keys = [key for key in self.data if self.data[key]['state'].get() == True]

        if len(keys) > 2:
            self.GEError = tk.Toplevel()
            self.GEError.title("Error -- too many datasets selected")
            tk.Label(self.GEError,text=u"Please select up to two datasets to plot gene expression on tSNE",fg="black",bg="white").grid(column=0, row=0)
            tk.Button(self.GEError, text="Ok", command=self.GEError.destroy).grid(column=0, row=1)
            self.wait_window(self.GEError)

        else:  
            self.getGeneSelection()
            if len(self.selectedGenes) < 1:
                print('Error: must select at least one gene')

            else:
                self.saveButton.config(state='normal')
                if self.scseq == True:
                    self.component_menu.config(state='disabled')
                    self.updateButton.config(state='disabled')
                    self.visMenu.entryconfig(6, state='disabled')
                else:
                    self.visMenu.entryconfig(4, state='disabled')

                self.resetCanvas()
                self.fig, self.ax = self.data[keys[0]]['scdata'].plot_gene_expression(self.selectedGenes, other_data=self.data[keys[1]]['scdata'] if len(keys) > 1 else None)
                for i in range(len(keys)):
                    plt.figtext(0.01, 0.75 - (0.5*i), keys[i], rotation='vertical')
                self.canvas = FigureCanvasTkAgg(self.fig, self)
                self.canvas.show()
                self.canvas.get_tk_widget().grid(column=1, row=1, rowspan=10, columnspan=4)
                self.currentPlot = '_'.join(self.selectedGenes) + '_tsne'
                self.geometry('950x550')

    def scatterGeneExpAgainstOther(self):
        keys = [key for key in self.data if self.data[key]['state'].get() == True]

        if len(keys) != 2:
            self.GEError = tk.Toplevel()
            self.GEError.title("Error -- must select two datasets")
            tk.Label(self.GEError,text=u"Please select exactly two datasets to compare gene expression",fg="black",bg="white").grid(column=0, row=0)
            tk.Button(self.GEError, text="Ok", command=self.GEError.destroy).grid(column=0, row=1)
            self.wait_window(self.GEError)    

        else:  
            self.getGeneSelection()
            if len(self.selectedGenes) < 1:
                print('Error: must select at least one gene')

            else:
                self.saveButton.config(state='normal')
                if self.scseq == True:
                    self.component_menu.config(state='disabled')
                    self.updateButton.config(state='disabled')
                    self.visMenu.entryconfig(6, state='disabled')
                else:
                    self.visMenu.entryconfig(4, state='disabled')

                self.resetCanvas()
                self.fig, self.axes = self.data[keys[0]]['scdata'].scatter_gene_expression_against_other_data(self.selectedGenes, 
                                                                                       other_data=self.data[keys[1]]['scdata'])
                for ax in self.axes:
                    ax.set_xlabel(keys[0])
                    ax.set_ylabel(keys[1])

                self.canvas = FigureCanvasTkAgg(self.fig, self)
                self.canvas.show()
                self.canvas.get_tk_widget().grid(column=1, row=1, rowspan=10, columnspan=4)
                self.currentPlot = '_'.join(self.selectedGenes) + '_compare_gene_exp'

    def getGeneSelection(self):
        #popup menu to get selected genes
        self.geneSelection = tk.Toplevel()
        self.geneSelection.title("Select Genes")
        tk.Label(self.geneSelection,text=u"Genes:",fg="black",bg="white").grid(row=0)

        keys = [key for key in self.data if self.data[key]['state'].get() == True]
        if len(keys) > 1:
            genes = reduce(np.intersect1d, tuple([self.data[key]['genes'] for key in keys])).tolist()
        else:
            genes = self.data[keys[0]]['genes'].tolist()

        self.geneInput = wishbone.autocomplete_entry.AutocompleteEntry(genes, self.geneSelection, listboxLength=6)
        self.geneInput.grid(row=1)
        self.geneInput.bind('<Return>', self.AddToSelected)

        self.geneSelectBox = tk.Listbox(self.geneSelection, selectmode=tk.EXTENDED)
        self.geneSelectBox.grid(row=2, rowspan=10)
        self.geneSelectBox.bind('<BackSpace>', self.DeleteSelected)
        self.selectedGenes = []

        tk.Button(self.geneSelection, text="Use selected genes", command=self.geneSelection.destroy).grid(row=12)
        tk.Button(self.geneSelection, text="Cancel", command=self.cancelGeneSelection).grid(row=13)
        self.wait_window(self.geneSelection)
    
    def cancelGeneSelection(self):
        self.selectedGenes = []
        self.geneSelection.destroy()

    def AddToSelected(self, event):
        self.selectedGenes.append(self.geneInput.get())
        self.geneSelectBox.insert(tk.END, self.selectedGenes[len(self.selectedGenes)-1])

    def DeleteSelected(self, event):
        selected = self.geneSelectBox.curselection()
        pos = 0
        for i in selected:
            idx = int(i) - pos
            self.geneSelectBox.delete( idx,idx )
            self.selectedGenes = self.selectedGenes[:idx] + self.selectedGenes[idx+1:]
            pos = pos + 1  

    def savePlot(self):
        self.plotFileName = filedialog.asksaveasfilename(title='Save Plot', defaultextension='.png', initialfile=self.fileNameEntryVar.get()+"_"+self.currentPlot)
        if self.plotFileName != None:
            self.fig.savefig(self.plotFileName)

    def setGate(self):
        keys = [key for key in self.data if self.data[key]['state'].get() == True]
        if len(keys) > 1:
            raise RuntimeError("Can only select a gate when one data set is plotted.")
        self.curKey = keys[0]
        #pop up for gate name
        self.gateOptions = tk.Toplevel()
        self.gateOptions.title(key + ": Create gate for start cells")
        tk.Label(self.gateOptions,text=u"Gate name:" ,fg="black",bg="white").grid(column=0, row=0)
        self.gateName = tk.StringVar()
        self.gateName.set('Gate ' + str(len(self.data[self.curKey]['gates']) + 1))
        tk.Entry(self.gateOptions, textvariable=self.gateName).grid(column=1,row=0)
        tk.Button(self.gateOptions, text="Select gate", command=self._setGate).grid(column=1, row=1)
        tk.Button(self.gateOptions, text="Cancel", command=self.gateOptions.destroy).grid(column=0, row=1)
        self.wait_window(self.gateOptions)

    def _setGate(self):
        self.gateOptions.destroy()
        self.buttonPress = self.canvas.mpl_connect('button_press_event', self._startGate)
        self.buttonRelease = self.canvas.mpl_connect('button_release_event', self._endGate)
        self.canvas.get_tk_widget().config(cursor='plus')
 
    def _startGate(self, event):
        self.start_x = event.xdata
        self.start_y = event.ydata

    def _endGate(self, event):
        #draw gate rectangle
        start_x = self.start_x if self.start_x < event.xdata else event.xdata
        start_y = self.start_y if self.start_y < event.ydata else event.ydata
        width = np.absolute(event.xdata-self.start_x)
        height = np.absolute(event.ydata-self.start_y)
        rect = Rectangle((start_x, start_y), width, height, 
                         fill=False, ec='black', alpha=1, lw=2)
        self.ax.add_patch(rect)
        self.canvas.draw()

        #disable mouse events
        self.canvas.mpl_disconnect(self.buttonPress)
        self.canvas.mpl_disconnect(self.buttonRelease)
        self.canvas.get_tk_widget().config(cursor='arrow')
 
        #save cell gate
        gate = Path([[start_x, start_y], 
                     [start_x + width, start_y],
                     [start_x + width, start_y + height], 
                     [start_x, start_y + height],
                     [start_x, start_y]])
        gated_cells = self.data[self.curKey]['scdata'].tsne.index[gate.contains_points(self.data[self.curKey]['scdata'].tsne)]
        self.data[self.curKey]['gates'][self.gateName.get()] = gated_cells

        #replot tSNE w gate colored
        self.fig.clf()
        plt.scatter(self.data[self.curKey]['scdata'].tsne['x'], 
                    self.data[self.curKey]['scdata'].tsne['y'], 
                    s=10, edgecolors='none', color='lightgrey')
        plt.scatter(self.data[self.curKey]['scdata'].tsne.ix[gated_cells, 'x'], 
                    self.data[self.curKey]['scdata'].tsne.ix[gated_cells, 'y'], 
                    s=10, edgecolors='none')
        self.canvas.draw()

        if self.scseq == True:
            self.visMenu.entryconfig(6, state='disabled')
        else:
            self.visMenu.entryconfig(4, state='disabled')

    def resetCanvas(self):
        self.fig.clf()
        if type(self.canvas) is FigureCanvasTkAgg:
            for item in self.canvas.get_tk_widget().find_all():
                self.canvas.get_tk_widget().delete(item)
        else:
            for item in self.canvas.find_all():
                self.canvas.delete(item)

    def quitWB(self):
        self.quit()
        self.destroy()

def launch():
    app = wishbone_gui(None)

    if platform.system() == 'Darwin':
        os.system('''/usr/bin/osascript -e 'tell app "Finder" to set frontmost of process "python" to true' ''')
    elif platform.system() == 'Windows':
        app.lift()
        app.call('wm', 'attributes', '.', '-topmost', True)
        app.after_idle(app.call, 'wm', 'attributes', '.', '-topmost', False)        
    elif platform.system() == 'Linux':
        app.focus_force()

    app.title('Wishbone')
    # app.mainloop()
    try:
        app.mainloop()
    except UnicodeDecodeError:
        pass

if __name__ == "__main__":
    launch()

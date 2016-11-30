#!/usr/local/bin/python3

import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle
from matplotlib.path import Path
from functools import reduce, partial
import magic
import os
import sys
import platform
import pandas as pd
import tkinter as tk
import numpy as np
from tkinter import filedialog, ttk
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

        #set up menu bar
        self.menubar = tk.Menu(self)
        self.fileMenu = tk.Menu(self.menubar, tearoff=0)
        self.menubar.add_cascade(label="File", menu=self.fileMenu)
        self.fileMenu.add_command(label="Load csv file", command=self.loadCSV)
        self.fileMenu.add_command(label="Load mtx file", command=self.loadMTX)
        self.fileMenu.add_command(label="Load saved session from pickle file", command=self.loadPickle)
        self.fileMenu.add_command(label="Save data", state='disabled', command=self.saveData)
        self.fileMenu.add_command(label="Exit", command=self.quitWB)

        self.analysisMenu = tk.Menu(self.menubar, tearoff=0)
        self.menubar.add_cascade(label="Analysis", menu=self.analysisMenu)
        self.analysisMenu.add_command(label="Principal component analysis", state='disabled', command=self.runPCA)
        self.analysisMenu.add_command(label="tSNE", state='disabled', command=self.runTSNE)
        self.analysisMenu.add_command(label="Diffusion map", state='disabled', command=self.runDM)
        self.analysisMenu.add_command(label="GSEA", state='disabled', command=self.runGSEA)
        self.analysisMenu.add_command(label="Wishbone", state='disabled', command=self.runWishbone)
        self.analysisMenu.add_command(label="MAGIC", state='disabled', command=self.runMagic)

        self.visMenu = tk.Menu(self.menubar, tearoff=0)
        self.menubar.add_cascade(label="Visualization", menu=self.visMenu)
        self.visMenu.add_command(label="Scatter plot", state='disabled', command=self.scatterPlot)
        
        self.config(menu=self.menubar)

        #intro screen
        tk.Label(self, text=u"MAGIC", font=('Helvetica', 48), fg="black", bg="white", padx=100, pady=50).grid(row=0)
        tk.Label(self, text=u"To get started, select a data file by clicking File > Load Data", fg="black", bg="white", padx=100, pady=25).grid(row=1)

        #update
        self.protocol('WM_DELETE_WINDOW', self.quitWB)
        self.grid_columnconfigure(0,weight=1)
        self.resizable(True,True)
        self.update()
        self.geometry(self.geometry())       
        self.focus_force()

    def loadCSV(self):
        self.dataFileName = filedialog.askopenfilename(title='Load data file', initialdir='~/.magic/data')
        if(self.dataFileName != ""):
            #pop up data options menu
            self.fileInfo = tk.Toplevel()
            self.fileInfo.title("Data options")
            tk.Label(self.fileInfo, text=u"File name: ").grid(column=0, row=0)
            tk.Label(self.fileInfo, text=self.dataFileName.split('/')[-1]).grid(column=1, row=0)

            tk.Label(self.fileInfo,text=u"Name:" ,fg="black",bg="white").grid(column=0, row=1)
            self.fileNameEntryVar = tk.StringVar()
            tk.Entry(self.fileInfo, textvariable=self.fileNameEntryVar).grid(column=1,row=1)

            tk.Button(self.fileInfo, text="Show distributions", command=partial(self.showRawDataDistributions, file_type='csv')).grid(column=0, row=3)

            #filter parameters
            self.filterCellMinVar = tk.StringVar()
            tk.Label(self.fileInfo,text=u"Filter by molecules per cell. Min:" ,fg="black",bg="white").grid(column=0, row=4)
            tk.Entry(self.fileInfo, textvariable=self.filterCellMinVar).grid(column=1,row=4)
            
            self.filterCellMaxVar = tk.StringVar()
            tk.Label(self.fileInfo, text=u" Max:" ,fg="black",bg="white").grid(column=2, row=4)
            tk.Entry(self.fileInfo, textvariable=self.filterCellMaxVar).grid(column=3,row=4)
            
            self.filterGeneNonzeroVar = tk.StringVar()
            tk.Label(self.fileInfo,text=u"Filter by nonzero cells per gene. Min:" ,fg="black",bg="white").grid(column=0, row=5)
            tk.Entry(self.fileInfo, textvariable=self.filterGeneNonzeroVar).grid(column=1,row=5)
            
            self.filterGeneMolsVar = tk.StringVar()
            tk.Label(self.fileInfo,text=u"Filter by molecules per gene. Min:" ,fg="black",bg="white").grid(column=0, row=6)
            tk.Entry(self.fileInfo, textvariable=self.filterGeneMolsVar).grid(column=1,row=6)

            #normalize
            self.normalizeVar = tk.BooleanVar()
            tk.Checkbutton(self.fileInfo, text=u"Normalize by library size", variable=self.normalizeVar).grid(column=0, row=7, columnspan=4)

            tk.Button(self.fileInfo, text="Cancel", command=self.fileInfo.destroy).grid(column=1, row=8)
            tk.Button(self.fileInfo, text="Load", command=partial(self.processData, file_type='csv')).grid(column=2, row=8)

            self.wait_window(self.fileInfo)

    def loadMTX(self):
        self.dataFileName = filedialog.askopenfilename(title='Load data file', initialdir='~/.magic/data')
        if(self.dataFileName != ""):
            #pop up data options menu
            self.fileInfo = tk.Toplevel()
            self.fileInfo.title("Data options")
            tk.Label(self.fileInfo, text=u"File name: ").grid(column=0, row=0)
            tk.Label(self.fileInfo, text=self.dataFileName.split('/')[-1]).grid(column=1, row=0)

            tk.Label(self.fileInfo,text=u"Name:" ,fg="black",bg="white").grid(column=0, row=1)
            self.fileNameEntryVar = tk.StringVar()
            tk.Entry(self.fileInfo, textvariable=self.fileNameEntryVar).grid(column=1,row=1)

            tk.Button(self.fileInfo, text="Select gene names file", command=self.getGeneNameFile).grid(column=0, row=2)

            tk.Button(self.fileInfo, text="Show distributions", command=partial(self.showRawDataDistributions, file_type='mtx')).grid(column=0, row=3)

            #filter parameters
            self.filterCellMinVar = tk.StringVar()
            tk.Label(self.fileInfo,text=u"Filter by molecules per cell. Min:" ,fg="black",bg="white").grid(column=0, row=4)
            tk.Entry(self.fileInfo, textvariable=self.filterCellMinVar).grid(column=1,row=4)
            
            self.filterCellMaxVar = tk.StringVar()
            tk.Label(self.fileInfo, text=u" Max:" ,fg="black",bg="white").grid(column=2, row=4)
            tk.Entry(self.fileInfo, textvariable=self.filterCellMaxVar).grid(column=3,row=4)
            
            self.filterGeneNonzeroVar = tk.StringVar()
            tk.Label(self.fileInfo,text=u"Filter by nonzero cells per gene. Min:" ,fg="black",bg="white").grid(column=0, row=5)
            tk.Entry(self.fileInfo, textvariable=self.filterGeneNonzeroVar).grid(column=1,row=5)
            
            self.filterGeneMolsVar = tk.StringVar()
            tk.Label(self.fileInfo,text=u"Filter by molecules per gene. Min:" ,fg="black",bg="white").grid(column=0, row=6)
            tk.Entry(self.fileInfo, textvariable=self.filterGeneMolsVar).grid(column=1,row=6)

            #normalize
            self.normalizeVar = tk.BooleanVar()
            tk.Checkbutton(self.fileInfo, text=u"Normalize by library size", variable=self.normalizeVar).grid(column=0, row=7, columnspan=4)

            tk.Button(self.fileInfo, text="Cancel", command=self.fileInfo.destroy).grid(column=1, row=8)
            tk.Button(self.fileInfo, text="Load", command=partial(self.processData, file_type='mtx')).grid(column=2, row=8)

            self.wait_window(self.fileInfo)

    def getGeneNameFile(self):
        self.geneNameFile = filedialog.askopenfilename(title='Select gene name file', initialdir='~/.magic/data')
        tk.Label(self.fileInfo,text=self.geneNameFile.split('/')[-1] ,fg="black",bg="white").grid(column=1, row=2)

    def loadPickle(self):
        if(self.dataFileName != ""):
            #pop up data options menu
            self.fileInfo = tk.Toplevel()
            self.fileInfo.title("Data options")
            tk.Label(self.fileInfo, text=u"File name: ").grid(column=0, row=0)
            tk.Label(self.fileInfo, text=self.dataFileName.split('/')[-1]).grid(column=1, row=0)

            tk.Label(self.fileInfo,text=u"Name:" ,fg="black",bg="white").grid(column=0, row=1)
            self.fileNameEntryVar = tk.StringVar()
            tk.Entry(self.fileInfo, textvariable=self.fileNameEntryVar).grid(column=1,row=1)

            tk.Button(self.fileInfo, text="Cancel", command=self.fileInfo.destroy).grid(column=0, row=2)
            tk.Button(self.fileInfo, text="Load", command=partial(self.processData, file_type='pickle')).grid(column=0, row=2)

            self.wait_window(self.fileInfo)

    def processData(self, file_type='csv'):

        if len(self.data) == 0:
            #clear intro screen
            for item in self.grid_slaves():
                item.grid_forget()

            self.data_list = ttk.Treeview()
            self.data_list.heading('#0', text='Data sets')
            self.data_list.grid(column=0, row=0, rowspan=11, sticky='NS')
            self.data_list.bind('<BackSpace>', self._deleteDataItem)

            #set up canvas for plots
            self.fig, self.ax = magic.mg.get_fig()
            self.canvas = FigureCanvasTkAgg(self.fig, self)
            self.canvas.show()
            self.canvas.get_tk_widget().grid(column=1, row=1, rowspan=10, columnspan=4,sticky='NSEW')


        if file_type == 'csv':    # sc-seq data
            scdata = magic.mg.SCData.from_csv(os.path.expanduser(self.dataFileName), 
                                                  data_type='sc-seq', normalize=False)
        elif file_type == 'mtx':   # sparse matrix
            scdata = magic.mg.SCData.from_mtx(os.path.expanduser(self.dataFileName), os.path.expanduser(self.geneNameFile))
            
        if len(self.filterCellMinVar.get()) > 0 or len(self.filterCellMaxVar.get()) > 0 or len(self.filterGeneNonzeroVar.get()) > 0 or len(self.filterGeneMolsVar.get()) > 0:
            scdata.filter_scseq_data(filter_cell_min=int(self.filterCellMinVar.get()) if len(self.filterCellMinVar.get()) > 0 else 0, 
                                     filter_cell_max=int(self.filterCellMaxVar.get()) if len(self.filterCellMaxVar.get()) > 0 else 0, 
                                     filter_gene_nonzero=int(self.filterGeneNonzeroVar.get()) if len(self.filterGeneNonzeroVar.get()) > 0 else 0, 
                                     filter_gene_mols=int(self.filterGeneMolsVar.get()) if len(self.filterGeneMolsVar.get()) > 0 else 0)

        if self.normalizeVar.get() == True:
            scdata = scdata.normalize_scseq_data() 
        wb = None

        if file_type == 'pickle':   # pickled Wishbone object
            wb = magic.mg.Wishbone.load(self.dataFileName)
            scdata = wb.scdata

        self.data[self.fileNameEntryVar.get()] = {'scdata' : scdata, 'wb' : wb, 'state' : tk.BooleanVar(),
                                                  'genes' : scdata.data.columns.values, 'gates' : {}}
        
        self.data_list.insert('', 'end', text=self.fileNameEntryVar.get() + ' (' + str(scdata.data.shape[0]) + ' x ' + str(scdata.data.shape[1]) + ')', open=True)
        
        #set up buttons
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
        self.analysisMenu.entryconfig(1, state='normal')
        if isinstance(scdata.tsne, pd.DataFrame):
            self.analysisMenu.entryconfig(2, state='normal')
        if isinstance(scdata.diffusion_eigenvectors, pd.DataFrame):
            self.analysisMenu.entryconfig(3, state='normal')
            self.analysisMenu.entryconfig(4, state='normal')

        #enable buttons
        self.analysisMenu.entryconfig(0, state='normal')
        self.fileMenu.entryconfig(1, state='normal')
        self.visMenu.entryconfig(0, state='normal')
        self.concatButton = tk.Button(self, text=u"Concatenate selected datasets", state='disabled', wraplength=80, command=self.concatenateData)
        self.concatButton.grid(column=0, row=10)
        if len(self.data) > 1:
            self.concatButton.config(state='normal')

        if wb:
            if isinstance(wb.trajectory, pd.Series):
                self.wishboneMenu.entryconfig(0, state='normal')
                self.wishboneMenu.entryconfig(1, state='normal')
                self.wishboneMenu.entryconfig(2, state='normal')

        self.geometry('1000x700')
        #destroy pop up menu
        self.fileInfo.destroy()

    def saveData(self):
        for key in self.data_list.selection():
            name = self.data_list.item(key)['text'].split(' (')[0]
            pickleFileName = filedialog.asksaveasfilename(title=name + ': save data', defaultextension='.p', initialfile=key)
            if pickleFileName != None:
                if self.data[name]['wb'] != None:
                    self.data[name]['wb'].save(pickleFileName)
                else:
                    self.data[name]['scdata'].save_as_wishbone(pickleFileName)

    def concatenateData(self):
        self.concatOptions = tk.Toplevel()
        self.concatOptions.title("Concatenate data sets")

        tk.Label(self.concatOptions, text=u"New data set name:", fg="black", bg="white").grid(column=0, row=0)
        self.nameVar = tk.StringVar()
        tk.Entry(self.concatOptions, textvariable=self.nameVar).grid(column=1, row=0)

        self.colVar = tk.IntVar()
        tk.Radiobutton(self.concatOptions, text='Concatenate columns', variable=self.colVar, value=0).grid(column=0, row=1)
        tk.Radiobutton(self.concatOptions, text='Concatenate rows', variable=self.colVar, value=1).grid(column=1, row=1)

        self.joinVar = tk.BooleanVar()
        self.joinVar.set(True)
        tk.Checkbutton(self.concatOptions, text=u"Outer join", variable=self.joinVar).grid(column=0, row=2, columnspan=2)

        tk.Button(self.concatOptions, text="Concatenate", command=self._concatenateData).grid(column=1, row=3)
        tk.Button(self.concatOptions, text="Cancel", command=self.concatOptions.destroy).grid(column=0, row=3)
        self.wait_window(self.concatOptions)

    def _concatenateData(self):
        to_concat = []
        for key in self.data_list.selection():
                to_concat.append(self.data[self.data_list.item(key)['text'].split(' (')[0]]['scdata'])

        scdata = to_concat[0].concatenate_data(to_concat[1:], axis=self.colVar.get(), join='outer' if self.joinVar.get() == True else 'inner',
                                               names=[self.data_list.item(key)['text'].split(' (')[0] for key in self.data_list.selection()])


        self.data[self.nameVar.get()] = {'scdata' : scdata, 'wb' : None, 'state' : tk.BooleanVar(),
                                         'genes' : scdata.data.columns.values, 'gates' : {}}
        self.data_list.insert('', 'end', text=self.nameVar.get() + ' (' + str(scdata.data.shape[0]) + ' x ' + str(scdata.data.shape[1]) + ')', open=True)

        self.concatOptions.destroy()

    def _deleteDataItem(self, event):
        for key in self.data_list.selection():
            name = self.data_list.item(key)['text'].split(' (')[0]
            if name in self.data:
                del self.data[name]
            else:
                data_set_name = self.data_list.item(self.data_list.parent(key))['text'].split(' (')[0]
                if 'Principal components' in name:
                    self.data[data_set_name]['scdata'].pca = None
                elif 'tSNE' in name:
                    self.data[data_set_name]['scdata'].tsne = None
                elif 'Diffusion components' in name:
                    self.data[data_set_name]['scdata'].diffusion_eigenvectors = None
                    self.data[data_set_name]['scdata'].diffusion_eigenvalues = None
                elif 'Wishbone' in name:
                    del self.data[data_set_name]['wb']
                elif 'MAGIC' in name:
                    del self.data[data_set_name + ' MAGIC']
                    self.data[data_set_name]['scdata'].magic = None
            self.data_list.delete(key)

    def showRawDataDistributions(self, file_type='csv'):
        if file_type == 'csv':    # sc-seq data
            scdata = magic.mg.SCData.from_csv(os.path.expanduser(self.dataFileName), 
                                              data_type='sc-seq', normalize=False)
        else:   # sparse matrix
            scdata = magic.mg.SCData.from_mtx(os.path.expanduser(self.dataFileName), 
                                              os.path.expanduser(self.geneNameFile))

        self.dataDistributions = tk.Toplevel()
        self.dataDistributions.title(self.fileNameEntryVar.get() + ": raw data distributions")

        fig, ax = scdata.plot_molecules_per_cell_and_gene()
        fig.tight_layout()

        canvas = FigureCanvasTkAgg(fig, self.dataDistributions)
        canvas.show()
        canvas.get_tk_widget().grid(column=0, row=0, rowspan=10, columnspan=4,sticky='NSEW')
        del scdata

        self.wait_window(self.dataDistributions)

    def runPCA(self):
        for key in self.data_list.selection():
            curKey = self.data_list.item(key)['text'].split(' (')[0]
            self.data[curKey]['scdata'].run_pca()

            self.data_list.insert(key, 'end', text=curKey + ' PCA' + 
                              ' (' + str(self.data[curKey]['scdata'].pca['loadings'].shape[0]) + 
                                ' x ' + str(self.data[curKey]['scdata'].pca['loadings'].shape[1]) + ')', open=True)

        #enable buttons
        self.analysisMenu.entryconfig(1, state='normal')

    def runTSNE(self):
        for key in self.data_list.selection():
            #pop up for # components
            self.tsneOptions = tk.Toplevel()
            self.curKey = key
            name = self.data_list.item(key)['text'].split(' (')[0]
            self.tsneOptions.title(name + ": tSNE options")
            if self.data[name]['scdata'].data_type == 'sc-seq':
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
        name = self.data_list.item(self.curKey)['text'].split(' (')[0]
        self.data[name]['scdata'].run_tsne(n_components=self.nCompVar.get(), perplexity=self.perplexityVar.get())
        self.data[name]['gates'] = {}

        #enable buttons
        self.analysisMenu.entryconfig(2, state='normal')

        self.data_list.insert(self.curKey, 'end', text=name + ' tSNE' + 
                              ' (' + str(self.data[name]['scdata'].tsne.shape[0]) + 
                                ' x ' + str(self.data[name]['scdata'].tsne.shape[1]) + ')', open=True)
        self.tsneOptions.destroy()

    def runDM(self):
        for key in self.data_list.selection():
            name = self.data_list.item(key)['text'].split(' (')[0]
            self.data[name]['scdata'].run_diffusion_map()
            self.data_list.insert(key, 'end', text=name + ' Diffusion components' +
                                  ' (' + str(self.data[name]['scdata'].diffusion_eigenvectors.shape[0]) + 
                                 ' x ' + str(self.data[name]['scdata'].diffusion_eigenvectors.shape[1]) + ')', open=True)

        #enable buttons
        self.analysisMenu.entryconfig(3, state='normal')
        self.analysisMenu.entryconfig(4, state='normal')

    def runGSEA(self):

        GSEAFileName = filedialog.askopenfilename(title='Select gmt File', initialdir='~/.magic/tools')

        if GSEAFileName != "":
            if 'mouse' in GSEAFileName:
                gmt_file_type = 'mouse'
            else:
                gmt_file_type = 'human'

            for key in self.data_list.selection():
                name = self.data_list.item(key)['text'].split(' (')[0]
                self.data[name]['scdata'].run_diffusion_map_correlations()
                self.data[name]['scdata'].data.columns = self.data[name]['scdata'].data.columns.str.upper()
                outputPrefix = filedialog.asksaveasfilename(title=name + ': input file prefix for saving output', initialdir='~/.magic/gsea')
                
                self.data[name]['gsea_reports'] = self.data[name]['scdata'].run_gsea(output_stem= os.path.expanduser(outputPrefix), 
                                                                                   gmt_file=(gmt_file_type, GSEAFileName.split('/')[-1]))


    def runWishbone(self):
        for key in self.data_list.selection():
            self.curKey = key
            self.curName = self.data_list.item(key)['text'].split(' (')[0]

            #popup menu for wishbone options
            self.wbOptions = tk.Toplevel()
            self.wbOptions.title(self.curName + ": Wishbone options")

            #s
            tk.Label(self.wbOptions,text=u"Start cell:",fg="black",bg="white").grid(column=0,row=0)
            self.start = tk.StringVar()
            tk.Entry(self.wbOptions, textvariable=self.start).grid(column=1,row=0)
            if(len(self.data[self.curName]['gates']) > 0):
                self.cell_gate = tk.StringVar()
                self.cell_gate.set('Use cell gate')
                self.gate_menu = tk.OptionMenu(self.wbOptions, self.cell_gate,
                                               *list(self.data[self.curName]['gates'].keys()))
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
        self.data[self.curName]['wb'] = magic.mg.Wishbone(self.data[self.curName]['scdata'])

        if self.cell_gate.get() == 'Use cell gate':
            self.data[self.curName]['wb'].run_wishbone(start_cell=self.start.get(), k=self.k.get(), 
                                                      components_list=[int(comp) for comp in self.compList.get().split(',')], 
                                                      num_waypoints=self.numWaypoints.get())
        else:
            #randomly select start cell in gate
            print('Using cell gate:')
            print(self.cell_gate.get())
            start_cell = random.sample(list(self.data[self.curName]['gates'][self.cell_gate.get()]), 1)[0]
            print('Start cell: ' + start_cell)
            self.data[self.curName]['wb'].run_wishbone(start_cell=start_cell, k=self.k.get(), 
                                                      components_list=[int(comp) for comp in self.compList.get().split(',')], 
                                                      num_waypoints=self.numWaypoints.get())
        
        #enable buttons
        self.wishboneMenu.entryconfig(0, state='normal')
        self.wishboneMenu.entryconfig(1, state='normal')
        self.wishboneMenu.entryconfig(2, state='normal')

        self.data_list.insert(self.curKey, 'end', text=self.curName +' Wishbone' +
                              ' (' + str(self.data[name]['wb'].trajectory.shape[0]) + 
                              ' x ' + str(self.data[name]['wb'].trajectory.shape[1]) + ')', open=True)
        self.wbOptions.destroy()

    def runMagic(self):
        for key in self.data_list.selection():
            #pop up for parameters
            self.magicOptions = tk.Toplevel()
            self.magicOptions.title(self.data_list.item(key)['text'].split(' (')[0] + ": MAGIC options")
            self.curKey = key

            tk.Label(self.magicOptions, text=u"Kernel type:", fg="black",bg="white").grid(column=0, row=0)
            self.kernelType = tk.StringVar()
            self.kernelType.set('gaussian')
            kernel_types = ['gaussian', 'tsne']
            self.kernel_menu = tk.OptionMenu(self.magicOptions, self.kernelType, *kernel_types)
            self.kernel_menu.grid(row=0, column=1)

            tk.Label(self.magicOptions,text=u"# of PCA components:" ,fg="black",bg="white").grid(column=0, row=1)
            self.nCompVar = tk.IntVar()
            self.nCompVar.set(20)
            tk.Entry(self.magicOptions, textvariable=self.nCompVar).grid(column=1,row=1)

            tk.Label(self.magicOptions,text=u"t:" ,fg="black",bg="white").grid(column=0, row=2)
            self.tVar = tk.IntVar()
            self.tVar.set(8)
            tk.Entry(self.magicOptions, textvariable=self.tVar).grid(column=1,row=2)

            tk.Label(self.magicOptions, text=u"Gaussian kernel only:", fg="black", bg="white").grid(column=0, row=3, columnspan=2)
            
            tk.Label(self.magicOptions,text=u"kNN:" ,fg="black",bg="white").grid(column=0, row=4)
            self.kNNVar = tk.IntVar()
            self.kNNVar.set(20)
            tk.Entry(self.magicOptions, textvariable=self.kNNVar).grid(column=1,row=4)

            tk.Label(self.magicOptions,text=u"kNN-autotune:" ,fg="black",bg="white").grid(column=0, row=5)
            self.autotuneVar = tk.IntVar()
            self.autotuneVar.set(0)
            tk.Entry(self.magicOptions, textvariable=self.autotuneVar).grid(column=1,row=5)

            tk.Label(self.magicOptions,text=u"Epsilon:" ,fg="black",bg="white").grid(column=0, row=6)
            self.epsilonVar = tk.IntVar()
            self.epsilonVar.set(0)
            tk.Entry(self.magicOptions, textvariable=self.epsilonVar).grid(column=1,row=6)

            tk.Label(self.magicOptions, text=u"tSNE kernel only:", fg='black', bg='white').grid(column=0, row=7, columnspan=2)

            tk.Label(self.magicOptions,text=u"Perplexity:" ,fg="black",bg="white").grid(column=0, row=8)
            self.perplexityVar = tk.IntVar()
            self.perplexityVar.set(30)
            tk.Entry(self.magicOptions, textvariable=self.perplexityVar).grid(column=1,row=8)

            tk.Label(self.magicOptions,text=u"k_kNN:" ,fg="black",bg="white").grid(column=0, row=9)
            self.k_kNNVar = tk.IntVar()
            self.k_kNNVar.set(100)
            tk.Entry(self.magicOptions, textvariable=self.k_kNNVar).grid(column=1,row=9)

            self.rescaleVar = tk.BooleanVar()
            tk.Checkbutton(self.magicOptions, text=u"Rescale data to 99th percentile", variable=self.rescaleVar).grid(column=0, row=10, columnspan=2)

            tk.Button(self.magicOptions, text="Cancel", command=self.magicOptions.destroy).grid(column=0, row=11)
            tk.Button(self.magicOptions, text="Run", command=self._runMagic).grid(column=1, row=11)
            self.wait_window(self.magicOptions)

    def _runMagic(self):
        name = self.data_list.item(self.curKey)['text'].split(' (')[0]
        self.data[name]['scdata'].run_magic(n_pca_components=self.nCompVar.get() if self.nCompVar.get() > 0 else None,
                                            t=self.tVar.get(), knn=self.kNNVar.get(), epsilon=self.epsilonVar.get(), 
                                            rescale=self.rescaleVar.get(), knn_autotune=self.autotuneVar.get(),
                                            perplexity=self.perplexityVar.get(), k_knn=self.k_kNNVar.get(),
                                            kernel=self.kernelType.get())
        
        self.data[name + ' MAGIC'] = {'scdata' : self.data[name]['scdata'].magic, 'wb' : None, 'state' : tk.BooleanVar(),
                                      'genes' : self.data[name]['scdata'].magic.data.columns.values, 'gates' : {}}
        
        self.data_list.insert(self.curKey, 'end', text=name + ' MAGIC' +
                              ' (' + str(self.data[name]['scdata'].magic.data.shape[0]) + 
                              ' x ' + str(self.data[name]['scdata'].magic.data.shape[1]) + ')', open=True)
        
        self.magicOptions.destroy()

    def plotPCA(self):
        self.saveButton.config(state='normal')
        self.component_menu.config(state='disabled')
        self.updateButton.config(state='disabled')

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
        keys = self.data_list.selection()
        if len(keys) > 1:
            self.fig = plt.figure(figsize=[8, 4 * int(np.ceil(len(keys)/2))])
            gs = gridspec.GridSpec(int(np.ceil(len(keys)/2)), 2)
            for i in range(len(keys)):
                name = self.data_list.item(keys[i])['text'].split(' (')[0]
                if 'PCA' in name:
                    name = name.split(' PCA')[0]
                self.ax = self.fig.add_subplot(gs[int(i/2), i%2])
                self.data[name]['scdata'].plot_pca_variance_explained(fig=self.fig, ax=self.ax, ylim=(0, self.yLimVar.get()), n_components=self.compVar.get())
                self.ax.set_title(name)
            gs.tight_layout(self.fig)

        else:
            name = self.data_list.item(keys[0])['text'].split(' (')[0]
            if 'PCA' in name:
                    name = name.split(' PCA')[0]
            self.fig, self.ax = self.data[name]['scdata'].plot_pca_variance_explained(ylim=(0, self.yLimVar.get()), n_components=self.compVar.get())
        
        self.canvas = FigureCanvasTkAgg(self.fig, self)
        self.canvas.show()
        self.canvas.get_tk_widget().grid(column=1, row=1, rowspan=10, columnspan=4) 
        self.currentPlot = 'pca'

        #enable buttons
        self.saveButton.config(state='normal')
        self.PCAOptions.destroy()

    def plotTSNE(self):
        keys = self.data_list.selection()
        self.getScatterSelection(plot_type='tsne')

        self.saveButton.config(state='normal')
        self.component_menu.config(state='disabled')
        self.updateButton.config(state='disabled')

        self.colorSelection = self.colorVar.get().split(', ')
        if (len(self.colorSelection) == 1 and len(self.colorSelection[0]) > 0) or len(self.colorSelection) > 1:
            self.resetCanvas()
            self.fig = plt.figure(figsize=[4*len(self.colorSelection), 4 * len(keys)])
            gs = gridspec.GridSpec(len(keys), len(self.colorSelection))
            for i in range(len(keys)):
                name = self.data_list.item(keys[i])['text'].split(' (')[0]
                if 'tSNE' in name:
                    name = name.split(' tSNE')[0]
                for j in range(len(self.colorSelection)):
                    self.ax = self.fig.add_subplot(gs[i, j])
                    if self.colorSelection[j] == 'wishbone trajectory':
                        self.fig, self.ax = self.data[name]['scdata'].plot_tsne(fig=self.fig, ax=self.ax,
                                                                                color=self.data[name]['wb'].trajectory)
                    elif self.colorSelection[j] == 'wishbone branches':
                        self.fig, self.ax = self.data[name]['scdata'].plot_tsne(fig=self.fig, ax=self.ax, 
                                                                                color=self.data[name]['wb'].branch)
                    elif 'diffusion component' in self.colorSelection[j]:
                        color = int(self.colorSelection[j].split('diffusion component ')[1])
                        self.fig, self.ax = self.data[name]['scdata'].plot_tsne(fig=self.fig, ax=self.ax, 
                                                                                color=self.data[name]['scdata'].diffusion_eigenvectors[color])
                    elif 'magic' in self.colorSelection[j]:
                        color = self.colorSelection[j].split(' magic')[0]
                        self.fig, self.ax = self.data[name]['scdata'].plot_tsne(fig=self.fig, ax=self.ax, 
                                                                                color=self.data[name]['scdata'].magic.data[color])
                    elif self.colorSelection[j] in self.data[name]['scdata'].data.columns:
                        self.fig, self.ax = self.data[name]['scdata'].plot_tsne(fig=self.fig, ax=self.ax, 
                                                                                color=self.data[name]['scdata'].data[self.colorSelection[j]])
                    elif self.colorSelection[j] == 'density':
                        self.data[name]['scdata'].plot_tsne(fig=self.fig, ax=self.ax, density=True)
                    else:
                        self.data[name]['scdata'].plot_tsne(fig=self.fig, ax=self.ax, color=self.colorSelection[j])
                    self.ax.set_title('Color: ' + self.colorSelection[j])
                    self.ax.set_xlabel(name + ' tsne_x')
                    self.ax.set_ylabel(name + ' tsne_y')
            gs.tight_layout(self.fig)

            self.canvas = FigureCanvasTkAgg(self.fig, self)
            self.canvas.show()
            self.canvas.get_tk_widget().grid(column=1, row=1, rowspan=10, columnspan=4) 
            self.currentPlot = 'tsne'
        
        # elif len(self.colorSelection) > 1:
        #     self.plotGeneExpOntSNE()

    def plotDM(self):

        keys = self.data_list.selection()
        if len(keys) > 2:
            self.DMError = tk.Toplevel()
            self.DMError.title("Error -- too many datasets selected")
            tk.Label(self.DMError,text=u"Please select up to two datasets to plot diffusion map components",fg="black",bg="white").grid(column=0, row=0)
            tk.Button(self.DMError, text="Ok", command=self.DMError.destroy).grid(column=0, row=1)
            self.wait_window(self.DMError)

        else:
            self.saveButton.config(state='normal')
            self.component_menu.config(state='disabled')
            self.updateButton.config(state='disabled')

            self.geometry('950x550')
            self.resetCanvas()

            self.fig, self.ax = self.data[self.data_list.item(keys[0])['text'].split(' (')[0]]['scdata'].plot_diffusion_components(other_data=self.data[self.data_list.item(keys[1])['text'].split(' (')[0]]['scdata'] if len(keys) > 1 else None)
            for i in range(len(keys)):
                plt.figtext(-0.05, 0.75 - (0.5*i), self.data_list.item(keys[i])['text'].split(' (')[0], rotation='vertical')

            self.canvas = FigureCanvasTkAgg(self.fig, self)
            self.canvas.show()
            self.canvas.get_tk_widget().grid(column=1, row=1, rowspan=10, columnspan=4,sticky='W') 
            self.currentPlot = 'dm_components'

    def showGSEAResults(self):
        keys = self.data_list.selection()
        if len(keys) != 1:
            self.GSEAError = tk.Toplevel()
            self.GSEAError.title("Error -- too many datasets selected")
            tk.Label(self.GSEAError,text=u"Please select exactly one dataset to view GSEA results.",fg="black",bg="white").grid(column=0, row=0)
            tk.Button(self.GSEAError, text="Ok", command=self.GSEAError.destroy).grid(column=0, row=1)
            self.wait_window(self.GSEAError)
        else:
            self.curName = self.data_list.item(keys[0])['text'].split(' (')[0]
            self.saveButton.config(state='disabled')
            self.component_menu.config(state='normal')
            self.updateButton.config(state='normal')

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
        pos_text = str(self.data[self.curName]['gsea_reports'][diff_component]['pos']).split('\n')
        pos_text = pos_text[1:len(pos_text)-1]
        pos_text = '\n'.join(pos_text)
        neg_text = str(self.data[self.curName]['gsea_reports'][diff_component]['neg']).split('\n')
        neg_text = neg_text[1:len(neg_text)-1]
        neg_text = '\n'.join(neg_text)
        self.canvas.create_text(5, 5, anchor='nw', text='Positive correlations:\n\n', font=('Helvetica', 16, 'bold'))
        self.canvas.create_text(5, 50, anchor='nw', text=pos_text)
        self.canvas.create_text(5, 150, anchor='nw', text='Negative correlations:\n\n', font=('Helvetica', 16, 'bold'))
        self.canvas.create_text(5, 200, anchor='nw', text=neg_text)

    def plotWBOnTsne(self):
        keys = self.data_list.selection()

        if len(keys) > 2:
            self.WBError = tk.Toplevel()
            self.WBError.title("Error -- too many datasets selected")
            tk.Label(self.WBError,text=u"Please select up to two datasets to plot their wishbone trajectories",fg="black",bg="white").grid(column=0, row=0)
            tk.Button(self.WBError, text="Ok", command=self.WBError.destroy).grid(column=0, row=1)
            self.wait_window(self.WBError)

        else:
            self.saveButton.config(state='normal')
            self.component_menu.config(state='disabled')
            self.updateButton.config(state='disabled')

            self.resetCanvas()

            self.fig, self.ax = self.data[self.data_list.item(keys[0])['text'].split(' (')[0]]['wb'].plot_wishbone_on_tsne(other_data=self.data[self.data_list.item(keys[1])['text'].split(' (')[0]]['wb'] if len(keys) > 1 else None)
            for i in range(len(keys)):
                    plt.figtext(0.05, 0.75 - (0.5*i), self.data_list.item(keys[i])['text'].split(' (')[0], rotation='vertical')
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
            self.component_menu.config(state='disabled')
            self.updateButton.config(state='disabled')

            self.resetCanvas()
            keys = self.data_list.selection()

            self.fig = plt.figure(figsize=[14, 4 * len(keys)])
            gs = gridspec.GridSpec(len(keys), 1)

            for i in range(len(keys)):
                name = self.data_list.item(keys[i])['text'].split(' (')[0]
                self.ax = self.fig.add_subplot(gs[i, 0])
                self.data[name]['wb'].plot_marker_trajectory(self.selectedGenes, fig=self.fig, ax=self.ax)
                self.ax.set_title(name)

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
        keys = self.data_list.selection()
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
                self.component_menu.config(state='disabled')
                self.updateButton.config(state='disabled')

                self.resetCanvas()

                name0 = self.data_list.item(keys[0])['text'].split(' (')[0]
                name1 = self.data_list.item(keys[1])['text'].split(' (')[0]
                vals1, tmp1, tmp2 = self.data[name0]['wb'].plot_marker_trajectory(self.selectedGenes)
                if len(keys) == 2:
                    vals2, tmp1, tmp2 = self.data[name1]['wb'].plot_marker_trajectory(self.selectedGenes)
                self.fig, self.ax = self.data[name0]['wb'].plot_marker_heatmap(vals1, other_data=[self.data[name1]['wb'], vals2])
                self.fig.set_size_inches(10, 4, forward=True)
                self.fig.tight_layout()
                for i in range(len(keys)):
                    name = self.data_list.item(keys[i])['text'].split(' (')[0]
                    plt.figtext(0.01, 0.75 - (0.5*i), name, rotation='vertical')

                self.canvas = FigureCanvasTkAgg(self.fig, self)
                self.canvas.show()
                self.canvas.get_tk_widget().grid(column=1, row=1, rowspan=10, columnspan=5)
                self.currentPlot = 'wishbone_marker_heatmap'

    def scatterPlot(self):
        keys = self.data_list.selection()
        if 'tSNE' in self.data_list.item(keys[0])['text']:
            self.plotTSNE()
        elif 'PCA' in self.data_list.item(keys[0])['text']:
            self.plotPCA()
        else:
            self.getScatterSelection()
            xSelection = self.xVar.get().split(', ')
            ySelection = self.yVar.get().split(', ')
            zSelection = self.zVar.get().split(', ')
            colorSelection = self.colorVar.get().split(', ')
            if len(xSelection[0]) > 0 and len(ySelection[0]) > 0 and len(xSelection) == len(ySelection):
                if len(colorSelection) == 1 and len(colorSelection) != len(xSelection):
                    colorSelection = np.repeat(colorSelection, len(xSelection))

                self.saveButton.config(state='normal')
                self.component_menu.config(state='disabled')
                self.updateButton.config(state='disabled')

                self.resetCanvas()
                keys = self.data_list.selection()
                self.fig = plt.figure(figsize=[4*len(xSelection), 4 * len(keys)])
                gs = gridspec.GridSpec(len(keys), len(xSelection))
                self.ax = []
                for i in range(len(keys)):
                    name = self.data_list.item(keys[i])['text'].split(' (')[0]
                    magic = False
                    if 'MAGIC' in name:
                        name = name.split(' MAGIC')[0]
                        magic = True

                    for j in range(len(xSelection)):
                        if len(zSelection[0]) > 0 and len(zSelection) == len(xSelection):
                            self.ax.append(self.fig.add_subplot(gs[i, j], projection='3d'))
                            genes = [xSelection[j], ySelection[j], zSelection[j]]
                        else:
                            self.ax.append(self.fig.add_subplot(gs[i, j]))
                            genes = [xSelection[j], ySelection[j]]

                        if colorSelection[j] == 'wishbone trajectory':
                            if magic:
                                self.data[name]['scdata'].magic.scatter_gene_expression(genes, fig=self.fig, ax=self.ax[len(self.ax)-1],
                                                                                                            color=self.data[name]['wb'].trajectory)
                            else:
                                self.data[name]['scdata'].scatter_gene_expression(genes, fig=self.fig, ax=self.ax[len(self.ax)-1],
                                                                                                            color=self.data[name]['wb'].trajectory)
                        elif colorSelection[j] == 'wishbone branches':
                            if magic:
                                self.data[name]['scdata'].magic.scatter_gene_expression(genes, fig=self.fig, ax=self.ax[len(self.ax)-1],
                                                                                                            color=self.data[name]['wb'].branch)
                            else:
                                self.data[name]['scdata'].scatter_gene_expression(genes, fig=self.fig, ax=self.ax[len(self.ax)-1],
                                                                                                            color=self.data[name]['wb'].branch)
                        elif 'diffusion component' in colorSelection[j]:
                            color = int(colorSelection[j].split('diffusion component ')[1])
                            if magic:
                                self.data[name]['scdata'].magic.scatter_gene_expression(genes, fig=self.fig, ax=self.ax[len(self.ax)-1],
                                                                                                            color=self.data[name]['scdata'].diffusion_eigenvectors[color])
                            else:
                                self.data[name]['scdata'].scatter_gene_expression(genes, fig=self.fig, ax=self.ax[len(self.ax)-1],
                                                                                                            color=self.data[name]['scdata'].diffusion_eigenvectors[color])
                        elif 'magic' in colorSelection[j]:
                            color = self.colorSelection[j].split(' magic')[0]
                            if magic:
                                self.data[name]['scdata'].magic.scatter_gene_expression(genes, fig=self.fig, ax=self.ax[len(self.ax)-1],
                                                                                                            color=self.data[name]['scdata'].magic.data[color])
                            else:
                                self.data[name]['scdata'].scatter_gene_expression(genes, fig=self.fig, ax=self.ax[len(self.ax)-1],
                                                                                                            color=self.data[name]['scdata'].magic.data[color])
                        elif colorSelection[j] in self.data[name]['scdata'].data.columns:
                            if magic:
                                self.data[name]['scdata'].magic.scatter_gene_expression(genes, fig=self.fig, ax=self.ax[len(self.ax)-1],
                                                                                                            color=self.data[name]['scdata'].magic.data[colorSelection[j]])
                            else:
                                self.data[name]['scdata'].scatter_gene_expression(genes, fig=self.fig, ax=self.ax[len(self.ax)-1],
                                                                                                            color=self.data[name]['scdata'].data[colorSelection[j]])
                        elif colorSelection[j] == 'density':
                            if magic:
                                self.data[name]['scdata'].magic.scatter_gene_expression(genes, fig=self.fig, ax=self.ax[len(self.ax)-1],
                                                                                                            density=True)
                            else:
                                self.data[name]['scdata'].scatter_gene_expression(genes, fig=self.fig, ax=self.ax[len(self.ax)-1],
                                                                                                            density=True)
                        else:
                            if magic:
                                self.data[name]['scdata'].magic.scatter_gene_expression(genes, fig=self.fig, ax=self.ax[len(self.ax)-1],
                                                                                                            color=colorSelection[j])
                            else:
                                self.data[name]['scdata'].scatter_gene_expression(genes, fig=self.fig, ax=self.ax[len(self.ax)-1],
                                                                                                            color=colorSelection[j])

                        self.ax[len(self.ax)-1].set_title('Color: ' + colorSelection[j])
                        self.ax[len(self.ax)-1].set_xlabel(name + ' ' + xSelection[j])
                        self.ax[len(self.ax)-1].set_ylabel(name + ' ' + ySelection[j])
                        if len(genes) == 3:
                            self.ax[len(self.ax)-1].set_zlabel(name + ' ' + zSelection[j])
                
                gs.tight_layout(self.fig, pad=1.2, w_pad=0.1)
                
                self.canvas = FigureCanvasTkAgg(self.fig, self)
                self.canvas.get_tk_widget().grid(column=1, row=1, rowspan=10, columnspan=4)
                self.canvas.show()
                if len(zSelection[0]) > 0 and len(zSelection) == len(xSelection):
                        for ax in self.ax:
                            ax.mouse_init()
                self.currentPlot = 'scatter'
            else:
                print('Error: must select at least one gene for x and y. x and y must also have the same number of selections.')

    def getScatterSelection(self, plot_type=''):
        #popup menu for scatter plot selections
        self.scatterSelection = tk.Toplevel()
        self.scatterSelection.title("Scatter plot options")

        #x
        if plot_type == 'tsne':
            tk.Label(self.scatterSelection, text=u"x: tsne_x", fg="black",bg="white").grid(row=0, column=0)
        else:
            tk.Label(self.scatterSelection, text=u"x:", fg="black",bg="white").grid(row=0, column=0)
            self.xVar = tk.StringVar()
            tk.Entry(self.scatterSelection, textvariable=self.xVar).grid(column=1,row=0)

        #y
        if plot_type == 'tsne':
            tk.Label(self.scatterSelection, text=u"y: tsne_y", fg="black",bg="white").grid(row=1, column=0)
        else:
            tk.Label(self.scatterSelection, text=u"y:", fg="black",bg="white").grid(row=1, column=0)
            self.yVar = tk.StringVar()
            tk.Entry(self.scatterSelection, textvariable=self.yVar).grid(column=1,row=1)

        #z
        if plot_type != 'tsne':
            tk.Label(self.scatterSelection, text=u"z:", fg="black",bg="white").grid(row=2, column=0)
            self.zVar = tk.StringVar()
            tk.Entry(self.scatterSelection, textvariable=self.zVar).grid(column=1,row=2)

        #color
        tk.Label(self.scatterSelection, text=u"color:", fg="black",bg="white").grid(row=3, column=0)
        self.colorVar = tk.StringVar()
        self.colorVar.set('blue')
        tk.Entry(self.scatterSelection, textvariable=self.colorVar).grid(column=1,row=3)

        tk.Button(self.scatterSelection, text="Plot", command=self.scatterSelection.destroy).grid(column=1, row=4)
        tk.Button(self.scatterSelection, text="Cancel", command=self._cancelScatter).grid(column=0, row=4)
        self.wait_window(self.scatterSelection)

    def _cancelScatter(self):
        self.colorVar.set('')
        self.scatterSelection.destroy()

    def savePlot(self):
        self.plotFileName = filedialog.asksaveasfilename(title='Save Plot', defaultextension='.png', initialfile=self.fileNameEntryVar.get()+"_"+self.currentPlot)
        if self.plotFileName != None:
            self.fig.savefig(self.plotFileName)

    def setGate(self):
        keys = self.data_list.selection()
        if len(keys) > 1:
            raise RuntimeError("Can only select a gate when one data set is plotted.")
        self.curKey = self.data_list.item(keys[0])
        #pop up for gate name
        self.gateOptions = tk.Toplevel()
        self.gateOptions.title(self.curKey + ": Create gate for start cells")
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
    print(platform.system())
    if platform.system() == 'Darwin':
        os.system('''/usr/bin/osascript -e 'tell app "Finder" to set frontmost of process "python" to true' ''')
    elif platform.system() == 'Windows':
        app.lift()
        app.call('wm', 'attributes', '.', '-topmost', True)
        app.after_idle(app.call, 'wm', 'attributes', '.', '-topmost', False)        
    elif platform.system() == 'Linux':
        app.focus_force()

    app.title('MAGIC')
    # app.mainloop()
    try:
        app.mainloop()
    except UnicodeDecodeError:
        pass

if __name__ == "__main__":
    launch()

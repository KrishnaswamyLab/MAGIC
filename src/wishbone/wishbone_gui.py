#!/usr/local/bin/python3

import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle
from matplotlib.path import Path
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

        #set up menu bar
        self.menubar = tk.Menu(self)
        self.fileMenu = tk.Menu(self.menubar, tearoff=0)
        self.menubar.add_cascade(label="File", menu=self.fileMenu)
        self.fileMenu.add_command(label="Load data", command=self.loadData)
        self.fileMenu.add_command(label="Exit Wishbone", command=self.quitWB)

        self.analysisMenu = tk.Menu(self.menubar, tearoff=0)
        self.menubar.add_cascade(label="Analysis", menu=self.analysisMenu)
        self.analysisMenu.add_command(label="Principal component analysis", state='disabled', command=self.runPCA)
        self.analysisMenu.add_command(label="tSNE", state='disabled', command=self.runTSNE)
        self.analysisMenu.add_command(label="Diffusion map", state='disabled', command=self.runDM)
        self.analysisMenu.add_command(label="GSEA", state='disabled', command=self.runGSEA)
        self.analysisMenu.add_command(label="Wishbone", state='disabled', command=self.runWishbone)

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
        self.visMenu.add_command(label="Gene expression", state='disabled', command=self.plotGeneExpOntSNE)
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

            self.normalizeVar = tk.BooleanVar()
            tk.Checkbutton(self.fileInfo, text=u"Normalize", variable=self.normalizeVar).grid(column=0, row=2, columnspan=2)
            tk.Label(self.fileInfo, text=u"The normalize parameter is used for correcting for library size among cells.").grid(column=0, row=3, columnspan=2)

            tk.Button(self.fileInfo, text="Cancel", command=self.fileInfo.destroy).grid(column=0, row=4)
            tk.Button(self.fileInfo, text="Load", command=self.processData).grid(column=1, row=4)

            self.wait_window(self.fileInfo)

    def processData(self):
        #clear intro screen
        for item in self.grid_slaves():
            item.grid_forget()

        #load data
        self.scdata = wishbone.wb.SCData.from_csv(os.path.expanduser(self.dataFileName), 
                data_type='sc-seq', normalize=self.normalizeVar.get())
        #get genes
        self.genes = self.scdata.data.columns.values

        #display file name
        tk.Label(self,text=u"File name: " + self.fileNameEntryVar.get(), fg="black",bg="white").grid(column=0,row=0)

        #set up canvas for plots
        self.fig, self.ax = wishbone.wb.get_fig()
        self.canvas = FigureCanvasTkAgg(self.fig, self)
        self.canvas.show()
        self.canvas.get_tk_widget().grid(column=1, row=1, rowspan=10, columnspan=4,sticky='NSEW')

        #set up visualization buttons
        tk.Label(self, text=u"Visualizations:", fg='black', bg='white').grid(column=0, row=1)
        self.PCAButton = tk.Button(self, text=u"PCA", state='disabled', command=self.plotPCA)
        self.PCAButton.grid(column=0, row=2)
        self.tSNEButton = tk.Button(self, text=u"tSNE", state='disabled', command=self.plotTSNE)
        self.tSNEButton.grid(column=0, row=3)
        self.DMButton = tk.Button(self, text=u"Diffusion map", state='disabled', command=self.plotDM)
        self.DMButton.grid(column=0, row=4)
        self.GSEAButton = tk.Button(self, text=u"GSEA Results", state='disabled', command=self.showGSEAResults)
        self.GSEAButton.grid(column=0, row=5)
        self.WBButton = tk.Button(self, text=u"Wishbone", state='disabled', command=self.plotWBOnTsne)
        self.WBButton.grid(column=0, row=6)
        self.geneExpButton = tk.Button(self, text=u"Gene expression", state='disabled', command=self.plotGeneExpOntSNE)
        self.geneExpButton.grid(column=0, row=7)
        self.setGateButton = tk.Button(self, text=u"Set gate", state='disabled', command=self.setGate)
        self.setGateButton.grid(column=0, row=8)
        self.saveButton = tk.Button(self, text=u"Save plot", state='disabled', command=self.savePlot)
        self.saveButton.grid(column = 4, row=0)
        self.diff_component = tk.StringVar()
        self.diff_component.set('Component 1')
        self.component_menu = tk.OptionMenu(self, self.diff_component,
                                            'Component 1', 'Component 2', 'Component 3',
                                            'Component 4', 'Component 5', 'Component 6', 
                                            'Component 7', 'Component 8', 'Component 9')
        self.component_menu.config(state='disabled')
        self.component_menu.grid(row=0, column=2)
        self.updateButton = tk.Button(self, text=u"Update component", command=self.updateComponent, state='disabled')
        self.updateButton.grid(column=3, row=0)

        #enable buttons
        self.analysisMenu.entryconfig(0, state='normal')
        self.geometry('800x550')
        #destroy pop up menu
        self.fileInfo.destroy()

    def runPCA(self):
        self.scdata.run_pca()

        #enable buttons
        self.analysisMenu.entryconfig(1, state='normal')
        self.visMenu.entryconfig(0, state='normal')
        self.PCAButton.config(state='normal')

    def runTSNE(self):
        #pop up for # components
        self.tsneOptions = tk.Toplevel()
        self.tsneOptions.title("tSNE options")
        tk.Label(self.tsneOptions,text=u"Number of components:" ,fg="black",bg="white").grid(column=0, row=0)
        self.nCompVar = tk.IntVar()
        self.nCompVar.set(5)
        tk.Entry(self.tsneOptions, textvariable=self.nCompVar).grid(column=1,row=0)
        tk.Button(self.tsneOptions, text="Run", command=self._runTSNE).grid(column=1, row=1)
        tk.Button(self.tsneOptions, text="Cancel", command=self.tsneOptions.destroy).grid(column=0, row=1)
        self.wait_window(self.tsneOptions)

    def _runTSNE(self):
        self.scdata.run_tsne(n_components=self.nCompVar.get())
        self.gates = {}

        #enable buttons
        self.analysisMenu.entryconfig(2, state='normal')
        self.visMenu.entryconfig(1, state='normal')
        self.visMenu.entryconfig(5, state='normal')
        self.tSNEButton.config(state='normal')
        self.geneExpButton.config(state='normal')
        self.tsneOptions.destroy()

    def runDM(self):
        self.scdata.run_diffusion_map()

        #enable buttons
        self.analysisMenu.entryconfig(3, state='normal')
        self.analysisMenu.entryconfig(4, state='normal')
        self.visMenu.entryconfig(2, state='normal')
        self.DMButton.config(state='normal')

    def runGSEA(self):
        self.GSEAFileName = filedialog.askopenfilename(title='Select gmt File', initialdir='~/.wishbone/tools')
        if self.GSEAFileName != "":
            self.scdata.run_diffusion_map_correlations()
            self.scdata.data.columns = self.scdata.data.columns.str.upper()
            self.outputPrefix = filedialog.asksaveasfilename(title='Input file prefix for saving output', initialdir='~/.wishbone/gsea')
            if 'mouse' in self.GSEAFileName:
                gmt_file_type = 'mouse'
            else:
                gmt_file_type = 'human'
            self.reports = self.scdata.run_gsea(output_stem= os.path.expanduser(self.outputPrefix), 
                     gmt_file=(gmt_file_type, self.GSEAFileName.split('/')[-1]))
            #enable buttons
            self.visMenu.entryconfig(3, state='normal')
            self.GSEAButton.config(state='normal')

    def runWishbone(self):
        #popup menu for wishbone options
        self.wbOptions = tk.Toplevel()
        self.wbOptions.title("Wishbone Options")

        #s
        tk.Label(self.wbOptions,text=u"Start cell:",fg="black",bg="white").grid(column=0,row=0)
        self.start = tk.StringVar()
        tk.Entry(self.wbOptions, textvariable=self.start).grid(column=1,row=0)
        if(len(self.gates) > 0):
            self.cell_gate = tk.StringVar()
            self.cell_gate.set('Use cell gate')
            self.gate_menu = tk.OptionMenu(self.wbOptions, self.cell_gate,
                                           *list(self.gates.keys()))
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
        self.wb = wishbone.wb.Wishbone(self.scdata)

        if self.cell_gate.get() == 'Use cell gate':
            self.wb.run_wishbone(start_cell=self.start.get(), k=self.k.get(), components_list=[int(comp) for comp in self.compList.get().split(',')], num_waypoints=self.numWaypoints.get())
        else:
            #randomly select start cell in gate
            print('Using cell gate:')
            print(self.cell_gate.get())
            start_cell = random.sample(list(self.gates[self.cell_gate.get()]), 1)[0]
            print(start_cell)
            self.wb.run_wishbone(start_cell=start_cell, k=self.k.get(), components_list=[int(comp) for comp in self.compList.get().split(',')], num_waypoints=self.numWaypoints.get())
        
        #enable buttons
        self.wishboneMenu.entryconfig(0, state='normal')
        self.wishboneMenu.entryconfig(1, state='normal')
        self.wishboneMenu.entryconfig(2, state='normal')
        self.WBButton.config(state='normal')
        self.wbOptions.destroy()

    def plotPCA(self):
        self.saveButton.config(state='normal')
        self.component_menu.config(state='disabled')
        self.updateButton.config(state='disabled')
        self.setGateButton.config(state='disabled')
        self.visMenu.entryconfig(6, state='disabled')

        #pop up for # components
        self.PCAOptions = tk.Toplevel()
        self.PCAOptions.title("PCA Plot Options")
        tk.Label(self.PCAOptions,text=u"Max variance explained (ylim):",fg="black",bg="white").grid(column=0, row=0)
        self.yLimVar = tk.DoubleVar()
        self.yLimVar.set(round(self.scdata.pca['eigenvalues'][0][0], 2))
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
        self.fig, self.ax = self.scdata.plot_pca_variance_explained(ylim=(0, self.yLimVar.get()), n_components=self.compVar.get())
        self.canvas = FigureCanvasTkAgg(self.fig, self)
        self.canvas.show()
        self.canvas.get_tk_widget().grid(column=1, row=1, rowspan=10, columnspan=4,sticky='NW') 
        self.currentPlot = 'pca'

        #enable buttons
        self.saveButton.config(state='normal')
        self.PCAOptions.destroy()

    def plotTSNE(self):
        self.saveButton.config(state='normal')
        self.component_menu.config(state='disabled')
        self.updateButton.config(state='disabled')
        self.setGateButton.config(state='normal')
        self.visMenu.entryconfig(6, state='normal')

        self.resetCanvas()
        self.fig, self.ax = self.scdata.plot_tsne()
        self.canvas = FigureCanvasTkAgg(self.fig, self)
        self.canvas.show()
        self.canvas.get_tk_widget().grid(column=1, row=1, rowspan=10, columnspan=4,sticky='NW') 
        self.currentPlot = 'tsne'

    def plotDM(self):
        self.saveButton.config(state='normal')
        self.component_menu.config(state='disabled')
        self.updateButton.config(state='disabled')
        self.setGateButton.config(state='disabled')
        self.visMenu.entryconfig(6, state='disabled')

        self.geometry('950x550')

        self.resetCanvas()
        self.fig, self.ax = self.scdata.plot_diffusion_components()
        self.canvas = FigureCanvasTkAgg(self.fig, self)
        self.canvas.show()
        self.canvas.get_tk_widget().grid(column=1, row=1, rowspan=10, columnspan=4,sticky='W') 
        self.currentPlot = 'dm_components'

    def showGSEAResults(self):
        self.saveButton.config(state='disabled')
        self.component_menu.config(state='normal')
        self.updateButton.config(state='normal')
        self.setGateButton.config(state='disabled')
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
        pos_text = str(self.reports[diff_component]['pos']).split('\n')
        pos_text = pos_text[1:len(pos_text)-1]
        pos_text = '\n'.join(pos_text)
        neg_text = str(self.reports[diff_component]['neg']).split('\n')
        neg_text = neg_text[1:len(neg_text)-1]
        neg_text = '\n'.join(neg_text)
        self.canvas.create_text(5, 5, anchor='nw', text='Positive correlations:\n\n', font=('Helvetica', 16, 'bold'))
        self.canvas.create_text(5, 50, anchor='nw', text=pos_text)
        self.canvas.create_text(5, 150, anchor='nw', text='Negative correlations:\n\n', font=('Helvetica', 16, 'bold'))
        self.canvas.create_text(5, 200, anchor='nw', text=neg_text)

    def plotWBOnTsne(self):
        self.saveButton.config(state='normal')
        self.component_menu.config(state='disabled')
        self.updateButton.config(state='disabled')
        self.setGateButton.config(state='disabled')
        self.visMenu.entryconfig(6, state='disabled')

        self.resetCanvas()
        self.fig, self.ax = self.wb.plot_wishbone_on_tsne()
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
            self.setGateButton.config(state='disabled')
            self.visMenu.entryconfig(6, state='disabled')

            self.resetCanvas()
            self.vals, self.fig, self.ax = self.wb.plot_marker_trajectory(self.selectedGenes)
            self.fig.set_size_inches(10, 4, forward=True)
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
        self.getGeneSelection()
        if len(self.selectedGenes) < 1:
            print('Error: must select at least one gene')

        else:
            self.saveButton.config(state='normal')
            self.component_menu.config(state='disabled')
            self.updateButton.config(state='disabled')
            self.setGateButton.config(state='disabled')
            self.visMenu.entryconfig(6, state='disabled')

            self.resetCanvas()
            self.vals, self.fig, self.ax = self.wb.plot_marker_trajectory(self.selectedGenes)
            self.fig, self.ax = self.wb.plot_marker_heatmap(self.vals)
            self.fig.set_size_inches(10, 4, forward=True)
            self.fig.tight_layout()
            self.canvas = FigureCanvasTkAgg(self.fig, self)
            self.canvas.show()
            self.canvas.get_tk_widget().grid(column=1, row=1, rowspan=10, columnspan=5, sticky='W')
            self.currentPlot = 'wishbone_marker_heatmap'

    def plotGeneExpOntSNE(self):
        self.getGeneSelection()
        if len(self.selectedGenes) < 1:
            print('Error: must select at least one gene')

        else:
            self.saveButton.config(state='normal')
            self.component_menu.config(state='disabled')
            self.updateButton.config(state='disabled')
            self.setGateButton.config(state='disabled')
            self.visMenu.entryconfig(6, state='disabled')

            self.resetCanvas()
            self.fig, self.ax = self.scdata.plot_gene_expression(self.selectedGenes)
            self.canvas = FigureCanvasTkAgg(self.fig, self)
            self.canvas.show()
            self.canvas.get_tk_widget().grid(column=1, row=1, rowspan=10, columnspan=4,sticky='W')
            self.currentPlot = 'gene_expression_tsne'
            self.geometry('950x550')

    def getGeneSelection(self):
        #popup menu to get selected genes
        self.geneSelection = tk.Toplevel()
        self.geneSelection.title("Select Genes")
        tk.Label(self.geneSelection,text=u"Genes:",fg="black",bg="white").grid(row=0)

        self.geneInput = wishbone.autocomplete_entry.AutocompleteEntry(self.genes.tolist(), self.geneSelection, listboxLength=6)
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
        #pop up for gate name
        self.gateOptions = tk.Toplevel()
        self.gateOptions.title("Create gate for start cells")
        tk.Label(self.gateOptions,text=u"Gate name:" ,fg="black",bg="white").grid(column=0, row=0)
        self.gateName = tk.StringVar()
        self.gateName.set('Gate ' + str(len(self.gates) + 1))
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
        gated_cells = self.scdata.tsne.index[gate.contains_points(self.scdata.tsne)]
        self.gates[self.gateName.get()] = gated_cells

        #replot tSNE w gate colored
        self.fig.clf()
        plt.scatter(self.scdata.tsne['x'], self.scdata.tsne['y'], s=10, edgecolors='none', color='lightgrey')
        plt.scatter(self.scdata.tsne.ix[gated_cells, 'x'], self.scdata.tsne.ix[gated_cells, 'y'], s=10, edgecolors='none')
        self.canvas.draw()

        self.setGateButton.config(state='disabled')
        self.visMenu.entryconfig(6, state='disabled')

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
    app.mainloop()

if __name__ == "__main__":
    launch()

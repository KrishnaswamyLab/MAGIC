import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
import wishbone
import os
import pandas as pd
from autocomplete_entry import AutocompleteEntry
import tkinter as tk
from tkinter import filedialog
import pickle

class wishbone_gui(tk.Tk):
    def __init__(self,parent):
        tk.Tk.__init__(self,parent)
        self.parent = parent
        self.initialize()

    def initialize(self):
        self.grid()
        self.vals = None

        #set up menu bar
        self.menubar = tk.Menu(self)
        self.fileMenu = tk.Menu(self.menubar, tearoff=0)
        self.menubar.add_cascade(label="File", menu=self.fileMenu)
        self.fileMenu.add_command(label="Load data", command=self.loadData)

        self.analysisMenu = tk.Menu(self.menubar, tearoff=0)
        self.menubar.add_cascade(label="Analysis", menu=self.analysisMenu)
        self.analysisMenu.add_command(label="Principal component analysis", state='disabled', command=self.runPCA)
        self.analysisMenu.add_command(label="tSNE", state='disabled', command=self.runTSNE)
        self.analysisMenu.add_command(label="Diffusion map", state='disabled', command=self.runDM)
        self.analysisMenu.add_command(label="Wishbone", state='disabled', command=self.runWishbone)

        self.visMenu = tk.Menu(self.menubar, tearoff=0)
        self.menubar.add_cascade(label="Visualization", menu=self.visMenu)
        self.visMenu.add_command(label="Principal component analysis", state='disabled', command=self.plotPCA)
        self.visMenu.add_command(label="tSNE", state='disabled', command=self.plotTSNE)
        self.visMenu.add_command(label="Diffusion map", state='disabled', command=self.plotDM)
        self.wishboneMenu = tk.Menu(self)
        self.visMenu.add_cascade(label="Wishbone", menu=self.wishboneMenu)
        self.wishboneMenu.add_command(label="On tSNE", state='disabled', command=self.plotWBOnTsne)
        self.wishboneMenu.add_command(label="Marker trajectory", state='disabled', command=self.plotWBMarkerTrajectory)
        self.wishboneMenu.add_command(label="Heat map", state='disabled', command=self.plotWBHeatMap)
        self.visMenu.add_command(label="Gene expression", state='disabled', command=self.plotGeneExpOntSNE)
        
        self.config(menu=self.menubar)

        #intro screen
        tk.Label(self, text=u"Wishbone", font=('Helvetica', 48), fg="black", bg="white", padx=100, pady=50).grid(row=0)
        tk.Label(self, text=u"To get started, select a data file by clicking File > Load Data", fg="black", bg="white", padx=100, pady=25).grid(row=1)

        #update
        self.grid_columnconfigure(0,weight=1)
        self.resizable(True,False)
        self.update()
        self.geometry(self.geometry())       

    def loadData(self):
        self.dataFileName = filedialog.askopenfilename()

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
        self.canvas.get_tk_widget().grid(column=1, row=1, rowspan=17, columnspan=4, sticky='NS')

        #set up visualization buttons
        tk.Label(self, text=u"Visualizations:", fg='black', bg='white').grid(column=0, row=1)
        self.PCAButton = tk.Button(self, text=u"PCA", state='disabled', command=self.plotPCA)
        self.PCAButton.grid(column=0, row=2)
        self.tSNEButton = tk.Button(self, text=u"tSNE", state='disabled', command=self.plotTSNE)
        self.tSNEButton.grid(column=0, row=3)
        self.DMButton = tk.Button(self, text=u"Diffusion map", state='disabled', command=self.plotDM)
        self.DMButton.grid(column=0, row=4)
        self.WBButton = tk.Button(self, text=u"Wishbone", state='disabled', command=self.plotWBOnTsne)
        self.WBButton.grid(column=0, row=5)
        self.geneExpButton = tk.Button(self, text=u"Gene expression", state='disabled', command=self.plotGeneExpOntSNE)
        self.geneExpButton.grid(column=0, row=6)

        #enable buttons
        self.analysisMenu.entryconfig(0, state='normal')

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
        tk.Entry(self.tsneOptions, textvariable=self.nCompVar).grid(column=1,row=0)
        tk.Button(self.tsneOptions, text=u"Run", command=self.tsneOptions.destroy).grid(column=0, columnspan=2, row=1)
        self.wait_window(self.tsneOptions)

        self.scdata.run_tsne(n_components=self.nCompVar.get())

        #enable buttons
        self.analysisMenu.entryconfig(2, state='normal')
        self.visMenu.entryconfig(1, state='normal')
        self.visMenu.entryconfig(3, state='normal')
        self.tSNEButton.config(state='normal')
        self.geneExpButton.config(state='normal')

    def runDM(self):
        self.scdata.run_diffusion_map()

        #enable buttons
        self.analysisMenu.entryconfig(3, state='normal')
        self.visMenu.entryconfig(2, state='normal')
        self.DMButton.config(state='normal')

    def runWishbone(self):
        #popup menu for wishbone options
        self.wbOptions = tk.Toplevel()
        self.wbOptions.title("Wishbone Options")

        #s
        tk.Label(self.wbOptions,text=u"Start cell:",fg="black",bg="white").grid(column=0,row=0)
        self.start = tk.StringVar()
        tk.Entry(self.wbOptions, textvariable=self.start).grid(column=1,row=0)

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
        tk.Checkbutton(self.wbOptions, text=u"Branch", variable=self.branch).grid(column=0, row=4, columnspan=2)

        tk.Button(self.wbOptions, text=u"Run", command=self.wbOptions.destroy).grid(column=0, columnspan=2, row=5)
        self.wait_window(self.wbOptions)

        self.wb = wishbone.wb.Wishbone(self.scdata)
        self.wb.run_wishbone(start_cell=self.start.get(), k=self.k.get(), components_list=[int(comp) for comp in self.compList.get().split(',')], num_waypoints=self.numWaypoints.get())

        #enable buttons
        self.wishboneMenu.entryconfig(0, state='normal')
        self.wishboneMenu.entryconfig(1, state='normal')
        self.WBButton.config(state='normal')

    def plotPCA(self):
        #pop up for # components
        self.PCAOptions = tk.Toplevel()
        self.PCAOptions.title("PCA Plot Options")
        tk.Label(self.PCAOptions,text=u"Max variance explained (ylim):",fg="black",bg="white").grid(column=0, row=0)
        self.yLimVar = tk.DoubleVar()
        tk.Entry(self.PCAOptions, textvariable=self.yLimVar).grid(column=1,row=0)
        tk.Label(self.PCAOptions, text=u"Number of components:", fg='black', bg='white').grid(column=0, row=1)
        self.compVar = tk.IntVar()
        tk.Entry(self.PCAOptions, textvariable=self.compVar).grid(column=1, row=1)
        tk.Button(self.PCAOptions, text=u"Plot", command=self.PCAOptions.destroy).grid(column=0, columnspan=2, row=2)
        self.wait_window(self.PCAOptions)

        self.fig.clf()
        self.fig, self.ax = self.scdata.plot_pca_variance_explained(ylim=(0, self.yLimVar.get()), n_components=self.compVar.get())
        self.canvas = FigureCanvasTkAgg(self.fig, self)
        self.canvas.show()
        self.canvas.get_tk_widget().grid(column=1, row=1, rowspan=17, columnspan=4, sticky='NS') 

    def plotTSNE(self):
        self.fig.clf()
        self.fig, self.ax = self.scdata.plot_tsne()
        self.canvas = FigureCanvasTkAgg(self.fig, self)
        self.canvas.show()
        self.canvas.get_tk_widget().grid(column=1, row=1, rowspan=17, columnspan=4, sticky='NS') 

    def plotDM(self):
        self.fig.clf()
        self.canvas.get_tk_widget().grid_forget()
        self.fig, self.ax = self.scdata.plot_diffusion_components()
        self.canvas = FigureCanvasTkAgg(self.fig, self)
        self.canvas.show()
        self.canvas.get_tk_widget().grid(column=1, row=1, rowspan=17, columnspan=4, sticky='NS') 

    def plotWBOnTsne(self):
        self.fig.clf()
        self.fig, self.ax = self.wb.plot_wishbone_on_tsne()
        self.canvas = FigureCanvasTkAgg(self.fig, self)
        self.canvas.show()
        self.canvas.get_tk_widget().grid(column=1, row=1, rowspan=17, columnspan=4, sticky='NS')

    def plotWBMarkerTrajectory(self):
        self.getGeneSelection()
        if len(self.selectedGenes) < 1:
            print('Error: must select at least one gene')
        self.fig.clf()
        self.vals, self.fig, self.ax = self.wb.plot_marker_trajectory(self.selectedGenes)
        self.canvas = FigureCanvasTkAgg(self.fig, self)
        self.canvas.show()
        self.canvas.get_tk_widget().grid(column=1, row=1, rowspan=17, columnspan=4, sticky='NS')
        
        #enable buttons
        self.wishboneMenu.entryconfig(2, state='normal')

    def plotWBHeatMap(self):
        self.fig.clf()
        self.fig, self.ax = self.wb.plot_marker_heatmap(self.vals)
        self.canvas = FigureCanvasTkAgg(self.fig, self)
        self.canvas.show()
        self.canvas.get_tk_widget().grid(column=1, row=1, rowspan=17, columnspan=4, sticky='NS')

    def plotGeneExpOntSNE(self):
        self.getGeneSelection()
        if len(self.selectedGenes) < 1:
            print('Error: must select at least one gene')
        self.fig.clf()
        self.fig, self.ax = self.scdata.plot_gene_expression(self.selectedGenes)
        self.canvas = FigureCanvasTkAgg(self.fig, self)
        self.canvas.show()
        self.canvas.get_tk_widget().grid(column=1, row=1, rowspan=17, columnspan=4, sticky='NS')

    def getGeneSelection(self):
        #popup menu to get selected genes
        self.geneSelection = tk.Toplevel()
        self.geneSelection.title("Select Genes")
        tk.Label(self.geneSelection,text=u"Genes:",fg="black",bg="white").grid(row=0)

        self.geneInput = AutocompleteEntry(self.genes.tolist(), self.geneSelection, listboxLength=6)
        self.geneInput.grid(row=1)
        self.geneInput.bind('<Return>', self.AddToSelected)

        self.geneSelectBox = tk.Listbox(self.geneSelection, selectmode=tk.EXTENDED)
        self.geneSelectBox.grid(row=2, rowspan=10)
        self.geneSelectBox.bind('<BackSpace>', self.DeleteSelected)
        self.selectedGenes = []

        tk.Button(self.geneSelection, text=u"Use selected genes", command=self.geneSelection.destroy).grid(row=12)
        self.wait_window(self.geneSelection)
    
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

if __name__ == "__main__":
    app = wishbone_gui(None)
    app.title('Wishbone')
    app.mainloop()

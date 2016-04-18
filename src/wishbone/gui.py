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

        #k
        label = tk.Label(self,text=u"k", anchor="w",fg="black",bg="white")
        label.grid(column=0,row=0,sticky='EW')

        self.kEntryVariable = tk.StringVar()
        self.kentry = tk.Entry(self, textvariable=self.kEntryVariable)
        self.kentry.grid(column=1,row=0,sticky='EW')
        self.kEntryVariable.set(u"15")

        #l
        label = tk.Label(self,text=u"l", anchor="w",fg="black",bg="white")
        label.grid(column=0,row=1,sticky='EW')

        self.lEntryVariable = tk.StringVar()
        self.lentry = tk.Entry(self, textvariable=self.lEntryVariable)
        self.lentry.grid(column=1,row=1,sticky='W')
        self.lEntryVariable.set(u"15")

        #s
        label = tk.Label(self,text=u"s", anchor="w",fg="black",bg="white")
        label.grid(column=0,row=2,sticky='EW')

        self.sEntryVariable = tk.StringVar()
        self.sentry = tk.Entry(self, textvariable=self.sEntryVariable)
        self.sentry.grid(column=1,row=2,sticky='W')
        self.sEntryVariable.set(u"1")

        #num graphs
        label = tk.Label(self,text=u"Number of graphs", anchor="w",fg="black",bg="white")
        label.grid(column=0,row=3,sticky='EW')

        self.gEntryVariable = tk.StringVar()
        self.gentry = tk.Entry(self, textvariable=self.gEntryVariable)
        self.gentry.grid(column=1,row=3,sticky='W')
        self.gEntryVariable.set(u"1")

        #num waypoints
        label = tk.Label(self,text=u"Number of waypoints", anchor="w",fg="black",bg="white")
        label.grid(column=0,row=4,sticky='EW')

        self.wEntryVariable = tk.StringVar()
        self.wentry = tk.Entry(self, textvariable=self.wEntryVariable)
        self.wentry.grid(column=1,row=4,sticky='W')
        self.wEntryVariable.set(u"150")


        #dm eigs
        label = tk.Label(self,text=u"# of eigen vectors (for diffusion map)", anchor="w",fg="black",bg="white")
        label.grid(column=0,row=5,sticky='EW')

        self.dmEntryVariable = tk.StringVar()
        self.dmentry = tk.Entry(self, textvariable=self.dmEntryVariable)
        self.dmentry.grid(column=1,row=5,sticky='W')
        self.dmEntryVariable.set(u"4")

        #verbose
        label = tk.Label(self,text=u"Verbose", anchor="w",fg="black",bg="white")
        label.grid(column=0,row=6,sticky='EW')

        self.vEntryVariable = tk.StringVar()
        self.ventry = tk.Entry(self, textvariable=self.vEntryVariable)
        self.ventry.grid(column=1,row=6,sticky='W')
        self.vEntryVariable.set(u"1")

        #metric
        label = tk.Label(self,text=u"Metric", anchor="w",fg="black",bg="white")
        label.grid(column=0,row=7,sticky='EW')

        self.mEntryVariable = tk.StringVar()
        self.mentry = tk.Entry(self, textvariable=self.mEntryVariable)
        self.mentry.grid(column=1,row=7,sticky='W')
        self.mEntryVariable.set(u"euclidean")

        #voting scheme
        label = tk.Label(self,text=u"Voting scheme", anchor="w",fg="black",bg="white")
        label.grid(column=0,row=8,sticky='EW')

        self.vsEntryVariable = tk.StringVar()
        self.vsentry = tk.Entry(self, textvariable=self.vsEntryVariable)
        self.vsentry.grid(column=1,row=8,sticky='W')
        self.vsEntryVariable.set(u"exponential")

        #branch
        label = tk.Label(self,text=u"Branch", anchor="w",fg="black",bg="white")
        label.grid(column=0,row=9,sticky='EW')

        self.bEntryVariable = tk.StringVar()
        self.bentry = tk.Entry(self, textvariable=self.bEntryVariable)
        self.bentry.grid(column=1,row=9,sticky='W')
        self.bEntryVariable.set(u"1")

        #band_sample
        label = tk.Label(self,text=u"Band sample", anchor="w",fg="black",bg="white")
        label.grid(column=0,row=10,sticky='EW')

        self.bsEntryVariable = tk.StringVar()
        self.bsentry = tk.Entry(self, textvariable=self.bsEntryVariable)
        self.bsentry.grid(column=1,row=10,sticky='W')
        self.bsEntryVariable.set(u"0")

        #partial order
        label = tk.Label(self,text=u"Partial order", anchor="w",fg="black",bg="white")
        label.grid(column=0,row=11,sticky='EW')

        self.poEntryVariable = tk.StringVar()
        self.poentry = tk.Entry(self, textvariable=self.poEntryVariable)
        self.poentry.grid(column=1,row=11,sticky='W')
        self.poEntryVariable.set(u"")

        #flock waypoints
        label = tk.Label(self,text=u"Flock waypoints", anchor="w",fg="black",bg="white")
        label.grid(column=0,row=12,sticky='EW')

        self.fwEntryVariable = tk.StringVar()
        self.fwentry = tk.Entry(self, textvariable=self.fwEntryVariable)
        self.fwentry.grid(column=1,row=12,sticky='W')
        self.fwEntryVariable.set(u"2")

        #plot data
        label = tk.Label(self,text=u"Plot data", anchor="w",fg="black",bg="white")
        label.grid(column=0,row=13,sticky='EW')

        self.pdEntryVariable = tk.StringVar()
        self.pdentry = tk.Entry(self, textvariable=self.pdEntryVariable)
        self.pdentry.grid(column=1,row=13,sticky='W')
        self.pdEntryVariable.set(u"")

        #search connected components
        label = tk.Label(self,text=u"Search connected components", anchor="w",fg="black",bg="white")
        label.grid(column=0,row=14,sticky='EW')

        self.sccEntryVariable = tk.StringVar()
        self.sccentry = tk.Entry(self, textvariable=self.sccEntryVariable)
        self.sccentry.grid(column=1,row=14,sticky='W')
        self.sccEntryVariable.set(u"1")


        #select data
        dataButton = tk.Button(self, text=u"Select data file", command=self.OnDataButtonClick)
        dataButton.grid(column=0, row=15)

        #select waypoints
        waypointsButton = tk.Button(self, text=u"Select landmarks file", command=self.OnWaypointsButtonClick)
        waypointsButton.grid(column=1, row=15)

        #calculate trajectory button
        button = tk.Button(self,text=u"Run wishbone", command=self.OnCalculateTrajectoryButtonClick)
        button.grid(column=0,row=17, columnspan=2, sticky='S')

        #Visualizations
        label = tk.Label(self,text=u"Visualize:", anchor="w",fg="black",bg="white")
        label.grid(column=2,row=17,sticky='EW')

        #diffusion map button
        dmButton = tk.Button(self, text=u"Diffusion map", command=self.ShowDiffusionMap)
        dmButton.grid(column=4, row=17)

        #trajectory button
        trajButton = tk.Button(self, text=u"Wishbone on tSNE", command=self.PlotWishboneOnTsne)
        trajButton.grid(column=3, row=18)

        #tSNE button
        tsneButton = tk.Button(self, text=u"tSNE", command=self.PlotTsne)
        tsneButton.grid(column=3, row=17)

        #diffusion eigen vectors
        evButton = tk.Button(self, text=u"Diffusion eigen vectors", command=self.PlotDiffusionEigenVectors)
        evButton.grid(column=5, row=17)

        #plot marker trajectory
        mtButton = tk.Button(self, text=u"Plot marker trajectory", command=self.PlotMarkerTrajectory)
        mtButton.grid(column=4, row=18)

        #plot marker heatmap
        mhButton = tk.Button(self, text=u"Plot marker heatmap", command=self.PlotMarkerHeatmap)
        mhButton.grid(column=5, row=18)

        #plot derivatives
        derivButton = tk.Button(self, text=u"Plot derivatives", command=self.PlotDerivatives)
        derivButton.grid(column=3, row=19)
        
        #set up canvas for plots
        self.fig, self.ax = wishbone.wb.get_fig()
        self.canvas = FigureCanvasTkAgg(self.fig, self)
        self.canvas.show()
        self.canvas.get_tk_widget().grid(column=2, row=0, rowspan=17, columnspan=4, sticky='NS')

        #set up listbox for genes
        label = tk.Label(self,text=u"Genes:", anchor="w",fg="black",bg="white")
        label.grid(column=6,row=0,sticky='NW')

        self.geneSelectBox = tk.Listbox(self, selectmode=tk.EXTENDED)
        self.geneSelectBox.grid(column=6, row=2, rowspan=16, sticky='NS')
        self.geneSelectBox.bind('<BackSpace>', self.DeleteSelected)
        self.selectedGenes = []
        gsButton = tk.Button(self, text=u"Plot gene expression", command=self.PlotSelectedGenes)
        gsButton.grid(column=6, row=17)
        gcButton = tk.Button(self, text=u"Compare expression of 2 markers", command=self.CompareGeneExpression)
        gcButton.grid(column=6, row=18)

        #update
        self.grid_columnconfigure(0,weight=1)
        self.resizable(True,False)
        self.update()
        self.geometry(self.geometry())       

    def OnDataButtonClick(self):
        self.dataFileName = filedialog.askopenfilename()
        # Load sample data
        self.scdata = wishbone.wb.SCData.from_fcs(os.path.expanduser('~/.wishbone/data/sample_masscyt.fcs'), 
            cofactor=None)
        self.genes = self.scdata.data.columns.values

        #gene selection setup
        self.geneInput = AutocompleteEntry(self.genes.tolist(), self, listboxLength=6)
        self.geneInput.grid(column=6, row=1, sticky='NS')
        self.geneInput.bind('<Return>', self.AddToSelected)

    
    def OnWaypointsButtonClick(self):
        self.waypointsFileName = filedialog.askopenfilename()
        self.waypoints = pd.DataFrame.from_csv(os.path.expanduser(self.waypointsFileName)).iloc[:, 0]
        self.waypoints = list(self.waypoints)

    def OnCalculateTrajectoryButtonClick(self):
        print(self.kEntryVariable.get())
        print(self.lEntryVariable.get())
        print(self.sEntryVariable.get())
        self.HideTrajectoryOptions()

        # Run tSNE
        # self.scdata.run_tsne()
        f = open('tsne', 'rb')
        self.scdata.tsne = pickle.load(f)
        # Run diffusion maps
        self.scdata.run_diffusion_map()

        # self.trajectory, self.waypoints, self.branches, self.bas = wishbone.wishbone(self.data, 
        #                                                                     k=int(self.kEntryVariable.get()), 
        #                                                                     l=int(self.lEntryVariable.get()), 
        #                                                                     s=int(self.sEntryVariable.get())-1, 
        #                                                                     num_graphs=int(self.gEntryVariable.get()), 
        #                                                                     num_waypoints=[int(self.wEntryVariable.get())] if self.wEntryVariable.get().isdigit() else self.landmarks, 
        #                                                                     dm_eigs=int(self.dmEntryVariable.get()),
        #                                                                     verbose=int(self.vEntryVariable.get()),
        #                                                                     metric=self.mEntryVariable.get(), 
        #                                                                     voting_scheme=self.vsEntryVariable.get(),
        #                                                                     branch=int(self.bEntryVariable.get()), 
        #                                                                     band_sample=int(self.bsEntryVariable.get()),
        #                                                                     partial_order= [int(i) for i in self.poEntryVariable.get().split(',')] if len(self.poEntryVariable.get()) > 0 else [],
        #                                                                     flock_waypoints=int(self.fwEntryVariable.get()),
        #                                                                     plot_data= [int(i) for i in self.pdEntryVariable.get().split(',')] if len(self.pdEntryVariable.get()) > 0 else [],
        #                                                                     search_connected_components=int(self.sccEntryVariable.get()))
        # self.wb = wishbone.wb.Wishbone(self.scdata)
        # self.wb.run_wishbone(int(self.sEntryVariable.get()), components_list=[1, 2, 3], num_waypoints=self.waypoints)
        f2 = open('wishb', 'rb')
        self.wb = pickle.load(f2)

        self.fig.clf()
        # self.canvas.get_tk_widget.delete('all')
        self.fig, self.ax = self.wb.plot_wishbone_on_tsne()
        self.canvas = FigureCanvasTkAgg(self.fig, self)
        self.canvas.show()
        self.canvas.get_tk_widget().grid(column=2, row=0, rowspan=17, columnspan=4, sticky='NS') 

    def ShowDiffusionMap(self):
        self.fig.clf()
        # self.canvas.get_tk_widget.delete('all')
        self.fig, self.ax = self.scdata.plot_diffusion_components()
        self.canvas = FigureCanvasTkAgg(self.fig, self)
        self.canvas.show()
        self.canvas.get_tk_widget().grid(column=2, row=0, rowspan=17, columnspan=4, sticky='NS') 

    def PlotWishboneOnTsne(self):
        self.fig.clf()
        # self.canvas.get_tk_widget.delete('all')
        self.fig, self.ax = self.wb.plot_wishbone_on_tsne()
        self.canvas = FigureCanvasTkAgg(self.fig, self)
        self.canvas.show()
        self.canvas.get_tk_widget().grid(column=2, row=0, rowspan=17, columnspan=4, sticky='NS') 

    def PlotTsne(self):
        self.fig.clf()
        # self.canvas.get_tk_widget.delete('all')
        self.fig, self.ax = self.scdata.plot_tsne()
        self.canvas = FigureCanvasTkAgg(self.fig, self)
        self.canvas.show()
        self.canvas.get_tk_widget().grid(column=2, row=0, rowspan=17, columnspan=4, sticky='NS') 

    def PlotDiffusionEigenVectors(self):
        self.fig.clf()
        # self.canvas.get_tk_widget.delete('all')
        self.fig, self.ax = self.scdata.plot_diffusion_eigen_vectors()
        self.canvas = FigureCanvasTkAgg(self.fig, self)
        self.canvas.show()
        self.canvas.get_tk_widget().grid(column=2, row=0, rowspan=17, columnspan=4, sticky='NS') 
        
    def CompareGeneExpression(self):
        selected = self.geneSelectBox.curselection()
        if len(selected) != 2:
            print('Error: select exactly two genes to plot')
        self.fig.clf()
        a = self.fig.add_subplot(111)
        a.scatter(self.scdata.data[self.selectedGenes[selected[0]]], self.scdata.data[self.selectedGenes[selected[1]]], s=10, edgecolors='none')
        a.set_xlim([0, 6])
        a.set_ylim([0, 6])
        a.set_xlabel(self.selectedGenes[selected[0]])
        a.set_ylabel(self.selectedGenes[selected[1]])
        self.canvas.draw()

    def PlotSelectedGenes(self):
        selected = self.geneSelectBox.curselection()
        if len(selected) < 1:
            print('Error: must select at least one gene')
        self.fig.clf()
        self.fig, self.ax = self.scdata.plot_gene_expression([self.selectedGenes[i] for i in selected])
        self.canvas = FigureCanvasTkAgg(self.fig, self)
        self.canvas.show()
        self.canvas.get_tk_widget().grid(column=2, row=0, rowspan=17, columnspan=4, sticky='NS')

    def PlotMarkerTrajectory(self):
        selected = self.geneSelectBox.curselection()
        if len(selected) < 1:
            print('Error: must select at least one gene')
        self.fig.clf()
        self.vals, self.fig, self.ax = self.wb.plot_marker_trajectory([self.selectedGenes[i] for i in selected], 
                    smoothing_factor=1.5, show_variance=True)
        self.canvas = FigureCanvasTkAgg(self.fig, self)
        self.canvas.show()
        self.canvas.get_tk_widget().grid(column=2, row=0, rowspan=17, columnspan=4, sticky='NS')

    def PlotMarkerHeatmap(self):
        if self.vals == None:
            print('Error: must plot marker trajectory first')
        self.fig.clf()
        self.fig, self.ax = self.wb.plot_marker_heatmap(self.vals)
        self.canvas = FigureCanvasTkAgg(self.fig, self)
        self.canvas.show()
        self.canvas.get_tk_widget().grid(column=2, row=0, rowspan=17, columnspan=4, sticky='NS')

    def PlotDerivatives(self):
        if self.vals == None:
            print('Error: must plot marker trajectory first')
        self.fig.clf()
        self.fig, self.ax = self.wb.plot_derivatives(self.vals)
        self.canvas = FigureCanvasTkAgg(self.fig, self)
        self.canvas.show()
        self.canvas.get_tk_widget().grid(column=2, row=0, rowspan=17, columnspan=4, sticky='NS')


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

    def HideTrajectoryOptions(self):
        for item in self.grid_slaves():
            if int(item.grid_info()['column']) < 2:
                item.grid_forget()
        self.grid_columnconfigure(0, minsize=0)
        self.grid_columnconfigure(1, minsize=0)


if __name__ == "__main__":
    app = wishbone_gui(None)
    app.title('Wishbone')
    app.mainloop()

#!/usr/bin/env python
# coding: utf-8


#transmon arch and runs
from qutip import *
import numpy as np
import matplotlib.pyplot as plt
import math
from tqdm import tqdm
from eigenfns import *


class transmon_arch:
    #init variables
    def __init__(self, N, Nlevel,error=0):
        #Energy scale: Ej = 20Ghz +- 10Ghz, Ec = 250MHz, T = (5-20)MHz 
        #ratio = 200/2.5/.1 
        self.N = N
        self.Nlevel = Nlevel
        self.Ej =  np.absolute(np.random.normal(200,50, N))
        self.Ec = np.absolute(np.random.normal(2.5,0.2, N))
        self.error = error #error width
        self.id = 0
        self.a = 0
        self.na = 0
        self.H_const = 0
        self.operators()
        self.H_c()
        
    def operators(self):
        Neye = []
        id = [] #tensored identiry operator
        a = [] #list of tensored destruction operators
        na = [] #list of number operators
        for i in  range(self.N):
            Neye.append(qeye(self.Nlevel))
        for i in range(self.N):
            temp = Neye.copy()
            temp[i] = destroy(self.Nlevel)
            a.append(tensor(temp))
        for i in range(self.N):
            na.append(a[i].dag()*a[i])
        id = tensor(Neye)

        self.id = id
        self.a = a
        self.na = na

    def H_c(self):
        #constant part of hamiltonian
        N = self.N
        #initialize v ,t
        v = np.zeros(N)
        for i in range(N):
            v[i] = math.sqrt(8*self.Ej[i]*self.Ec[i])
        
        H1 = 0
        H2 = 0
        for i in range(N):
            e = abs(np.random.normal(1, self.error))
            H1 = H1 + self.na[i]*v[i]*e
            H2 = H2 + self.Ec[i]*(self.na[i]*(self.na[i]+1))*e
            
        H2 = -0.5*H2        
        self.H_const = H1 + H2
    
    def H_int(self,T,RWA = False,t_matrix = 0):
        N = self.N
        Ec = self.Ec
        Ej = self.Ej
        if t_matrix == 0:
            t = np.zeros((N,N))
            for i in range(N-1):
                t[i][i+1] = (T/(4*(2*((Ec[i]*Ec[i+1])**1/float(2)))**(1/float(2))))*((Ej[i]*Ej[i+1])**(1/float(4)))
        else:
            t = t_matrix
        
        H3 = 0
        if RWA:
            for i in range(N):
                for j in range(N):
                    H3 = H3 + t[i][j]*(self.a[i]*self.a[j].dag()+self.a[i].dag()*self.a[j])
            
        if not RWA:
            for i in range(N):
                for j in range(N):
                    H3 = H3 + t[i][j]*(self.a[i] + self.a[i].dag())*(self.a[j].dag()+self.a[j])
        return H3
    
   
        
class run_data:
    def __init__(self, t_arch,tlist,RWA=False, update_error = False, global_error = 0):
        self.t_arch = t_arch
        self.tlist = tlist
        self.spectrum_arr = 0
        self.cspectrum_arr = 0
        self.evec = []
        self.c_evec = []
        self.RWA = RWA
        self.update_error = update_error
        self.global_error = global_error

    def spectrum_cal(self):
        spectrum_arr = []
        for t in tqdm(self.tlist):
            H3 = self.t_arch.H_int(t,RWA = self.RWA)
            if self.update_error:
                self.t_arch.H_c()
            H = abs(np.random.normal(1, self.global_error))*self.t_arch.H_const + H3
            evals = H.eigenenergies()
            spectrum_arr.append([t,evals])

        self.spectrum_arr = spectrum_arr

    def evec_cal(self):
        evec_arr = []
        spectrum_arr = []
        for t in tqdm(self.tlist):
            H3 = self.t_arch.H_int(t,RWA = self.RWA)
            if self.update_error:
                self.t_arch.H_c()
            H = abs(np.random.normal(1, self.global_error))*self.t_arch.H_const + H3
            eigen_vals, eigen_vecs = H.eigenstates()
            spectrum_arr.append([t,eigen_vals])
            evec_arr.append(eigen_vecs)

        self.evec = evec_arr
        self.spectrum_arr = spectrum_arr

        
    def cspectrum_cal(self):
        cspectrum_arr = []
        c_evec_arr = []
        for t in tqdm(self.tlist):
            H3 = self.t_arch.H_int(t,RWA = self.RWA)
            if self.update_error:
                self.t_arch.H_c()
            H = abs(np.random.normal(1, self.global_error))*self.t_arch.H_const + H3
            evals,evecs, labels = cstates_label(H, self.t_arch.na)
            cspectrum_arr.append([t, evals, labels])
            c_evec_arr.append(evecs)
        
        self.cspectrum_arr = cspectrum_arr
        self.c_evec = c_evec_arr
        
    def cal_all(self):
        self.evec_cal()
        self.cspectrum_cal()
        return
        

class plot_spectrum:    
    def __init__(self, N ,Nlevel,spectrum_arr):
            self.title = str(str(N) + " tranmons " + str(Nlevel)+ " level systems")
            self.x_label = str("T value")
            self.y_label = str("eigenvalues")
            self.size = 10
            self.c = 'b'
            self.x = 0
            self.y = 0
            self.label = 0
            self.spectrum_arr = spectrum_arr
            self.spectrum_data()

    def spectrum_data(self):
        x = []
        y = []
        label = []
        fig = 0
        if len(self.spectrum_arr[0]) == 2:
            for i in self.spectrum_arr:
                for e in i[1]:
                    x.append(i[0])
                    y.append(e)
            self.x = x
            self.y = y
        elif len(self.spectrum_arr[0]) == 3:
            for i in self.spectrum_arr:
                for e,l in zip(i[1],i[2]):
                    x.append(i[0])
                    y.append(e)
                    label.append(l)
            self.x = x
            self.y = y
            self.label = label
    def label_t(self,t=0):
        label_t_arr = []
        label_t_y = []
        for x_temp,y_temp,l in zip(self.x,self.y,self.label):
            if x_temp == t:
                label_t_arr.append(l)
                label_t_y.append(y_temp)
        return label_t_arr , label_t_y

    def view_fig(self, y_min = -1, y_max = -1, use_label=False, title = -1):
        self.fig = plt.figure(figsize=(self.size,self.size))
        if not((y_min ==-1 or y_max ==-1) or (y_min > y_max) ):
            plt.ylim(y_min, y_max)
        plt.scatter(self.x,self.y,s = 0.5, color = self.c)
        if use_label:
            label_t_arr , label_t_y = self.label_t()
            for l,y in zip(label_t_arr , label_t_y):
                plt.text(0.0,y,str(l))
        if title == -1:
            title = self.title
        plt.title(title)
        plt.xlabel(self.x_label)
        plt.ylabel(self.y_label)
        
def mult_plot(*plots, size = 10, y_min = -1, y_max = -1, use_label = False, title = "overlapped",x_label = "T value", y_label = "eigenvalues"):
    fig = plt.figure(figsize=(size,size))
    if not((y_min ==-1 or y_max ==-1) or (y_min > y_max) ):
        plt.ylim(y_min, y_max)
    
    for p in plots:
        plt.scatter(p.x,p.y,s = 0.5, color = p.c, alpha=0.5)
        if use_label:
            label_t_arr , label_t_y = p.label_t()
            for l,y in zip(label_t_arr , label_t_y):
                plt.text(0.0,y,str(l))
    plt.title(str(title))
    #assuming all have same x,y labels
    plt.xlabel(x_label)
    plt.ylabel(y_label)

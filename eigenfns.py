#!/usr/bin/env python
# coding: utf-8

from qutip import *

#functions on eigenvalues
def check_val(expect_val):
    #function for checking if a^\dag a lies in a certain range 
    if temp > 1.1 or (temp <0.9 and temp >0.1):
        return 1
    else:
        retun 0
    
def cstates(H,na):
    #returns eigenvecs and eigenvalues of states that are computational
    #i.e have at max 1 excitation of transmons
    eigen_vals, eigen_vecs = H.eigenstates()
    c_val = []
    c_vec = []
    for i in range(len(eigen_vals)):
        c = 1
        vec = eigen_vecs[i]
        for nai in na:
            if expect(nai,vec) > 1.1:
                c = 0
                break
        if c==1 :
            c_val.append(eigen_vals[i])
            c_vec.append(eigen_vecs[i])

    return c_val,c_vec

def cstates_label(H,na):
    #returns evals, estates and occupation string of all cstates
    labels= []
    eigen_vals, eigen_vecs = H.eigenstates()
    c_val = []
    c_vec = []
    for i in range(len(eigen_vals)):
        if len(c_val) > 2**len(na):
            break
        c = 1
        blabel = []
        vec = eigen_vecs[i]
        for nai in na:
            temp = expect(nai,vec)
            if temp > 1.1:
                c = 0
                break
            blabel.append(temp)
        
        if c==1 :
            c_val.append(eigen_vals[i])
            c_vec.append(eigen_vecs[i])
            labels.append(blabel)

    return c_val,c_vec ,labels

def band_states(H,na,band_no):
    #returns all eigenvectors in given band no
    eigen_vals, eigen_vecs = H.eigenstates()
    b_val = []
    b_vec = []
    total_na = 0
    error_tol = 0.2

    for nai in na:
        total_na = total_na + nai    

    for i in range(len(eigen_vals)):
        vec = eigen_vecs[i]
        c = expect(total_na,vec)
        if c < (band_no + error_tol) and c > (band_no - error_tol)  :
            b_val.append(eigen_vals[i])
            b_vec.append(eigen_vecs[i])

    return b_val,b_vec


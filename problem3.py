# -*- coding: utf-8 -*-
"""
Created on Sat Mar 28 11:19:08 2015

@author: kevinhe
"""

import math
import random
import pylab

def rg(n):
    return range(2*n+1)
    
def gene0(n):
    m=[] 
    for i in rg(n):
        t=[]
        for j in rg(n):
            if (i-n)**2+(j-n)**2>n**2:
                t+=[-1]
            else:
                t+=[0]
        m+=[t]
    m[n][n]=1
    return m

def bpoints(m,n):
    ps=[]
    for i in range(1,2*n):
        for j in range(1,2*n):
            if m[i][j]==0:
                if -1 in adjvs(m,(i,j)):
                    ps+=[(i,j)]
    ps+=[(0,n),(n,0),(n,2*n),(2*n,n)]
    return ps

def adj(m,loc):
    (i1,j1)=loc
    n=(len(m)-1)/2
    vs=[(i1+1,j1),(i1-1,j1),(i1,j1+1),(i1,j1-1)]
    vr=[]
    for (p,q) in vs:
        if (p in rg(n)) and (q in rg(n)):
            if m[p][q]!=-1:
                vr+=[(p,q)]
    return vr

def adjvs(m,loc):
    return [m[p][q] for (p,q) in adj(m,loc)]
    
def notjoin(m,loc):
    if 1 in adjvs(m,loc):
        return False
    return True
    
def main(d):
    n=100/d
    m=gene0(n)
    ps=bpoints(m,n)
    while notcomp(m,n):
        (i,j)=random.choice(ps)
        m[i][j]=1
        while notjoin(m,(i,j)):
            (ip,jp)=random.choice(adj(m,(i,j)))
            (m[i][j],m[ip][jp])=(m[ip][jp],m[i][j])
            (i,j)=(ip,jp)
        pplot(m)
    return m
    
def notcomp(m,n):
    for (p,q) in bpoints(m,n):
        if m[p][q]==1:
            return False
    return True

def pplot(m):
    pylab.contour(m)
    pylab.show()
    
main(25)

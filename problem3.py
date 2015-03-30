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
    '''
    generate the initial grid 
    with -1 representing points beyond boundary
    and 0 the empty grid sites
    and 1 the occupied sites
    '''
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
    '''
    find all boundary points
    '''
    ps=[]
    for i in range(1,2*n):
        for j in range(1,2*n):
            if m[i][j]!=-1:
                if ((i-n)**2+(j+1-n)**2)>n**2 or ((i+1-n)**2+(j-n)**2)>n**2 or ((i-n)**2+(j-1-n)**2)>n**2 or ((i-1-n)**2+(j-n)**2)>n**2 :
                    ps+=[(i,j)]
    ps+=[(0,n),(n,0),(n,2*n),(2*n,n)]
    return ps

def adj(m,loc):
    '''
    loc is the tuple storing the coordinate
    return the coordinates of all adjacent non-boundary points around loc
    '''
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
    '''
    values of all adjacent points 
    '''
    return [m[p][q] for (p,q) in adj(m,loc)]
    
def notjoin(m,loc):
    '''
    check if the moving monomer has joined the cluster
    '''
    if 1 in adjvs(m,loc):
        return False
    return True
    
def main(d):
    '''
    main function
    takes d, which is the grid size, as argument 
    '''
    n=100/d
    m=gene0(n)
    maxr=0
    maxrs=[]
    mass=1
    masses=[]
    while notcomp(m,n):
        (i,j)=random.choice(bpoints(m,n))
        m[i][j]=1
        
        while notjoin(m,(i,j)):
            (ip,jp)=random.choice(adj(m,(i,j)))
            (m[i][j],m[ip][jp])=(m[ip][jp],m[i][j])
            (i,j)=(ip,jp)
           
        if ((i-n)**2+(j-n)**2)*d**2>maxr**2:
            maxr=math.sqrt(((i-n)**2+(j-n)**2))*d
            
        mass=mass+1
        
        maxrs=maxrs+[math.log(maxr)]
        masses+=[math.log(mass)]
        

    pylab.plot(masses,maxrs)
    pylab.xlabel('number of molecules')
    pylab.ylabel('maximum radius')
    
    #pylab.close()
    #pplot(m)
    return m
    
def notcomp(m,n):
    '''
    check if the cluster is complete, which means the cluster has touched the boundary
    '''
    for (p,q) in bpoints(m,n):
        if m[p][q]==1:
            return False
    return True

def pplot(m):
    pylab.contour(m,levels=[-0.5,0.5])
    pylab.show()

for l in range(10):
    main(5)
pylab.show()



# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

#################### IMPORTING ALL PACKAGES WE USED ###########################

import math
from random import *
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

################## LIST OF MATERIALS' PROPERITES ##############################

Mwater=18 * 10 **(-3)
Dwater=1000
sigma_abs_Water = 0.665 * 10**-28
sigma_sct_Water = 103 * 10**-28

Mlead=207* 10 **(-3)
Dlead=11340
sigma_abs_Lead = 0.158 * 10**-28
sigma_sct_Lead = 11.2 * 10**-28

Mgraphite=12* 10 **(-3)
Dgraphite=2266
sigma_abs_Graphite = 0.0045 * 10**-28
sigma_sct_Graphite = 4.74 * 10**-28

Mbe=9*10**-3#Beryllium
Dbe=1850
sigma_sct_Be=7.63*10**-28
sigma_abs_Be=0.0076*10**-28

################### DEFINING ALL CONSTANTS ####################################

M=Mgraphite
D=Dgraphite
Na=6.02214086 * 10**23 #/mol
A_abs=sigma_abs_Graphite #m^2
A_sct=sigma_sct_Graphite #m^2
A=A_abs+A_sct
mfp=M/(D*Na*A)

M2=Mwater
D2=Dwater
A2_abs=sigma_abs_Water
A2_sct=sigma_sct_Water
A2=A2_abs+A2_sct
mfp2=M2/(D2*Na*A2)

############ DEFINING PROBABILITY OF ABSORPTION AND MATERIAL SIZE(S) ##########

woodcock = False

Pabs1=A_abs/A
x0=10.0
y0=10.0
z0=100

Pabs2=A2_abs/A2
x2=10
y2=10
z2=10000


################## FUNCTION TO CALCULATE THE PATH DISTANCE ####################

def calDistance():
    rand = random()
    r = (-mfp)*math.log(1-rand)
    return r

################ OLD UNIFORM SCATTERING ANLGE FUNCTION ########################
    
#def genAngle():
#    phi = np.random.uniform(0,np.pi*2)
#    costheta = np.random.uniform(-1,1)
#    theta = np.arccos( costheta )
#    return (theta,phi)

#################### HENYEY-GREENSTEIN PHASE FUNCTION #########################
g=0.9
 
def HGangle():
    rand1=random()
    costheta2=(1/(2*g))*(1+(g**2)-(((1-g**2)/(1-g+2*g*rand1))**2))
    theta1=np.arccos(costheta2)
    phi=2*np.pi*random()
    return(theta1, phi)
    
###################### FUNCTION TO CALCULATE NEW COORDINATE ###################
    
def NewCord(x1,y1,z1,r,theta,phi):
    x1=x1+r*math.sin(theta)*math.cos(phi)
    y1=y1+r*math.sin(theta)*math.sin(phi)
    z1=z1+r*math.cos(theta)
    return (x1,y1,z1)

###################### FUNCTION TO CHECK IF COORD OUTSIDE MATERIAL ############
    
def checkCord(x1,y1,z1,refed,trans,z0):
    if woodcock == True:
        if z1 > z2:
            trans = True
        elif z1 < 0 or x1 > x0 or x1 < 0 or y1 > y0 or y1 < 0:
            refed = True
    else:
         if z1 > z0:
            trans = True
         elif z1 < 0 or x1 > x0 or x1 < 0 or y1 > y0 or y1 < 0:
            refed = True
    return (refed,trans)

###################### FUNCTION TO MAKE A LIST ################################ 
     
def makeList(x_list,y_list,z_list,x,y,z):
    x_list.append(x)
    y_list.append(y)
    z_list.append(z)
    return (x_list,y_list,z_list)
    
##################### FUNCTION TO INTEGRATE ABSORPTION ########################    

def proceed(x1,y1,z1,r,theta,phi):
    (x1,y1,z1)=NewCord(x1,y1,z1,r,theta,phi)

    absed = False
    rand=random()
    if woodcock == True: 
        if 0 < x1 < x0 or 0 < y1 < y0 or 0 < z1 < z0:
            if rand < Pabs1:
                absed = True
        elif 0 < x1 < x2 or 0 < y1 < y2 and  z0 < z1 < z2:  
            if rand < Pabs2:
                absed = True
    else:   
        if rand < Pabs1:
            absed =  True            
    return(x1,y1,z1,absed)
        
##################### FUNCTION TO PLOT 1 PATH #################################
    
def plotPath(x_start,x_end,y_start,y_end,z_start,z_end):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    for i in range (len(x_start)):
       ax.plot([x_start[i],x_end[i]],[y_start[i],y_end[i]],[z_start[i],z_end[i]])
    
##################### FUNCTION TO FIND PATH OF 1 NEUTRON ######################
       
def neutronPath(z0):
    
    trans = False
    refed = False
    absed = False
    x1=x0/2
    y1=y0/2
    z1=0
    x_start=[]
    y_start=[]
    z_start=[]
    x_end=[]
    y_end=[]
    z_end=[]
    thetai=0.0
    phii=0.0
    xsection1 = A
    xsection2 = A2
    if woodcock == True:
        if xsection1 > xsection2:
            maxxsection=xsection1
        elif xsection1 < xsection2: 
            maxxsection=xsection2
     
        if 0 < x1 < x0 or 0 < y1 < y0 or 0 < z1 < z0:
            Fict=xsection1/maxxsection
        elif 0 < x1 < x2 or 0 < y1 < y2 and  z0 < z1 < z2:
            Fict=xsection2/maxxsection
        
        u=random()
        r0=calDistance()
        if Fict < u: 
            extra=mfp2-mfp
            r=r0+extra   
        else:
            r=r0
    else:
        r=calDistance()        
    (x_start,y_start,z_start)=makeList(x_start,y_start,z_start,x1,y1,z1)
    (x1,y1,z1,absed)=proceed(x1,y1,z1,r,thetai,phii)
    (x_end,y_end,z_end)=makeList(x_end,y_end,z_end,x1,y1,z1)
    (refed,trans)=checkCord(x1,y1,z1,refed,trans,z0)
    
    while refed == False and trans == False and absed == False:
       (theta,phi)=HGangle()
       if woodcock == True:
            if xsection1 > xsection2:
                maxxsection=xsection1
            elif xsection1 < xsection2: 
                maxxsection=xsection2
         
            if 0 < x1 < x0 or 0 < y1 < y0 or 0 < z1 < z0:
                Fict=xsection1/maxxsection
            elif 0 < x1 < x2 or 0 < y1 < y2 and  z0 < z1 < z2:
                Fict=xsection2/maxxsection
            
            u=random()
            r0=calDistance()
            if Fict < u: 
                extra=mfp2-mfp
                r=r0+extra   
            else:
                r=r0
       else:
            r=calDistance()        
       (x_start,y_start,z_start)=makeList(x_start,y_start,z_start,x1,y1,z1)
       (x1,y1,z1,absed)=proceed(x1,y1,z1,r,theta,phi)
       (x_end,y_end,z_end)=makeList(x_end,y_end,z_end,x1,y1,z1)
       (refed,trans)=checkCord(x1,y1,z1,refed,trans,z0)
    return (x_start,x_end,y_start,y_end,z_start,z_end,refed,trans,absed)
    
##################### FUNCTION TO FIND PROBABILITIES ##########################
    
def findProb(n,z0):
    refed_counter=0
    trans_counter=0
    absed_counter=0
    for i in range(n):
        (x_start,x_end,y_start,y_end,z_start,z_end,refed,trans,absed)=neutronPath(z0)
        if refed == True:
            refed_counter += 1
        elif trans == True:
            trans_counter += 1
        elif absed == True:
            absed_counter += 1
    refed_prob = float(refed_counter)/n
    trans_prob = float(trans_counter)/n
    absed_prob = float(absed_counter)/n
    trans_percent=trans_prob*100
    transvar=((trans_percent*(n-trans_percent))/n)**0.5
    return(trans_prob)

########## FUNCTION TO PLOT LOG(THICKNESS) AGAINST TRANSMISSION ###############
    
translist=[]
thickness=[]

def Graph(z0):    
    for i in range(200):
        z0=z0+0.2
        z0=(float((z0)))
        thickness.append(z0)
        (trans_problog)=(np.log(findProb(100,z0)))
        #print trans_problog
        translist.append(trans_problog)
    plt.plot(thickness, translist)
    plt.xlabel('Thickness(m)')
    plt.ylabel('log(Transmission Probability)')
    plt.title('Graphite')
#from datetime import datetime
#start_time = datetime.now()
#Graph(0.0) 
#end_time = datetime.now()
#print('Duration: {}'.format(end_time - start_time))   
#print translist
Graph(10)
#(trans_prob,transvar)=findProb(10000,4)
#print trans_prob
#print transvar

######################## PLOTTING 1 NEUTRON PATH ##############################
#
#(x_start,x_end,y_start,y_end,z_start,z_end,refed,trans,absed)=neutronPath(z0)
#plotPath(x_start,x_end,y_start,y_end,z_start,z_end)
#    
#plt.savefig('Grpahite2 3d plot')
#
###################### HISTOGRAM ################################################
#
#x123=[]
#for i in range(10**4):
#    rand = random()
#    r =(-mfp)*math.log(1-rand)
#    x123.append(r)
#
#s = x123 
#
#plt.hist(s, 10**2, normed=True, facecolor='green', alpha=0.5)
#plt.xlabel('Path distance')
#plt.ylabel('Frequency Density')
##plt.title('Histogram of free path distribution')
#plt.savefig('Path proof1')
##plt.show()
#
###################### OLD RANDOM ANGLE FUCNTION ###############################
#
##this is a function that has a uniform distribution
#
#from mpl_toolkits.mplot3d import Axes3D
#x=[]
#y=[]
#z=[]
#n=1000
#for i in range(n):
#    phi = np.random.uniform(0,np.pi*2)
#    costheta = np.random.uniform(-1,1)
#    theta = np.arccos( costheta )
#    x.append(np.sin( theta) * np.cos( phi ))
#    y.append(np.sin( theta) * np.sin( phi ))
#    z.append(np.cos( theta ))
#fig = plt.figure(2)
#ax = fig.add_subplot(111, projection='3d')
#ax.set_aspect("equal")
#ax.scatter(x,y,z)
#
######################### HG FUNCTION TEST #####################################
#
##FUNCTION TO PLOT SCATTERING ANGLES AS SPHERE USING HENYEY-GREENSTEIN
#
#g=0.9
#
#x=[]
#y=[]
#z=[]
#
#thetalist=[]
#
#for i in range(1000):
#    rand1=random()
#    costheta2=(1/(2*g))*(1+(g**2)-(((1-g**2)/(1-g+2*g*rand1))**2))
#    theta1=np.arccos(costheta2)
#    thetalist.append(theta1)
#    phi=2*np.pi*random()
#    x.append(np.sin( theta1) * np.cos( phi ))
#    y.append(np.sin( theta1) * np.sin( phi ))
#    z.append(np.cos( theta1))
#fig = plt.figure(2)
#ax = fig.add_subplot(111, projection='3d')
#ax.set_aspect("equal")
#ax.scatter(x,y,z)
#ax.set_xlabel('X axis')
#ax.set_ylabel('Y axis')
#ax.set_zlabel('Z axis')
#fig = plt.figure(3)
##ax.text2D(0.05, 0.95, "Henyey-Greenstein Distribution", transform=ax.transAxes)
#ax.dist = 11








        
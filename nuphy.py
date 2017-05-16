#!/usr/bin/python3

print('You entered NuPhy.py, merger of calc and srimbuster');
print('---------------------------------------------------');

#http://stackoverflow.com/questions/7427101/dead-simple-argparse-example-wanted-1-argument-3-results
import argparse

import numpy as np
import pandas as pd
import NuPhyPy.db.ReadNubase  as db
import NuPhyPy.Reactions.Kinematics as kin

################# SRIM
#import NuPhyPy.srim as sr
import NuPhyPy.srim.srim as sr

#################### PARSE 
parser=argparse.ArgumentParser(description='NuPhy.py');

parser.add_argument('mode' )
parser.add_argument('-i','--incomming',  default="h2,c12" )
parser.add_argument('-o','--outgoing',  default="h2,c12" )
parser.add_argument('-e','--energy', required=True)
parser.add_argument('-t','--theta',  default="15")
args=parser.parse_args() 

if args.mode=='react':
    print(args.mode,' incomming=',args.incomming, 'outgoing=',args.outgoing)
    ip=args.incomming.split(',')
    op=args.outgoing.split(',')
    nu1=db.isotope( ip[0] )
    nu2=db.isotope( ip[1] )
    nu3=db.isotope( op[0] )
    if len(op)==1:
        print( nu1.A +  nu2.A - nu3.A )
        print( nu1.Z +  nu2.Z - nu3.Z )
        nu4=db.isotope( nu1.A +  nu2.A - nu3.A, nu1.Z +  nu2.Z - nu3.Z )   
    else:
        nu4=db.isotope( op[1] )
    TKE=float( args.energy)
    th=float( args.theta)
    res=kin.react( nu1,nu2,nu3,nu4, TKE=TKE, theta=th,silent=0)

###############################################################  SRIM    
if args.mode=='srim':
    print(args.mode)
    ipath=sr.CheckSrimFiles()
        
    h1=db.isotope('h1')
    Eini=5.800
    rho=2.253
    mgcm2=sr.get_mgcm2( 25.5,  rho )
    print( mgcm2,  sr.get_um( mgcm2, rho))
    ipath=sr.CheckSrimFiles()
    TRIMIN=sr.PrepSrimFile( ion=h1, energy=Eini, angle=0., number=300 ,
                            mater='c', thick=mgcm2, dens=rho  )

    
    tmpp=sr.run_srim(ipath, TRIMIN)

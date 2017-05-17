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

from NuPhyPy.db.ReadNubase import gas,densities,elements
#import NuPhyPy.db.ReadNubase as db
#from NuPhyPy.db.ReadNubase import gas,densities,elements



def material_density( matnam ):
    print('i... testing if ',matnam,'is in elements ')
    if matnam.title() in elements:
        CC=matnam.title()
        print(CC,type(CC))
        zzz=elements.index(  matnam.title() )
        #print(zzz)
        dens=densities[ elements.index( matnam.title() ) ]
        print('i... element ',matnam.title() ,'found, density is set to:', dens)
    else:
        print('i... element NOT found, maybe it is an isotope?')
        isotope=db.isotope( matnam.title() )
        dens=isotope.isodensity
        print('i...  isotope:',isotope.name,'found;  density is set to:',dens)
    return dens




            
#################### PARSE 
parser=argparse.ArgumentParser(description='NuPhy.py');

parser.add_argument('mode' )
parser.add_argument('-i','--incomming',  default="h2,c12" , help='REACT')
parser.add_argument('-o','--outgoing',  default="h2,c12"  , help='REACT')
parser.add_argument('-a','--angle',  default="15" , help='REACT')

parser.add_argument('-e','--energy', required=True , help='REACT/SRIM')

parser.add_argument('-t','--thickness',  default="4" , help='SRIM')
parser.add_argument('-m','--material',   default="c12",help='SRIM')
parser.add_argument('-n','--number',  default="100" , help='SRIM')
parser.add_argument('-d','--density',  default="0" , help='SRIM')

parser.add_argument('-s','--silent', action="store_true",   help='SRIM')

#parser.add_argument('-t','--thickness',  default="4" , help='SRIM')
args=parser.parse_args() 

if args.mode=='react':
    print(args.mode,' incomming=',args.incomming, 'outgoing=',args.outgoing)
    ip=args.incomming.split(',')
    op=args.outgoing.split(',')
    nu1=db.isotope( ip[0] )
    nu2=db.isotope( ip[1] )
    nu3=db.isotope( op[0] )
    if len(op)==1:
        print( 'i...',nu1.A +  nu2.A - nu3.A )
        print('i...', nu1.Z +  nu2.Z - nu3.Z )
        nu4=db.isotope( nu1.A +  nu2.A - nu3.A, nu1.Z +  nu2.Z - nu3.Z )   
    else:
        nu4=db.isotope( op[1] )
    TKE=float( args.energy)
    th=float( args.angle)
    res=kin.react( nu1,nu2,nu3,nu4, TKE=TKE, theta=th,silent=0)

###############################################################  SRIM    
if args.mode=='srim':
    print('i...',args.mode)
    ipath=sr.CheckSrimFiles()
        
    n0=db.isotope(args.incomming )
    material=args.material
    Eini=float( args.energy )
    number=int(args.number)

    ##### MY UNIT WILL BE mg/cm2
    rho=float(args.density)
    if rho==0:
        rho=material_density(material)
    thick=args.thickness
    thickok=False
    if thick.find('ug')>0:
        print('D... ug/cm2')
        thick=float(thick.split('ug')[0])/1000
        thickok=True
    elif thick.find('mg')>0:
        print('D... mg/cm2')
        thick=float( thick.split('mg')[0] )
        thickok=True
    elif thick.find('um')>0:
        print('!... um ! I need rho')
        thick=float(thick.split('um')[0])
        mgcm2=sr.get_mgcm2( thick,  rho ) # um in, mgcm2 out
        thickok=True
        quit()
    if not(thickok):
        print('!...  thicknesses must be in ug,mg or um')
        quit()
    print('i...',material,'thickness',thick,'mg/cm2','for rho=',rho)
    #    print( material_density(material) )
    #print('D...',thick)

    TRIMIN=sr.PrepSrimFile( ion=n0, energy=Eini, angle=0., number=number ,
                            mater=material, thick=thick, dens=rho  )

    # RUN ############################
    if args.silent:
        tmpp=sr.run_srim(ipath, TRIMIN,  silent=True)
    else:
        tmpp=sr.run_srim(ipath, TRIMIN,  silent=False)
    print(tmpp[:5])
    deint=tmpp['e'].max()-tmpp['e'].min()
    sigma=tmpp['e'].std()
    mean=tmpp['e'].mean()
    median=tmpp['e'].median()
    print()
    print( "{:.3f} MeV (median {:.3f}) +- {:.3f}  hi-lo={:.3f}  Eloss={:.3} MeV".format( mean, median,  sigma , deint ,  Eini-mean)   )
    #print( tmpp['e'].max(),tmpp['e'].min(), de  )
    #plt.hist( tmpp['e'], 20 , facecolor='red', alpha=0.25)
#    print("R...    E mean +- std")
#    print(tmpp['e'].mean(), '  ' ,tmpp['e'].std() )
 #   print(tmpp['e'].mean(), '  ' ,tmpp['e'].std() )

    

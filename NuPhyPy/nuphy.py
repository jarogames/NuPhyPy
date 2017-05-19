#!/usr/bin/python3

#print('You entered NuPhy.py, merger of calc and srimbuster');
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

import NuPhyPy.srim.compounds as srcomp

def material_density( matnam ):
    print('i... testing if ',matnam,'is in compounds ')
    #print( 'D...  keys=',srcomp.material_text.keys()  )
    if matnam.title() in srcomp.material_text.keys():
        print('i... FOUND in compounds')
        dens=srcomp.get_density( matnam.title()  )
        print('D... density from compounds=',dens)
    elif matnam.title() in elements: 
        print('i... testing if ',matnam,'is in elements ')
        CC=matnam.title()
        #print(CC,type(CC))
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

parser.add_argument('mode' ,help='react; srim; store ')
parser.add_argument('-i','--incomming',  default="h2,c12" , help='REACT')
parser.add_argument('-o','--outgoing',  default="h2,c12"  , help='REACT')
parser.add_argument('-a','--angle',  default=0 , help='REACT')

parser.add_argument('-e','--energy', default=5.804 , help='REACT+SRIM')

parser.add_argument('-t','--thickness',  default="100ug" , help='SRIM')
parser.add_argument('-m','--material',   default="c12",help='SRIM')
parser.add_argument('-n','--number',  default=100 , help='SRIM')
parser.add_argument('-d','--density',  default=0 , help='SRIM')
parser.add_argument('-P','--Pressure',  default=1013.25e+3 ,type=float, help='SRIM')
parser.add_argument('-T','--Temperature',  default=273.15 , type=float,help='SRIM')

parser.add_argument('-s','--silent', action="store_true",   help='SRIM')

parser.add_argument('-S','--Store', default="",   help='SRIM + STORE')

parser.add_argument('-f','--fwhm', default=0.,   help='STORE - convolute with a detector resolution (FWHM) in [MeV]')

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
    rho=float(args.density)

    if material.title() in srcomp.material_gas:
        """
        trim assumes the target as STP before 1982
        T=273.15 K
        p=1013.25 kPa
        SRIMHelp/helpTR23.rtf
        """
        rho=material_density(material)
        R=1013.25e+3/273.15/rho
        rho2=args.Pressure/R/args.Temperature
        print('i...rho at STD (0deg,1013kPa)=',rho,' now=',rho2)
        rho=rho2

    
    ##### MY UNIT WILL BE mg/cm2
    if rho==0:
        rho=material_density(material)
    thick=args.thickness
    thickok=False
    if thick.find('ug')>0:
        #print('D... ug/cm2')
        thick=float(thick.split('ug')[0])/1000
        thickok=True
    elif thick.find('mg')>0:
        #print('D... mg/cm2')
        thick=float( thick.split('mg')[0] )
        thickok=True
    elif thick.find('um')>0:
        print('!... um ! I need rho=',rho)
        thick=float(thick.split('um')[0])
        mgcm2=sr.get_mgcm2( thick,  rho ) # um in, mgcm2 out
        thickok=True
        #thick=sr.get_um( thick,rho )
        thick=mgcm2
        #quit()
    if not(thickok):
        print('!...  thicknesses must be in ug,mg or um')
        quit()
    print('i...',material,'thickness',thick,'mg/cm2','for rho=',rho,
          '...', 1000*thick/rho/1e-2  ,'A  ',
           1000*thick/rho/1e+2,'um')
    #    print( material_density(material) )
    #print('D...',thick)

    TRIMIN=sr.PrepSrimFile( ion=n0, energy=Eini, angle=0., number=number ,
                            mater=material, thick=thick, dens=rho  )

    # RUN ############################
    if args.silent:
        tmpp=sr.run_srim(ipath, TRIMIN,  silent=True, nmax=number)
    else:
        tmpp=sr.run_srim(ipath, TRIMIN,  silent=False, nmax=number)
    print(tmpp[-5:])
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

    ### MAYBE - I WILL NUMBER ALREADY FROM HERE
    if args.Store!="":
        store = pd.HDFStore( args.Store )
        existing=len(store.keys())
        print('D... already existing is', existing )
        #print(store)
        if material.title() in srcomp.material_gas:
            pt='P{}_T{}'.format( args.Pressure, args.Temperature )
        else:
            pt=""
            
        fname='{:03d}_{}_in_{:_<6s}_w{:_<6s}_r{:3.1f}_{}_n{:04d}_e{:06.3f}_ef{:06.3f}'.format( existing,
                            args.incomming, args.material, args.thickness, float(args.density),
                             pt, int(args.number), float(args.energy), mean )
        fname=fname.replace('.','_')
        store[fname] = tmpp
        print(store)
        store.close()

    

###############################################################  STORE
if args.mode=='store':
    import matplotlib.pyplot as plt
    if args.Store!="":
        stor=args.Store.split(',')
        store = pd.HDFStore( stor[0] )
        if len(stor)>1:
            plots=[]
            for inx in range( len(stor)-1 ):
                dfname=sorted(store.keys())[ int(stor[inx+1])]
                print('o... openning: ', dfname )
                df=store[dfname]
                if float(args.fwhm)>0.: ########## GENERATE GAUSS ############ fwhm=2.355sigma
                    print('i...  mean before convolution: {:.3f} {:.4f}'.format(df['e'].mean(),df['e'].std() ))
                    df['fwhm']=np.random.normal( 0.0, float(args.fwhm)/2.355 ,  len(df) )
                    df['e']=df['e']+df['fwhm']
                    print('i... mean with    convolution: {:.3f} {:.4f}'.format(df['e'].mean(),df['e'].std() ))
                plt.hist( df['e'], 20, ec='k',alpha=0.3,label=dfname) 
            plt.legend( loc=4 , fontsize="x-small" )
            plt.show()
        else:
            for i,v in enumerate(sorted(store.keys())):print("{:03d} {}".format(i,v) )
        
        store.close()
    else:
        print('!... filename missing: use -S; open an item by specifying  line number after comma')

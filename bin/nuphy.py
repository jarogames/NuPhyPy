#!/usr/bin/python3

print('---------------------------------------------------');
print('You entered NuPhy.py, merger of calc-react and srimbuster');
print("""
Three MODES:
  nuphy.py  react  -i h2,c12  -o h1  -e 19  -a 11
  nuphy.py  srim  -i h1 -m si -e 28 -w 5500um -n 100
  nuphy.py  srim  -i he4 -m c  -e 5.804 -w 100ug -n 100
STORE in file:
  nuphy.py  srim  -i he4 -m c  -e 5.804 -w 100ug -n 100 -S ~/srim_h1h2
 LIST:
  nuphy.py  store -S ~/srim_h1h2
 PLOT
  nuphy.py  store -S ~/srim_h1h2,0,1 
  nuphy.py  store -S ~/srim_h1h2,0,1 -f yz
  nuphy.py  store -S ~/srim_h1h2,0,1 -f x
  nuphy.py  store -S ~/srim_h1h2,0,1 -f cos
  nuphy.py  store -S ~/srim_h1h2,0,1 -f cosy
  nuphy.py  store -S ~/srim_h1h2,0,1 -f 0.100
 RRATE:
  nuphy.py  rrate -n 10nA  -t 10ug          -m c   -xs 0.1
  nuphy.py  rrate -n 10uA  -t 10um   -d 1.8 -m c   -xs 0.1
  nuphy.py  rrate -n 1uA   -t 100ug         -m h2  -xs 0.1   
  ...  ... current Nb/t;  rho_S;  Na/Mm ;  100mb=> R
""")
print('---------------------------------------------------');
#time.sleep(1)
#http://stackoverflow.com/questions/7427101/dead-simple-argparse-example-wanted-1-argument-3-results
import argparse

import numpy as np
import pandas as pd
import NuPhyPy.db.ReadNubase  as db
import NuPhyPy.Reactions.Kinematics as kin

import re
################# SRIM
#import NuPhyPy.srim as sr
import NuPhyPy.srim.srim as sr

from NuPhyPy.db.ReadNubase import gas,densities,elements,molarweights
#import NuPhyPy.db.ReadNubase as db
#from NuPhyPy.db.ReadNubase import gas,densities,elements

import NuPhyPy.srim.compounds as srcomp


######################################
# look for compounds/elements/isotopes : return tabulated density
#
######################################
def material_density( matnam ):
    #print('i... testing if ',matnam,'is in compounds ')
    #print( 'D...  keys=',srcomp.material_text.keys()  )
    mm_amu=1.0  # NONSENSE
    if matnam.title() in srcomp.material_text.keys():
        print('F... FOUND in compounds')
        dens=srcomp.get_density( matnam.title()  )
        print('i... FOUND density from compounds=',dens)  # mm_amu ???
        # mm_amu determined probably somewhere in SRIM.... # i dont need
    elif matnam.title() in elements: 
        print('F...  ',matnam,' FOUND in elements ')
        CC=matnam.title()
        #print(CC,type(CC))
        zzz=elements.index(  matnam.title() )
        #print(zzz)
        dens=densities[ elements.index( matnam.title() ) ]
        mm_amu=molarweights[  elements.index( matnam.title() ) ]  # Mm for mix of isotopes
        print('i... element ',matnam.title() ,'found, density is set to:', dens)
    else:
        print('i... element NOT found, maybe it is an isotope?')
        isotope=db.isotope( matnam.title() )
        dens=isotope.isodensity
        mm_amu=isotope.amu   # This is for cross sections  Mm pure
        print('i...  isotope: {} found;  density is set to: {:.4f} g/cm3'.format(isotope.name,dens) )
    return dens,mm_amu




#           
#################### PARSE


parser=argparse.ArgumentParser(description="""NuPhy.py

10um of Si:

./nuphy.py srim  -i he4 -m si -e 5.804 -n 100 -t 10um  

with density:
./nuphy.py srim  -i he4 -m si -e 5.804 -n 100 -t 10um  -d 1.2 

gas:
./nuphy.py srim  -i he4 -m air -e 5.804 -n 100 -t 10000um  -P 150 -T 293 

store:
./nuphy.py store


""");

parser.add_argument('mode' ,help='react; srim; store')
parser.add_argument('-i','--incomming',  default="h2,c12" , help='REACT')
parser.add_argument('-o','--outgoing',  default="h2,c12"  , help='REACT')
parser.add_argument('-a','--angle',  default=0 , help='REACT')

parser.add_argument('-e','--energy', default=5.804 , help='REACT+SRIM')

parser.add_argument('-t','--thickness',  default="100ug" , help='SRIM')
parser.add_argument('-m','--material',   default="c12",help='SRIM')
parser.add_argument('-n','--number',  default=100 , help='SRIM')
parser.add_argument('-d','--density',  default="0" , help='SRIM')
parser.add_argument('-P','--Pressure',  default=1013.25e+3 ,type=float, help='pressure in Pascals e.g.1013250')
parser.add_argument('-T','--Temperature',  default=273.15 , type=float,help='SRIM')

parser.add_argument('-s','--silent', action="store_true",   help='SRIM')

parser.add_argument('-xs','--xsection', default=1.0,   help='XS in barns')

parser.add_argument('-S','--Store', default="",   help='Determine the STORE file; When reading- tell spectrum # -S file,0,1')

parser.add_argument('-f','--fwhm', default=0.,   help='with STORE:  convolute with a detector resolution (FWHM) in [MeV]; BUT ALSO  -f yz ... plot Y vs. Z scatter;  -f cos ... plot radians scatter; -f view ... view')

#parser.add_argument('-t','--thickness',  default="4" , help='SRIM')
args=parser.parse_args() 





############################################
#
#  reaction
#
##########################################
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
#
#  srim    -  first functions -
#
###############################


##################
#  question on GAS = we use 2x
#################
def isitGas( material ):
    if material.title() in srcomp.material_gas:
        return 1
    if material.title() in elements:
        eZ=elements.index(material.title())
        if gas[ eZ ]==1:
            return 1
    try:
        isotope=db.isotope( material , silent=True  )
    except:
        return 0
    if isotope==None:
        eZ=isotope.Z
        if gas[ eZ ]==1:
            return 1
    return 0
        





#####################################
# get_thick_rho  ......  I extract the part
#         return   thick [ug/cm2] and rho g/cm3
#####################################
def get_thick_rho( material,  thick, rho ):   #convert properly um and find rho
    print('D... in  get_thick_rho ::: DUPLICITE CODE:', material, thick, rho, '-------------')
    rho=float(rho)
    if rho==0:
        rho,mm_amu=material_density(material) # compound/element/isot = ALL
    # GASEOUS DENSITY
    # 1/ if compound - rho from function
    # 2/ element:
    rho2=rho # to keep if solid phase
    #========test gas phase =========
    if material.title() in srcomp.material_gas:
        """
        trim assumes the target as STP before 1982
        T=273.15 K
        p=1013.25 kPa
        SRIMHelp/helpTR23.rtf
        """
        R=1013.25e+3/273.15/rho
        rho2=args.Pressure/R/args.Temperature
        print('i...rho at STD (0deg,1013kPa)=',rho,' NEW now=',rho2)
    elif material.title() in elements:
        #print("D... rho - elements")
        eZ=elements.index(material.title())
        print("D... rho - elements  eZ={:d}".format(eZ))
        #print(gas)
        if gas[ eZ ]==1:
            R=1013.25e+3/273.15/rho
            rho2=args.Pressure/R/args.Temperature
            print('i...rho at STD (0deg,1013kPa)=',rho,' NEW now=',rho2)
    else: # could be also a gaseous  isotope
        print('D... maybe gaseous isotope?')
        try:
            isotope=db.isotope( material , silent=True )
        except:
            isotope=None
        if isotope==None:
            print('D... not a gaseous isotope')
        else:
            eZ=isotope.Z
            if gas[ eZ ]==1:
                R=1013.25e+3/273.15/rho
                rho2=args.Pressure/R/args.Temperature
                print('i...rho at STD (0deg,1013kPa)=',rho,' NEW now=',rho2)
    rho=rho2
    print("D... rho WAS deduced")
    # THICKNESS TO mgcm2:    ##### MY UNIT WILL BE mg/cm2
    thickok=False
    if thick.find('ug')>0:
        thick=float(thick.split('ug')[0])/1000
        thickok=True
    elif thick.find('mg')>0:
        thick=float( thick.split('mg')[0] )
        thickok=True
    elif thick.find('um')>0:
        print('!... um ! I need rho=',rho)
        thick=float(thick.split('um')[0])
        thick=sr.get_mgcm2( thick,  rho ) # um in, mgcm2 out
        thickok=True
    elif thick.find('mm')>0:
        print('!... mm ! I need rho=',rho)
        thick=float(thick.split('mm')[0])
        thick=sr.get_mgcm2( thick*1000,  rho ) # um in, mgcm2 out
        thickok=True
    if not(thickok):
        print('!...  thicknesses must be in ug,mg or um,mm')
        quit()
    print('i... {} thickness {:.6f} mg/cm2 for rho={:.3f} ... {:.0f} A = {:.2f}um'.format( material.title(),thick,
                                                rho,1000*thick/rho/1e-2  ,   1000*thick/rho/1e+2 ) )
    return thick, rho, mm_amu
    






######################################
#PREPARE TRIMIN
#  prepare single layer, return one TRIM.IN line
########################################
def prepare_trimin( material,  thick,  rho  ):
    '''
    Here I prepare materials:  single layers
    '''
    print('D... preparing TRIMIN:', material, thick, rho, '-------------')
    #print('D... PV/T density:')
    rho=float(rho)
    if rho==0:
        rho,mm_amu=material_density(material) # compound/element/isot = ALL
    # GASEOUS DENSITY
    # 1/ if compound - rho from function
    # 2/ element:
    rho2=rho # to keep if solid phase
    if material.title() in srcomp.material_gas:
        """
        trim assumes the target as STP before 1982
        T=273.15 K
        p=1013.25 kPa
        SRIMHelp/helpTR23.rtf
        """
        R=1013.25e+3/273.15/rho
        rho2=args.Pressure/R/args.Temperature
        print('i...rho at STD (0deg,1013kPa)=',rho,' NEW now=',rho2)
    elif material.title() in elements:
        print("D... rho - elements")
        eZ=elements.index(material.title())
        print("D... rho - elements  eZ=",eZ)
        if gas[ eZ ]==1:
            R=1013.25e+3/273.15/rho
            rho2=args.Pressure/R/args.Temperature
            print('i...rho at STD (0deg,1013kPa)=',rho,' now=',rho2)
    else: # could be also a gaseous  isotope
        print('D... maybe gaseous isotope?')
        try:
            isotope=db.isotope( material , silent=True )
        except:
            isotope=None
        if isotope==None:
            print('D... not a gaseous isotope')
        else:
            eZ=isotope.Z
            if gas[ eZ ]==1:
                R=1013.25e+3/273.15/rho
                rho2=args.Pressure/R/args.Temperature
                print('i...rho at STD (0deg,1013kPa)=',rho,' now=',rho2)
    rho=rho2
    print("D... rho deduced")
    # THICKNESS TO mgcm2:    ##### MY UNIT WILL BE mg/cm2
    thickok=False
    if thick.find('ug')>0:
        thick=float(thick.split('ug')[0])/1000
        thickok=True
    elif thick.find('mg')>0:
        thick=float( thick.split('mg')[0] )
        thickok=True
    elif thick.find('um')>0:
        print('!... um ! I need rho=',rho)
        thick=float(thick.split('um')[0])
        thick=sr.get_mgcm2( thick,  rho ) # um in, mgcm2 out
        thickok=True
    elif thick.find('mm')>0:
        print('!... mm ! I need rho=',rho)
        thick=float(thick.split('mm')[0])
        thick=sr.get_mgcm2( thick*1000,  rho ) # um in, mgcm2 out
        thickok=True
    if not(thickok):
        print('!...  thicknesses must be in ug,mg or um, mm')
        quit()
    print('i... {} thickness {:.6f} mg/cm2 for rho={:.3f} ... {:.0f} A = {:.2f}um'.format( material.title(),thick,
                                                rho,1000*thick/rho/1e-2  ,   1000*thick/rho/1e+2 ) )

    # AT THIN MOMENT I HAVE A GOOD rho an thick in mgcm2

    #print('D... goto TRIMIN MULTI')
    TRIMIN=sr.PrepSrimFile( ion=n0, energy=Eini, angle=0., number=number,
                            mater=material, thick=thick, dens=rho  )
    
    return TRIMIN
############# END OF PREPARE TRIMIN











###############################
#  SRIM MODE
#        here i try to prepare also multilayer  materials
#        
##############################

if args.mode=='srim':
    print('i...',args.mode)
    ipath=sr.CheckSrimFiles()
        
    n0=db.isotope(args.incomming )
    Eini=float( args.energy )
    number=int(args.number)
    rho=args.density
    thick=args.thickness
    material=args.material.title()  # this will get complicated with layers
    nmats=len(material.split(','))
    print("D... counting number of layers",nmats)
    if nmats>1:
        print('!... ',nmats,'materials - TEST REGIME:')
        if nmats!=len(thick.split(',')):
            print('!... NOT THE SAME NUMBER OF THICKNESSES')
            quit()
        if nmats!=len(rho.split(',')):
            if float(rho)!=0.:
                print('!... NOT THE SAME NUMBER OF densities')
                quit()
            else:
                rho=','.join( map(str,[0]*nmats) )
        print('i... PREPARING ANALYSIS  mat, thick, rho:')
        TRILIST=[]
        for imat in range(nmats):
            print(imat,'... =========', material.split(',')[imat].title(), thick.split(',')[imat],  rho.split(',')[imat],
                  "================================"  )
            TRILIST.append( prepare_trimin(  material.split(',')[imat].title(), thick.split(',')[imat],  rho.split(',')[imat]   )   )
        print('i... I GOT ALL TRIM.IN files. Now somebody must merge....')
        ############################   MERGING  LAYERS ##############
        print('D--------------')
        for xi,xii in enumerate(TRILIST):  print( xi, xii,'\n')          # PRINT
        TRIMIN="==> SRIM-2013.00 This file controls TRIM Calculations.\r\n"
        #TRIMIN= TRIMIN + TRILIST[0].split("\n")[0].rstrip()+'\r\n'
        TRIMIN= TRIMIN + TRILIST[0].split("\n")[1].rstrip()+'\r\n'
        TRIMIN= TRIMIN + TRILIST[0].split("\n")[2].rstrip()+'\r\n'
        TRIMIN= TRIMIN + TRILIST[0].split("\n")[3].rstrip()+'\r\n'
        TRIMIN= TRIMIN + TRILIST[0].split("\n")[4].rstrip()+'\r\n'
        TRIMIN= TRIMIN + TRILIST[0].split("\n")[5].rstrip()+'\r\n'
        TRIMIN= TRIMIN + TRILIST[0].split("\n")[6].rstrip()+'\r\n'
        TRIMIN= TRIMIN + TRILIST[0].split("\n")[7].rstrip()+'\r\n'   # target material  num ele and layers
        #TRIMIN=TRIMIN+"            8 \r\n"

        
        #print(TRIMIN)
        #quit()
        #############################################################
#        with open('TRIM.IN.ALL','w') as f:
#            for imat in range(nmats):
#                f.write(TRILIST[imat])
        line8=[]      # Target material+1
        TRILISTTOT=[] # totallist for Atom
        lineLay=[]    # Layer Layer long line
        layerList=[]
        for imat in range(nmats):
            listOfLines=TRILIST[imat].split('\n')
            if listOfLines[7].find('Target material')<0:
                print('!... Target material line not detected...quit')
                quit()
            # define line8 (7: Target material)
            line8.append(   re.findall(r'".+"|[\w\.]+',  listOfLines[8] )  )
            TRILISTTOT.extend(  listOfLines  )
            # lineLay...
            liLa=[ i for i in listOfLines if re.match('^Layer\s+Layer\s+Name',i)  ]
            lineLay.extend(liLa)  # line with columns for layers.
            #---- duplicate: but i find #line of Layer Layer +2
            for j,v in enumerate( listOfLines ): 
                if v.find('Layer')>=0 and v.find('Density')>0:
                    print(imat,'LAYER LINE= #',j+2+1) # starts with 0
                    layerList.append( re.findall(r'".+"|[\w\.\-]+', listOfLines[j+2] ) )
        #print(line8)
        layname='...'.join( [ i[0].strip("\"").rstrip() for i in line8]  )
        layname='"'+layname+'"'
        nelems=sum( map(int, [i[1] for i in line8] ) )
        nlayers=sum( map(int, [i[2] for i in line8] ) )
        print(layname,nelems,nlayers)
        TRIMIN=TRIMIN+"{} {} {}\r\n".format(layname,nelems,nlayers)
        #TRIMIN= TRIMIN + TRILIST[0].split("\n")[8].rstrip()+'\r\n'  # He 4  into C  - nedavat
        TRIMIN= TRIMIN + TRILIST[0].split("\n")[9].rstrip()+'\r\n'
        TRIMIN= TRIMIN + TRILIST[0].split("\n")[10].rstrip()+'\r\n' 
        TRIMIN= TRIMIN + TRILIST[0].split("\n")[11].rstrip()+'\r\n' 
        
        # NOW I need to get all Atom .. = ... = lines
        #print(TRILISTTOT)
        #regex = re.compile('^Atom\s\d+\s=\s')
        atoms=[i for i in TRILISTTOT if re.match('^Atom\s\d+\s=\s',i) ]
        if len(atoms)!=nelems:
            print('!... Atoms lines and # elements differ!')
            quit()
        #        atoms=[m.group(1) for l in TRILISTTOT for m in [regex.search(l)] if m]
        #atoms=re.findall(r'Atom\s\d+\s=\s', TRILISTTOT)
        for i,v in enumerate(atoms):
            atoms[i]=re.sub( 'Atom\s(\d+)\s','Atom '+str(i+1)+' ' , v ).rstrip()
        print( '\n'.join( atoms ) )
        TRIMIN=TRIMIN+ '\r\n'.join(atoms)+'\r\n'
        # NOW I need Layer Layer
        preLayer=re.sub(r'^(.+Density).+$',r'\1', lineLay[0] )
        for i,v in enumerate(lineLay):
            lineLay[i]=re.sub('^.+?Density ','', v).strip()
        print( preLayer,'  '.join(lineLay) )
        TRIMIN=TRIMIN+preLayer+' '+'  '.join(lineLay)+'\r\n'
        # NOW I need the same with next "stoich stoich stoich...."
        lineStoich="Numb.   Description                (Ang) (g/cm3)   "+" Stoich "*nelems
        print(lineStoich)
        TRIMIN=TRIMIN+lineStoich+'\r\n'
        # NOW I need 1   "Layer 1";  2  "Layer 2"
        zeroes=0
        for i,v in enumerate(layerList):
            zstring=' 0   '*zeroes
            zpost=' 0   '*(nelems-zeroes -  int(line8[i][1]) )
            prepart='  '.join(v[1:4])
            postpart='  '.join(v[4:])
            print( "{} {} {} {} {}".format(i+1,prepart,zstring,postpart, zpost) )
            TRIMIN=TRIMIN+ "{} {} {} {} {}\r\n".format(i+1,prepart,zstring,postpart, zpost) 
            zeroes=zeroes+int(line8[i][1])
            # 16 fields;
        # NOW I need to copy lines with GAS.....
        lineGas="0  Target layer phases (0=Solid, 1=Gas)"
        print(lineGas)
        TRIMIN=TRIMIN+lineGas+'\r\n'
        for i,v in enumerate(layerList):
            ###print( 'material====',v )
            print(  isitGas(v[1]),' '   , end=' ')
            TRIMIN=TRIMIN+ str(isitGas(v[1]) )+' '
        print()
        TRIMIN=TRIMIN+'\r\n'
        ######## RAGG
        lineBragg="Target Compound Corrections (Bragg)"
        print(lineBragg)  # i want to put all lines  12 ( i think)
        TRIMIN=TRIMIN+lineBragg+'\r\n'
        for imat in range(nmats):
            listOfLines=TRILIST[imat].split('\n')
            for jmat in range( len(listOfLines) ):
                if listOfLines[jmat].find('Target Compound Corrections')>=0:
                    #print( listOfLines[jmat].rstrip() )
                    print( listOfLines[jmat+1].rstrip()  , end=' ')
                    TRIMIN=TRIMIN+ listOfLines[jmat+1].rstrip() +' '
                    break
            if jmat==len(listOfLines)-1:
                print('!... Bragg line not found')
                quit()
        print()
        TRIMIN=TRIMIN+'\r\n'
        #Individual target atom displacement energies (eV)
        lineAtomDisp="Individual target atom displacement energies (eV)"
        print(lineAtomDisp)
        TRIMIN=TRIMIN+lineAtomDisp+'\r\n'
        for imat in range(nmats):
            listOfLines=TRILIST[imat].split('\n')
            for jmat in range( len(listOfLines) ):
                if listOfLines[jmat].find('Individual target atom displacement energies')>=0:
                    #print( listOfLines[jmat].rstrip() )
                    print( listOfLines[jmat+1].rstrip()  , end=' ')
                    TRIMIN=TRIMIN+ listOfLines[jmat+1].rstrip() +' '
                    break
            if jmat==len(listOfLines)-1:
                print('!... Displacement line not found')
                quit()
        print()
        TRIMIN=TRIMIN+'\r\n'

        #Individual target atom displacement energies (eV)
        lineAtomDisp="Individual target lattice binding energies (eV)"
        print(lineAtomDisp)
        TRIMIN=TRIMIN+lineAtomDisp+'\r\n'
        for imat in range(nmats):
            listOfLines=TRILIST[imat].split('\n')
            for jmat in range( len(listOfLines) ):
                if listOfLines[jmat].find('Individual target atom lattice binding energies')>=0:
                    #print( listOfLines[jmat].rstrip() )
                    print( listOfLines[jmat+1].rstrip()  , end=' ')
                    TRIMIN=TRIMIN+ listOfLines[jmat+1].rstrip() +' '
                    break
            if jmat==len(listOfLines)-1:
                print('!... Lattice binding line not found')
                quit()
        print()
        TRIMIN=TRIMIN+'\r\n'
        #Individual target atom displacement energies (eV)
        lineAtomDisp="Individual target atom surface binding energies (eV)"
        print(lineAtomDisp)
        TRIMIN=TRIMIN+lineAtomDisp+'\r\n'
        for imat in range(nmats):
            listOfLines=TRILIST[imat].split('\n')
            for jmat in range( len(listOfLines) ):
                if listOfLines[jmat].find('Individual target atom surface binding energies')>=0:
                    #print( listOfLines[jmat].rstrip() )
                    print( listOfLines[jmat+1].rstrip()  , end=' ')
                    TRIMIN=TRIMIN+ listOfLines[jmat+1].rstrip() +' '
                    break
            if jmat==len(listOfLines)-1:
                print('!... Surface binding line not found')
                quit()
        print()
        TRIMIN=TRIMIN+'\r\n'
        print('Stopping Power Version (1=2011, 0=2011)')
        TRIMIN=TRIMIN+'Stopping Power Version (1=2011, 0=2011)\r\n'
        print(' 0 ')
        TRIMIN=TRIMIN+' 0\r\n'
        
        print('\n\n\n',TRIMIN,"\n\n\n")
        #quit()

    else:
        print("D... one material - preparing TRIMIN")
        TRIMIN=prepare_trimin(  material, thick,  rho   ) 
    #########################################
    # PREPARE FILE
    ##########################################
#    print('D... goto TRIMIN')
#    TRIMIN=sr.PrepSrimFile( ion=n0, energy=Eini, angle=0., number=number,
#                            mater=material, thick=thick, dens=rho  )

    # RUN ############################
    if args.silent:
        tmpp=sr.run_srim(ipath, TRIMIN,  silent=True, nmax=number)
    else:
        tmpp=sr.run_srim(ipath, TRIMIN,  silent=False, nmax=number)
    print(tmpp[-5:])
    if 'e' in tmpp:
        deint=tmpp['e'].max()-tmpp['e'].min()
        sigma=tmpp['e'].std()
        mean=tmpp['e'].mean()
        median=tmpp['e'].median()
        print()
        print( "{:.3f} MeV (median {:.3f}) +- {:.3f}  hi-lo={:.3f}  Eloss={:.3} MeV".format( mean, median,  sigma , deint ,  Eini-mean)   )
    else:
        #print(tmpp)
        print()
        print('{:.3f} +- {:.4f} um implanted depth'.format( tmpp['x'].mean(), tmpp['x'].std() )  )
        mean=0.0
        median=0.0
        deint=0.0
        sigma=0.0
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
            
        fname='n{:03d}_{}_in_{:_<6s}_w{:_<6s}_r{:3.1f}_{}_n{:04d}_e{:06.3f}_ef{:06.3f}'.format( existing,
                            args.incomming, args.material, args.thickness, float(args.density),
                             pt, int(args.number), float(args.energy), mean )
        fname=fname.replace('.','_')
        store[fname] = tmpp
        print(store)
        store.close()

# RUN SRIM ENDED






        
###############################################################  
#
#   Reaction Rate         
#   rate=Nb/t; Nt/A ; sigma==cross section==xs
#####################################
if args.mode=='rrate':
    Nb=args.number
    if Nb[-2:]=='nA':   Nb=float(Nb[:-2])*1e-9/1.602176e-19
    elif Nb[-2:]=='uA': Nb=float(Nb[:-2])*1e-6/1.602176e-19
    elif Nb[-2:]=='mA': Nb=float(Nb[:-2])*1e-3/1.602176e-19
    else: Nb=float(Nb)

    sigma=  float(args.xsection)  # give it in barns= cm2

    rho_S,just_rho,Mm_amu=get_thick_rho( args.material , args.thickness, args.density )
    rho_S = rho_S  /1e+3 # mg -> g/cm2
    Na=6.24+23   # i think - Avogadro at/ mol   or 1e+26/kmol
    Mm=Mm_amu    # molar weight -  isotopic should be
    
    rrate=Nb * rho_S *Na /Mm *sigma
    print("Beam {:.4g} cps; XS={:.4g} barns; rho_S={:.4g} mg/cm2  Mm={:.3f} g/mol".format( Nb, sigma, 1000*rho_S, Mm) )
    print("Reaction rate: {:.4g} cps".format(rrate))
    print("Reaction rate: {:.4g} cps/sr".format(rrate/4/3.1415926 ))





    

    
###############################################################  STORE
#
#   (UN) store  and   plot
#
#####################################
if args.mode=='store':
    import matplotlib.pyplot as plt
    if args.Store!="":
        stor=args.Store.split(',')
        store = pd.HDFStore( stor[0] )
        if len(stor)>1:
            plots=[]
            #========= i want to plot dE histogram with convolution fwhm
            for inx in range( len(stor)-1 ):
                dfname=sorted(store.keys())[ int(stor[inx+1])]
                print('o... openning: ', dfname )
                df=store[dfname]
                plotme=False
                #==================  if fwhm .... CONVOLUTE;  yz ... scatter;  cosy cosz ...angles
                if args.fwhm=="x": #======= implantation plot
                    plt.hist( df['x'], 20, ec='k',alpha=0.3,label=dfname+' implant [um]')
                    plotme=True
                if args.fwhm=="y": #======= implantation plot
                    plt.hist( df['y'], 20, ec='k',alpha=0.3,label=dfname+' implant [um]')
                    plotme=True
                if args.fwhm=="z": #======= implantation plot
                    plt.hist( df['z'], 20, ec='k',alpha=0.3,label=dfname+' implant [um]')
                    plotme=True
                                        
                    #plt.scatter( df['z']/1e+4, df['y']/1e+4,alpha=0.3,label=dfname+' [um]' ) 
                if args.fwhm=="yz": #======= scatter plot  Y Z 
                    plt.scatter( df['z'], df['y'],alpha=0.3,label=dfname+' [um]' )
                    plotme=True
                if args.fwhm=="cos" and  'e' in df.keys(): #======= scatter plot cosy cosz
                    plt.scatter( df['cosz'], df['cosy'], alpha=0.3, label=dfname+' [rad]' )
                    plotme=True

                if args.fwhm=="cosy" and  'e' in df.keys(): #======= scatter plot cosy cosz
                    plt.hist( df['cosy'], 20, ec='k',alpha=0.3,label=dfname+' implant [rad]')
                    plotme=True
                if args.fwhm=="cosz" and  'e' in df.keys(): #======= scatter plot cosy cosz
                    plt.hist( df['cosz'], 20, ec='k',alpha=0.3,label=dfname+' implant [rad]')
                    plotme=True
                if args.fwhm=="view": #======= VIEW DATA
                    print(df)
                try:
                    floatfwh=float(args.fwhm)
                except:
                    floatfwh=-1
                if floatfwh>=0.: #====== GENERATE GAUSS ======== fwhm=2.355sigma
                    if 'e' in df.keys():
                        print("...... MEAN ENERGY MODE")
                        print('i...  mean before convolution: {:.3f} {:.4f}'.format(df['e'].mean(),df['e'].std() ))
                        df['fwhm']=np.random.normal( 0.0, float(args.fwhm)/2.355 ,  len(df) )
                        df['e']=df['e']+df['fwhm']
                        print('i...  mean with   convolution: {:.3f} {:.4f}'.format(df['e'].mean(),df['e'].std() ))
                        plt.hist( df['e'], 20, ec='k',alpha=0.3,label=dfname)
                        plotme=True
                    elif floatfwh>0:
                        print("...... MEAN ENERGY MODE FOR IMPLANT")
                        df['fwhm']=np.random.normal( 0.0, float(args.fwhm)/2.355 ,  len(df) )
                        ENE=dfname.split('_')[-4]+'.'+dfname.split('_')[-3]
                        ENE=float(ENE[1:])
                        df['e']=ENE+df['fwhm']
                        print('i...  mean with   convolution: {:.3f} {:.4f}'.format(df['e'].mean(),df['e'].std() ))
                        plt.hist( df['e'], 20, ec='k',alpha=0.3,label=dfname)
                        plotme=True

                    else:
                        print("....... IMPLANT MODE")
                        plt.hist( df['x'], 20, ec='k',alpha=0.3,label=dfname+' [um]')
                        plotme=True
                    
            if plotme:
                plt.legend( loc=4 , fontsize="x-small" )
                plt.show()
        else:
            for i,v in enumerate(sorted(store.keys())):print("{:03d} {}".format(i,v) )
        
        store.close()
    else:
        print('!... filename missing: use -S; open an item by specifying  line number after comma')

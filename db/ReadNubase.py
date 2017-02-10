#!/home/ojr/anaconda3/bin/python
import sys # stdout.write
import re
import numpy as np
from time import sleep
import os

def isfloat(value):
        try:
                float(value)
                return True
        except ValueError:
                return False

############ GLOBAL STUFF:  MASSTABLE DATA AND LIST
masstable=''
masslist=[]
massnp=np.zeros( shape=(300,150)  )
massnp=massnp-999999.  # invalid value



elements=['n','H','He','Li','Be','B','C','N','O','F','Ne',
	 		'Na','Mg','Al','Si','P','S','Cl','Ar',
 			'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',  'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',
 			'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe',
 			'Cs', 'Ba',
 			'La','Ce','Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu','Hf',
 			'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi',
 			'Po', 'At', 'Rn', 'Fr','Ra', 'Ac',  'Th', 'Pa', 'U',
 			'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 
 			'Fm', 'Md', 'No', 'Lr','Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 
 			'Ds', 'Rg', 'Cn', '113', 'Fl', '115', 'Lv', '117', '118'];


class isotope:
        massline=''   # MY LINE
        mex=0.0
        dmex=0.0
        spin=0
        halflife=None
        IS=0.0
        name=''
        A=0
        Z=0
        N=0
        spin=9999.
        parity=0.0
        amu=0.0
        def __init__(self, A,Z):
                global masstable
                global masslist
                global massnp
                global elements
####################################### LOAD MASSTABLE FIRST
#                print( os.path.dirname( os.path.abspath(__file__) ) )
                databasepath= os.path.dirname( os.path.abspath(__file__) )
                if ( len(masstable)<=0):
                        with open(databasepath+'/nubase.mas12') as f:
                                masstable=f.read()
                                masslist=masstable.split('\n')
                        print('masses loaded ... ', len(masslist),'lines')
                        for li in masslist:
                                if (len(li)>0):
                                        #print( li,  int(li[0:3]), int(li[4:7]) ,  )
                                        if ( li[7] is '0' ):
                                                flo=li[17:24]+'.'+li[25:29]
                                                #print(flo, int(li[0:3]), int(li[4:7])  , li)
                                                if isfloat( flo ):
                                                        massnp[ int(li[0:3]), int(li[4:7]) ] = float(  flo )
                                        #print( int(li[0:3]), int(li[4:7]),    )
                                        #sleep(0.2)
                        # i want a map ... numpy [a,z]
######################## MASSTABLE LOADED ######## READ Z A 
                if (type(Z) is str):
                        sys.stdout.write("element "+Z+' has Z=')
                        try:
                                Z=elements.index(Z)
                                print(Z)
                        except:
                                print('\n!... NO SUCH ELEMENT ',Z)
                                return None
                else:
                        print('Z=',Z,' is element',elements[Z])
                if (A>=massnp.shape[0] or Z>=massnp.shape[1]):
                        print('No isotope',A,Z)
                        return None
                #print('alloc ...',A,Z)
                Anu='%03d' % A
                Znu='%03d0' % Z
                rese='^'+Anu+' '+Znu+' '
                #print( '... going to search ',rese )
                if ( massnp[A,Z]>-999999.):
                        for li in masslist:
                                if ( re.search( rese , li) ):
                                        self.massline=li.rstrip()
                                        #print(li)
                                        #break
                                        if (len(self.massline)>0):
                                                self.name=self.massline[10:16].strip()
                                                self.mex=massnp[A,Z]
                                                self.dmex=   float( self.massline[30:33]+'.'+self.massline[34:38] )
                                                halv=self.massline[61:64]+'.'+self.massline[65:68]
                                                halv=halv.replace('<', '')
                                                halv=halv.replace('>', '')
                                                halv=halv.replace('stb.', 'None')  # precisely put 0. for stables
                                                #print( halv )
                                                if ( isfloat(halv) ):
                                                        self.halflife=float( halv )
                                                        #print(  self.halflife )
                                                else:
                                                        self.halflife= None
                                                        #print(  'NOT FLOAT ',self.halflive )

                                                uniti=self.massline[69:71]
                                                if (uniti.strip()  == 'Zy'):
                                                        self.halflife=self.halflife*365*24*3600*1e+24
                                                if (uniti.strip()  == 'Zy'):
                                                        self.halflife=self.halflife*365*24*3600*1e+21
                                                if (uniti.strip()  == 'Ey'):
                                                        self.halflife=self.halflife*365*24*3600*1e+18
                                                if (uniti.strip()  == 'Py'):
                                                        self.halflife=self.halflife*365*24*3600*1e+15
                                                if (uniti.strip()  == 'Ty'):
                                                        self.halflife=self.halflife*365*24*3600*1e+12
                                                if (uniti.strip()  == 'Gy'):
                                                        self.halflife=self.halflife*365*24*3600*1e+9
                                                if (uniti.strip()  == 'My'):
                                                        self.halflife=self.halflife*365*24*3600*1e+6
                                                if (uniti.strip()  == 'ky'):
                                                        self.halflife=self.halflife*365*24*3600*1e+3
                                                if (uniti.strip()  == 'y'):
                                                        self.halflife=self.halflife*365*24*3600
                                                if (uniti.strip()  == 'd'):
                                                        self.halflife=self.halflife*24*3600
                                                if (uniti.strip()  == 'h'):
                                                        self.halflife=self.halflife*3600
                                                if (uniti.strip()  == 'm'):
                                                        self.halflife=self.halflife*60
                                                if (uniti.strip()  == 'ms'):
                                                        self.halflife=self.halflife*1e-3
                                                if (uniti.strip()  == 'us'):
                                                        self.halflife=self.halflife*1e-6
                                                if (uniti.strip()  == 'ns'):
                                                        self.halflife=self.halflife*1e-9
                                                if (uniti.strip()  == 'ps'):
                                                        self.halflife=self.halflife*1e-12
                                                if (uniti.strip()  == 'fs'):
                                                        self.halflife=self.halflife*1e-15
                                                if (uniti.strip()  == 'as'):
                                                        self.halflife=self.halflife*1e-18
                                                if (uniti.strip()  == 'zs'):
                                                        self.halflife=self.halflife*1e-21
                                                if (uniti.strip()  == 'ys'):
                                                        self.halflife=self.halflife*1e-24
                                                #print( 'mex %f dmex %f' % (self.mex,  self.dmex ) )

                                                IS=self.massline[110:].split(' ')[0]
                                                if (self.halflife==None):
                                                        self.IS=float( IS[3:])
                                                else:
                                                        self.IS=0.0
                                                self.A=A
                                                self.Z=Z
                                                self.N=A-Z
                                                spin=self.massline[79:93]
                                                self.spinfull=spin
                                                parity=None
                                                if (spin.find('+')>=0):
                                                        parity=+1
                                                if (spin.find('-')>=0):
                                                        parity=-1
                                                spin=spin.replace('+','')
                                                spin=spin.replace('-','')
                                                
                                                spin=spin.replace('(','')
                                                spin=spin.replace(')','')
                                                spin=spin.replace('#','')
                                                if (spin.split() ): spin=spin.split()[0]
                                                if (spin.split(',') ): spin=spin.split(',')[0]
                                                #spin=spin.split(',')[0]
                                                #print( 'SPIN==', spin)
                                                spinh=re.findall( r'(\d+)/2', spin)
                                                if ( len(spinh)>0): 
                                                        #print('spinh:',spinh)
                                                        spin=int(spinh[0])/2.
                                                if (isfloat(spin)):
                                                        self.spin=float(spin)
                                                else:
                                                        self.spin=-9999.
                                                self.parity=parity  
                                                self.amu= (self.mex/1000.  + self.A * 931.49403)/931.49403;
                                                #print( self.name , self.mex, self.dmex, ' ...  ', self.halflife , '?'+uniti.strip()+'?' , IS, self.IS )
                                                break

                else:
                        return None
#########################################################################
###########     ISOTOPE CREATED
#
#
#########################################################################
#        def mex(self):
#                return self.mex
#        def dmex(self):
#                return self.dmex
#        def halflife(self):
#                return self.halflife
#        def IS(self):
#                return self.IS
#        def name(self):
#                return self.name




###########################
#             MAIN 
####################################################################
if __name__ == "__main__":
    print("ReadNuBase.py is being run directly...")
    n14=isotope( 14,'N' )
    n15=isotope( 15,7 )
    o15=isotope( 15,8 )

    print( 'AMU',n14.name,n14.amu )
    print( 'AMU',n15.name,n15.amu )
#    print( 'n14 half, spin, parity', n14.halflife , n14.spin , n14.parity , n14.IS )
#    print( 'n15 half, spin, parity', n15.halflife , n15.spin , n15.parity , n15.IS )
    print( 'o15 half, spin, parity, mex, abund', o15.halflife , o15.spin , o15.parity , o15.mex, o15.IS)
#for a in range(260):
#        for z in range(200):
#                he3=isotope(a,z)
                ###print( he3.mex   )
                ###print( he3.mex, he3.dmex, he3.halflife )
                ###if (he3.name):print( he3.name, he3.mex, he3.dmex, he3.halflife, he3.IS , he3.A, he3.Z, he3.N , he3.spinfull, he3.spin, '*')
#                if (he3.name):print( he3.name, he3.spinfull, he3.spin, he3.parity, '*')
                ###print( he3.mex, he3.dmex, he3.halflife, he3.IS()  )
else:
    print("ReadNubase is being imported into another module")

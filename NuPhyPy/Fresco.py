#!/usr/bin/env python3

#!/home/ojr/anaconda3/bin/python
import os # search directory
import re 
import numpy as np
import matplotlib.pyplot as plt
####
import subprocess  ## to run fresco with pip out
#
#https://t2.lanl.gov/nis/data/endf/    decay and neutronxs data
#PyNE======================================
#http://pyne.io/examples/endf_reader.html
#http://pyne.io/examples/half_life.html
#================= masstable .......pip install masstable
#
#https://pypi.python.org/pypi/masstable/0.2.1
#
def isfloat(value):
    try:
        float(value)
        return True
    except ValueError:
        return False

############################## i want to call import NuPhyPy.Fresco as fr
########################                       o16fr=fr.newtask(  'coulomb.inp' , silent=True )

class newtask:
    '''
    SECTION contains sections  &name ..... /
    some values must ne inegers  .. integs
    '''
    SECTION={}
    npartitions=0
    npartitionsta=0
    pot=0
#    couplings=0
    LABEL=''
    integs=['iter', 'treneg', 'chans', 'smats', 'xstabl', 'zp', 'zt', 
        'nex', 'cpot', 'ptyp', 'ptyt', 'bandt', 'bandp', 'kp', 'type', 'shape', 'kn1',
        'ic1', 'ic2','in', 'kind', 'nn', 'l', 'kbpot', 'isc', 'ipc',
        'in', 'ib','ia','kn',
        'icfrom', 'icto', 'kind', 'ip1', 'ip2', 'ip3' ]
#    SECTION['FRESCO']={'hcm':0.05, 'rmatch':30}  #dict  of dicts

    lastinput=''

    def __init__(self,filename, silent=True):      ######## READ CONFIG TO SECTION
        '''
        MAIN PROC
        init: opens .inp file
              extracts a LABEL
              detects NAMELIST
              removes comments / blablabla
              and anaylses the SECTION  &start .... /  in identsection
                  after he has SECTION keys
        '''
        self.SECTION={} # 
        #frescopath= os.path.dirname( os.path.abspath(__file__) )
        #print("FREPATH", frescopath)
        self.lastinput=filename
        #print( self.SECTION['FRESCO']['hcm'])
        with open( filename ,'r') as out:
            line1=out.readline().rstrip()
            self.LABEL=line1
            print( 'FRESCO /',filename,'/ ... LABEL: ',line1)
            if ( out.readline().rstrip() != 'NAMELIST' ):
                print('problem with line2 - not NAMELIST')
                exit(1)
            section=''    
            for line in out:
                line=line.rstrip()+' '#re.sub( '\n' , ' ', line) ADD SPACE to compensate \n
                section=section+line
                if (section.count('/')>0 ):
                    section=re.sub( '/.+' , '', section)   # clear all after the slash
                    self.identsection( section )
                    section=''
            if not silent:
                print('--------GETTING FILE TO MEMORY----------------------------------------')
                for i in sorted(self.SECTION.keys()):
                    print(i, self.SECTION[i] )
                print('--------OK FILE IN MEMORY---------------------------------------------')




#######################################################
#
#
#  ident section ---- called from   constructor
#
#
#################################################
    def identsection(self, section ):
        '''
        2nd to MAIN
        This procedure takes one section and compares with
              FRESCO,  PARTITION  STATES  POT  Overlap   Coupling   CFP
        With the input it can call:
               filldictFRE
        or
               filldict
        the problems were several:
               read parameter lists, 
               elab(1)
               POT kp=i   complex sections
        '''
        sobfresco=re.search( r'&FRESCO [\w\d=\s\.\-\()\:]+',  section , re.IGNORECASE )
        if (sobfresco):
        #    print( '======FRESCOG:\n',sobfresco.group(0) )
#            self.filldict(  sobfresco.group(0) ,'FRESCO')
            self.filldictFRE(  sobfresco.group(0) )
#            print( self.SECTION)
        
        sobfresco=re.search( r'&PARTITION [\w\d=\'\"\w\s\.\-\()\:]+',  section , re.IGNORECASE )
        if (sobfresco):
            #print( 'PARTIIT:',sobfresco.group(0) )
            self.npartitions+=1
            self.filldict(  sobfresco.group(0) ,'PARTITION'+str(self.npartitions) )
            

        sobfresco=re.search( r'&STATES [\w\d=\'\"\w\s\.\-\()\:]+',  section , re.IGNORECASE )  ### INSIDE PARTITION
        if (sobfresco):
            #print( 'STATIS:',sobfresco.group(0) )
            self.filldict(  sobfresco.group(0)  ,'PARTITION'+str(self.npartitions) , 'STATES')
            

        sobfresco=re.search( r'&POT [\w\d=\'\"\w\s\.\-\()\:]+',  section , re.IGNORECASE   )  ### 
        if (sobfresco):
            kp=re.findall( 'kp=\s*(\d+)',  section ) 
            typep=re.findall( 'type=\s*(\d+)',  section ) #...... and by type (0=Coulomb)
            if (len(typep)==0):
                typep=[0]
            #print( 'POT:'+str( kp[0] ) ,sobfresco.group(0) )
            self.filldict(  sobfresco.group(0)  ,'POT'+str( kp[0] ) +'_'+str(typep[0]) )

        sobfresco=re.search( r'&Overlap [\w\d=\'\"\w\s\.\-\()\:]+',  section , re.IGNORECASE )  ### KN1
        if (sobfresco):
            kp=re.findall( 'kn1=\s*(\d+)',  section ) 
            #print( 'Overlap:'+str( kp[0] ) ,sobfresco.group(0) )
            self.filldict(  sobfresco.group(0)  ,'Overlap'+str( kp[0] ) )
            
        sobfresco=re.search( r'&Coupling [\w\d=\'\"\w\s\.\-\()\:]+',  section , re.IGNORECASE )  ### 
        if (sobfresco):
            #print( 'Coupling:' ,sobfresco.group(0) )
            #self.couplings+=1
            self.filldict(  sobfresco.group(0)  ,'Coupling' )

        sobfresco=re.search( r'&CFP [\w\d=\'\"\w\s\.\-\()\:]+',  section , re.IGNORECASE )  ### CFP should corresp.2 Overlap KN=KN1
        if (sobfresco):
            kn=re.findall( 'kn=\s*(\d+)',  section )      #by kn
            #print( 'CFP:',sobfresco.group(0) )
            self.filldict(  sobfresco.group(0)  ,'CFP'+str(kn[0])   )
            


            


###########################################################
#
#         FRESCO doesnt have subsetions or complications.
#
#
############################################

    def filldictFRE( self , group  ):
        #print('=====', group ,'=====\n')
        if (  not 'FRESCO' in  self.SECTION ):
            self.SECTION['FRESCO']={}
        # new strategy:   
        #       1/ remove all \s  before  =
        #       2/ rm all \s  inside '  '
        #       3/ rm all \s  inside ()
        group = re.sub( r'\s+=\s*', '=', group )
        group = re.sub( r'\=\s+', '=', group )
        group = re.sub( r'(\(\S*)\s+(\S*\))', r'\1\2', group )
        group = re.sub( r'(\(\S*)\s+(\S*\))', r'\1\2', group )
        group = re.sub( r'(\'.*)\s+(.*\')', r'\1\2', group )

        #print('=====', group ,'=====\n')
        items=[]
        # if : do one thing if not : do other
        itran=re.findall( r'(\w[\w\d\(\)\:]+)=',  group )   # not for elab(1)
        if ( len(itran)>0):
            #print( 'RANGE () ',itran )
            dic1={}
            dic2=[]
            for i in itran:  # i have identified te NAMES, now i get prams for each
                j=i.replace( '(', '\(')
                j=j.replace( ')', '\)')
                j=j.replace( ':', '\:')
                j=j+'='
                g=re.search( j+'[\d\.\-\s]+' , group )
                dic2.append( g.start() )
                dic1[i]=len( dic2 )-1
                #print( j,  g.start()   )  
            dic2.append( len(group))
            for i in dic1.keys():
                a=dic2[ dic1[i] ]+len(i)+1 # without name and=
                b=dic2[ dic1[i]+1 ] 
                #print('/',i, '/', group[a:b] )
                self.SECTION['FRESCO'][ i ]=group[a:b]
                if isfloat(group[a:b]):
                    self.SECTION['FRESCO'][ i ]=float(group[a:b])
#                print('/',i, dic1[i],'.',  a, '-',  b , group[a:b])

##            print( 'itran3', itran[0][3] )
            #items1=re.findall( r'([\d\.\-\+]+)',   itran[0][3]  )   # numbers......
            #ii=0
            #for i in range( int( itran[0][1] ),  int( itran[0][2] )+1  ):
                #ip123=itran[0][0]+str( i )
                ##print('i',i,  ip123 ,'=', items1[ii] )
                #self.SECTION['FRESCO'][ ip123 ]= float( items1[ii] )
                #ii=ii+1




                
###################################################################
#
# every section tries to decode its parameters;
# 1/ POT  -    i need to breakdown  p1 p2 p3 p4 p5 p6 - WORKS
#
# 2/ FRESCO   elab  nlab -  i need to keep elab(1:3)  nlab(1:2)
######################################



    def filldict( self , group , secname, subsec=''  ):
        #print('=====', group )
        #print('=====')
        # i added possibility a=...0.0
        #   and  elab(1)
        if (  not secname in  self.SECTION ):
            self.SECTION[secname]={}
        else:
            print('SEC', secname, 'exists:', self.SECTION[secname] )
        items=[]
        itemsbis=[]
        # if : do one thing if not : do other
        itran=re.findall( r'([\w\d]+)\((\d)\:(\d)\)=\s*([\-\.\+\d\s]+)',  group )   # not for elab(1)
        if ( len(itran)>0):
            #print( 'RANGE () ',itran )
#            print( 'itran3', itran[0][3] )
            items1=re.findall( r'([\d\.\-\+]+)',   itran[0][3]  )   # numbers......
            ii=0
            for i in range( int( itran[0][1] ),  int( itran[0][2] )+1  ):
                ip123=itran[0][0]+str( i )
                #print('i',i,  ip123 ,'=', items1[ii] )
                self.SECTION[secname][ ip123 ]= float( items1[ii] )
                ii=ii+1

#            print( 'items1', items1 )
            #print( self.SECTION )
        items=re.findall( r'([\w\d\()\:]+)=\s*([\-\.\+\d\'\w]+)',  group )   #  for elab(1)
        #items=items + itemsbis
        
        #print(' items = ',items )
        if (subsec==''):
            for i in items:
#                print('Q',i)
                if (  not(':' in i[0] ) ):
#                    print('P',i)
                    i1=i[1]
                    #i1=re.sub( "\'" , '', i1)    # remove 'Proton' --- Proton
                    #print(i,'->',i1)
                    if ( isfloat( i1 ) ):        # isfloat means it is a number
                        if (  i[0] in self.integs ):
                            self.SECTION[secname][i[0]]=int( i1 )
                        else:
                            self.SECTION[secname][i[0]]=float( i1 )
                    else:
                        self.SECTION[secname][i[0]]=i1
        if (subsec!=''):
            if (  not subsec in  self.SECTION ):
                self.SECTION[secname][subsec]={}
            for i in items:
#                print('Q',i)
                if (  not(':' in i[0] ) ):
#                    print('P',i)
                    i1=i[1]
                    #i1=re.sub( "\'" , '', i1)    # remove 'Proton' --- Proton
                    #print(i,'->',i1)
                    if ( isfloat( i1 ) ):        # isfloat means it is a number
                        if (  i[0] in self.integs ):
                            self.SECTION[secname][subsec][i[0]]=int( i1 )
                        else:
                            self.SECTION[secname][subsec][i[0]]=float( i1 )
                    else:
                        self.SECTION[secname][subsec][i[0]]=i1
        return 0





    

#############################################################################
#
#
#   helpers for  WRITE
#
#
########################################################


    def decomtuples(self, dicti, nopat='xxxxxxxx'):
        '''
        Used in WRITE
        '''
        #print('DECOM TUPLE',dicti)
        wrsec=''
        states={}
        for i in sorted( dicti.keys() ):
            p123=re.search( nopat, i )
            if ( p123 ):
                print('XXXXXXXX pattern found', i, nopat)
                continue
            #if isinstance( dicti[i] , dict ):
            if (i == 'STATES'):
                #print('i see instance tuple')
                states=dicti[i]
            else:
                #print('i see instance NON tuple')
                if ( i in self.integs):
                    wrsec=wrsec+' '+i+'='+str(int(dicti[i]))+' '
                else:
                    wrsec=wrsec+' '+i+'='+str(dicti[i])+' '
        if ( len(states)>0 ):
                wrsec=wrsec+ '/\n &STATES '+self.decomtuples( states )
        return wrsec    



    def writesection(self, section , end=1 ):
        wrsec=''
        if ( section in self.SECTION):
            rsection=re.sub( '[\d_]*$', '',section )   # remove pot1, cfp2 ...
            wrsec=wrsec+' &'+rsection+ self.decomtuples( self.SECTION[section]  )
            if (end==1):
                wrsec=wrsec + '\n &'+rsection+' / \n'
            else:
                wrsec=wrsec + ' /\n'
            return wrsec
        else:
            return ''



    def writeinput(self, filename ,silent=True ):
        self.lastinput=filename
        wrfi=self.LABEL+'\n'
        wrfi=wrfi+'NAMELIST\n'
        wrfi=wrfi+self.writesection('FRESCO')
        for i in range(9):
            wrfi=wrfi+self.writesection('PARTITION'+str(i) ,0 )
        wrfi=wrfi+' &PARTITION /\n'
        for i in range(9):
            for j in range(9):
                wrfi=wrfi+self.writesection('POT'+str(i)+'_'+str(j), 0) # but i need to strip 1_1
        wrfi=wrfi+' &POT /\n'
        for i in range(9):
            wrfi=wrfi+self.writesection('Overlap'+str(i),0)
        wrfi=wrfi+ ' &Overlap /\n'
        wrfi=wrfi+self.writesection('Coupling', 0)
        for i in range(9):
            wrfi=wrfi+self.writesection('CFP'+str(i), 0)
        wrfi=wrfi+' &CFP /\n'
        wrfi=wrfi+' &Coupling /\n'
        realPath = os.path.realpath(__file__)
        if not silent:
            print( 'i... going to write to', filename,'at',  os.getcwd()  )
        with open(  filename , 'w') as fout:
            fout.write(wrfi);
        if not silent:
            print('--------WRITTEN FILE FROM MEMORY--------------------------------------')
            for i in sorted(self.SECTION.keys()):
                print(i, self.SECTION[i] )
            print('--------OK -----------------------------------------------------------')
          



            
####################################################
#
#
#              RUN
#
###########################################################

    def run(self,  inputfile='' , silent=True ):
        frescopath= os.path.dirname( os.path.abspath(__file__) )
#        print("FREPATH=", frescopath, 'current is ',os.getcwd())
        fname=inputfile
        if ( len(fname)==0):
            fname=self.lastinput
        if not silent:
            print('i... trying to run with',fname)
        self.lastinput=fname
        ######################## BINARIES ##
        runpath=frescopath+'/fresco_binaries/'
        runname="fresco-winxp.exe"
        runname="runfresco"
        runok=0
        for file in os.listdir( runpath ):
            if (file==runname):
                if not silent:
                    print('+...',runname,'found at ',runpath)
                runok=1
        if runok==0:
            print("!...",runname,'not found at',runpath)
            return -1
#        for file in os.listdir(os.environ['HOME']+"/bin"):
#            if (file == runname ):
#                runpath=os.environ['HOME']+"/bin"
#        for file in os.listdir(os.environ['HOME']+"/bin"):
#            if (file == runname ):
#                runpath=os.environ['HOME']+"/bin"
#            if len(runpath)<3:
#            print("!... sorry, runfresco  NOT FOUND at",runpath)
#            return -1
#        print( '+...',runname ,'found at {',runpath,'}')
        ##################### search for INPUT ########
        inputok=0
        for file in os.listdir(os.getcwd()):
            if (file==fname):
                if not silent:
                    print('+...inputfile=',fname,'found at ',os.getcwd())
                inputok=1
        if inputok==0:
            print("!... sorry, inputfile=",fname,"not found at ",os.getcwd())
            return -1
        CMD=runpath+'/'+runname+' '+ fname
        if not silent:
            print('i... running ',CMD)
        pf=subprocess.Popen(CMD.split(),stdout=subprocess.PIPE )
        while pf.poll() is None:
            l = pf.stdout.readline() # This blocks until it receives a newline.
            if not silent:
                print('  ',l.rstrip().decode('utf-8') )
            # When the subprocess terminates there might be unconsumed output 
            # that still needs to be processed.
        if not silent:
            print( pf.stdout.read().rstrip().decode('utf-8') )
            print("X... exiting fresco .run")
#        if (os.system(CMD) ):
#            print('!... Error with',CMD)







###################################################
#
#
#             READ  OUTFILE
#
#     keyword is name of projectile e.g. in partition
#
#      1/check number of energies
##
#####################################################


    def readout(self, keyword='----------',  filename=''  ):
        '''
        I read  -  input.out  file by default
        the outgoing projectile of the partition
            and only the line given by PARTITION1 namep; projectile name
            with string: CROSS SECTIONS FOR OUTGOING
        It is ok for Rutherford.
        '''
        print('extract XS PARTITION', keyword)
        fname=filename
        if (len(fname)==0):
            fname=self.lastinput  # CRAZY_GETS .ipn file
            fname=fname.replace( '.inp', '.out')
            print('i... trying to blindly read assumed output file: /', fname,'/')
        else:
            print('i... reading file: /', fname,'/')
        ang=[]
        val=[]
        ang2=[]
        val2=[]
        energies=[]
        with open( fname ) as f:  # I will parse by PARTITION
            for i in f:
                if ('ENERGY' in  i)and('=' in i)and('MeV' in i):
                    g=re.findall( 'ENERGY\s*=\s*([\d\.]+)\s+MeV', i)
                    if ( len(g)>0 ):
                        #print(i.rstrip() , g )
                        if isfloat(g[0]):
                            energies.append(  float(g[0])  )
        f.close()
        print('i... detected energies during read:', energies,'MeV')
        with open( fname ) as f:  # I will parse by PARTITION
            for i in f:
                #print('kuku:', i.rstrip()  )
                if ((  'CROSS SECTIONS FOR OUTGOING' in i  ) and (keyword.strip("'") in i ) ):  # FIND IT
                    #print('JUK :', i.rstrip()  )
                    f.readline()
        #            f.readline()
                    break
            for i in f:
                if (( 'Finished' in i ) or ( len(i.strip() )==0 ) or ( 'CROSS SECTIONS FOR OUTGOING 'in i )):  ## FINISH
                    break
                rg=re.findall( r'([\d\.]+) deg.: X-S =\s+([\d\.\+\-E]+) ', i.rstrip() )
                if (len(rg)>0):
                    ang.append(rg[0][0] ) # i want one occurence =[0] with two values=in tuple
                    val.append(rg[0][1] ) # i want one occurence =[0] with two values=in tuple
                else:
                    rg=re.findall( r'/R\s+=\s+([\d\.\+\-E]+)\s*', i )
                    if (len(rg)>0):
                        #ang2.append( ang[-1] ) # i want one occurence =[0] with two values=in tuple
                        val2.append( rg[0] ) # i want one occurence =[0] with two values=in tuple
                        #print( rg[0] )
#    #################### SECOND ROUND #####################
#            ang2=[]
#            val2=[]
#            #for i in f:
#                #print('1st',i)    
#            #if     ( 'CROSS SECTIONS FOR OUTGOING 'in i ):
#                    #break
#            for i in f:
#                if ( 'Finished' in i ):
#                    break
#                #print('2nd',i.rstrip() )    
#                rg=re.findall( r'([\d\.]+) deg.: X-S =\s+([\d\.\+\-E]+) ', i.rstrip() )
#                if (len(rg)>0):
#                    ang2.append(rg[0][0] ) # i want one occurence =[0] with two values=in tuple
#                    val2.append(rg[0][1] ) # i want one occurence =[0] with two values=in tuple
    
    
        npang=np.array( ang )
        npval=np.array( val )
        npval2=np.array( val2 )

        npang=npang.astype('float') 
        npval=npval.astype('float') 
        npval2=npval2.astype('float') 

        return    npang,npval,npval2
#        npang2=np.array( ang2 )
#    
#        plt.plot( npang, npval )
#        plt.plot( npang2, npval2  )
#    
#        plt.gca().set_xlim(5., 50.)
#        plt.xlim(5,50)# OK
#        plt.ylim(0.01 , 50000)# OK
#        #plt.semilogy() #OK
#        #plt.set_xlim(2.,50)
#        plt.grid(True)
#        plt.yscale('log')
#        #plt.axis('tight')
#        plt.show()







###################################################
#
#
#             READ  DATA
#
#####################################################

    def dataread(self,  filename ):
        '''
         nothing?   angle value dvalue
        '''
        an=[]
        val=[]
        dval=[]
        with open( filename ) as f:
            for i in f:
                sss=i.strip().split()
                an.append( sss[0] )
                val.append( sss[1] )
                dval.append( sss[2] )
                #print(i.rstrip() )
        npan=np.array( an ).astype('float')
        npval=np.array( val ).astype('float')
        npdval=np.array( dval ).astype('float')
        return npan, npval,npdval


#fr=frespar( 'n14dp_jm.inp' )
#fr.writedict()
# runfresco 
# check errors...
# plot basics ... matplotlib
# combine with nastro + readout of nudat??? ... automatic partitions with ground states.

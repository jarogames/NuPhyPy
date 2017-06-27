#https://www.nist.gov/programs-projects/ansiieee-n4242-standard
import xml.etree.ElementTree as etree    
import numpy as  np
#import matplotlib.pyplot as plt
#import mpld3
import re

###################### DECODE  REAL TIME AND LIVE TIME FROM STRING
def decode_RT(datum):
        #print(datum)
#        y=re.search(r'P(\d\d)Y(\d\d)M(\d\d)DT(\d\d)H(\d\d)M([\d\.])S', datum)
        y=re.search(r'P(\d\d)Y(\d\d)M(\d\d)DT(\d\d)H(\d\d)M(.*)S', datum)
        #print( y.group(1),  y.group(2), y.group(3),  y.group(4),  y.group(5) ,  y.group(6) )
        return int(y.group(1))*365*86400+int(y.group(2))*30*86400+int(y.group(3))*86400+int(y.group(4))*3600+int(y.group(5))*60+float(y.group(6))


###################### FIND LAST REASONABLE BIN, SKIP SATURATION EVENTS peak @ end
import sys
def last_real_bin(y, limit=10, artfwid=10):
    artf=0
    for i in range(len(y)-1,-1,-1):
        if (y[i]<limit and artf<=0): artf=-i
        if (y[i]>limit and artf<=0): artf=i;print(y[i],artf)
        if (y[i]>limit and i>artf-artfwid):
            #print('artifact found at'+str(i))
            sys.stdout.write(str(i)+' ')

        elif (y[i]>limit):
            print('end bin found at',i)
            return i
    return i


######################################### SEARCH FOR PEAKS - method from TOP
# 1. split spectrum to 100 vertical belts from 0 to MAXY  and 
# 1a.  create a histogram:  what part of spe i over this Y value
#      x: Yvalue=limit;   y: part of the spectrum over the limit
# 2. this spectrum (from hi to low!) rises with steps (when Compton is found)
#    I SEARCH FOR DERIVATIONS:
#           rise in derivations means that some compton appeared
# 3. plot of derivations (looking from high to low) has high values at low
#    I Search xhigh to xlow (which is Yvalue) and
#       at   1/20 limit of maximum I define the threshold (x=Yvalues)
#      step  =  -1% of max
def get_integder(y):
    mx=int(max(y))
    nonul=[x for x in y if x>0]
    intg=[]
    intgd=[]
    intgx=[]
    step=-int(mx/100) ### DIVIDING to 100 steps, reasobably fast
    if step==0: step=-1
    for i in range(mx,-1,step):
        num=0
        tot=0
        for j in range(len(y)):
            tot=tot+1
            #print( y[j] )
            if y[j]>i:
                 num=num+1
        #print(i, num)
        if (num>0):
                intg.append(num/len(nonul))
                if len(intg)>1: intgd.append( intg[-1]-intg[-2] )
                intgx.append(i)

    #SEARCH FOR y THRESHOLD    
    maxd100=max(intgd)/20  # this is sensitive to compton...
    print('searching derivation',maxd100,'out of ',max(intgd))
    for xh in range(len(intgd)):  
        if intgd[xh]>maxd100:  
            trsh=intgx[xh]    #intgx contains Y values 
            #print(' ',xh, '.bin in deriv x yvalue;  current derivation',intgd[xh], 'x val',intgx[xh],'trsh',trsh)
            break
    #return intg,intgd,intgx
#intg,intgd,intgx=get_integder(y23)                
#fig, ax = plt.subplots(figsize=(12, 5))
#plt.subplot(211)
##ax.set_xscale("log")
#ax.set_yscale("log")
#plt.plot( intgx, intg )
#plt.subplot(212)
#plt.plot( intgx[:-1], intgd )
#ax.set_yscale("log")
#plt.show()
    #print('compare to mx',mx,'derivaton rises at ', trsh)
    return trsh


#######################################FIND PEAKS
#
#  two lists are appended with just maximum and x
#  Spectra should be ZEROES and PEAKS only
#             symplified
#
def find_peaks(y,  setx,  sety):
    #print('finding',len(y))
    maxy=0
    xmamx=0
    for i in range(len(y)):
        if (y[i]==0 and maxy>0):
            setx.append( xmax )
            sety.append( maxy )
            maxy=0
            xmax=0
            #print( setx[-1], sety[-1] )
        if (y[i]>maxy):
            maxy=y[i]
            xmax=i


###########################################################
#
#
#    search backwards with findpeaks by chunk, automatic labi
#
#
# i look at derivations to find a streep rise. Cmp to max rise
def iterate_chunks_backwards(y, spelen,labi=0):
    spelen=1000
    peaksx=[]
    peaksy=[]
    if labi==0:
        labi=last_real_bin( y ,10, 10 )
    #fig2, ax2 = plt.subplots(figsize=(12, 5))
    #ax2.set_yscale("log", nonposy='clip')
    #plt.plot(x,y,  drawstyle='steps')
    for i in range(labi-spelen,-1,-spelen):
        trsh=get_integder(y[i:labi])  # trick - i go until end of spectrum
        #print(i,'-',labi,'derivaton rises at ', trsh)
        y2=[q if q>trsh else 0 for q in y]
        y2=[ y2[q] if q>i else 0 for q in range(len(y2))]
        y2=[ y2[q] if q<i+spelen else 0 for q in range(len(y2))]
        #plt.plot(x,y2, 'r' )
        ###plt.plot(x,y2,  drawstyle='steps')
        find_peaks(y2, peaksx,peaksy)
        #ax.set_yscale("log")
###    peakse=[x[o] for o in peaksx]    
    #print(  list(zip(peakse,peaksy) ) )
    #plt.show()
    return peaksx,peaksy


    
################################### get_better_peaks
#
#    2 lists of peaks:  x,y
#       log-lin fit. what is above the line I accept and return
#
import numpy as np
from math import log
import datetime as dt

def get_better_peaks( peakse, peaksy ):
        peaksyl=[ log(q) if q>0 else 1 for q in peaksy]
        reg=np.polyfit( peakse, peaksyl,1)
        
        pkx =[peakse[i] for i in range( len(peakse) )       if peaksyl[i]>=peakse[i]*reg[0] + reg[1] ]
        pky =[peaksy[i] for i in range( len(peakse) )       if peaksyl[i]>=peakse[i]*reg[0] + reg[1] ]
#        pkyl =[log(peaksy[i]) for i in range( len(peakse) ) if peaksyl[i]>=peakse[i]*reg[0] + reg[1] ]
#        fig3, ax3 = plt.subplots(figsize=(12, 5))
#        r_x, r_y = zip(*((i, i*reg[0] + reg[1]) for i in range( int(x[labi])   )))
#        plt.plot(r_x, r_y, "r-")

#        plt.plot(peakse,peaksyl, '.')
#        ax3.set_yscale("log")
#        plt.plot(pkx,pkyl, 'r+')
#        plt.show()
        return pkx,pky








######################################### READ n42 SPECTRUM #################
def read(name, emin=10,emax=3000, verbose=True ):
        print('reading',name)
#        mpld3.enable_notebook()
        tree = etree.parse(name)  
        root = tree.getroot()                    
        #root
        spes=[]
        cv=[]
        for c in root:
            #print(c.tag)
            for d in c:
                if ( c.tag.find('EnergyCalibration')>0 and d.tag.find('CoefficientValues')>0 ):
                    if verbose: print('===============  CV',d.text)
                    cv=[float(x) for x in d.text.split()]
                    if verbose:print(cv)
                if ( c.tag.find('RadMeasurement')>0 and d.tag.find('RealTimeDuration')>0 ):
                    if verbose: print('===============  RT',d.text)
                    rt=decode_RT( d.text )
                    duration=dt.timedelta(seconds=float(rt))
                    if verbose: print(rt)
                if ( c.tag.find('RadMeasurement')>0 and d.tag.find('StartDateTime')>0 ):
                    if verbose: print('===============  SD',d.text)
                    #t = dt.datetime()
                    startdate=d.text.split('T')[0].split('-')
                    starttime=d.text.split('T')[1].split(':')
                    startdate=list(map( int, startdate ))
                    starttime=list(map( int, starttime ))
                    #SD 2017-6-23T14:32:20  not RT P00Y00M00DT00H01M52.68S 
                    #dt=decode_RT( d.text )
                    start=dt.datetime( *startdate, *starttime )
                    if verbose: print(  start )
                if ( c.tag.find('RadMeasurement')>0 and d.tag.find('Spectrum')>0 ):
                    for e in d:
                        if ( e.tag.find('LiveTimeDuration')>0 ):
                            if verbose: print('===============  LT',e.text)
                            lt=decode_RT( e.text )
                            if verbose: print(lt)
                        if ( e.tag.find('ChannelData')>0 ):
                            if verbose: print('===============  LEN=',len(e.text))
                            spes=(e.text).replace('\n',' ').split()
                #print('  ',d.tag,'  ',d.text)
#        print( type(spes) , len(spes))
        spes=[float(i) for i in spes]
        ss=np.array( spes )
        CPS=ss.sum()/duration.total_seconds() 
        if verbose: print( 'Spectrum SUM={:.0f}   CPS={:.2f}'.format( ss.sum() , ss.sum()/duration.total_seconds() ) )
        x=np.linspace(0, cv[0]+cv[1]*len(spes), len(spes))
        xli=[ int(1000*q)/1000 for q in x]
        x=np.array(xli)
        deadt=1.0-lt/rt
        if verbose: print('deadtime = {:.2f}%'.format(100*deadt) )
        return (x,ss, deadt,start,duration, CPS)

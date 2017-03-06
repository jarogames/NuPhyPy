#!/usr/bin/env python3
#from sympy import *
from math import sin,cos,tan,pi

def closest(iz=2, iZ2=6, iTKE=1, theta=180):
#    TKE,z,Z2=symbols('TKE z Z2')
    dpar=197.3 * iz*iZ2/137/ iTKE
    print('closest approach',dpar,'fm')
#    dpar=bb.subs(z,2).subs(Z2,6).subs(TKE,1)
    #print(iz,iZ2,iTKE) # 2 79 1 ... 227.52 fm
    if (theta!=180):
        dpar=dpar/2/tan(theta/180*3.1415926/2)
        rmin=dpar*cos( theta/180*3.1415926/2)/(1-sin(theta/180*3.1415926/2))
        print('closest approach at',theta,'deg.:',rmin,'fm')
        return rmin
    return(dpar)

def impact_par(iz=2, iZ2=6, iTKE=1, theta=45):
    dpar=197.3 * iz*iZ2/137/ iTKE
    b=dpar/2/tan(theta/180*3.1415926/2)
    print('impact parameter=',b,'fm')
    return b

def xs_above(z=2, Z2=6, TKE=1, theta=45):
    theta2=theta/180*3.1415926
    sig=z*z/4*pi*Z2*Z2*(197.3/137/TKE)**2*( (1+cos(theta2))/(1-cos(theta2)))
    sig=sig/100
    print('xs           for angles > ',theta,'deg.:',sig,'barns')
    return sig


def xs_dsdo(z=2, Z2=79, TKE=1, theta=45, silent=0):
    theta2=theta/180*3.1415926
    dpar=197.3 * z*Z2/137/ TKE
    dsdo=dpar**2/16/sin(theta2/2)**4
    dsdo=dsdo/100
    if silent==0:
        print('differential xs at {:7.3f}'.format(theta),'deg.      : {:13.3f}'.format(dsdo),'barns/srad')
    return dsdo

from scipy import integrate as nintegrate
from math import sin,cos,tan,pi

def dsdo_integr(theta, z, Z2, TKE ): # theta in RAD just here
    res=xs_dsdo(z,Z2,TKE,theta/pi*180, silent=1)  # make silent from here
    return res*sin(theta)*2*pi

def xs_above_num(z=2, Z2=6, TKE=1, theta=45):
    result = nintegrate.quad( dsdo_integr, theta/180*pi,pi, args=(z,Z2,TKE) )
    # dsdo is already in barns
    xsin=result[0]
    print('numerical xs for angles > ',theta,'deg.:',xsin,'barns')
    return xsin

def get_dOmega(dist=183, dx=1, dy=3):
    sang= dx*dy/4/pi/dist/dist
    print('collimator {}x{} at {}     dOmega={:.3e} srad'.format(dx,dy,dist,sang) )
    return sang

if __name__ == "__main__":
    theta=45
    xs_dsdo(z=2, Z2=79, TKE=1, theta= theta)
    xs_above(2,79,1,     theta)
    xs_above_num(2,79,1, theta)
    get_dOmega()
######    dsdo_integr(theta, z, Z2, TKE ):


"""
Cosmological functions
Author: Arta Khosravi
Started: 
Last update: Jan. 2024
"""

import numpy as np
import scipy as sp
from scipy.integrate import quad
import scipy.integrate as integrate
import math
from constants_Acc import *


def Mdot(Lbol): 
    return (Lbol/((epsilon_const)*(c**2))) #https://www.aanda.org/articles/aa/pdf/2020/10/aa36649-19.pdf
def Mdotch(Lbol,eps): 
    return (Lbol/((eps)*(c**2)))
def BHARpos(SFR):
    return ((10**(-3.72+0.51))*((SFR)**(1.05+0.33)))*(6.30276*1e+22)
def BHARneg(SFR):
    return ((10**(-3.72-0.51))*((SFR)**(1.05-0.33)))*(6.30276*1e+22)
def BHARLX(Lx):
    return ((1-ep)*kbol*(Lx))/(ep*(c**2)) #https://uhra.herts.ac.uk/bitstream/handle/2299/21348/JS_AAM_1.pdf?sequence=1&isAllowed=y
def epsilonpos(SFR,Lbol):
    return (Lbol)/((c**2)*(((10**(-3.72+0.51))*((SFR)**(1.05+0.33)))*(6.30276*1e+22)))
def epsilonneg(SFR,Lbol):
    return (Lbol)/((c**2)*(((10**(-3.72-0.51))*((SFR)**(1.05-0.33)))*(6.30276*1e+22)))
def epX(Lx,BHAR):
    return (kbol*(Lx))/(BHAR*(c**2))
def epBH(Lbol,BHAR):
    return (Lbol/((BHAR)*(c**2)))

####-----Eddington Ratio-----####
def Ledd(M):
    return 1.26*1e+44*(M/1e+6)
def eddration(M,Lbol):
    return Lbol/(1.26*1e+44*(M/1e+6))
def Leddar(Lbol,M):
    return Lbol/(1.3e+38*M)
def eddratio(M,Lbol):
    return 1.26e+44*(M/1e+6)/Lbol
def MBHM(Lopt,MdotL):
    return 2.6*(1e+8/MdotL)*((Lopt/1e+45)**(1.5))*((cosi)**-1.5)

####-----Accretion Disk efficiency-----####
def Epsilon5100(Lbol,M,L5100): #eq6
    return 0.105*(Lbol/1e+46)*(M/1e+8)*((L5100/1e+45)**(-1.5))*((cosi)**1.5)
def EpsilonL(Lbol,M,L_Lambda,Lambda): #https://arxiv.org/pdf/2309.12944.pdf
    return  0.073*(Lbol/1e+46)*(M/1e+8)*((L_Lambda/1e+45)**(-1.5))*((Lambda/5100)**(-2))*((cosi)**1.5)
def Epsilonopt(Lbol,M,L_Lambda,Lambda): #https://academic.oup.com/mnras/article/419/3/2529/1069495?login=false
    return  0.063*((Lbol/1e+46)**0.99)*((M/1e+8)**0.89)*((L_Lambda/1e+45)**(-1.5))
def Epsilonopt_new(Lbol,M): #from the fit
    return 0.063*((Lbol/1e+46)**0.99)*((M/1e+8)**0.89)*(((0.990495699129281*Lbol-0.137138386847888) /1e+45)**(-1.5))
 
def Mcrit(mbh,epsilon):#https://arxiv.org/pdf/2201.05300.pdf
    return 1.4e+18*1e-3*mbh*((epsilon)**-1/10)
def Mdotog(Lbol,epsilon):
    return Lbol/(epsilon*c**2)
def mdot(Lbol,mbh,epsilon):
    return Mdotog(Lbol,epsilon)/Mcrit(mbh,epsilon)
def Mdotacc(Mbh,Lopt):
    return 3.5*Msun*(yr**(-1))*((Mbh*1e-8)**(-0.89))*((Lopt*1e-45)**(1.5))


def MBHD11(Lopt,Mdot):
    return 2.6*np.power((10**Lopt)*1e-45/(0.8),3/2)*1e+8/(10**Mdot)
def EpsilonKopt(Lbol,M,Lopt): #https://academic.oup.com/mnras/article/419/3/2529/1069495?login=false
    return  0.063*((Lbol/Lopt)**0.99)*((M/1e+8)**0.89)*((Lopt/1e+45)**(-0.51))

def tau(t):
    return t/tedd
def m(MH):
    return (MH)/MH0
def eta(MH):
    return eta0*m(MH)*np.exp(-m(MH)/mcr)
def L_new(MH,Ledd):
    return (1/(alpha))*(ep0)*(Ledd)*np.log(1+alpha*eta(MH))
def epsilon_n(MH):
    return (ep0/alpha*eta(MH))*np.log(1+alpha*eta(MH))
def dm(dtau,MH):
    return dtau*m(MH)*(alpha*eta(MH)-ep0*np.log(1+alpha*eta(MH)))
def Llimlow(MH,Ledd):
    return (ep0)*(Ledd)*eta(MH)
def Llimhigh(MH,Ledd):
    return (1/(alpha))*(ep0)*(Ledd)*np.log(alpha*eta(MH))

 
def Mdotopt(M,Lopt):
    return 3.5*(6.30276*1e+22)*((M/1e+8)**(-0.89))*((Lopt/1e+45)**1.5)
def Epsilon(Lbol,M,Lopt): #eq6
    return 0.063*((Lbol/1e+46)**(0.99))*((M/1e+8)**(0.89))*((Lopt/1e+45)**(-1.5))
def Epsiloneq4(Lbol,M,Lopt): #eq4
    return Lbol/(Mdotopt(M,Lopt)*(9e+16))
#https://academic.oup.com/mnras/article/419/3/2529/1069495?view=extract



def RISCOn(Lbol,M,L5100):
    return G*M/((c**2)*2*(0.105*(Lbol/1e+46)*(M/1e+8)*((L5100/1e+45)**(-1.5))*((cosi)**1.5)))
def Z_1(a):
    return 1+(((1-a**2)**1/3)*((((1+a)**1/3))+(((1-a)**1/3))))
def Z_2(a):
    return (3*(a**2)+((Z_1(a),2)))**(1/2)
def R_Isco_plus(a):
    return 3+Z_2(a)+(((3-Z_1(a))*(3+Z_1(a)+2*Z_2(a)))**(1/2))
def R_Isco_neg(a):
    return 3+Z_2(a)-(((3-Z_1(a))*(3+Z_1(a)+2*Z_2(a)))**(1/2))
def Epsilona_plus(a):
    return 1-(((R_Isco_plus(a))**3/2) -(2*((R_Isco_plus(a))**1/2))+a)/((((R_Isco_plus(a))**3/4))*((((((R_Isco_plus(a))**3/2))-3*(((R_Isco_plus(a))**1/2))+2*(a))**1/2)))
def Epsilona_neg(a):
    return 1-((np.power(R_Isco_neg(a),3/2) -(2*np.power(R_Isco_neg(a),1/2))+a)/((np.power(R_Isco_neg(a),3/4))*(np.power((np.power(R_Isco_neg(a),3/2))-3*(np.power(R_Isco_neg(a),1/2))+2*(a),1/2))))

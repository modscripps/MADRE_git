#!/usr/bin/env python
#TODO add timestamp flag and numchannel in the header from the board
#     so we can use it to set the reader up


import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt


global fid
global blocks
global SBEcal

## local library
#  reads and apply calibration to the conductivity data
def get_CalSBE(filename='0133.cal'):
    class SBEcal_class:
        pass
    SBEcal=SBEcal_class
    fid=open('/Users/aleboyer/ARNAUD/SCRIPPS/SHEAR_PROBE/SBE49/' + filename)
    line=fid.readline()
    SBEcal.SN=line[-5:-1]
   ## Temperature Cal 
    line=fid.readline()
    SBEcal.TempCal_date=line[-10:-1]
    line=fid.readline()
    SBEcal.ta0=float(line[-14:-1])
    line=fid.readline()
    SBEcal.ta1=float(line[-14:-1])
    line=fid.readline()
    SBEcal.ta2=float(line[-14:-1])
    line=fid.readline()
    SBEcal.ta3=float(line[-14:-1])
    line=fid.readline()
    SBEcal.toffset=float(line[-14:-1])
   ## Conductivity Cal 
    line=fid.readline()
    SBEcal.CondCal_date=line[-10:-1]
    line=fid.readline()
    SBEcal.g=float(line[-14:-1])
    line=fid.readline()
    SBEcal.h=float(line[-14:-1])
    line=fid.readline()
    SBEcal.i=float(line[-14:-1])
    line=fid.readline()
    SBEcal.j=float(line[-14:-1])
    line=fid.readline()
    SBEcal.tcor=float(line[-14:-1])
    line=fid.readline()
    SBEcal.pcor=float(line[-14:-1])
    line=fid.readline()
    SBEcal.cslope=float(line[-14:-1])
   ## Pressure Cal
    line=fid.readline()
    SBEcal.PresCal_date=line[-10:-1]
    line=fid.readline()
    SBEcal.pa0=float(line[-14:-1])
    line=fid.readline()
    SBEcal.pa1=float(line[-14:-1])
    line=fid.readline()
    SBEcal.pa2=float(line[-14:-1])
    line=fid.readline()
    SBEcal.ptca0=float(line[-14:-1])
    line=fid.readline()
    SBEcal.ptca1=float(line[-14:-1])
    line=fid.readline()
    SBEcal.ptca2=float(line[-14:-1])
    line=fid.readline()
    SBEcal.ptcb0=float(line[-14:-1])
    line=fid.readline()
    SBEcal.ptcb1=float(line[-14:-1])
    line=fid.readline()
    SBEcal.ptcb2=float(line[-14:-1])
    line=fid.readline()
    SBEcal.ptempa0=float(line[-14:-1])
    line=fid.readline()
    SBEcal.ptempa1=float(line[-14:-1])
    line=fid.readline()
    SBEcal.ptempa2=float(line[-14:-1])
    line=fid.readline()
    SBEcal.poffset=float(line[-14:-1])
  
    return SBEcal 

 


def EPSI_temp(sample):
    T1  = Count2Volt_unipol(int(sample[:6],16))
    T2  = Count2Volt_unipol(int(sample[6:12],16))
    return T1,T2
def EPSI_shear(sample):
    S1  = Count2Volt_unipol(int(sample[12:18],16))
    S2  = Count2Volt_unipol(int(sample[18:24],16))
    return S1,S2 
#def EPSI_cond(sample):
#    C  = int(sample[24:30],16)
#    return C 
#def EPSI_accel(sample):
#    A1  = int(sample[30:36],16)
#    A2  = int(sample[36:42],16)
#    A3  = int(sample[42:48],16)
#    return A1,A2,A3 
def EPSI_accel(sample):
    A1  = Volt2g(Count2Volt_unipol(int(sample[24:30],16)))
    A2  = Volt2g(Count2Volt_unipol(int(sample[30:36],16)))
    A3  = Volt2g(Count2Volt_unipol(int(sample[36:42],16)))
    return A1,A2,A3 



def SBE_temp(SBEcal,SBEsample):
    a0 = SBEcal.ta0;
    a1 = SBEcal.ta1;
    a2 = SBEcal.ta2;
    a3 = SBEcal.ta3;
    
    rawT = int(SBEsample[:6],16);
    mv = (rawT-524288)/1.6e7;
    r = (mv*2.295e10 + 9.216e8)/(6.144e4-mv*5.3e5);
    temperature = a0+a1*np.log(r)+a2*np.log(r)**2+a3*np.log(r)**3;
    temperature = 1/temperature - 273.15;
    return temperature

def SBE_cond(SBEcal,SBEsample,pressure,temperature):
    g    = SBEcal.g;
    h    = SBEcal.h;
    i    = SBEcal.i;
    j    = SBEcal.j;
    tcor = SBEcal.tcor;
    pcor = SBEcal.pcor;
     
    f = int(SBEsample[6:12],16)/256/1000;
    conductivity = (g+h*f**2+i*f**3+j*f**4)/(1+tcor*temperature+pcor*pressure);
 
    return conductivity
 
#  reads and apply calibration to the pressure data
def SBE_Pres(SBEcal,SBEsample):
    pa0     = SBEcal.pa0;
    pa1     = SBEcal.pa1;
    pa2     = SBEcal.pa2;
    ptempa0 = SBEcal.ptempa0;
    ptempa1 = SBEcal.ptempa1;
    ptempa2 = SBEcal.ptempa2;
    ptca0   = SBEcal.ptca0;
    ptca1   = SBEcal.ptca1;
    ptca2   = SBEcal.ptca2;
    ptcb0   = SBEcal.ptcb0;
    ptcb1   = SBEcal.ptcb1;
    ptcb2   = SBEcal.ptcb2;
    
    rawP = int(SBEsample[12:18],16);
    y    = int(SBEsample[18:22],16)/13107;
 
    t = ptempa0+ptempa1*y+ptempa2*y**2;
    x = rawP-ptca0-ptca1*t-ptca2*t**2;
    n = x*ptcb0/(ptcb0+ptcb1*t+ptcb2*t**2);
     
    pressure = (pa0+pa1*n+pa2*n**2-14.7)*0.689476;
     
    return pressure


def Count2Volt_unipol(count,N=24,Full_Range=2.5,Gain=1):
    Vin=Full_Range*np.array(count)/2**N/Gain;
    return Vin;

def Count2Volt_bipol(count,N=24,Full_Range=2.5,Gain=1):
    Vin=Full_Range/Gain*(np.array(count)/2**(N-1)-1);
    return Vin;

def Volt2g(V,offset=1.65):
    g_in=(np.array(V)-offset)/.66;
    return g_in

## end local library TO DO: move it somewhere else

## local library for plotting


def read_CTDsample(SBEcal,SBEsample):
    P     = SBE_Pres(SBEcal,SBEsample,)
    T     = SBE_temp(SBEcal,SBEsample)
    C     = SBE_cond(SBEcal,SBEsample,P,T)
    return P,T,C
    
def read_EPSIblock(EPSIblock):
    Temp    = np.asarray(EPSI_temp(EPSIblock))
    Shear   = np.asarray(EPSI_shear(EPSIblock))
#    Cond    = EPSI_cond(sample)
    Accel   = np.asarray(EPSI_accel(EPSIblock))

    t1=Temp[0]
    t2=Temp[1]
    s1=Shear[0]
    s2=Shear[1]
    a1=Accel[0]
    a2=Accel[1]
    a3=Accel[2]
    
    return t1,t2,s1,s2,a1,a2,a3


def open_datafile(filename='SBE_EPSI_hex_1.dat'):
    fid=open('/Users/aleboyer/ARNAUD/SCRIPPS/SHEAR_PROBE/PROCESS_EPSI/data/' + filename,'rb')
    eof=fid.seek(0,2)
    fid.seek(0) 
    print('reading file ... NB: if too long change to mmap?')
    lines=fid.read() 
    blocks=lines.split(b'##')
    return fid,eof,blocks
    
    
def read_datablock(blocks):
    lenblock=[len(i) for i in blocks]
    flag_block=[i==995 for i in lenblock]
    OKblocks=[i.split(b'\r\n') for i in blocks if len(i)==995]
    
    
    timeaxis =np.zeros(len(OKblocks))
    P =np.zeros(len(OKblocks))
    T =np.zeros(len(OKblocks))
    C =np.zeros(len(OKblocks))
    SBEn =np.zeros(len(OKblocks))
    timeEPSI = []
    EPSI=[]
    starttime=float(OKblocks[0][1])
    for ind,i in enumerate(OKblocks):
        SBEn[ind]=int(i[0].split(b',')[0],16)
        timeaxis[ind]= starttime+SBEn[ind]/16
        [P[ind],T[ind],C[ind]]= read_CTDsample(SBEcal,i[2])
        if i[-2][1]==33:
            EPSI+=i[3:-2]
            timeEPSI+=(timeaxis[ind]+np.arange(len(i[3:-2]))/len(i[3:-2])/16).tolist()
        else:
            EPSI+=i[3:-1]
            timeEPSI+=(timeaxis[ind]+np.arange(len(i[3:-1]))/len(i[3:-1])/16).tolist()
    
    t1 =np.zeros(len(EPSI))
    t2 =np.zeros(len(EPSI))
    s1 =np.zeros(len(EPSI))
    s2 =np.zeros(len(EPSI))
    a1 =np.zeros(len(EPSI))
    a2 =np.zeros(len(EPSI))
    a3 =np.zeros(len(EPSI))
    
    for ind,i in enumerate(EPSI):
        [t1[ind],t2[ind],s1[ind],s2[ind],a1[ind],a2[ind],a3[ind]]=read_EPSIblock(i)

    return timeEPSI,t1,t2,s1,s2,a1,a2,a3,P,T,C,timeaxis    
    

fid,eof,blocks=open_datafile('SPROULnovember/SBE_EPSI_hex_0.dat')
SBEcal=get_CalSBE()
timeEPSI,t1,t2,s1,s2,a1,a2,a3,P,T,C,timeCTD=read_datablock(blocks)

EPSImat={'time':timeEPSI, \
         'Sensor1':t1,  \
         'Sensor2':t2,  \
         'Sensor3':s1,  \
         'Sensor4':s2,  \
         'Sensor5':[0], \
         'Sensor6':a1,  \
         'Sensor7':a2,  \
         'Sensor8':a3}
sio.savemat('../data/EPSITEST_sproul_epsi1.mat',EPSImat)

SBEmat={'time':timeCTD,\
        'T':T,\
        'C':C,\
        'P':P}
sio.savemat('../data/EPSITEST_sproul_sbe1.mat',SBEmat)




#!/usr/bin/env python

##  only look at shear and temp from madre
import numpy as np
import time
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.widgets import RadioButtons,Slider
from matplotlib import style

# style of graphes
style.use('fivethirtyeight')

global SBEsample
global EPSIsample
global fid
global start_time
global start_Epsicount
global RTCM1
global timestamp0
global byte_per_sample
global ax0
global epsisample_per_block
global aux1sample_per_block
global channel

## local library
#  reads and apply calibration to the conductivity data
def get_CalSBE(filename='SBE49_4935239-0058_cal.dat'):
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
    T1  = int(sample[:6],16)
    T2  = int(sample[6:12],16)
    return T1,T2
def EPSI_shear(sample):
    S1  = int(sample[12:18],16)
    S2  = int(sample[18:24],16)
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
    A1  = int(sample[24:30],16)
    A2  = int(sample[30:36],16)
    A3  = int(sample[36:42],16)
    return A1,A2,A3 



def SBE_temp(SBEcal,SBEsample):
    a0 = SBEcal.ta0;
    a1 = SBEcal.ta1;
    a2 = SBEcal.ta2;
    a3 = SBEcal.ta3;
    
    rawT = int(SBEsample.raw[:6],16);
    mv = (rawT-524288)/1.6e7;
    r = (mv*2.295e10 + 9.216e8)/(6.144e4-mv*5.3e5);
    SBEsample.temperature = a0+a1*np.log(r)+a2*np.log(r)**2+a3*np.log(r)**3;
    SBEsample.temperature = 1/SBEsample.temperature - 273.15;
    return SBEsample 

def SBE_cond(SBEcal,SBEsample):
    g    = SBEcal.g;
    h    = SBEcal.h;
    i    = SBEcal.i;
    j    = SBEcal.j;
    tcor = SBEcal.tcor;
    pcor = SBEcal.pcor;
     
    f = int(SBEsample.raw[6:12],16)/256/1000;
    SBEsample.conductivity = (g+h*f**2+i*f**3+j*f**4)/(1+tcor*SBEsample.temperature+pcor*SBEsample.pressure);
 
    return SBEsample
 
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
    
    rawP = int(SBEsample.raw[12:18],16);
    y    = int(SBEsample.raw[18:22],16)/13107;
 
    t = ptempa0+ptempa1*y+ptempa2*y**2;
    x = rawP-ptca0-ptca1*t-ptca2*t**2;
    n = x*ptcb0/(ptcb0+ptcb1*t+ptcb2*t**2);
     
    SBEsample.pressure = (pa0+pa1*n+pa2*n**2-14.7)*0.689476;
     
    return SBEsample


def Count2Volt_unipol(count,N=24,Full_Range=2.5,Gain=1):
    Vin=Full_Range*np.array(count)/2**N/Gain;
    return Vin;
def Volt2Count_unipol(Vin,N=24,Full_Range=2.5,Gain=1):
    counts=(2**N * Vin * Gain) / Full_Range;
    return counts;

def Count2Volt_bipol(count,N=24,Full_Range=2.5,Gain=1):
    Vin=Full_Range/Gain*(np.array(count)/2**(N-1)-1);
    return Vin;
def Volt2count_bipol(Vin,N=24,Full_Range=2.5,Gain=1):
    counts=2**N * ( (Vin * Gain) / Full_Range +1 );
    return counts;

def Volt2g(V,offset=1.65):
    g_in=(np.array(V)-offset)/.66;
    return g_in

## end local library TO DO: move it somewhere else

## local library for plotting

def init_figure(num_sub):
    fig = plt.figure(0)
    #fig1 = plt.figure(1)
    #fig.add_subplot(1,1,1,autoscale_on=True)
    return fig

def init_axes(fig):
    ax0=plt.axes([.3,.4,.5,.5])
    ax1=plt.axes([.1,.1,.15,.2])
    ax2=plt.axes([.3,.01,.5,.05])
    ax3=plt.axes([.3,.1,.5,.05])
    ax4=plt.axes([.1,.4,.15,.5])
    # set up subplot 
    ax0.cla()
    lineF0, = ax0.plot([],[], '.-', alpha=0.8, color="gray", markerfacecolor="red")
    lineV1, = ax0.plot([],[], '.-', alpha=0.8, color="gray", markerfacecolor="green")
    lineV2, = ax0.plot([],[], '.-', alpha=0.8, color="gray", markerfacecolor="cyan")
    radio = RadioButtons(ax1, ('Counts', 'Volts')) 
    radio1 = RadioButtons(ax4, ('t1', 't2','s1','s2','a1','a2','a3'))
    Sl_ymin = Slider(ax2, 'Ymin', 0, 1, valinit=.5)
    Sl_ymax = Slider(ax3, 'Ymax', 0, 1, valinit=.5)
    ax0.set_xlabel('Second')
    ax0.set_ylabel('Counts')
    return lineF0,lineV1,lineV2,radio,radio1,Sl_ymin,Sl_ymax,ax0



def open_datafile(filename='MADRE2.1.dat'):
    fid=open('../data/' + filename,'br')
    eof=fid.seek(0,2)
    fid.seek(eof-8*1050) 
    line=fid.readline()
#    while(len(line)!=29):
    while(line[:6]!=b'$MADRE'):
          line=fid.readline()
          print('coucou open')
    fid.seek(fid.tell()-61)
    #fid.seek(0)
    return fid,eof
    
    
def read_header(fid):
    class header_class:
          pass
    header=header_class
    header.raw=fid.readline()
    header.EpsiStamp      = int(header.raw[6:6+8],16)
    header.TimeStamp      = int(header.raw[15:15+8],16)
    header.Voltage        = int(header.raw[24:24+8],16)
    header.aux1_checksum  = int(header.raw[33:33+8],16)
    header.aux2_checksum  = int(header.raw[42:42+8],16)
    header.block_checksum = int(header.raw[51:51+8],16)

    header.local_timestamp= time.time()
    return header    

def read_block(fid):   
    block=[fid.readline() for i in range(epsisample_per_block)] # 24 is for the length of the SBE49 sample
    return block
def read_aux1sample(fid):   
    block=[fid.readline() for i in range(aux1sample_per_block)] # 24 is for the length of the SBE49 sample
    return block


def Aux1block2sample(SBEcal,Auxsample,Auxblock):
    Auxsample.raw = Auxblock[:22]
    Auxsample     = SBE_Pres(SBEcal,SBEsample)
    Auxsample     = SBE_temp(SBEcal,SBEsample)
    Auxsample     = SBE_cond(SBEcal,SBEsample)
    return SBEsample

 
 
def Epsiblock2sample(EPSIsample,Epsiblock):

    EPSIsample.raw = Epsiblock
    epsisamples    = EPSIsample.raw  ## TODO change because it only works when madre send hex
    
    EPSIsample.T1 = [];
    EPSIsample.T2 = [];
    EPSIsample.S1 = [];
    EPSIsample.S2 = [];
    EPSIsample.C  = [];
    EPSIsample.A1 = [];
    EPSIsample.A2 = [];
    EPSIsample.A3 = [];
    for (i,sample) in enumerate(epsisamples): 
           t1,t2    = EPSI_temp(sample[:-2])
           s1,s2    = EPSI_shear(sample[:-2])
#           c        = EPSI_cond(sample[:-2])
           a1,a2,a3 = EPSI_accel(sample[:-2])
           EPSIsample.T1.append(t1);   EPSIsample.T2.append(t2);           
           EPSIsample.S1.append(s1);   EPSIsample.S2.append(s2);          
#           EPSIsample.C.append(c);           
           EPSIsample.A1.append(a1);EPSIsample.A2.append(a2);           
           EPSIsample.A3.append(a3);           
           
    return EPSIsample

def update_sample(fid,EPSIsample,SBEsample):
    
    header    = read_header(fid) 
    Epsiblock = read_block (fid) 
    #timestamp=header.local_timestamp
    EPSIsample.RTC=np.append(EPSIsample.RTC[1:].tolist(),header.TimeStamp)
    EPSI_freq=32768.0*epsisample_per_block/np.mean(np.diff(EPSIsample.RTC))
    ##time axis
    local_EPSItime=EPSIsample.timeaxis[-1]+np.arange(1,epsisample_per_block+1)/EPSI_freq
    

    EPSIsample.timeaxis=np.array(EPSIsample.timeaxis[epsisample_per_block:].tolist() \
                                 + local_EPSItime.tolist())
    
    ## get samples
    EPSIsample=Epsiblock2sample(EPSIsample,Epsiblock)
    
    ## concatanate EPSI sample 
    # fpo7
    local_var=EPSIsample.T1
    EPSIsample.temp1=np.append(EPSIsample.temp1[epsisample_per_block:],local_var)
    local_var=EPSIsample.T2
    EPSIsample.temp2=np.append(EPSIsample.temp2[epsisample_per_block:],local_var)

    # shear            
    local_var=EPSIsample.S1
    EPSIsample.shear1=np.append(EPSIsample.shear1[epsisample_per_block:],local_var)
    local_var=EPSIsample.S2
    EPSIsample.shear2=np.append(EPSIsample.shear2[epsisample_per_block:],local_var)

    # micro conductivity
    local_var=Count2Volt_bipol(EPSIsample.C)
    EPSIsample.conductivity=np.array(EPSIsample.conductivity[epsisample_per_block:].tolist()+ \
                                     local_var.tolist())

    # acelleration
    local_var=np.copy(EPSIsample.A1)
    EPSIsample.accelx=np.append(EPSIsample.accelx[epsisample_per_block:],local_var)
    local_var=EPSIsample.A2
    EPSIsample.accely=np.append(EPSIsample.accely[epsisample_per_block:],local_var)
    local_var=EPSIsample.A3
    EPSIsample.accelz=np.append(EPSIsample.accelz[epsisample_per_block:],local_var)

    ##concatanate SBE samples
#    SBEsample.T = np.array(SBEsample.T[1:].tolist()+[SBEsample.temperature])
#    SBEsample.C = np.array(SBEsample.C[1:].tolist()+[SBEsample.conductivity])
#    SBEsample.P = np.array(SBEsample.P[1:].tolist()+[SBEsample.pressure])

    return EPSIsample,SBEsample

def animate(i,radio,radio1,Sl_ymin,Sl_ymax,ax0,index1):

    if(radio1.value_selected=='t1'):
       data=EPSIsample.temp1
    if(radio1.value_selected=='t2'):
       data=EPSIsample.temp2
    if(radio1.value_selected=='s1'):
       data=EPSIsample.shear1
    if(radio1.value_selected=='s2'):
       data=EPSIsample.shear2
    if(radio1.value_selected=='a1'):
       data=EPSIsample.accelx
    if(radio1.value_selected=='a2'):
       data=EPSIsample.accely
    if(radio1.value_selected=='a3'):
       data=EPSIsample.accelz

    ax0xmax=0xffffff
    ax0xmin=0
    if(radio.value_selected=='Volts'):
       data=Count2Volt_unipol(data)
       ax0xmax=3.3
       ax0xmin=0
    
    ax0xmin=Sl_ymin.val*0xffffff 
    ax0xmax=Sl_ymax.val*0xffffff 
    if(radio.value_selected=='Volts'):
    	ax0xmin=Count2Volt_unipol(Sl_ymin.val*0xffffff) 
    	ax0xmax=Count2Volt_unipol(Sl_ymax.val*0xffffff) 
    	ax0.set_ylabel('Volts')
    else: 
        ax0.set_ylabel('Counts')
    
    posi=fid.tell()
    eof=fid.seek(0,2)
    fid.seek(posi)

    if(i==0):
        index1=0
        
    if (eof-posi)>2*EPSIWordlength:
        update=True
        print('update')
        index1=0
    else:
        update=False
        print('trou du cul')
        index1+=1

    print(i)
    if update:
        update_sample(fid,EPSIsample,SBEsample)

    lineF0.set_data(EPSIsample.timeaxis[index1:index1+Lplot],data[index1:index1+Lplot])
    lineV1.set_data(EPSIsample.timeaxis[index1:index1+Lplot],EPSIsample.timeaxis[index1:index1+Lplot]*0+max(data[index1:index1+Lplot]))
    lineV2.set_data(EPSIsample.timeaxis[index1:index1+Lplot],EPSIsample.timeaxis[index1:index1+Lplot]*0+min(data[index1:index1+Lplot]))
     
    tmin1=min(EPSIsample.timeaxis[index1:index1+Lplot])
    tmax1=max(EPSIsample.timeaxis[index1:index1+Lplot])

    fig.axes[0].set_xlim([tmin1,tmax1])   
    fig.axes[0].set_ylim([ax0xmin,ax0xmax]) 
    #fig.axes[0].text(tmin1+.9*(tmax1-tmin1),ax0xmin+.9*(ax0xmax-ax0xmin),'rms=%3.3f' % (np.sqrt(np.mean(data**2))))  
    fig.axes[0].legend(['rms=%3.3f' % np.sqrt(np.mean(data**2))],loc=2)  
    
    return lineF0,lineV1,lineV2,index1


    

#####################################################################



# define classe SBE sample and epsi sample
# define classe SBE sample and epsi sample
class SBEsample_class:
    pass
class EPSIsample_class:
    pass

EPSIsample  = EPSIsample_class
SBEsample   = SBEsample_class


#SBEcal      = get_CalSBE()

fid,eof=open_datafile()


Aux1WordLength = 0
ADCWordlength  = 3
number_of_sensor  = 7
EpsisampleWordLength= ADCWordlength*number_of_sensor
epsisample_per_block  = 160 
aux1sample_per_block=0
EPSIWordlength = ADCWordlength *  number_of_sensor * epsisample_per_block

header    = read_header(fid) # reader first header to get time stamp,
Epsiblock = read_block (fid) 

RTCm1 =header.TimeStamp

print('header')
print(header.raw)

start_time=header.local_timestamp
start_Epsicount=header.EpsiStamp
LBuffer=320 # number of point in buffer 
Lplot  =160 # number of point in the plot 

index1=0

#SBEsample.timeaxis  = np.zeros(Lplot)+header.timestamp
#SBEsample.T  = np.zeros(Lplot)
#SBEsample.C  = np.zeros(Lplot)
#SBEsample.P  = np.zeros(Lplot)

EPSIsample.timeaxis  = np.zeros(LBuffer)
EPSIsample.RTC       = np.zeros(3)
EPSIsample.temp1     = np.zeros(LBuffer)
EPSIsample.temp2     = np.zeros(LBuffer)
EPSIsample.shear1    = np.zeros(LBuffer)
EPSIsample.shear2    = np.zeros(LBuffer)
EPSIsample.accelx    = np.zeros(LBuffer)
EPSIsample.accely    = np.zeros(LBuffer)
EPSIsample.accelz    = np.zeros(LBuffer)
EPSIsample.conductivity    = np.zeros(LBuffer)
  


fig=init_figure(1)
lineF0,lineV1,lineV2,radio,radio1,Sl_ymin,Sl_ymax,ax0=init_axes(fig)
#ani = animation.FuncAnimation(fig, animate,frames=100,interval=10, blit=True)
ani = animation.FuncAnimation(fig, animate,fargs=(radio,radio1,Sl_ymin,Sl_ymax,ax0,index1),interval=1)
plt.show()
   




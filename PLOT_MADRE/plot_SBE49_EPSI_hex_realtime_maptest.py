#!/usr/bin/env python

##  only look at shear and temp from madre
import numpy as np
import sys
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
global byte_per_sample
global ax0

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



def open_datafile(filename='SBE_EPSI_hex.dat'):
    fid=open('../data/' + filename,'br')
    eof=fid.seek(0,2)
    fid.seek(eof-8*1050) 
    line=fid.readline()
#    while(len(line)!=29):
    while(len(line)!=30):
          line=fid.readline()
          print('coucou open')
    fid.seek(fid.tell()-30)
    #fid.seek(0)
    return fid,eof
    
    
def read_datablock(fid):
    class header_class:
          pass
    header=header_class
    header.raw=fid.readline()
    header.byte_per_block = int(header.raw[-9:-1],16)
    header.SBEcount       = int(header.raw[2:10],16)
    header.TX_block_size  = int(header.raw[12:19],16)
    
    block=fid.read(header.TX_block_size+24) # 24 is for the length of the SBE49 sample
    # 50 is the number of byte per sample
    line=fid.readline()
    if (len(line)==30 and line[:2]==b'##'):
       valid_block=1
       fid.seek(fid.tell()-30)
       header.lengthBlock=fid.tell()
       header.timestamp=header.SBEcount
    else:
       valid_block=0
       print('coucou bad block')
       print(line)
       print(header.TX_block_size)
       print(header.raw)
       print(block)
#       while(len(line)!=29):
       while(len(line)!=30):
          line=fid.readline()
       fid.seek(fid.tell()-30)
    return valid_block,block,header    
 
 
def block2sample(SBEcal,SBEsample,EPSIsample,block,header):
    SBEsample.raw = block[:22]
    SBEsample     = SBE_Pres(SBEcal,SBEsample)
    SBEsample     = SBE_temp(SBEcal,SBEsample)
    SBEsample     = SBE_cond(SBEcal,SBEsample)

    EPSIsample.raw = block[23:]
    epsisamples    = EPSIsample.raw.split(b'\n')  ## TODO change because it only works when madre send hex
    
    EPSIsample.T1 = [];
    EPSIsample.T2 = [];
    EPSIsample.S1 = [];
    EPSIsample.S2 = [];
    EPSIsample.C  = [];
    EPSIsample.A1 = [];
    EPSIsample.A2 = [];
    EPSIsample.A3 = [];
    for (i,sample) in enumerate(epsisamples[1:-1]): 
        if (i*byte_per_sample<header.byte_per_block):
           t1,t2    = EPSI_temp(sample)
           s1,s2    = EPSI_shear(sample)
#           c        = EPSI_cond(sample)
           a1,a2,a3 = EPSI_accel(sample)
           EPSIsample.T1.append(t1);   EPSIsample.T2.append(t2);           
           EPSIsample.S1.append(s1);   EPSIsample.S2.append(s2);          
#           EPSIsample.C.append(c);           
           EPSIsample.A1.append(a1);EPSIsample.A2.append(a2);           
           EPSIsample.A3.append(a3);           
           
    return EPSIsample,SBEsample

def update_sample(fid,EPSIsample,SBEsample):
    valid_block,block,header = read_datablock(fid) 
    if valid_block==1: # good block to read
    
        # create timestamp here using the number of SBE sample
        # TODO create time stamp on the board
        # rollover is used when the number of sbe samples 
        # is higher than 0xffffffff (uint32_t)
        nb_EPSIsample=int(header.byte_per_block/byte_per_sample)
        timestamp0=SBEsample.timeaxis[-1]
        timestamp=header.timestamp
        if timestamp0==0:
            EPSI_freq=325.0
        else:
            EPSI_freq=(header.byte_per_block/byte_per_sample)/(timestamp-timestamp0)
        
        ##time axis
        SBEsample.timeaxis=np.array(SBEsample.timeaxis[1:].tolist()+[timestamp])
        local_EPSItime=timestamp+np.arange(1,header.byte_per_block/byte_per_sample+1)/EPSI_freq
        EPSIsample.timeaxis=np.array(EPSIsample.timeaxis[nb_EPSIsample:].tolist() + local_EPSItime.tolist())
        
        ## get samples
        EPSIsample,SBEsample=block2sample(SBEcal,SBEsample,EPSIsample,block,header)
        
        ##concatanate SBE samples
        SBEsample.T = np.array(SBEsample.T[1:].tolist()+[SBEsample.temperature])
        SBEsample.C = np.array(SBEsample.C[1:].tolist()+[SBEsample.conductivity])
        SBEsample.P = np.array(SBEsample.P[1:].tolist()+[SBEsample.pressure])
        ## concatanate EPSI sample 
        # fpo7
        local_var=EPSIsample.T1
        #EPSIsample.temp1=np.array(EPSIsample.temp1[nb_EPSIsample:].tolist()+local_var.tolist())
        EPSIsample.temp1=np.append(EPSIsample.temp1[nb_EPSIsample:],local_var)
        local_var=EPSIsample.T2
        #EPSIsample.temp2=np.array(EPSIsample.temp2[nb_EPSIsample:].tolist()+local_var.tolist())
        EPSIsample.temp2=np.append(EPSIsample.temp2[nb_EPSIsample:],local_var)
        # shear            
        local_var=EPSIsample.S1
        #EPSIsample.shear1=np.array(EPSIsample.shear1[nb_EPSIsample:].tolist()+local_var.tolist())
        EPSIsample.shear1=np.append(EPSIsample.shear1[nb_EPSIsample:],local_var)
        local_var=EPSIsample.S2
        #EPSIsample.shear2=np.array(EPSIsample.shear2[nb_EPSIsample:].tolist()+local_var.tolist())
        EPSIsample.shear2=np.append(EPSIsample.shear2[nb_EPSIsample:],local_var)
        # micro conductivity
        local_var=Count2Volt_bipol(EPSIsample.C)
        EPSIsample.conductivity=np.array(EPSIsample.conductivity[nb_EPSIsample:].tolist()+local_var.tolist())
        # acelleration
        #local_var=Volt2g(Count2Volt_unipol(EPSIsample.A1))
        local_var=EPSIsample.A1
        EPSIsample.accelx=np.append(EPSIsample.accelx[nb_EPSIsample:],local_var)
        #local_var=Volt2g(Count2Volt_unipol(EPSIsample.A2))
        local_var=EPSIsample.A2
        EPSIsample.accely=np.append(EPSIsample.accely[nb_EPSIsample:],local_var)
        #local_var=Volt2g(Count2Volt_unipol(EPSIsample.A3))
        local_var=EPSIsample.A3
        EPSIsample.accelz=np.append(EPSIsample.accelz[nb_EPSIsample:],local_var)
    return EPSIsample,SBEsample,timestamp

def animate(i,radio,radio1,Sl_ymin,Sl_ymax,ax0):

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
    while (eof-posi)>4*1050:
        update_sample(fid,EPSIsample,SBEsample)
        posi=fid.tell()
    

    lineF0.set_data(EPSIsample.timeaxis-start_time,data)
    lineV1.set_data(EPSIsample.timeaxis-start_time,EPSIsample.timeaxis*0+max(data))
    lineV2.set_data(EPSIsample.timeaxis-start_time,EPSIsample.timeaxis*0+min(data))
    
     
    tmin1=min(EPSIsample.timeaxis-start_time)
    tmax1=max(EPSIsample.timeaxis-start_time)


    
    fig.axes[0].set_xlim([tmin1,tmax1])   
    fig.axes[0].set_ylim([ax0xmin,ax0xmax]) 
    #fig.axes[0].text(tmin1+.9*(tmax1-tmin1),ax0xmin+.9*(ax0xmax-ax0xmin),'rms=%3.3f' % (np.sqrt(np.mean(data**2))))  
    fig.axes[0].legend(['rms=%3.3f' % np.sqrt(np.mean(data**2))],loc=2)  
    return lineF0,lineV1,lineV2


    

#####################################################################



# define classe SBE sample and epsi sample
# define classe SBE sample and epsi sample
class SBEsample_class:
    pass
class EPSIsample_class:
    pass


SBEsample   = SBEsample_class
EPSIsample  = EPSIsample_class
byte_per_sample=44


SBEcal      = get_CalSBE()

fid,eof=open_datafile()
valid_block,block,header = read_datablock(fid) # reader first header to get time stamp, 
print('header')
print(header.raw)

start_time=header.timestamp

#fid.seek(0)# go back to begin of the file
nb_of_sample=(eof-(eof%header.lengthBlock))/header.lengthBlock
Lplot=500 # number of point in the plot 
#Lplot=64 # number of point in the plot 




SBEsample.timeaxis  = np.zeros(Lplot)+header.timestamp
SBEsample.T  = np.zeros(Lplot)
SBEsample.C  = np.zeros(Lplot)
SBEsample.P  = np.zeros(Lplot)

EPSIsample.timeaxis  = np.zeros(Lplot)+header.timestamp
EPSIsample.temp1     = np.zeros(Lplot)
EPSIsample.temp2     = np.zeros(Lplot)
EPSIsample.shear1    = np.zeros(Lplot)
EPSIsample.shear2    = np.zeros(Lplot)
EPSIsample.accelx    = np.zeros(Lplot)
EPSIsample.accely    = np.zeros(Lplot)
EPSIsample.accelz    = np.zeros(Lplot)
EPSIsample.conductivity    = np.zeros(Lplot)
  


fig=init_figure(1)
lineF0,lineV1,lineV2,radio,radio1,Sl_ymin,Sl_ymax,ax0=init_axes(fig)
#ani = animation.FuncAnimation(fig, animate,frames=100,interval=10, blit=True)
ani = animation.FuncAnimation(fig, animate,fargs=(radio,radio1,Sl_ymin,Sl_ymax,ax0,),interval=1)
plt.show()
   




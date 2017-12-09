#!/usr/bin/env python

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import style

# style of graphes
style.use('fivethirtyeight')

global SBEsample
global EPSIsample
global fid
global start_time
global byte_per_sample

## local library
#  reads and apply calibration to the conductivity data
def get_CalSBE(filename='SBE49_4935239-0058_cal.dat'):
    class SBEcal_class:
        pass
    SBEcal=SBEcal_class
    fid=open('../SBE49/' + filename)
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

def Count2Volt_bipol(count,N=24,Full_Range=2.5,Gain=1):
    Vin=Full_Range/Gain*(np.array(count)/2**(N-1)-1);
    return Vin;

def Volt2g(V,offset=1.65):
    g_in=(np.array(V)-offset)/.66;
    return g_in

## end local library TO DO: move it somewhere else

## local library for plotting

def init_figure(num_sub):
    fig = plt.figure(0)
#    xlim_min=[0,0,-.3]
#    xlim_max=[30,1e-3,0]
#    ylim_min=[]
#    ylim_max=[]
    for i in range(num_sub):
#        fig.add_subplot(1,num_sub, i+1,autoscale_on=False,xlim(xlim_min[i],xlim_max[i]))
        fig.add_subplot(1,num_sub, i+1,autoscale_on=True)
    return fig

def init_axes(fig):
    ax0=fig.axes[0]
    ax1=fig.axes[1]
    ax2=fig.axes[2]
    ax3=fig.axes[3]
    ax4=fig.axes[4]
    ax5=fig.axes[5]
    # set up subplot 
    ax0.cla()
    ax0.set_title('SBE T')
    ax0.set_ylabel('Time (s)')
    ax0.set_xlabel('C')
    line0, = ax0.plot([],[], '.-', alpha=0.8, color="gray", markerfacecolor="red")
    ax1.cla()
    ax1.set_title('SBE C')
    ax1.set_xlabel('')
    ax1.set_yticklabels('')
    line1, = ax1.plot([],[], '.-', alpha=0.8, color="gray", markerfacecolor="red")
    
    ax2.cla()
    ax2.set_title('SBE Pressure')
    ax2.set_xlabel('dB')
    ax2.set_yticklabels('')
    line2, = ax2.plot([],[], '.-', alpha=0.8, color="gray", markerfacecolor="red")

    ax3.cla()
    ax3.set_title('FPO7')
    ax3.set_xlabel('Volt')
    ax3.set_yticklabels('')
    lineF1, = ax3.plot([],[], '.-', alpha=0.8, color="gray", markerfacecolor="red")
    lineF2, = ax3.plot([],[], '.-', alpha=0.8, color="gray", markerfacecolor="green")
    
    ax4.cla()
    ax4.set_title('Shear')
    ax4.set_xlabel('Volt')
    ax4.set_yticklabels('')
    lineS1, = ax4.plot([],[], '.-', alpha=0.8, color="gray", markerfacecolor="red")
    lineS2, = ax4.plot([],[], '.-', alpha=0.8, color="gray", markerfacecolor="green")
    
    ax5.cla()
    ax5.set_title('Acceleration ')
    ax5.set_xlabel('G')
    ax5.set_yticklabels('')
    lineA1, = ax5.plot([],[], '.-', alpha=0.8, color="gray", markerfacecolor="red")
    lineA2, = ax5.plot([],[], '.-', alpha=0.8, color="gray", markerfacecolor="green")
    lineA3, = ax5.plot([],[], '.-', alpha=0.8, color="gray", markerfacecolor="blue")
    
    return line0,line1,line2,lineF1,lineF2,lineS1,lineS2,lineA1,lineA2,lineA3



def open_datafile(filename='SBE_EPSI_hex.dat'):
    fid=open('../data/' + filename,'r')
    eof=fid.seek(0,2)
    fid.seek(eof-8*1050) 
    line=fid.readline()
    while(len(line)!=29):
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
    tempo_timestamp       = fid.readline()
    
    block=fid.read(header.TX_block_size+23-int(header.TX_block_size/byte_per_sample)) # 24 is for the length of the SBE49 sample
    # 50 is the number of byte per sample
    line=fid.readline()
    print(len(line))
    if (len(line)==29 and line[:2]=='##'):
       valid_block=1
       fid.seek(fid.tell()-30)
       header.lengthBlock=fid.tell()
       header.timestamp=float(tempo_timestamp[:-1])
    else:
       valid_block=0
       print('coucou bad block')
       print(line)
       print(header.TX_block_size)
       print(header.raw)
       print(block)
       while(len(line)!=29):
          line=fid.readline()
       fid.seek(fid.tell()-30)
    return valid_block,block,header    
 
 
def block2sample(SBEcal,SBEsample,EPSIsample,block,header):
    SBEsample.raw = block[:22]
    SBEsample     = SBE_Pres(SBEcal,SBEsample)
    SBEsample     = SBE_temp(SBEcal,SBEsample)
    SBEsample     = SBE_cond(SBEcal,SBEsample)

    EPSIsample.raw = block[23:]
    epsisamples    = EPSIsample.raw.split('\n')  ## TODO change because it only works when madre send hex
    
    EPSIsample.T1 = [];
    EPSIsample.T2 = [];
    EPSIsample.S1 = [];
    EPSIsample.S2 = [];
    EPSIsample.C  = [];
    EPSIsample.A1 = [];
    EPSIsample.A2 = [];
    EPSIsample.A3 = [];
    for (i,sample) in enumerate(epsisamples): 
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
        local_var=Count2Volt_bipol(EPSIsample.T1)
        EPSIsample.temp1=np.array(EPSIsample.temp1[nb_EPSIsample:].tolist()+local_var.tolist())
        local_var=Count2Volt_bipol(EPSIsample.T2)
        EPSIsample.temp2=np.array(EPSIsample.temp2[nb_EPSIsample:].tolist()+local_var.tolist())
        # shear            
        local_var=Count2Volt_bipol(EPSIsample.S1)
        EPSIsample.shear1=np.array(EPSIsample.shear1[nb_EPSIsample:].tolist()+local_var.tolist())
        local_var=Count2Volt_bipol(EPSIsample.S2)
        EPSIsample.shear2=np.array(EPSIsample.shear2[nb_EPSIsample:].tolist()+local_var.tolist())
        # micro conductivity
        local_var=Count2Volt_bipol(EPSIsample.C)
        EPSIsample.conductivity=np.array(EPSIsample.conductivity[nb_EPSIsample:].tolist()+local_var.tolist())
        # acelleration
        local_var=Volt2g(Count2Volt_unipol(EPSIsample.A1))
        EPSIsample.accelx=np.array(EPSIsample.accelx[nb_EPSIsample:].tolist()+local_var.tolist())
        local_var=Volt2g(Count2Volt_unipol(EPSIsample.A2))
        EPSIsample.accely=np.array(EPSIsample.accely[nb_EPSIsample:].tolist()+local_var.tolist())
        local_var=Volt2g(Count2Volt_unipol(EPSIsample.A3))
        EPSIsample.accelz=np.array(EPSIsample.accelz[nb_EPSIsample:].tolist()+local_var.tolist())
    return EPSIsample,SBEsample,timestamp

    
def animate(i):
    posi=fid.tell()
    eof=fid.seek(0,2)
    fid.seek(posi)
    while (eof-posi)>4*1050:
        update_sample(fid,EPSIsample,SBEsample)
        posi=fid.tell()
            
    line0.set_data(SBEsample.T,SBEsample.timeaxis-start_time)
    line1.set_data(SBEsample.C,SBEsample.timeaxis-start_time)
    line2.set_data(SBEsample.P,SBEsample.timeaxis-start_time)
    
    print(fid.tell())
    print(eof)

    lineF1.set_data(EPSIsample.temp1,EPSIsample.timeaxis-start_time)
    lineF2.set_data(EPSIsample.temp2,EPSIsample.timeaxis-start_time)
    lineS1.set_data(EPSIsample.shear1,EPSIsample.timeaxis-start_time)
    lineS2.set_data(EPSIsample.shear2,EPSIsample.timeaxis-start_time)
    lineA1.set_data(EPSIsample.accelx,EPSIsample.timeaxis-start_time)
    lineA2.set_data(EPSIsample.accely,EPSIsample.timeaxis-start_time)
    lineA3.set_data(EPSIsample.accelz,EPSIsample.timeaxis-start_time)
    
    
    #ax0xmax=max(SBEsample.T)
    #ax0xmin=min(SBEsample.T)
    ax0xmax=26
    ax0xmin=22

    ax1xmax=max(SBEsample.C)
    ax1xmin=min(SBEsample.C)

    ax2xmax=max(SBEsample.P)
    ax2xmin=min(SBEsample.P)

    tmin1=min(EPSIsample.timeaxis-start_time)
    tmax1=max(EPSIsample.timeaxis-start_time)


    ax3xmax=max([max(EPSIsample.temp1),max(EPSIsample.temp2)])
    ax3xmin=min([min(EPSIsample.temp1),min(EPSIsample.temp2)])

    ax4xmax=max([max(EPSIsample.shear1),max(EPSIsample.shear2)])
    ax4xmin=min([min(EPSIsample.shear1),min(EPSIsample.shear2)])

    ax5xmax=max([max(EPSIsample.accelx),max(EPSIsample.accely),max(EPSIsample.accelz)])
    ax5xmin=min([min(EPSIsample.accelx),min(EPSIsample.accely),min(EPSIsample.accelz)])
    
    
    
    fig.axes[0].set_ylim([tmin1,tmax1])   
    fig.axes[0].set_xlim([ax0xmin,1.1*ax0xmax])   

    fig.axes[1].set_ylim([tmin1,tmax1])   
    fig.axes[1].set_xlim([ax1xmin,ax1xmax])   

    fig.axes[2].set_ylim([tmin1,tmax1])   
    fig.axes[2].set_xlim([ax2xmin,ax2xmax])   


    fig.axes[3].set_ylim([tmin1,tmax1])
    fig.axes[3].set_xlim([ax3xmin,ax3xmax])   
    fig.axes[4].set_ylim([tmin1,tmax1])   
    fig.axes[4].set_xlim([ax4xmin,ax4xmax])   
    fig.axes[5].set_ylim([tmin1,tmax1])   
    fig.axes[5].set_xlim([ax5xmin,ax5xmax])   

    return line0,line1,line2,lineF1,lineF2,lineS1,lineS2,lineA1,lineA2,lineA3


    

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
  


fig=init_figure(6)
line0,line1,line2,lineF1,lineF2,lineS1,lineS2,lineA1,lineA2,lineA3=init_axes(fig)
#ani = animation.FuncAnimation(fig, animate,frames=100,interval=10, blit=True)
ani = animation.FuncAnimation(fig, animate,interval=1)
plt.show()
   




#!/usr/bin/env python

import serial
import time
import sys
import glob

serport=glob.glob('/dev/tty.usbserial-*')

ser = serial.Serial(serport[0],460800)  # open serial port
print(ser.name)         # check which port was really used

#def open_datafile(filename='../data/MADRE2.1_' +time.strftime("%m%d%Y_%H%M%S")+ '.dat'):
def open_datafile(filename='../data/MADRE2.1.dat'):
    fid=open(filename,'wb+')
    return fid 



start_time = time.time()

State=0
compt=0
rollover=0
ser.flushInput()


if len(sys.argv)>2:
   filename=sys.argv[2]
   fid=open_datafile(filename)
else:
    fid=open_datafile()

Aux1WordLength = 0
#ADCWordlength  = 6
ADCWordlength  = 3
number_of_sensor  = 7
EpsisampleWordLength= ADCWordlength*number_of_sensor
epsisample_per_block  = 160 

EPSIWordlength = ADCWordlength *  number_of_sensor * epsisample_per_block

State=0
count=0

while True:
    if (State==0):
#      while(ser.in_waiting>0):
#          ser.readline()
#          print('wait')
          
          time.sleep(.00001)
          line=ser.readline()
          if line[:6]==b'$MADRE':
             #State=1
             len_header    = len(line)
             EpsiStamp     = int(line[6:6+8],16)
             TimeStamp     = int(line[15:15+8],16)
             Voltage       = int(line[24:24+8],16)
             Checksum_aux1 = int(line[33:33+8],16)
             Checksum_aux2 = int(line[42:42+8],16)
             Checksum_map  = int(line[51:51+8],16)
                      
             newHeader=ser.read(5)
             if newHeader==b'$AUX1':
                 print('AUXheader='+newHeader)
                 newHeader=ser.read(5)
                 
             if newHeader==b'$EPSI':
                 block=ser.read(EPSIWordlength)
                 endblock=ser.read(2)
                 if endblock==b'\r\n':
                     count+=1
                     if count==5:
                         State=1
                         print('Lock up')
                         print('Start recording')
                     else:
                         print(5-count)    
    
    if (State==1):

         line=ser.readline()
         print(line) 
         if (line[:6]==b'$MADRE')==False:
            State=0
            count=0
         else:
        
             fid.write(line)
             fid.flush() 
             EpsiStamp     = int(line[6:6+8],16)
             #TimeStamp     = int(line[15:15+8],16)
             #Voltage       = int(line[24:24+8],16)
             #Checksum_aux1 = int(line[33:33+8],16)
             #Checksum_aux2 = int(line[42:42+8],16)
             #Checksum_map  = int(line[51:51+8],16)
             newHeader=ser.read(5)
             if Checksum_aux1>0:
                 print('AUXheader='+newHeader)
                 newHeader=ser.read(5)
             if Checksum_aux2>0:
                 print('AUXheader='+newHeader)
                 newHeader=ser.read(5)
                     
             block=ser.read(EPSIWordlength) # 
             epsisamples=[ block[i*EpsisampleWordLength:(i+1)*EpsisampleWordLength] \
                           for i in range(epsisample_per_block) ]
             endblock=ser.read(2) # 2 is for /r/n
        #     fid.write(line)
             if(ADCWordlength==3): 
                [fid.write(str.encode(samples.hex() + '\r\n')) for samples in epsisamples]
             else:
                [fid.write(samples + b'\r\n') for samples in epsisamples]
                
             fid.flush() 

    


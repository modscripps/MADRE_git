#!/usr/bin/env python

import serial
import time
import sys

ser = serial.Serial('/dev/tty.usbserial-FTYVXOQA',460800)  # open serial port
#ser = serial.Serial('/dev/tty.usbserial-FTYVXOI5',460800)  # open serial port
#ser = serial.Serial('/dev/tty.usbserial-FTYVXOH6',460800)  # open serial port
print(ser.name)         # check which port was really used


def open_datafile(filename='SBE_EPSI_hex.dat'):
    fid=open('../data/' + filename,'wb+')
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

while (State==0):
      while(ser.in_waiting>0):
          ser.readline()
          print('wait')
          
      print('Empty buffer')
      time.sleep(.00001)
      line=ser.readline()
      if len(line)==30:
         byte_per_block = int(line[-10:-2],16)
         SBE_samples    = int(line[2:10],16)
         TX_block_size  = int(line[12:19],16)
         block=ser.read(TX_block_size+24) # 24 is for the length of the SBE49 sample
         line1=ser.readline()
         print("block\r\n")
         print(block)
         print("line1\r\n")
         print(line1)
         if (line1[:2]==b'##'): 
            State=1
            fid.write(line1)
            time_stamp = int(line[2:10],16)
            offset_t   = time_stamp
            timestamp=start_time+(time_stamp-offset_t) 
            str_timestamp=str(timestamp)+'\r\n'
            b=bytearray()
            b.extend(map(ord,str_timestamp))
            fid.write(b)
            timestamp0=timestamp

while (State==1):
      while(ser.in_waiting<100):
        test=1	
        strcoucou='coucou\n'
        b=bytearray()
        b.extend(map(ord,strcoucou))
        time.sleep(.00001)


      line=ser.readline()
      fid.write(line)
      if len(line)==30:
          time_stamp     = int(line[2:10],16)
          timestamp=start_time+(time_stamp-offset_t)+rollover*0xffffffff 
          if timestamp<timestamp0:
              rollover+=1
          timestamp0=timestamp
          timestamp=start_time+(time_stamp-offset_t)+rollover*0xffffffff 
          str_timestamp=str(timestamp)+'\r\n'
          b=bytearray()
          b.extend(map(ord,str_timestamp))
          fid.write(b)

      fid.flush() 
      print(line) 
      compt+=1
      if compt<0:
         State=0
         ser.close()


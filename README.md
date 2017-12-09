# MADRE_git

Pre-requisit: 
  - download anaconda3 for python
  - in a terminal "conda install pyserial"


1/ Connect MADRE and your laptop with an FTDI serial device. The python reader will an issue if you have more than one FTDI device. It is possible to set the name of the device in read_MADRE2.1.py. TODO: give the user a choice if more than 1 device is available

2/ in a terminal, go into READ_MADRE: python read_MADRE2.1.py

3/ in  another terminal: python plot_MADRE2.1_realtime_maptest.py. There is a bug in this routine, soon to be corrected  

4/ in  another terminal: python plot_MADRE2.1_realtime_spectrum_maptest.py. Not implemented yet. As soon as the bug above is fixed this routine will be available


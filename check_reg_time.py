# PRINT THE FIRST 100 TIME STAMPS TO ENSURE THAT THE DATA ARE ON A REGULAR GRID
# I had to do this because DP3 phaseshift does weird things by adding in new rows with irregularly spaced time axis
# NOTE: CHANGE nant BELOW TO YOUR CASE (ELSE THE CODE WILL GIVE WRONG RESULTS)
# Originally written to check uGMRT data that had been passes through DP3 phaseshift
# Incidentally, I learnt that DP3 addes in new timeslots with weird spacings
# Maybe because it gets confused with irregular timespacing due to interleaved calibration scans
#
import casacore.tables as tab
import numpy as np
import sys
#
nant=30
ntime = 100
ms = sys.argv[1]
t = tab.table(ms)
nbase=int(nant*(nant-1)/2)
#
a1 = np.zeros((ntime,))
a2 = np.zeros((ntime,))
time = np.zeros((ntime,))
#
for i in range(ntime):
    time[i]=t[i*nbase]["TIME"]
    a1[i]=t[i*nbase]["ANTENNA1"]
    a2[i]=t[i*nbase]["ANTENNA2"]
# Print antenna 1 and antenna 2 to make sure we have the current nbase
print ("FIRST %d TIMESLOTS"%ntime)
print ("ANTENNA1:  ")
print (a1)

print ("ANTENNA2:  ")
print (a2)
print ("DELTA_TIME ")
print (np.diff(time))

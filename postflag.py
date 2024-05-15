# Code for an additional round of outlier flagging after the LoTSS pipeline. 
# This flagging is sometimes necessary because of large outliers in the data
# These outliers are rare and as such do not affect the final 8-hr synthesis image significantly
# But they can create problems in snapshot images made for variable and transient source detection
# At the same time, we do not want to flag bright transients
# So the flagging here has two steps (Note: It is based on Stokes V component of vis)
# Step - 1: Reject all data that have abs(StokesV vis)> 100 Jy as these transients would have been detected in the
# standard 8-hr images as Stokes V sources.
# Step - 2: Flag all data with abs(Stokes V vis)>20*MAD for each timestep
# I.E the MAD is taken over all baselines and spectral channels.
# This ensures that genuine bright browd band transients are not rejected but RFI is as it is often narrowband
# and often afflicts a restructed set of baselines
# There is also a uv filter hardcoded between 200m and 90km which is the same filter that is used in the snapshot imaging
#
#
import numpy as np
import pyrap.tables as pt
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import sys
import configparser
import os
#
#
#
def pflag(msname):
    #
    fs = open("postflag.log","a")
    fs.write("--- Working on %s ---\n"%msname) 
    print ("--- Working on MS %s ---"%msname)
    #
    t = pt.table(msname,readonly=False)
    nrows = t.nrows()
    # READ IN THE FLAG AND ANTENNA TABLES
    f = t.getcol("FLAG")
    np.savez("%s.oldflags.npz"%msname[:-9])
    a21 = t.getcol("ANTENNA2") - t.getcol("ANTENNA1")
    amin,amax = np.amin(a21),np.amax(a21)
    #
    # I DONT TRUST THE ANTENNA TABLES AS IT SOMETIMES HAS MORE ANTENNAS THAN THE MS ITSELF
    # ALSO SOMETIMES THE MS DOES NOT HAVE THE AUTOCORRELATIONS AND SOMETIMES IT DOES
    # SO I NEED A WAY TO FIND OUT ALL OF THE ABOVE WITH THE MS ITSELF
    #
    nant = np.amax(a21)+1
    if amin==0:             # DOES IT HAVE AUTOCORRELATIONS?
        nbase = int(nant*(nant-1)/2+nant)
    else:                   # NO AUTOCORRELATIONS
        nbase = int(nant*(nant-1)/2)
    ntime = int(nrows/nbase)
    #
    #print ("Nant = %d"%nant)
    #print ("Nbase = %d"%nbase)
    #print ("Nrows = %d"%nrows)
    #
    #
    # READ THE DATA FLAGS AND UVW VALUES
    #
    d = t.getcol("DATA")                            # VISIBILITY DATA (NTIME X NCHAN X NPOL)
    f = t.getcol("FLAG")                            # FLAGS (NTIME X NCHAN X NPOL)
    uvw = t.getcol("UVW")
    blen = (np.sum(uvw**2,axis=1))**0.5             # BASELINE LENGTHS IN METRES (NTIMEX1) 
    #
    # FLAG ALL STOKES V OVER 100 JY
    #
    vraw = np.absolute(-1j*d[:,:,1]+1j*d[:,:,2])    # STOKES V ABS. VISIBILITY
    init_flagfrac=np.mean(f)                        # INITIAL FLAG FRACTION
    I=vraw>100                                      # CLIP AT 100 JY
    for pol in range(4):
        f[I,pol]=True
    #
    clip_flagfrac = np.mean(f)                      # FLAG FRACTION AFTER CLIPPING
    fs.write ("Initial flagfrac = %f\n"%init_flagfrac)
    fs.write ("Post clipping flagfrac = %f\n"%clip_flagfrac)
    #
    # FLAG  MASK TO ONLY TAKE BASELINES BETWEEN 200M and 90KM (DONT APPLY TO THE DATA!)
    #
    blen_flag = np.zeros(d[:,:,1].shape,dtype=bool)
    for i in range(len(d[0,:,0])):
        blen_flag[:,i] = np.logical_or(blen<200.0,blen>90000)
    #
    fraw = np.logical_or(f[:,:,1],f[:,:,2])
    #
    #
    for i in range(ntime):                      # FOR EACH INTEGRATION
        vt = vraw[i*nbase:(i+1)*nbase,:]       # VISIBILITIES FOR THIS TIMESLOT (NBASE X NCHAN)
        ft = fraw[i*nbase:(i+1)*nbase,:]
        vt_noflag = vt[~ft]
        # FLAG ENTIRE TIMESLOT IF LESS THAN 10% OF THE DATA ARE UNFLAGGED IN THAT SLOT
        if np.mean(f[i*nbase:(i+1)*nbase,:,:])>0.9:
            f[i*nbase:(i+1)*nbase,:,:]=True
        else:
            mean = np.mean(vt_noflag)              # MEDIAN OF ABS VIS
            std = np.median(np.absolute(vt_noflag-mean)) # MEDIAN ABSOLUTE DEVIATION OF ABS VIS
            newflag = np.logical_or((vt-mean)>20*std,(vt-mean)<-20*std) # NEW FLAGS (NBASE X NCHAN)
            for pol in range(4):                # OR GATE THE NEW FLAGS INTO THE OLD FLAGS
                f[i*nbase:(i+1)*nbase,:,pol] |= newflag
    #
    final_flagfrac = np.mean(f)                 # FINAL FLAGGED FRACTRION
    fs.write ("Final flagfrac = %f\n"%final_flagfrac)
    #
    # IF FINAL FLAGFRAC EXCEEDS 85% THEN FLAG THE WHOLE DATASET
    #
    if final_flagfrac>0.85:
        f[:,:,:]=True
    np.savez("%s.newflags.npz"%msname[:-9])
    #
    t.putcol("FLAG",f)
    t.close()
    fs.write("---       END       ---\n")
    fs.close()
#
#
# --- MAIN BEGINS HERE ---
#
filelist=os.listdir()
for name in filelist:
	if name.endswith('.ms.archive'):
		pflag(name)
#

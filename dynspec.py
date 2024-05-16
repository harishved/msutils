import numpy as np
import casacore.tables as tab
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
#
#
def rotmat(RA0,DEC0,RA1,DEC1):
   # Calculate rotation matricies for UVW coordinate recalculation
   # during phaseshifting visibilities
   # Done using an intermediate coordinate system XYZ
   # Let XYZ be a system where Z points to NCP, X to RA=6h and Y to 12h 
   # In the UVW convention we have:
   # W points to a source, UV are perpendicular to W and V,W, and NCP
   # lie in a plane and U points in direction of increating RA
   # Consider the rotation XYZ --> UVW
   # This can be done by:
   # First rotate around Z axis by an angle RA
   # Then rotate along new X axis by angle pi/2-DEC
   # Below, M1 is the matric to get the first rotation and M2 for second
   # Then to go from XYZ--> UVW left multiply by M2M1
   # to go from UVW--> XYZ left multiply by (M2M1)^-1 = (M2M1)^T
   #
   # Inputs: 
   # RA0,DEC0 --> existing RA and DEC in radian
   # This should be typically the value in the PHASE_DIR of the FIELD subtable
   # RA1,DEC1 --> New RA and DEC in radian
   # This is the direction in which ms should be phase rotated
   #
   cosRA0 = np.cos(RA0)
   sinRA0 = np.sin(RA0)
   cosDE0 = np.cos(DEC0)
   sinDE0 = np.sin(DEC0)
   cosRA1 = np.cos(RA1)
   sinRA1 = np.sin(RA1)
   cosDE1 = np.cos(DEC1)
   sinDE1 = np.sin(DEC1)

   M1_0 = np.array([ [ cosRA0, sinRA0, 0],\
                     [-sinRA0, cosRA0, 0],\
                     [0,       0,      1]])

   M1_1 = np.array([ [ cosRA1, sinRA1, 0],\
                     [-sinRA1, cosRA1, 0],\
                     [0,       0,      1]])

   M2_0 = np.array([ [1,      0,    0   ],\
                     [0,  sinDE0, cosDE0],\
                     [0, -cosDE0, sinDE0]])

   M2_1 = np.array([ [1,      0,    0   ],\
                     [0,  sinDE1, cosDE1],\
                     [0, -cosDE1, sinDE1]])
   
   r1 = np.transpose(np.matmul(M2_0,M1_0))
   r2 = np.matmul(M2_1,M1_1)
   rotmat = np.matmul(r2,r1)
   return rotmat
#
#
def shift_ms(msname,RA1,DEC1,colname=["DATA"],update_ms=False):
   # Phase shift the input measurement set in situ
   # msname = Input measurement set name
   # RA1, DEC1 = old phase centre
   # RA2, DEC2 = New phase centre
   # colname = list of column names that needs to be phase shifted
   # Ensure that RA and DEC are in radians
   #
   # Read in freq vector and ra dec of phase centre
   t = tab.table(msname+"/SPECTRAL_WINDOW")
   freq = t.getcol("CHAN_FREQ")
   t.close()
   t = tab.table(msname+"/FIELD")
   RA0,DEC0 = t.getcol("PHASE_DIR").flatten()
   print ("Old phase centre = %f, %f rad"%(RA0,DEC0))
   print ("New phase centre = %f, %f rad"%(RA1,DEC1))

   # Read in old uvw data
   if update_ms:
      t = tab.table(msname,readonly=False)
   else:
      t = tab.table(msname)

   old_uvw = t.getcol("UVW")

   # Calculate new uvw
   r = rotmat(RA0,DEC0,RA1,DEC1); 
   new_uvw = (np.matmul(r,old_uvw.T)).T

   # Calculate phase shift factor
   lam = 2.99792458e8/freq
   diff_w = np.reshape(new_uvw[:,2]-old_uvw[:,2],(old_uvw.shape[0],1))
   phasor = np.exp(-2*np.pi*1j*np.matmul( diff_w, 1.0/lam))

   # Write new data into ms
   if update_ms:
      t.putcol("UVW",new_uvw)

   result=[]
   
   for c in list(colname):
      print ("Shifting column %s"%c)
      d = t.getcol(c)
      dnew = np.zeros(d.shape,dtype=d.dtype)
      for i in range(d.shape[2]):
          dnew[:,:,i]=d[:,:,i]*phasor
      if update_ms:    
         t.putcol(c,dnew)
      else:
         result.append(dnew)

   t.close()
   
   if update_ms:
      return
   else:
      return result
#
#

def nant_nbase(msname):
    # Calculate the number of antennas and baselines in a ms
    # ASsume a max of 100 antennas possible to speed up
    t = tab.table(msname)
    a1 = t.getcol("ANTENNA1",0,int(100*99/2))
    a2 = t.getcol("ANTENNA2",0,int(100*99/2))
    a21 = a2-a1
    a21min = np.amin(a21)
    nant = np.amax(a21)+1
    if a21min==0:
        nbase = int(nant*(nant-1)/2+nant)
    else:
        nbase = int(nant*(nant-1)/2)
    return nant,nbase
#

def calc_weights(msname,colname,SEFD=300,update_ms=False):
   # Calculate the noise in the MS using time differencing and 
   # use that to update the WEIGHT_SPECTRUM and SIGMA_SPECTRUM
   # SIGMA_SPECTRUM is the rms noise per visibility = SEFD/(2*BW*t)^1/2
   # WEIGHT_SPECTRUM is just the weight (ideally given by 2*BW*t
   # msname = Input measurement set
   # colname = Column name to use to calculate noise values 

   nant,nbase = nant_nbase(msname)
   if update_ms:
      t = tab.table(msname,readonly=False)
   else:
      t = tab.table(msname)
   #
   # Check if SIGMA and WEIGHT columns exist; id not then add them
   nchan,npol = t[0][colname].shape
   allcols = t.colnames()
   for tp in ["SIGMA_SPECTRUM","WEIGHT_SPECTRUM"]:
      if tp not in allcols:
         print ("Adding column %s"%tp)
         dminfo = {'TYPE': 'TiledShapeStMan','SPEC': {'DEFAULTTILESHAPE': [npol, nchan, 128]}}
         dminfo["NAME"] = tp
         cd = t.getcoldesc(colname)
         cd["NAME"] = tp
         cd["valueType"]='float'
         t.addcols(tab.maketabdesc(tab.makecoldesc(tp,cd)),dminfo)

   alld = t.getcol(colname)
   allf = t.getcol("FLAG")

   ntime = int(t.nrows()/nbase)

   noise_std = np.zeros((nbase,nchan,npol))

   for i in range(nbase):
      f = allf[i::nbase,:,:] 
      d = np.ma.array(alld[i::nbase,:,:],mask=f)
      dd = np.ma.diff(d,axis=0)
      noise_std[i,:,:] = np.ma.std(dd[:,:,:],axis=0)/2**0.5
  
   noise_std[noise_std==0] = 1e10

   weights = (SEFD/noise_std)**2
   if update_ms:
      for i in range(ntime):
         t.putcol("SIGMA_SPECTRUM",noise_std,i*nbase,nbase)
         t.putcol("WEIGHT_SPECTRUM",weights,i*nbase,nbase)
      t.close()
      return

   t.close()
   return noise_std, weights
#
#

def calc_dynspec(msname,ra,dec,method="simple",colname="DATA"):
   #
   # method = simple assumes that the sigma and weight column have the correct entries
   # If they do not then run calc_weights with update_ms=True before
   # 
   # Shift the visibilities to target
   d = shift_ms(msname,RA1,DEC1,colname=["DATA"],update_ms=False)[0]
   
   t = tab.table(msname+"/SPECTRAL_WINDOW")
   fvec = t.getcol("CHAN_FREQ").flatten()

   nant,nbase = nant_nbase(msname) # Num and ant and blines in ms
   t = tab.table(msname)
   nchan,npol = t[0][colname].shape
   ntime = int(t.nrows()/nbase)
   tvec = t.getcol("TIME",0,t.nrows(),nbase).flatten()

   dspec = np.zeros((ntime,nchan,npol)) # Output dynamic spectrum
   
   w = t.getcol("WEIGHT_SPECTRUM")
   s = t.getcol("SIGMA_SPECTRUM")
   f = t.getcol("FLAG")

   t.close()

   dm = np.ma.array(d.real,mask=f) # Input visibility data (shifted)
   wm = np.ma.array(w,mask=f) # Weights
   sm = np.ma.array(s,mask=f) # Noise per visibility

   dspec_noise = np.zeros(dspec.shape) # Error on dynspec values per pixel

   for i in range(ntime):
      for j in range(nchan):
         for k in range(npol):
            dspec[i,j,k] = np.ma.sum(dm[i*nbase:i*nbase+nbase,j,k]*w[i*nbase:i*nbase+nbase,j,k])/np.ma.sum(w[i*nbase:i*nbase+nbase,j,k])
            dspec_noise[i,j,k] = (1.0/np.ma.sum(1/sm[i*nbase:i*nbase+nbase,j,k]**2))**0.5

   return dspec,dspec_noise,tvec,fvec

def plot_dynspec(d,n,t,f,opfname="dynspec.pdf"):
   # d = Dynamic spec data (ntime,nchan,npol)
   # n = noise per pixel (ntime,nchan,npol)
   # t = Time stamps in epoch
   # f = Freqcienies in Hz
   # 
   t0 = t[0]/24/3600.0 # Starting time in MJD days
   t-=t[0]
   t/=60.0 # Now t is in minutes

   f/=1e6 # Now f is in MHz

   d*=1e3
   n*=1e3 # Now data and noise are in mJy

   n[n==0]=np.nan
   w = 1/n**2
   dt = np.nansum(d*w,axis=1)/np.nansum(w,axis=1)
   df = np.nansum(d*w,axis=0)/np.nansum(w,axis=0)

   nt,nf,npol = d.shape

   with PdfPages(opfname) as pdf:
      for i in range(npol):
         med = np.nanmedian(d[:,:,i])
         mad = np.nanmedian(np.absolute(d[:,:,i]-med))
         fig, axes=plt.subplots(2,2, sharex='col', sharey='row', width_ratios=[4,1], height_ratios=[1,4])
         img=axes[1,0].imshow(d[:,:,i].T, aspect='auto',interpolation='none', cmap='cividis', vmin=med-10*mad, vmax=med+10*mad,extent=[t[0],t[-1],f[0],f[-1]], origin='lower')
         axes[1,0].set_xlabel('Time (min)')
         axes[1,0].set_ylabel('Frequency (MHz)')
         axes[0,0].plot(t,dt[:,i])
         axes[1,1].plot(df[:,i],f)

         plt.subplots_adjust(wspace=0, hspace=0)
         plt.colorbar(img, label='Flux density (mJy)',ax=axes.ravel().tolist())
         fig.delaxes(axes[0,1])

         pdf.savefig()
         plt.close()

#
# TESTING
msname = "temp.ms"
RA1 = 2.7076476017675555
DEC1 = 1.1493846140091937
#
# If the sigma_spectrum and weight_spectrum columns do not have 
# meaningful values, then we must first calculate them using 
# ime differencing
# This is because dynspec code expects to find the values here that will
# allow in inverse variance weighting
calc_weights(msname, colname="DATA",update_ms=True)
# Calculate the dynspec
d,n,t,f = calc_dynspec(msname,RA1,DEC1,method="simple",colname="DATA")
# plot the dynspace
plot_dynspec(d,n,t,f)


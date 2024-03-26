# Convert GMRT RR LL into format that DP3 likes (with 4 POL PRODUCTS)
# Should be able to work on any data though
# If using on linear feeds then make sure that the POLARIZATION table values are changed
# Right now circ. feeds have been hardcoded in the code
# Harish Vedantham
# 20 March 2024, ASTRON
# Need to change the following
# FLAG, DATA, WEIGHT, SIGMA, WEIGHT_SPECTRUM values (2,n) --> (4,n) arrays
# and their column descriptors
# NUM_CORR, CORR_TYPE and CORR_PRODUCT values in the POLARIZATION TABLE
#
#
import sys
import casacore.tables as tab
import numpy as np

msin=sys.argv[1]

t=tab.table(msin,readonly=False)
'''
for col in ["FLAG","DATA","WEIGHT_SPECTRUM"]:
    ip_desc = t.getcoldesc(col)
    print (col, ip_desc)

for col in ["WEIGHT","SIGMA"]:
    ip_desc = t.getcoldesc(col)
    print (col, ip_desc)
'''
nchan=t[0]["FLAG"].shape[0]
dminfo = {'TYPE': 'TiledShapeStMan','SPEC': {'DEFAULTTILESHAPE': [4, nchan, 128]}}

for col in ["FLAG","DATA","WEIGHT_SPECTRUM"]:
   cd = t.getcoldesc(col)
   cd["shape"][1]=4
   cd["NAME"]=col
   dminfo["NAME"]=col
   t.renamecol(col,col+"1")
   t.addcols(tab.maketabdesc(tab.makecoldesc(col,cd)),dminfo)


dminfo = {'TYPE': 'TiledShapeStMan','SPEC': {'DEFAULTTILESHAPE': [4, 65520]}}

for col in ["WEIGHT","SIGMA"]:
   cd = t.getcoldesc(col)
   cd["shape"][0]=4
   cd["NAME"]=col
   dminfo["NAME"]=col
   t.renamecol(col,col+"1")
   t.addcols(tab.maketabdesc(tab.makecoldesc(col,cd)),dminfo)

nrows=t.nrows()
nbase=128
nblock=int(np.floor(float(nrows)/float(nbase)))
for i in range(nblock+1):
    id_start=min(nrows-1,i*nbase)
    id_end=min(nrows,i*nbase+nbase)
    print ("Writing rows %d - %d"%(id_start,id_end))
    for col in ["FLAG","DATA","WEIGHT_SPECTRUM"]:
        dat=t.getcol(col+"1",id_start,id_end-id_start)
        datnew = np.zeros((id_end-id_start,dat.shape[1],4),dtype=type(dat[0,0,0]))
        datnew[:,:,0]=np.fliplr(dat[:,:,0])
        datnew[:,:,3]=np.fliplr(dat[:,:,1])
        t.putcol(col,datnew,id_start,id_end-id_start)
        
for col in ["FLAG","DATA","WEIGHT_SPECTRUM","WEIGHT","SIGMA"]:
    t.removecols([col+"1"])
t.flush()
t.close()

t = tab.table(sys.argv[1]+"/POLARIZATION",readonly=False)
t.putcol("CORR_TYPE",np.array([[5,6,7,8]],dtype=np.int32))
t.putcol("CORR_PRODUCT",np.array([[[0,0],[0,1],[1,0],[1,1]]],dtype=np.int32))
t.putcol("NUM_CORR",np.array([[4]],dtype=np.int32))
t.flush()
t.close()

t = tab.table(sys.argv[1]+"/SPECTRAL_WINDOW",readonly=False)
d = t.getcol("CHAN_FREQ")
t.putcol("CHAN_FREQ",np.fliplr(d))
d = t.getcol("CHAN_WIDTH")
t.putcol("CHAN_WIDTH",np.absolute(d))
t.flush()
t.close()
#
tab.addImagingColumns(sys.argv[1])
#t = tab.table(sys.argv[1]+"/ANTENNA",readonly=False)
#t.removerows([31,30])
#t.flush()
#t.close()

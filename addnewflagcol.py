# Copy FLAGS to an OLD_FLAG column and make a new FLAG column
import sys
import casacore.tables as tab
import numpy as np
#
msin=sys.argv[1]
t=tab.table(msin,readonly=False)
npol = t[0]["FLAG"].shape[1]
nchan = t[0]["FLAG"].shape[0]
old_flags = t.getcol("FLAG") # Existing flag data
coldesc = t.getcoldesc("FLAG") # Column descriptor
dminfo = {'TYPE': 'TiledShapeStMan','SPEC': {'DEFAULTTILESHAPE': [npol, nchan, 128]}} # Table info (minimalist)
t.renamecol("FLAG","OLD_FLAG") # Rename existing flag column
t.addcols(tab.maketabdesc(tab.makesoldesc("FLAG",coldesc)),dminfo) # Add new flag column
t.putcol("FLAG",old_flags) # copy old flags into new flag column
t.flush()
t.close()
#

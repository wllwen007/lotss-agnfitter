import sys
import os
from astropy.table import Table, Column
import astropy.units as au
import numpy as np
import utils.plot_util as pp

#field = 'Bootes'
#field = 'LH'
#field = 'EN1'
field = sys.argv[1]

if field=='all':
    fields = ['Bootes','EN1','LH']
else:
    fields = [field]


for field in fields:
    
    rcat = Table.read(field+'_agnfitter_v0.6.fits')  
    ocat = Table.read(field+'_opt_agnfitter.fits')
    
    nr = len(rcat)
    no = len(ocat)
    
    print('{field}: radio - {nr}\n{field}: optical - {no}'.format(field=field, nr=nr, no=no))

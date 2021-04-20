import glob
import linecache
import os
from astropy.table import Table, Column
import numpy as np
#from get_running import get_running

clobber=True

def get_range(vals):
    if len(vals) == 0:
        return ''
    vals.sort()
    srange = ''
    s = vals[0]
    srange += str(s)
    for i in range(1,len(vals)):
        if vals[i]-vals[i-1] == 1:
            if srange[-1] != '-': srange+='-'
        else:
            if srange[-1] == '-': 
                srange+=str(vals[i-1])
            srange+=','+str(vals[i])
    if srange[-1] == '-': srange+=str(vals[-1])
    return srange

#lookup = Table.read('array_id_lookup.fits')
#zsamples_cat = 'Bootes_v2014_agnfitter_zsamples_v2.fits'
#zcat = Table.read(zsamples_cat)
cat = Table.read('data/bootes_cat_agnfitter_in.fits')
cat.add_column(Column(np.arange(len(cat)), 'array_id'))


##zsel = np.array([ii in cat['id'] for ii in zcat['id']])
#zsel = np.array([ii in lookup['id'] for ii in zcat['id']])
#zcat = zcat[zsel]

#single = (zcat['zselect'] == 0.) | (zcat['z_spec'] > 0)
#mult = ~single
sid = cat['ID']
#mid = zcat['id'][mult]



path = '/beegfs/lofar/wwilliams/lofar_surveys/deep/agnfitter/agnfitter_out_photz/agnfitter_out/'


outfiles0 = [path+'{id:d}/parameter_outvalues_{id:d}.txt'.format(id=t) for t in sid ] 



## updated - add GA 1-30 mu
'''
# Output for source 1783740
Rows are: 2.5, 16, 50, 84, 97.5 percentiles # and max like
-----------------------------------------------------
tau age Nh irlum SB BB GA TO EBVbbb EBVgal log Mstar SFR_opt LIR(8-1000) Lbb(0.1-1) Lbbdered(0.1-1) Lga(0.1-1) Ltor(1-30) Lsb(1-30) SFR_IR -ln_like
'''
# 16th and 84th percentiles correspond to -/+1sigma (at least for gauss dist)

#cmd = "sed '6q;d' OUTPUT/1783740/parameter_outvalues_1783740.txt"
linehead = linecache.getline(outfiles0[0], 4).strip() + ' '
linehead = linehead.replace('log Mstar','log_Mstar')
linehead = linehead.replace('-ln_like','neg_ln_like')
header = '# id '+linehead.replace(' ','_med ')+linehead.replace(' ','_errlow ')+linehead.replace(' ','_errhigh ') #+linehead.replace(' ','_ml ')

lines0 = [header+'\n']

missing0 = []
for fi,filen in enumerate(outfiles0):
    #print fi
    if os.path.exists(filen):
    
        #idn = filen.split('/')[-1].replace('parameter_outvalues_','').replace('.txt','')
        idn = filen.split('/')[-2]

        medianline = linecache.getline(filen, 7).strip()
        errlowline = linecache.getline(filen, 6).strip()
        errhighline = linecache.getline(filen, 8).strip()
        mlline = linecache.getline(filen, 10).strip()
        
        newline = idn+' '+medianline+' '+errlowline+' '+errhighline+' '+mlline
        
        lines0.append(newline+'\n')
    else:
        missing0.append(filen.split('/')[-2])
        
#print '0' , ','.join(missing0)

outfile = path+'/agnfitter_out.txt'

if os.path.exists(outfile.replace('.txt','.fits')): 
    if clobber:
        os.remove(outfile.replace('.txt','.fits'))
with open(outfile,'w') as f:

    f.writelines(lines0)


# convert to fits
os.system('stilts tcopy in={inf:s} ifmt=ascii ofmt=fits out={out:s}'.format(inf=outfile,out=outfile.replace('.txt','.fits')))

missingid0 = [cat['array_id'][cat['ID']==int(m)][0] for m in missing0]
#missingid0 = missing0
print 'total rows missing: ',len(missingid0)
print 'missing array_ids: ' , get_range(missingid0)
print

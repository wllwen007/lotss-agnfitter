
import glob
import os
from astropy.table import Table
import sys
import numpy as np

field=sys.argv[1]
print 'field is',field


if field=='EN1':
    mask='/beegfs/lofar/deepfields/ELAIS_N1_optical/radio_optical_overlap_masks/image_full_ampphase_di_m.NS_shift.int.facetRestored-scaled.pybdsm.rmsd_I_spmask.fits'
    optcat='/beegfs/lofar/deepfields/ELAIS_N1_optical/catalogues/correct_merging/add_uncat/EN1_MASTER_opt_spitzer_merged_cedit_apcorr_adduncat_lite.fits'
    src='/beegfs/lofar/deepfields/science_ready_catalogs/EN1_opt_spitzer_merged_vac_opt3as_irac4as_all_hpx_masses_public.fits'
elif field=='Bootes':
    mask='/beegfs/lofar/deepfields/Bootes_merged_optical/radio_optical_overlap_masks/image_full_ampphase_di_m.NS_shift.blanked.scaled.rms_spmask.fits'
    optcat='/beegfs/lofar/deepfields/Bootes_merged_optical/add_uncat/Bootes_MASTER_opt_spitzer_merged_adduncat_lite.fits'
    src='/beegfs/lofar/deepfields/science_ready_catalogs/Bootes_opt_spitzer_merged_vac_opt3as_irac4as_all_hpx_masses_public.fits'
elif field=='LH':
    mask='/beegfs/lofar/deepfields/Lockman_edited_cats/radio_optical_overlap_masks/image_full_ampphase_di_m.NS_shift.blanked.scaled.rms_spmask.fits'
    optcat='/beegfs/lofar/deepfields/Lockman_edited_cats/optical/add_uncat/LH_MASTER_opt_spitzer_merged_cedit_apcorr_adduncat_lite.fits'
    src='/beegfs/lofar/deepfields/science_ready_catalogs/LH_opt_spitzer_merged_vac_opt3as_irac4as_all_hpx_masses_public.fits'
else:
    raise RuntimeError('Field not supported!')
    
#g=sorted(glob.glob('sources-v*.fits'))

#infile=g[-1]
#print 'Processing',infile
#outfile=infile.replace('sources','filtered')

#filter_table(infile,mask,outfile)

#flagname=outfile.replace('filtered','flagged')
#print 'Applying final flags'
#final_flag(field,outfile,flagname)

#mergeout=flagname.replace('flagged','merged')

#os.system('/soft/topcat/stilts tskymatch2 ra1=optRA dec1=optDec ra2=ALPHA_J2000 dec2=DELTA_J2000 error=1.5 join=all1 in1=%s in2=%s out=%s' % (flagname,optcat,mergeout))

#mergeout2=mergeout.replace('merged','merged_src')

#if field=='bootes':
#    print 'Fixing ID column!'
#    t=Table.read(mergeout)
#    t['NUMBER']=t['ID'].astype('int')
#    t['NUMBER'].mask=np.isnan(t['ID'])
#    t.write(mergeout,overwrite=True)

mergeout = field+'_full_optcat.fits'

os.system('/soft/topcat/stilts tmatch2  join=all1 values1=ID values2=ID matcher=exact in1=%s in2=%s out=%s' % (optcat,src,mergeout))

t=Table.read(mergeout)
finalname = field+'_full_optcat.fits'

print 'Remove unnecessary columns'

cols=['RA_2','DEC_2','FLAG_OVERLAP_2','FLAG_CLEAN_2','id', 'ID_OPTICAL_2', 'ID_SPITZER_2','ID_2','EBV_2','FLAG_DEEP_2']
for fcol in t.colnames:
    if fcol.endswith("_fluxerr"):
        cols.append(fcol)
        cols.append(fcol[:-3])

for c in cols:
    if c in t.colnames:
        print c,
        sys.stdout.flush()
        t.remove_columns(c)

print
print 'Rename columns'
#t['Separation_1'].name='Separation'

for column in ['ID','RA','DEC','FLAG_OVERLAP','flag_clean', 'ID_OPTICAL', 'ID_SPITZER','EBV','FLAG_DEEP']:

    oldcol=column+'_1'
    if oldcol in t.colnames:
        t[oldcol].name=column
        print oldcol,'->',column,':',

print
print 'Remove whitespace padding:',
sys.stdout.flush()
stringcols=[n for (n,ty) in t.dtype.descr if 'S' in ty]
for c in stringcols:
    print c,
    sys.stdout.flush()
    t[c]=[s.rstrip() for s in t[c]]
print

print 'Remove -99s:',
dblcols=[n for (n,ty) in t.dtype.descr if ('f8' in ty or 'f4' in ty)]
for c in dblcols:
    print c,
    sys.stdout.flush()
    t[c]=np.where(t[c]==-99,np.nan,t[c])
print


print 'Sorting'
t.sort('RA')

print 'Writing to disk'
t.write(finalname,overwrite=True)

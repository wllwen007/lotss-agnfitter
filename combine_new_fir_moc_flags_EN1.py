from astropy.table import Table, Column
import numpy as np
import os

### new XID+ fluxes

tnew = Table.read('/beegfs/lofar/deepfields/ELAIS_N1_FIR_prelim/XID+_lofar_ELAIS-N1_v0.5_20200113.fits')
t = Table.read('/beegfs/lofar/deepfields/data_release/en1/final_cross_match_catalogue-v0.5.fits')


# add flag

#badmips = (tnew['F_MIPS_24']==1e20)
### add a flag for these sources...
#newflag = np.ones(len(tnew), dtype=int)
#newflag[badmips] = 2
#tnew.add_column(Column(name='newXIDp', data=newflag))


newflag = np.zeros(len(t), dtype=int)
t.add_column(Column(name='newXIDp', data=newflag))
t.add_column(Column(name='XID+_rerun_mips', data=newflag))
t.add_column(Column(name='XID+_rerun_pacs', data=newflag))
t.add_column(Column(name='XID+_rerun_SPIRE', data=newflag))


# annoying!!!
t.rename_column('flag_spire_250','flag_SPIRE_250')
t.rename_column('flag_spire_350','flag_SPIRE_350')
t.rename_column('flag_spire_500','flag_SPIRE_500')

##

colsnew = tnew.colnames
cols = t.colnames


#don't replace 
#['Source_Name',  'RA_1',  'Dec_1','RA_2',  'Dec_2',  'newXIDp']


mipscols = ['F_MIPS_24',  'FErr_MIPS_24_u',  'FErr_MIPS_24_l',  'Bkg_MIPS_24',  'Sig_conf_MIPS_24',  'Rhat_MIPS_24',  'n_eff_MIPS_24',  'Pval_res_24',  'XID+_rerun_mips']

pacscols = ['F_PACS_100',  'FErr_PACS_100_u',  'FErr_PACS_100_l',  'F_PACS_160',  'FErr_PACS_160_u',  'FErr_PACS_160_l',  'Bkg_PACS_100',  'Bkg_PACS_160',  'Sig_conf_PACS_100',  'Sig_conf_PACS_160',  'Rhat_PACS_100',  'Rhat_PACS_160',  'n_eff_PACS_100',  'n_eff_PACS_160',  'Pval_res_100',  'Pval_res_160',  'flag_PACS_100',  'flag_PACS_160',  'XID+_rerun_pacs']

spirecols = ['F_SPIRE_250',  'FErr_SPIRE_250_u',  'FErr_SPIRE_250_l',  'F_SPIRE_350',  'FErr_SPIRE_350_u',  'FErr_SPIRE_350_l',  'F_SPIRE_500',  'FErr_SPIRE_500_u',  'FErr_SPIRE_500_l',  'Bkg_SPIRE_250',  'Bkg_SPIRE_350',  'Bkg_SPIRE_500',  'Sig_conf_SPIRE_250',  'Sig_conf_SPIRE_350',  'Sig_conf_SPIRE_500',  'Rhat_SPIRE_250',  'Rhat_SPIRE_350',  'Rhat_SPIRE_500',  'n_eff_SPIRE_250',  'n_eff_SPIRE_500',  'n_eff_SPIRE_350',  'Pval_res_250',  'Pval_res_350',  'Pval_res_500',  'flag_SPIRE_250',  'flag_SPIRE_350',  'flag_SPIRE_500',  'XID+_rerun_SPIRE']

print 'replacing mips'
cc = mipscols[0]
replace = np.isnan(t[cc])
t['newXIDp'][replace] += 1
minds = []
for ts in t[replace]['Source_Name']:
    minds.append(np.where(ts == tnew['Source_Name'])[0][0])
for i in range(len(mipscols)):
    cc = mipscols[i]
    t[cc][replace] = tnew[cc][minds]
    
print 'replacing pacs'
cc = pacscols[0]
replace = np.isnan(t[cc])
t['newXIDp'][replace] += 2
minds = []
for ts in t[replace]['Source_Name']:
    minds.append(np.where(ts == tnew['Source_Name'])[0][0])
for i in range(len(pacscols)):
    cc = pacscols[i]
    t[cc][replace] = tnew[cc][minds]
    
print 'replacing spire'
cc = spirecols[0]
replace = np.isnan(t[cc])
t['newXIDp'][replace] += 4
minds = []
for ts in t[replace]['Source_Name']:
    minds.append(np.where(ts == tnew['Source_Name'])[0][0])
for i in range(len(spirecols)):
    cc = spirecols[i]
    t[cc][replace] = tnew[cc][minds]
        
    

## Ian added RA, Dec cols ??? some, and some F24 flux vales are 1e20 ... manually hack these
#copycols = []
#for cc in colsnew:
    #if 'Source_Name' in cc: continue
    #if 'RA' in cc: continue    ## don't copy these...
    #if 'DEC' in cc: continue    ## don't copy these...
    #if 'Dec' in cc: continue    ## don't copy these...
    #if cc in cols:
        #copycols.append(cc)

#for ti in tnew:
    #mind = np.where(t['Source_Name'] == ti['Source_Name'])[0][0]
    #for cc in copycols:
        #if np.isnan(t[mind][cc]) or ('flag' in cc) or ('newXIDp'):
            #t[mind][cc] = ti[cc]
        #else:
            #print 'value already exists'
            #sys.exit(1)
            ## this should not happen



## for the sources in tnew - replace the columns
#minds = []
#for ti in tnew:
    #mind = np.where(t['Source_Name'] == ti['Source_Name'])[0][0]
    #minds.append(mind)
    #print ti['Source_Name'],mind
#for cc in copycols:
    #if cc not in ['newXIDp']:
        #tnew[cc][tnew[cc]==1e20] = np.nan   # fix the 1e20 aargh
    #t[cc][minds] = tnew[cc]


# save
t.write('/beegfs/lofar/wwilliams/lofar_surveys/deep/agnfitter/EN1_final_cross_match_catalogue-v0.5_mergenew.fits',overwrite=True)


# MOC flag


#t['ID'][t['ID'].mask] =-1
#tmoc['ID'][tmoc['ID'].mask] =-2 ## there shouldn't be any in here...
#tt = join(t, tmoc, keys=['ID','ALPHA_J2000','DELTA_J2000'])
#tt['ID'][tt['ID']==-1] = np.nan


# for the sources in tmoc - add the column
#minds = []
#mindst = []
#for ti in t:
    #if (ti['ID'] in tmoc['ID']):
        #mind = np.where(ti['ID'] == tmoc['ID'])[0][0]
        #minds.append(mind)
        #mindst.append(True)
        ##print ti['Source_Name'],mind
    #else:
        #mindst.append(False)
    

ft = '/beegfs/lofar/wwilliams/lofar_surveys/deep/agnfitter/EN1_final_cross_match_catalogue-v0.5_mergenew.fits'

# MOC flags
ftmoc = '/beegfs/lofar/deepfields/ELAIS_N1_optical/catalogues/correct_merging/add_uncat/EN1_MASTER_opt_spitzer_merged_cedit_apcorr_adduncat_extra_run2_ecorr_dr-1_lite_fservs.fits'
ftmatchmoc = 'temp_merge_moc.fits'
    
os.system('stilts tmatch2 in1={in1} in2={in2} out={out} values1=ID values2=ID join=all1 find=best suffix2=moc matcher=exact'.format(in1=ft, in2=ftmoc, out=ftmatchmoc))

#t.add_column(Column(name='FLAG_SERVS', data=np.zeros(len(t),dtype=bool)))
#t['FLAG_SERVS'] = tmoc['FLAG_SERVS'][minds]

# Ken's outlier flags
ftoutlier = '/beegfs/lofar/deepfields/science_ready_catalogs/EN1_opt_spitzer_merged_vac_opt3as_irac4as_all_hpx_public.fits'
ftmatchout = 'temp_merge_out.fits'

os.system('stilts tmatch2 in1={in1} in2={in2} out={out} values1=ID values2=id join=all1 find=best suffix2=outlier matcher=exact'.format(in1=ft, in2=ftoutlier, out=ftmatchout))


t = Table.read(ft)
toutlier = Table.read(ftmatchout)
tmoc = Table.read(ftmatchmoc)

for col in t.colnames:
    if 'flux_corr' in col:
        col1 = col.replace('_corr','')
        Table([t[col], toutlier[col1]])
        print col1, col, np.sum(toutlier[col1] == -90),np.sum(np.isnan(t[col])),
        t[col][toutlier[col1] == -90] = np.nan
        print np.sum(np.isnan(t[col]))

t.add_column(tmoc['FLAG_SERVS'])

t.write('/beegfs/lofar/wwilliams/lofar_surveys/deep/agnfitter/EN1_final_cross_match_catalogue-v0.5_mergenew_outlier.fits',overwrite=True)


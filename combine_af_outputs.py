import sys
import os
import glob
from astropy.table import Table, vstack, join

field = sys.argv[1]

version = 'v0.7_newtest'
version ='v0.8'
version = 'v1.0'

if field=='all':
    fields = ['Bootes','EN1','LH']
else:
    fields = [field]
    
for field in fields:


    sourcecat = Table.read(field+'_agnfitter'+ '_' + version+'.fits')

    #sourcecat.keep_columns(['Source_Name', 'radioID', 'ID', 'z1_median', 'Z_BEST', 'newXIDp'])
    sourcecat.keep_columns(['Source_Name', 'radioID', 'ID', 'z1_median', 'Z_BEST', 'XID+_rerun_mips','XID+_rerun_pacs','XID+_rerun_SPIRE', 'CHANGE_FLAG_ID', 'CHANGE_FLAG_ZBEST', 'CHANGE_FLAG_FIR', 'CHANGE_FLAG_DEEP', 'FLAG_GOOD'])


    outtables = []

    print('getting outlists')
    outlist = glob.glob('{field}_{v}/agnfitter_out/*/*.ascii'.format(field=field,v=version))
    nalist = glob.glob('{field}_{v}/agnfitter_out/*/na*'.format(field=field,v=version))

    fitslist = glob.glob('{field}_{v}/agnfitter_out/*/*.txt'.format(field=field,v=version))

    #for f in fitslist:
        #t  =Table.read(f)
        ##print f
        ###t.write(f.replace('.fits','.ascii'),format='ascii',overwrite=True)
        #outtables.append(Table.read(f))
    print( len(outlist),'tables to append')

    #with open('fits.list','w') as f:
        #for ff in fitslist:
            #f.write(ff+'\n')
        #cmd = 'stilts tcat ifmt=fits lazy=true in=@fits.list out=tjoin.fits' 


    print('concatenating')
    with open('tout.ascii','w') as f:
        f.write('radioID metal_med tau_med age_med EBVgal_med Tdust_med fracPAH_med Nh_med EBVbbb_med GA_med SB_med TO_med BB_med log_Mstar_med SFR_opt_med LIR_8_1000_med Lbb_0p1_1_med Lbbdered_0p1_1_med Lga_0p1_1_med Lga_1_100_med Ltor_1_100_med Lsb_1_100_med Lbb_1_30_med Lbbdered_1_30_med Lga_1_30_med Ltor_1_30_med Lsb_1_30_med Lga_1_8_med Ltor_1_8_med Lsb_1_8_med SFR_IR_med neg_ln_like_med metal_p16 tau_p16 age_p16 EBVgal_p16 Tdust_p16 fracPAH_p16 Nh_p16 EBVbbb_p16 GA_p16 SB_p16 TO_p16 BB_p16 log_Mstar_p16 SFR_opt_p16 LIR_8_1000_p16 Lbb_0p1_1_p16 Lbbdered_0p1_1_p16 Lga_0p1_1_p16 Lga_1_100_p16 Ltor_1_100_p16 Lsb_1_100_p16 Lbb_1_30_p16 Lbbdered_1_30_p16 Lga_1_30_p16 Ltor_1_30_p16 Lsb_1_30_p16 Lga_1_8_p16 Ltor_1_8_p16 Lsb_1_8_p16 SFR_IR_p16 neg_ln_like_p16 metal_p84 tau_p84 age_p84 EBVgal_p84 Tdust_p84 fracPAH_p84 Nh_p84 EBVbbb_p84 GA_p84 SB_p84 TO_p84 BB_p84 log_Mstar_p84 SFR_opt_p84 LIR_8_1000_p84 Lbb_0p1_1_p84 Lbbdered_0p1_1_p84 Lga_0p1_1_p84 Lga_1_100_p84 Ltor_1_100_p84 Lsb_1_100_p84 Lbb_1_30_p84 Lbbdered_1_30_p84 Lga_1_30_p84 Ltor_1_30_p84 Lsb_1_30_p84 Lga_1_8_p84 Ltor_1_8_p84 Lsb_1_8_p84 SFR_IR_p84 neg_ln_like_p84 metal_ml tau_ml age_ml EBVgal_ml Tdust_ml fracPAH_ml Nh_ml EBVbbb_ml GA_ml SB_ml TO_ml BB_ml log_Mstar_ml SFR_opt_ml LIR_8_1000_ml Lbb_0p1_1_ml Lbbdered_0p1_1_ml Lga_0p1_1_ml Lga_1_100_ml Ltor_1_100_ml Lsb_1_100_ml Lbb_1_30_ml Lbbdered_1_30_ml Lga_1_30_ml Ltor_1_30_ml Lsb_1_30_ml Lga_1_8_ml Ltor_1_8_ml Lsb_1_8_ml SFR_IR_ml neg_ln_like_ml'+'\n')
    # this is now done in run_subset_agnfitter
    #os.system("tail -n1 -q `ls {field}/agnfitter_out/*/*.ascii` >> tout.ascii".format(field=field))
    os.system("cat {field}_{v}/agnfitter_out/*/*.ascii >> tout.ascii".format(field=field,v=version))

    outtable = Table.read('tout.ascii', format='ascii')

    # add flag where AF couldn't be run
    natab = '# radioID afFLAG\n'
    for f in nalist:
        natab += '{id} 1\n'.format(id=f.split('_')[-1])
    nacat = Table.read(natab, format='ascii')


    #for i in range(len(sourcecat)):
        #radioID = sourcecat['radioID'][i]
        #jobdir = '{field}/agnfitter_out/{rid}'.format(field=field, rid=radioID)
        #outfits = '{jobdir}/parameter_outvalues_{job_ind}.fits'.format(job_ind=radioID, jobdir=jobdir)
        
        #if outfits in fitslist:
            #outtables.append(Table.read(outfits))
            
    #outtable = vstack(outtables, join_type='exact')

    #outtable = Table.read('tjoin.fits')

    fullout = join(sourcecat, outtable, keys='radioID', join_type='left')
    fullout = join(fullout, nacat, keys='radioID', join_type='left')
    fullout['afFLAG'][fullout['afFLAG']!=1] = 0

    fullout.write(field+ '_' + version+'_agnfitter_out.fits', overwrite=True)

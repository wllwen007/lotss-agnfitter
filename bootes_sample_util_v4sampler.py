import LF_util
import numpy as np
from utils.fits_util import load_fits
import astropy.coordinates as ac
import k3match
import matplotlib.pyplot as plt
import utils.plot_util as pp
import os

from t_calc_kcor import calc_kcor


#CATPATH = '/home/wwilliams/'
#CATPATH = '/home/wendy/'
#CATPATH = '/net/laak/data2/wwilliams/phd/multi_bootes/'

if 'uhppc58' in os.uname()[1]:
    local = True
else:
    local = False
if local:
    CATPATH = '/local/wwilliams/phd/multi_bootes/'
else:
    CATPATH = '/data/lofar/wwilliams/phd/multi_bootes/'

#fast_mass_offset = +0.45
fast_mass_offset = 0.

maglim_bright = 10.  # quick check - this is lower than all mags in cat
#maglim_faint = 24.5
maglim_faint = 24.
#maglim_faint = 23.5
#maglim_faint = 23.3
fmaglim_bright =  1.e10*10**(-0.4*maglim_bright)  # AB mag  - 25 zp
fmaglim_faint =  1.e10*10**(-0.4*maglim_faint)    # AB mag  - 25 zp
fmaglim_faint =  2412.50*10**(-0.4*maglim_faint)    # ? mag  - 25 zp

radiofluxlim = 0.5e-3               #very open - only excludes 6 sources
#radiofluxlim = 0.1e-3               #very open - only excludes 6 sources
#radiofluxlim = 1.e-3               #high flux density limit (changed from 100)
#radiofluxlim = 1.5e-3               #high flux density limit (changed from 100)
#radiofluxlim = 2.e-3               #high flux density limit (changed from 100)


areaBOOTES = 19.2*(np.pi/180.)**2.      # 
areaBOOTES1 = 7.*(np.pi/180.)**2.      #  should be the combined area
areaBOOTES = 8.98*(np.pi/180.)**2.      #  should be the combined area
areaBOOTES = 9.247*(np.pi/180.)**2.      #  should be the combined area? area of rms map, masked by ndwfs coverage... but ndwfs coverage is less??
## calculated to be the finite area of the masked rmsmap used below
#rmsmapBOOTES = CATPATH+'/LOFAR150/rms.radio_opt_masked.fits'
rmsmapBOOTES = '/local/wwilliams/projects/radio_imaging/bootes_hba_cycle2_obs4/mosaic_v4/rms_masked.fits'
if local:
    completenessBOOTES = '/local/wwilliams/projects/radio_imaging/bootes_hba_cycle2_obs4/mosaic_v4/CR/LOFAR_BOOTES_Detection_fraction_int_i.npy'
else:
    completenessBOOTES = '/data/lofar/wwilliams/projects/radio_imaging/bootes_hba_cycle2_obs4/mosaic_v4/CR/LOFAR_BOOTES_Detection_fraction_int_i.npy'
    

#areaBOOTES = areaBOOTES*10.

# area calc from /net/laak/data2/wwilliams/phd/multi_bootes/NDWFS/data/mosaic/ndwfs_sdwfs_bootes.source.radiomask.fits
#areaBOOTES_optsamp = 0.01*areaBOOTES  # optical subsample 1
#areaBOOTES_optsamp = 0.1*areaBOOTES  # optical subsample 10
areaBOOTES_optsamp = 0.9*areaBOOTES  # optical subsample 10


#gabyTOsplit = 0.8
#gabyTOsplit = 0.25
gabyTOsplit = 0.15
gabyTOsplit = 0.25
gabyTOsplit_err = 0.2
gabyTOsplit2 = 0.15
gabyTOsplit1 = 0.35
gabyTOsplit_err = 0.2

def plot_seds_EAZY_bin(sed_name, eazydir, opt_ind, savedir,  zrange=(0,6), lrange=[3000,8.e4]):
    
    import threedhst.eazyPy as eazy
    
    #tempfilt, coeffs, temp_seds, pz = readEazyBinary(OUTPUT_DIRECTORY=WD+'OUTPUT/')

    zoutcat = eazydir+'_Iband_sample.zout.fits'
    
    zcat = load_fits(zoutcat)
    if 'z_peak' in zcat.dtype.names:
        zname = 'z_peak'
    elif 'z_a' in zcat.dtype.names:
        zname = 'z_a'
    else:
        zname = 'z_1'
    if 'chi_p' in zcat.dtype.names:
        chiname = 'chi_p'
    elif 'chi_1' in zcat.dtype.names:
        chiname = 'chi_1'
    else:
        chiname = 'chi_a'
    
    Fs = 1.  #plot Jy Hz
    cl = 2.998e19  # angs/s
    
    
    #f1 = pl.figure()
    #ax1 = f1.add_subplot(111)
    #ax1.set_xscale('log')
    
    i = np.where(zcat.id == opt_ind)[0][0]
    
        
    ids = str(opt_ind)
    
    
    lambdaz, temp_sed, lci, obs_sed, fobs, efobs = \
        eazy.getEazySED(i, OUTPUT_DIRECTORY=eazydir+'/OUTPUT')
    
    zgrid, pz = eazy.getEazyPz(i, OUTPUT_DIRECTORY=eazydir+'/OUTPUT' )

    ##### start plot
    f = plt.figure(figsize=[8,5],dpi=100)
    
    #### Plot parameters
    plotsize=20
    alph=0.7
    
        
    ax = f.add_subplot(111)
    f.subplots_adjust(wspace=0.18, hspace=0.0,left=0.09,bottom=0.15,right=0.98,top=0.98)
    ax.set_yscale('log')
    ax.set_xscale('log')
        
    ax.plot(lambdaz, temp_sed, linewidth=1.0, color='blue',alpha=alph)
    
    #### template fluxes integrated through the filters
    ax.scatter(lci,obs_sed,
            c='red',marker='o',s=plotsize,alpha=alph)

    #### Observed fluxes w/ errors
    ax.errorbar(lci,fobs,yerr=efobs,ecolor=None,
            color='black',fmt='o',alpha=alph)
    
    #### Set axis range and titles
    ax.semilogx()
    ax.set_xlim(lrange[0],lrange[1])
    ax.set_ylim(-0.05*max(obs_sed),1.1*max(obs_sed))
    ax.set_xlabel(r'$\lambda$ [$\AA$]')
    ax.set_ylabel(r'$f_\lambda$')
    
    ##### P(z)
    if pz is not None:
        axp = f.add_subplot(444)
            
        axp.plot(zgrid, pz, linewidth=1.0, color='orange',alpha=alph)
        axp.fill_between(zgrid,pz,np.zeros(zgrid.size),color='yellow')

        if zcat['z_spec'][i] > 0:
            axp.plot(zcat['z_spec'][i]*np.ones(2), np.array([0,1e6]),color='red',alpha=0.4)

        #### Set axis range and titles
        axp.set_xlim(0,np.ceil(np.max(zgrid)))
        axp.set_xlim(zrange[0], zrange[1])
        axp.set_ylim(0,1.1*max(pz))
        axp.set_xlabel(r'$z$')
        axp.set_ylabel(r'$p(z)$')
        

    f.savefig(savedir+'/'+sed_name+'.png')
    
    plt.close(f)
    
    
    
    return


def SMenvelope_bootes(z=None, m=None):
    p = np.array([ 8.83868081,  2.93526103,  4.06156463])
    if z is not None:
        return p[0]-p[1]*np.exp(-p[2]*z)
    else:
        t = -1.*np.log((p[0]-m)/p[1])/p[2]
        t[np.isnan(t)] = 100.
        return t

def SMenvelope(z=None, m=None):
    if z is not None:
        zenv,lgMenv = np.loadtxt("/local/wwilliams/projects/LF_by_mass//Mstar_redshift_completeness_emp_uvista_v4.1_95.dat").transpose()
        zenv = np.append(zenv, 100.)
        lgMenv = np.append(lgMenv, 10.991006)
        return np.interp(z, zenv, lgMenv)
    elif m is not None:
        zenv,lgMenv = np.loadtxt("/local/wwilliams/projects/LF_by_mass//Mstar_redshift_completeness_emp_uvista_v4.1_95.dat").transpose()
        zenv = np.append(zenv, 100.)
        lgMenv = np.append(lgMenv, 10.991006)
        return np.interp(m, lgMenv, zenv)
    else:
        return

def SMenvelope_cosmos(z):
    zenv,lgMenv = np.loadtxt("/local/wwilliams/projects/LF_by_mass//Mstar_redshift_completeness_emp_uvista_v4.1_95.dat").transpose()
    zenv = np.append(zenv, 100.)
    lgMenv = np.append(lgMenv, 10.991006)
    return np.interp(z, zenv, lgMenv)


def SMenvelope_cosmos_z(M):
    zenv,lgMenv = np.loadtxt("/local/wwilliams/projects/LF_by_mass//Mstar_redshift_completeness_emp_uvista_v4.1_95.dat").transpose()
    zenv = np.append(zenv, 100.)
    lgMenv = np.append(lgMenv, 10.991006)
    return np.interp(M, lgMenv, zenv)

def get_areas(incat,rmsfits=CATPATH+'/LOFAR150/rms.radio_opt_masked.fits', savefile='radio_opt_areas',corFlux=1):
    import pyfits as pf
    import os
    
    if os.path.exists(savefile+'.npy'):
        areadat = np.load(savefile+'.npy')
        areas = areadat #['areal']
        
        if len(areas) == len(incat):
            return areas
        else:
            print "saved file has different length"
            print len(areas), len(incat)
    else:
        print "file does not exist: "+savefile+".npy"
    
    print 'getting areas'
    #rmsdat = 0.79*pf.getdata(rmsfits)
    rmsdat = pf.getdata(rmsfits)
    rmsdat = rmsdat[np.isfinite(rmsdat)]
    rmsdat = rmsdat[rmsdat>0]
    rmsdat = rmsdat.flatten()
    rmsdat = rmsdat/corFlux
    #rmshead = pf.getheader(rmsfits)
    #ps = rmshead.get('CDELT2')
    pixels = np.sum(np.isfinite(rmsdat))
    
    areas = -1.*np.ones(len(incat))
    #pixels = -1.*np.ones(len(incat))
    for i in range(len(incat)):
        fp = incat['Peak_flux'][i]
        rmslim = fp/4.5
        #rmslim = incat['Resid_Isl_rms'][i]
        #rmslim = incat['Isl_rms'][i]
        pixel = np.sum(rmsdat < rmslim)
        areaf = 1.0*pixel/pixels
        areas[i] = areaf
        #pixels[i] = pixel
    #areas[pixels==0] = np.nan
    #pixels[pixels==0] = np.nan

    areas[areas==0] = 0.25  ### NB force bad ones

    np.save(savefile, areas)

    #incat = append_field(incat, 'Pixels', pixels)
    #incat = append_field(incat, 'Area', areas)
    #pf.writeto(cat_areafile,incat)
    return areas



def get_area_for_flux(flux, snr_limit=5):
    '''calculate area corrections for non-uniform VLA rms'''
    vla_areas = np.load('/net/laak/data2/wwilliams/projects/multi_cosmos/vla/deep_vla_areal.npz')
    indices = [np.sum(fluxi >=snr_limit*vla_areas['flux'])-1 for fluxi in flux]
    areal = vla_areas['areal'][indices]
    return areal


def match_indices(ind_in1, ind_in2, ind_array1, ind_array2):
    idx1, d1 = k3match.nearest_cartesian(ind_array1, np.ones_like(ind_array1), np.ones_like(ind_array1), ind_in1, np.ones_like(ind_in1), np.ones_like(ind_in1))
    idx2, d2 = k3match.nearest_cartesian(ind_array2, np.ones_like(ind_array2), np.ones_like(ind_array2), ind_in2, np.ones_like(ind_in2), np.ones_like(ind_in2))
    indsel = (d1==0) & (d2==0)
    idx1 = idx1[indsel]
    idx2 = idx2[indsel]
    if sum(~indsel) != 0:
        print "not all returned"
        print sum(~indsel)
    return idx1, idx2, indsel


def RL_q24(z):
    return -0.15*(1+z)

def get_q24(S24, S150, lim=None):
    #S24 = 3631*10**(-0.4*m24)  # in Jy
    #S24
    S1400 = S150*(1400./150)**-0.7  #in Jy
    q24 = np.log10(S24/S1400)
    
    MipsGoodInd = ((S24 >0) & (S24<90))
    MipsLimInd = ~MipsGoodInd
    q24[MipsLimInd] = np.nan
    
    q24_limits = np.zeros_like(q24)
    q24_limits[MipsLimInd] = np.log10(3631*10**(-0.4*lim)/S1400[MipsLimInd])
    q24_limits[MipsGoodInd] = np.nan

    return q24, q24_limits

def get_RV(mopt, S150):
    Sr = 3631*10**(-0.4*mopt)  # in Jy
    S1400 = S150*(1400./150)**-0.7  #in Jy
    RV = np.log10(S1400/Sr)
    return RV




def select_good_radio_sample(good=True, withSED=False, withArea=True, plot=True, zsampleind=None, do_kcor=False):
    
    global maglim_bright
    global maglim_faint
    global radiofluxlim
    
    global  gabyTOsplit
    
    match = 'LR'
    
    print "selecting good radio sample...",

    ## NOTE TEMPORARY FILES - TO BE UPDATED ##

    # opt cats with all extra info...
    eazydir = 'eazyfast_irf_b2014_lofar_v4_lr_no_ch3_ch4'
    sample = load_fits(CATPATH+'run_eazy_lofar/{eazydir}_Iband_sample.cat.fits'.format(eazydir=eazydir))
    #zsample = load_fits(CATPATH+'run_eazy_lofar/{eazydir}_Iband_sample.zout.kenz.fits'.format(eazydir=eazydir))
    if zsampleind is not None:
        zsample = load_fits(CATPATH+'run_eazy_lofar/{eazydir}_Iband_sample.zout.kenz.samples.v2.fits'.format(eazydir=eazydir))
        #zsample = t['zsamples'][:,zsampleind]
    
    fsample = load_fits(CATPATH+'run_eazy_lofar/{eazydir}_Iband_sample.fout.fits'.format(eazydir=eazydir))

    # in cat - this is what is input to run_EAZY
    opticalsample = load_fits(CATPATH+'MikeBrown/v2014/Bootes_merged_Icorr_2014a_all_ap4_mags_LOFARv5_LR.fits')
    
    #opticalsample.sort(order='INDEX')
    #import ipdb; ipdb.set_trace()


    rcat = '/local/wwilliams/projects/radio_imaging/bootes_hba_cycle2_obs4/catalogue_v4/'+'mosaic_v4.pbcor.pybdsm.corpos.corflux.optflagmerge.fits'
    lrcat = '/local/wwilliams/projects/radio_imaging/bootes_hba_cycle2_obs4/catalogue_v4/LOFAR_v5_BROWNv2/'+'LOFAR_v5_BROWNv2_matches.fits'
    ocat = CATPATH+'MikeBrown/v2014/Bootes_merged_Icorr_2014a_all_ap4_mags_hmap.anom.fits'

    if not os.path.isfile(CATPATH+'MikeBrown/'+'lrcatsimplerz.fits'):
        print 'matching...'
        print CATPATH+'MikeBrown/'+'lrcatsimplerz.fits'
        s='stilts tmatch2 matcher=exact in2={catr} in1={catlr} values2=LOFAR_id values1=radio_ID join=all1 find=all out={catout}'.format(catr=rcat, catlr=lrcat, catout=CATPATH+'MikeBrown/'+'t.fits')
        os.system(s)
        s='stilts tmatch2 matcher=exact in1={cat1} in2={cato} values1=optical_ID values2=id join=all1 suffix2="_LR" find=best out={catout}'.format(cat1=CATPATH+'MikeBrown/'+'t.fits', cato=ocat, catout=CATPATH+'MikeBrown/'+'lrcat.fits')
        os.system(s)
        
        
        matchrad = 5.  # 2arcsec
        s='stilts tmatch2 matcher=sky in1={cat1} in2={cato} values1="radio_RA radio_DEC" values2="ALPHA_J2000 DELTA_J2000" join=all1 find=best params={radius:f} suffix1="_LR" suffix2="_S"  out={catout}'.format(cat1=CATPATH+'MikeBrown/'+'lrcat.fits', cato=ocat, catout=CATPATH+'MikeBrown/'+'lrcatsimpler.fits', radius=matchrad)
        os.system(s)
        
        zcat1 = CATPATH+'run_eazy_lofar/{eazydir}_Iband_sample.zout.fits'.format(eazydir=eazydir)
        s='stilts tmatch2 matcher=exact in1={cat1} in2={cato} values1=id_LR values2=id join=all1 suffix2="_LR" find=best out={catout}'.format(cat1=CATPATH+'MikeBrown/'+'lrcatsimpler.fits', cato=zcat1, catout=CATPATH+'MikeBrown/'+'lrcatsimplerz.fits')
        os.system(s)
    
    #sample = load_fits(ocat)
    matchdatf = load_fits(CATPATH+'MikeBrown/'+'lrcatsimplerz.fits')
    print len(matchdatf)
    matchdat = matchdatf[matchdatf.Prank==0]
    print len(matchdat)
    Pcut = 0.0
    Pcut = 0.6
    matchdat = matchdat[matchdat.P>Pcut]
    print len(matchdat)

    m_radiosample = matchdat
    m_opticalsample = matchdat
    
    #import ipdb; ipdb.set_trace()
    ii = np.array([np.where(opticalsample.id == t)[0][0] for t in matchdat.id_LR])
    
    
    # this indices should all be the same opticalsample is fed into eazy...
    zii = np.array([np.where(zsample.id == t)[0][0] for t in matchdat.id_LR])
    #fii = np.array([np.where(fsample.id == t)[0][0] for t in matchdat.INDEX_LR])

    m_opticalsample = opticalsample[ii]
    m_zsample = zsample[zii]
    m_fsample = fsample[zii]

    if good:
        
        print "(all {n})...".format(n=len(m_radiosample)) ,
        

        #### FLUX LIMIT ###
        fluxlimind = np.where(m_radiosample.Total_flux > radiofluxlim)[0]   # radio flux limit
        
        print "good (r {n})...".format(n=len(fluxlimind)) ,
        
        m_radiosample = m_radiosample[fluxlimind]
        m_opticalsample = m_opticalsample[fluxlimind]
        m_zsample = m_zsample[fluxlimind]
        m_fsample = m_fsample[fluxlimind]
            

        ##### MAGNITUDE LIMITS ###
        maglimind = np.where((m_opticalsample.MAG_AUTO >= maglim_bright) & (m_opticalsample.MAG_AUTO <= maglim_faint))[0]   # magnitude cut
        
        print "good (o {n})...".format(n=len(maglimind)) ,
        
        m_radiosample = m_radiosample[maglimind]
        m_opticalsample = m_opticalsample[maglimind]
        m_zsample = m_zsample[maglimind]
        m_fsample = m_fsample[maglimind]

    # Redshift
    zspec = m_zsample['z_spec']
    if zsampleind is not None:
        zz = m_zsample['zsamples'][:,zsampleind]
    else:
        #zz = m_zsample['z_peak']
        zz = m_zsample['z1_median']
        
    #has_no_z2 = ~(m_zsample['z2_median'] > 0)
    #z1_q = (m_zsample['z1_max']-m_zsample['z1_min'])/(1+m_zsample['z1_median']) < 0.2
    #keep_z = has_no_z2 & z1_q
    keep_z = np.ones(len(zz),dtype=bool)  # temp keep everything
    #keep_z = (zspec>0) | ((zspec>-90)&(zspec<0))  # temp keep only zspecs
    print 
    print 'throwing away {n:.0f} bad phot-zs'.format(n=np.sum(~keep_z))
    print 'keeping {n:.0f} good phot-zs'.format(n=np.sum(keep_z))
    zz[~keep_z] = np.nan
    
    ## use zspec where available ##
    zz[zspec>0] = zspec[zspec>0]
    zz[(zspec>-90)&(zspec<0)] = -1*zspec[(zspec>-90)&(zspec<0)]  # include also tha 'bad' -ve zspec values
    
    has_spec = (zspec>0) | ((zspec>-90)&(zspec<0))
    
    print '{n:.0f} bad phot-zs have spec-zs'.format(n=np.sum((~keep_z) & has_spec))
    print '{n:.0f} good phot-zs have spec-zs'.format(n=np.sum((keep_z) & has_spec))
    zz[zz<-90] = np.nan
    
    #zphot = m_zsample['z_peak']
    zphot = m_zsample['z1_median']

    # Radio flux
    ff = m_radiosample.Total_flux   #in Jy
    power = np.log10(LF_util.RadioPower(ff, zz, alpha=0.7))  # assumed spec ind
    
    #XX
    rad_id = m_radiosample.LOFAR_id
    opt_id = m_opticalsample.id
    
    # Magnitudes
    Imag = m_opticalsample.MAG_AUTO 
    #Imag = m_opticalsample.I_MAG   
    fnu_Imag = 2412.50*10.**(-0.4*Imag)    # zp  - in Jy
    Ilum = LF_util.OpticalLuminosity(fnu_Imag, zz)  # in  W/Hz
    
    ##### TESTING ### 
    ## use rest-frame colour
    ##mi = -2.5*np.log10(m_colsample.RF_F4) + 8.9
    #mi = -2.5*np.log10(m_colsample.RF_F3) + 8.9
    #fnu_Imag_0 = 2412.50*10.**(-0.4*mi)    # zp  - in Jy
    #Ilum_0 = LF_util.OpticalLuminosity(fnu_Imag_0, zz)  # in  W/Hz
    
    #if plot:
        #f,ax = pp.paper_single_ax()
        #c = ax.scatter(Ilum_0, Ilum,c=zz,edgecolor='none',vmin=0,vmax=3)
        #cbar = plt.colorbar(c)
        #ax.plot([1e46, 1e52], [1e46, 1e52],'k')
        #ax.loglog()
        #cbar.set_label('$z$')
        #ax.set_ylabel('$L_I^0$ (rest)')
        #ax.set_xlabel('$L_I$')
        #pp.fig_save_many(f,'bootes_sample_select_good_radio_sample_opt_lum_kcor.png')
        #plt.close()
        #f,ax = pp.paper_single_ax()
        #ax.plot(zz, -2.5*np.log10(Ilum/Ilum_0),'k.',alpha=0.25)
        #ax.plot([0, 3], [0, 0],'k')
        #ax.set_xlim(0,3)
        ##ax.loglog()
        #ax.set_ylabel('$-2.5 \log L_I/L_I^0$')
        #ax.set_xlabel('$z$')
    
    #dM0 = -2.5*np.log10(Ilum/Ilum_0)
    #print 'using median dM for sources with no RF Imag'
    #nanind = np.where(np.isnan(Ilum_0))[0]
    #for i in nanind:
        #zi = zz[i]
        #Ilumi = Ilum[i]
        #dz = np.max((0.1, 0.1*zi))
        #imed = np.where((zz>zi-dz)&(zz<zi+dz))[0]
        #dMmed = np.nanmedian(dM0[imed])
        #Ilum_0[i] = Ilum[i]*10**(0.4*dMmed)

    #if plot:
        #ax.plot(zz[nanind], -2.5*np.log10(Ilum/Ilum_0)[nanind],'r.',alpha=0.25)
        
        #pp.fig_save_many(f,'bootes_sample_select_good_radio_sample_opt_lum_kcor2.png')
        #plt.close()
        
    ## use rest frame
    if do_kcor:
        Ilum = Ilum_0
    
    #BI = m_opticalsample['magBw'] - m_opticalsample['magI']
    #B = m_opticalsample['magBw']
    #kc = calc_kcor('B',zz, 'B - Ic', BI)
    #print kc
    #mi = m_opticalsample.MAG_AUTO  + kc
    #fnu_Imag = 2412.50*10.**(-0.4*mi)    # zp  - in Jy
    #Ilum = LF_util.OpticalLuminosity(fnu_Imag, zz)  # in  W/Hz
    
    # MIR mags and LUM
    #mirmag = m_opticalsample['24_flux'] 
    ##mirmag2 = m_opticalsample['ch4_MAG']
    #mirmag2 = m_opticalsample['ch4_flux']
    fnu_mirmag = m_opticalsample['24_flux'] 
    fnu_mirmag2 = m_opticalsample['ch4_flux']
    #fnu_mirmag = 3631.*10.**(-0.4*mirmag)    # MIPS
    #fnu_mirmag[(mirmag<0)] = np.nan
    #fnu_mirmag[(mirmag>90)] = np.nan
    #fnu_mirmag2 = 63.10*10.**(-0.4*mirmag2)    # IRAC4
    #fnu_mirmag2[(mirmag2<0)] = np.nan
    #fnu_mirmag2[(mirmag2>90)] = np.nan
    miralpha = np.log10(fnu_mirmag/fnu_mirmag2)/np.log10(24./8.)
    mirlum = LF_util.OpticalLuminosity2(fnu_Imag, zz, miralpha)

    #import ipdb ;  ipdb.set_trace()
    
    Imag_AB = -2.5*np.log10(fnu_Imag/3631.)
    IMag_AB = LF_util.OpticalMag(Imag_AB, zz)


    #cmag_I = m_opticalsample.I_MAG
    #cmag_Bw = m_opticalsample.Bw_MAG 
    #fnu_cmag_I = 2412.50*10.**(-0.4*cmag_I)    # zp
    #fnu_cmag_Bw = 3627.50*10.**(-0.4*cmag_Bw)    # zp
    fnu_cmag_I = m_opticalsample.I_flux
    fnu_cmag_Bw = m_opticalsample.Bw_flux
    cmag_I_AB = -2.5*np.log10(fnu_cmag_I/3631.)
    cmag_Bw_AB = -2.5*np.log10(fnu_cmag_Bw/3631.)
    cMag_I_AB = LF_util.OpticalMag(cmag_I_AB, zz)
    cMag_Bw_AB = LF_util.OpticalMag(cmag_Bw_AB, zz)
    cMag_BI = cMag_Bw_AB-cMag_I_AB


    smass = m_fsample.lmass + fast_mass_offset
    sfr = m_fsample.lsfr
    
    #XX
    
    if withArea:
        if good:
            areafl = get_areas(m_radiosample, savefile='radio_opt_good_areas')
        else:
            areafl = get_areas(m_radiosample, savefile='radio_opt_areas')
        print areafl


    if withSED:
        ## updated 10 aug
        agnfitterout = load_fits(CATPATH+"/agnfitter/OUTPUT_kenz_old/agnfitter_out.fits") 
        ## updated 27 jul 2017
        
        
        if zsampleind is not None:
            agnfitterout = load_fits(CATPATH+"/agnfitter/OUTPUT_v2_zsamples_v2/agnfitter_out_{i:03d}.fits".format(i=zsampleind)) ##

        m_optical_hasgaby = np.array([ src in agnfitterout.id for src in m_opticalsample.id])

        m_optical_agnfitterout_ind = np.array([ np.where(agnfitterout.id == src)[0][0] for src in m_opticalsample.id if  src in agnfitterout.id ])
        

        l10 = np.log(10.)

        #import ipdb; ipdb.set_trace()
        # in W
        LIR_e = 1e-7*10**np.float64(agnfitterout['LIR_8_1000_med'])
        LIR_em = LIR_e*l10*np.float64(agnfitterout['LIR_8_1000_med']-agnfitterout['LIR_8_1000_errlow'])
        LIR_ep = LIR_e*l10*np.float64(agnfitterout['LIR_8_1000_errhigh']-agnfitterout['LIR_8_1000_med'])
        
        Slam = '_1_30_'
        #Slam = '_1_100_'
        
        TO_e = 10**np.float64(agnfitterout['Ltor'+Slam+'med'])
        BB_e = 10**np.float64(agnfitterout['Lbb'+Slam+'med'])
        SB_e = 10**np.float64(agnfitterout['Lsb'+Slam+'med'])
        GA_e = 10**np.float64(agnfitterout['Lga'+Slam+'med'])  ##NB: the L_galaxy in the same freq regime as the torus and starburst
        TO_em = TO_e*l10*np.float64(agnfitterout['Ltor'+Slam+'med']-agnfitterout['Ltor'+Slam+'errlow'])
        BB_em = TO_e*l10*np.float64(agnfitterout['Lbb'+Slam+'med']-agnfitterout['Lbb'+Slam+'errlow'])
        SB_em = SB_e*l10*np.float64(agnfitterout['Lsb'+Slam+'med']-agnfitterout['Lsb'+Slam+'errlow'])
        GA_em = GA_e*l10*np.float64(agnfitterout['Lga'+Slam+'med']-agnfitterout['Lga'+Slam+'errlow'])
        TO_ep = TO_e*l10*np.float64(agnfitterout['Ltor'+Slam+'errhigh']-agnfitterout['Ltor'+Slam+'med'])
        BB_ep = TO_e*l10*np.float64(agnfitterout['Lbb'+Slam+'errhigh']-agnfitterout['Lbb'+Slam+'med'])
        SB_ep = SB_e*l10*np.float64(agnfitterout['Lsb'+Slam+'errhigh']-agnfitterout['Lsb'+Slam+'med'])
        GA_ep = GA_e*l10*np.float64(agnfitterout['Lga'+Slam+'errhigh']-agnfitterout['Lga'+Slam+'med'])
        
        
        fSB_e = SB_e/(GA_e+SB_e)
        
        #Tot_e =  GA_e+SB_e+TO_e
        #fSB1_e = SB_e/Tot_e
        #fTO_e = TO_e/Tot_e
        
        #new fTO
        Tot1_e =  GA_e+BB_e+TO_e
        fTO_e = (TO_e+BB_e)/Tot1_e
        
        fTO_em = ((GA_e**2.*(TO_em**2.+BB_em**2.) + (TO_e+GA_e)**2.*GA_em**2.)**0.5)/(Tot1_e**2.) 
        fTO_ep = ((GA_e**2.*(TO_ep**2.+BB_ep**2.) + (TO_e+GA_e)**2.*GA_ep**2.)**0.5)/(Tot1_e**2.) 
        
        
        #fTO_e = TO_e/(TO_e+GA_e)
        ## errors #
        #fTO_em = fTO_e**2.*(GA_e/TO_e)*np.sqrt( (TO_em/TO_e)**2. + (GA_em/GA_e)**2. )
        #fTO_ep = fTO_e**2.*(GA_e/TO_e)*np.sqrt( (TO_ep/TO_e)**2. + (GA_ep/GA_e)**2. )
        
        fSB_em = fSB_e**2.*(GA_e/SB_e)*np.sqrt( (SB_em/SB_e)**2. + (GA_em/GA_e)**2. )
        fSB_ep = fSB_e**2.*(GA_e/SB_e)*np.sqrt( (SB_ep/SB_e)**2. + (GA_ep/GA_e)**2. )
        
        
        Tot_em = np.sqrt( (TO_em/TO_e)**2. + (GA_em/GA_e)**2. +  (SB_em/SB_e)**2.)
        Tot_ep = np.sqrt( (TO_ep/TO_e)**2. + (GA_ep/GA_e)**2. +  (SB_ep/SB_e)**2.)
        #fTO1_em = fTO1_e*np.sqrt( (TO_em/TO_e)**2. + (Tot_em/Tot_e)**2. )
        #fTO1_ep = fTO1_e*np.sqrt( (TO_ep/TO_e)**2. + (Tot_ep/Tot_e)**2. )
        #fSB1_em = fSB1_e*np.sqrt( (SB_em/SB_e)**2. + (Tot_em/Tot_e)**2. )
        #fSB1_ep = fSB1_e*np.sqrt( (SB_ep/SB_e)**2. + (Tot_ep/Tot_e)**2. )
        
          
        ln_l = np.nan*np.zeros(len(m_opticalsample))
        ln_l[m_optical_hasgaby] = agnfitterout['neg_ln_like_med'][m_optical_agnfitterout_ind]
        
        
        LIR = np.nan*np.zeros(len(m_opticalsample))
        LIR[m_optical_hasgaby] = LIR_e[m_optical_agnfitterout_ind]
        LIR_e_u = np.nan*np.zeros(len(m_opticalsample))
        LIR_e_u[m_optical_hasgaby] = LIR_ep[m_optical_agnfitterout_ind]
        LIR_e_l = np.nan*np.zeros(len(m_opticalsample))
        LIR_e_l[m_optical_hasgaby] = LIR_em[m_optical_agnfitterout_ind]
        
        fTO = np.nan*np.zeros(len(m_opticalsample))
        fTO[m_optical_hasgaby] = fTO_e[m_optical_agnfitterout_ind]
        fSB = np.nan*np.zeros(len(m_opticalsample))
        fSB[m_optical_hasgaby] = fSB_e[m_optical_agnfitterout_ind]
        
        fTO_e_u = np.nan*np.zeros(len(m_opticalsample))
        fTO_e_u[m_optical_hasgaby] = fTO_ep[m_optical_agnfitterout_ind]
        fTO_e_l = np.nan*np.zeros(len(m_opticalsample))
        fTO_e_l[m_optical_hasgaby] = fTO_em[m_optical_agnfitterout_ind]
        fSB_e_u = np.nan*np.zeros(len(m_opticalsample))
        fSB_e_u[m_optical_hasgaby] = fSB_ep[m_optical_agnfitterout_ind]
        fSB_e_l = np.nan*np.zeros(len(m_opticalsample))
        fSB_e_l[m_optical_hasgaby] = fSB_em[m_optical_agnfitterout_ind]
        
        
        L150 = 10**power
        q150 =  np.log10((LIR / 3.75e12) / L150)
        sigSF = 3.5
        sedSF = 1.*(q150 >= 1.62-sigSF*0.24)  # Ivison+ 2010 - alpha = -0.8
        
        sigSF = 2.5
        sedSF = (q150 >= (1.72*(1+zz)**-0.22)-sigSF*0.529)  # Ivison+ 2010
        
        sigSF = 1.5
        sedSF = (q150 >= (1.72*(1+zz)**-0.22)-sigSF*0.529) & (power < 26.)  # Ivison+ 2010
        
        sedHERG = 1.*((fTO >= gabyTOsplit) & ~sedSF)
        sedLERG = 1.*((fTO <  gabyTOsplit) & ~sedSF)
        sedHERGc = 1.*((fTO >= gabyTOsplit1) & ~sedSF)
        sedLERGc = 1.*((fTO <  gabyTOsplit2) & ~sedSF)
        
        sedHERGu = 1.*((fTO >= gabyTOsplit1-gabyTOsplit_err) & ~sedSF)
        sedLERGu = 1.*((fTO <  gabyTOsplit2+gabyTOsplit_err) & ~sedSF)
        
        sedSF = 1.0*sedSF
        sedSF[np.isnan(LIR)] = np.nan
        sedHERG[np.isnan(fTO)] = np.nan
        sedLERG[np.isnan(fTO)] = np.nan
        sedHERGc[np.isnan(fTO)] = np.nan
        sedLERGc[np.isnan(fTO)] = np.nan
        sedHERGu[np.isnan(fTO)] = np.nan
        sedLERGu[np.isnan(fTO)] = np.nan
        
        print 
        print "agnfitterout: {f:d} of {f2:d} bad fits".format(f=np.sum(ln_l<-100),f2=len(ln_l<-100))
        #print "agnfitterout: {f:d} of {f2:d}".format(f=len(m_optical_hasgaby), f2=len(m_opticalsample))
        
        fTO[ln_l<-100] = np.nan
        fTO[ln_l<-100] = np.nan
        fSB[ln_l<-100] = np.nan
        fSB[ln_l<-100] = np.nan
        
        
        sedL = ln_l
        
    

    
    savearr = np.array(( opt_id, rad_id, zz, zspec, zphot, ff, power, Imag, Imag_AB, IMag_AB, cMag_BI, Ilum, mirlum, smass, sfr ))
    savenames = 'opt_id, rad_id, z, zspec, zphot, radio_flux, power, opt_mag, opt_mag_AB, opt_Mag_AB, opt_Mag_BI_AB, opt_lum, mir_lum, smass, sfr'
    savefmts = 'i8,i8, f8, f8, f8, f8, f8, f8, f8, f8, f8, f8, f8, f8, f8'
    if withArea:
        savearr = np.vstack((savearr, np.array(areafl)))
        savenames = savenames+ ',areal'
        savefmts = savefmts+',f8'
        
    if withSED:
        savearr = np.vstack((savearr, np.array((LIR,LIR_e_u,LIR_e_l,fTO,fTO_e_u,fTO_e_l,fSB,fSB_e_u,fSB_e_l, sedHERG,sedLERG, sedHERGc,sedLERGc, sedHERGu,sedLERGu,sedSF,sedL))))
        savenames = savenames+ ',LIR,LIR_eu,LIR_el,fTO,fTO_eu,fTO_el,fSB,fSB_eu,fSB_el,sedHERG,sedLERG,sedHERGc,sedLERGc,sedHERGu,sedLERGu,sedSF,sedL'
        savefmts = savefmts+',f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8'
        
        
        
    lf_cat = np.core.records.fromarrays(savearr, names=savenames, formats = savefmts)
    
    print "done"
    

    return m_opticalsample, lf_cat



def select_good_sample(good=True, withMIR=False, withcol=True, do_kcor=False, plot=False):
    
    global maglim_bright
    global maglim_faint
    global radiofluxlim
    
    print "selecting good sample...",

    ## NOTE TEMPORARY FILES - TO BE UPDATED ##

    # opt cats with all extra info...
    eazydir = 'eazyfast_irf_bright_b2014_no_ch3_ch4'
    sample = load_fits(CATPATH+'run_eazy_lofar/{eazydir}_Iband_sample.cat.fits'.format(eazydir=eazydir))
    zsample = load_fits(CATPATH+'run_eazy_lofar/{eazydir}_Iband_sample.zout.fits'.format(eazydir=eazydir))
    

    # in cat - this is what is input to run_EAZY
    opticalsample = load_fits(CATPATH+'MikeBrown/v2014/Bootes_merged_Icorr_2014a_all_ap4_mags_hmap.anom.fits')

    sel_opt_ind = np.searchsorted(opticalsample.id, sample.id)
    opticalsample = opticalsample[sel_opt_ind]



    Noptical = len(opticalsample)
    Nsrcs = Noptical


    #import ipdb ; ipdb.set_trace()
    if good:
        ##### MAGNITUDE LIMITS ###
        maglimind = np.where((opticalsample.MAG_AUTO >= maglim_bright) & (opticalsample.MAG_AUTO <= maglim_faint))   # magnitude cut
        
        m_opticalsample = opticalsample[maglimind]
        m_zsample = zsample[maglimind]
        
        ## GOOD DEEP ##
        deepind = np.where(m_opticalsample.FLAG_DEEP == 1)   # magnitude cut
        
        m_opticalsample = opticalsample[deepind]
        m_zsample = m_zsample[deepind]
        
    else:
        m_opticalsample = opticalsample
        m_zsample = zsample

    Noptical = len(m_opticalsample)
    Nsrcs = Noptical


    # matched samples:



    # get the MIR sources - for all optical sources
    if withMIR:
        print "adding MIR...",
        
        #m1 = m_opticalsample.ch1_ma
        #m2 = m_opticalsample.ch2_ma
        #m3 = m_opticalsample.ch3_ma
        #m4 = m_opticalsample.ch4_ma
        #miragn1 = pp.stern_mask(m1, m2, m3, m4)
        #miragn2 = pp.donley_mask(m1, m2, m3, m4)
        m1 = -2.5*np.log10(m_opticalsample['ch1_flux']/3631.)
        m2 = -2.5*np.log10(m_opticalsample['ch2_flux']/3631.)
        m3 = -2.5*np.log10(m_opticalsample['ch3_flux']/3631.)
        m4 = -2.5*np.log10(m_opticalsample['ch4_flux']/3631.)
        #m1 = m_opticalsample.ch1_ma
        #m2 = m_opticalsample.ch2_ma
        #m3 = m_opticalsample.ch3_ma
        #m4 = m_opticalsample.ch4_ma
        miragn1 = pp.stern_mask(m1, m2, m3, m4)
        miragn2 = pp.donley_mask(m1, m2, m3, m4)


    # Redshift
    zspec = m_zsample['z_spec']
    zz = m_zsample['z_peak']
    ## use zspec where available ##
    zz[zspec>0] = zspec[zspec>0]
    zz[(zspec>-90)&(zspec<0)] = -1*zspec[(zspec>-90)&(zspec<0)]  # include also tha 'bad' -ve zspec values
    
    zz[zz<-90] = np.nan
    
    zphot = m_zsample['z_peak']
    
    #import ipdb ; ipdb.set_trace()
    # Magnitudes
    Imag = m_opticalsample.MAG_AUTO 
    #Imag = m_opticalsample.I_MAG   
    fnu_Imag = 2412.50*10.**(-0.4*Imag)    # zp
    Ilum = LF_util.OpticalLuminosity(fnu_Imag, zz)
    
        #### TESTING ### 
    if do_kcor:
        # use rest-frame colour
        #mi = -2.5*np.log10(m_colsample.RF_F4) + 8.9
        mi = -2.5*np.log10(m_colsample.RF_F3) + 8.9
        fnu_Imag_0 = 2412.50*10.**(-0.4*mi)    # zp  - in Jy
        Ilum_0 = LF_util.OpticalLuminosity(fnu_Imag_0, zz)  # in  W/Hz
        
        if plot:
            f,ax = pp.paper_single_ax()
            c = ax.scatter(Ilum_0, Ilum,c=zz,edgecolor='none',vmin=0,vmax=3)
            cbar = plt.colorbar(c)
            ax.plot([1e46, 1e52], [1e46, 1e52],'k')
            ax.loglog()
            cbar.set_label('$z$')
            ax.set_ylabel('$L_I^0$ (rest)')
            ax.set_xlabel('$L_I$')
            pp.fig_save_many(f,'bootes_sample_select_good_sample_opt_lum_kcor.png')
            plt.close()
            f,ax = pp.paper_single_ax()
            ax.plot(zz, -2.5*np.log10(Ilum/Ilum_0),'k.',alpha=0.25)
            ax.plot([0, 3], [0, 0],'k')
            ax.set_xlim(0,3)
            #ax.loglog()
            ax.set_ylabel('$-2.5 \log L_I/L_I^0$')
            ax.set_xlabel('$z$')
        
        dM0 = -2.5*np.log10(Ilum/Ilum_0)
        print 'using median dM for sources with no RF Imag'
        nanind = np.where(np.isnan(Ilum_0))[0]
        for i in nanind:
            zi = zz[i]
            Ilumi = Ilum[i]
            dz = np.max((0.1, 0.1*zi))
            imed = np.where((zz>zi-dz)&(zz<zi+dz))[0]
            dMmed = np.nanmedian(dM0[imed])
            Ilum_0[i] = Ilum[i]*10**(0.4*dMmed)

        if plot:
            ax.plot(zz[nanind], -2.5*np.log10(Ilum/Ilum_0)[nanind],'r.',alpha=0.25)
            
            pp.fig_save_many(f,'bootes_sample_select_good_sample_opt_lum_kcor2.png')
            plt.close()
        
        
        Ilum = Ilum_0
    
    # MIR mags and LUM
    #mirmag = m_opticalsample['24_MAG'] 
    #mirmag2 = m_opticalsample['ch4_MAG']
    #fnu_mirmag = 3631.*10.**(-0.4*mirmag)    # MIPS
    #fnu_mirmag[(mirmag<0)] = np.nan
    #fnu_mirmag[(mirmag>90)] = np.nan
    #fnu_mirmag2 = 63.10*10.**(-0.4*mirmag2)    # IRAC4
    #fnu_mirmag2[(mirmag2<0)] = np.nan
    #fnu_mirmag2[(mirmag2>90)] = np.nan
    #miralpha = np.log10(fnu_mirmag/fnu_mirmag2)/np.log10(24./8.)
    #mirlum = LF_util.OpticalLuminosity2(fnu_Imag, zz, miralpha)


    fnu_mirmag = m_opticalsample['24_flux'] 
    fnu_mirmag2 = m_opticalsample['ch4_flux']
    #fnu_mirmag = 3631.*10.**(-0.4*mirmag)    # MIPS
    #fnu_mirmag[(mirmag<0)] = np.nan
    #fnu_mirmag[(mirmag>90)] = np.nan
    #fnu_mirmag2 = 63.10*10.**(-0.4*mirmag2)    # IRAC4
    #fnu_mirmag2[(mirmag2<0)] = np.nan
    #fnu_mirmag2[(mirmag2>90)] = np.nan
    miralpha = np.log10(fnu_mirmag/fnu_mirmag2)/np.log10(24./8.)
    mirlum = LF_util.OpticalLuminosity2(fnu_Imag, zz, miralpha)


    Imag_AB = -2.5*np.log10(fnu_Imag/3631.)
    IMag_AB = LF_util.OpticalMag(Imag_AB, zz)

    #
    #cmag_I = m_opticalsample.I_MAG
    #cmag_Bw = m_opticalsample.Bw_MAG 
    #fnu_cmag_I = 2412.50*10.**(-0.4*cmag_I)    # zp
    #fnu_cmag_Bw = 3627.50*10.**(-0.4*cmag_Bw)    # zp
    fnu_cmag_I = m_opticalsample.I_flux
    fnu_cmag_Bw = m_opticalsample.Bw_flux
    cmag_I_AB = -2.5*np.log10(fnu_cmag_I/3631.)
    cmag_Bw_AB = -2.5*np.log10(fnu_cmag_Bw/3631.)
    cMag_I_AB = LF_util.OpticalMag(cmag_I_AB, zz)
    cMag_Bw_AB = LF_util.OpticalMag(cmag_Bw_AB, zz)
    cMag_BI = cMag_Bw_AB-cMag_I_AB

    smass = m_fsample.lmass + fast_mass_offset
    sfr = m_fsample.lsfr
    
    

    savearr = np.array((zz, zspec, zphot, Imag, Imag_AB, IMag_AB, cMag_BI, Ilum, mirlum, smass, sfr))
    savenames = 'z, zspec, zphot, opt_mag, opt_mag_AB, opt_Mag_AB, opt_Mag_BI_AB, opt_lum, mir_lum, smass, sfr'
    #savearr = np.array((zz, zspec, Imag, Ilum, smass, sfr))
    #savenames = 'z, zspec, opt_mag, opt_lum, smass, sfr'
    savefmts = 'f8, f8, f8, f8, f8, f8, f8, f8, f8, f8, f8'
    print savearr.shape
    print savenames
    print savefmts
    lf_cat = np.core.records.fromarrays(savearr, names=savenames, formats = savefmts)

    print "done"

    return m_opticalsample, lf_cat


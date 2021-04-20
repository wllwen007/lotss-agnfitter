'''
AGNfitter setting file:

required:
CATALOG_settings
FILTERS_settings
MCMC_settings
OUTPUT_settings

For default use (test example with 2 redshifts and default filter set)

Change only the functions which state 
***USER INPUT NEEDED***.
'''


def CATALOG_settings():

    """==================================
    ***USER INPUT NEEDED***

    Set the right values to be able to read your catalog's format.
    FITS option is not available yet.
    =================================="""


    cat = dict()


    ##GENERAL
    cat['path'] ='/beegfs/lofar/wwilliams/lofar_surveys/deep/AGNfitter/'

    cat['outpath'] ='/beegfs/lofar/wwilliams/lofar_surveys/deep/agnfitter/##FIELD##_##VERSION##/'
    
    cat['filename'] = '/beegfs/lofar/wwilliams/lofar_surveys/deep/agnfitter/##FIELD##_agnfitter_##VERSION##.fits'

    cat['use_central_wavelength'] = False # Option to use central wavelength if no wavelengths in table

    ##OTHERS (user input ***not*** needed)
    cat['output_folder'] =  cat['outpath']+'/agnfitter_out/'#if no special OUTPUT folder, leave default
    #cat['output_folder'] =  cat['outpath']+'/agnfitter_out_zbest/'#if no special OUTPUT folder, leave default
    cat['dict_path'] = cat['outpath']+'/MODELSDICT_##FIELD##'


    cat['filetype'] = 'FITS' ## catalog file type: 'ASCII' or 'FITS'. 
    cat['name'] = 'radioID' #'ID'            ## If ASCII: Column index (int) of source IDs
    cat['redshift'] = 'Z_BEST' #'z'              ## If ASCII:  Column index(int) of redshift
 
    ##FREQUENCIES/WAVELENGTHS 
    cat['freq/wl_suffix'] = '_wl'
    cat['freq/wl_format'] = 'wavelength' ## Gives your catalog *observed*
    cat['freq/wl_unit'] = u.Angstrom       ## Astropy unit of freq or wavelength

    ##FLUXES  
    cat['flux_unit'] = u.Jy             ## Astropy unit of *flux* (astropy-units)
    cat['flux_suffix'] = '_f'
    cat['fluxerr_suffix'] = '_e'
    cat['ndflag_bool'] = False          ## Does you catalog has columns with flags 1(0) for 
    #cat['ndflag_suffix'] = '_nd'         ## If ASCII: List of column indexes (int)

    return cat


def FILTERS_settings():

    """==================================
    Set the photometric bands included in your catalog,
    in order to integrate the models over their response curves.
    =================================="""

    filters = dict()
    
    #filters['dict_zarray'] = zz  # The grid of redshifts needed to fit your catalog
    #filters['dict_zarray'] = np.arange(0,7.00001,0.01)  # The grid of redshifts needed to fit your catalog
    filters['filterset'] = 'filterset_file' # OPTIONS: 
                                           # 'BANDSET_default' (for testing)
                                           # 'BANDSET_settings' (choosing relevant filters below, as given by your catalog)
                                           # if your filter is not included, go to DICTIONARIES_AGNfitter to add.

    #filters['file'] = '/beegfs/lofar/wwilliams/lofar_surveys/deep/agnfitter/##FIELD##_agnfitter_filters.lst'
    filters['path'] = '/beegfs/lofar/wwilliams/lofar_surveys/deep/agnfitter/##FIELD##_##VERSION##/'  
    filters['file'] = '/beegfs/lofar/wwilliams/lofar_surveys/deep/agnfitter/##FIELD##_##VERSION##_agnfitter_filters.lst' 
    
    filters['add_filters']= True # If 'True' please add them below in ADD FILTERS
    
    
    return filters

def MODELS_settings():

    """==================================
    Work in progress
    =================================="""


    models = dict()
    models['path'] = 'models/' 
    models['modelset'] = 'modelsv1'


    models['GALAXY'] = 'BC03_metal'   ### Current options:
                                ### 'BC03' (Bruzual & Charlot 2003)
                                ### 'BC03_metal' (Bruzual & Charlot 2003), with metallicities
    models['STARBURST'] = 'S17_newmodel' ### Current options:
                                ### 'DH02_CE01' (Dale & Helou 2002 + Chary & Elbaz 2001)
                                ### 'S07' (Schreiber et al. 2017 (submitted))

    models['BBB'] ='R06' ### Current options:
                         ### 'R06' (Richards et al. 2006) ## Needs 2 manual changes in PARAMETERSPACE_AGNfitter.py
                         ### 'SN12' (Slone&Netzer 2012)
                         ### 'D12_S' (Done et al. 2012) for Schwarzschild BH, with x-ray predictions
                         ### 'D12_K' (Done et al. 2012) for Kerr BH, with x-ray predictions

    models['TORUS'] ='S04' ### Current options:
                           ### 'S04' (Silva et al. 2004)

    models['XRAYS'] = False ### If X-ray data is available and informative for the fit

    models['RADIO'] = False ### If radio data is available and informative for the fit

    models['PRIOR_energy_balance'] = True ### Default:True
                                          ### True: Sets a lower limit to the dust emission luminosity ('starburst' model)
                                          ### as given by the observed attenuation in the stellar component SED.
    models['PRIOR_AGNfraction'] = True  ### Default: True
                                        ### True: - *IF* blue/UV bands (around 1500 Angstrom) are 10 times higher than expected by the galaxy luminosity function by Parsa, Dunlop et al. 2014. 
                                        ###         this option rejects AGN-to-GAL ratios lower than 1 (log =0). It then applies a Gaussian prior probability with log ratio=2, with a sigma of 2.
                                        ###       - In this cases it also applies a Gaussian prior on the galaxy normalization, i.e. stellar mass (usually unconstrained in these cases) to 
                                        ###         populate physically expected ranges for QSO hosts -> 10^9 - 10^11. 
                                        ###       - *ELSE IF* blue/UV bands (around 1500 Angstrom) are below 10 times the expected value by Parsa, Dunlop et al. 2014. 
                                        ###         this option gives preference to galaxy contribution in the optical UV, with Gaussian prior probability centered on AGN to GALAXY log ratios of -1. 
                                        ###          and sigma 1, i.e. accretion disk is disfavoured at least the data strongly prefers it.
                                        ### False:- Non-informative prior
    models['PRIOR_galaxy_only'] = False ### Default:False 
                                        ### True: sets all AGN contribution to 0.ß
    return models

def MCMC_settings():

    """==================================
    Set your preferences for the MCMC sampling.
    =================================="""

    mc = dict()

    mc['Nwalkers'] = 100  ## number of walkers 
    mc['Nburnsets']= 2   ## number of burn-in sets
    mc['Nburn'] = 4000 ## length of each burn-in sets
    mc['Nmcmc'] = 10000  ## length of each burn-in sets
    mc['iprint'] = 1000 ## show progress in terminal in steps of this many samples

    return mc

def OUTPUT_settings():

    """==================================
    Set your preferences for the production of OUTPUT files. 
    =================================="""

    out = dict()

    out['plot_format'] = 'png'
    out['plot_residuals'] = True

    #CHAIN TRACES
    out['plot_tracesburn-in'] = True    
    out['plot_tracesmcmc'] = True

    #BASIC OUTPUT
    out['Nsample'] = 1000 ## out['Nsample'] * out['Nthinning'] <= out['Nmcmc']
    out['Nthinning'] = 10 ## This describes thinning of the chain to sample
    out['writepar_meanwitherrors'] = True ##Write output values for all parameters in a file.
    out['plot_posteriortriangle'] = True ##Plot triangle with all parameters' PDFs?
    
    out['save_chains'] = False  # saves space for now (may want to not save burnins)

    #INTEGRATED LUMINOSITIES
    out['calc_intlum'] = True  
    out['realizations2int'] = 1000 #This process is very time consuming.
                                #Around 100-1000 is recomendend for computational reasons.
                                #If you want to plot posterior triangles of 
                                #the integrated luminosities, should be > 1000.
    out['plot_posteriortrianglewithluminosities'] = True  # requires out['calc_intlum']=True

    out['save_posterior_luminosities'] = True 

    #INTEGRATION RANGES
    out['intlum_models'] = ['sb','bbb', 'bbbdered', 'gal', 'gal', 'tor','sb','bbb', 'bbbdered', 'gal', 'tor','sb', 'gal', 'tor','sb']  #leave 'sb' always 
                                                                        #as first element
    out['intlum_freqranges_unit'] = u.micron   #Astropy unit 
    out['intlum_freqranges'] = np.array([[8.,1000.],[0.1,1.],[0.1,1.],[0.1,1.],[1.,100.],[1.,100.],[1.,100.],[1.,30.],[1.,30.],[1.,30.],[1.,30.],[1.,30.],[1.,8.],[1.,8.],[1.,8.]])
    out['intlum_names'] = ['LIR_8_1000','Lbb_0p1_1', 'Lbbdered_0p1_1', 'Lga_0p1_1', 'Lga_1_100', 'Ltor_1_100','Lsb_1_100','Lbb_1_30', 'Lbbdered_1_30', 'Lga_1_30', 'Ltor_1_30','Lsb_1_30', 'Lga_1_8', 'Ltor_1_8','Lsb_1_8']

    #SED PLOTTING
    out['realizations2plot'] = 10

    out['plotSEDrealizations'] = True

    return out

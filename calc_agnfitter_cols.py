import numpy as np
from astropy.table import Table, MaskedColumn, Column, join
#import LF_util
from radio_lf import util as LF_util  
import sys

#AGNsplit = 0.15
AGNsplit = 0.25
#AGNsplit_err = 0.2
AGNsplit2 = 0.15
AGNsplit1 = 0.35
AGNsplit_err = 0.2

#clobber = True

Slam = '_1_30_'



field = sys.argv[1]

if field=='all':
    fields = ['Bootes','EN1','LH']
else:
    fields = [field]


version = 'v1.0'

zpapply = 'atlas'


for field in fields:
    
    print (field)
    
    if field == 'EN1':
        infileradio = '/beegfs/lofar/deepfields/lgz/en1/final-v1.0.fits'
        infile = 'EN1_v1.0_agnfitter_out.fits'
        outfile = 'EN1_v1.0_agnfitter_out_calc.fits'
    elif field == 'LH':
        infileradio = '/beegfs/lofar/deepfields/lgz/lockman/final-v1.0.fits'
        infile = 'LH_v1.0_agnfitter_out.fits'
        outfile = 'LH_v1.0_agnfitter_out_calc.fits'
    elif field == 'Bootes':
        infileradio = '/beegfs/lofar/deepfields/lgz/bootes/final-v1.0.fits'
        infile = 'Bootes_v1.0_agnfitter_out.fits'
        outfile = 'Bootes_v1.0_agnfitter_out_calc.fits'
    else:
        print('not implemented:',field)
        sys.exit(1)


    agnfitterout = Table.read(infile)
    radiocat = Table.read(infileradio)

    cat = join(radiocat, agnfitterout, keys='Source_Name')

    print(len(agnfitterout), len(radiocat), len(cat))

    

    l10 = np.log(10.)

    # needed for qIR and sedSF which needs L150
    ## in W
    LIR = 1e-7*10**np.float64(cat['LIR_8_1000_med'])
    LIR_e_l = LIR*l10*np.float64(cat['LIR_8_1000_med']-cat['LIR_8_1000_p16'])
    LIR_e_u = LIR*l10*np.float64(cat['LIR_8_1000_p84']-cat['LIR_8_1000_med'])
    
    #Slam = '_1_100_'
    
    TO = 10**np.float64(cat['Ltor'+Slam+'med'])
    BB = 10**np.float64(cat['Lbb'+Slam+'med'])
    SB = 10**np.float64(cat['Lsb'+Slam+'med'])
    GA = 10**np.float64(cat['Lga'+Slam+'med'])  ##NB: the L_galaxy in the same freq regime as the torus and starburst
    TO_em = TO*l10*np.float64(cat['Ltor'+Slam+'med']-cat['Ltor'+Slam+'p16'])
    BB_em = TO*l10*np.float64(cat['Lbb'+Slam+'med']-cat['Lbb'+Slam+'p16'])
    SB_em = SB*l10*np.float64(cat['Lsb'+Slam+'med']-cat['Lsb'+Slam+'p16'])
    GA_em = GA*l10*np.float64(cat['Lga'+Slam+'med']-cat['Lga'+Slam+'p16'])
    TO_ep = TO*l10*np.float64(cat['Ltor'+Slam+'p84']-cat['Ltor'+Slam+'med'])
    BB_ep = TO*l10*np.float64(cat['Lbb'+Slam+'p84']-cat['Lbb'+Slam+'med'])
    SB_ep = SB*l10*np.float64(cat['Lsb'+Slam+'p84']-cat['Lsb'+Slam+'med'])
    GA_ep = GA*l10*np.float64(cat['Lga'+Slam+'p84']-cat['Lga'+Slam+'med'])
    
    
    fSB = SB/(GA + SB)
    
    #new fTO
    Tot =  GA + BB + TO
    fAGN = (TO + BB)/Tot
    
    fAGN_em = ((GA**2.*(TO_em**2.+BB_em**2.) + (TO+GA)**2.*GA_em**2.)**0.5)/(Tot**2.) 
    fAGN_ep = ((GA**2.*(TO_ep**2.+BB_ep**2.) + (TO+GA)**2.*GA_ep**2.)**0.5)/(Tot**2.) 
    
    
    ## errors #
    fSB_em = fSB**2.*(GA/SB)*np.sqrt( (SB_em/SB)**2. + (GA_em/GA)**2. )
    fSB_ep = fSB**2.*(GA/SB)*np.sqrt( (SB_ep/SB)**2. + (GA_ep/GA)**2. )
    
    Tot_em = np.sqrt( (TO_em/TO)**2. + (GA_em/GA)**2. +  (SB_em/SB)**2.)
    Tot_ep = np.sqrt( (TO_ep/TO)**2. + (GA_ep/GA)**2. +  (SB_ep/SB)**2.)
        
    ln_l = cat['neg_ln_like_ml']
    
    cat.add_column(MaskedColumn(data=Tot, name='Tot'))
    cat.add_column(MaskedColumn(data=Tot_em, name='Tot_em'))
    cat.add_column(MaskedColumn(data=Tot_ep, name='Tot_ep'))
    cat.add_column(MaskedColumn(data=fAGN, name='fAGN'))
    cat.add_column(MaskedColumn(data=fAGN_em, name='fAGN_em'))
    cat.add_column(MaskedColumn(data=fAGN_ep, name='fAGN_ep'))
    cat.add_column(MaskedColumn(data=fSB, name='fSB'))
    cat.add_column(MaskedColumn(data=fSB_em, name='fSB_em'))
    cat.add_column(MaskedColumn(data=fSB_ep, name='fSB_ep'))
    cat.add_column(MaskedColumn(data=ln_l, name='ln_l'))
    
    
    ff = cat['Total_flux']   #in Jy
    cat.remove_column('Z_BEST_2')
    cat.rename_column('Z_BEST_1', 'Z_BEST')
    zz = cat['Z_BEST']   #in Jy

    power = np.log10(LF_util.RadioPower(ff, zz, alpha=-0.7))  # assumed spec ind
    
    cat.add_column(MaskedColumn(data=power, name='logL150'))
    
    L150 = 10**power
    q150 =  np.log10((LIR / 3.75e12) / L150)
    cat.add_column(MaskedColumn(data=LIR, name='LIR'))
    cat.add_column(MaskedColumn(data=LIR_e_l, name='LIR_em'))
    cat.add_column(MaskedColumn(data=LIR_e_u, name='LIR_ep'))
    cat.add_column(MaskedColumn(data=q150, name='q150'))
    
    sigSF = 1.5
    sedSF = (q150 >= (1.72*(1+zz)**-0.22)-sigSF*0.529) & (power < 26.)  # Ivison+ 2010
    
    
    sedAGN = ~sedSF
    
    sedHERG = (fAGN >= AGNsplit) & ~sedSF
    sedLERG = (fAGN <  AGNsplit) & ~sedSF
    sedHERGc = (fAGN >= AGNsplit1) & ~sedSF
    sedLERGc = (fAGN <  AGNsplit2) & ~sedSF
    sedHERGu = (fAGN >= AGNsplit2) & ~sedSF
    sedLERGu = (fAGN <  AGNsplit1) & ~sedSF
    
    
    goodAF = (ln_l>-30)
    
    cat.add_column(MaskedColumn(data=goodAF, name='goodAF'))
    cat.add_column(MaskedColumn(data=sedSF, name='sedSF'))
    cat.add_column(MaskedColumn(data=sedAGN, name='sedAGN'))
    cat.add_column(MaskedColumn(data=sedHERG, name='sedHERG'))
    cat.add_column(MaskedColumn(data=sedHERGc, name='sedHERGc'))
    cat.add_column(MaskedColumn(data=sedHERGu, name='sedHERGu'))
    cat.add_column(MaskedColumn(data=sedLERG, name='sedLERG'))
    cat.add_column(MaskedColumn(data=sedLERGc, name='sedLERGc'))
    cat.add_column(MaskedColumn(data=sedLERGu, name='sedLERGu'))
    
    
    cat.write(outfile, overwrite=True)

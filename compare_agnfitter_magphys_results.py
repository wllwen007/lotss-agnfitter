import numpy as np
from radio_lf import util as LF_util  
from utils import plotting as pp
from astropy.table import Table, join
import matplotlib.pyplot as plt
import astropy.units as u

mcat = Table.read('/data2/wwilliams/projects/lofar_surveys/deep/magphys/ELAIS_N1_lowzmagphys_20200219.fits')
acat = Table.read('EN1_agnfitter_out_calc.fits')

mcat.rename_column('galaxy_id', 'Source_Name')

cat = join(mcat, acat)


badm = cat['chi2_99_flag']==1
bada = cat['neg_ln_like_med'] < -20

print(np.sum(badm), 'bad magphys fits')
print(np.sum(bada), 'bad agnfitter fits')
print(np.sum(badm & bada), 'bad magphys&agnfitter fits')
print(np.sum(badm | bada), 'bad magphys|agnfitter fits')

checkz = 0
if checkz:
    f,ax = pp.paper_single_ax()
    ax.scatter(cat['Z_BEST'], cat['z'],marker='.',s=4,alpha=0.2)
    pp.set_attrib(ax,xlabel='$z_{AF}$', ylabel='$z_{M}$')
    pp.fig_save_many(f,'plots/magtoaf_zs')  
    ## ok


f,ax = pp.paper_single_ax()
ax.scatter(cat['neg_ln_like_med'], cat['chi2'],marker='.',c='C0',s=4,alpha=0.2)
ax.scatter(cat['neg_ln_like_med'][badm], cat['chi2'][badm],marker='.',c='C1',s=4,alpha=0.4)
ax.scatter(cat['neg_ln_like_med'][bada], cat['chi2'][bada],marker='.',c='C2',s=4,alpha=0.4)
pp.set_attrib(ax,xlabel='$\log L_{AF}$', ylabel='$\chi_{M}$')
pp.fig_save_many(f,'plots/magtoaf_gofs')


cat_aLIR = np.log10(10**np.float64(cat['LIR_8_1000_med'])/u.Lsun.to('erg/s'))
f,ax = pp.paper_single_ax()
ax.scatter(cat_aLIR, np.log10(cat['Ldust_50']),marker='.',s=4,alpha=0.2, c='C0')
ax.scatter(cat_aLIR[badm], np.log10(cat['Ldust_50'][badm]),marker='.',s=4,alpha=0.2, c='C1')
ax.scatter(cat_aLIR[bada], np.log10(cat['Ldust_50'][bada]),marker='.',s=4,alpha=0.2, c='C2')
pp.set_attrib(ax,xlabel='$ L_{IR,AF}$ $(L_\odot)$', ylabel='$L_{dust,M}$ $(L_\odot)$')
pp.plot_equal(ax)
pp.fig_save_many(f,'plots/magtoaf_LIRs')


f,ax = pp.paper_single_ax()
ax.scatter(cat['z'], np.log10(cat['Ldust_50']),marker='.',s=4,alpha=0.2, c='C0')
ax.scatter(cat['z'][badm], np.log10(cat['Ldust_50'][badm]),marker='.',s=4,alpha=0.2, c='C1')
#ax.scatter(cat['z'][bada], np.log10(cat['Ldust_50'][bada]),marker='.',s=4,alpha=0.2, c='C2')
pp.set_attrib(ax,xlabel='$z$', ylabel='$L_{dust,M}$ $(L_\odot)$')
pp.fig_save_many(f,'plots/mag_LIR_z')



f,ax = pp.paper_single_ax()
ax.scatter(cat['log_Mstar_med'], np.log10(cat['Mstar_50']),marker='.',s=4,alpha=0.2, c='C0')
ax.scatter(cat['log_Mstar_med'][badm], np.log10(cat['Mstar_50'][badm]),marker='.',s=4,alpha=0.2, c='C1')
ax.scatter(cat['log_Mstar_med'][bada], np.log10(cat['Mstar_50'][bada]),marker='.',s=4,alpha=0.2, c='C2')
pp.set_attrib(ax,xlabel='$M_{AF}$ $(M_\odot)$', ylabel='$M_{M}$ $(M_\odot)$')
pp.plot_equal(ax)
pp.fig_save_many(f,'plots/magtoaf_masss')

l_msfr = np.log10(cat['SFR_50'])
l_asfr1 = np.log10(cat['SFR_opt_med'])
l_asfr2 = np.log10(cat['SFR_IR_med'])
f,ax = pp.paper_single_ax()
ax.scatter(l_asfr1, l_msfr,marker='.',s=4,alpha=0.2, c='C0')
ax.scatter(l_asfr1[badm], l_msfr[badm],marker='.',s=4,alpha=0.2, c='C1')
ax.scatter(l_asfr1[bada], l_msfr[bada],marker='.',s=4,alpha=0.2, c='C2')
pp.set_attrib(ax,xlabel='$\log SFR_{opt,AF}$', ylabel='$\log SFR_{M}$ $(M_\odot/yr)$')
pp.plot_equal(ax,zorder=-10)
pp.fig_save_many(f,'plots/magtoaf_sfrs')
f,ax = pp.paper_single_ax()
ax.scatter(l_asfr2, l_msfr,marker='.',s=4,alpha=0.2, c='C0')
ax.scatter(l_asfr2[badm], l_msfr[badm],marker='.',s=4,alpha=0.2, c='C1')
ax.scatter(l_asfr2[bada], l_msfr[bada],marker='.',s=4,alpha=0.2, c='C2')
pp.set_attrib(ax,xlabel='$\log SFR_{IR,AF}$', ylabel='$\log SFR_{M}$ $(M_\odot/yr)$')
pp.plot_equal(ax,zorder=-10)
pp.fig_save_many(f,'plots/magtoaf_sfrs')


sigSF = 1.5
zlim = np.arange(0,7,0.01)
sedSFlim = (1.72*(1+zlim)**-0.22)-sigSF*0.529  # Ivison+ 2010

zz = cat['Z_BEST']
ff = cat['Total_flux']
power = np.log10(LF_util.RadioPower(ff, zz, alpha=-0.7))  # assumed spec ind
L150 = 10**power
mLIR = (cat['Ldust_50']*u.Lsun.to('W')/ 3.75e12)
mq150 =  np.log10((mLIR ) / L150)
f,ax = pp.paper_single_ax()
ax.scatter(cat['Z_BEST'], mq150, marker='.',s=4,alpha=0.2)
ax.scatter(cat['Z_BEST'][badm], mq150[badm], marker='.',s=4,alpha=0.2,c='C1')
ax.plot(zlim, sedSFlim, 'gray')
pp.set_attrib(ax,xlabel='$z$', ylabel='$q_{{150}}$',ylim=(-4,3))
pp.fig_save_many(f,'plots/mag_z_q150')  


#aLIR = (((10**np.float64(cat['LIR_8_1000_med'])/3.75e12))*u.erg/u.s).to('W')
aLIR = 1e-7*10**np.float64(cat['LIR_8_1000_med'])/3.75e12
aq150 =  np.log10((aLIR ) / L150)
f,ax = pp.paper_single_ax()
ax.scatter(cat['Z_BEST'], aq150, marker='.',s=4,alpha=0.2)
ax.scatter(cat['Z_BEST'][bada], aq150[bada], marker='.',s=4,alpha=0.2,c='C2')
ax.plot(zlim, sedSFlim, 'gray')
pp.set_attrib(ax,xlabel='$z$', ylabel='$q_{{150}}$',ylim=(-4,3))
pp.fig_save_many(f,'plots/af_z_q150')  



#f,ax = pp.paper_single_ax()
#ax.scatter(cat['Z_BEST'], cat['q150'], marker='.',s=4,alpha=0.2)
##ax.scatter(cat['Z_BEST'], qq, marker='.',s=4,alpha=0.2)
#pp.set_attrib(ax,xlabel='$z$', ylabel='$q_{{150}}$',ylim=(-4,3))
#pp.fig_save_many(f,'plots/af_z_q150_b')  

#Mstar_50  # _best/_bayes

#f,ax = pp.paper_single_ax()
#ax.scatter(acat['Z_BEST'], acat['logL150'],marker='.',s=4,alpha=0.2)
#pp.set_attrib(ax,xlabel='$z$', ylabel='$\log L_{{150}}$')
#pp.fig_save_many(f,'plots/z_power')  

#f,ax = pp.paper_single_ax()
#ax.scatter(acat['Z_BEST'], acat['q150'],marker='.',s=4,alpha=0.2)
#pp.set_attrib(ax,xlabel='$z$', ylabel='$q_{{150}}$')
#pp.fig_save_many(f,'plots/z_q150')  

#f,ax = pp.paper_single_ax()
#ax.scatter(acat['Z_BEST'], np.log10(acat['LIR']),marker='.',s=4,alpha=0.2)
#pp.set_attrib(ax,xlabel='$z$', ylabel='$\log L_{{IR}}$')
#pp.fig_save_many(f,'plots/z_LIR')  

import numpy as np
from utils import plotting as pp
from astropy.table import Table
import matplotlib.pyplot as plt


cat = Table.read('EN1_agnfitter_out_calc.fits')

f,ax = pp.paper_single_ax()
ax.scatter(cat['Z_BEST'], cat['logL150'],marker='.',s=4,alpha=0.2)
pp.set_attrib(ax,xlabel='$z$', ylabel='$\log L_{{150}}$')
pp.fig_save_many(f,'plots/z_power')  

f,ax = pp.paper_single_ax()
ax.scatter(cat['Z_BEST'], cat['q150'],marker='.',s=4,alpha=0.2)
pp.set_attrib(ax,xlabel='$z$', ylabel='$q_{{150}}$')
pp.fig_save_many(f,'plots/z_q150')  

f,ax = pp.paper_single_ax()
ax.scatter(cat['Z_BEST'], np.log10(cat['LIR']),marker='.',s=4,alpha=0.2)
pp.set_attrib(ax,xlabel='$z$', ylabel='$\log L_{{IR}}$')
pp.fig_save_many(f,'plots/z_LIR')  


sedLERG = cat['sedLERG'] & cat['goodAF']
sedHERG = cat['sedHERG'] & cat['goodAF']
sedSF = cat['sedSF'] & cat['goodAF']

f,ax = pp.paper_single_ax()
ax.scatter(cat['fTO'], cat['fSB'],marker='.',s=4,alpha=0.2)
ax.scatter(cat['fTO'][sedSF], cat['fSB'][sedSF],marker='.',s=4,alpha=0.2,label='SF')
ax.scatter(cat['fTO'][sedLERG], cat['fSB'][sedLERG],marker='.',s=4,alpha=0.2,label='LERG')
ax.scatter(cat['fTO'][sedHERG], cat['fSB'][sedHERG],marker='.',s=4,alpha=0.2,label='HERG')
pp.set_attrib(ax,xlabel='$f_{TO}$', ylabel='$f_{SB}$')
pp.fig_save_many(f,'plots/af_ratios')  

ax.set_xscale('log')
ax.set_yscale('log')
pp.fig_save_many(f,'plots/af_ratios_log')  


from mocpy import MOC, World2ScreenMPL
from astropy.coordinates import Angle, SkyCoord
import astropy.units as u
import matplotlib.pyplot as plt
import sys

field = sys.argv[1]

if field.lower() == 'bootes':
    path = 'Bootes_FIR/'
    ra = 217.5
    dec = 34
    ext = '_bootes'
elif field.lower() == 'lh':
    path = 'Lockman_FIR/'
    ra = 161.25
    dec = 58
    ext = '_Lockman'
elif field.lower() == 'en1':
    path = 'ELAIS_N1_FIR_prelim/'
    ra = 243.75
    dec = 55
    ext = ''

# Load a MOC
mocmips = MOC.from_fits(path+'mocs/mips_moc{}.fits'.format(ext))
mocpacs = MOC.from_fits(path+'mocs/pacs_moc{}.fits'.format(ext))
mocspire = MOC.from_fits(path+'mocs/spire_moc{}.fits'.format(ext))

# Plot the MOC using matplotlib
fig = plt.figure(figsize=(15, 10))
# Define a astropy WCS easily
with World2ScreenMPL(fig, 
        fov=10 * u.deg,
        center=SkyCoord(ra, dec, unit='deg', frame='icrs'),
        coordsys="icrs",
        rotation=Angle(0, u.degree),
        projection="AIT") as wcs:
    ax = fig.add_subplot(1, 1, 1, projection=wcs)
    # Call fill with a matplotlib axe and the `~astropy.wcs.WCS` wcs object.
    mocmips.fill(ax=ax, wcs=wcs, alpha=0.5, fill=True, color="green", label='MIPS')
    mocmips.border(ax=ax, wcs=wcs, alpha=0.5, color="black")
    mocpacs.fill(ax=ax, wcs=wcs, alpha=0.5, fill=True, color="red", label='PACS')
    mocspire.fill(ax=ax, wcs=wcs, alpha=0.5, fill=True, color="blue", label='SPIRE')
plt.legend()
plt.xlabel('ra')
plt.ylabel('dec')
plt.title('Coverage of MIPS')
plt.grid(color="black", linestyle="dotted")
plt.savefig(path+'mocs/moc.png')

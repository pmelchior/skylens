#!/bin/env python

from os import system, popen, environ
import fitsio
import numpy as np
from glob import glob
from random import randint

def kappa2alpha(kappafile, alphafile):
    fits = fitsio.FITS(kappafile)
    kappa = fits[0].read()
    h = fits[0].read_header()
    fits.close()
    sidel = h['SIDEL']
    Dl = h['DL']
    Ds = h['DS']
    Dls = h['DLS']

    # smooth roll-off on the edge
    L = kappa.shape[0]
    roll = 10
    window = np.outer(np.sin(np.pi/2 * np.linspace(0,1,roll)), np.ones(L))
    kappa[0:roll, :] *= window
    kappa[:, 0:roll] *= window.T
    window = np.flipud(window)
    kappa[L-roll:, :] *= window
    kappa[:, L-roll:] *= window.T

    # mechanics of making gradient field in FFT
    mu = kappa.mean()
    kappa_ = np.fft.fft2(kappa)
    k = np.fft.fftfreq(L, d=sidel/L) # spacing now in arcsec
    Ka = np.add.outer(k, 1j*k)
    with np.errstate(divide="ignore", invalid="ignore"):
        alpha_ = kappa_ / Ka
    alpha_[Ka == 0] = 0
    alpha = np.fft.ifft2(alpha_);
    
    # since kappa has no-zero mean, we need to add a linear function
    #tilt = np.outer(np.ones(L), np.arange(L) - L/2)*mu
    #alpha.real += tilt
    #alpha.imag += tilt.T
    alpha *= 1. / 3600 / 180 * np.pi / Ds * Dls # now in radians

    # write alpha to FITS in SkyLens++ format
    fits = fitsio.FITS(alphafile, 'rw')
    sidel2 = h['SIDEL2']
    h['SIDEL'] = sidel2 # doesn't update comments. Be it...
    h['SIDEL2'] = sidel
    fits.write(np.array(np.real(alpha), dtype='float32'), header=h)
    fits.write(np.array(np.imag(alpha), dtype='float32'))
    fits.close()

def getMstar(z):
    data = np.loadtxt('mstar_z.cmbf.dat')
    return np.interp(z, data[:,0], data[:,1])

def getSeed():
    return randint(0, 4294967295) # whole unsiged int range

def getFieldSize(M):
    return 60*2*(M/1e14)**0.5 # arcsec

def createMainConfig(outdir):
    fp = open(outdir + "/shear_accuracy.conf", "w")
    fp.write("PROJECT\tS\tSHEAR_ACC\n")
    fp.write("OUTFILE\tS\tSHEAR_ACC.fits\n")
    fp.write("TELESCOPE\tS\ttelescope.conf\n")
    fp.write("FILTER\tS\ttelescopes/Subaru/I.fits\n")
    fp.write("EXPTIME\tI\t2000\n")
    fp.write("SKY\tD\t22.7\n")
    fp.write("ATMOSPHERE\tD\t0.\n")
    fp.write("AIRMASS\tD\t1.0\n")
    fp.write("LENSES\tVS\tmoka_lens.conf\n")
    fp.write("OVERSAMPLING\tI\t1\n")
    fp.write("NOISE\tI\t0\n")
    fp.close()

def createLensConfig(zl, mokafile, outdir):
    fp = open(outdir + "/moka_lens.conf", "w")
    fp.write("REDSHIFT\tD\t%.2f\t# lens redshift\n" % zl)
    fp.write("ANGLEFILE\tS\t%s\t# FITS file with deflection angles\n" %  mokafile)
    fp.close()
    

def createTelescopeConfig(M, outdir):
    fp = open(outdir + "/telescope.conf", "w")
    fp.write("DIAMETER\tD\t1.0\n")
    fp.write("FLAT-ACC\tD\t0.0\n")
    fp.write("GAIN\tD\t1\n")
    fp.write("PIXSIZE\tD\t0.05\n")
    fp.write("RON\tD\t10.0\n")
    fp.write("MIRROR\tD\t1\n")
    fp.write("OPTICS\tD\t1\n")
    fp.write("CCD\tD\t1\n")
    fp.write("FOV_X\tD\t%.1f\n" % getFieldSize(M))
    fp.write("FOV_Y\tD\t%.1f\n" % getFieldSize(M))
    fp.close()

def createMOKAInput(M, zl, zs, seed, outdir):
    fp = open(outdir + "/INPUT", "w")
    fp.write("\n!...Omega0\n0.3\n!...OmegaL\n0.7\n!...h0\n0.7\n!...w_dark_energy\n-1\n")
    fp.write("!...lens_redshift\n%.2f\n" % zl)
    fp.write("!...halo_mass\n%.4e\n" % M)
    fp.write("!...source_redshift\n%.2f\n" % zs)
    fp.write("!...field_of_view\n%.2f\n" % (getFieldSize(M)/60))
    fp.write("!...JAFFE_HERNQUIST_or_DM_for_the_BCG_or_not\nNO\n")
    fp.write("!...number_of_haloes\n1\n")
    fp.write("!...npixels\n2048\n")
    fp.write("!...write_fits_for_what_component\nNO\n")
    fp.write("!...write_SkyLens_file_FORTRAN_or_CPP\nCPP\n")
    fp.write("!...input_file_mass_concentration\nNO\n") 
    fp.write("!...scatter_in_concentration_log_normal\n0.25\n")
    fp.write("!...mass_resolution_subhalo_mass_function\n1.e10\n")
    fp.write("!...SISc_or_NFWc_for_satellitess_trunkated_SIS_and_NFW_not_trunkated\nNFWc\n")
    fp.write("!...satellite_distribution_model_NFW_or_GAO\nGAO\n")
    fp.write("!...mstar_zl\n%.5e\n" % getMstar(zl))
    fp.write("!..elliptical_or_spherical_halo\n1\n")
    fp.write("!...n_pix_smooth_fourier\n0\n")
    fp.write("!...ADC\nYES\n")
    fp.write("!...zero_padding_size\n2\n")
    fp.write("!..bcg_velocity_disp_prof_and_lensing_comps_prof\nNO\n")
    fp.write("!..RNG_seed\n%d" % seed)
    fp.close()

import sys
import tempfile
import sqlite3

if (len(sys.argv) < 4):
    print "usage: " + sys.argv[0] + " <Mass> <Redshift> <Repetitions>"
    exit(0)

M = float(sys.argv[1])
zl = float(sys.argv[2])
repetition = int(sys.argv[3])

for repeat in xrange(repetition):
    # create tempdir
    tmpdir = tempfile.mkdtemp(prefix='moka_batch_', dir='.')

    # run MOKA for halo with given parameters
    seed = getSeed()
    createMOKAInput(M, zl, 1.5, seed, tmpdir)
    print "running MOKA", M, zl, repeat
    popen("cd " + tmpdir + " && " + environ["MOKABIN"])
    mokafile = "/".join(glob(tmpdir + "/fits/0SkyLens*.fits")[0].split("/")[2:])
    # update SkyLens config files with new lens parameters
    createMainConfig(tmpdir)
    createLensConfig(zl, mokafile, tmpdir)
    createTelescopeConfig(M, tmpdir)

    print "running SkyLens++"
    # run range of source redshifts per lens
    for zs in [0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 2, 2.5, 3, 4]:
        if zl < zs:
            popen("cd " + tmpdir + " && ../../bin/shear_accuracy -o -c shear_accuracy.conf -N 1000 -z %1.2f -s %d > SHEAR_ACC_%1.2f.log" % (zs, seed, zs))

    # gather the outputs and delete tmpdir...
    


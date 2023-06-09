#!/usr/bin/python

import argparse
import numpy as np
from tqdm import tqdm

import matplotlib
import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.style.ROOT)
    
# in order: bb, tautau, mumu, cc, tt, gg, gammagamma, Zgamma, WW, ZZ
w_sm = {20:     sum((5.244E-4, 3.949E-5, 1.465E-7, 2.936E-5, 0.,       1.895E-5, 2.861E-8, 0.000E0, 3.560E-11, 3.100E-12)),
        30:     sum((7.709E-4, 6.086E-5, 2.198E-7, 3.921E-5, 0.,       1.522E-5, 9.586E-8, 0.000E0, 6.432E-10, 1.586E-10)),
        40:     sum((9.752E-4, 8.193E-5, 2.932E-7, 4.851E-5, 0.,       1.939E-5, 2.246E-7, 0.000E0, 5.186E-9, 1.407E-9)),
        50:     sum((1.163E-3, 1.028E-4, 3.663E-7, 5.738E-5, 0.,       2.879E-5, 4.380E-7, 0.000E0, 2.720E-8, 7.734E-9)),
        60:     sum((1.342E-3, 1.237E-4, 4.397E-7, 6.601E-5, 0.,       4.305E-5, 7.642E-7, 0.000E0, 1.103E-7, 3.129E-8)),
        70:     sum((1.516E-3, 1.445E-4, 5.129E-7, 7.442E-5, 0.,       6.318E-5, 1.237E-6, 0.000E0, 3.839E-7, 1.080E-7)),
        80:     sum((1.684E-3, 1.653E-4, 5.862E-7, 8.258E-5, 0.,       9.002E-5, 1.898E-6, 0.000E0, 1.227E-6, 3.257E-7)),
        90:     sum((1.849E-3, 1.860E-4, 6.592E-7, 9.063E-5, 0.,       1.245E-4, 2.805E-6, 0.000E0, 4.509E-6, 9.188E-7)),
        100:    sum((2.011E-3, 2.068E-4, 7.325E-7, 9.852E-5, 0.,       1.674E-4, 4.035E-6, 1.226E-7, 2.67E-5, 2.775E-6)),
        110:    sum((2.170E-3, 2.275E-4, 8.058E-7, 1.063E-4, 0.,       2.198E-4, 5.699E-6, 1.114E-6, 1.327E-4, 1.234E-5)),
        120:    sum((2.328E-3, 2.483E-4, 8.791E-7, 1.140E-4, 0.,       2.827E-4, 7.977E-6, 3.881E-6, 4.826E-4, 5.536E-5)),
        125:    sum((2.406E-3, 2.587E-4, 9.159E-7, 1.178E-4, 0.,       3.185E-4, 9.428E-6, 6.266E-6, 8.489E-4, 1.066E-4)),
        125.09: sum((2.407E-3, 2.589E-4, 9.166E-7, 1.179E-4, 0.,       3.191E-4, 9.460E-6, 6.313E-6, 8.573E-4, 1.078E-4)),
        130:    sum((2.483E-3, 2.690E-4, 9.526E-7, 1.216E-4, 0.,       3.572E-4, 1.116E-5, 9.555E-6, 1.442E-3, 1.933E-4)),
        140:    sum((2.637E-3, 2.898E-4, 1.025E-6, 1.290E-4, 0.,       4.445E-4, 1.584E-5, 2.007E-5, 3.963E-3, 5.531E-4)),
        150:    sum((2.787E-3, 3.105E-4, 1.099E-6, 1.364E-4, 0.,       5.459E-4, 2.347E-5, 3.998E-5, 1.159E-2, 1.403E-3)),
        160:    sum((2.938E-3, 3.313E-4, 1.172E-6, 1.438E-4, 0.,       6.629E-4, 4.336E-5, 9.587E-5, 7.105E-2, 3.402E-3)),
        170:    sum((3.088E-3, 3.521E-4, 1.246E-6, 1.511E-4, 0.,       7.975E-4, 5.812E-5, 1.518E-4, 3.466E-1, 8.791E-3)),
        180:    sum((3.233E-3, 3.728E-4, 1.319E-6, 1.582E-4, 0.,       9.509E-4, 6.441E-5, 1.864E-4, 5.603E-1, 3.694E-2)),
        190:    sum((3.376E-3, 3.935E-4, 1.392E-6, 1.652E-4, 0.,       1.125E-3, 7.011E-5, 2.195E-4, 7.758E-1, 2.093E-1)),
        200:    sum((3.519E-3, 4.144E-4, 1.466E-6, 1.723E-4, 0.,       1.325E-3, 7.535E-5, 2.512E-4, 1.013E0, 3.521E-1)),
        210:    sum((3.659E-3, 4.349E-4, 1.539E-6, 1.790E-4, 0.,       1.549E-3, 8.011E-5, 2.810E-4, 1.280E0, 4.879E-1)),
        220:    sum((3.798E-3, 4.557E-4, 1.612E-6, 1.859E-4, 0.,       1.804E-3, 8.450E-5, 3.093E-4, 1.584E0, 6.319E-1)),
        230:    sum((3.938E-3, 4.765E-4, 1.685E-6, 1.927E-4, 0.,       2.093E-3, 8.855E-5, 3.360E-4, 1.928E0, 7.904E-1)),
        240:    sum((4.636E-3, 5.657E-4, 2.000E-6, 2.269E-4, 0.,       2.754E-3, 1.050E-4, 4.108E-4, 2.317E0, 9.673E-1)),
        250 :   sum((4.214E-3, 5.178E-4, 1.832E-6, 2.062E-4, 0.,       2.792E-3, 9.572E-5, 3.849E-4, 2.750E0, 1.166E0)),
        260 :   sum((4.351E-3, 5.384E-4, 1.905E-6, 2.129E-4, 1.946E-7, 3.215E-3, 9.890E-5, 4.072E-4, 3.232E0, 1.386E0)),
        270 :   sum((4.487E-3, 5.592E-4, 1.978E-6, 2.195E-4, 1.235E-5, 3.699E-3, 1.018E-4, 4.282E-4, 3.765E0, 1.632E0)),
        280 :   sum((4.689E-3, 5.885E-4, 2.081E-6, 2.294E-4, 7.015E-5, 4.315E-3, 1.061E-4, 4.545E-4, 4.351E0, 1.904E0)),
        290 :   sum((4.756E-3, 6.008E-4, 2.125E-6, 2.328E-4, 2.247E-4, 4.894E-3, 1.071E-4, 4.665E-4, 4.992E0, 2.202E0)),
        300 :   sum((4.890E-3, 6.215E-4, 2.198E-6, 2.393E-4, 5.767E-4, 5.644E-3, 1.094E-4, 4.840E-4, 5.690E0, 2.529E0)),
        320 :   sum((5.155E-3, 6.629E-4, 2.345E-6, 2.522E-4, 2.862E-3, 7.597E-3, 1.137E-4, 5.157E-4, 7.264E0, 3.271E0)),
        350 :   sum((5.549E-3, 7.250E-4, 2.564E-6, 2.715E-4, 2.420E-1, 1.363E-2, 1.156E-4, 5.537E-4, 1.010E1, 4.623E0)),
        400 :   sum((6.195E-3, 8.288E-4, 2.931E-6, 3.032E-4, 4.315E0,  2.573E-2, 7.226E-5, 5.545E-4, 1.623E1, 7.574E0)),
        450 :   sum((6.829E-3, 9.322E-4, 3.297E-6, 3.342E-4, 8.701E0,  3.420E-2, 3.885E-5, 5.369E-4, 2.432E1, 1.150E1)),
        500 :   sum((7.452E-3, 1.036E-3, 3.663E-6, 3.648E-4, 1.278E1,  4.052E-2, 1.802E-5, 5.151E-4, 3.460E1, 1.654E1)),
        550 :   sum((8.070E-3, 1.139E-3, 4.029E-6, 3.948E-4, 1.651E1,  4.539E-2, 7.930E-6, 4.925E-4, 4.734E1, 2.280E1)),
        600 :   sum((8.677E-3, 1.243E-3, 4.395E-6, 4.245E-4, 1.994E1,  4.921E-2, 7.031E-6, 4.707E-4, 6.274E1, 3.039E1)),
        650 :   sum((9.281E-3, 1.347E-3, 4.764E-6, 4.542E-4, 2.314E1,  5.243E-2, 1.429E-5, 4.508E-4, 8.110E1, 3.945E1)),
        700 :   sum((9.874E-3, 1.450E-3, 5.128E-6, 4.831E-4, 2.613E1,  5.497E-2, 2.905E-5, 4.325E-4, 1.026E2, 5.010E1)),
        750 :   sum((1.046E-2, 1.554E-3, 5.493E-6, 5.118E-4, 2.898E1,  5.711E-2, 5.099E-5, 4.173E-4, 1.275E2, 6.245E1)),
        800 :   sum((1.105E-2, 1.658E-3, 5.861E-6, 5.406E-4, 3.169E1,  5.896E-2, 7.996E-5, 4.051E-4, 1.561E2, 7.663E1)),
        850 :   sum((1.163E-2, 1.761E-3, 6.229E-6, 5.691E-4, 3.430E1,  6.050E-2, 1.159E-4, 3.961E-4, 1.886E2, 9.276E1)),
        900 :   sum((1.221E-2, 1.865E-3, 6.594E-6, 5.971E-4, 3.681E1,  6.179E-2, 1.587E-4, 3.911E-4, 2.252E2, 1.110E2)),
        950 :   sum((1.278E-2, 1.969E-3, 6.962E-6, 6.250E-4, 3.926E1,  6.293E-2, 2.087E-4, 3.899E-4, 2.662E2, 1.314E2)),
        1000:   sum((1.334E-2, 2.072E-3, 7.326E-6, 6.528E-4, 4.163E1,  6.392E-2, 2.660E-4, 3.932E-4, 3.118E2, 1.541E2))}

def plot_ratio_vs_sin_lambdas(masses, lambdas112, sthetas):

    for mx in tqdm(masses):
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        ax.set_title("$M_{{X}}$={}GeV".format(mx))

        for lbd in tqdm(lambdas112):
            ratio = [calc_width(mx, st, lbd)/mx for st in sthetas]
            ax.plot(sthetas, ratio, "-o", label="$\lambda_{{112}}$={}GeV".format(lbd))
        
        ax.set_xlabel(r'$\sin( \theta)$')
        ax.set_ylabel(r'$\frac{\Gamma_{X}}{M_{X}}$')
        plt.legend(loc="upper left")
        plt.tight_layout()
        plt.savefig('plots/ratio_vs_sin_lambdas_{}.png'.format(str(mx).replace('.', 'p').replace('-', 'm')))

def plot_ratio_vs_sin_masses(masses, lambdas112, sthetas):

    for lbd in tqdm(lambdas112):
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        ax.set_title("$\lambda_{{112}}$={}GeV".format(lbd))

        for mx in masses:
            ratio = [calc_width(mx, st, lbd)/mx for st in sthetas]
            ax.plot(sthetas, ratio, "-o", label="$M_{{X}}$={}GeV".format(mx))
        
        ax.set_xlabel(r'$\sin( \theta)$')
        ax.set_ylabel(r'$\frac{\Gamma_{X}}{M_{X}}$')
        plt.legend(loc="upper left")
        plt.tight_layout()
        plt.savefig('plots/ratio_vs_sin_masses_{}.png'.format(str(lbd).replace('.', 'p').replace('-', 'm')))

def calc_width(mres, stheta, lambda112):
    """
    Calculate width of singlet scalar as given by https://arxiv.org/pdf/2010.00597.pdf
    in eqs. A.1 and A.3, where "1" refers to SM Higgs and "2" refers to the singlet scalar.
    """        
    m1 = 125. # GeV, recommended by the Yellow Report Section 4
    m2 = mres
    width_2_to_11 = (lambda112**2 * np.sqrt(1-4*m1**2/m2**2)) / (8*np.pi*m2)
    width_2 = stheta**2*w_sm[m2] + width_2_to_11
    return np.round(width_2, 6) 

def main(opt):
    sthetas = np.arange(0.,1.00001,.05)
    
    plot_ratio_vs_sin_lambdas((300, 600, 1000), np.arange(0,301,100), sthetas)
    plot_ratio_vs_sin_masses((250, 400, 700, 1000), np.arange(0,301,100), sthetas)

if __name__=='__main__':
    example = 'python ' + __file__
    parser = argparse.ArgumentParser(description="Plot width/mass ratios. \nExample: " + example)
    FLAGS = parser.parse_args()
    main(FLAGS)

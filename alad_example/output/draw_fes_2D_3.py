#!/usr/bin/env python3
import matplotlib
import os
import argparse
import numpy as np
from matplotlib.figure import figaspect
from scipy.interpolate import interp1d
import matplotlib.ticker as ticker
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
import pandas as pd
import matplotlib.pyplot as plt
# matplotlib.use("pgf")
plt.rcParams.update({
    "pgf.texsystem": "lualatex",
    "font.family": "FreeSans",  # use serif/main font for text elements
    "mathtext.fontset": "stix",
    "text.usetex": False,     # use inline math for ticks
    "pgf.rcfonts": False,    # don't setup fonts from rc parameters
    "axes.labelsize": 28,
    "axes.linewidth": 2.0,
    "font.size": 24,
    'axes.unicode_minus': False,
    "pgf.preamble": '\n'.join([
         "\\usepackage{units}",          # load additional packages
         "\\usepackage{metalogo}",
         "\\usepackage{unicode-math}",   # unicode math setup
         r"\setmathfont{MathJax_Math}",
         r"\setmainfont{FreeSans}",  # serif font via preamble
         ])
})


# parser = argparse.ArgumentParser()
# parser.add_argument("pmf", help = "specify the PMF file")
# parser.add_argument("-o", "--output", default = "output.png", help = "specify the PNG output image file")
# parser.add_argument("--xtitle", default = "CV1", help = "title along X axis")
# parser.add_argument("--ytitle", default = "CV2", help = "title along Y axis")
# parser.add_argument("--levels", default = 25, type = int, help = "number of levels")
# args = parser.parse_args()

def plotfes(pmffilename, pngfilename, xtitle, ytitle, title, level=31):
    x, y, z = np.genfromtxt(pmffilename, unpack=True)
    z = np.clip(z, 0, 16.0)
    binx = len(set(x))
    biny = len(set(y))
    xi = x.reshape(binx, biny)
    yi = y.reshape(binx, biny)
    zi = z.reshape(binx, biny)
    fig = plt.figure()
    cf = plt.contourf(xi, yi, zi, np.linspace(0, 15.0, level), cmap='turbo')
    plt.xlabel(xtitle)
    plt.ylabel(ytitle)
    plt.title(title)
    ax = plt.gca()
    ax.set_xlim(-180, 180)
    ax.set_ylim(-180, 180)
    ax.tick_params(direction='in', which='major', length=6.0, width=1.0, top=True, right=True, pad=6.0)
    ax.tick_params(direction='in', which='minor', length=3.0, width=1.0, top=True, right=True, pad=6.0)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.xaxis.set_major_locator(plt.MultipleLocator(90))
    ax.yaxis.set_major_locator(plt.MultipleLocator(90))
    clb = plt.colorbar(cf)
    clb.ax.set_title("kcal/mol", fontsize=20)
    clb.ax.xaxis.get_major_formatter()._usetex = False
    clb.ax.yaxis.get_major_formatter()._usetex = False
    # clbticks = [float(i.get_position()[1]) for i in clb.ax.get_yticklabels()]
    # clbticksstr = ['{:.1f}'.format(i) for i in clbticks]
    # print(clbticksstr)
    # clb.ax.set_yticklabels(clbticksstr, fontsize = 20)
    plt.savefig(pngfilename, dpi=400, bbox_inches='tight', transparent=False)
    plt.close(fig)
    return

# plotfes(args.pmf, args.output, xtitle = args.xtitle, ytitle = args.ytitle, level = args.levels)
#plotfes('trialanine_apath_cpp+1.czar.pmf', 'apath_2d.png', r'$\xi (s)$', r'$\xi (z)$')
#plotfes('trialanine_gpath_2d.czar.pmf', 'gpath_2d.png', r'$\xi (s)$', r'$\xi (z)$')


if __name__ == '__main__':
    stride = 1
    for t in range(0+stride, 32+stride, stride):
        pmffilename = f'out_{t:03d}.0.pmf'
        pngfilename = f'out_{t:03d}.png'
        title = f'{t} ns'
        plotfes(pmffilename, pngfilename, r'$\phi$ (°)', r'$\psi$ (°)', title)
    # plotfes('out_final.pmf', 'out_final.png', r'$\phi$ (°)', r'$\psi$ (°)', 'NAMD')

#!/usr/bin/env python3
from histogram import HistogramScalar
from histogram import Axis
from reweight import convert_probability_to_pmf
import numpy as np
import math

if __name__ == '__main__':
    ax1 = Axis(-180.0, 180.0, 72, True)
    ax2 = Axis(-180.0, 180.0, 72, True)
    hist = HistogramScalar(ax=[ax1, ax2])
    kbt = 300.0 * 0.0019872041
    with open('alad_opes.colvars.opes_metad1.misc.traj', 'r') as f_data:
        for line in f_data:
            if line.startswith('#'):
                continue
            fields = line.split()
            if len(fields) == 8:
                time_ns = int(round(float(fields[0]))) / 1e3
                # print(time_ns)
                phi = float(fields[1])
                psi = float(fields[2])
                bias = float(fields[3])
                prob = np.exp(bias / kbt)
                hist[[phi, psi]] += prob
                if time_ns > 0 and math.isclose(time_ns % 1, 0):
                    hist_pmf = convert_probability_to_pmf(hist, kbt)
                    # hist_pmf.data /= 4.184
                    with open(f'to_ffmpeg/out_{int(time_ns):03d}.pmf', 'w') as f_pmf:
                        hist_pmf.write_to_stream(f_pmf)
    with open('out_final.pmf', 'w') as f_pmf:
        hist_pmf = convert_probability_to_pmf(hist, kbt)
        # hist_pmf.data /= 4.184
        hist_pmf.write_to_stream(f_pmf)

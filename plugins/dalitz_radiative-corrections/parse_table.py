#!/usr/bin/env python3

import re
import numpy as np

'''
The intention of this script is to parse the ascii files provided in the source of the paper arxiv:1711.11001
and dump them in a formatted way to a file which can be included via the TGraph2D constructor within the
PDalitzCorrections Plugin. The correction factors are applied as a weight to the pseudo-scalar meson Dalitz decays.
'''

def read_file(path):
    decimal = re.compile(r"\d?\.\d*")
    with open(path) as f:
        corrections, corrs = [], {}
        for line in f:
            if line.startswith('#') and 'y = ' in line:
                y = decimal.search(line).group(0)
                y = float(y)
                corrs[y] = []
            elif not line.startswith('#') and line.strip():
                vals = line.split()
                x = float(vals[0])
                vals = [vals[i] for i in [3, 4, 8]]
                corr = sum(map(float, vals))
                corrs[y].append(x)
                corrections.append((x, y, corr))

    print('parsed y values:', corrs.keys())
    print('individual data points per y value:')
    for y in corrs.keys():
        print('{:.2f}: {:3d}'.format(y, len(corrs[y])))

    return corrections

def write_table(corrections, path):
    with open(path, 'w+') as f:
        for val in corrections:
            f.write('%g %g %g\n' % val)

def write_header(corrections, path, name='corrections'):
    with open(path, 'w+') as f:
        f.write('#include <vector>\n\n')
        f.write('struct correction {\n')
        f.write('\tdouble x;\n')
        f.write('\tdouble y;\n')
        f.write('\tdouble c;\n')
        f.write('};\n\n')
        f.write('static const std::vector<correction> %s = {' % name)
        struct = "{%g, %g, %g}"
        f.write(', '.join(map(lambda x : struct % x, corrections)))
        f.write('};\n')
    print('C++ header written to file:', path)

def main():
    # eta -> e+ e- g
    corr = read_file('/home/sascha/1711.11001/anc/eta_e')
    write_header(corr, 'eta_dilepton_radiative_corrections.h', 'eta_ee_corrections')
    # eta -> mu+ mu- g
    corr = read_file('/home/sascha/1711.11001/anc/eta_mu')
    write_header(corr, 'eta_dimuon_radiative_corrections.h', 'eta_mumu_corrections')
    # eta' -> e+ e- g
    corr = read_file('/home/sascha/1711.11001/anc/etap_e')
    write_header(corr, 'etap_dilepton_radiative_corrections.h', 'etap_ee_corrections')
    # eta' -> mu+ mu- g
    corr = read_file('/home/sascha/1711.11001/anc/etap_mu')
    write_header(corr, 'etap_dimuon_radiative_corrections.h', 'etap_mumu_corrections')

if __name__ == '__main__':
    main()

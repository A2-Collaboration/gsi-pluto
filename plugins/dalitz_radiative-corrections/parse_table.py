#!/usr/bin/env python3

import re
import numpy as np

'''
The intention of this script is to parse the ascii tables provided in the source of the paper arxiv:1711.11001
and dump them in a formatted way to a header file which can be included in the Pluto PDalitzCorrections plugin
and apply the correction factors from it as an additional weight to the pseudo-scalar meson Dalitz decays.
'''

def read_table(path):
    data = np.loadtxt(path)

def read_file(path):
    decimal = re.compile(r"\d?\.\d*")
    with open(path) as f:
        corrections = {}
        for line in f:
            if line.startswith('#') and 'y = ' in line:
                y = decimal.search(line).group(0)
                y = float(y)
                corrections[y] = []
            elif not line.startswith('#') and line.strip():
                vals = line.split()
                x = float(vals[0])
                vals = [vals[i] for i in [3, 4, 8]]
                corr = sum(map(float, vals))
                corrections[y].append((x, corr))

    print('parsed y values:', corrections.keys())
    print('individual data points per y value:')
    for y in corrections.keys():
        print('{:.2f}: {:3d}'.format(y, len(corrections[y])))

    return corrections

def write_header(corrections, path, namespace=''):
    y_vals, corrs = [], []
    # there are 999 x values for y = 0, but this makes things much more complicated
    # for an easier implementation reduce the 999 values to just the 66 used for every other value
    # therefore pick the x values from y = 0.1 and just store the corrections for the exsiting x values
    x_vals, _ = zip(*corrections[0.1])
    for y in corrections:
        if y == 0.:
            _corr = [c for x, c in corrections[y] if x in x_vals]
            if len(_corr) is not len(x_vals):
                print('Something went wrong, only %d corrections found for %d reference x values' % (len(_corr), len(x_vals)))
                return
        else:
            _x, _corr = zip(*corrections[y])
        y_vals.append(y)
        if not x_vals:
            x_vals = _x
        corrs.extend(_corr)
    with open(path, 'w+') as f:
        f.write('#include <vector>\n\n')
        indent = ''
        if namespace:
            indent = ' '*4
            f.write('namespace %s {\n' % namespace)
        f.write('%sstatic const std::vector<double> x_values = {%s};\n'
                % (indent, ', '.join(map(str, x_vals))))
        f.write('%sstatic const std::vector<double> y_values = {%s};\n'
                % (indent, ', '.join(map(str, y_vals))))
        f.write('%sstatic const std::vector<double> corr_values = {%s};\n'
                % (indent, ', '.join(map(str, corrs))))
        if namespace:
            f.write('}')
    print('C++ header written to file:', path)

def main():
    corr = read_file('/home/sascha/1711.11001/anc/etap_e')
    write_header(corr, 'radiative_corrections_etap.h', 'etap')

if __name__ == '__main__':
    main()

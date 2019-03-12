#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
from datetime import datetime
from nefis.dataset import Nefis
import tables
import numpy as np
import libwaq


NEFIS_INFO = True  # PRINT NEFIS INFO

# CONSTANTS
atom = tables.Float32Atom()
fname = '../VossoEutro.def'
fname = '../VossoEutroCI2.def'
fnamelga = '../com.lga'
fnamecco = '../com.cco'

# Open lga and cco to get grid coordinates
ind, parsLGA = libwaq.readLGA(fnamelga) # ind[m, n]
X, Y, parsCCO = libwaq.readCCO(fnamecco)  # X[m, n], Y[m, n]

dt = datetime.now()
print('Processing coordinates:')
# Get and Verify Coordinates
m, n = ind.shape
k = parsLGA['kmax']

# To Do: Assert that all m and n from different sources are equal
print("shape: (m={}, n={}, z={})".format(m, n, k))

# Transform Y matrix in a Y vector (This only work for rectangular grid)
y = np.zeros(n)
for i in range(m):
    vec = Y[i, :]
    y[vec > 0] = vec[vec > 0]

# Transform X matrix in a X vector (This only work for rectangular grid)
x = np.zeros(m)
for j in range(n):
    vec = X[:, j]
    x[vec > 0] = vec[vec > 0]

coords = {'x': x, 'y': y}

coordmn = np.full((parsLGA['nmnw'], 2), np.nan, dtype=np.int32)
for i in range(m):
    for j in range(n):
        if ind[i, j] > 0:
            coordmn[ind[i, j]-1, 0] = i
            coordmn[ind[i, j]-1, 1] = j
print('Time Elapsed: ', datetime.now()-dt)

print('Get nefis info')
dt = datetime.now()
# Open WAQ Nefis
nf = Nefis(fname)

# Print NEFIS infos
if NEFIS_INFO:
    print(nf.groups)
    for key, item in nf.variables.items():
        info = item.flat()
        print(item.group, key, info['attributes']['description'])
        print('\tshape: ', info['shape'])
        print('\tunits: ', info['attributes']['units'])
exit()
# Get Name of Substances
varname = nf.get_data("DELWAQ_PARAMS", "SUBST_NAMES")
sizes = nf.get_data("DELWAQ_PARAMS", "SIZES")
varnames = [varname[i*20:i*20+20].strip() for i in range(sizes[0])]
print(varnames, sizes, len(varnames))
print('Time Elapsed: ', datetime.now()-dt)

# Choose Vars
# vars = [b'NO3', b'PO4']#, b'Chlo', b'cTR1', b'dTR1']
vars = [sys.argv[1].encode('UTF-8')]
id = varnames.index(vars[0])
print(id)
print(vars)
print('Converting to Hdf:')
varnames = vars
for nvar, varname in enumerate(varnames):
    nt = nf.groups['DELWAQ_RESULTS']['group_size']
    varname = "SUBST_{:03d}".format(id+1)
    nloc = nf.variables[varname].shape[0]
    fnameout = vars[nvar].decode('utf-8') + 'CI2.h5'
    print(fnameout)
    # print('NT: {}'.format(nt))
    # print('shape: {}'.format(nf.variables[varname].shape))
    d0 = datetime.now()
    with tables.open_file(fnameout, mode='w', title=vars[nvar].decode('utf8')) as f:
        f.create_array(f.root, "x", coords['x'], "X Coordinate")
        f.create_array(f.root, "y", coords['y'], "Y Coordinate")
        array_e = f.create_earray(
            f.root, vars[nvar].decode('utf8'), atom=atom, shape=(0, k, m, n))
        for t in range(nt):
            print('{}/{}\r'.format(t, nt), sep='', end='', flush=True)
            tmp = nf.get_data(
                'DELWAQ_RESULTS', varname, t=t)
            noseg = int(len(tmp)/k)
            matrix = np.full((1, k, m, n), np.nan, dtype=np.float32)
            for z in range(k):
                for seg in range(noseg):
                    i, j = coordmn[seg]
                    matrix[0, k-z-1, i, j] = tmp[z*noseg + seg]
            array_e.append(matrix)
    print(datetime.now() - d0)

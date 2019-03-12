#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Delft3d-Flow results to hdf5
"""
from datetime import datetime
import numpy as np
import tables
from nefis.dataset import Nefis
import matplotlib.pyplot as plt


def getx(matrix):
    m, n = matrix.shape
    x = np.zeros(m)
    for j in range(n):
        vec = matrix[:, j]
        mask = ~np.isnan(vec)
        x[mask] = vec[mask]
    return x


def gety(matrix):
    m, n = matrix.shape
    y = np.zeros(n)
    for i in range(m):
        vec = matrix[i, :]
        mask = ~np.isnan(vec)
        y[mask] = vec[mask]
    return y


# CONSTANTS
atom = tables.Float32Atom()
fname = '../trim-vossoroca20zVert.def'

# What to do?
DO_COORD = True
DO_SPATIAL = True
DO_TIMESERIES = False

# Nefis
nf = Nefis(fname)

# Print nefis infos
print(nf.groups)
for key, item in nf.variables.items():
    info = item.flat()
    print(item.group, key, info['attributes']['description'])
    print('\tshape: ', info['shape'])
    print('\tunits: ', info['attributes']['units'])
exit()
# Get Coordinates
if DO_COORD == True:
    varnames = ['ZK']
    for varname in varnames:
        fname = varname + '.h5'
        print(fname)
        var = nf.variables[varname]
        group = var.group
        shape = var.shape
        output = nf.get_data(group, varname)
        # print(group, shape)
        # print(output.shape)
        with tables.open_file(fname, mode='w', title=varname) as f:
            f.create_array(f.root, varname, output, varname)
    varnames = ['XCOR', 'YCOR']
    for varname in varnames:
        fname = varname + '.h5'
        print(fname)
        var = nf.variables[varname]
        group = var.group
        shape = var.shape
        output = nf.get_data(group, varname)
        output[output <= 0] = np.nan
        if varname == 'XCOR':
            vecx = getx(output)
        if varname == 'YCOR':
            vecy = gety(output)
        with tables.open_file(fname, mode='w', title=varname) as f:
            f.create_array(f.root, varname, output, varname)
# Common Output Variables
varnames = ['R1', 'U1', 'V1', 'WPHY', 'VICUV']
varnames = ['S1']
varnames = ['DICWW', 'R1', 'U1', 'V1', 'WPHY', 'VICUV']
varnames = ['WINDU', 'WINDV']
exit()
if DO_SPATIAL == True:
    for varname in varnames:
        print(varname)
        nt = nf.groups['map-series']['group_size']
        varshape = nf.variables[varname].shape
        lenshape = len(varshape)
        if lenshape == 4:
            m, n, z, nv = nf.variables[varname].shape
            shape0 = (0, z, n, m)
            shape1 = (1, z, n, m)
            shape2 = shape3 = (0, ...)
        elif lenshape == 3:
            m, n, z = nf.variables[varname].shape
            shape0 = (0, z, n, m)
            shape1 = (1, z, n, m)
            shape2 = (0, ...)
            shape3 = (...)
        elif lenshape == 2:
            m, n = nf.variables[varname].shape
            shape0 = (0, n, m)
            shape1 = (1, n, m)
            shape2 = (0, ...)
            shape3 = (...)
        else:
            raise ValueError

        d0 = datetime.now()

        # Write matrix in Order for Spatial Analysis
        fnameout = varname + '.h5'
        with tables.open_file(fnameout, mode='w', title=varname) as f:
            f.create_array(f.root, "x", vecx, "X Coordinate")
            f.create_array(f.root, "y", vecy, "Y Coordinate")
            array_e = f.create_earray(
                f.root, varname, atom=atom, shape=shape0)
            for i in range(nt):
                print(str(i)+'\r', sep='', end='', flush=True)
                output = np.empty(shape1, dtype=np.float32)
                output[shape2] = nf.get_data(
                    'map-series', varname, t=i)[shape3]
                array_e.append(output)
        print('spatial: ', datetime.now() - d0)

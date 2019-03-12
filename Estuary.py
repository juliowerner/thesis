#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Estimate Horizontal and Vertical Estuary
"""
import numpy as np
import tables
import matplotlib.pyplot as plt
atom = tables.Float32Atom()

fnameU = 'HorMagVel.h5'
fnameD = 'VICUV.h5'
fnameout = 'Estuary.h5'
varname = "Estuary"
with tables.open_file(fnameU, mode='r') as fu, \
        tables.open_file(fnameD, mode='r') as fd, \
        tables.open_file(fnameout, mode='w', title=varname) as fout:
    print(fu)
    print(fd)
    nt, z, n, m = fu.root.MHV.shape
    array_e = fout.create_earray(
            fout.root, varname, atom=atom,
            shape=(0, z, n, m), expectedrows=nt)
    for i in range(nt):
        print(str(i)+'\r', sep='', end='', flush=True)
        U = fu.root.MHV[i, :, :, :]
        U = np.ma.array(U, mask=U == -999.0)
        D = fd.root.VICUV[i, :, :, :]
        D = np.ma.array(D, mask=D == -999.0)
        mag = D/(U*U)
        output = np.empty((1, z, n, m), dtype=np.float32)
        output[0, :, :, :] = mag
        array_e.append(output)


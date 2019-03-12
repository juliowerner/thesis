#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
"""
import numpy as np
import tables
from datetime import datetime as dt

atom = tables.Float32Atom()

Clim = 0.01 # mg/l

fname = "./PO4CI2.h5"
varname = "PO4"
varnameout = 'PO4E1'
with tables.open_file(fname, mode='r+') as f:
    try:
        f.remove_node('/'+varnameout)
    except tables.exceptions.NoSuchNodeError:
        pass
    tbl1 = f.get_node(f.root, 'PO4SWP')
    z, n, m, nt = tbl1.shape
    output = np.empty((z, n, m), dtype=np.float32)
    print('total:', z*n*m)
    for k in range(z):
        print(k)
        d0 = dt.now()
        tmp = tbl1[k, :, :, :]
        for j in range(n):
            for i in range(m):
                output[k, j, i] = np.sum(tmp[j, i, :] >= Clim)/nt
        print(dt.now()-d0)
    array = f.create_array(f.root, varnameout, output)

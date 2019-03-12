#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
"""
import numpy as np
import tables
from datetime import datetime as dt

atom = tables.Float32Atom()

fname = "./PO4.h5"
varname = "PO4"
d0 = dt.now()
with tables.open_file(fname, mode='r+') as f:
    tbl1  = f.get_node(f.root, varname)
    nt, z, n, m = tbl1.shape
    array_e = f.create_earray(
            f.root, varname+'SWP', atom=atom,
            shape=(z, n, m, 0), expectedrows=nt)
    output = np.empty((z, n, m, 1), dtype=np.float32)
    for i in range(nt):
        print(str(i)+'\r', sep='', end='', flush=True)
        t = tbl1[i,:,:,:]
        output[:,:,:,0] = t
        array_e.append(output)
print(dt.now()-d0)

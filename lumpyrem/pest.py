from lumpyrem import lumprem, run
import os
import pandas as pd

class Pst():
    """ A Pest setup class. Facilities generating PEST control, template and instruction files from an ennsemble of LUMPREM models.
    """
    def __init__(self, controlfile='temp.pst', models=None):
        self.controlfile = controlfile
        self.models = models


def write_pst(controlfile, models, params='all'):
    # get number of parameters
    param_grp = 12
    num_param = 12*len(models)
    for m in models:
        if type(m.irrigfile) == tuple:
            num_param += 1
            if param_grp < 13:
                param_grp += 1
        if type(m.vegfile) == tuple:
            num_param += 2
            if param_grp < 14:
                param_grp += 2

    with open(controlfile, 'w+') as f:
        f.write('pcf $\n')
        f.write("* control data\n")
        f.write("restart estimation\n")
        f.write("{0}\t{1}\t{2}\t{3}\t{4}{5}".format(num_param, 0, param_grp,0,0,'\n'))
        f.write("{0}\t{1}{2}".format(len(models), "0  single  point  1   0   0",'\n'))
        f.write("10.0  -3.0  0.3  0.03  10  999\n")
        f.write("10.0   10.0    0.001\n")
        f.write("0.1   boundscale\n")
        f.write("50  0.005  4  4  0.005  4\n")
        f.write("1  1  1\n")
        f.write("* singular value decomposition\n")
        f.write("1\n")
        f.write("10000  5.0e-7\n")
        f.write("0\n")
        f.write("* parameter groups\n")
        f.write("maxvol     relative  0.015   0.00005    switch   2.0  parabolic\n")
        f.write("irigvf     relative  0.015   0.0001     switch   2.0  parabolic\n")
        f.write("rdelay     relative  0.015   0.0001     switch   2.0  parabolic\n")
        f.write("mdelay     relative  0.015   0.0001     switch   2.0  parabolic\n")
        f.write("ks         relative  0.015   0.0        switch   2.0  parabolic\n")
        f.write("m          relative  0.015   0.0        switch   2.0  parabolic\n")
        f.write("l          relative  0.015   0.0        switch   2.0  parabolic\n")
        f.write("mfmax      relative  0.015   0.0001     switch   2.0  parabolic\n")
        f.write("offset     absolute  0.1     0.0        switch   2.0  parabolic\n")
        f.write("f1         relative  0.015   0.0        switch   2.0  parabolic\n")
        f.write("f2         relative  0.015   0.0        switch   2.0  parabolic\n")
        f.write("power      absolute  0.015   0.0        switch   2.0  parabolic\n")
        if param_grp >= 13:
            f.write("gwirfr     relative  0.015   0.0001     switch   2.0  parabolic\n")
        if param_grp >= 14:
            f.write("crfac      relative  0.015   0.0001     switch   2.0  parabolic\n")
            f.write("gamma      relative  0.015   0.0        switch   2.0  parabolic\n")
        f.write("* parameter data\n")
        for m in models:
            name = m.lumprem_model_name
            f.write("{0:<10}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}{10}".format("maxvol_"+name, 'log', 'factor', m.maxvol ,0.001,10.0,'maxvol',1.0,0.0,1,'\n'))
            f.write("{0:<10}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}{10}".format("irigvf_"+name, 'log', 'factor', m.irrigvolfrac ,0.001,1.0,'irigvf',1.0,0.0,1,'\n'))
            f.write("{0:<10}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}{10}".format("rdelay_"+name, 'log', 'factor', m.rdelay ,0.001,50.0,'rdelay',1.0,0.0,1,'\n'))
            f.write("{0:<10}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}{10}".format("mdelay_"+name, 'log', 'factor', m.mdelay ,0.001,10.0,'mdelay',1.0,0.0,1,'\n'))
            f.write("{0:<10}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}{10}".format("ks_"+name, 'log', 'factor', m.mdelay ,1e-5,100.0,'ks',1.0,0.0,1,'\n'))
            f.write("{0:<10}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}{10}".format("m_"+name, 'log', 'factor', m.M ,1e-2,10.0,'m',1.0,0.0,1,'\n'))
            f.write("{0:<10}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}{10}".format("l_"+name, 'log', 'factor', m.L ,0.1,1.0,'l',1.0,0.0,1,'\n'))
            f.write("{0:<10}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}{10}".format("mfmax_"+name, 'log', 'factor', m.mflowmax ,0.1,1.0,'mfmax',1.0,0.0,1,'\n'))
            f.write("{0:<10}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}{10}".format("offset_"+name, 'none', 'relative', m.offset ,-1000.0,10000.0,'offset',1.0,0.0,1,'\n'))
            f.write("{0:<10}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}{10}".format("f1_"+name, 'log', 'factor', m.mflowmax ,1e-5,1e4,'fac1',1.0,0.0,1,'\n'))
            f.write("{0:<10}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}{10}".format("f2_"+name, 'log', 'factor', m.mflowmax ,1e-5,1e4,'fac2',1.0,0.0,1,'\n'))
            f.write("{0:<10}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}{10}".format("power_"+name, 'none', 'relative', m.power ,-3.0,3.0,'power',1.0,0.0,1,'\n'))
            if param_grp >= 13:
                f.write("{0:<10}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}{10}".format("gwirfr_"+name, 'log', 'factor', m.irrigvolfrac ,0.001,1.0,'gwirfr',1.0,0.0,1,'\n'))
            if param_grp >= 14:
                f.write("{0:<10}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}{10}".format("crfac_"+name, 'log', 'factor', m.vegfile[0] ,0.001,3.0,'crfac',1.0,0.0,1,'\n'))
                f.write("{0:<10}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}{10}".format("gamma_"+name, 'log', 'factor', m.vegfile[1] ,0.001,100.0,'gamma',1.0,0.0,1,'\n'))
        f.write("* model input/output\n")
        for m in models:
            f.write("lr_"+m.lumprem_model_name+".tpl  lr_"+m.lumprem_model_name+".in\n")
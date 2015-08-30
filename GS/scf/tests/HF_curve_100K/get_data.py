#!/usr/bin/env python

import os
import re
import subprocess
import shutil
import time
import numpy as np

input_template_file = "input_template_file"
maindir = os.getcwd()

r_vec = np.arange(0.5,3.0,.1)

# substitute multiple keywords
def multigsub(subs,str):
    for k,v in subs.items():
        str = re.sub(k,v,str)
    return str



os.chdir(maindir)

summary = open(maindir + '/summary.txt', "w+")
summary.write('r   SCF   FTHF\n')

for r_val in r_vec:

    r_str = "%1.5f" % r_val
    job_name = r_str
    direc = maindir + "/" + job_name
    file_open = open(direc + "/" + "output.dat")
    summary.write('%s' % job_name)
    rhf = []
    fthfa = []
    fail  = []
    for line in file_open:
        if(re.search('CALCULATION FAILED', line)):
            fail.append(':-(')
        else:
            if(re.search('@RHF Final Energy', line)):
                print 'RHF match'
                rhf = re.findall(r'[+-]\d+.\d+', line)
                rhf.append(rhf[0])
            if(re.search('FTHFEnergy', line)):
                print 'FTHF match'
                fthf = re.findall(r'[+-]\d+.\d+', line)
                fthfa.append(fthf[0])
    if(not len(fail)):
        summary.write('  %s   %s\n' % (rhf[-1], fthfa[-1]))
    else:
        summary.write(' %s \n' % fail[0])
    
    

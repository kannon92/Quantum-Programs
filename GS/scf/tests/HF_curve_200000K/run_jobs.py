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

dir_search = os.listdir(maindir)
for dirs in dir_search:
    if os.path.isdir(dirs):
        subprocess.call('rm -rf %s' % dirs, shell=True)

# substitute multiple keywords
def multigsub(subs,str):
    for k,v in subs.items():
        str = re.sub(k,v,str)
    return str



for r_val in r_vec:

    r_str = "%1.5f" % r_val
    job_name = r_str
    
    if os.path.exists(job_name):
    
        subprocess.call('rm -rf %s' % job_name, shell=True)
        os.mkdir(job_name)
    else:
        os.mkdir(job_name)
    

    input_template = open(input_template_file,"r").read()
    subs = {"Revals" : job_name}
    print subs
    input = multigsub(subs,input_template)
    os.chdir(maindir + "/" + job_name)
    input_fname = "input.dat"
    input_file = open(input_fname,"w+")
    input_file.write(input)
    input_file.close()

            #write script file
    subprocess.call('psi4', shell=True)
    time.sleep(0.01)
    os.chdir(maindir)
os.chdir(maindir)


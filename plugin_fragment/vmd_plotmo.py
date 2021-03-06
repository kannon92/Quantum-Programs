#!/usr/bin/env python

import sys
import re
import subprocess
#import multigsub

cubefile = sys.argv[1]
def multigsub(subs,str):
    for k,v in subs.items():
        str = re.sub(k,v,str)
    return str
# Check if the user wants to reverse the phase of the MOs
isovalue = 0.02
if (len(sys.argv) == 3):
    isovalue = float(sys.argv[2])

vmd_exe = "/Applications/VMD\ 1.9.2.app/Contents/MacOS/startup.command.csh"
#vmd_exe = "/Volumes/VMD-1.9.1/VMD\ 1.9.1.app/Contents/MacOS/startup.command.csh"

vmd_script_name = "vmd_mo_script.vmd"

vmd_template = """#
# VMD script to plot MOs from cube files
#

# Load the molecule and change the atom style/colors
mol load cube CUBEFILE
mol modcolor 0 0 Element
mol modstyle 0 0 Licorice 0.150000 50.000000 100.000000
#mol modstyle 0 0 CPK 1.100000 0.600000 50.000000 17.000000
color Element C silver
color Element H white
material change ambient Opaque 0.310000
material change diffuse Opaque 0.720000
material change specular Opaque 0.500000
material change shininess Opaque 0.480000
material change opacity Opaque 1.000000
material change outline Opaque 0.000000
material change outlinewidth Opaque 0.000000
material change transmode Opaque 0.000000

material change ambient EdgyShiny 0.310000
material change diffuse EdgyShiny 0.720000
material change specular Opaque 0.750000
material change shininess EdgyShiny 1.0000

light 1 off
light 0 rot y  30.0
light 0 rot x -30.0
#rotate z by 0
#rotate x by 30
#rotate y by 90

# Eliminate the axis and perfect the view
axes location Off
display projection Orthographic
display depthcue off
color Display Background white"""

vmd_template_surface = """#
# Add the surfaces
mol color ColorID 3
mol representation Isosurface ISOVALUE1 0 0 0 1 1
mol selection all
mol material Glass2
mol addrep 0
mol color ColorID 7
mol representation Isosurface ISOVALUE2 0 0 0 1 1
mol selection all
mol material Glass2
mol addrep 0

scale by 0.80000

# Save a TGA snapshot of the scene
# Optionally make the surfaces transparent
render snapshot CUBEFILE.tga /usr/bin/open %s
mol modmaterial 2 0 Transparent
mol modmaterial 1 0 Transparent

quit
"""

vmd_script_surface = multigsub({"ISOVALUE1" : "%.5f" % isovalue,"ISOVALUE2" : "%.5f" % (-isovalue),"CUBEFILE" : cubefile},vmd_template_surface)

vmd_script = open(vmd_script_name,"w+")

vmd_template = re.sub("CUBEFILE",cubefile,vmd_template)
vmd_script.write(vmd_template + "\n" + vmd_script_surface)

vmd_script.close()

subprocess.call(("%s -e %s" % (vmd_exe,vmd_script_name)), shell=True)


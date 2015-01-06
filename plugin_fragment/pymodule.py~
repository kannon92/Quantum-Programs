import psi4
import re
import os
import inputparser
import math
import warnings
from driver import *
from wrappers import *
from molutil import *
import p4util
from p4xcpt import *


def run_plugin_fragment(name, **kwargs):
    r"""Function encoding sequence of PSI module and plugin calls so that
    plugin_fragment can be called via :py:func:`~driver.energy`. For post-scf plugins.

    >>> energy('plugin_fragment')

    """
    lowername = name.lower()
    kwargs = p4util.kwargs_lower(kwargs)

    # Your plugin's psi4 run sequence goes here
    psi4.plugin('plugin_fragment.so')

# Integration with driver routines
procedures['energy']['plugin_fragment'] = run_plugin_fragment


def localize():
    psi4.plugin('plugin_fragment.so')
    pass

import psi4
import re
import os
import math
import warnings
import aliases
from driver import *
from wrappers import *
from molutil import *


def run_prop(name, **kwargs):
    r"""Function encoding sequence of PSI module and plugin calls so that
    blah can be called via :py:func:`~driver.energy`.

    >>> energy('blah')

    """
    returnvalue = psi4.plugin('/Users/kevinhannon/Programs/PsiPlugins_psi4public/plugin_prop/plugin_prop.so')

    return returnvalue


# Integration with driver routines
procedures['energy']['prop'] = run_prop
psi4.plugin_load('/Users/kevinhannon/Programs/Quantum_Programs/Quantum-Programs/plugin_prop/plugin_prop.so')


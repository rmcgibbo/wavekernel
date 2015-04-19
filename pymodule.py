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


def sample_descriptors(n=100):
    psi4.set_local_option('WAVEKERNEL', 'MODE', 'SAMPLE_DESCRIPTORS')
    psi4.set_local_option('WAVEKERNEL', 'NUM_SAMPLE_DESCRIPTORS', n)
    psi4.plugin('wavekernel.so')

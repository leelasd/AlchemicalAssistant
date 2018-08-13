from __future__ import division, print_function

import sys

# OpenMM Imports
import simtk.openmm as mm
import simtk.openmm.app as app

# ParmEd Imports
from parmed.charmm import CharmmPsfFile, CharmmCrdFile, CharmmParameterSet
from parmed.openmm import StateDataReporter
from parmed import unit as u

# Load the CHARMM files
print('Loading CHARMM files...')
params = CharmmParameterSet('BNZ_2_PRD.prm')
ala5_gas = CharmmPsfFile('zero.psf')
ala5_crds = app.PDBFile('A2B_gas.pdb')

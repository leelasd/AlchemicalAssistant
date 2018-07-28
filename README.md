# My project's README

Alchemical Assistant is based on 
(LigParGen)Ligand Parameter Generator Based on 
Jorgensen group's OPLS-AA/CM1A(-LBCC) FF
Created on Mon Feb 15 15:40:05 2016

@author: Leela Sriram Dodda leela.dodda@yale.edu
@author: Matt R. Robinson leela.dodda@yale.edu
@author: William L. Jorgensen Lab 

##FF formats provided : 
--------------------
CHARMM/NAMD  .prm & .rtf  
GROMACS      .itp & .gro 
BOSS/MCPRO   .z
TINKER       .key & .xyz

##Input Files supported : 
--------------------

BOSS Z-matrix

################################################## 

if using Z-matrix file 
Usage: LigParGen -A INI.z -B FIN.z  


## REQUIREMENTS:

BOSS (need to set BOSSdir in bashrc and cshrc)
Preferably Anaconda python with following modules
pandas 
argparse
numpy
openbabel
rdkit

Please cite following references: 

1. LigParGen web server: an automatic OPLS-AA parameter generator for organic ligands  
   Leela S. Dodda  Israel Cabeza de Vaca  Julian Tirado-Rives William L. Jorgensen 
   Nucleic Acids Research, Volume 45, Issue W1, 3 July 2017, Pages W331–W336
   
2. 1.14**CM1A-LBCC: Localized Bond-Charge Corrected CM1A Charges for Condensed-Phase Simulations
   Leela S. Dodda, Jonah Z. Vilseck, Julian Tirado-Rives , and William L. Jorgensen 
   Department of Chemistry, Yale University, New Haven, Connecticut 06520-8107, United States
   J. Phys. Chem. B, 2017, 121 (15), pp 3864–3870
   
3. Accuracy of free energies of hydration using CM1 and CM3 atomic charges.
   Udier–Blagović, M., Morales De Tirado, P., Pearlman, S. A. and Jorgensen, W. L. 
   J. Comput. Chem., 2004, 25,1322–1332. doi:10.1002/jcc.20059


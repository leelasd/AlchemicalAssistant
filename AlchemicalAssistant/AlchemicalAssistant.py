import argparse
import os
from AlchemicalAssistant.MakeFEPZmat import Alchemify


#def main(**kwargs):
def main():
    parser = argparse.ArgumentParser(
        prog='AlchemicalAssistant',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""
	SCRIPT TO CREATE FEP TOPOLOGIES FROM LigParGen generated 
	BOSS ZMATRIX FOR GMX, NAMD, BOSS and MCPRO
	Created on Mon Feb 15 15:40:05 2016
	@author: Leela S. Dodda leela.dodda@yale.edu
	@author: William L. Jorgensen Lab 
	
	Usage: python Alchemist.py -A BNZ.z -B PRD.z 
        THIS CREATES a BNZ_2_PRD.z Z-Matrix 

	REQUIREMENTS:
	BOSS (need to set BOSSdir in bashrc and cshrc)
	Preferably Anaconda python with following modules
	pandas 
	argparse
	numpy
	RDKit for MCCS Algorithm  
        openbabel 
	"""
    )
    parser.add_argument(
        "-A", "--Astate", help="Paste SMILES code from CHEMSPIDER or PubChem", type=str)
    parser.add_argument(
        "-B", "--Bstate", help="Submit MOL file from CHEMSPIDER or PubChem", type=str)
    args = parser.parse_args()
    ref_z = args.Astate
    ali_z = args.Bstate
    Alchemify(ref_z=ref_z, ali_z=ali_z)
    os.system('/bin/rm sum out plt.pdb log optzmat AMOL.p BMOL.p aligned.mol LLN')

if __name__ == "__main__":
    main()

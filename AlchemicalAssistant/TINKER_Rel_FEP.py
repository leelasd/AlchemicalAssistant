import pickle 
import pandas as pd
from Vector_algebra import AtomNum2Mass 
from FEPBOSSReader import tor_cent

def xyz_prep(atoms, coos, bonds, f_df, resid='A2B',pdbname='COMBO'):
    opdb = open(pdbname+'_tinker.xyz', 'w+')
    opdb.write('%6d %s LigParGen generated OPLS-AA/CM1A Parameters\n'%(len(atoms),resid))
    for i,r in f_df.iterrows():
        atom_blist = list(bonds[bonds.cl1==i+3]['cl2'])+list(bonds[bonds.cl2==i+3]['cl1'])
        line_atom_blist = ''.join(['%6d'%(l-2) for l in atom_blist])
        opdb.write('%6d%3s%14.6f%12.6f%12.6f%6d%s\n'%(i+1,r.SYM,r.X,r.Y,r.Z,r.OPLSN,line_atom_blist))
    opdb.close()
    return None

def tinker_prm(complete_atom_df,all_dfs, pdbname='COMBO'):
#    at_df = pd.read_csv('csv_tinker.csv') 
    at_df = complete_atom_df
    df_bnds = all_dfs['BONDS']
    ang_df = all_dfs['ANGLE']
    tor_df = all_dfs['DIHED']
    uid2onum = at_df.set_index('UID')['OPLSN'].to_dict()
    prm = open(pdbname + '.key', 'w+')
    prm.write(
'''

      ##############################
      ##                          ##
      ##  Force Field Definition  ##
      ##                          ##
      ##############################


forcefield              OPLS-AA

vdwindex                TYPE
vdwtype                 LENNARD-JONES
radiusrule              GEOMETRIC
radiustype              SIGMA
radiussize              DIAMETER
epsilonrule             GEOMETRIC
torsionunit             1.0
imptorunit              1.0
vdw-14-scale            2.0
chg-14-scale            2.0
electric                332.06
dielectric              1.0


      #############################
      ##                         ##
      ##  Atom Type Definitions  ##
      ##                         ##
      #############################


''')

    for i,r in at_df.iterrows():
        numbs = list(df_bnds[df_bnds.cl1==i+3]['cl2'])+list(df_bnds[df_bnds.cl2==i+3]['cl1'] )
        prm.write('atom %10d %4d %5s %8s %10d %10.3f %5d\n' %
            (r.OPLSN,r.OPLSN,r.at_ty,'\"'+r.at_symb+'\"',r.at_num,AtomNum2Mass(r.at_num),len(numbs)))
    prm.write(
'''


      ################################
      ##                            ##
      ##  Van der Waals Parameters  ##
      ##                            ##
      ################################


''')
    for i,r in at_df.iterrows():
        prm.write('vdw %11d %16.4f %8.4f \n' %
                (r.OPLSN, r.SIG, r.EPS))
    prm.write(
'''


      ##################################
      ##                              ##
      ##  Bond Stretching Parameters  ##
      ##                              ##
      ##################################


'''
)
    # ask about this one
    for index, row in df_bnds.iterrows():
        atom1_type = uid2onum[row.Ncl1]
        atom2_type = uid2onum[row.Ncl2]
        R = row['RIJ']
        K = row['KIJ']

        prm.write('bond %10d %4d %16.2f %8.4f \n' %
                    (atom1_type, atom2_type, K, R))

    prm.write(
'''


      ################################
      ##                            ##
      ##  Angle Bending Parameters  ##
      ##                            ##
      ################################


''')
    for index, row in ang_df.iterrows():
        atom1_type = uid2onum[row['Ncl1']]
        atom2_type = uid2onum[row['Ncl2']]
        atom3_type = uid2onum[row['Ncl3']]
#        R = float(row['R'])
#        K = float(row['K'])

        prm.write('angle %9d %4d %4d %8.2f %8.2f \n' %
                    (atom1_type, atom2_type, atom3_type, row.K, row.R))
    prm.write(
'''


      ################################
      ##                            ##
      ##   Urey-Bradley Parameters  ##
      ##                            ##
      ################################


ureybrad     35   34   35      38.25     1.5537


      #####################################
      ##                                 ##
      ##  Improper Torsional Parameters  ##
      ##                                 ##
      #####################################



''')
    for index, row in tor_df.iterrows():
        if row['IMP'] == True:
            cen_nums = tor_cent([row.I,row.J,row.K,row.L],list(df_bnds.UID))
#            atom1_type = int(num2typ2symb[cen_nums[1]][1][5:])#int[int(row['I'])][1].strip('_opls'))
            atom1_type = at_df.ix[cen_nums[1]-3]['OPLSN'] 
#            atom2_type = int(num2typ2symb[cen_nums[2]][1][5:])#int[int(row['I'])][1].strip('_opls'))
            atom2_type = at_df.ix[cen_nums[2]-3]['OPLSN'] 
#            atom3_central_type = int(num2typ2symb[cen_nums[0]][1][5:]) #int(types[int(row['J'])][1].strip('_opls'))
            atom3_central_type = at_df.ix[cen_nums[0]-3]['OPLSN'] 
            #atom4_type = int(num2typ2symb[cen_nums[3]][1][5:])
            atom4_type =  at_df.ix[cen_nums[3]-3]['OPLSN'] 

            V2 = float(row['V2'])
            gamma = 180.0
            n = 2

            # ordering for this is weird
            # see https://ryanmrichard.github.io/ForceManII/tinkerformat.html
            prm.write('imptors %7d %4d %4d %4d %12.3f %4.1f %2d \n' %
                    (atom1_type, atom2_type, atom3_central_type, atom4_type, V2, gamma, n))
    prm.write(
'''


      ############################
      ##                        ##
      ##  Torsional Parameters  ##
      ##                        ##
      ############################


''')
    for index, row in tor_df.iterrows():
        if row['IMP'] == False:
            atom1_type = at_df.ix[row.I-3]['OPLSN']#int(types[int(row['I'])][1].strip('_opls'))
            atom2_type = at_df.ix[row.J-3]['OPLSN']#int(types[int(row['J'])][1].strip('_opls'))
            atom3_type = at_df.ix[row.K-3]['OPLSN']#int(types[int(row['K'])][1].strip('_opls'))
            atom4_type = at_df.ix[row.L-3]['OPLSN']#int(types[int(row['L'])][1].strip('_opls'))

            V1 = float(row['V1'])
            gamma1 = 0.0
            n1 = 1

            V2 = float(row['V2'])
            gamma2 = 180.0
            n2 = 2

            V3 = float(row['V3'])
            gamma3 = 0.0
            n3 = 3

            prm.write('torsion %7d %4d %4d %4d %12.3f %4.1f %2d %6.3f %4.1f %2d %6.3f %4.1f %2d \n' %
                    (atom1_type, atom2_type, atom3_type, atom4_type, V1, gamma1, n1, V2, gamma2, n2, V3, gamma3, n3))
    prm.write(
'''
torsion       0    0    0    0        0.000  0.0  1  0.000 180.0  2  0.000  0.0  3

      ########################################
      ##                                    ##
      ##  Atomic Partial Charge Parameters  ##
      ##                                    ##
      ########################################


''')
#    types_idx = 0
    for i,r in at_df.iterrows():
        prm.write('charge %11d %16.4f \n' %
                (r.OPLSN, r.Q))

    prm.close()
    return None 

#all_dfs = pickle.load(open("BNZ_2_PRD.p", "rb")) 
#tinker_prm(all_dfs, resid='A2B')

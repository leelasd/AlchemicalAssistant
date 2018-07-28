import pandas as pd
import networkx as nx
import numpy as np
from MolReaders import impr_id, tor_id
from Vector_algebra import Distance, angle, dihedral, pairing_func, ang_id
import collections


def CENT2IMPty(cent):
    ce2ty = {'C': 160, 'N': 161, 'CA': 162, 'CM': 221}
#    print(cent)
    if cent in list(ce2ty.keys()):
        return ce2ty[cent]
    elif cent[0] == 'C': return 162
    #elif cent[0] == 'N': return 161 
    else:    
        return 0


def Get_Add_Int(mol_icords, Z_BONDS, Z_ANGLES, Z_TORSIONS):
    all_bonds_mol, all_angles_mol, all_torsions_mol = mol_icords[
        'BONDS'], mol_icords['ANGLES'], mol_icords['TORSIONS']
    Z_B = {pairing_func(i[0] - 3, i[1] - 3): [i[0] - 3, i[1] - 3]
           for i in Z_BONDS.values()}
    Z_A = {ang_id([i[0] - 3, i[1] - 3, i[2] - 3]): [i[0] - 3,
                                                    i[1] - 3, i[2] - 3] for i in Z_ANGLES.values()}
    Z_T = {tor_id([i[0] - 3, i[1] - 3, i[2] - 3, i[3] - 3],False,list(Z_B.keys())): [i[0] - 3,
                                                              i[1] - 3, i[2] - 3, i[3] - 3] for i in Z_TORSIONS.values()}
    Z_Ad_B, Z_Ad_A, Z_Ad_T = collections.OrderedDict(
    ), collections.OrderedDict(), collections.OrderedDict()
    Z_Ad_I = collections.OrderedDict()
    for b_ij in all_bonds_mol.values():
        uid_b_ij = pairing_func(b_ij[0], b_ij[1])
        if uid_b_ij not in list(Z_B.keys()):
            Z_Ad_B[uid_b_ij] = [b_ij[0], b_ij[1]]
    for a_ij in all_angles_mol.keys():
        if a_ij not in list(Z_A.keys()):
            Z_Ad_A[a_ij] = [i for i in all_angles_mol[a_ij]]
    for t_ij in all_torsions_mol.keys():
        if t_ij not in list(Z_T.keys()):
            Z_Ad_T[t_ij] = [i + 3 for i in all_torsions_mol[t_ij]]
    for c in mol_icords['IMPROPERS'].values():
        Z_Ad_I[impr_id([c[0], c[2], c[3]])] = [i for i in c]
    return(Z_Ad_B, Z_Ad_A, Z_Ad_T, Z_Ad_I)


def print_FEPZMAT(atoms, G_mol, mol_icords, coos, umatchA, umatchB, extra, G_vars, num2typ, zmat_name='TEST_A2B.z', resid='A2B'):
    Z_ATOMS = {1: 'X', 2: 'X'}
    Z_INI = {1: -1, 2: -1}
    Z_FIN = {1: -1, 2: -1}
    Z_BONDS = {1: (1, 0, 0.000), 2: (2, 1, 1.00), 3: (3, 2, 1.00)}
    Z_ANGLES = {1: (1, 0, 0, 0.000), 2: (2, 1, 0, 0.000),
                3: (3, 2, 1, 90.00), 4: (4, 3, 2, 90.0)}
    Z_TORSIONS = {1: (1, 0, 0, 0, 0.00), 2: (2, 1, 0, 0, 0.00), 3: (
        3, 2, 1, 0, 0.00), 4: (4, 3, 2, 1, 0.00), 5: (5, 4, 3, 2, 90.0)}
    for i in range(len(atoms)):
        Z_ATOMS[i + 3] = atoms[i]
    for i in range(len(atoms)):
        Z_INI[i + 3] = G_mol.node[i]['init_no']
        Z_FIN[i + 3] = G_mol.node[i]['finl_no']
        #print(i,Z_INI[i + 3],Z_FIN[i + 3])
    ######## collect all the bonds ######################
    n_ats = 0
    B_LINK = {}
    for i in range(1, len(G_mol.nodes())):
        # if n_ats > 0:
        neigs = np.sort(list(G_mol.neighbors(i)))
        B_LINK[i] = neigs[0]
        Z_BONDS[i + 3] = (i + 3, neigs[0] + 3,
                          Distance(coos[i], coos[neigs[0]]))
        n_ats += 1
    ######## Collect all the angles #####################
    n_ats = 0
    A_LINK = {}
    for i in G_mol.nodes():
        if n_ats > 1:
            neigs = np.sort(list(G_mol.neighbors(B_LINK[i])))
            neigs = neigs[neigs < i]
            if i > (len(coos) - len(umatchB) - 1):
                neigs = np.array([j for j in neigs if j not in umatchA])
            A_LINK[i] = neigs[0]
            ang = angle(coos[i], coos[B_LINK[i]], coos[neigs[0]])
            Z_ANGLES[i + 3] = (i + 3, B_LINK[i] + 3, neigs[0] + 3, ang)
        n_ats += 1
    ######## Collect all the IC torsions/impropers ######
    n_ats = 0
    for i in G_mol.nodes():
        if n_ats > 2:
            neigs = list(G_mol.neighbors(A_LINK[i]))
            neigs = np.array(
                [j for j in neigs if j not in [i, B_LINK[i], A_LINK[i]]])
            neigs = np.sort(neigs)
            neigs = neigs[neigs < i]
            if len(neigs) < 1:
                neigs = [j for j in list(G_mol.neighbors(
                    B_LINK[i])) if j not in [i, A_LINK[i]]]
                if impr_id([i, A_LINK[i], neigs[0]]) in mol_icords['IMPROPERS'].keys():
                    del mol_icords['IMPROPERS'][
                        impr_id([i, A_LINK[i], neigs[0]])]
            [ti, tj, tk, tl] = [i, B_LINK[i], A_LINK[i], neigs[0]]
            dihed = dihedral(coos[ti], coos[tj], coos[tk], coos[tl])
            Z_TORSIONS[i + 3] = (ti + 3, tj + 3, tk + 3, tl + 3, dihed)
        n_ats += 1
######## Removing the repeated IC using Additional ICs ##############
    UID_2_IB = G_vars.set_index('UID')['IB'].to_dict()
    UID_2_FB = G_vars.set_index('UID')['FB'].to_dict()

    GG_V = {i: [UID_2_IB[pairing_func(Z_BONDS[i][0] - 3, Z_BONDS[i][1] - 3)], UID_2_FB[pairing_func(
        Z_BONDS[i][0] - 3, Z_BONDS[i][1] - 3)]] for i in range(4, len(G_mol.nodes()) + 3)}
    Z_Ad_B, Z_Ad_A, Z_Ad_T, Z_Ad_I = Get_Add_Int(
        mol_icords, Z_BONDS, Z_ANGLES, Z_TORSIONS)
   #Z_BONDS = {1: (1, 0, 0.000), 2: (2, 1, 1.00), 3: (3, 2, 1.00)}
    GG_V[1] = [0.0, 0.0]
    GG_V[2] = [1.0, 1.0]
    GG_V[3] = [1.0, 1.0]
# PRINTING ACTUAL Z-MATRIX
    ofile = open(zmat_name, 'w+')
    ofile.write('MCPRO/BOSS FEP Z-Matrix by LSD\n')
    for i in range(1, len(atoms) + 3):
        ofile.write('%4d %-3s%5d%5d%5d%12.6f%4d%12.6f%4d%12.6f%4s%5d\n'
                    % (i, Z_ATOMS[i], Z_INI[i], Z_FIN[i], Z_BONDS[i][1], GG_V[i][0], Z_ANGLES[i][-2], Z_ANGLES[i][-1], Z_TORSIONS[i][-2], Z_TORSIONS[i][-1], resid[0:3], 1)
                    )
    ofile.write(
        '                    Geometry Variations follow    (2I4,F12.6)\n')
    for k in range(4, len(G_mol.nodes()) + 3):
        ofile.write('%4d%4d%12.6f\n' % (k, 1, GG_V[k][1]))
    ofile.write('                     Variable Bonds follow         (I4)\n')
    for i in range(4, len(atoms) + 3):
        ofile.write('%4d\n' % i)
    ofile.write('                    Additional Bonds follow       (2I4)\n')
    if len(Z_Ad_B) > 0:
        for i in Z_Ad_B.values():
            ofile.write('%4d%4d\n' % (i[0] + 3, i[1] + 3))
    # CREATE A FUNCTION TO DEFINE ADDITIONAL BONDS IN CASE OF RINGS
    ofile.write('''                    Harmonic Constraints follow   (2I4,4F10.4)
                    Variable Bond Angles follow   (I4)\n''')
    for i in range(5, len(atoms) + 3):
        ofile.write('%4d\n' % i)
    ofile.write('                    Additional Bond Angles follow (3I4)\n')
    if len(Z_Ad_A) > 0:
        for i in Z_Ad_A.values():
            ofile.write('%4d%4d%4d\n' % (i[0] + 3, i[1] + 3, i[2] + 3))
    # CREATE A FUNCTION TO DEFINE ADDITIONAL BONDS IN CASE OF RINGS
    ofile.write(
        '                    Variable Dihedrals follow     (3I4,F12.6)\n')
    for i in range(6, len(atoms) + 3):
        ofile.write('%4d%4d%4d%12.6f\n' % (i, -1, -1, 0.000))
    ofile.write('                    Additional Dihedrals follow   (6I4)\n')
    if len(Z_Ad_T) > 0:
        for k in Z_Ad_T.keys():
            torsion = Z_Ad_T[k]
            ofile.write('%4d%4d%4d%4d%4d%4d\n' %
                        (torsion[0], torsion[1], torsion[2], torsion[3], -1, -1))
    add_imp_lines = []
    if len(Z_Ad_I) > 0:
        for k in Z_Ad_I.keys():
            imp = Z_Ad_I[k]
            ui = Z_INI[imp[1] + 3]
            uf = Z_FIN[imp[1] + 3]
            imp_ity = CENT2IMPty(num2typ[ui])
            imp_fty = CENT2IMPty(num2typ[uf])
            #print(num2typ[ui],num2typ[uf],imp_ity,imp_fty)
            #print(ui,uf,imp_ity,imp_fty)
            add_imp_lines.append('%4d%4d%4d%4d%4d%4d\n' % (
                imp[0] + 3, imp[1] + 3, imp[2] + 3, imp[3] + 3, imp_ity, imp_fty))
#            ofile.write('%4d%4d%4d%4d%4d%4d\n'%(imp[0]+3,imp[1]+3,imp[2]+3,imp[3]+3,imp_ity,imp_fty))
    ofile.write('                    Domain Definitions follow     (4I4)\n')
    for i in umatchA:
        for j in umatchB:
            ofile.write('%4d%4d%4d%4d\n' % (i + 3, i + 3, j + 3, j + 3))
    ofile.write(
        '''                    Conformational Search (2I4,2F12.6)
                    Local Heating Residues follow (I4 or I4-I4)
                    Final blank line\n
 Final Non-Bonded Parameters for QM (AM1 CM1Ax1.14) Atoms:\n
''')
    for i in extra:
        ofile.write('%s' % i)
    ofile.close()
    return add_imp_lines


def GetGeomVars(ALL_BNDS, A_UID, B_UID, A_UID_RIJ, B_UID_RIJ, FEP_ITY, FEP_FTY):
    geom_var = pd.DataFrame(ALL_BNDS).T
    geom_var.columns = ['BI', 'BJ']
    geom_var['UID'] = [pairing_func(i, j)
                       for i, j in zip(geom_var.BI, geom_var.BJ)]
    geom_var['INI_BT'] = ['%s-%s' %
                          (FEP_ITY[i], FEP_ITY[j]) for i, j in zip(geom_var.BI, geom_var.BJ)]
    geom_var['FIN_BT'] = ['%s-%s' %
                          (FEP_FTY[i], FEP_FTY[j]) for i, j in zip(geom_var.BI, geom_var.BJ)]
    test_A = []
    test_B = []
    BL_A = []
    BL_B = []
    # print(geom_var)
    for rb in geom_var.iterrows():
        i, r = rb
        if r.UID in A_UID:
            test_A.append(1)
            BL_A.append(A_UID_RIJ[r.UID])
        if r.UID in B_UID:
            test_B.append(1)
            BL_B.append(B_UID_RIJ[r.UID])
        if (r.UID not in A_UID):
            test_A.append(0)
            if ('100' in r.INI_BT) or ('4894' in r.INI_BT):
                BL_A.append(0.300)
            else:
                BL_A.append(B_UID_RIJ[i])
        if (r.UID not in B_UID):
            test_B.append(0)
            if ('100' in r.FIN_BT) or ('4894' in r.FIN_BT):
                BL_B.append(0.300)
            else:
                BL_B.append(A_UID_RIJ[i])
    geom_var['INIT'] = test_A
    geom_var['FINAL'] = test_B
    geom_var['IB'] = BL_A
    geom_var['FB'] = BL_B
    return(geom_var)

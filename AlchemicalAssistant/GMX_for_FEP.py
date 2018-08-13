import pandas as pd
import numpy as np
from AlchemicalAssistant.Vector_algebra import AtomNum2Symb, pairing_func

DUM_ANG_CONST = 10.0
DUM_BND_CONST = 300.0
DUM_THT_CONST = 0.3

def merge_two_dicts(x, y):
    """Given two dicts, merge them into a new dict as a shallow copy."""
    z = x.copy()
    z.update(y)
    return z

def Refine_PDB_file(fname):
    flines = open(fname, 'r+').readlines()
    pdb_lines = []
    for line in flines:
        if ('ATOM' in line) or ('HETATM' in line):
            line = line.rstrip()
            line = line.lstrip()
            ws = line.split()
            if not 'DUM' in line:
                pdb_lines.append([ws[2], ws[-3], ws[-2], ws[-1]])
    pdb_df = pd.DataFrame(pdb_lines)
    pdb_df.columns = ['AT', 'X', 'Y', 'Z']
    pdb_df[['X', 'Y', 'Z']] = pdb_df[['X', 'Y', 'Z']].apply(pd.to_numeric)
    return pdb_df


def printGRO(gro_name, resid='A2B'):
    import os
    os.system('cp plt.pdb %s.pdb' % gro_name)
    combo = Refine_PDB_file('plt.pdb')
    atoms = combo.AT
    coos = combo[['X', 'Y', 'Z']].as_matrix()
    pdb2gro(atoms, coos, gro_name, resid='A2B')
    return atoms


def pdb2gro(atoms, coos, gro_name, resid):
    grof = open(gro_name + '.gro', 'w+')
    grof.write('LIGPARGEN GENERATED GRO FILE\n')
    grof.write('%5d\n' % len(atoms))
    num = 0
    for (i, j) in zip(atoms, coos):
        num += 1
        grof.write('%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n' %
                   (1, resid, i, num, j[0] * 0.1, j[1] * 0.1, j[2] * 0.1))
    grof.write('%10.5f%10.5f%10.5f\n' % (1.00000, 1.00000, 1.00000))
    grof.close()
    return None


def gmxDihed(df, state='A'):
    [f1, f2, f3, f4] = [df['%sV1' % state], df['%sV2' %
                                               state], df['%sV3' % state], df['%sV4' % state]]
    cdat = [f2 + (f1 + f3) * 0.5, 1.5 * f3 - 0.5 * f1, 4.0 *
            f4 - f2, -2.0 * f3, -4.0 * f4, 0.0, 0.00]
    cdat = np.array(cdat) * 4.184
    return cdat


def ReturnInts(dfA, dfB):
    dfA.columns = ['A' + i for i in dfA.columns[:-1]] + [dfA.columns[-1]]
    dfB.columns = ['B' + i for i in dfB.columns[:-1]] + [dfB.columns[-1]]
    comT = list(set(dfA.UID) & set(dfB.UID))
    dfA = dfA.sort_values('UID')
    dfA.index = dfA.UID
    dfB = dfB.sort_values('UID')
    dfB.index = dfB.UID
    tAB = pd.concat([dfA[dfA.UID.isin(comT)], dfB[
                    dfB.UID.isin(comT)]], axis=1, join='inner')
    tA = dfA[~dfA.UID.isin(comT)]
    tB = dfB[~dfB.UID.isin(comT)]
    # tAB.to_csv('TAB.csv')
    # tA.to_csv('TA.csv')
    # tB.to_csv('TB.csv')
    return (tAB, tA, tB)


def PrintBNDS(bAB, bA, bB, itp):
    itp.write(
        '''
[ bonds ]
;          atomnrs   funct <<<<<<<<<<<<<<<<parametersA <<<<<<<<<<<<<<<<<<<<<<<<parametersB\n'''
    )
    for row in bAB.iterrows():
        i, r = row
        itp.write('%5d %5d %5d %11.4f %10.3f %11.4f %10.3f ; PERTURB\n' % (
            r.Acl1, r.Acl2, 1, r.ARIJ * 0.1, r.AKIJ * 836.8, r.BRIJ * 0.1, r.BKIJ * 836.8))
    for row in bA.iterrows():
        i, r = row
        itp.write('%5d %5d %5d %11.4f %10.3f %11.4f %10.3f ; ANNIHILATE\n' % (
            r.Acl1, r.Acl2, 1, r.ARIJ * 0.1, r.AKIJ * 836.8, 0.09, DUM_BND_CONST * 836.8))
    for row in bB.iterrows():
        i, r = row
        itp.write('%5d %5d %5d %11.4f %10.3f %11.4f %10.3f ; CREATE\n' % (
            r.Bcl1, r.Bcl2, 1, 0.09, DUM_BND_CONST * 836.8, r.BRIJ * 0.1, r.BKIJ * 836.8))
    return None
# ANGLES PRINTING


def PrintANGS(aAB, aA, aB, itp):
    itp.write(
        '''
[ angles ]
;          atomnrs   funct <<<<<<<<<<<<<<<<parametersA <<<<<<<<<<<<<<<<<<<<<<<<parametersB\n'''
    )
    for row in aAB.iterrows():
        i, r = row
        itp.write('%5d %5d %5d %5d %10.3f %10.3f %10.3f %10.3f ; PERTURB\n' % (
            r.Acl1, r.Acl2, r.Acl3, 1, r.AR, r.AK * 8.368, r.BR, r.BK * 8.368))
    for row in aA.iterrows():
        i, r = row
#        itp.write('%5d %5d %5d %5d %10.3f %10.3f %10.3f %10.3f ; ANNIHILATE\n'%(r.Acl1,r.Acl2,r.Acl3,1,r.AR,r.AK*8.368,r.AR,DUM_ANG_CONST*8.368))
        itp.write('%5d %5d %5d %5d %10.3f %10.3f %10.3f %10.3f ; ANNIHILATE\n' % (
            r.Acl1, r.Acl2, r.Acl3, 1, r.AR, r.AK * 8.368, r.AR, r.AK * 8.368))
    for row in aB.iterrows():
        i, r = row
#        itp.write('%5d %5d %5d %5d %10.3f %10.3f %10.3f %10.3f ; CREATE\n'%(r.Bcl1,r.Bcl2,r.Bcl3,1,r.BR,DUM_ANG_CONST*8.368,r.BR,r.BK*8.368))
        itp.write('%5d %5d %5d %5d %10.3f %10.3f %10.3f %10.3f ; CREATE\n' % (
            r.Bcl1, r.Bcl2, r.Bcl3, 1, r.BR, r.BK * 8.368, r.BR, r.BK * 8.368))
    return None


def PrintDIHS(tAB, tA, tB, itp):
    itp.write(
        '''
[ dihedrals ]
;                atomnrs   funct <<<<<<<<<<<<<<<<parametersA <<<<<<<<<<<<<<<<<<<<<<<<parametersB\n'''
    )
    for row in tAB.iterrows():
        i, r = row
        at = gmxDihed(r, 'A')
        bt = gmxDihed(r, 'B')
        itp.write('%5d%5d%5d%5d        3      %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f ; PERTURB\n' % (
            r.AI, r.AJ, r.AK, r.AL, at[0], at[1], at[2], at[3], at[4], 0.00, bt[0], bt[1], bt[2], bt[3], bt[4], 0.00))
    for row in tA.iterrows():
        i, r = row
        at = gmxDihed(r, 'A')
        bt = [0.0, 0.0, 0.0, 0.0, 0.0]
        itp.write('%5d%5d%5d%5d        3      %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f ; ANNIHILATE\n' % (
            r.AI, r.AJ, r.AK, r.AL, at[0], at[1], at[2], at[3], at[4], 0.00, bt[0], bt[1], bt[2], bt[3], bt[4], 0.00))
    for row in tB.iterrows():
        i, r = row
        at = [0.0, 0.0, 0.0, 0.0, 0.0]
        bt = gmxDihed(r, 'B')
        itp.write('%5d%5d%5d%5d        3      %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f ; CREATE\n' % (
            r.BI, r.BJ, r.BK, r.BL, at[0], at[1], at[2], at[3], at[4], 0.00, bt[0], bt[1], bt[2], bt[3], bt[4], 0.00))
    return None


def PrintPairs(tordf_A, tordf_B, itp):
    '''
    tordf_A: DF with proper and improper torsions of State A
    tordf_B: DF with proper and improper torsions of State B
    itp    : File object to write 
    '''
    pA = tordf_A[tordf_A.ADIH_TYP == False][['AI', 'AL']]
    pB = tordf_B[tordf_B.BDIH_TYP == False][['BI', 'BL']]
    pA['UID'] = [pairing_func(i, j) for i, j in zip(pA.AI, pA.AL)]
    pB['UID'] = [pairing_func(i, j) for i, j in zip(pB.BI, pB.BL)]
    a_dict = pA.set_index('UID')[['AI', 'AL']].to_dict()
    b_dict = pB.set_index('UID')[['BI', 'BL']].to_dict()
#    pALL_I = {**a_dict['AI'], **b_dict['BI']}
#    pALL_L = {**a_dict['AL'], **b_dict['BL']}
    pALL_I = merge_two_dicts(a_dict['AI'], b_dict['BI'])
    pALL_L = merge_two_dicts(a_dict['AL'], b_dict['BL'])
    itp.write('''[ pairs ]\n;   ai    aj   funct\n''')
    for k in list(pALL_I.keys()):
        itp.write('%5d %5d %5d\n' % (pALL_I[k], pALL_L[k], 1))
    return None


def PrintAtoms(amap, df_qlj, num2name, comAT, itp):
    qlj = df_qlj.set_index(
        'OPLSN')[['AN', 'AT', 'Q', 'EPS', 'SIG', 'AMASS']].to_dict()
    # print(qlj)
    itp.write('''; GENERATED BY: 
; Alchemical Assistant @ LigParGen Server
; Writted by Leela S. Dodda leela.dodda@yale.edu
; Jorgensen Lab @ Yale University

[ atomtypes ]\n''')
    for row in df_qlj.iterrows():
        i, r = row
        itp.write('%10s %5s %10.4f     0.000    A    %10.5E   %10.5E\n' % (
            'opls_%d' % r.OPLSN, r.AT, r.AMASS, r.SIG * 0.1, r.EPS * 4.184))
    itp.write('''
[ moleculetype ]
; Name               nrexcl
A2B                   3

[ atoms ]
;   nr       type  resnr    res   atom   cgnr     charge       mass      typeB    chargeB      massB comments\n''')
    for row in amap.iterrows():
        i, r = row
        a, b = r.ITY, r.FTY
#        print(a,b)
        if a not in list(qlj['AN'].keys()):  # a = b
            itp.write(' %5d %10s %6d %6s %6s %6d %10.4f %10.4f %10s %10.4f %10.4f\n' % (i + 1, 'opls_%d' % a, 1, 'A2B',
                                                                                        comAT[i], ((i + 1) / 30.0) + 1, 0.000, qlj['AMASS'][b], 'opls_%d' % b, qlj['Q'][b], qlj['AMASS'][b]))
        elif b not in list(qlj['AN'].keys()):  # a = b
            itp.write(' %5d %10s %6d %6s %6s %6d %10.4f %10.4f %10s %10.4f %10.4f\n' % (i + 1, 'opls_%d' % a, 1, 'A2B',
                                                                                        comAT[i], ((i + 1) / 30.0) + 1, qlj['Q'][a], qlj['AMASS'][a], 'opls_%d' % b, 0.00, qlj['AMASS'][a]))
        else:
            itp.write(' %5d %10s %6d %6s %6s %6d %10.4f %10.4f %10s %10.4f %10.4f\n' % (i + 1, 'opls_%d' % a, 1, 'A2B',
                                                                                        comAT[i], ((i + 1) / 30.0) + 1, qlj['Q'][a], qlj['AMASS'][a], 'opls_%d' % b, qlj['Q'][b], qlj['AMASS'][b]))
    return None


def WriteGmxTop(molA, molB, amap, df_qlj, rel_itp_name):
    tordf_A, tordf_B, angdf_A, angdf_B, bnddf_A, bnddf_B = molA['ALL_DIHEDS'], molB[
        'ALL_DIHEDS'], molA['ANGLES'], molB['ANGLES'], molA['BONDS'], molB['BONDS']
    df_qlj[['AN', 'Q', 'EPS', 'SIG']] = df_qlj[
        ['AN', 'Q', 'EPS', 'SIG']].apply(pd.to_numeric)
    bnddf_A[['cl1', 'cl2']] = bnddf_A[
        ['cl1', 'cl2']].apply(lambda x: x + 1, axis=0)
    tordf_A[['I', 'J', 'K', 'L']] = tordf_A[
        ['I', 'J', 'K', 'L']].apply(lambda x: x + 1, axis=0)
    tordf_B[['I', 'J', 'K', 'L']] = tordf_B[
        ['I', 'J', 'K', 'L']].apply(lambda x: x + 1, axis=0)
    # MANIPULATE ANGLES
    angdf_A[['cl1', 'cl2', 'cl3']] = angdf_A[
        ['cl1', 'cl2', 'cl3']].apply(lambda x: x + 1, axis=0)
    angdf_B[['cl1', 'cl2', 'cl3']] = angdf_B[
        ['cl1', 'cl2', 'cl3']].apply(lambda x: x + 1, axis=0)
    # MANIPULATE BONDS
    bnddf_B[['cl1', 'cl2']] = bnddf_B[
        ['cl1', 'cl2']].apply(lambda x: x + 1, axis=0)
    extra_rows = []
    for i in list(set(amap.ITY) - set(df_qlj.OPLSN)):
        extra_rows.append([i, 1, 'DA', 0.000000, 0.000000, 0.000000, 1.008])
    for i in list(set(amap.FTY) - set(df_qlj.OPLSN)):
        extra_rows.append([i, 1, 'DA', 0.000000, 0.000000, 0.000000, 1.008])
    df_qlj = df_qlj.append(pd.DataFrame(extra_rows, columns=[
                           'OPLSN', 'AN', 'AT', 'Q', 'EPS', 'SIG', 'AMASS']))
    combo_atoms = printGRO(rel_itp_name)
    tAB, tA, tB = ReturnInts(tordf_A, tordf_B)
    aAB, aA, aB = ReturnInts(angdf_A, angdf_B)
    bAB, bA, bB = ReturnInts(bnddf_A, bnddf_B)
    #print(molA['ZMAT'])
    za_df = molA['ZMAT']
#    za_df = pd.read_csv('ZM_A.csv')
    za_df = za_df[za_df.ITY > 0]
    zb_df = molB['ZMAT']
#    zb_df = pd.read_csv('ZM_B.csv')
    zb_df = zb_df[zb_df.ITY > 0]
    zb_df.ITY = zb_df.ITY + 9000
#    num2name = {**za_df.set_index('ITY')['ATOM'].to_dict(), **zb_df.set_index('ITY')['ATOM'].to_dict()}
    num2name = merge_two_dicts(za_df.set_index('ITY')['ATOM'].to_dict(), 
                               zb_df.set_index('ITY')['ATOM'].to_dict())
    itp = open('%s.itp' % rel_itp_name, 'w+')
    PrintAtoms(amap, df_qlj, num2name, combo_atoms, itp)
    PrintBNDS(bAB, bA, bB, itp)
    PrintANGS(aAB, aA, aB, itp)
    PrintDIHS(tAB, tA, tB, itp)
    PrintPairs(tordf_A, tordf_B, itp)
    itp.close()
    return None

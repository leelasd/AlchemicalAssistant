from MolReaders import ReadMol2File, make_graphs, ReadMolFile, SPZmat
from MolReaders import bossElement2Num, GetMol2Tripos, GetDummyType, EXPZmat, tor_id
from FEPBOSSReader import BOSSReader
from FEP_ZMAT import GetGeomVars, print_FEPZMAT
from Vector_algebra import pairing_func, ang_id, IsImproper, AtomNum2Mass
from GMX_for_FEP import WriteGmxTop
import pandas as pd
import sys
import pickle
from ReadZmat import AddImpLinestoZmat
from NAMD_Rel_FEP import MakeA2B

def merge_two_dicts(x, y):
    z = x.copy()   # start with x's keys and values
    z.update(y)    # modifies z with y's keys and values & returns None
    return z

def CallBOSSDualTopology(zname):
    '''
    Like Mentioned in Alchemify, This is all the info that it needs for now
    pickle files for both the molecules 
    ''' 
    amol = pickle.load(open("AMOL.p", "rb"))
    bmol = pickle.load(open("BMOL.p", "rb"))
    amol = amol.MolData
    bmol = bmol.MolData
    MakeA2B(amol,bmol,zname[:-2])
    return None


def WriteFEPZmat(ref_zmat, ali_zmat, IQV, FQV, QLJ_df, molA, molB):
    from rdkit import Chem
    from rdkit.Chem import rdMolAlign
    from rdkit.Chem import rdFMCS
    from rdkit.Chem.rdFMCS import FindMCS, AtomCompare, BondCompare
    import numpy as np
    import networkx as nx
    ref_mol = ref_zmat[:-2] + '.mol2'  # 'toluene.mol2'
    ali_mol = ali_zmat[:-2] + '.mol2'  # 'chlorobenzene.mol2'
    '''
    Do not sanitize the molecules, RDKit will freak out and give errors
    And All we want is to do MCSS, we dont care much about health of molecule
    '''
    mol1 = Chem.MolFromMol2File(ref_mol, removeHs=False, sanitize=False)
    mol2 = Chem.MolFromMol2File(ali_mol, removeHs=False, sanitize=False)
    tri_ref = GetMol2Tripos(ref_mol)
    tri_ali = GetMol2Tripos(ali_mol)
    # FIRST ONE
    _fmcs_params = dict(maximizeBonds=False, threshold=1.0, timeout=60,
                        verbose=False, matchValences=True,
                        ringMatchesRingOnly=True, completeRingsOnly=True,
                        atomCompare=AtomCompare.CompareAny,
                        bondCompare=BondCompare.CompareAny)
    try:
        mcs = rdFMCS.FindMCS([mol1, mol2], **_fmcs_params)
    except ValueError:
        print('\n Max Common Substructure calculation \n failed for this molecule!! \n Please be judicious ')
        sys.exit()
    core = Chem.MolFromSmarts(mcs.smartsString)
    match1 = mol1.GetSubstructMatch(core)
    match2 = mol2.GetSubstructMatch(core)
    from rdkit.Chem import AllChem
    AllChem.AlignMol(mol2, mol1, atomMap=list(zip(match2, match1)))
    #Chem.AlignMol(mol2, mol1, atomMap=list(zip(match2, match1)))
    Chem.MolToMolFile(mol2, 'aligned.mol', kekulize=False)
    cs_ref, ats_ref, bs_ref = ReadMol2File(ref_mol)
    cs_ali, ats_ali, bs_ali = ReadMolFile('aligned.mol')
    BtoA = {j: i for i, j in zip(match1, match2)}
    umatch1 = [i for i in range(mol1.GetNumAtoms()) if i not in match1]
    umatch2 = [i for i in range(mol2.GetNumAtoms()) if i not in match2]
    score = mol1.GetNumHeavyAtoms() + mol2.GetNumHeavyAtoms() - \
        2 * core.GetNumHeavyAtoms()
    score = 10**(-1.0 * score)
    print('SCORE FOR THIS TRANSFORMATION %3.3f' % (score))
    print(umatch1, umatch2)
    iR_1 = {at.GetIdx(): at.IsInRing() for at in mol1.GetAtoms()}
    iR_2 = {at.GetIdx(): at.IsInRing() for at in mol2.GetAtoms()}
    FEP_ITY = {i: int(IQV[i]) for i in list(IQV.keys())}
    GMX_FEP_ITY = {i: int(IQV[i]) for i in list(IQV.keys())}
    maxITY = np.array([int(i) for i in FEP_ITY.values()]).min()
    for ei, ej in zip(umatch2, range(mol1.GetNumAtoms(), mol1.GetNumAtoms() + len(umatch2))):
        BtoA[ei] = ej
        GMX_FEP_ITY[ej] = maxITY + ej
        if iR_2[ei]:
            FEP_ITY[ej] = maxITY + ej
        else:
            FEP_ITY[ej] = GetDummyType(tri_ali[ei])
    FEP_FTY = {}
    GMX_FEP_FTY = {}
    for i in list(FQV.keys()):
        FEP_FTY[BtoA[i]] = int(FQV[i])
        GMX_FEP_FTY[BtoA[i]] = int(FQV[i])
    maxFTY = np.array([int(i) for i in FEP_FTY.values()]).min()
    for i, j in zip(umatch1, range(mol2.GetNumAtoms(), mol2.GetNumAtoms() + len(umatch1))):
        GMX_FEP_FTY[i] = maxFTY + j
        if iR_1[i]:
            FEP_FTY[i] = maxFTY + j
        else:
            FEP_FTY[i] = GetDummyType(tri_ref[i])
##################
    ALL_AT = QLJ_df.set_index('OPLSN')['AT'].to_dict()
    ALL_AN = QLJ_df.set_index('OPLSN')['AN'].to_dict()
    ALL_QS = QLJ_df.set_index('OPLSN')['Q'].to_dict()
    ALL_SI = QLJ_df.set_index('OPLSN')['SIG'].to_dict()
    ALL_EP = QLJ_df.set_index('OPLSN')['EPS'].to_dict()
    FINAL_Q_SIG_EPS = {}
    assert(len(FEP_FTY.keys()) == len(FEP_ITY.keys()))
    for NN in range(len(FEP_FTY.keys())):
        if str(FEP_FTY[NN]) in list(QLJ_df.OPLSN):  # print('AS IS')
            va = str(FEP_FTY[NN])
            FINAL_Q_SIG_EPS[int(va)] = [ALL_AN[va], ALL_AT[va], ALL_QS[
                va], ALL_SI[va], ALL_EP[va]]
        elif str(FEP_ITY[NN]) in list(QLJ_df.OPLSN):  # print('MODIFIED')
            va = str(FEP_FTY[NN])
            vm = str(FEP_ITY[NN])
            FINAL_Q_SIG_EPS[int(va)] = [ALL_AN[vm], ALL_AT[
                vm], '0.000000', '0.000000', '0.000000']
    for NN in range(len(FEP_ITY.keys())):
        if str(FEP_ITY[NN]) in list(QLJ_df.OPLSN):  # print('AS IS')
            va = str(FEP_ITY[NN])
            FINAL_Q_SIG_EPS[int(va)] = [ALL_AN[va], ALL_AT[va], ALL_QS[
                va], ALL_SI[va], ALL_EP[va]]
        elif str(FEP_FTY[NN]) in list(QLJ_df.OPLSN):  # print('MODIFIED')
            vm = str(FEP_FTY[NN])
            va = str(FEP_ITY[NN])
            FINAL_Q_SIG_EPS[int(va)] = [ALL_AN[vm], ALL_AT[
                vm], '0.000000', '0.000000', '0.000000']
    exx_lines = []
    for LN in np.sort(list(FINAL_Q_SIG_EPS.keys())):
        if (LN != 100) and (LN != 4894):
            w = FINAL_Q_SIG_EPS[LN]
            exx_lines.append('{:>4d}{:>3s} {:<2s} {:>10s}{:>10s}{:>10s}\n'.format(
                LN, w[0], w[1], w[2], w[4], w[3]))
    f_qlj = pd.DataFrame(FINAL_Q_SIG_EPS).T
    f_qlj.columns = ['AN', 'AT', 'Q', 'EPS', 'SIG']
    f_qlj['OPLSN'] = f_qlj.index
    f_qlj = f_qlj[f_qlj.OPLSN > 100]
    f_qlj = f_qlj[['OPLSN', 'AN', 'AT', 'Q', 'EPS', 'SIG']]
    f_qlj['AMASS'] = [AtomNum2Mass(n) for n in f_qlj.AN]
    # GETTING INFO ABOUT DUMMIES
    typ_df = pd.DataFrame({
        'ITY': [GMX_FEP_ITY[i] for i in np.sort(list(GMX_FEP_ITY.keys()))],
        'FTY': [GMX_FEP_FTY[i] for i in np.sort(list(GMX_FEP_FTY.keys()))]})
    typ_df['NUM'] = [i for i in typ_df.index]
    typ_df['NUM_B'] = typ_df.FTY - 9800
    typ_df['NUM_A'] = typ_df.ITY - 800
#    f_qlj.to_csv('QLJ.csv')
    molA, molB = TranslateAliMol(molA, molB, BtoA)
    ##### THIS IS FOR GROMACS ####
    new_umatch2 = [BtoA[i] for i in umatch2]
    new_bs_ali = {'BI': [BtoA[i] for i in bs_ali['BI']], 'BJ': [
        BtoA[i] for i in bs_ali['BJ']], 'RIJ': [i for i in bs_ali['RIJ']]}
    new_ats_ali = {BtoA[i]: ats_ali[i] for i in ats_ali.keys()}
    new_cs_ali = {BtoA[i]: cs_ali[i] for i in cs_ali.keys()}
    G_A, MC_A = make_graphs(ats_ref, cs_ref, bs_ref)
    G_B, MC_B = make_graphs(new_ats_ali, new_cs_ali, new_bs_ali)
#    try: 
#        ALL_TORS = {**MC_A['TORSIONS'], **MC_B['TORSIONS']}
#        ALL_IMPS = {**MC_A['IMPROPERS'], **MC_B['IMPROPERS']}
#        ALL_ANGS = {**MC_A['ANGLES'], **MC_B['ANGLES']}
#    except SyntaxError: 
    ALL_TORS = merge_two_dicts(MC_A['TORSIONS'],  MC_B['TORSIONS'] )
    ALL_IMPS = merge_two_dicts(MC_A['IMPROPERS'], MC_B['IMPROPERS'])
    ALL_ANGS = merge_two_dicts(MC_A['ANGLES'],    MC_B['ANGLES']   )
# GET ALL BONDS
    ALL_BNDS = {}
    A_UID = []
    A_UID_RIJ = {}
    for i in MC_A['BONDS']:
        ALL_BNDS[pairing_func(i[0], i[1])] = i
        A_UID.append(pairing_func(i[0], i[1]))
        A_UID_RIJ[pairing_func(i[0], i[1])] = G_A[i[0]][i[1]]['distance']
    B_UID = []
    B_UID_RIJ = {}
    for i in MC_B['BONDS']:
        ALL_BNDS[pairing_func(i[0], i[1])] = i
        B_UID.append(pairing_func(i[0], i[1]))
        B_UID_RIJ[pairing_func(i[0], i[1])] = G_B[i[0]][i[1]]['distance']
    geom_var = GetGeomVars(ALL_BNDS, A_UID, B_UID,
                           A_UID_RIJ, B_UID_RIJ, FEP_ITY, FEP_FTY)
    #ALL_COOS = {**cs_ref, **new_cs_ali}
    #ALL_ATYS = {**ats_ref, **new_ats_ali}
    ALL_COOS = merge_two_dicts(cs_ref, new_cs_ali)
    ALL_ATYS = merge_two_dicts(ats_ref, new_ats_ali)
######
    G_combo = nx.DiGraph()
    for i in ALL_COOS.keys():
        G_combo.add_node(i, XYZ=ALL_COOS[i], elem=ALL_ATYS[i], init_no=FEP_ITY[i], finl_no=FEP_FTY[i],
                         atno=bossElement2Num(ALL_ATYS[i]))
    for k in ALL_BNDS.keys():
        [i, j] = ALL_BNDS[k]
        G_combo.add_edge(i, j)
        G_combo.add_edge(j, i)
    COMBO_MOL_ICOORDS = {'BONDS': ALL_BNDS, 'ANGLES': ALL_ANGS,
                         'TORSIONS': ALL_TORS, 'IMPROPERS': ALL_IMPS}
    rel_zmat_name = '%s_2_%s.z' % (ref_zmat[:-2], ali_zmat[:-2])
    oplsnum2aty = f_qlj.set_index(['OPLSN'])['AT'].to_dict()
    add_imp_lines = print_FEPZMAT(ALL_ATYS, G_combo, COMBO_MOL_ICOORDS, ALL_COOS, umatch1, new_umatch2,
                                  extra=exx_lines, G_vars=geom_var, zmat_name=rel_zmat_name, resid='A2B', num2typ=oplsnum2aty)
    EXPZmat(rel_zmat_name)
    AddImpLinestoZmat(rel_zmat_name, add_imp_lines)
    SPZmat(rel_zmat_name)
    rel_itp_name = rel_zmat_name[:-2]
    WriteGmxTop(molA, molB, typ_df, f_qlj, rel_itp_name)
    return rel_zmat_name


def TranslateAliMol(molA, molB, BtoA):
    molB['BONDS']['cl1'] = [BtoA[i] for i in molB['BONDS']['cl1']]
    molB['BONDS']['cl2'] = [BtoA[i] for i in molB['BONDS']['cl2']]
    molB['ANGLES']['cl1'] = [BtoA[i] for i in molB['ANGLES']['cl1']]
    molB['ANGLES']['cl2'] = [BtoA[i] for i in molB['ANGLES']['cl2']]
    molB['ANGLES']['cl3'] = [BtoA[i] for i in molB['ANGLES']['cl3']]
    molB['ALL_DIHEDS']['I'] = [BtoA[i] for i in molB['ALL_DIHEDS']['I']]
    molB['ALL_DIHEDS']['J'] = [BtoA[i] for i in molB['ALL_DIHEDS']['J']]
    molB['ALL_DIHEDS']['K'] = [BtoA[i] for i in molB['ALL_DIHEDS']['K']]
    molB['ALL_DIHEDS']['L'] = [BtoA[i] for i in molB['ALL_DIHEDS']['L']]
    molB['BONDS']['UID'] = [pairing_func(i, j) for i, j in zip(
        molB['BONDS']['cl1'], molB['BONDS']['cl2'])]
    molB['ANGLES']['UID'] = [ang_id([i, j, k]) for i, j, k in zip(
        molB['ANGLES']['cl1'], molB['ANGLES']['cl2'], molB['ANGLES']['cl3'])]
    molB['ALL_DIHEDS']['DIH_TYP'] = [
        IsImproper(t) for t in molB['ALL_DIHEDS']['TY']]
    molA['BONDS']['UID'] = [pairing_func(i, j) for i, j in zip(
        molA['BONDS']['cl1'], molA['BONDS']['cl2'])]
    molA['ANGLES']['UID'] = [ang_id([i, j, k]) for i, j, k in zip(
        molA['ANGLES']['cl1'], molA['ANGLES']['cl2'], molA['ANGLES']['cl3'])]
    molA['ALL_DIHEDS']['DIH_TYP'] = [
        IsImproper(t) for t in molA['ALL_DIHEDS']['TY']]
    bndlist = list(molB['BONDS']['UID'])+list(molA['BONDS']['UID'])
    molA['ALL_DIHEDS']['UID'] = [tor_id([i, j, k, l],SIS,bndlist) for i, j, k, l,SIS in
                                 zip(molA['ALL_DIHEDS']['I'], molA['ALL_DIHEDS']['J'],
                                     molA['ALL_DIHEDS']['K'], molA['ALL_DIHEDS']['L'],molA['ALL_DIHEDS']['DIH_TYP'])]
    molB['ALL_DIHEDS']['UID'] = [tor_id([i, j, k, l],SIS,bndlist) for i, j, k, l,SIS in 
                                 zip(molB['ALL_DIHEDS']['I'], molB['ALL_DIHEDS']['J'],
                                     molB['ALL_DIHEDS']['K'], molB['ALL_DIHEDS']['L'],molB['ALL_DIHEDS']['DIH_TYP'])]
    #molA['ALL_DIHEDS'].to_csv('DIH_A.csv', index=False)
    #molB['ALL_DIHEDS'].to_csv('DIH_B.csv', index=False)
    return (molA, molB)


def Alchemify(ref_z, ali_z):
    '''
    1. BOSSReader gets all the molecule data and writes it to pickle file
    '''
    import pandas as pd
    molA = BOSSReader(ref_z, optim=0, charge=0, lbcc=False)
    molA.cleanup()
    molB = BOSSReader(ali_z, optim=0, charge=0, lbcc=False)
    molB.cleanup()
    pickle.dump(molA, open("AMOL.p", "wb"))
    pickle.dump(molB, open("BMOL.p", "wb"))
    molA = molA.MolData
    molB = molB.MolData
    '''
    2.  Save The Z-matrix and Q-LJ data for future.  
    '''
    I_QLJ = molA['Q_LJ']
    F_QLJ = molB['Q_LJ']
    F_QLJ = {k: '9' + F_QLJ[k][1:] for k in F_QLJ.keys()}
    I_TY = {d: I_QLJ[d].split()[0] for d in I_QLJ.keys()}
    F_TY = {d: F_QLJ[d].split()[0] for d in F_QLJ.keys()}
    arrays_qlj = {i.strip().split()[0]: i.strip().split()
                  for i in I_QLJ.values()}
    for i in F_QLJ.values():
        word = i.strip().split()
        arrays_qlj[word[0]] = word
    df_QLJ = (pd.DataFrame(arrays_qlj).T)
    df_QLJ.columns = ['OPLSN', 'AN', 'AT', 'Q', 'EPS', 'SIG']
    '''
    3. WriteFEPZmat Does MCSS and Creates BOSS/MCPRO and GMX Stuff
    4. For CallBOSSDualTopology which also does NAMD needs just this info 
       and This Probably will change 
    '''
    print('Writing BOSS Single-Topology Z-matrix')
    r_zmat = WriteFEPZmat(ref_z, ali_z, I_TY, F_TY, df_QLJ, molA, molB)
    print('Writing BOSS Dual-Topology Z-matrix')
    CallBOSSDualTopology(r_zmat)
    return(None)

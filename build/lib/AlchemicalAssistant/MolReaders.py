import numpy as np
import networkx as nx
from AlchemicalAssistant.ReadZmat import GetZmatQLJDat
from AlchemicalAssistant.FEPBOSSReader import tor_cent
from AlchemicalAssistant.Vector_algebra import Distance, pairing_func
import os


def EXPZmat(zmat):
    '''
    Warning: 
    Need to have BOSSdir in the path
    zmat: string zmat file name with extension
    '''
    assert ('BOSSdir' in os.environ) and os.path.isfile((os.environ[
        'BOSSdir'] + '/scripts/xSPM')), 'Please Make sure $BOSSdir is defined \n xZmol and related files are in scripts directory of BOSS'
    execfile = os.environ['BOSSdir'] + '/scripts/xSPM > /tmp/olog'
    cmd = execfile + ' ' + zmat[:-2]
    os.system(cmd)
    os.system('cp sum %s.z' % (zmat[:-2]))
#    QLJ_DAT = GetZmatQLJDat(zmat)
    return None


def SPZmat(zmat):
    '''
    Warning: 
    Need to have BOSSdir in the path
    zmat: string zmat file name with extension
    '''
    assert ('BOSSdir' in os.environ) and os.path.isfile((os.environ[
        'BOSSdir'] + '/scripts/xSPM')), 'Please Make sure $BOSSdir is defined \n xZmol and related files are in scripts directory of BOSS'
    execfile = os.environ['BOSSdir'] + '/scripts/xSP > /tmp/olog'
    cmd = execfile + ' ' + zmat[:-2]
    os.system(cmd)
    os.system('cp sum %s.z' % (zmat[:-2]))
    return None


def ZMAT2MOL(zmat):
    '''
    Warning: 
    Need to have BOSSdir in the path
    zmat: string zmat file name with extension
    '''
    assert ('BOSSdir' in os.environ) and os.path.isfile((os.environ[
        'BOSSdir'] + '/scripts/xZmol')), 'Please Make sure $BOSSdir is defined \n xZmol and related files are in scripts directory of BOSS'
    execfile = os.environ['BOSSdir'] + '/scripts/xSPM > /tmp/olog'
    cmd = execfile + ' ' + zmat[:-2]
    os.system(cmd)
    cmd = 'mv plt.pdb %s.pdb' % zmat[:-2]
    os.system(cmd)
    os.system('babel -ipdb %s.pdb -omol2 %s.mol2 ---errorlevel 1 -b > LLN 2>&1' %
              (zmat[:-2], zmat[:-2]))
    QLJ_DAT = GetZmatQLJDat(zmat)
    return QLJ_DAT


def bossElement2Num(elem):
    symb2mass = {
        'H':  1,
        'B':  5,
        'C':  6,
        'N':  7,
        'O':  8,
        'F':  9,
        'Si': 14,
        'P': 15,
        'S': 16,
        'Cl': 17,
        'CL': 17,
        'Br': 35,
        'I': 53,
    }
    import re
    word = re.split('(\d+)', elem)
    elem = word[0]
    try:
        res = symb2mass[elem]
    except NameError:
        print("Mass for atom %s is not available \n add it to symb2mass dictionary")
    return res


def CheckExclusions(t, excl):
    combo = np.array([len(set(t) & set(e)) for e in excl])
    combo = combo[combo >= 2]
#    print(combo)
    if len(combo) > 0:
        return(True)
    else:
        return(False)


def ReadMol2File(fname):
    mollines = open(fname, 'r').readlines()
    [nats, nbonds] = map(int, mollines[2].split()[0:2])
    cooslines = mollines[7:7 + nats]
    coos = {}
    atypes = {}
    for i in range(nats):
        els = cooslines[i].split()
        coos[i] = [float(e) for e in els[2:5]]
        atypes[i] = els[1]
    bondlines = mollines[8 + nats:8 + nats + nbonds]
    bonds = {'BI': [], 'BJ': [], 'RIJ': []}
    for line in bondlines:
        [bi, bj] = map(int, line.split()[1:3])
        bonds['BI'].append(bi - 1)
        bonds['BJ'].append(bj - 1)
        bonds['RIJ'].append(Distance(coos[bi - 1], coos[bj - 1]))
    return (coos, atypes, bonds)


def ReadMolFile(fname):
    mollines = open(fname, 'r').readlines()
    [nats, nbonds] = map(int, mollines[3].split()[0:2])
    cooslines = mollines[4:4 + nats]
    coos = {}
    atypes = {}
    for i in range(nats):
        els = cooslines[i].split()
        coos[i] = [float(e) for e in els[0:3]]
        atypes[i] = els[3]
    bondlines = mollines[4 + nats:4 + nats + nbonds]
    bonds = {'BI': [], 'BJ': [], 'RIJ': []}
    for line in bondlines:
        [bi, bj] = map(int, line.split()[0:2])
        bonds['BI'].append(bi - 1)
        bonds['BJ'].append(bj - 1)
        bonds['RIJ'].append(Distance(coos[bi - 1], coos[bj - 1]))
    return (coos, atypes, bonds)


def tor_id(a,isImp,bndlist):
    bond = pairing_func(a[1], a[2])
    ends = pairing_func(a[0], a[3])
    if isImp: 
        ndata = tor_cent(a, bndlist)
        aimp = np.sort(ndata[1:])
        return(ang_id(aimp))
    else: return(pairing_func(bond, ends))


def ang_id(a):
    bond_a = pairing_func(a[0], a[1])
    bond_b = pairing_func(a[1], a[2])
    return(pairing_func(bond_a, bond_b))


def impr_id(a):
    a = np.sort(np.array(a))
    bond_a = pairing_func(a[0], a[1])
    bond_b = pairing_func(a[1], a[2])
    return(pairing_func(bond_a, bond_b))


def make_graphs(atoms, coos, bonds):
    G = nx.DiGraph()
    # ADD NODES USING ATOM TYPES AND COORDINATES
    for i in coos.keys():
        G.add_node(i, XYZ=coos[i], elem=atoms[i],
                   atno=bossElement2Num(atoms[i]))
    for (i, j, rij) in zip(bonds['BI'], bonds['BJ'], bonds['RIJ']):
        G.add_edge(i, j, distance=rij)
        G.add_edge(j, i, distance=rij)
#    all_ps = dict(nx.algorithms.all_pairs_shortest_path_length(G))
    all_ps = dict(nx.algorithms.all_pairs_shortest_path_length(G,cutoff=3))
    all_paths = []
    for s in all_ps.keys():
        for e in all_ps[s].keys():
            if all_ps[s][e] == 1:
                all_paths += list(nx.algorithms.all_simple_paths(G,s,e,cutoff=1)) 
#                all_paths += list(nx.algorithms.shortest_simple_paths(G, s, e))
            elif all_ps[s][e] == 2:
                all_paths += list(nx.algorithms.all_simple_paths(G,s,e,cutoff=2)) 
#                all_paths += list(nx.algorithms.shortest_simple_paths(G, s, e))
            elif all_ps[s][e] == 3:
                all_paths += list(nx.algorithms.all_simple_paths(G,s,e,cutoff=3)) 
#                all_paths += list(nx.algorithms.shortest_simple_paths(G, s, e))
    all_bonds = [p for p in all_paths if len(set(p)) == 2]
    new_angs = [p for p in all_paths if len(set(p)) == 3]
    new_tors = [p for p in all_paths if len(set(p)) == 4]
    bndlist = [pairing_func(a[0],a[1]) for a in all_bonds]
    dict_new_tors = {tor_id(t,False,bndlist): t for t in new_tors}
    dict_new_angs = {ang_id(t): t for t in new_angs}
    imp_keys = [n for n in G.nodes() if G.degree(n) / 2 == 3]
    all_imps = {}
    for i in imp_keys:
        nei = np.sort(list(G.neighbors(i)))
        if G.node[i]['atno'] == 6:
            all_imps[impr_id([nei[0], nei[1], nei[2]])] = [
                nei[0], i, nei[1], nei[2]]
    # print(all_imps)
    MOL_ICOORDS = {'BONDS': all_bonds,
                   'ANGLES': dict_new_angs, 'TORSIONS': dict_new_tors, 'IMPROPERS': all_imps}
    return(G, MOL_ICOORDS)


def GetMol2Tripos(mol):
    lines = open(mol, 'r').readlines()
    nats = int(lines[2].split()[0])
    coos = [lines[7 + aa].split()[5] for aa in range(0, nats)]
    return coos


def GetDummyType(typ):
    typ2dum = {'C.3': 100, 'H': 100, 'C.cat': 100, 'N.3': 100, 'C.2': 4894, 'C.1': 4894, 'C.ar': 4894, 'N.2': 4894, 'N.1': 4894, 'N.ar': 4894, 'N.am': 4894,
               'N.pl3': 4894, 'O.2': 4894, 'S.2': 4894, 'S.O': 4894, 'S.O2': 4894, 'N.4': 100, 'O.3': 100, 'S.3': 100, 'P.3': 100, 'F': 100, 'Cl': 100, 'Br': 100}
    #print(typ, typ2dum[typ])
    if typ in typ2dum.keys():
        return typ2dum[typ]
    else:
        return 100

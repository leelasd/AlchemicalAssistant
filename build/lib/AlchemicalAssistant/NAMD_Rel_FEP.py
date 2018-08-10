import pickle 
import numpy as np 
import pandas as pd
from AlchemicalAssistant.FEPBOSSReader import bossPdbAtom2Element,ucomb,tor_cent
from AlchemicalAssistant.Vector_algebra import pairing_func,AtomNum2Symb,AtomNum2Mass
from AlchemicalAssistant.MolReaders import ang_id,tor_id
from AlchemicalAssistant.TINKER_Rel_FEP import xyz_prep,tinker_prm

def pdb_prep(atoms, coos, resid='A2B',pdbname='COMBO'):
    opdb = open(pdbname+'_NAMD.pdb', 'w+')
    opdb.write('REMARK LIGPARGEN GENERATED PDB FILE\n')
    num = 0
    for (i, j) in zip(atoms, coos):
        num += 1
        opdb.write('%-6s%5d %4s %3s  %4d    %8.3f%8.3f%8.3f\n' %
                   ('ATOM', num, i, resid, 1, j[0], j[1], j[2]))
    opdb.write('END\n')
    opdb.close()
    return None

###
def MapMolecules(map_dict,num): 
    if num in map_dict.keys(): return map_dict[num] 
    else: return num+1

def TranslateICs(amol,bmol,map_dict,zdf):
    amol['BONDS'][['cl1','cl2']] = amol['BONDS'][['cl1','cl2']].apply(lambda x: x+3)
    bmol['BONDS'][['cl1','cl2']] = bmol['BONDS'][['cl1','cl2']].applymap(lambda x: MapMolecules(map_dict,x+3))
    all_bonds = pd.concat([amol['BONDS'],bmol['BONDS']],axis=0) 
    all_bonds['UID'] = [pairing_func(i,j) for i,j in zip(all_bonds.cl1,all_bonds.cl2)] 
    amol['ANGLES'][['cl1','cl2','cl3']] = amol['ANGLES'][['cl1','cl2','cl3']].apply(lambda x: x+3)
    bmol['ANGLES'][['cl1','cl2','cl3']] = bmol['ANGLES'][['cl1','cl2','cl3']].applymap(lambda x: MapMolecules(map_dict,x+3))
    all_angles = pd.concat([amol['ANGLES'],bmol['ANGLES']],axis=0) 
    all_angles['UID'] = [ang_id([i,j,k]) for i,j,k in zip(all_angles.cl1,all_angles.cl2,all_angles.cl3)]
    amol['ALL_DIHEDS'][['I', 'J', 'K', 'L']] = amol['ALL_DIHEDS'][['I', 'J', 'K', 'L']].apply(lambda x: x+3)
    bmol['ALL_DIHEDS'][['I', 'J', 'K', 'L']] = bmol['ALL_DIHEDS'][['I', 'J', 'K', 'L']].applymap(lambda x: MapMolecules(map_dict,x+3))
    amol['ALL_DIHEDS']['STATE'] = ['INI' for i in amol['ALL_DIHEDS']['I']]
    bmol['ALL_DIHEDS']['STATE'] = ['FIN' for i in bmol['ALL_DIHEDS']['I']]
    all_diheds = pd.concat([amol['ALL_DIHEDS'],bmol['ALL_DIHEDS']],axis=0) 
    imp_tys = [160,161,162,221]
    is_imp = []
    for i in all_diheds.TY: 
        if i in imp_tys: is_imp.append(True)
        else: is_imp.append(False)
    all_diheds['IMP'] = is_imp 
    all_diheds['UID'] = [tor_id([r.I,r.J,r.K,r.L],r.IMP,list(all_bonds.UID)) for i,r in all_diheds.iterrows()] 

    bnd_dat = zdf[zdf.NJ>2][['NUM','NJ']]
    bnd_dat['UID'] = [pairing_func(r.NUM,r.NJ) for i,r in bnd_dat.iterrows()]
    ang_dat = zdf[zdf.NK>2][['NUM','NJ','NK']]
    ang_dat['UID'] = [ang_id([r.NUM,r.NJ,r.NK]) for i,r in ang_dat.iterrows()]
    tor_dat = zdf[zdf.NL>2][['NUM','NJ','NK','NL']]
    tor_dat['IMP'] = [False if ucomb([r.NUM,r.NJ,r.NK,r.NL], list(all_bonds.UID)) == 3 else True for n,r in tor_dat.iterrows()]
    tor_dat['UID'] = [tor_id([r.NUM,r.NJ,r.NK,r.NL],r.IMP,list(all_bonds.UID)) for i,r in tor_dat.iterrows()]
    add_diheds = all_diheds[~all_diheds.UID.isin(tor_dat.UID)]
    add_angles = all_angles[~all_angles.UID.isin(ang_dat.UID)] 
    add_bonds = all_bonds[~all_bonds.UID.isin(bnd_dat.UID)] 
    I_QLJ = amol['Q_LJ']
    F_QLJ = bmol['Q_LJ']
    F_QLJ = {k: '9' + F_QLJ[k][1:] for k in F_QLJ.keys()}
    #print(I_QLJ,F_QLJ)
    all_qljs = list(I_QLJ.values())+ list(F_QLJ.values())
    #print(add_diheds)
    comboIC = {
            'ADD_BONDS': add_bonds,
            'ADD_ANGLE': add_angles,
            'ADD_DIHED': add_diheds,
            'BONDS': all_bonds,
            'ANGLE': all_angles,
            'DIHED': all_diheds,
            'QQLJS': all_qljs, 
            }
    return comboIC

def MakeA2B(amol,bmol,zname,shrink=False):
    com_coos = pd.concat([ amol['XYZ'],bmol['XYZ']],axis=0)
    com_coos.index = [i for i in range(3,len(com_coos.X)+3)]
    pertz_bmol = bmol['ZMAT'][2:]
    azmat = amol['ZMAT']
    counter = azmat.NUM.count()
    final_index = [ i for i in range(1,counter+pertz_bmol.NUM.count()+1)]
    a2b_zdf = pd.concat([azmat,pertz_bmol],axis=0)
    a2b_zdf.index = final_index
    ### GET DATA BEFORE YOU CHANGE THE BOND LENGTHS 
    common = (a2b_zdf[(a2b_zdf.index>3)]['RIJ'].to_dict())
    ### CHANGE THE ATOM NUMBERS AND OTHER STUFF NOW 
    map_dict = {}
    for i, r in a2b_zdf[counter:].iterrows():  map_dict[r.NUM] = i
    pertz_bmol.NUM = [MapMolecules(map_dict,num) for num in pertz_bmol.NUM]
    pertz_bmol.NJ  = [MapMolecules(map_dict,num) for num in pertz_bmol.NJ]
    pertz_bmol.NK  = [MapMolecules(map_dict,num) for num in pertz_bmol.NK]
    pertz_bmol.NL  = [MapMolecules(map_dict,num) for num in pertz_bmol.NL]
# DO NOT DELETE THESE LINES THESE ARE IMPT FOR SHRINKING THE MOLECULE
    if shrink: 
        pert_init_bond = [] 
        for i,j in zip(pertz_bmol.NUM,pertz_bmol.RIJ): #print(i)
            if i > counter+1: pert_init_bond.append(0.3)
            else: pert_init_bond.append(j)
        pertz_bmol.RIJ = [0.05]+ pert_init_bond[1:]
        pertz_bmol.ITY = [4894 for num in pertz_bmol.ITY]
        azmat.FTY = [4894 if num>0 else -1 for num in azmat.ITY]
    else: 
        pertz_bmol.RIJ = [0.01]+ list(pertz_bmol.RIJ)[1:]#pert_init_bond[1:]
        pertz_bmol.FTY = [num+9000 for num in pertz_bmol.FTY]
        pertz_bmol.ITY = [num-100 for num in pertz_bmol.FTY]
        azmat.FTY = [num+8800 if num>0 else -1 for num in azmat.ITY]
    a2b_zdf = pd.concat([azmat,pertz_bmol],axis=0)
    a2b_zdf.index = list(a2b_zdf.NUM) 
    a2b_zdf[['RIJ','TIJK','PIJKL']] = a2b_zdf[['RIJ','TIJK','PIJKL']].apply(pd.to_numeric)
    fin_geomvars = {}
    for i,r in a2b_zdf[a2b_zdf.index>3].iterrows(): #print(i,r.RIJ)  
        if r.RIJ == 0.05: fin_geomvars[i]=0.05
        elif r.RIJ!=0.3: fin_geomvars[i]=0.3
        else: fin_geomvars[i] = float(common[i])
    a2b_geomvars = a2b_zdf[a2b_zdf.index>3][['NUM','RIJ']] 
    a2b_geomvars['FIN_RIJ'] = [fin_geomvars[i] for i in a2b_geomvars.NUM]
    if shrink: 
        a2b_geomvars = a2b_geomvars[a2b_geomvars.FIN_RIJ>0.05]
    else: a2b_geomvars=None
    a2b_zdf.ATOM = [bossPdbAtom2Element(attype) for attype in a2b_zdf.ATOM]
    com_coos.at_symb = list(a2b_zdf.ATOM[2:])
    comboIC = TranslateICs(amol,bmol,map_dict,a2b_zdf)
    WriteNamdRtfAndPrm(amol,bmol,comboIC,zname)
    WriteDualTopZmat(a2b_zdf,a2b_geomvars,comboIC,counter,zname)
    return None
    #return(a2b_zdf,a2b_geomvars,com_coos,comboIC)
######### MAIN PART

def WriteNamdRtfAndPrm(amol,bmol,comboIC,zname):
    all_xyz = pd.concat([amol['XYZ'],bmol['XYZ']],axis=0)
    all_xyz['SYM']=[AtomNum2Symb(i) for i in all_xyz.at_num]
    all_xyz.index = range(0,len(all_xyz.X))
    all_xyz['UID']=['%s%02d'%(i,j) for i,j in zip(all_xyz.SYM,all_xyz.index) ]
    num2symb = {i+3:j for i,j in zip(all_xyz.index,all_xyz.UID)}
    comboIC['BONDS']['Ncl1'] = [num2symb[i] for i in comboIC['BONDS']['cl1']]
    comboIC['BONDS']['Ncl2'] = [num2symb[i] for i in comboIC['BONDS']['cl2']]
    comboIC['ANGLE']['Ncl1'] = [num2symb[i] for i in comboIC['ANGLE']['cl1']]
    comboIC['ANGLE']['Ncl2'] = [num2symb[i] for i in comboIC['ANGLE']['cl2']]
    comboIC['ANGLE']['Ncl3'] = [num2symb[i] for i in comboIC['ANGLE']['cl3']]
    dataQLJ = [i.split() for i in comboIC['QQLJS']]
    dfQLJ = pd.DataFrame(dataQLJ)
    dfQLJ.columns = ['OPLSN','atom_num','at_ty','Q','SIG','EPS']
    final_df = (pd.concat([all_xyz,dfQLJ],axis=1))
    final_df[['OPLSN','at_num','Q','SIG','EPS']]=final_df[['OPLSN','at_num','Q','SIG','EPS']].apply(pd.to_numeric)
    final_df['NAMD']=['%s%3d'%(s,num-8900) if num>9000 else '%s%3d'%(s,num) for s,num in zip(final_df.SYM,final_df.OPLSN)]
    num2namd = {i+3:j for i,j in zip(final_df.index,final_df.NAMD)}
    #print(final_df)
    rtf = open(zname+'.rtf','w+')
    rtf.write('! Relative FEP RTF file for NAMD/CHARMM by LSD\n')
    for i,r in final_df.iterrows(): rtf.write('MASS %d %s %10.6f %s\n'%(i+1,r.NAMD,AtomNum2Mass(r.at_num),r.SYM))
    rtf.write('AUTO ANGLE DIHE\n')
    rtf.write('RESI A2B %6.3f\n'%(final_df.Q.sum()))
    for i,r in final_df.iterrows(): rtf.write('ATOM %5s %5s %12.6f\n'%(r.UID,r.NAMD,r.Q))
    for i,r in comboIC['BONDS'].iterrows(): rtf.write('BOND %5s %5s\n'%(r.Ncl1,r.Ncl2))
#    imp_names={} 
    for i,r in comboIC['DIHED'][comboIC['DIHED'].IMP==True].iterrows(): #rtf.write('IMPR %5s %5s %5s %5s %5s'%())
        imp_order = tor_cent([r.I,r.J,r.K,r.L],list(comboIC['BONDS'].UID))
#        imp_names[i] = [num2namd[i] for i in imp_order]
        ans = [num2symb[i] for i in imp_order]
        rtf.write('IMPR %5s %5s %5s %5s\n'%(ans[1],ans[0],ans[2],ans[3]))
    rtf.write('PATCH FIRST NONE LAST NONE\nEND\n')
    rtf.close()
    coos = all_xyz[['X','Y','Z']].as_matrix()
    pdb_prep(list(all_xyz.UID), coos,pdbname=zname)
    xyz_prep(list(all_xyz.UID), coos,comboIC['BONDS'],final_df, pdbname=zname)
    tinker_prm(final_df,comboIC, pdbname=zname)
#    pickle.dump(comboIC, open(zname + ".p", "wb"))
#    tinker_prm(comboIC, resid=zname)
    ### RTF FILE AND PDB FILE WORKS PERFECTLY FINE ## NOW JUST THE PARAMETER FILE :D 
#    print(final_df)
    prm = open(zname+'.prm','w+')
    prm.write('BOND\n')
    for i,r in comboIC['BONDS'].iterrows():prm.write('%-5s %5s %12.6f %12.6f\n'%(num2namd[r.cl1],num2namd[r.cl2],r.KIJ,r.RIJ))
    prm.write('\nANGLE\n')
    for i,r in comboIC['ANGLE'].iterrows():prm.write('%-5s %5s %5s %12.6f %12.6f\n'%(num2namd[r.cl1],num2namd[r.cl2],num2namd[r.cl3],r.K,r.R))
    #print( comboIC['DIHED'])
    prm.write('\nDIHEDRAL \n')
    for i,r in comboIC['DIHED'][comboIC['DIHED'].IMP==False].iterrows(): 
       torname = '-'.join([num2namd[i] for i in [r.I,r.J,r.K,r.L]])
       for pot in range(1, 5):
           prm.write('%s %4.5f %d %4.5f \n' % (torname.replace("-", " "), r['V' + str(pot)], pot, 180.00 * abs(pot % 2 - 1)))
    prm.write('\nIMPROPER \n')
    for i,r in comboIC['DIHED'][comboIC['DIHED'].IMP==True].iterrows(): 
       pot = 2
       imp_order = tor_cent([r.I,r.J,r.K,r.L],list(comboIC['BONDS'].UID))
       imprname = [num2namd[i] for i in imp_order]
#       imprname = imp_names[i] 
       nameImpr = '-'.join([imprname[1],imprname[0],imprname[2],imprname[3]])
       prm.write('%s %4.5f %d %4.5f \n' % (nameImpr.replace(
           "-", " "), r['V' + str(pot)], pot, 180.00 * abs(pot % 2 - 1)))
    prm.write(
        '\nNONBONDED nbxmod 5 atom cdiel switch vatom vdistance vswitch - \ncutnb 14.0 ctofnb 12.0 ctonnb 11.5 eps 1.0 e14fac 0.5  geom\n')
    for i,r in final_df.iterrows():
        prm.write('%s 0.00 %3.6f %3.6f 0.00 %3.6f %3.6f \n'%(r.NAMD,r.EPS*-1.0,r.SIG*0.561231,r.EPS*-0.5,r.SIG*0.561231))
    prm.close()
    return None

############# FOR BOSS/MCPRO Dual-top FEPS ######################
def WriteDualTopZmat(zdf,a2b_geomvars,comboIC,counter,zname):
    ofile = open(zname+'_DT.z', 'w+')
    ofile.write('BOSS/MCPRO DualTopology Z-Matrix by LSD\n')
    for i,r in zdf.iterrows():
        ofile.write('%4d %-3s%5d%5d%5d%12.6f%4d%12.6f%4d%12.6f%4s%5d\n'%(r.NUM,r.ATOM,r.ITY,r.FTY,r.NJ,r.RIJ,r.NK,r.TIJK,r.NL,r.PIJKL,'A2B',1))
    ofile.write(
        '                    Geometry Variations follow    (2I4,F12.6)\n')
    if a2b_geomvars is not None: 
       for k,r in a2b_geomvars.iterrows(): ofile.write('%4d%4d%12.6f\n' % (r.NUM, 1,r.FIN_RIJ))
    ofile.write('                     Variable Bonds follow         (I4)\n')
    for k,r in zdf[zdf.NJ>2].iterrows(): 
        if pairing_func(r.NUM,r.NJ) in list(comboIC['BONDS'].UID): ofile.write('%4d\n' % (r.NUM))
    #for k,r in a2b_geomvars.iterrows(): ofile.write('%4d\n' % (r.NUM))
    ofile.write('                    Additional Bonds follow       (2I4)\n')
    bonds = comboIC['ADD_BONDS']
    for k,r in bonds.iterrows(): ofile.write('%4d%4d\n'%(r.cl1,r.cl2))
    ofile.write('''                    Harmonic Constraints follow   (2I4,4F10.4)
                    Variable Bond Angles follow   (I4)\n''')
    for k,r in zdf[zdf.NK>2].iterrows(): 
        if ang_id([r.NUM,r.NJ,r.NK]) in list(comboIC['ANGLE'].UID): ofile.write('%4d\n' % (r.NUM))
    angles = comboIC['ADD_ANGLE']
    ofile.write('                    Additional Bond Angles follow (3I4)\n')
    for k,r in angles.iterrows(): ofile.write('%4d%4d%4d\n'%(r.cl1,r.cl2,r.cl3)) 
    ofile.write(
        '                    Variable Dihedrals follow     (3I4,F12.6)\n')
    for k,r in zdf[zdf.NL>2].iterrows(): 
        arr = np.array([r.NUM,r.NJ,r.NK,r.NL])
        arr_A = arr[arr<=counter]
        arr_B = arr[arr>counter]
        if (len(arr_A)==4 or len(arr_B)==4): #ofile.write('%4d%4d%4d%12.6f\n' % (r.NUM,-1, -1, 0.000)) 
            if len(arr_A)==4: ofile.write('%4d%4d%4d%12.6f\n' % (r.NUM,-1, 100, 0.000))
            elif len(arr_B)==4: ofile.write('%4d%4d%4d%12.6f\n' % (r.NUM,100, -1, 0.000))
    ofile.write('                    Additional Dihedrals follow   (6I4)\n')
    torsions = comboIC['ADD_DIHED']
    for k,r in torsions[torsions.STATE == 'INI'].iterrows(): ofile.write('%4d%4d%4d%4d%4d%4d\n' %(r.I,r.J,r.K,r.L,r.TY,100))
    for k,r in torsions[torsions.STATE == 'FIN'].iterrows(): ofile.write('%4d%4d%4d%4d%4d%4d\n' %(r.I,r.J,r.K,r.L,100,r.TY))
    ofile.write('                    Domain Definitions follow     (4I4)\n')
    ofile.write('%4d%4d%4d%4d\n' % (3, counter,counter+1,zdf.NUM.count()))
    ofile.write(
        '''                    Conformational Search (2I4,2F12.6)
                    Local Heating Residues follow (I4 or I4-I4)
                    Final blank line\n
 Final Non-Bonded Parameters for QM (AM1 CM1Ax1.14) Atoms:\n
''')
    all_qljlines = {int(line.lstrip()[0:4]): line for line in comboIC['QQLJS']}
    #for l in comboIC['QQLJS']: ofile.write('%s'%l) 
    if a2b_geomvars is None: 
        new_lines = []
        for line in comboIC['QQLJS']:  
            if int(line.lstrip()[0:4]) < 9000: all_qljlines[int(line[0:4])+8800] = ('%4d'%(int(line[0:4])+8800)+line[4:13]+'0.000000  0.000000  0.000000\n') 
            elif int(line.lstrip()[0:4])>=9800: all_qljlines[int(line[0:4])-100] = ('%4d'%(int(line[0:4])-100 )+line[4:13]+'0.000000  0.000000  0.000000\n')
    #print(all_qljlines)
    sorted_keys = np.sort(list(all_qljlines.keys()))
    for k in sorted_keys: ofile.write('%s'%all_qljlines[k])
    ofile.close()
    return None

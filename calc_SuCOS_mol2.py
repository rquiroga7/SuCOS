# Copyright <2019> <University of Oxford>
# This code is licensed under MIT license (see LICENSE.txt for details)

import argparse, os, gzip
import numpy as np
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem, rdShapeHelpers
from rdkit.Chem.FeatMaps import FeatMaps
from rdkit import RDConfig
from rdkit.Chem import rdMolAlign
from rdkit import RDLogger
from rdkit.Chem.MolStandardize import rdMolStandardize





#################################################
#### Setting up the features to use in FeatureMap
fdef = AllChem.BuildFeatureFactory(os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef'))
#    keep = ('Donor','Acceptor','NegIonizable','PosIonizable','Aromatic')

fmParams = {}
for k in fdef.GetFeatureFamilies():
    fparams = FeatMaps.FeatMapParams()
    fmParams[k] = fparams

keep = ('Donor', 'Acceptor', 'NegIonizable', 'PosIonizable', 'ZnBinder',
        'Aromatic', 'Hydrophobe', 'LumpedHydrophobe')
#################################################
def _FragIndicesToMol(oMol, indices):
    em = Chem.rdchem.EditableMol(Chem.rdchem.Mol())

    newIndices = {}
    for i, idx in enumerate(indices):
        em.AddAtom(oMol.GetAtomWithIdx(idx))
        newIndices[idx] = i

    for i, idx in enumerate(indices):
        at = oMol.GetAtomWithIdx(idx)
        for bond in at.GetBonds():
            if bond.GetBeginAtomIdx() == idx:
                oidx = bond.GetEndAtomIdx()
            else:
                oidx = bond.GetBeginAtomIdx()
            # make sure every bond only gets added once:
            if oidx < idx:
                continue
            em.AddBond(newIndices[idx], newIndices[oidx], bond.GetBondType())
    res = em.GetMol()
    res.ClearComputedProps()
    Chem.rdmolops.GetSymmSSSR(res)
    res.UpdatePropertyCache(False)
    res._idxMap = newIndices
    return res

"""
sanifix4.py
Original code from rdkit [James Davidson]
"""
import logging

from rdkit import Chem


def _FragIndicesToMol(oMol, indices):
    em = Chem.EditableMol(Chem.Mol())

    newIndices = {}
    for i, idx in enumerate(indices):
        em.AddAtom(oMol.GetAtomWithIdx(idx))
        newIndices[idx] = i

    for i, idx in enumerate(indices):
        at = oMol.GetAtomWithIdx(idx)
        for bond in at.GetBonds():
            if bond.GetBeginAtomIdx() == idx:
                oidx = bond.GetEndAtomIdx()
            else:
                oidx = bond.GetBeginAtomIdx()
            # make sure every bond only gets added once:
            if oidx < idx:
                continue
            em.AddBond(newIndices[idx], newIndices[oidx], bond.GetBondType())
    res = em.GetMol()
    res.ClearComputedProps()
    Chem.GetSymmSSSR(res)
    res.UpdatePropertyCache(False)
    res._idxMap = newIndices
    return res


def _recursivelyModifyNs(mol, matches, indices=None):
    if indices is None:
        indices = []
    res = None
    while len(matches) and res is None:
        tIndices = indices[:]
        nextIdx = matches.pop(0)
        tIndices.append(nextIdx)
        nm = Chem.Mol(mol.ToBinary())
        nm.GetAtomWithIdx(nextIdx).SetNoImplicit(True)
        nm.GetAtomWithIdx(nextIdx).SetNumExplicitHs(1)
        cp = Chem.Mol(nm.ToBinary())
        try:
            Chem.SanitizeMol(cp)
        except ValueError:
            res, indices = _recursivelyModifyNs(nm, matches, indices=tIndices)
        else:
            indices = tIndices
            res = cp
    return res, indices


def AdjustAromaticNs(m, nitrogenPattern="[n&D2&H0;r5,r6]"):
    """
    default nitrogen pattern matches Ns in 5 rings and 6 rings in order to be able
    to fix: O=c1ccncc1
    """
    Chem.GetSymmSSSR(m)
    m.UpdatePropertyCache(False)

    # break non-ring bonds linking rings:
    em = Chem.EditableMol(m)
    linkers = m.GetSubstructMatches(Chem.MolFromSmarts("[r]!@[r]"))
    plsFix = set()
    for a, b in linkers:
        em.RemoveBond(a, b)
        plsFix.add(a)
        plsFix.add(b)
    nm = em.GetMol()
    for at in plsFix:
        at = nm.GetAtomWithIdx(at)
        if at.GetIsAromatic() and at.GetAtomicNum() == 7:
            at.SetNumExplicitHs(1)
            at.SetNoImplicit(True)

    # build molecules from the fragments:
    fragLists = Chem.GetMolFrags(nm)
    frags = [_FragIndicesToMol(nm, x) for x in fragLists]

    # loop through the fragments in turn and try to aromatize them:
    ok = True
    for i, frag in enumerate(frags):
        cp = Chem.Mol(frag)
        try:
            Chem.SanitizeMol(cp)
        except ValueError:
            matches = [x[0] for x in frag.GetSubstructMatches(Chem.MolFromSmarts(nitrogenPattern))]
            lres, indices = _recursivelyModifyNs(frag, matches)
            if not lres:
                # print 'frag %d failed (%s)'%(i,str(fragLists[i]))
                ok = False
                break
            else:
                revMap = {}
                for k, v in frag._idxMap.items():
                    revMap[v] = k
                for idx in indices:
                    oatom = m.GetAtomWithIdx(revMap[idx])
                    oatom.SetNoImplicit(True)
                    oatom.SetNumExplicitHs(1)
    if not ok:
        return None
    return m


def sanifix(m):
    if m is None:
        return None
    try:
        m.UpdatePropertyCache(False)
        cp = Chem.Mol(m.ToBinary())
        Chem.SanitizeMol(cp)
        return cp
    except ValueError as e:
        logging.debug(e, Chem.MolToSmiles(m))
        try:
            m = AdjustAromaticNs(m)
            if m is not None:
                Chem.SanitizeMol(m)
            return m
        except Exception as ee:
            logging.debug(ee, Chem.MolToSmiles(m))
            return None
    except RuntimeError as e:
        logging.debug(e, Chem.MolToSmiles(m))
        logging.info("The faulty smiles is: {}".format(Chem.MolToSmiles(m)))
        raise e
		

def get_FeatureMapScore(small_m, large_m, score_mode=FeatMaps.FeatMapScoreMode.Best):
    featLists = []
    for m in [small_m, large_m]:
        rawFeats = fdef.GetFeaturesForMol(m)
        # filter that list down to only include the ones we're intereted in
        featLists.append([f for f in rawFeats if f.GetFamily() in keep])
    fms = [FeatMaps.FeatMap(feats=x, weights=[1] * len(x), params=fmParams) for x in featLists]
    fms[0].scoreMode = score_mode
    fm_score = fms[0].ScoreFeats(featLists[1]) / min(fms[0].GetNumFeatures(), len(featLists[1]))
    return fm_score

def neutralize_aliphatic_nitrogens(m):
    m2 = Chem.EditableMol(m)
    pattern = Chem.MolFromSmarts("[N+1]")
    at_matches = m.GetSubstructMatches(pattern)
    at_matches_list = [y[0] for y in at_matches]
    #print(at_matches_list)
    if len(at_matches_list) > 0:
        for at_idx in at_matches_list:
            atom = m.GetAtomWithIdx(at_idx)
            chg = atom.GetFormalCharge()
            hcount = atom.GetTotalNumHs(includeNeighbors=True)
            #atom.SetFormalCharge(0)
            to_remove=min(hcount,chg)
            removed=0
            #atom.SetNumExplicitHs(hcount - chg)
            for a in atom.GetNeighbors():
                if a.GetAtomicNum() == 1 and removed<to_remove:
                    print(a.GetIdx())
                    m2.RemoveAtom(a.GetIdx())
                    removed=removed+1
        m3=m2.GetMol()
        for at_idx in at_matches_list:
            atom = m3.GetAtomWithIdx(at_idx)
            atom.SetFormalCharge(0)
            atom.UpdatePropertyCache()
    Chem.SanitizeMol(m3)
    return m3

def neutralize_all_nitrogens(m):
    m2 = Chem.EditableMol(m)
    pattern = Chem.MolFromSmarts("[n+1]")
    at_matches = m.GetSubstructMatches(pattern)
    at_matches_list = [y[0] for y in at_matches]
    at_matches_list.sort(reverse=True)
    #at_matches_list
    #print(at_matches_list)
    if len(at_matches_list) > 0:
        for at_idx in at_matches_list:
            atom = m.GetAtomWithIdx(at_idx)
            chg = atom.GetFormalCharge()
            hcount = atom.GetTotalNumHs(includeNeighbors=True)
            #atom.SetFormalCharge(0)
            to_remove=min(hcount,chg)
            removed=0
            #atom.SetNumExplicitHs(hcount - chg)
            if removed<to_remove:
                for a in atom.GetNeighbors():
                    if a.GetAtomicNum() == 1:
                        print(a.GetIdx())
                        m2.RemoveAtom(a.GetIdx())
                        removed=removed+1
        m3=m2.GetMol()
        for at_idx in at_matches_list:
            atom = m3.GetAtomWithIdx(at_idx)
            atom.SetFormalCharge(0)
            atom.UpdatePropertyCache()
    m4 = Chem.EditableMol(m3)
    pattern = Chem.MolFromSmarts("[N+1]")
    at_matches = m3.GetSubstructMatches(pattern)
    at_matches_list = [y[0] for y in at_matches]
    at_matches_list.sort(reverse=True)
    #print(at_matches_list)
    if len(at_matches_list) > 0:
        for at_idx in at_matches_list:
            atom = m3.GetAtomWithIdx(at_idx)
            chg = atom.GetFormalCharge()
            hcount = atom.GetTotalNumHs(includeNeighbors=True)
            #atom.SetFormalCharge(0)
            to_remove=min(hcount,chg)
            removed=0
            #atom.SetNumExplicitHs(hcount - chg)
            for a in atom.GetNeighbors():
                if a.GetAtomicNum() == 1 and removed<to_remove:
                    print(a.GetIdx())
                    m4.RemoveAtom(a.GetIdx())
                    removed=removed+1
        m5=m4.GetMol()
        for at_idx in at_matches_list:
            atom = m5.GetAtomWithIdx(at_idx)
            atom.SetFormalCharge(0)
            atom.UpdatePropertyCache()
    Chem.SanitizeMol(m5)
    return m5

def neutralize_aromatic_nitrogens(m):
    pattern = Chem.MolFromSmarts("[n+1]")
    at_matches = m.GetSubstructMatches(pattern)
    at_matches_list = [y[0] for y in at_matches]
    at_matches_list.sort(reverse=True)
    #at_matches_list
    #print(at_matches_list)
    if len(at_matches_list) > 0:
        for at_idx in at_matches_list:
            atom = m.GetAtomWithIdx(at_idx)
            hcount = atom.GetTotalNumHs(includeNeighbors=True)
            chg = atom.GetFormalCharge()
            if hcount>0 and chg==1:
                atom.SetFormalCharge(0)
                atom.UpdatePropertyCache()
    try:
        Chem.SanitizeMol(m)
    except:
        pass
    return m

def protonate_aromatic_nitrogens(m):
    pattern = Chem.MolFromSmarts("[n+1]")
    at_matches = m.GetSubstructMatches(pattern)
    at_matches_list = [y[0] for y in at_matches]
    at_matches_list.sort(reverse=True)
    #at_matches_list
    #print(at_matches_list)
    if len(at_matches_list) > 0:
        for at_idx in at_matches_list:
            atom = m.GetAtomWithIdx(at_idx)
            hcount = atom.GetTotalNumHs(includeNeighbors=True)
            chg = atom.GetFormalCharge()
            dg = atom.GetDegree()
            if hcount==0 and chg==1 and dg==3:
                atom.SetNumExplicitHs(1)
                atom.UpdatePropertyCache()
    try:
        Chem.SanitizeMol(m)
    except:
        pass
    return m



def add_nitrogen_charges(mol):
    m = Chem.MolFromMol2File(mol,sanitize=False)
    m.UpdatePropertyCache(strict=False)
    ps = Chem.DetectChemistryProblems(m)
    if not ps:
        Chem.SanitizeMol(m)
        return m
    else:
        for p in ps:
            if p.GetType()=='AtomValenceException':
                at = m.GetAtomWithIdx(p.GetAtomIdx())
                if at.GetAtomicNum()==7 and at.GetFormalCharge()==0 and at.GetExplicitValence()==4:
                    at.SetFormalCharge(1)
            if p.GetType()=='KekulizeException':
                protonate_aromatic_nitrogens(m)
                AdjustAromaticNs(m)
                ps2 = Chem.DetectChemistryProblems(m)
                for p2 in ps2:
                    if p2.GetType()=='KekulizeException':
                        neutralize_aromatic_nitrogens(m)
                        AdjustAromaticNs(m)
#               print( 'fixed:',Chem.MolToSmiles(m))
    Chem.SanitizeMol(m)
    return m

def main(ref_file, prb_file, write=True, return_all=False,          score_mode=FeatMaps.FeatMapScoreMode.Best):
    lg = RDLogger.logger()
    lg.setLevel(RDLogger.CRITICAL)
    if type(ref_file) == str:
        if os.path.splitext(ref_file)[-1] == '.sdf':
            #reflig = Chem.MolFromMolFile(ref_file, sanitize=False)
            reflig = add_nitrogen_charges(ref_file)
        elif os.path.splitext(ref_file)[-1] == '.mol2':
            #reflig = Chem.MolFromMol2File(ref_file, sanitize=False)
            reflig = add_nitrogen_charges(ref_file)
    elif type(ref_file) == rdkit.Chem.rdchem.Mol:
        reflig = ref_file

    if type(prb_file) == str:
        #if os.path.splitext(prb_file)[-1] == '.sdf':
        #    prb_mols = Chem.SDMolSupplier(prb_file, sanitize=True)
        if os.path.splitext(prb_file)[-1] == '.mol2':
            #prb_mol = Chem.MolFromMol2File(prb_file, sanitize=False)
            prb_mol = add_nitrogen_charges(prb_file)
        elif os.path.splitext(prb_file)[-1] == '.gz':
            tmp = os.path.splitext(prb_file)[0]
            if os.path.splitext(tmp)[-1] == '.sdf':
                inf = gzip.open(prb_file)
                prb_mols = Chem.ForwardSDMolSupplier(inf, sanitize=False)
    elif type(prb_file) == rdkit.Chem.rdchem.Mol:
        prb_mols = [prb_file]

    try: reflig
    except NameError:
        raise ValueError("Incorrect file format for ref lig" )
    try: prb_mol
    except NameError:
        raise ValueError("Incorrect file format for prb lig" )

    if write: w = Chem.SDWriter("%s_SuCOS_score.sdf" % os.path.splitext(prb_file)[0])
    ######prb_mols = [x for x in prb_mols if x]

    #####for prb_mol in prb_mols:

    ##############################################
    ####### Feature map
    ##############################################
    fm_score = get_FeatureMapScore(reflig, prb_mol, score_mode)
    #fm_score = np.clip(fm_score, 0, 1) #Commented out, no longer necessary if using score_mode=Best (~Marc)
    ##############################################

    #tversky_ind = rdShapeHelpers.ShapeTverskyIndex(reflig, prb_mol, 1.0, 0.0)
    #SuCOS_score = 0.5*fm_score + 0.5*tversky_ind

    protrude_dist = rdShapeHelpers.ShapeProtrudeDist(reflig, prb_mol,
            allowReordering=False)
    protrude_dist = np.clip(protrude_dist, 0, 1)
    SuCOS_score = 0.5*fm_score + 0.5*(1 - protrude_dist)

    #print ("********************************")
    #print ("SuCOS score:\t%f" % SuCOS_score)
    print ("%f" % SuCOS_score)
    #print ("********************************")
    prb_mol.SetProp("SuCOS_score", str(SuCOS_score))
    prb_mol.SetProp("Volume_score", str(1 - protrude_dist))
    prb_mol.SetProp("Feature_score", str(fm_score))
    if write:
        w.write(prb_mol)
    if return_all:
        return SuCOS_score, fm_score, (1 - protrude_dist)
    else:
        return SuCOS_score

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="run SuCOS")
    parser.add_argument('--lig1', help='the smaller/reference ligand, in sdf or mol2 file\
                        format.')
    parser.add_argument('--lig2', help='the larger/query ligand(s), in sdf or .sdf.gz\
                        file format.')
    parser.add_argument('--write', action='store_true', default=False,
                        help='writes the SuCOS score into a sdf file with\
                        suffix _SuCOS_score.sdf')
    parser.add_argument('--return_all', action='store_true', default=False)
    parser.add_argument('--score_mode', choices=['all', 'closest', 'best'],
                        help='choose the scoring mode for the feature map,\
                        default is all.')

    args = parser.parse_args()

    ref_file = args.lig1
    prb_file = args.lig2

    if args.score_mode:
        #print("mode specified:")
        if args.score_mode == 'all':
            #print ("Feature maps scoring by all. WARNING: THIS IS NOT NORMALIZED")
            score_mode = FeatMaps.FeatMapScoreMode.All
        elif args.score_mode == 'closest':
            #print ("Feature maps scoring by closest. WARNING: THIS IS NOT NORMALIZED")
            score_mode = FeatMaps.FeatMapScoreMode.Closest
        elif args.score_mode == 'best': #FeatMapScoreMode.Best normalizes SuCOS ##MARC 
            #print ("Feature maps scoring by best. This is standard.")
            score_mode = FeatMaps.FeatMapScoreMode.Best
        #else:
            #print ("This is not an option")
    else:
        #print("standard mode utilized:")
        #print ("Feature maps scoring by best. This is standard.")
        score_mode = FeatMaps.FeatMapScoreMode.Best

    main(ref_file, prb_file, score_mode, args.write, args.return_all)

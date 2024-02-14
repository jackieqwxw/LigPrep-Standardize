#!/usr/bin/env python3

import pandas as pd
import sys
import os
import rdkit
from rdkit import Chem
from rdkit.Chem.rdmolfiles import SDWriter
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem import AllChem

def standardize(mol): 
    # removeHs, disconnect metal atoms, normalize the molecule, reionize the molecule
    clean_mol = rdMolStandardize.Cleanup(mol) 
    # if many fragments, get the "parent" (the actual mol we are interested in) 
    parent_clean_mol = rdMolStandardize.FragmentParent(clean_mol)
    # try to neutralize molecule
    uncharger = rdMolStandardize.Uncharger()
    uncharged_parent_clean_mol = uncharger.uncharge(parent_clean_mol)    
    # canonical tautomers
    enumerator = rdMolStandardize.TautomerEnumerator()
    canon = enumerator.Canonicalize(uncharged_parent_clean_mol)
    return canon

def neutralize_atoms(mol):
    pattern = Chem.MolFromSmarts("[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]")
    at_matches = mol.GetSubstructMatches(pattern)
    at_matches_list = [y[0] for y in at_matches]
    if len(at_matches_list) > 0:
        for at_idx in at_matches_list:
            atom = mol.GetAtomWithIdx(at_idx)
            chg = atom.GetFormalCharge()
            hcount = atom.GetTotalNumHs()
            atom.SetFormalCharge(0)
            atom.SetNumExplicitHs(hcount - chg)
            atom.UpdatePropertyCache()
    return mol

def protonate_ali_amino(mol):#charging sp3 Nitrogens with +1 
    for at in mol.GetAtoms(): 
            if at.GetAtomicNum() == 7 and at.GetHybridization()==rdkit.Chem.rdchem.HybridizationType.SP3 and at.GetFormalCharge()==0:
                at.SetFormalCharge(1)
    return mol

def deprotonation_cooh(mol):
    deprotonate_cooh  =  AllChem.ReactionFromSmarts("[C:1](=[O:2])-[OH1:3]>>[C:1](=[O:2])-[O-H0:3]")
    m_deprot  =  deprotonate_cooh.RunReactants((mol,))
    return m_deprot[0][0]  if  m_deprot  else  mol 

## Headers
print("##############################################################")
print("#  Ligand preparation / standardization for virtual screening\n#")
print("#  Computational Biophysics & Drug Design (Kireev lab)\n#")
print("#  Developed by Akila Mettu and Xiaowen Wang")
print("##############################################################\n")

## command line arguments setting
if len(sys.argv) < 3:
    print("Usage: python3 standardize.py input_file output_file")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

## read input files from smiles
read_smiles = Chem.SmilesMolSupplier(input_file,titleLine=False,sanitize=True)
mols = [x for x in read_smiles if x is not None]
legends = [x.GetProp("_Name") for x in read_smiles if x is not None]
com = list(zip(mols,legends))

## write molecules to smiles
failmol = 0
modnum = 0
with open(output_file, "w") as fout:
    for i,j in enumerate(com):
        mol, label = com[i]
        try:
            st=standardize(mol)
            st.UpdatePropertyCache()
            neu=neutralize_atoms(st)
            pro=protonate_ali_amino(neu)
            dep=deprotonation_cooh(pro) 
            Chem.SanitizeMol(dep)
            org = Chem.MolToSmiles(mol)
            mod = Chem.MolToSmiles(dep,isomericSmiles=True)
            if org != mod:
               modnum += 1     
            dep.SetProp('_Name',label)
            ## preserve the stereo if specified 
            fout.write(Chem.MolToSmiles(dep, isomericSmiles=True))
            fout.write(' '+label+'\n')
        except Chem.KekulizeException:
            #print(f"Unable to kekulize molecule with name {label}. Skipping...")
            failmol += 1

current_directory = os.getcwd()
print("%s" %format(i+1,','),"ligands read from %s" %('['+current_directory+'/'+input_file+']'))
print("%s" %format(modnum,','),"ligands modified")
print("%s" %format(i+1-failmol,','),"ligands written to %s" %('['+current_directory+'/'+output_file+']'))
#print("%s" %format(failmol,','),"ligands failed")
fout.close()

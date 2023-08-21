from rdkit import Chem
from rdkit.Chem import Descriptors

def find_prop(smile):
  mol = Chem.MolFromSmiles(smile)

  # Calculate molecular properties
  formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
  molar_weight = (Descriptors.MolWt(mol))
  num_h_bonds = (Descriptors.NumHDonors(mol))
  boiling_pt = (Descriptors.TPSA(mol))  # Approximation for boiling point
  melting_pt = (Descriptors.MolLogP(mol) ) # Approximation for melting point
  log_p = (Descriptors.MolLogP(mol))
  num_h_bond_acceptors = (Descriptors.NumHAcceptors(mol))
  num_h_bond_donors = ( Descriptors.NumHDonors(mol))

  ans = True
  if(molar_weight <= 500 and log_p <= 5 and num_h_bond_donors <= 5 and num_h_bond_acceptors <= 10):
    ans = True
  else:
    ans = False

  return {"formula":formula, "molar_weight":"{:.4f}".format(molar_weight), "num_h_bonds":"{:.4f}".format(num_h_bonds),"boiling_pt":"{:.4f}".format(boiling_pt),"melting_pt":"{:.4f}".format(melting_pt), "log_p":"{:.4f}".format(log_p), "num_h_bond_acceptors":"{:.4f}".format(num_h_bond_acceptors), "num_h_bond_donors":"{:.4f}".format(num_h_bond_donors), "Lipinski":ans}

def find_efficacy(value):
  if(value >= 0.55):
    return True
  return False

def find_effectiveness(pic50, receptors):

  ic50_res = ""
  if(pic50 >= 6):
      ic50_res = "Acceptable and Can be effective for Diabeties"
  else:
      ic50_res = "May not be as effective for Diabeties"

  map = {"IC50_res": ic50_res}

  score = 0
  
  receptor_names = ['NR.AhR', 'NR.AR', 'NR.AR.LBD', 'NR.Aromatase', 'NR.ER','NR.ER.LBD','NR.PPAR.gamma','SR.ARE','SR.ATAD5','SR.HSE','SR.MMP','SR.p53']
  recep = []
  for i in range(12):
    if( i != 6):
      if(receptor_names[i] == 0):
        map[f'rep{i+1}_res'] = "Acceptable"
      else:
        score += 1
        recep.append(receptors[i])
        map[f'rep{i+1}_res'] = "Not Acceptable"
    else:
      score += 1
      map[f'rep{i+1}_res'] = "Acceptable"
   
  effects = ""
  recep_str = [str(item) for item in recep]
  if(score >= 2):
    effects = "Drug may not be as effective as it is acting on different receptors" 
  else:
    effects = "None and is Acceptable for Diabeties"

  return {"IC50_res": ic50_res, "Effects": effects}

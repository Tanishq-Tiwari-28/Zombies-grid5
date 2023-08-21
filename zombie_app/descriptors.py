from rdkit import Chem
from rdkit.Chem import Descriptors3D, AllChem
from rdkit.ML.Descriptors.MoleculeDescriptors import MolecularDescriptorCalculator
import pandas as pd
from rdkit import Chem
from mordred import Calculator, descriptors
from rdkit.Chem import GetMolFrags
import mordred
from IPython.display import display
from PIL import Image
from rdkit.Chem import Draw
def more_descriptors(smile):
  mapping = {}

  mol = Chem.MolFromSmiles(smile)

  chosen_descriptors = ['BalabanJ', 'BertzCT', 'Chi0', 'Chi0n', 'Chi0v', 'Chi1', 'Chi1n', 'Chi1v', 'Chi2n', 'Chi2v', 'Chi3n', 'Chi3v', 'Chi4n', 'Chi4v', 'EState_VSA1', 'EState_VSA10', 'EState_VSA11', 'EState_VSA2', 'EState_VSA3', 'EState_VSA4', 'EState_VSA5', 'EState_VSA6', 'EState_VSA7', 'EState_VSA8', 'EState_VSA9', 'ExactMolWt', 'FpDensityMorgan1', 'FpDensityMorgan2', 'FpDensityMorgan3', 'FractionCSP3', 'HallKierAlpha', 'HeavyAtomCount', 'HeavyAtomMolWt', 'Ipc', 'Kappa1', 'Kappa2', 'Kappa3', 'LabuteASA', 'MaxAbsEStateIndex', 'MaxAbsPartialCharge', 'MaxEStateIndex', 'MaxPartialCharge', 'MinAbsEStateIndex', 'MinAbsPartialCharge', 'MinEStateIndex', 'MinPartialCharge', 'MolLogP', 'MolMR', 'MolWt', 'NHOHCount', 'NOCount', 'NumAliphaticCarbocycles', 'NumAliphaticHeterocycles', 'NumAliphaticRings', 'NumAromaticCarbocycles', 'NumAromaticHeterocycles', 'NumAromaticRings', 'NumHAcceptors', 'NumHDonors', 'NumHeteroatoms', 'NumRadicalElectrons', 'NumRotatableBonds', 'NumSaturatedCarbocycles', 'NumSaturatedHeterocycles', 'NumSaturatedRings', 'NumValenceElectrons', 'PEOE_VSA1', 'PEOE_VSA10', 'PEOE_VSA11', 'PEOE_VSA12', 'PEOE_VSA13', 'PEOE_VSA14', 'PEOE_VSA2', 'PEOE_VSA3', 'PEOE_VSA4', 'PEOE_VSA5', 'PEOE_VSA6', 'PEOE_VSA7', 'PEOE_VSA8', 'PEOE_VSA9', 'RingCount', 'SMR_VSA1', 'SMR_VSA10', 'SMR_VSA2', 'SMR_VSA3', 'SMR_VSA4', 'SMR_VSA5', 'SMR_VSA6', 'SMR_VSA7', 'SMR_VSA8', 'SMR_VSA9', 'SlogP_VSA1', 'SlogP_VSA10', 'SlogP_VSA11', 'SlogP_VSA12', 'SlogP_VSA2', 'SlogP_VSA3', 'SlogP_VSA4', 'SlogP_VSA5', 'SlogP_VSA6', 'SlogP_VSA7', 'SlogP_VSA8', 'SlogP_VSA9', 'VSA_EState1', 'VSA_EState10', 'VSA_EState2', 'VSA_EState3', 'VSA_EState4', 'VSA_EState5', 'VSA_EState6', 'VSA_EState7', 'VSA_EState8', 'VSA_EState9', 'fr_Al_COO', 'fr_Al_OH', 'fr_Al_OH_noTert', 'fr_ArN', 'fr_Ar_COO', 'fr_Ar_N', 'fr_Ar_NH', 'fr_Ar_OH', 'fr_COO', 'fr_COO2', 'fr_C_O', 'fr_C_O_noCOO', 'fr_C_S', 'fr_HOCCN', 'fr_Imine', 'fr_NH0', 'fr_NH1', 'fr_NH2', 'fr_N_O', 'fr_Ndealkylation1', 'fr_Ndealkylation2', 'fr_Nhpyrrole', 'fr_SH', 'fr_aldehyde', 'fr_alkyl_carbamate', 'fr_alkyl_halide', 'fr_allylic_oxid', 'fr_amide', 'fr_amidine', 'fr_aniline', 'fr_aryl_methyl', 'fr_azide', 'fr_azo', 'fr_barbitur', 'fr_benzene', 'fr_benzodiazepine', 'fr_bicyclic', 'fr_diazo', 'fr_dihydropyridine', 'fr_epoxide', 'fr_ester', 'fr_ether', 'fr_furan', 'fr_guanido', 'fr_halogen', 'fr_hdrzine', 'fr_hdrzone', 'fr_imidazole', 'fr_imide', 'fr_isocyan', 'fr_isothiocyan', 'fr_ketone', 'fr_ketone_Topliss', 'fr_lactam', 'fr_lactone', 'fr_methoxy', 'fr_morpholine', 'fr_nitrile', 'fr_nitro', 'fr_nitro_arom', 'fr_nitro_arom_nonortho', 'fr_nitroso', 'fr_oxazole', 'fr_oxime', 'fr_para_hydroxylation', 'fr_phenol', 'fr_phenol_noOrthoHbond', 'fr_phos_acid', 'fr_phos_ester', 'fr_piperdine', 'fr_piperzine', 'fr_priamide', 'fr_prisulfonamd', 'fr_pyridine', 'fr_quatN', 'fr_sulfide', 'fr_sulfonamd', 'fr_sulfone', 'fr_term_acetylene', 'fr_tetrazole', 'fr_thiazole', 'fr_thiocyan', 'fr_thiophene', 'fr_unbrch_alkane', 'fr_urea', 'qed']
  mol_descriptor_calculator = MolecularDescriptorCalculator(chosen_descriptors)
  list_of_descriptor_vals = list(mol_descriptor_calculator.CalcDescriptors(mol))

  for i in range(len(list_of_descriptor_vals)):
    mapping[chosen_descriptors[i]] = list_of_descriptor_vals[i]

  mol = Chem.AddHs(mol)
  AllChem.EmbedMolecule(mol, maxAttempts=500)

  descriptor_names = ['Asphericity', 'Eccentricity', 'InertialShapeFactor', 'NPR1', 'NPR2', 'PMI1', 'PMI2', 'PMI3', 'RadiusOfGyration', 'SpherocityIndex']

  for descriptor_name in descriptor_names:
    try:
        descriptor_func = getattr(Descriptors3D, descriptor_name)
        descriptor_value = descriptor_func(mol)
        mapping[descriptor_name] = descriptor_value
    except:
        print(f"Descriptor '{descriptor_name}' not found in RDKit.")

  descriptor_names = ['CalcNumBridgeheadAtoms', 'CalcNumHBA', 'CalcNumHBD', 'CalcNumHeavyAtoms', 'CalcNumHeterocycles', 'CalcNumLipinskiHBA', 'CalcNumLipinskiHBD', 'CalcNumRings', 'CalcNumSpiroAtoms', 'CalcPBF', 'CalcPhi', ]

  for descriptor_name in descriptor_names:
    try:
        descriptor_func = getattr(Descriptors3D.rdMolDescriptors, descriptor_name)
        descriptor_value = descriptor_func(mol)
        mapping[descriptor_name[4:]] = descriptor_value
    except:
        print(f"Descriptor '{descriptor_name}' not found in RDKit.")


  for i in range(5, 10):
    mapping['Chi'+str(i)+'v']=Descriptors3D.rdMolDescriptors.CalcChiNv(mol,i)
    mapping['Chi'+str(i)+'n']=Descriptors3D.rdMolDescriptors.CalcChiNn(mol,i)


  try:
    crippen = Descriptors3D.rdMolDescriptors.CalcCrippenDescriptors(mol)
    for i in range(1, len(crippen)+1):
      mapping[f'crippen{i}'] = crippen[i-1]

    Morse = Descriptors3D.rdMolDescriptors.CalcMORSE(mol)
    for i in range(1, len(Morse)+1):
      mapping[f'Morse{i}'] = Morse[i-1]

    RDF = Descriptors3D.rdMolDescriptors.CalcRDF(mol)
    for i in range(1, len(RDF)+1):
      mapping[f'RDF{i}'] = RDF[i-1]

    WHIM = Descriptors3D.rdMolDescriptors.CalcWHIM(mol)
    for i in range(1, len(WHIM)+1):
      mapping[f'WHIM{i}'] = WHIM[i-1]
  except:
    print("hello")

  return mapping

def smiles_des(smile):
    print("in smilesdes")
    calc = Calculator(descriptors, ignore_3D=False)
    
    mols = []
    mol = Chem.MolFromSmiles(smile)

    if mol is not None:
            fragments = GetMolFrags(mol)
            if len(fragments) == 1:  
                mols.append(mol)
    else:
        return None
                    
    
    descriptors_df = calc.pandas(mols)
    mapping = more_descriptors(smile)
    
    similar_keys = []
    for key in descriptors_df.columns:
        if(key in mapping):
            similar_keys.append(key)
            

    for key in similar_keys:
        if key in mapping:
            del mapping[key]
            
    
    
    # appending columns to df
    for column_name, value in mapping.items():
        descriptors_df[column_name] = value
    # print(descriptors_df)
    # print("hsj;gsfh;shf;wsifh;sgh;fgh" , type(descriptors_df))
    # descriptors_df = pd.concat([descriptors_df, pd.DataFrame(mapping)], axis=1)
    # print("df" , descriptors_df)
    row_dict = descriptors_df.to_dict(orient='records')[0]


    for key in row_dict:
        if isinstance(row_dict[key], mordred.error.Error) or isinstance(row_dict[key], mordred.error.Missing):  # Replace with the actual type
            row_dict[key] = 0
    

    return row_dict

def des_main(descriptors_df, file_path):
    
    with open(file_path, 'r') as f:
        array = f.read()
        f.close()

    array = array.split()

    values_to_pass = []
    try:
      for i in array:
          values_to_pass.append(descriptors_df[i])
    except:
       return "Molecule NOT FOUND"
    return values_to_pass

from io import BytesIO   

def draw_structure(smiles):
    mol = Chem.MolFromSmiles(smiles)
    print("Cbchdsbvjkd", smiles)
    img = Draw.MolToImage(mol)
    img_path = "C:/Users/tanis/OneDrive/Desktop/Tanishq/Flipkart/Zombies/static/img/molecule.png"
    img.save(img_path)

  
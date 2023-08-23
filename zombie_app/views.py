from django.shortcuts import render, redirect , HttpResponse
import joblib
# model = joblib.load('./randomForest.joblib')
import pandas as pd
import numpy as np
from rdkit import Chem
from mordred import Calculator, descriptors
from rdkit.Chem import GetMolFrags
# from chembl_webresource_client.new_client import new_client
# import chembl_webresource_client
# from .chembl import ChEMBLClientSingleton
from .descriptors import draw_structure
from .result import find_prop , find_efficacy , find_effectiveness
import os
from reportlab.lib.pagesizes import letter
from reportlab.pdfgen import canvas
import io
from urllib.parse import unquote
from reportlab.pdfgen import canvas
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfbase import pdfmetrics
from reportlab.lib import colors
from PIL import Image as PILImage

from django.http import FileResponse
import time
from django.core.cache import cache


# def load_ml_model():
#     # Try to get the model from the cache
#     model = cache.get('cached_ml_model')

#     if model is None:
#         # Model is not in the cache, load it from joblib file
#         model = joblib.load('./randomForest.joblib')
#         # Store the model in the cache for future requests
#         cache.set('cached_ml_model', model)
#     else:
#         # Use the cached model
#         pass

#     return model

def smiles_descriptors(smile):
    print('in smiles_des')
    from .descriptors import smiles_des
    mapping = smiles_des(smile)
    mapping['isomeric_smile'] = smile
    # print("mapping"  ,mapping)
    # print(des)
    return mapping


def chembl_smiles(id):
    print('in chembl_smiles')
    print(id)
    # print(mols)
    try:
        from .chembl import molecule
        res = molecule.filter(chembl_id=id).only(['molecule_chembl_id', 'pref_name', 'molecule_structures'])
        print(res)
        smile = res[0]['molecule_structures']['canonical_smiles']
        mol = Chem.MolFromSmiles(smile)
        isomeric=''
        if mol:
            isomeric = Chem.MolToSmiles(mol, isomericSmiles=True)
        # print(smile)
        df = smiles_descriptors(isomeric)
        print(df)
        df['isomeric_smile'] = isomeric
        return df
        # if df != None:
            
        # else:
        #     return ''

    except:
        return ''
  
def pubchem_smiles(compound_name):
    print("in pubchem")
    import pubchempy as pcp
    try:
        compound = pcp.get_compounds(compound_name, 'name')[0]
        print(compound)
        smiles = compound.to_dict(properties=['isomeric_smiles'])['isomeric_smiles']
        df =  smiles_descriptors(smiles)
        df['isomeric_smile'] = smiles
        return df
    except IndexError:
        return "Molecule with the given name doesn't exist"

    
def name_chembl(name):
    print('in names_chembl')

    try:
        from .chembl import molecule
        print("Trying to query ChEMBL with name:", name)
        mols = molecule.filter(molecule_synonyms__molecule_synonym=name).only('molecule_chembl_id')
        print("Query result:", mols)
        if(not mols):
            pubchem_smiles(name)
        if len(mols) >= 1:
            print("Found molecule, getting smiles...")
            res = chembl_smiles(mols[0]['molecule_chembl_id'])
            # print("Smiles result:", res)
            return res
            # if res != None:
            #     return res
            # else:
            #     return "Smiles not found for the molecule"  
        else:
            return "Molecule with the given name doesn't exist"
    
    except Exception as e:
        return "Error: " + str(e)
          

# print("ame function" , name_chembl("Sildenafil"))



def home(request):
    
    result = None
    pdf_data = None
    # pdf_viewed = request.session.get('pdf_viewed', False)
    try:

        print('in view')
        if request.method == 'POST':
            # if 'resubmit' in request.POST:
            #     # Clear PDF viewed status and redirect to the form page
            #     request.session.pop('pdf_viewed', None)
            #     return redirect('home')
            animation = True
            input_data = request.POST.get('input')
            input_type = request.POST.get('input-type')
            if input_type == 'drug_name':
                result = name_chembl(input_data)
    
                draw_structure(result['isomeric_smile'])
                
            elif input_type == 'chembl_id':
                result = chembl_smiles(input_data)
                draw_structure(result['isomeric_smile'])
                # print(result)
            elif input_type == 'smiles':
                result = smiles_descriptors(input_data)
                # print(result)
                draw_structure(input_data)
            # result_dict = result.to_dict()
            # print(type(result_dict))
            # print(result_dict)
            request.session['result'] = result
            request.session['input_data'] = input_data
            request.session['input_type'] = input_type
            # request.session['image'] = image
            # request.session.save()
            print("to resultview")
            # time.sleep(2)       
            return redirect('result') 
    except:
        error  = 1
        return render(request , 'index.html' , {'error':error})
            
    return render(request, 'index.html') #'pdf_viewed': pdf_viewedx 



def result(request):
    error = 0
    from .descriptors import des_main
    model1 = joblib.load('.\\joblib\\Efficacy_Model.joblib')
    model2 = joblib.load('.\\joblib\\NR.AhR.joblib')
    model3 = joblib.load('.\\joblib\\NR.AR.joblib')
    model4 = joblib.load('.\\joblib\\NR.AR.LBD.joblib')
    model5 = joblib.load('.\\joblib\\NR.Aromatase.joblib')
    model6 = joblib.load('.\\joblib\\NR.ER.joblib')
    model7 = joblib.load('.\\joblib\\NR.ER.LBD.joblib')
    model8 = joblib.load('.\\joblib\\NR.PPAR.gamma.joblib') #diabeties
    model9 = joblib.load('.\\joblib\\SR.ARE.joblib')
    model10 = joblib.load('.\\joblib\\SR.ATAD5.joblib')
    model11 = joblib.load('.\\joblib\\SR.HSE.joblib')
    model12 = joblib.load('.\\joblib\\SR.MMP.joblib')
    model13 = joblib.load('.\\joblib\\SR.p53.joblib')
    model14 = joblib.load('.\\joblib\\Pic50_regressor.joblib')
    print("model8" , model8)
    result = request.session.get('result')
    print("fsdfgsdgs" , type(result))   
    path1 = ".\\modelspath\\Efficacy_Model.txt"
    path2 = ".\\modelspath\\NR.AhR.txt"
    path3 = ".\\modelspath\\NR.AR.txt"
    path4 = ".\\modelspath\\NR.AR.LBD.txt"
    path5 = ".\\modelspath\\NR.Aromatase.txt"
    path6 = ".\\modelspath\\NR.ER.txt"
    path7 = ".\\modelspath\\NR.ER.LBD.txt"
    path8 = ".\\modelspath\\NR.PPAR.gamma.txt"
    path9 = ".\\modelspath\\SR.ARE.txt"
    path10 = ".\\modelspath\\SR.ATAD5.txt"
    path11 = ".\\modelspath\\SR.HSE.txt"
    path12 = ".\\modelspath\\SR.MMP.txt"
    path13 = ".\\modelspath\\SR.p53.txt"
    path14  = ".\\modelspath\\Pic50_features.txt"

    try:
        des1 = des_main(result, path1)
        des2 = des_main(result, path2)
        des3 = des_main(result, path3)
        des4 = des_main(result, path4)
        des5 = des_main(result, path5)
        des6 = des_main(result, path6)
        des7 = des_main(result, path7)
        des8 = des_main(result, path8)
        print("des8" , des8)
        des9 = des_main(result, path9)
        des10 = des_main(result, path10)
        des11 = des_main(result, path11)
        des12 = des_main(result, path12)
        des13 = des_main(result, path13)
        des14 = des_main(result , path14)


        input_data_des1 = np.array(des1).reshape(1, len(des1))
        input_data_des2 = np.array(des2).reshape(1, len(des2))
        input_data_des3 = np.array(des3).reshape(1, len(des3))
        input_data_des4 = np.array(des4).reshape(1, len(des4))
        input_data_des5 = np.array(des5).reshape(1, len(des5))
        input_data_des6 = np.array(des6).reshape(1, len(des6))
        input_data_des7 = np.array(des7).reshape(1, len(des7))
        input_data_des8 = np.array(des8).reshape(1, len(des8))
        print("input data 8" , input_data_des1)
        input_data_des9 = np.array(des9).reshape(1, len(des9))
        input_data_des10 = np.array(des10).reshape(1, len(des10))
        input_data_des11 = np.array(des11).reshape(1, len(des11))
        input_data_des12 = np.array(des12).reshape(1, len(des12))
        input_data_des13 = np.array(des13).reshape(1, len(des13))    
        input_data_des14 = np.array(des14).reshape(1, len(des14))  
        # print(input_data_des14)

        predictions_des1 = model1.predict_proba(input_data_des1)
        predictions_des2 = model2.predict(input_data_des2)
        predictions_des3 = model3.predict(input_data_des3)
        predictions_des4 = model4.predict(input_data_des4)
        predictions_des5 = model5.predict(input_data_des5)
        predictions_des6 = model6.predict(input_data_des6)
        predictions_des7 = model7.predict(input_data_des7)
        predictions_des8 = model8.predict(input_data_des8)
        predictions_des9 = model9.predict(input_data_des9)
        predictions_des10 = model10.predict(input_data_des10)
        predictions_des11 = model11.predict(input_data_des11)
        predictions_des12 = model12.predict(input_data_des12)
        predictions_des13 = model13.predict(input_data_des13)    
        predictions_des14 = model14.predict(input_data_des14)
        
        print("Predictions for des1:", predictions_des1)
        print("Predictions for des2:", predictions_des2)
        print("Predictions for des3:", predictions_des3)
        print("Predictions for des4:", predictions_des4)
        print("Predictions for des5:", predictions_des5)
        print("Predictions for des6:", predictions_des6)
        print("Predictions for des7:", predictions_des7)
        print("Predictions for des8:", predictions_des8)
        print("Predictions for des9:", predictions_des9)
        print("Predictions for des10:", predictions_des10)
        print("Predictions for des11:", predictions_des11)
        print("Predictions for des12:", predictions_des12)
        print("Predictions for des13:", predictions_des13)
        print("Predictions for des14:", predictions_des14)
        # Retrieve data from session
        result = request.session.get('result')
        input_data = request.session.get('input_data')
        input_type = request.session.get('input_type')
        # image = request.session.get('image')

        final_list = [
            predictions_des1,
            predictions_des2,
            predictions_des3,
            predictions_des4,
            predictions_des5,
            predictions_des6,
            predictions_des7,
            predictions_des8,
            predictions_des9,
            predictions_des10,
            predictions_des11,
            predictions_des12,
            predictions_des13,
            predictions_des14

        ]
        final_list = [predictions.tolist() for predictions in final_list]

        if(input_type == 'smiles'):
            prop = find_prop(input_data)
        prop = find_prop(result['isomeric_smile'])   
        efficacy = find_efficacy(predictions_des1[0][0])
        if(efficacy == True):
            efficacy = "Drug Can be effacacious with confidence score of " + str(predictions_des1[0][0])
        else:
            efficacy = "Drug Is Not effacious as the confidence score is " + str(predictions_des1[0][0])
        
        effectiveness = find_effectiveness(predictions_des14[0] , final_list)
        print(type(effectiveness))
        print((effectiveness))
        print((efficacy))
        print(prop)
        print("in result")
        print("fsdfgsdgs" , input_data)
        # Clear session data after retrieving it
        request.session['prop']  = prop
        request.session['efficacy']  = efficacy
        request.session['effectiveness']  = effectiveness
        request.session['predictions']  = final_list
        
    # request.session.pop('image', None)

    except:
        error = 1
        return render(request , 'home.html' , {'error':error})
    return render(request, 'result.html', {
        'result': result,
        'input_data': input_data,
        'input_type': input_type,
        'predictions' : final_list,
        # 'image' : image,
        'efficacy' : efficacy,
        'effectiveness' : effectiveness,
        'prop' : prop
    })

from django.conf import settings

def download_pdf(request):
    print("in download")

   
    story = []
    prop = request.session.get('prop')
    efficacy = request.session.get('efficacy')
    effectiveness = request.session.get('effectiveness')
    stored_predictions = request.session.get('predictions', [])
    predictions = [np.array(pred) for pred in stored_predictions]
    
    request.session.pop('result', None)
    request.session.pop('input_data', None)
    request.session.pop('prop' , None)
    request.session.pop('efficacy' , None)
    request.session.pop('effectiveness' , None)
    request.session.pop('predictions' , None)


    buffer = io.BytesIO()
    p = canvas.Canvas(buffer, pagesize=letter)

    documentTitle = 'Drug Quality Testing Report'
    title = 'Drug Quality Testing'
    subTitle = 'Test Your Drug, Be Safe'
    
    image_path = (os.path.join(settings.BASE_DIR, 'static/img/molecule.png'))
    image = PILImage.open(image_path)
    p.drawImage(image_path, 350, 500, width=200, height=200)


    textLines1 = [
        f"{key}: {value}" for key, value in prop.items()
    ]
    textLines2 = efficacy.split('\n')
    textLines3 = [
        f"{key}: {value}" for key, value in effectiveness.items()
    ]
    receptors = ['Efficacy_Model', 'NR.AhR', 'NR.AR', 'NR.AR.LBD', 'NR.Aromatase', 'NR.ER', 'NR.ER.LBD', 'NR.PPAR.gamma', 'SR.ARE', 'SR.ATAD5', 'SR.HSE', 'SR.MMP', 'SR.p53', 'Pic50_value']

    textLines4 = [
    f"{receptors[i]}: {', '.join(map(str, predictions[i]))}"
    for i in range(len(predictions))
    ]

    p.setTitle(documentTitle)
    p.setFont("Helvetica-Bold", 18)
    p.drawString(100, 750, title)
    p.setFont("Helvetica", 12)
    p.drawString(100, 720, subTitle)
    p.line(100, 710, 500, 710)
    
    y = 680
    for line in textLines1:
        p.drawString(100, y, line)
        y -= 20
    for line in textLines2:
        p.drawString(100, y, line)
        y -= 20
    for line in textLines3:
        p.drawString(100, y, line)
        y -= 20
    for line in textLines4:
        p.drawString(100, y, line)
        y -= 20
    p.save()

    buffer.seek(0)
    return FileResponse(buffer, as_attachment=True, filename="drug_quality_report.pdf")
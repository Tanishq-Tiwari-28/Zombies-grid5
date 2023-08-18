from django.shortcuts import render, redirect , HttpResponse
import joblib
# model = joblib.load('./randomForest.joblib')
import pandas as pd
from rdkit import Chem
from mordred import Calculator, descriptors
from rdkit.Chem import GetMolFrags
# from chembl_webresource_client.new_client import new_client
# import chembl_webresource_client
# from .chembl import ChEMBLClientSingleton

from reportlab.lib.pagesizes import letter
from reportlab.pdfgen import canvas
import io
from urllib.parse import unquote
from reportlab.pdfgen import canvas
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfbase import pdfmetrics
from reportlab.lib import colors
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
    calc = Calculator(descriptors, ignore_3D=False)
    mols = []
    mol = Chem.MolFromSmiles(smile)
    if mol != None:
            fragments = GetMolFrags(mol)
            if len(fragments) == 1:  # Exclude molecules with multiple fragments
                mols.append(mol)
    else:
        return ''
    descriptors_df = calc.pandas(mols)
    print(descriptors_df)
    return descriptors_df


def chembl_smiles(id):
    print(id)
    try:
        from .chembl import molecule
        res = molecule.filter(chembl_id=id).only(['molecule_chembl_id', 'pref_name', 'molecule_structures'])
        smile = res[0]['molecule_structures']['canonical_smiles']
        print(smile)
        df = smiles_descriptors(smile)
        print(df)
        if df != '':
            return df
        else:
            return ''

    except:
        return ''
  

def name_chembl(name):
    try:
        from .chembl import molecule
        # print(molecule)
        mols = molecule.filter(molecule_synonyms__molecule_synonym__iexact=name).only('molecule_chembl_id')
        print(mols)
        if len(mols) >= 1:
            res = chembl_smiles(mols[0]['molecule_chembl_id'])
            print(res)
            if (res == ''):
                return res
            else:
                return "Smiles not found for the molecule"  
        else:
            return "Molecule with the given name doesn't exist"

    except Exception as e:
        return "Error: " + str(e)

          

# name_chembl("Sildenafil")


def home(request):
    result = None
    pdf_data = None
    # pdf_viewed = request.session.get('pdf_viewed', False)
    print('in view')
    animation = False
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
            print(result)
        elif input_type == 'chembl_id':
            result = chembl_smiles(input_data)
            print(result)
        elif input_type == 'smiles':
            result = smiles_descriptors(input_data)
            print(result)
        animation = False
        request.session['result'] = result
        request.session['input_data'] = input_data
        request.session['input_type'] = input_type
        request.session['animation'] = animation
        # time.sleep(2)
        return redirect('result')   
    return render(request, 'index.html') #'pdf_viewed': pdf_viewedx 


def download_pdf(request):
    model = joblib.load('./randomForest.joblib')
    print(model)
    result = request.GET.get('result')
    input_data = request.GET.get('input_data')
    input_type= request.GET.get('input_type')
    print(result , input_data , input_type)
    if result:
        # Generate PDF
        buffer = io.BytesIO()
        p = canvas.Canvas(buffer, pagesize=letter)

        documentTitle = 'Drug Quality Testing Report'
        title = 'Drug Quality Testing'
        subTitle = 'Test Your Drug, Be Safe'
        textLines = [
            f"Input Type: {input_type}",
            f"Input Data: {input_data}",
        ]
        p.setTitle(documentTitle)
        p.setFont("Helvetica-Bold", 18)
        p.drawString(100, 750, title)
        p.setFont("Helvetica", 12)
        p.drawString(100, 720, subTitle)
        p.line(100, 710, 500, 710)
        y = 680
        for line in textLines:
            p.drawString(100, y, line)
            y -= 20
        p.save()
        buffer.seek(0)
        return FileResponse(buffer, as_attachment=True, filename="hello.pdf")

def result(request):
    # Retrieve data from session
    result = request.session.get('result')
    input_data = request.session.get('input_data')
    input_type = request.session.get('input_type')
    animation = request.session.get('animation')
    print("in result")
    print(result)
    # Clear session data after retrieving it
    request.session.pop('result', None)
    request.session.pop('input_data', None)
    request.session.pop('input_type', None)
    request.session.pop('animation', None)

    return render(request, 'result.html', {
        'result': result,
        'input_data': input_data,
        'input_type': input_type,
        'animation': animation,
    })

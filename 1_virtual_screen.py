import streamlit as st


import pandas as pd
import base64
import numpy as np 
import pickle
import streamlit as st
from PIL import Image
import matplotlib.pyplot as plt

#-------Deep Learning-------
from rdkit import Chem
from rdkit.Chem import MACCSkeys
from rdkit.Chem.AllChem import GetMorganFingerprintAsBitVect
import json 

with open(f'id2nm.json', 'r') as f:
    id2name = json.load(f)
    
with open(f'nm2id.json', 'r') as f:
    name2id = json.load(f)
    
    
s2m = lambda s: Chem.MolFromSmiles(s)
m2s = lambda m: Chem.MolToSmiles(m)
canonicalise = lambda s: m2s(s2m(s)) 
get_dict_query = lambda q : data[str(int(q))]
get_maccs = lambda smiles: np.array(MACCSkeys.GenMACCSKeys(s2m(smiles)))
get_mrgn2 = lambda smiles: np.array(GetMorganFingerprintAsBitVect(s2m(smiles), 2, nBits=2048))
get_mrgn3 = lambda smiles: np.array(GetMorganFingerprintAsBitVect(s2m(smiles), 3, nBits=2048))



CUDA_VISIBLE_DEVICES=""



def build_isactive_model(dataset, target,balanced):
    
    print('is active model building...')
    
    #prediction using maccs fingerprint
    if balanced=="unbalance":
        # Apply RF to make predictions
        load_model_1 = pickle.load(open(f'models/{balanced}/{target}/{target}_RF_maccs_is_active.pkl', 'rb'))
    else:
        load_model_1 = pickle.load(open(f'models/{balanced}/{target}/{target}_RF_opt_balanced_maccs_is_active.pkl', 'rb'))

    prediction_1 = load_model_1.predict(dataset.mccs_fp.tolist())
    _RF_1 = f"{id2name[target]}_RF_1"
    predicted_active_1 = pd.DataFrame(prediction_1, columns=[_RF_1])
    
    if balanced=="unbalance":
        # Apply ANN to make predictions
        load_model_3 = pickle.load(open(f'models/{balanced}/{target}/{target}_ANN_maccs_is_active.pkl', 'rb'))
    else:
        load_model_3 = pickle.load(open(f'models/{balanced}/{target}/{target}_ANN_opt_balanced_maccs_is_active.pkl', 'rb'))

    
    prediction_3 = load_model_3.predict(dataset.mccs_fp.tolist())
    
    _ANN_1 = f"{id2name[target]}_ANN_1"
    predicted_active_3 = pd.DataFrame(prediction_3, columns=[_ANN_1])
    predicted_active_df_3 = predicted_active_1.join(predicted_active_3)
    
    #Apply model using Morgan fingerprint
    if balanced=="unbalance":
        # Apply RF to make predictions
        load_model_4 = pickle.load(open(f'models/{balanced}/{target}/{target}_RF_morgan2_is_active.pkl', 'rb'))
    else:
        load_model_4 = pickle.load(open(f'models/{balanced}/{target}/{target}_RF_opt_balanced_morgan2_is_active.pkl', 'rb'))
    prediction_4 = load_model_4.predict(dataset.mrgn_fp.tolist())
    _RF_2 = f"{id2name[target]}_RF_2"
    predicted_active_4 = pd.DataFrame(prediction_4, columns=[_RF_2])
    predicted_active_df_4 = predicted_active_df_3.join(predicted_active_4)

    if balanced=="unbalance":
        # Apply ANN to make predictions
        load_model_5 = pickle.load(open(f'models/{balanced}/{target}/{target}_ANN_morgan2_is_active.pkl', 'rb'))
    else:
        load_model_5 = pickle.load(open(f'models/{balanced}/{target}/{target}_ANN_opt_balanced_morgan2_is_active.pkl', 'rb'))

    prediction_5 = load_model_5.predict(dataset.mrgn_fp.tolist())
    _ANN_2=  f"{id2name[target]}_ANN_2"                                                       
    predicted_active_5 = pd.DataFrame(prediction_5, columns=[_ANN_2])
    predicted_active_df_5 = predicted_active_df_4.join(predicted_active_5)
    
    return predicted_active_df_5



# Molecular descriptor calculator
def desc_calc(input_data):
    # Performs the descriptor calculation
   input_data["mccs_fp"] = input_data["canonical_smiles"].apply(get_maccs)
   input_data["mrgn_fp"]= input_data["canonical_smiles"].apply(get_mrgn2)
   return input_data


def Predict(uploaded_file,target,balanced="unbalance"):
    # load_data.to_csv('molecule.smi', sep = '\t', header = False, index = False)
    # load_data = uploaded_file
    
    # load_data["canonical_smiles"] = load_data["SMILES"].apply(canocalis)
    # prepared_data = load_data

    dataset = desc_calc(uploaded_file)
    # Apply trained model to make prediction on query compounds
    all_df = build_isactive_model(dataset,target,balanced)
    return all_df


def add_target_columns(combo_df):
    rename = lambda x : x.split('_')[:2] if len(x.split('_')) >=4 else x.split('_')[:1]
    models = {}
    for index in range(len(combo_df)):
        d = combo_df.T
        d = d[index]
        models[index] = ','.join(list(set([rename(x)[0] for x in list(d[d==1].index)])))
        df = pd.DataFrame([models]).T
        df= df.rename(columns={
            0:'Target'
        })

        df['target_number'] = df.Target.apply(lambda x : 0 if ''.join(x.split(',')) =='' else len(x.split(',')))
        final_df = combo_df.join(df)
    return final_df



# File download
def filedownload(df):
    csv = df.to_csv(index=False)
    b64 = base64.b64encode(csv.encode()).decode()  # strings <-> bytes conversions
    href = f'<a href="data:file/csv;base64,{b64}" download="Prediction.csv">Download Predictions</a>'
    return href


st.title("minKLIFSAI Virtual sceening page")


# Page title
st.markdown("""
# Predict kinase inhibitor activity \n
---
""")


sorted_sector_unique = sorted(list(id2name.values()))
selected_sector = st.sidebar.multiselect('Select Target for virtual screen', sorted_sector_unique)


balanced_state = st.sidebar.radio('Choose type of model:',["imbalance","balanced"], index=0)
if balanced_state== "balanced":
    balanced_label = "balance_opt"
else: 
    balanced_label = "unbalance"


if balanced_label== "balance_opt":
   # Logo image
   image = Image.open('balanced.png')
   st.image(image, use_column_width=True)

if balanced_label== "unbalance":
   # Logo image
   image = Image.open('imbalance.png')
   st.image(image, use_column_width=True)

         
# Sidebar
with st.sidebar.header('1. Upload your CSV data'):
    uploaded_file = st.sidebar.file_uploader("Upload your input file", type=['txt'])
    st.sidebar.markdown("""
[Example input file (same order)](./example.txt)
""")

# id2name_= lambda x : name2id[x]

options = st.multiselect('Uniprot to name', sorted_sector_unique)
dic = {name2id[x]:x for x in options}
st.write(dic)


if st.sidebar.button('Predict'):
#     try:
    try:
        file = pd.read_table(uploaded_file, sep=' ', header=0, names= ["SMILES","Name"])
        load_data = file
        load_data["canonical_smiles"] = load_data["SMILES"].apply(canonicalise)
        loaded_data_only = load_data[:100]
    except:
        print("Invalid Input..!")
    prepared_data = loaded_data_only
    df_list =[]
    try:
        targets = [name2id[x] for x in selected_sector]

    except:
        targets = sorted(list(id2name.keys()))
    if not targets:
        targets = sorted(list(id2name.keys()))
    my_bar = st.progress(0)
    contr = 0
    for target in targets:
        df = Predict(prepared_data,target,balanced_label)
        loaded_data_only = loaded_data_only.join(df)

        contr+= int(1/len(targets)*100)

        my_bar.progress(contr)

    combo_df = loaded_data_only.drop(['mccs_fp','mrgn_fp','canonical_smiles'], axis=1)
    st.write(combo_df)
    st.write('Download Prediction data...')
    st.markdown(filedownload(combo_df), unsafe_allow_html=True)
    combo_df = loaded_data_only.drop(['mccs_fp','mrgn_fp','canonical_smiles'], axis=1)
    final_df = add_target_columns(combo_df)
    final_df = final_df[['Target','target_number']]
    st.write(final_df)
    st.write('Download targets in details...')
    st.markdown(filedownload(final_df), unsafe_allow_html=True)        
#     except:
#         st.write("Something went wrong...!")
else:
    st.info('Upload input data in the sidebar to start!')
 

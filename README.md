# minKLIFSAI
minimum KLIFSAI (a machine learning model to predict Kinase activity)

To start the app download the repo and extract the files in the folder write:

`pip install -r requirements.txt`

then run:

`streamlit run 1_virtual_screen.py`

This was the initial step toward providing a model for every kinase as a tool to discover selective kinase inhibitor.

Repo is a streamlit app that can be used to make a virtual screen for chemical data to predict kinase activity of 480 target id according to Uniprot.
Data model trained on the data from PubChem. 1920 prediction was produced, you can customize the prediction that you want by choosing target you want from sidebar menu and upload the data
and click predict. Model are consisting of two part, one trained on imbalanced data which characterize by high specificity. Another balance data
which train on a data that Inactive compound was less or equal to active compound.
Analysis of all model found in paper supplementary files, here:


Notes: This app is not published yet and all copyright reserved to the owner untill paper published, after that it will be for Non-commercial and for academic use, non-profit,
or research who want to provide a cancer treatment for free.

Metric of each target model found in the paper.

Target name will found in the file: `id2name.json` file


Project are open for contribution, if you want to contribute:

1- Develop a model better than provided

2- Rename it like the model you want to replace

3- Add metric that prove the result like AUC or ROC plot or a reference for a paper and data used if possible.

4- Write you licence either want to be open source, MIT, copyrighted whatever. The default licence here is a MIT.


If you used this tool, please cite:

[1] 1. Abd Elaleem M. Src kinase app: valid inhibitor generation and prediction with explanation using predictive model and selfies. ChemRxiv. 2022; doi:10.26434/chemrxiv-2022-090pf-v2

[2] [10.5281/zenodo.13370507](https://zenodo.org/records/13370507).

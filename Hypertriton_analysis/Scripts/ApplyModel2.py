import os

os.environ["MKL_NUM_THREADS"] = "1" 
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ['OMP_NUM_THREADS'] = "1"


import xgboost as xgb
import numpy as np
import json
import pandas as pd

from hipe4ml.tree_handler import TreeHandler
from hipe4ml.model_handler import ModelHandler


def ApplyModel(pt_min, pt_max, dataH, MultipleModels, model_file_list, slice_list):
    """
    Applies the models which have been trained by ML_Hypertriton.py to the data
    ---
    pt_min: minimum pt
    pt_max: maximum pt
    dataH: dataHandler to which the model is applied to
    MultipleModels: whether multiple models are used or not
    """

    # if single model file is used
    if MultipleModels == False:
        model_file = f'../Models/Hypertriton_model_{pt_min}<pt<{pt_max}'
        model_clf = xgb.XGBClassifier(n_jobs = 2)
        model_hdl = ModelHandler(model_clf,  {'n_jobs' : 2})
        model_hdl.load_model_handler(model_file)
        modellist = model_hdl

    # if multiple models are used
    elif MultipleModels == True:
        modellist=[]
        for model in model_file_list:
            model_clf = xgb.XGBClassifier(n_jobs = 2)
            model_hdl = ModelHandler(model_clf,  {'n_jobs' : 2})
            model_hdl.load_model_handler(model)
            modellist.append(model_hdl)
    
        dataH.slice_data_frame("pt", slice_list, delete_original_df = False)



    #apply model
    dataH.apply_model_handler(modellist, output_margin=True,  column_name="model_output")
    
    if MultipleModels:
        df = pd.concat([dataH.get_slice(i) for i in range(len(slice_list))], ignore_index = True)
        dataH.set_data_frame(df)

    return dataH
    BDT_cut = df["model_output"][i]
    print(BDT_cut)
    selected_data = dataH.get_subset(f'model_output>{BDT_cut}')

    before=dataH.get_n_cand()
    after=selected_data.get_n_cand()

    return after/before



def Procedure_for_one_pt_range(pt_min, pt_max, dataH, MultipleModels, OutputFile):


    print(f'{pt_min} < pt < {pt_max}')

    # apply model to datahandler
    dataHnew =  ApplyModel(pt_min, pt_max, dataH, MultipleModels, model_file_list, slice_list)

    # save dataHandler in root file
    dataHnew.write_df_to_root_files(OutputFile, "df")
    print("Dataframe saved to .root file.")



if __name__ == "__main__":
    # pt ranges:
    pt_min = 3
    pt_max = 6

    # Input files
    MC_file = '../Data/SignalTable_B_20g7.root'
    bkg_file = '../Data/DataTable_18LS_pass3.root'
    data_file = '../Data/DataTable_18_pass3_EPangle2.root'

    # Output file
    OutputFile = "Output.root"
    
    # Multiple models?
    MultipleModels = True

    # if MultipleModels:
    model_file_list = [f'../Models/Hypertriton_model_3.0<pt<3.5', f'../Models/Hypertriton_model_3.5<pt<4.0']
    slice_list = [[3,3.5], [3.5,4]]


    dataH = TreeHandler(library = "np")
    dataH.get_handler_from_large_file(
    file_name = data_file, 
    tree_name = "DataTable", 
    preselection = f'{pt_min} < pt < {pt_max} and 1 < ct < 35',
    max_workers = 2, 
    )

    Procedure_for_one_pt_range(pt_min, pt_max, dataH, MultipleModels, model_file_list, slice_list, OutputFile)



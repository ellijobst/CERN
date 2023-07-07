import numpy as np
import pandas as pd
import xgboost as xgb
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from hipe4ml.model_handler import ModelHandler
from hipe4ml.tree_handler import TreeHandler
from hipe4ml import plot_utils
from matplotlib.backends.backend_pdf import PdfPages

def save_output_as_pdf(filename):
    '''
    this functions saves the output plots from the maschine learning to a .pdf file
    ---
    filename: str, Path/name of the file
    '''
    
    p = PdfPages(filename)
      
    # get number of figures
    fig_nums = plt.get_fignums()  
    figs = [plt.figure(n) for n in fig_nums]
      
    # iterating over the numbers in list
    for fig in figs: 
        
        # and saving the files
        fig.savefig(p, format='pdf') 
        plt.close('all')
    # close the object
    p.close()  

def apply_ML_model(model_file, input_file, output_file):
    # create Model
    model_clf = xgb.XGBClassifier()
    model_hdl = ModelHandler(model_clf)

    #load model
    model_hdl.load_model_handler(model_file)

    #load data
    dataH = TreeHandler()
    dataH.get_handler_from_large_file(
        file_name=input_file, 
        tree_name='DataTable', 
        preselection =f'{pt_min} < pt < {pt_max} and 1 < ct < 35',
        )

    #apply model
    dataH.apply_model_handler(model_hdl, False)

    selected_data_hndl = dataH.get_subset('model_output>0.90')

    #plot
    labels_list = ["after selection","before selection"]
    colors_list = ['orangered', 'cornflowerblue']

    plot_utils.plot_distr(
        [selected_data_hndl, dataH], 
        column='m', 
        bins=200, 
        labels=labels_list, 
        colors=colors_list, 
        density=True,
        fill=True, 
        histtype='step', 
        alpha=0.5,
        )

    ax = plt.gca()
    ax.set_xlabel(r'm($^3_\Lambda H- ^3He - \pi^-$) (GeV/$c^2$)')#TODO: check!
    ax.xaxis.set_label_coords(0.9, -0.075)
    ax.set_xlim([2.9479616, 3.1]) #TODO: check!
    # ax.set_yscale('log')

    # save output as pdf
    save_output_as_pdf(output_file)  
    print('output has been saved as:', output_file)

if __name__ == "__main__":

    pt_min = 3
    pt_max = 4

    model_file = f'../Models/Hypertriton_model_{pt_min}<pt<{pt_max}'
    input_file = '../Data/DataTable_18_pass3.root'
    output_file = f"../Output/ML_Hypertriton_output_{pt_min}<pt<{pt_max}_all.pdf"

    apply_ML_model(model_file, input_file, output_file)
#!/usr/bin/python3
import numpy as np
import pandas as pd
import xgboost as xgb
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from hipe4ml.model_handler import ModelHandler
from hipe4ml.tree_handler import TreeHandler
from hipe4ml.analysis_utils import train_test_generator
from hipe4ml import plot_utils
from matplotlib.backends.backend_pdf import PdfPages


#TODO: add logger

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


def ML_Hypertriton(dataH,bkgH,promptH,filename, pt_min, pt_max):
    '''
    this function trains a model on the input training signal and bkg data, and applies it afterwards
    '''
    # define a test/train dataset
    train_test_data = train_test_generator([promptH, bkgH], [1,0], test_size=0.5, random_state=42)

    #select variables for plots
    vars_to_draw = [
        'V0CosPA',
        'pt', 
        'ProngsDCA', 
        'PiProngPvDCAXY',
        'He3ProngPvDCAXY',
        'He3ProngPvDCA', 
        'PiProngPvDCA',
        'NpidClustersHe3',
        'NitsClustersHe3',
        'TPCnSigmaHe3',
        'TPCnSigmaPi',
        'm',
        # 'ct',
    ]

    #plotting
    leg_labels = ['background', 'signal']

    # plot_utils.plot_distr([bkgH, promptH], vars_to_draw, bins=100, labels=leg_labels, \
    #                     log=True, density=True, figsize=(12, 7), alpha=0.3, grid=False)
    
    # plt.subplots_adjust(left=0.06, bottom=0.06, right=0.99, top=0.96, hspace=0.55, wspace=0.55)
    # plot_utils.plot_corr([bkgH, promptH], vars_to_draw, leg_labels)

    # select variables for training, remove those that are very different for bkg and signal events
    features_for_train = vars_to_draw.copy()
    features_for_train.remove('m')
    features_for_train.remove('pt')

    print('Features used for Training:', features_for_train)

    # create Model
    model_clf = xgb.XGBClassifier()
    model_hdl = ModelHandler(model_clf, features_for_train)

    # pt_min = 4
    # pt_max = 5

    # model_file = f'../Models/Hypertriton_model_{pt_min}<pt<{pt_max}'

    # model_hdl.load_model_handler(model_file)

    # define the shape of the model (how deep, how many trees etc.)
    hyper_pars_ranges = {
        'n_estimators': (50, 500),
        'max_depth': (5, 20), 
        'learning_rate': (0.01, 0.3), 
        'gamma': (0.3, 1.1), 
        'min_child_weight': (1,12),
        'subsample':(0.5, 0.9),
        'colsample_bytree':(0.5, 0.9),
        }
    
    # OBACHT: Do not run on alicecerno2 with n_jobs=-1 !!
    model_hdl.optimize_params_optuna(train_test_data, hyper_pars_ranges, cross_val_scoring='roc_auc',\
                                     timeout=120, n_jobs=8, n_trials=100, direction='maximize')

    model_hdl.train_test_model(train_test_data)

    #saving model
    # model_hdl.dump_model_handler(f'../Models/Hypertriton_model_{pt_min}<pt<{pt_max}')
    # print('Model saved as:', f'Hypertriton_model_{pt_min}<pt<{pt_max}')

    print('model has been trained')
    y_pred_train = model_hdl.predict(train_test_data[0], True)
    y_pred_test = model_hdl.predict(train_test_data[2], True)

    # print(y_pred_train)
    # print(y_pred_test)

    model_hdl.dump_model_handler(f'../Models/Hypertriton_model_{pt_min}<pt<{pt_max}')
    print('Model saved as:', f'Hypertriton_model_{pt_min}<pt<{pt_max}')

    

    #compare trainig and testing set, in case the classifier picked up some features from the training dataset
    plt.rcParams["figure.figsize"] = (10, 7)

    ml_out_fig = plot_utils.plot_output_train_test(model_hdl, train_test_data, 100, 
                                                False, leg_labels, True, density=True)
    #compare efficiencies for signal and bkg
    roc_train_test_fig = plot_utils.plot_roc_train_test(train_test_data[3], y_pred_test,
        train_test_data[1], y_pred_train, None, leg_labels)

    #apply model & plot everything
    dataH.apply_model_handler(model_hdl, False)
    selected_data_hndl = dataH.get_subset('model_output>0.90')

    labels_list = ["after selection","before selection"]
    colors_list = ['orangered', 'cornflowerblue']
    plot_utils.plot_distr([selected_data_hndl, dataH], column='m', bins=200, labels=labels_list, colors=colors_list, density=True,fill=True, histtype='step', alpha=0.5)
    ax = plt.gca()
    ax.set_xlabel(r'm($^3_\lambda H- ^3He - \pi^-$) (GeV/$c^2$)')#TODO: check!
    ax.xaxis.set_label_coords(0.9, -0.075)
    ax.set_xlim([2.9479616, 3.1])
    # ax.set_yscale('log')

    # save output as pdf
    save_output_as_pdf(filename)  
    print('output has been saved as:', filename)


if __name__ == "__main__":
    #specify the pt ranges we want to take a look at
    pt = [2,3,4,5,6,9]

    for i in range(len(pt)-1):
        # specify which dataset to look at
        pt_min = pt[i]
        pt_max = pt[i+1]
        print(f'{pt_min} < pt < {pt_max}')

        # get files
        dataH = TreeHandler()
        promptH = TreeHandler()
        bkgH = TreeHandler()

        # Since Hyptertriton decays into two differently charged particles,
        # when two Like Sign particles are measured it has to be background signal
        dataH.get_handler_from_large_file(file_name='../Data/DataTable_18_pass3.root',tree_name='DataTable', 
                            preselection =f'{pt_min} < pt < {pt_max} and 1 < ct < 35')
        dataH = dataH.get_subset(size=170001)

        promptH.get_handler_from_large_file(file_name='../Data/SignalTable_20g7.root',tree_name='SignalTable', 
                            preselection =f'{pt_min} < pt < {pt_max} and 1 < ct < 35')
        promptH = promptH.get_subset(size=5667)
        
        bkgH.get_handler_from_large_file(file_name='../Data/DataTable_18LS_pass3.root',tree_name='DataTable', 
                            preselection =f'{pt_min} < pt < {pt_max} and 1 < ct < 35')
        
        # df = dataH.get_data_frame()
        # print(df[:10])
        # print(np.amin(df['m']), np.amax(df['m']))
        # break

        # alternatively: select bkg events via the invariant mass
        # bkgH = dataH.get_subset('m < 2.991-3*0.0017 or m > 2.991+3*0.0017', size=promptH.get_n_cand()*3)
        # limits are taken from 
        # https://www.researchgate.net/publication/334719549_Highlights_of_the_production_of_anti-hyper-nuclei_and_exotica_with_ALICE_at_the_LHC/download
        # bkgH = bkgH.get_subset(size=170001)

        filename=f"../Output/ML_Hypertriton_output_{pt_min}_pt_{pt_max}.pdf"  
        
        ML_Hypertriton(dataH, bkgH, promptH, filename, pt_min, pt_max)
        
    
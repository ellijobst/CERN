import os

os.environ["MKL_NUM_THREADS"] = "1" 
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ['OMP_NUM_THREADS'] = "1"


import xgboost as xgb
import numpy as np
import json

from apply_model import get_efficiency, get_significance, apply_ML_model
from HypertritonRestframeBoost import get_df
from hipe4ml.tree_handler import TreeHandler
from hipe4ml.model_handler import ModelHandler
from hipe4ml.analysis_utils import train_test_generator


def TrainModel(bkgH, promptH, pt_min, pt_max, optuna=False):
    # define a test/train dataset
    train_test_data = train_test_generator([promptH, bkgH], [1,0], test_size=0.5, random_state=42)

    # select variables for training, remove those that are very different for bkg and signal events
    features_for_train = [
        'V0CosPA',
        'ProngsDCA', 
        'PiProngPvDCAXY',
        'He3ProngPvDCAXY',
        'He3ProngPvDCA', 
        'PiProngPvDCA',
        'NpidClustersHe3',
        'NitsClustersHe3',
        'TPCnSigmaHe3',
        'TPCnSigmaPi',
    ]

    print('Features used for Training:', features_for_train)

    # create Model
    model_clf = xgb.XGBClassifier(n_jobs=2) #for some reason this is not working
    model_hdl = ModelHandler(model_clf, features_for_train, {'n_jobs' : 2})

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
    
    # ATTENTION: Do not run on alicecerno2 with n_jobs=-1 !!
    if optuna == True:
        model_hdl.optimize_params_optuna(train_test_data, hyper_pars_ranges, cross_val_scoring='roc_auc',\
                                     timeout=120, n_jobs=2, n_trials=100, direction='maximize')

    model_hdl.train_test_model(train_test_data, {'n_jobs' : 2})

    print('model has been trained')

    # save model
    model_hdl.dump_model_handler(f'../Models/Hypertriton_model_{pt_min}<pt<{pt_max}')
    print('Model saved as:', f'Hypertriton_model_{pt_min}<pt<{pt_max}')


def OptimizeBDT(pt_min, pt_max, MC_file, bkg_file):
    # load model and output file
    model_file = f'../Models/Hypertriton_model_{pt_min}<pt<{pt_max}'

    # BDT cuts
    BDT_cuts = np.linspace(-15,10,100)

    # calculate efficiency for different BDT values
    efficiency = get_efficiency(MC_file, model_file, BDT_cuts, pt_min, pt_max)
    with open(f"../Output/SignEff/efficiency_{pt_min}_pt_{pt_max}.json", 'w') as f:
        json.dump(efficiency, f, indent=2) 

    print('efficiency saved as:', f"../Output/SignEff/efficiency_{pt_min}_pt_{pt_max}.json")

    # calculate significance for different BDT values
    # TODO: check if this is really the MC or the real data that should be used here? (Adjust: tree_name, data_file)
    significance = get_significance(efficiency, BDT_cuts, data_file=MC_file, bkg_file=bkg_file, tree_name='SignalTable', pt_cut=(pt_min, pt_max), model_file=model_file)
    
    with open(f"../Output/SignEff/significance_{pt_min}_pt_{pt_max}.json", 'w') as f:
        json.dump(significance, f, indent=2) 

    print('significance saved as:', f"significance_{pt_min}_pt_{pt_max}.json")

    # find best BDT cut
    signxeff = [x*y for x,y in zip(significance, efficiency)]

    with open(f"../Output/SignEff/signxeff_{pt_min}_pt_{pt_max}.json", 'w') as f:
        json.dump(signxeff, f, indent=2) 
    
    print('SignificancexEfficiency saved as:', f"../Output/SignEff/signxeff_{pt_min}_pt_{pt_max}.json")
    
    Best_BDT = BDT_cuts[signxeff.index(max(signxeff))]
    Best_BDT_Sign = significance[signxeff.index(max(signxeff))]
    print("Best BDT Cut: ", Best_BDT)

    return Best_BDT, Best_BDT_Sign

# not used atm
def ApplyModel(pt_min, pt_max, Best_BDT):

    MC = "MC"
    data_file = '../Data/SignalTable_B_20g7.root'
    data_tree = 'SignalTable'

    # load model
    model_file = f'../Models/Hypertriton_model_{pt_min}<pt<{pt_max}'

    # apply model to data (cut on ct is included in this function)
    selected_data_hndl, dataH =  apply_ML_model(model_file, data_file, tree_name=data_tree, output_file=None, BDT_cut=Best_BDT, pt_cut=(pt_min, pt_max))

    before=dataH.get_n_cand()
    after=selected_data_hndl.get_n_cand()

    print("BDT CUT EFFICIENCY:", after/before)

    # apply hypertriton restframe boost to get cos(theta*) columns
    df_with_cos_theta = get_df(selected_data_hndl)

    # save df as .root fule
    selected_data_hndl.set_data_frame(df_with_cos_theta)
    selected_data_hndl.write_df_to_root_files(f"SelectedDataFrame{MC}_{pt_min}_pt_{pt_max}", "df")
    print("df saved to .root file.")


def Procedure_for_one_pt_range(pt_min, pt_max):
    print(f'{pt_min} < pt < {pt_max}')
    #-----------------------LOAD TREES----------------------
    dataH = TreeHandler(library = "np")
    promptH = TreeHandler(library = "np")
    bkgH = TreeHandler(library = "np")

    data_file = '../Data/DataTable_18_pass3.root'
    MC_file = f'../Data/SignalTable_20g7.root'
    bkg_file = f'../Data/DataTable_18LS_pass3.root'
    model_file = f'../Models/Hypertriton_model_{pt_min}<pt<{pt_max}'

    dataH.get_handler_from_large_file(
        file_name = data_file,
        tree_name = 'DataTable', 
        preselection = f'{pt_min} < pt < {pt_max} and 1 < ct < 35',
        max_workers = 2,
        )
    dataH = dataH.get_subset(size=170001)

    promptH.get_handler_from_large_file(
        file_name = MC_file,
        tree_name ='SignalTable', 
        preselection = f'{pt_min} < pt < {pt_max} and 1 < ct < 35',
        max_workers = 2,
        )
    promptH = promptH.get_subset(size=5667)
    
    bkgH.get_handler_from_large_file(
        file_name = bkg_file,
        tree_name = 'DataTable', 
        preselection = f'{pt_min} < pt < {pt_max} and 1 < ct < 35',
        max_workers = 2,
        )
    bkgH = bkgH.get_subset(size=170001)

    #----------------------TRAIN MODEL----------------------
    # TODO: set optuna to true when running over weekend
    TrainModel(bkgH, promptH, pt_min, pt_max, optuna=False)
    
    #------------------BDT OPTIMIZATION---------------------
    BDTOpt = OptimizeBDT(pt_min, pt_max, MC_file, bkg_file)
    Best_BDT = BDTOpt[0]
    Significance = BDTOpt[1]

    #--------------APPLY MODEL WITH BDT CUT-----------------
    # ApplyModel(pt_min, pt_max, Best_BDT)

    return Significance


if __name__ == "__main__":
    # pt ranges:
    max_int_size = 7#GeV/c
    min_int_size = 2#GeV/c #before: 1.6
    step_size = 1#GeV/c
    maximumPt=9#GeV/c
    minimumPt=2#GeV/c

    pt_ranges=[]
    significances=[]

    for i in range(int((max_int_size-min_int_size)/step_size)+1):
        int_size = max_int_size - i*step_size
        print("pt interval size", int_size)
        # if i < 3:
        #     continue
        for j in range(i+1):
            pt_min = minimumPt+step_size*j
            pt_max = pt_min + int_size

            Significance = Procedure_for_one_pt_range(pt_min, pt_max)

            pt_ranges.append((pt_min, pt_max))
            significances.append(Significance)

    with open(f"../Output/ptSignificances{step_size}GeV.json", 'w') as f:
        json.dump(significances, f, indent=2) 

    with open(f"../Output/ptIntervals{step_size}GeV.json", 'w') as f:
        json.dump(pt_ranges, f, indent=2) 

    Best_pT = pt_ranges[significances.index(max(significances))]
    print("best pt interval:", Best_pT)

            

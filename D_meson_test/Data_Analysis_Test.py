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


#MC data
promptH = TreeHandler('prompt.root','treeMLDplus')

#real data
dataH = TreeHandler('data.root','treeMLDplus')

#background, cut on invariant mass, at least 3 times the counts for background as we have for the signal are needed
bkgH = dataH.get_subset('inv_mass < 1.82 or 1.92 < inv_mass < 2.00', size=promptH.get_n_cand()*3)
print('number of candidates:', promptH.get_n_cand())
# define a test/train dataset
train_test_data = train_test_generator([promptH, bkgH], [1,0], test_size=0.5, random_state=42)

#select variables
vars_to_draw = promptH.get_var_names()


leg_labels = ['background', 'signal']

#plotting
plot_utils.plot_distr([bkgH, promptH], vars_to_draw, bins=100, labels=leg_labels, log=True, density=True, figsize=(12, 7), alpha=0.3, grid=False)
plt.subplots_adjust(left=0.06, bottom=0.06, right=0.99, top=0.96, hspace=0.55, wspace=0.55)
plot_utils.plot_corr([bkgH, promptH], vars_to_draw, leg_labels)

# plt.show()

features_for_train = vars_to_draw.copy()

#remove invariant mass, since it is known, and transverse momentum, since it can differ in the MCs from the real data
features_for_train.remove('inv_mass')
features_for_train.remove('pt_cand')
print('features for training:', features_for_train)

model_clf = xgb.XGBClassifier()
model_hdl = ModelHandler(model_clf, features_for_train)

# define the shape of the model (how deep, how many trees etc.)
hyper_pars_ranges = {'n_estimators': (200, 1000), 'max_depth': (
    2, 4), 'learning_rate': (0.01, 0.1)}
model_hdl.optimize_params_optuna(train_test_data, hyper_pars_ranges, cross_val_scoring='roc_auc', timeout=120,
                                 n_jobs=-1, n_trials=100, direction='maximize')

model_hdl.train_test_model(train_test_data)

y_pred_train = model_hdl.predict(train_test_data[0], False)
y_pred_test = model_hdl.predict(train_test_data[2], False)

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
plot_utils.plot_distr([selected_data_hndl, dataH], column='inv_mass', bins=200, labels=labels_list, colors=colors_list, density=True,fill=True, histtype='step', alpha=0.5)
ax = plt.gca()
ax.set_xlabel(r'm(K$^-\pi^+\pi^+$) (GeV/$c^2$)')
ax.margins(x=0)
ax.xaxis.set_label_coords(0.9, -0.075)

# plt.show()

def save_image(filename):
    
    # PdfPages is a wrapper around pdf 
    # file so there is no clash and
    # create files with no error.
    p = PdfPages(filename)
      
    # get_fignums Return list of existing
    # figure numbers
    fig_nums = plt.get_fignums()  
    figs = [plt.figure(n) for n in fig_nums]
      
    # iterating over the numbers in list
    for fig in figs: 
        
        # and saving the files
        fig.savefig(p, format='pdf') 
          
    # close the object
    p.close()  

# name your Pdf file
filename = "Data_Analysis_test_ouput.pdf"  
  
# call the function
save_image(filename)  
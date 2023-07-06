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


pt_min = 2
pt_max = 3

modelH = load_model_handler(f'Models/Hypertriton_model_{pt_min}<pt<{pt_max}')



#apply model & plot everything
dataH = TreeHandler('DataTable_18_pass3.root','DataTable')
dataH.apply_model_handler(model_hdl, False)

selected_data_hndl = dataH.get_subset('model_output>0.90')

labels_list = ["after selection","before selection"]
colors_list = ['orangered', 'cornflowerblue']

# TODO: make this work!
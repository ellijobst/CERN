import numpy as np
import pandas as pd
import xgboost as xgb
import matplotlib.pyplot as plt
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

def apply_ML_model(model_file, data_file, tree_name, output_file, BDT_cut, pt_cut=(3,4)):
    '''
    This function applies a Model to dataH.
    ---
    model_file: file in which the model is saved
    data_file: file to which to model should be applied
    output_file: path to which the output plots should be saved as .pdf file
    BDT_cut: the BDT_cut applied to the data
    ---
    returns: [dataH after the BDT cut, dataH]
    '''
    # set limits for pt
    pt_min = pt_cut[0]
    pt_max = pt_cut[1]

    # create Model
    model_clf = xgb.XGBClassifier()
    model_hdl = ModelHandler(model_clf)

    #load model
    model_hdl.load_model_handler(model_file)

    #load data and cut on pt
    dataH = TreeHandler()
    dataH.get_handler_from_large_file(
        file_name=data_file, 
        tree_name=tree_name, 
        preselection =f'{pt_min} < pt < {pt_max} and 1 < ct < 35',
        )

    #apply model
    dataH.apply_model_handler(model_hdl, False)

    # cut on BDT
    selected_data_hndl = dataH.get_subset(f'model_output>{BDT_cut}')

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
    plt.close('all')
    # save output as pdf
    if output_file!=None:
        save_output_as_pdf(output_file)  
        print('output has been saved as:', output_file)

    return [selected_data_hndl, dataH]

def calculate_efficiency(model_file, MC_file, BDT_cut, generatedH, recoH, output):
    # apply model+BDTcut to MC data
    recoBDTH = apply_ML_model(model_file, MC_file, tree_name='SignalTable', output_file=None, BDT_cut=BDT_cut)[0]

    df_recoBDTH = recoBDTH.get_data_frame()
    df_generatedH = generatedH.get_data_frame()

    if output != False:
        print('Len Reco + HDF:', len(df_recoBDTH))
        print('# generated Hypertritons:', len(df_generatedH))
        plot_utils.plot_distr(
            [generatedH, recoH, recoBDTH],
            column='pt', 
            bins=200, 
            labels=['generated Hypertritons', 'reconstructed Hypertritons', 'reco Hypertritons after BDT_cut'],
            colors = ['orangered', 'cornflowerblue', 'limegreen'],
            density=False,#normalize to total number of counts
            fill=True, 
            histtype='step', 
            alpha=0.5,
            )
        save_output_as_pdf(f'../Output/pt_disitribution_gen_reco_BDT>{BDT_cut}_{pt_min}<pt<{pt_max}.pdf')
        print('Output saved as:',f'../Output/pt_disitribution_gen_reco_BDT>{BDT_cut}_{pt_min}<pt<{pt_max}.pdf')


    # number of events in the corresponding pt range
    n_gen = len(df_generatedH)
    n_reco = len(df_recoBDTH)

    # calculate efficiency
    eff = n_reco/n_gen #TODO: check! Es könnte auch sein, dass man durch die recoH teilen muss, dann würde es mehr sinn ergeben
    print(BDT_cut, eff)
    return eff

def get_efficiency(MC_file, model_file):
    '''
    This function calculates the efficiency depending on the BDT cut. For the efficiency we only use MC files
    TODO: there is a whole function for this!
    ---
    MC_file: input Monte Carlo file, with generated and reconstructed Hypertritons
    model_file: model that is applied to the reconstructed Hypertritons
    BDT_cut: BDT_cut that is applied
    '''
    # number of expected hypertritons = 2*278971416*2.6e-5*0.25*efficiency
    # number of signal reco is the number of entries in the 3sigma region, apply pt cut and BDT cut
    # number of bkg is the number of bkg entries in the 3sigma region
    # to find the region fit a gaussian to the peak, for bgk we use the like sign file with a cut on the mass
    # get significance S = nbkg/sqrt(nbkg+nsignal)
    # plot significance over bdt cut
    # plot efficiency*significance over bdt cut -> find maximum
    
    # load MC events

    # generated Hypertritons
    generatedH = TreeHandler()
    generatedH.get_handler_from_large_file(
        file_name=MC_file, 
        tree_name='GenTable', #TODO: check!
        preselection =f'{pt_min} < pt < {pt_max}', #only apply pt cut
        )

    # reconstructed Hypertritons 
    recoH = TreeHandler()
    recoH.get_handler_from_large_file(
        file_name=MC_file, 
        tree_name='SignalTable', 
        preselection =f'{pt_min} < pt < {pt_max} and 1 < ct < 35', #apply all cuts
        )
    
    # apply model to reconstructed Hypertritons
    # pt cut is already considered here!


    BDT_cuts = np.linspace(0,0.99,100)

    efficiency = [calculate_efficiency(MC_file=MC_file, model_file=model_file, BDT_cut=BDT_cut, generatedH=generatedH, recoH=recoH, output=False) for BDT_cut in BDT_cuts]
    
    return [BDT_cuts, efficiency]
    
def get_significance():
    pass
    #TODO
    # N_hyp = 2*278971416*2.6e-5*0.25*eff
    

if __name__ == "__main__":

    pt_min = 3
    pt_max = 4

    model_file = f'../Models/Hypertriton_model_{pt_min}<pt<{pt_max}'
    data_file = '../Data/DataTable_18_pass3.root'
    output_file = f"../Output/ML_Hypertriton_output_{pt_min}<pt<{pt_max}_all.pdf"
    MC_file = f'../Data/SignalTable_20g7.root'

    # calculate efficiency for different BDT Cut values in the range of 0 to 1 in steps of 0.01
    BDT_cuts, efficiency = get_efficiency(MC_file, model_file)
    # plot efficency over BDT cut
    # TODO: find better way(for sure there is a hipe4ml function!)
    plt.plot(
        BDT_cuts, efficiency,
        # label=r'$\varepsilon = Number of ^3_\Lambda H_{reco} / Number of ^3_\Lambda H_{gen}$',
        color='orangered',
    )
    plt.xlabel('BDT cut')
    plt.ylabel(r'Efficiency $\varepsilon $')
    # plt.legend()

    save_output_as_pdf(f'../Output/efficiency_BDT_{pt_min}<pt<{pt_max}.pdf')
    plt.show()
    #TODO: missing other steps

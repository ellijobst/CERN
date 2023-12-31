#!/usr/bin/python3
import numpy as np
import pandas as pd
import xgboost as xgb
import matplotlib.pyplot as plt
import json
import scipy.optimize as opt
import scipy.stats as stats
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

    fig_nums = plt.get_fignums()  
    figs = [plt.figure(n) for n in fig_nums]

    for fig in figs: 
        fig.savefig(p, format='pdf') 
        # plt.close('all')
    # close the object
    p.close()  

def norm(m, loc, scale, A):
    return  A*stats.norm.pdf(m, loc=loc, scale=scale)

def StartParamsNorm(df):
    count, bins = np.histogram(df['m'], bins=50)
    centered_bins = (bins[:-1] + bins[1:]) / 2
    # count, bins = np.histogram(df, bins=50)
    count_max = np.amax(count)

    #calculate the first 2 moments from the ditribution
    m = np.mean(np.array(df))
    s = np.std(np.array(df))

    #calculate params for the distribution
    loc = m
    scale = s

    #calculate amplitude
    dist = np.array(stats.norm.pdf(x=centered_bins,loc=loc,scale=scale))
    dist_max = np.amax(dist)

    A = count_max/dist_max
    return [loc, scale, A]

def apply_ML_model(model_file, data_file, tree_name, output_file, BDT_cut, pt_cut):
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
    model_clf = xgb.XGBClassifier(n_jobs = 2)
    model_hdl = ModelHandler(model_clf,  {'n_jobs' : 2})

    #load model
    model_hdl.load_model_handler(model_file)

    #load data and cut on pt
    dataH = TreeHandler(library = "np")
    dataH.get_handler_from_large_file(
        file_name = data_file, 
        tree_name = tree_name, 
        preselection = f'{pt_min} < pt < {pt_max} and 1 < ct < 35',
        max_workers = 2, 
        )

    #apply model
    dataH.apply_model_handler(model_hdl, output_margin=True)# note: this needs to be set to true!!

    # cut on BDT
    selected_data_hndl = dataH.get_subset(f'model_output>{BDT_cut}')

    # without BDT cut
    # selected_data_hndl = dataH

    #plot BDT Value distribution to figure out the BDT Value bin limits
    '''
    plt.figure()
    labels_list = ['BDT value distribution']
    colors_list = ['orangered']

    plot_utils.plot_distr(
        [dataH], 
        column='model_output', 
        bins=200, 
        labels=labels_list, 
        colors=colors_list, 
        density=True,
        fill=True, 
        histtype='step', 
        alpha=0.5,
        )

    ax = plt.gca()
    plt.show()
    # ax.set_xlabel(r'm($^3$He $\pi^-$) (GeV/$c^2$)')#TODO: check!
    # ax.xaxis.set_label_coords(0.9, -0.075)
    # ax.set_xlim([2.9479616, 3.1]) #TODO: check!
    # ax.set_yscale('log')
    plt.close()
    '''

    # save output as pdf
    if output_file!=None:
        save_output_as_pdf(output_file)  
        print('output has been saved as:', output_file)

    return [selected_data_hndl, dataH]

def calculate_efficiency(model_file, MC_file, BDT_cut, generatedH, recoH, output, pt_min, pt_max):
    # apply model+BDTcut to MC data
    recoBDTH = apply_ML_model(model_file, MC_file, tree_name='SignalTable', output_file=None, BDT_cut=BDT_cut, pt_cut=(pt_min, pt_max))[0]

    n_reco = recoBDTH.get_n_cand() 
    n_gen = generatedH.get_n_cand()

    # note: BDT absolute value is saved as 'model_output' column

    # plot distribution
    if output == True:
        counts = plot_utils.plot_distr(
            [recoBDTH],
            column='model_output', 
            bins=10, 
            labels=['reco Hypertritons'],
            colors = ['orangered'],
            density=False,#normalize to total number of counts
            fill=True, 
            histtype='step', 
            alpha=0.5,
            )
        
        # save_output_as_pdf(f'../Output/pt_disitribution_gen_reco_BDT>{BDT_cut}_{pt_min}<pt<{pt_max}.pdf')
        plt.show()
        print(counts)
        # print('Output saved as:',f'../Output/pt_disitribution_gen_reco_BDT>{BDT_cut}_{pt_min}<pt<{pt_max}.pdf')

    # calculate efficiency

    eff = n_reco/n_gen 
    print(BDT_cut, eff)
    return eff

def get_efficiency(MC_file, model_file, BDT_cuts, pt_min, pt_max):
    '''
    This function calculates the efficiency depending on the BDT cut. For the efficiency we only use MC files
    TODO: there is a whole function for this!
    ---
    MC_file: input Monte Carlo file, with generated and reconstructed Hypertritons
    model_file: model that is applied to the reconstructed Hypertritons
    BDT_cut: BDT_cut that is applied
    '''

    # generated Hypertritons
    generatedH = TreeHandler(library = "np")
    generatedH.get_handler_from_large_file(
        file_name=MC_file, 
        tree_name='GenTable', 
        preselection =f'{pt_min} < pt < {pt_max}', #only apply pt cut
        max_workers = 2,
        )

    # reconstructed Hypertritons 
    recoH = TreeHandler(library = "np")
    recoH.get_handler_from_large_file(
        file_name=MC_file, 
        tree_name='SignalTable', 
        preselection =f'{pt_min} < pt < {pt_max} and 1 < ct < 35', #apply all cuts
        max_workers = 2,
        )
    
    print(recoH.get_n_cand(), generatedH.get_n_cand())
    
    # apply model to reconstructed Hypertritons
    # pt cut is already considered here!

    efficiency = [calculate_efficiency(MC_file=MC_file, model_file=model_file, BDT_cut=BDT_cut, generatedH=generatedH, recoH=recoH, output=False, pt_min = pt_min, pt_max=pt_max) for BDT_cut in BDT_cuts]
    
    return efficiency
    
def get_significance(efficiency, BDT_cuts, data_file, bkg_file, tree_name, pt_cut, model_file):
    # pt_min = pt_cut[0]
    # pt_max = pt_cut[1]
   
    # get number of hypertritons
    N_hyp = [2*278971416*2.6e-5*0.25*eff for eff in efficiency]
    
    # plt Number of Hypertritons
    '''
    plt.plot(
        BDT_cuts, N_hyp,
        color='orangered',
    )
    plt.xlabel('BDT cut')
    plt.ylabel(r'$N_{Hyp}$')
    # plt.legend()
    plt.show()
    save_output_as_pdf(f'../Output/Nhyp_BDT_{pt_min}<pt<{pt_max}.pdf')
    print('Output saved as:', f'../Output/Nhyp_BDT_{pt_min}<pt<{pt_max}.pdf')
    plt.close()
    '''

    # determine cut on invariant mass with gaussian fit
    selected_data_hndl, dataH =  apply_ML_model(model_file, data_file, tree_name, output_file=None, BDT_cut=0.7, pt_cut=pt_cut)
   
    df = selected_data_hndl.get_data_frame()
    count, bins = np.histogram(df['m'], bins=50)
    centered_bins = (bins[:-1] + bins[1:]) / 2

    m_Hyp = 2.991
    p0=[m_Hyp, 3*0.0017, 150]
    # p0=StartParamsNorm(df)
    optimizedParameters, pcov = opt.curve_fit(norm, centered_bins, count, p0=p0) #TODO:vlt zu root fit ändern

    # plt.figure()
    # plot_utils.plot_distr(selected_data_hndl, column=['m'], bins=50)
    # plt.plot(np.linspace(bins[0], bins[-1], 100), norm(np.linspace(bins[0], bins[-1], 100), *optimizedParameters), label='Gaussian Fit')
    # plt.show()
    # plt.close()

    m_Hyp_fit = optimizedParameters[0]
    std = optimizedParameters[1]

    m_min = m_Hyp_fit-3*std
    m_max = m_Hyp_fit+3*std

    # get number of background events

    N_bkg = [calculate_number_of_background_events(model_file, 
                                                   data_file, 
                                                   bkg_file, 
                                                   tree_name, 
                                                   output_file=None, 
                                                   BDT_cut=BDT_cut, 
                                                   m_cut = (m_min, m_max),
                                                   pt_cut=pt_cut) for BDT_cut in BDT_cuts]
    # N_bkg = 119166
    significance = [n_hyp/np.sqrt(n_hyp+n_bkg) for n_hyp,n_bkg in zip(N_hyp,N_bkg)]
    
    # plot significance
    '''
    plt.figure()
    plt.plot(BDT_cuts, significance, color='orangered')
    plt.xlabel('BDT cut')
    plt.ylabel(r'Significance')
    plt.show()
    save_output_as_pdf(f'../Output/Sign_BDT_{pt_min}<pt<{pt_max}.pdf')
    print('Output saved as:', f'../Output/Sign_BDT_{pt_min}<pt<{pt_max}.pdf')
    plt.close()
    '''
    return significance

def calculate_number_of_background_events(model_file, data_file, bkg_file, tree_name, output_file, m_cut, BDT_cut, pt_cut):
    # set limits for preselection
    pt_min = pt_cut[0]
    pt_max = pt_cut[1]

    m_min = m_cut[0]
    m_max = m_cut[1]

    # load model
    model_clf = xgb.XGBClassifier(n_jobs = 2)
    model_hdl = ModelHandler(model_clf, {"n_jobs" : 2})
    model_hdl.load_model_handler(model_file)

    # load bkg file and apply model
    bkgH = TreeHandler(library = "np")
    bkgH.get_handler_from_large_file(
        file_name=bkg_file,
        tree_name='DataTable', 
        preselection =f'{pt_min} < pt < {pt_max} and 1 < ct < 35 and {m_min} < m < {m_max}',
        max_workers = 2,
    )
    bkgH.apply_model_handler(model_hdl, output_margin=True)
    
    # N_before = bkgH.get_n_cand()

    # cut on BDT
    selected_bkg_hndl = bkgH.get_subset(f'model_output>{BDT_cut}')
    N_bkg = selected_bkg_hndl.get_n_cand()
    # N_bkg_after = N_before - N_bkg
    print(N_bkg)
    return N_bkg


if __name__ == "__main__":
    
    # input data files
    data_file = '../Data/DataTable_18_pass3.root'#NOTE: This file is outdated!!
    MC_file = f'../Data/SignalTable_20g7.root'
    bkg_file = f'../Data/DataTable_18LS_pass3.root'
    
    pt = [4,5]

    for i in range(len(pt)-1):
        # specify which dataset to look at

        print("Pt_min:")
        pt_min = input()
        print("pt_max:")
        pt_max = input()



        print(f'{pt_min} < pt < {pt_max}')

        # model and output file
        model_file = f'../Models/Hypertriton_model_{pt_min}<pt<{pt_max}'
        output_file = f"../Output/ML_Hypertriton_output_{pt_min}<pt<{pt_max}_all.pdf"

   
        # apply_ML_model(model_file, data_file, tree_name='DataTable', output_file='generic_output_file.pdf', pt_cut=(pt_min, pt_max), BDT_cut=0.8)

        # TODO: put this part in a function
        # TODO: add single scripts for efficiency, significance. etc.


        BDT_cuts = np.linspace(-15,10,100)

        # calculate efficiency for different BDT values

        # efficiency = get_efficiency(MC_file, model_file, BDT_cuts, pt_min, pt_max)
        # with open(f"../Output/SignEff/efficiency_{pt_min}_pt_{pt_max}.json", 'w') as f:
        #     json.dump(efficiency, f, indent=2) 

        # print('efficiency saved as:', f"efficiency_{pt_min}_pt_{pt_max}.json")

        with open(f"../Output/SignEff/efficiency_{pt_min}_pt_{pt_max}.json", 'r') as f:
            efficiency = json.load(f)


        # calculate significance for different BDT values

        # significance = get_significance(efficiency, BDT_cuts, data_file, bkg_file, tree_name='DataTable', pt_cut=(pt_min, pt_max), model_file =model_file)
        
        # with open(f"../Output/SignEff/significance_{pt_min}_pt_{pt_max}.json", 'w') as f:
        #     json.dump(significance, f, indent=2) 

        print('significance saved as:', f"significance_{pt_min}_pt_{pt_max}.json")

        with open(f"../Output/SignEff/significance_{pt_min}_pt_{pt_max}.json", 'r') as f:
            significance = json.load(f)
        
        
        # plot efficency over BDT cut
        # TODO: find better way(for sure there is a hipe4ml function!)  
        

        plt.figure()

        signxeff = [x*y for x,y in zip(significance, efficiency)]

        with open(f"../Output/SignEff/signxeff_{pt_min}_pt_{pt_max}.json", 'w') as f:
            json.dump(signxeff, f, indent=2) 
        
        print('SignificancexEfficiency saved as:', f"../Output/SignEff/signxeff_{pt_min}_pt_{pt_max}.json")
        
        Best_BDT = BDT_cuts[signxeff.index(max(signxeff))]
        print(Best_BDT)
        # hier ist es noch der richtige!!

        plt.figure()
        plt.plot(
            BDT_cuts, efficiency, 
            label=f"BestBDTEff = {efficiency[signxeff.index(max(signxeff))]}",
            color='orangered',
        )
        plt.xlabel('BDT cut')
        plt.axvline(Best_BDT)
        plt.legend()
        plt.ylabel(r'Efficiency $\varepsilon $')
        print('Efficiency saved.')


        plt.figure()
        plt.plot(
            BDT_cuts, significance,
            color='orangered',
        )
        plt.xlabel('BDT cut')
        plt.ylabel(r'Significance')
        print('Significance saved.')

        
        plt.plot(
            BDT_cuts, signxeff,
            color='orangered',
            label=f'Maximum at {Best_BDT}.'
        )
        plt.axvline(Best_BDT)
        plt.xlabel('BDT cut')
        plt.ylabel(r'Significance x Efficiency')
        plt.legend()


        save_output_as_pdf(f'../Output/SignEff/SignEff_output_{pt_min}_pt_{pt_max}.pdf')
        plt.close('all')

        break
        print('Done!')
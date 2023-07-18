#!/usr/bin/python3
import numpy as np
import pandas as pd
import xgboost as xgb
import matplotlib.pyplot as plt
from hipe4ml.model_handler import ModelHandler
from hipe4ml.tree_handler import TreeHandler
from hipe4ml import plot_utils
from matplotlib.backends.backend_pdf import PdfPages



def get_gamma(m,px,py,pz):
    E = get_energy(m,px,py,pz)
    gamma = 1/np.sqrt(1-(px**2 + py**2 + pz**2)/E**2)
    return gamma

#pz is along beam 

#calculate Energy
def get_energy(m,px,py,pz):
    return np.sqrt(m**2+px**2+py**2+pz**2)

def get_pT(px,py):
    return np.sqrt(px**2+py**2)

#cuts 
#pseudorapidity < 0.8
def get_peudorapidity(px,py,pz):
    p = np.sqrt(px**2+py**2+pz**2)
    eta = 1/2* np.log((p+pz)/(p-pz))
    return eta

def get_m_Hyp(px_He,py_He,pz_He, px_Pi,py_Pi,pz_Pi, px_Hyp, py_Hyp,pz_Hyp):
    # get energy of daughters
    E_He = get_energy(m_He, px_He, py_He, pz_He)
    E_Pi = get_energy(m_Pi, px_Pi, py_Pi, pz_Pi)
    E_Hyp = E_Pi + E_He

    p_Hyp2 = px_Hyp**2+py_Hyp**2+pz_Hyp**2

    m_Hyp = np.sqrt(E_Hyp**2-p_Hyp2)

    return m_Hyp



def get_CMSFourVector(m,px,py,pz,m_Hyp, px_Hyp, py_Hyp, pz_Hyp):
    
    # calculate gamma and beta factors of Hypertriton
    E_Hyp = get_energy(m_Hyp,px_Hyp,py_Hyp,pz_Hyp)
    gamma = get_gamma(m_Hyp,px_Hyp,py_Hyp,pz_Hyp)
    betax = px_Hyp/E_Hyp
    betay = py_Hyp/E_Hyp
    betaz = pz_Hyp/E_Hyp
    beta = np.sqrt(1-1/gamma**2)

    #calculate Energy of daughter (momenta are already in file)
    E = get_energy(m, px, py, pz)
   
    # transform into the rest frame of the Hypertriton
    E_CMS = gamma*(E-(betax*px+betay*py+betaz*pz))
    px_CMS = -gamma*betax*E+px+(gamma-1)/beta**2 *(betax**2*px+betax*betay*py+betax*betaz*pz)
    py_CMS = -gamma*betay*E+py+(gamma-1)/beta**2 *(betax*betay*px+betay**2*py+betay*betaz*pz)
    pz_CMS = -gamma*betaz*E+pz+(gamma-1)/beta**2 *(betax*betaz*px+betay*betaz*py+betaz*betaz*pz)

    return [E_CMS, px_CMS, py_CMS, pz_CMS]

def boost(df, i):
    px_Pi = df['pxPi'][i]
    py_Pi = df['pyPi'][i]
    pz_Pi = df['pzPi'][i]

    px_He = df['pxHe3'][i]
    py_He = df['pyHe3'][i]
    pz_He = df['pzHe3'][i]

    px_Hyp = px_Pi+px_He
    py_Hyp = py_Pi+py_He
    pz_Hyp = pz_Pi+pz_He

    m_Hyp = get_m_Hyp(px_He, py_He, pz_He, px_Pi, py_Pi, pz_Pi, px_Hyp, py_Hyp, pz_Hyp) #GeV

    FourVector_Hyp = get_CMSFourVector(m_Hyp,px_Hyp, py_Hyp, pz_Hyp, m_Hyp, px_Hyp, py_Hyp, pz_Hyp )
    FourVector_He3 = get_CMSFourVector(m_He,px_He, py_He, pz_He, m_Hyp, px_Hyp, py_Hyp, pz_Hyp )
    FourVector_Pi = get_CMSFourVector(m_Pi,px_Pi, py_Pi, pz_Pi, m_Hyp, px_Hyp, py_Hyp, pz_Hyp )

if __name__ == "__main__":
    # Simulated Hypertritons
    labFrameH = TreeHandler()
    labFrameH.get_handler_from_large_file(
        file_name="../Data/SignalTable_B_20g7.root", 
        tree_name='GenTable',
    )   

    m_He = 2.80839 #GeV
    m_Pi = 0.139570 #GeV
    
    df=labFrameH.get_data_frame()
    
        


    # TODO: cut on pt, pseudorapidity, and transversemomentum
    # TODO: calculate cos theta/theta* for two different axis, once beam and once hypertriton momentum
    # TODO: plot ?
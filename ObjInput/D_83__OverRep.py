# -*- coding: utf-8 -*-
###############################################################################
# --- D_83__OverRep.py --------------------------------------------------------
###############################################################################
import Core.C_00__GenConstants as GC

# --- general -----------------------------------------------------------------
sNmSpec = GC.S_OV_REP_L

# --- data specific (general) -------------------------------------------------
cSep = ';'
sMCorrectL = GC.S_M_CORR_BON    # None / S_M_CORR_BON
sMCorrectS = sMCorrectL[:3]

# --- data specific (Pho data) ------------------------------------------------
nMinPD, nMaxPD = 1, 832
dSrtPD = {'Idx': {GC.S_MN_CONC: {'Asc': False}},
          'Col': {'Srt': ('float', 0), 'Asc': True}}
lElCPD = [GC.S_BIN_C_2]

# --- data specific ("BinOps" from basic data) --------------------------------
nMinBO, nMaxBO = 1, round(42*3330/10)
# nMinBO, nMaxBO = 1, 1000
dSrtBO = {'Idx': {GC.S_PNCI: {'Asc': False},
                  GC.S_PCI: {'Asc': False},
                  GC.S_SPEAR_V: {'Asc': False},
                  GC.S_CORR_V: {'Asc': False}},
          'Col': {'Srt': ('float', 0), 'Asc': True}}
lElCBO = [GC.S_MET_D]     # columns to evaluate
                            # None: eval. MetD and PhoD columns of BO file
                            # [GC.S_MET_D]: eval. MetD column of BO file
                            # [GC.S_PHO_D]: eval. PhoD column of BO file
                            # [GC.S_BIN_C_2]: eval. BinCode2 column of BO file

# --- data specific ("BinOps" from cluster result data) -----------------------
nMinCR, nMaxCR = 1, 10
dSrtCR = {'Idx': {GC.S_CORR_V: {'Asc': False}},
          'Col': {'Srt': ('int', 3), 'Asc': True}}
lElCCR = None

# --- dictionaries containing profile-type specific input ---------------------
lTpX = [GC.S_PD, GC.S_BO]
lTpY = [GC.S_N_OCC, GC.S_OVER_REP, GC.S_UNDER_REP]
lSXAx = ['Top $\it{n}$ of the most frequent phosphopeptides',
         'Top $\it{n}$ of the strongest interactions']
lSYAx = [GC.S_Y_N_OCC, GC.S_Y_P_VAL, GC.S_Y_P_VAL]
lNDigRndYAx = [0, 2, 2]
lDoPYAx = [False, True, True]

# --- names and paths of files and dirs ---------------------------------------
dISPr = {'No': {0: ''}, 'Tr': {0: ''}, 'Dv': {0: ''}}
sNOcc = GC.S_N_OCC_ABS
sPValOv = GC.S_P_VAL_OV
sPValUn = GC.S_P_VAL_UN
sPOf = GC.S_P_OF

# --- graphics parameters -----------------------------------------------------
nmPlt_Prof = 'Profile'  # name prefix of the plot
thrProf = 0.05          # plot threshold for profiles
sComp = '>='            # comparison string (value with threshold)
szFontLeg = 'small'     # font size of legend
iIncr = 1               # increase of file number
jIncr = 10              # number of entities (e.g. metabolites) per plot
coordAnchorBox = (1.1, 0.5)         # coordinates of the legend anchor box
dClrSpec = {'1.1': (0.12, 0.47, 0.71),
            '1.2': (1.0, 0.5, 0.05),
            '1.3': (0.17, 0.63, 0.17),
            '2.1': (0.84, 0.15, 0.16),
            '2.2': (0.58, 0.4, 0.74),
            '3.1': (0.55, 0.34, 0.29),
            '3.2': (0.89, 0.47, 0.76),
            '3.4': (0.5, 0.5, 0.5),
            '3.5': (0.74, 0.74, 0.13),
            '3.6': (0.09, 0.75, 0.81),
            '3.99': (0.12, 0.47, 0.71),
            '4.1': (1.0, 0.5, 0.05),
            '4.2': (0.17, 0.63, 0.17),
            '4.3': (0.84, 0.15, 0.16),
            '5.3': (0.58, 0.4, 0.74),
            '6.4': (0.55, 0.34, 0.29),
            '6.5': (0.89, 0.47, 0.76),
            '7.1': (0.5, 0.5, 0.5),
            '8.1': (0.74, 0.74, 0.13),
            '8.2': (0.09, 0.75, 0.81),
            '8.3': (0.12, 0.47, 0.71),
            '10.1': (1.0, 0.5, 0.05),
            '10.2': (0.17, 0.63, 0.17),
            '10.3': (0.84, 0.15, 0.16),
            '10.5': (0.58, 0.4, 0.74),
            '10.8': (0.55, 0.34, 0.29),
            '11.1': (0.89, 0.47, 0.76),
            '11.3': (0.5, 0.5, 0.5),
            '11.8': (0.74, 0.74, 0.13),
            '11.9': (0.09, 0.75, 0.81),
            '12.1': (0.12, 0.47, 0.71),
            '12.2': (1.0, 0.5, 0.05),
            '12.4': (0.17, 0.63, 0.17),
            '13.1': (0.84, 0.15, 0.16),
            '13.2': (0.58, 0.4, 0.74),
            '15': (0.55, 0.34, 0.29),
            '15.1': (0.89, 0.47, 0.76),
            '15.2': (0.5, 0.5, 0.5),
            '16.2': (0.74, 0.74, 0.13),
            '16.5': (0.09, 0.75, 0.81),
            '16.8': (0.12, 0.47, 0.71),
            '17.1': (1.0, 0.5, 0.05),
            '17.2': (0.17, 0.63, 0.17),
            '17.3': (0.84, 0.15, 0.16),
            '17.5': (0.58, 0.4, 0.74),
            '17.6': (0.55, 0.34, 0.29),
            '18.2': (0.89, 0.47, 0.76),
            '18.4': (0.5, 0.5, 0.5),
            '20': (0.74, 0.74, 0.13),
            '20.1': (0.09, 0.75, 0.81),
            '20.2': (0.12, 0.47, 0.71),
            '21.1': (1.0, 0.5, 0.05),
            '21.2': (0.17, 0.63, 0.17),
            '21.4': (0.84, 0.15, 0.16),
            '21.5': (0.58, 0.4, 0.74),
            '21.6': (0.55, 0.34, 0.29),
            '23.1': (0.89, 0.47, 0.76),
            '23.2': (0.5, 0.5, 0.5),
            '23.3': (0.74, 0.74, 0.13),
            '23.4': (0.09, 0.75, 0.81),
            '24.1': (0.12, 0.47, 0.71),
            '25': (1.0, 0.5, 0.05),
            '25.6': (0.17, 0.63, 0.17),
            '26.13': (0.84, 0.15, 0.16),
            '26.16': (0.58, 0.4, 0.74),
            '26.17': (0.55, 0.34, 0.29),
            '26.22': (0.89, 0.47, 0.76),
            '26.3': (0.5, 0.5, 0.5),
            '26.4': (0.74, 0.74, 0.13),
            '26.5': (0.09, 0.75, 0.81),
            '26.6': (0.12, 0.47, 0.71),
            '26.7': (1.0, 0.5, 0.05),
            '26.8': (0.17, 0.63, 0.17),
            '26.9': (0.84, 0.15, 0.16),
            '27.1': (0.58, 0.4, 0.74),
            '27.2': (0.55, 0.34, 0.29),
            '27.3': (0.89, 0.47, 0.76),
            '27.4': (0.5, 0.5, 0.5),
            '28.1': (0.74, 0.74, 0.13),
            '28.2': (0.09, 0.75, 0.81),
            '28.99': (0.12, 0.47, 0.71),
            '29.1': (1.0, 0.5, 0.05),
            '29.2': (0.17, 0.63, 0.17),
            '29.3': (0.84, 0.15, 0.16),
            '29.4': (0.58, 0.4, 0.74),
            '29.5': (0.55, 0.34, 0.29),
            '29.6': (0.89, 0.47, 0.76),
            '30.1': (0.5, 0.5, 0.5),
            '30.11': (0.74, 0.74, 0.13),
            '30.2': (0.09, 0.75, 0.81),
            '30.3': (0.12, 0.47, 0.71),
            '30.4': (1.0, 0.5, 0.05),
            '30.5': (0.17, 0.63, 0.17),
            '30.6': (0.84, 0.15, 0.16),
            '30.7': (0.58, 0.4, 0.74),
            '30.8': (0.55, 0.34, 0.29),
            '31.1': (0.89, 0.47, 0.76),
            '31.2': (0.5, 0.5, 0.5),
            '31.3': (0.74, 0.74, 0.13),
            '31.4': (0.09, 0.75, 0.81),
            '31.5': (0.12, 0.47, 0.71),
            '32': (1.0, 0.5, 0.05),
            '33.1': (0.17, 0.63, 0.17),
            '33.99': (0.84, 0.15, 0.16),
            '34.1': (0.58, 0.4, 0.74),
            '34.11': (0.55, 0.34, 0.29),
            '34.12': (0.89, 0.47, 0.76),
            '34.13': (0.5, 0.5, 0.5),
            '34.14': (0.74, 0.74, 0.13),
            '34.15': (0.09, 0.75, 0.81),
            '34.16': (0.12, 0.47, 0.71),
            '34.17': (1.0, 0.5, 0.05),
            '34.18': (0.17, 0.63, 0.17),
            '34.19': (0.84, 0.15, 0.16),
            '34.2': (0.58, 0.4, 0.74),
            '34.21': (0.55, 0.34, 0.29),
            '34.3': (0.89, 0.47, 0.76),
            '34.4': (0.5, 0.5, 0.5),
            '34.5': (0.74, 0.74, 0.13),
            '34.7': (0.09, 0.75, 0.81),
            '34.8': (0.12, 0.47, 0.71),
            '34.9': (1.0, 0.5, 0.05),
            '34.98': (0.17, 0.63, 0.17),
            '34.99': (0.84, 0.15, 0.16),
            '35.1': (0.58, 0.4, 0.74),
            '35.2': (0.55, 0.34, 0.29)}

dClrBinC = {'2.1': (0.12, 0.47, 0.71),
            '4.1': (1.0, 0.5, 0.05),
            '17.2': (0.17, 0.63, 0.17),
            '29.2': (0.84, 0.15, 0.16),
            '31.4': (0.58, 0.4, 0.74),
            '33.99': (0.55, 0.34, 0.29)}

# --- assertions --------------------------------------------------------------
assert lElCBO is None or type(lElCBO) == list
if lElCBO is not None:
    for s in lElCBO:
        assert s in [GC.S_MET_D, GC.S_PHO_D] + GC.L_S_PHO_CL
if lElCBO == None or lElCBO == [GC.S_MET_D] or lElCBO == [GC.S_PHO_D]:
    dSrtBO['Col']['Srt'] = ('str', 0)
elif [s in [GC.S_BIN_C_2, GC.S_BIN_C_1] for s in lElCBO] == [True]*len(lElCBO):
    dSrtBO['Col']['Srt'] = ('float', 0)

assert len(lSXAx) == len(lTpX)
assert (len(lSYAx) == len(lTpY) and len(lNDigRndYAx) == len(lTpY) and
        len(lDoPYAx) == len(lTpY))

# --- create input dictionary -------------------------------------------------
dIO = {# --- general
       'sNmSpec': sNmSpec,
       # --- data specific (general)
       'cSep': cSep,
       'sMCorrectL': sMCorrectL,
       'sMCorrectS': sMCorrectS,
       # --- data specific (Pho data)
       'nMinPD': nMinPD,
       'nMaxPD': nMaxPD,
       'dSrtPD': dSrtPD,
       'lElCPD': lElCPD,
       # --- data specific ("BinOps" from basic data)
       'nMinBO': nMinBO,
       'nMaxBO': nMaxBO,
       'dSrtBO': dSrtBO,
       'lElCBO': lElCBO,
       # --- data specific ("BinOps" from cluster result data)
       'nMinCR': nMinCR,
       'nMaxCR': nMaxCR,
       'dSrtCR': dSrtCR,
       'lElCCR': lElCCR,
       # --- dictionaries containing profile-type specific input
       'lTpX': lTpX,
       'lTpY': lTpY,
       'dTpX': {lTpX[k]: lSXAx[k] for k in range(len(lTpX))},
       'dTpY': {lTpY[k]: (lSYAx[k], lNDigRndYAx[k], lDoPYAx[k]) for k in
                range(len(lTpY))},
       # --- names and paths of files and dirs
       'dISPr': dISPr,
       'sNOcc': sNOcc,
       'sPValOv': sPValOv,
       'sPValUn': sPValUn,
       'sPOf': sPOf,
       'sPValOvPOf': sPValOv + '_' + sPOf,
       'sPValUnPOf': sPValUn + '_' + sPOf,
       # --- graphics parameters
       'nmPlt_Prof': nmPlt_Prof,
       'thrProf': thrProf,
       'sComp': sComp,
       'szFontLeg': szFontLeg,
       'iIncr': iIncr,
       'jIncr': jIncr,
       'coordAnchorBox': coordAnchorBox,
       'dClrSpec': dClrSpec,
       'dClrBinC': dClrBinC}

###############################################################################

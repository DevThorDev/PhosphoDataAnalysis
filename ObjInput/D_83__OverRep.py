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
dSrtBO = {'Idx': {GC.S_IPDC: {'Asc': False},
                  GC.S_SPDS: {'Asc': False},
                  GC.S_CORR_V: {'Asc': False}},
          'Col': {'Srt': ('float', 0), 'Asc': True}}
lElCBO = [GC.S_BIN_C_2]     # columns to evaluate
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
       'coordAnchorBox': coordAnchorBox}

###############################################################################

# -*- coding: utf-8 -*-
###############################################################################
# --- S__PlotDerivIC.py -------------------------------------------------------
# plot the derivation of the concordance index IC from the measurements
# in three steps
###############################################################################
import os, time, math, itertools

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# --- CONSTANTS ---------------------------------------------------------------
S_SPACE = ' '
S_USC = '_'
S_BAR = '|'
S_SLASH = '/'
S_DOT = '.'
S_CSV = 'csv'
S_PDF = 'pdf'
S_SPSLSP = S_SPACE + S_SLASH + S_SPACE

S_M = 'M'
S_P = 'P'
S_I = 'I'
S_C = 'C'
L_S_SPC = [S_M, S_P, S_I, S_C]

S_0 = '0'
S_1 = '1'
S_2 = '2'
S_3 = '3'
S_4 = '4'
S_5 = '5'
S_6 = '6'
L_S_0_6 = [S_0, S_1, S_2, S_3, S_4, S_5, S_6]

S_MET = 'Metabolite'
S_METS = S_MET + 's'
S_PHO = 'Phosphopeptide'
S_PHOS = S_PHO + 's'

S_IC = 'IC'
S_IC_MEX = '$I_C$'
S_IC_M = S_IC + '(-)'
S_IC_P = S_IC + '(+)'
L_S_IC_ALL = [S_IC_M, S_IC_P, S_IC]
S_IC_M_MEX = '$I_{C(-)}$'
S_IC_P_MEX = '$I_{C(+)}$'
S_IC_ALL_MEX = S_IC_M_MEX + S_SPSLSP + S_IC_P_MEX + S_SPSLSP + S_IC_MEX
L_S_IC_ALL_MEX = [S_IC_M_MEX, S_IC_P_MEX, S_IC_MEX]

S_GT = 'GT'
S_GT0 = S_GT + S_0
S_GT1 = S_GT + S_1
S_GT5 = S_GT + S_5
L_S_GT = [S_GT0, S_GT1, S_GT5]

NM_GT0 = 'wild type'
NM_GT1 = '$\it{pgm}$ mutant'
NM_GT5 = '$\it{sweet11/12}$ mutant'
D_NM_GT = {S_GT0: NM_GT0, S_GT1: NM_GT1, S_GT5: NM_GT5}

S_FT1 = 'DR'
S_FT2 = 'DS'
S_FT3 = 'NR'
S_FT4 = 'NS'
L_S_FT = [S_FT1, S_FT2, S_FT3, S_FT4]
L_S_FT_CHG = [L_S_FT[i] + S_BAR + L_S_FT[j] for j in range(len(L_S_FT))
              for i in range(len(L_S_FT)) if i != j]
D_S_FT_CHG = {s: [S_USC.join([t, s]) for t in L_S_FT_CHG] for s in L_S_SPC}

S_CLR0 = S_C + S_0
S_CLR1 = S_C + S_1
S_CLR2 = S_C + S_2
S_CLR3 = S_C + S_3
S_CLR4 = S_C + S_4
S_CLR5 = S_C + S_5
L_S_CLR = [S_CLR0, S_CLR1, S_CLR2, S_CLR3, S_CLR4, S_CLR5]

S_DEV_SD = 'DeviationSD'
S_DERIV = 'Deriv'
S_IC_DERIV = S_IC + S_DERIV
S_CMB_TO_IC = 'CombineToIC'

S_CMB_C = 'Combined'
S_DEV_S = 'deviation'
S_DEVS_C = 'Deviations between features'
S_CMP_S = 'component'
S_CMB_DEV = S_CMB_C + S_SPACE + S_DEV_S
S_IC_CMP = S_IC_MEX + S_SPACE + S_CMP_S

S_YLBL_DEV_SD_PLT = 'Difference as multiples of standard deviation'
S_YLBL_IC_DERIV_PLT = S_DEVS_C + S_SPSLSP + S_IC_CMP
S_YLBL1_CMB_TO_IC_PLT = S_IC_CMP

S_BASE_CL = 'BaseClass'
S_INP_DATA = 'InputData'
S_ROOT_CL = 'RootClass'
S_PLTR = 'Plotter'
S_DEV_SD_PLTR = S_DEV_SD + S_PLTR
S_IC_DERIV_PLTR = S_IC_DERIV + S_PLTR
S_CMB_TO_IC_PLTR = S_CMB_TO_IC + S_PLTR

S_PLT = 'Plot'
NM_DEV_SD_PLT = S_1 + S_USC + S_DEV_SD + S_PLT
NM_IC_DERIV_PLT = S_2 + S_USC + S_IC_DERIV + S_PLT
NM_CMB_TO_IC_PLT1 = S_3 + S_USC + S_1 + S_USC + S_CMB_TO_IC + S_PLT
NM_CMB_TO_IC_PLT2 = S_3 + S_USC + S_2 + S_USC + S_CMB_TO_IC + S_PLT

STY_LN_SOL = '-'
STY_LN_DSH = '--'
STY_LN_DDS = '-.'

EPS = .01
R04 = 4

# --- INPUT -------------------------------------------------------------------
# --- flow control ------------------------------------------------------------
doPlot_DevSD = False             # True / False
doPlot_ICDrv = False             # True / False
doPlot_CmbToIC = True           # True / False

# --- data specific input -----------------------------------------------------
maxSLen = 20
sSep = ';'

# --- graphics parameters / all plots -----------------------------------------
clrDef = 'k'                    # default colour
szFontLeg = 'small'             # font size of legend
nCharDsp = 60                   # number of chars displayed for legend item
posLegXY = (.5, 1.02)           # coordinates of the legend anchor box

# --- graphics parameters / deviation SD plot ---------------------------------
dTupIn_DevSD = {'A': (S_GT0, 'Alanine', S_FT2 + S_BAR + S_FT1),
                'B': (S_GT0, 'Alanine', S_FT2 + S_BAR + S_FT4),
                'C': (S_GT1, 'Docosanoic_acid', S_FT4 + S_BAR + S_FT1),
                'D': (S_GT1, 'Docosanoic_acid', S_FT1 + S_BAR + S_FT4),
                'E': (S_GT5, 'Putrescine', S_FT2 + S_BAR + S_FT1),
                'F': (S_GT5, 'Putrescine', S_FT3 + S_BAR + S_FT2),

                'N': (S_GT0, 'Docosanoic_acid', S_FT4 + S_BAR + S_FT3),
                'O': (S_GT1, 'Alanine', S_FT4 + S_BAR + S_FT1),
                'P': (S_GT1, 'Docosanoic_acid', S_FT4 + S_BAR + S_FT1),

                'X_S': (S_GT5, 'Sucrose', S_FT2 + S_BAR + S_FT1)}

nmPlt_DevSD = NM_DEV_SD_PLT         # name prefix of the deviation SD plot
tFigSz_DevSD = (2., 5.)             # (width, height): figure size [inches]
symMark_DevSD = 'x'                 # marker symbol
szMark_DevSD = 25                   # marker size
lenMnBar_DevSD = .3                 # length of mean bar
lwdMnBar_DevSD = 1.5                # line width of mean bar
dltLwdMnBar_DevSD = .5              # delta of width of black line to col. line
lenWMnBar_DevSD = .05               # length of mean bar whiskers
lwdWMnBar_DevSD = .6                # line width of mean bar whiskers
lwdLnSD_DevSD = .8                  # line width of mean bar
dltLwdLnSD_DevSD = .25              # delta of width of black line to col. line
lenWLnSD_DevSD = .1                 # length of mean bar whiskers
lwdWLnSD_DevSD = .4                 # line width of mean bar whiskers
lstCnLnSD_DevSD = STY_LN_DSH        # line style of connecting line
lwdCnLnSD_DevSD = .8                # line width of connecting line
dPosFt_DevSD = {0: .75, 1: .25}     # dictionary of positions of features 0|1
axXLim_DevSD = (0., 1.)             # limits for the x-axis, or None
axYLim_DevSD = None                 # limits for the y-axis, or None
adaptYLim_DevSD = True              # adapt y-limits using values to plot?
sXLbl_DevSD = None                  # x-label string, or None
sYLbl_DevSD = S_YLBL_DEV_SD_PLT     # y-label string, or None
lblXTck_DevSD = None                # labels of x-ticks
stepYTck_DevSD = 1                  # step distance for y-ticks, or None
plotTtl_DevSD = True                # add title to plot?
dFontTtl_DevSD = {'size': 10}       # font dictionary of the title
locTtl_DevSD = 'left'               # location of the title
padTtl_DevSD = 10.                  # padding of the title
degRotXLbl_DevSD = 90               # degree rotation of x-labels

# --- graphics parameters / IC derivation plot --------------------------------
dTupIn_ICDrv = {'A': (S_GT0, 'Alanine', 'ESLS(1)PGQQHVSQNTAVKPEGR'),
                'B': (S_GT0, 'Alanine', 'T(0.001)PS(0.996)QRPS(0.003)TSSSSGGFNIGK'),
                'C': (S_GT1, 'Glutamic_acid', 'TADS(1)DGES(1)GEIKFEDNNLPPLGEK'),
                'D': (S_GT5, 'Putrescine', 'S(0.999)FS(0.001)VADFPR'),
                'E': (S_GT5, 'Fructose', 'TEEDENDDEDHEHKDDKT(0.854)S(0.144)PDS(0.001)IVMVEAK'),

                'N': (S_GT0, 'Docosanoic_acid', 'SLEELS(1)GEAEVS(1)HDEK'),
                'O': (S_GT1, 'Alanine', 'TFDELS(1)DGEVYEDS(1)D'),
                'P': (S_GT1, 'Docosanoic_acid', 'S(0.003)PS(0.997)YKEVALAPPGSIAK'),

                'X_S1': (S_GT5, 'Sucrose', 'S(1)DLQTPLVRPK'),
                'X_S2': (S_GT5, 'Sucrose', 'NFANS(1)FGRK'),

                'Y_GT0_1': (S_GT0, 'Nonanoic_acid', 'T(0.011)FDELS(0.759)DT(0.23)EVYEDS(1)D'),
                'Y_GT0_2': (S_GT0, 'myo-Inositol', 'TFDELS(1)DTEVYEDS(1)D'),
                'Y_GT1_1': (S_GT1, 'Phenylalanine', 'S(0.825)RS(0.137)VDES(0.039)FANSFSPR'),
                'Y_GT1_2': (S_GT1, 'Phenylalanine', 'S(0.001)GRT(0.004)S(0.996)EPNS(1)EDEAAGVGK'),
                'Y_GT5_1': (S_GT5, 'Aspartic_acid', 'IGS(0.999)S(0.001)EMLIEGEDVR'),
                'Y_GT5_2': (S_GT5, 'Beta-alanine', 'DIS(1)PTAAGLGLPVTGGK'),
                'Y_GT5_3': (S_GT5, 'Docosanoic_acid', 'LSRPGS(1)GS(1)VSGLASQR'),
                'Y_GT5_4': (S_GT5, 'Ornithine', 'TDSEVTSLAAS(0.024)S(0.976)PARS(1)PR')}

nmPlt_ICDrv = NM_IC_DERIV_PLT       # name prefix of the IC derivation plot
tFigSz_ICDrv = (6., 4.)             # (width, height): figure size [inches]
lwdPlt_ICDrv = .75                  # line width in plot
lClr_ICDrv  = L_S_CLR[:len(L_S_SPC)]    # list of colours for barplot
axXLim_ICDrv = None                 # limits for the x-axis, or None
axYLim_ICDrv = (-11, 11)            # limits for the y-axis, or None
plotTtl_ICDrv = False               # add title to plot?
dFontTtl_ICDrv = {'size': 10}       # font dictionary of the title
locTtl_ICDrv = 'left'               # location of the title
padTtl_ICDrv = 70.                  # padding of the title
locLegend_ICDrv = 'lower center'    # location of the legend
wdthGrp_ICDrv = .8                  # width of bar group
wdthBar_ICDrv = wdthGrp_ICDrv/len(D_S_FT_CHG) - EPS    # width of single bars
sXLbl_ICDrv = None                  # x-label string, or None
sYLbl_ICDrv = S_YLBL_IC_DERIV_PLT   # y-label string, or None
lblXTck_ICDrv = L_S_FT_CHG          # labels of x-ticks
stepYTck_ICDrv = 2                  # step distance for y-ticks, or None
degRotXLbl_ICDrv = 90               # degree rotation of x-labels
plotVLines_ICDrv = True             # plot vertical lines between groups?
lwdVLine_ICDrv = .2                 # width of vertical line
clrVLine_ICDrv = clrDef             # colour of vertical line

# --- graphics parameters / combine to IC plot B(ase) -------------------------
dTupIn_CmbToIC = dTupIn_ICDrv
lwdPlt_CmbToIC = .75                # line width in plot
dClr_CmbToIC = {S_IC_M: (.557, .153, .706),    # colour dictionary
                S_IC_P: (.839, .153, .471),
                S_IC: S_CLR3}
dFontTtl_CmbToIC = {'size': 12}     # font dictionary of the title
locTtl_CmbToIC = 'left'             # location of the title
padTtl_CmbToIC = 10.                # padding of the title
sXLbl_CmbToIC = None                # x-label string, or None
lblXTck_CmbToIC = None              # labels of x-ticks
degRotXLbl_CmbToIC = 90             # degree rotation of x-labels
lwdVLine_CmbToIC = .2               # width of vertical line
clrVLine_CmbToIC = clrDef           # colour of vertical line

# --- graphics parameters / combine to IC plot 1 ------------------------------
nmPlt_CmbToIC1 = NM_CMB_TO_IC_PLT1  # name prefix of the IC derivation plot 1
tFigSz_CmbToIC1 = (5., 4.)          # (width, height): figure size [inches]
alphaClr_CmbToIC1 = .8              # alpha of colours defined via tuples
axXLim_CmbToIC1 = None              # limits for the x-axis, or None
axYLim_CmbToIC1 = (-2, 2)           # limits for the y-axis, or None
plotTtl_CmbToIC1 = False             # add title to plot?
sYLbl_CmbToIC1 = S_YLBL1_CMB_TO_IC_PLT    # y-label string, or None
wdthBar_CmbToIC1 = .6 - EPS         # width of bars
lblXTck_CmbToIC1 = L_S_FT_CHG       # labels of x-ticks
stepYTck_CmbToIC1 = .5              # step distance for y-ticks, or None
plotVLines_CmbToIC1 = True          # plot vertical lines between groups?

# --- graphics parameters / combine to IC plot 2 ------------------------------
nmPlt_CmbToIC2 = NM_CMB_TO_IC_PLT2  # name prefix of the IC derivation plot 2
# tFigSz_CmbToIC2 = (3., 4.)          # (width, height): figure size [inches]
tFigSz_CmbToIC2 = (1.5, 1.5)        # (width, height): figure size [inches]
lwdBarIC_CmbToIC2 = 1.5             # line width of final IC bar
hatBarIC_CmbToIC2 = 'xx'            # hatch of final IC bar
alphaClr_CmbToIC2 = .8              # alpha of colours defined via tuples
alphaClrLow_CmbToIC2 = .4           # alpha of colours defined via tuples
axXLim_CmbToIC2 = None              # limits for the x-axis, or None
axYLim_CmbToIC2 = None              # limits for the y-axis, or None
plotTtl_CmbToIC2 = False            # add title to plot?
sYLbl_CmbToIC2 = None               # y-label string, or None
wdthBar_CmbToIC2 = .6 - EPS         # width of bars
lblXTck_CmbToIC2 = L_S_IC_ALL_MEX   # labels of x-ticks
stepYTck_CmbToIC2 = None            # step distance for y-ticks, or None
plotVLines_CmbToIC2 = True          # plot vertical lines between groups?

# --- names and paths of files and dirs ---------------------------------------
sMs, sPs, s1, s2, s3 = S_METS, S_PHOS, 'Means', 'SDs', 'Deviations'
dSFIn_Ms = {(sMs, sGT): S_USC.join([sMs, s1, s2, s3, sGT]) for sGT in L_S_GT}
dSFIn_Ps = {(sPs, sGT): S_USC.join([sPs, s1, s2, s3, sGT]) for sGT in L_S_GT}
dSFIn_DevSD = dSFIn_Ms | dSFIn_Ps
s1, s2, s3, s4, s5, s6 = 'Corr', 'BinOp', 'MetD', 'PhoD', 'DvSD', 'AllD'
dSFIn_ICDrv = {sGT: (S_IC_DERIV + S_USC*2 + s1 + S_USC*2 +
                     S_USC.join([s2, s3, s5, sGT, s6, s4, s5, sGT, s6]))
               for sGT in L_S_GT}
sFOut_DevSD = nmPlt_DevSD
sFOut_ICDrv = nmPlt_ICDrv
sF1Out_CmbToIC = nmPlt_CmbToIC1
sF2Out_CmbToIC = nmPlt_CmbToIC2

sDirInCSV = 'InpData'
sDirOutPDF = 'OutPDF'

pBaseIn = os.path.join('..', '..', '..', '25_Papers', '01_FirstAuthor',
                       '04_SysBio_DataAnalysis', '15_Figures',
                       'A__IC_Derivation')
pBaseOut = os.path.join('..', '..', '..', '25_Papers', '01_FirstAuthor',
                        '04_SysBio_DataAnalysis', '15_Figures',
                        'A__IC_Derivation')
# pBaseOut = os.path.join('..', '..', '..', '..', '..', '..', 'V')
pInCSV = os.path.join(pBaseIn, sDirInCSV)
pOutPDF = os.path.join(pBaseOut, sDirOutPDF)

# --- derived values ----------------------------------------------------------
dPFIn_DevSD = {tK: os.path.join(pInCSV, dSFIn_DevSD[tK] + S_DOT + S_CSV)
               for tK in dSFIn_DevSD}
dPFIn_ICDrv = {sGT: os.path.join(pInCSV, dSFIn_ICDrv[sGT] + S_DOT + S_CSV)
               for sGT in L_S_GT}

# --- assertions --------------------------------------------------------------
assert len(L_S_SPC) == 4 and len(D_S_FT_CHG) == 4

# --- INPUT DICTIONARY --------------------------------------------------------
dInput = {# --- constants
          'lSGT': L_S_GT,
          'sBase': S_BASE_CL,
          'sInpDat': S_INP_DATA,
          'sRoot': S_ROOT_CL,
          'sPltr': S_PLTR,
          'sPltr_DevSD': S_DEV_SD_PLTR,
          'sPltr_ICDrv': S_IC_DERIV_PLTR,
          'sPltr_CmbToIC': S_CMB_TO_IC_PLTR,
          'sMet': S_MET,
          'sPho': S_PHO,
          'R04': R04,
          # --- flow control
          'doPlot_DevSD': doPlot_DevSD,
          'doPlot_ICDrv': doPlot_ICDrv,
          'doPlot_CmbToIC': doPlot_CmbToIC,
          # --- data specific input
          'maxSLen': maxSLen,
          'sSep': sSep,
          # --- graphics parameters / all plots
          'plot_Gen': {'clrDef': clrDef,
                       'szFontLeg': szFontLeg,
                       'nCharDsp': nCharDsp,
                       'posLegXY': posLegXY},
          # --- graphics parameters / deviation SD plot
          'plot_DevSD': {'dTupIn': dTupIn_DevSD,
                         'nmPlt': nmPlt_DevSD,
                         'tFigSz': tFigSz_DevSD,
                         'symMark': symMark_DevSD,
                         'szMark': szMark_DevSD,
                         'lenMnBar': lenMnBar_DevSD,
                         'lenHfMnBar': lenMnBar_DevSD/2.,
                         'lwdMnBar': lwdMnBar_DevSD,
                         'dltLwdMnBar': dltLwdMnBar_DevSD,
                         'lenWMnBar': lenWMnBar_DevSD,
                         'lwdWMnBar': lwdWMnBar_DevSD,
                         'lwdLnSD': lwdLnSD_DevSD,
                         'dltLwdLnSD': dltLwdLnSD_DevSD,
                         'lenWLnSD': lenWLnSD_DevSD,
                         'lwdWLnSD': lwdWLnSD_DevSD,
                         'lstCnLnSD': lstCnLnSD_DevSD,
                         'lwdCnLnSD': lwdCnLnSD_DevSD,
                         'dPosFt': dPosFt_DevSD,
                         'axXLim': axXLim_DevSD,
                         'axYLim': axYLim_DevSD,
                         'sXLbl': sXLbl_DevSD,
                         'sYLbl': sYLbl_DevSD,
                         'lblXTck': lblXTck_DevSD,
                         'stepYTck': stepYTck_DevSD,
                         'adaptYLim': adaptYLim_DevSD,
                         'plotTtl': plotTtl_DevSD,
                         'dFontTtl': dFontTtl_DevSD,
                         'locTtl': locTtl_DevSD,
                         'padTtl': padTtl_DevSD,
                         'degRotXLbl': degRotXLbl_DevSD},
          # --- graphics parameters / IC derivation plot
          'plot_ICDrv': {'dTupIn': dTupIn_ICDrv,
                         'nmPlt': nmPlt_ICDrv,
                         'tFigSz': tFigSz_ICDrv,
                         'lwdPlt': lwdPlt_ICDrv,
                         'lClr': lClr_ICDrv,
                         'axXLim': axXLim_ICDrv,
                         'axYLim': axYLim_ICDrv,
                         'plotTtl': plotTtl_ICDrv,
                         'dFontTtl': dFontTtl_ICDrv,
                         'locTtl': locTtl_ICDrv,
                         'padTtl': padTtl_ICDrv,
                         'locLegend': locLegend_ICDrv,
                         'wdthGrp': wdthGrp_ICDrv,
                         'wdthBar': wdthBar_ICDrv,
                         'sXLbl': sXLbl_ICDrv,
                         'sYLbl': sYLbl_ICDrv,
                         'lblXTck': lblXTck_ICDrv,
                         'stepYTck': stepYTck_ICDrv,
                         'degRotXLbl': degRotXLbl_ICDrv,
                         'plotVLines': plotVLines_ICDrv,
                         'lwdVLine': lwdVLine_ICDrv,
                         'clrVLine': clrVLine_ICDrv},
          # --- graphics parameters / combine to IC plot B(ase)
          'plotB_CmbToIC': {'dTupIn': dTupIn_ICDrv,
                            'lwdPlt': lwdPlt_CmbToIC,
                            'dClr': dClr_CmbToIC,
                            'dFontTtl': dFontTtl_CmbToIC,
                            'locTtl': locTtl_CmbToIC,
                            'padTtl': padTtl_CmbToIC,
                            'sXLbl': sXLbl_CmbToIC,
                            'lblXTck': lblXTck_CmbToIC,
                            'degRotXLbl': degRotXLbl_CmbToIC,
                            'lwdVLine': lwdVLine_CmbToIC,
                            'clrVLine': clrVLine_CmbToIC},
          # --- graphics parameters / combine to IC plot 1
          'plot1_CmbToIC': {'nmPlt': nmPlt_CmbToIC1,
                            'tFigSz': tFigSz_CmbToIC1,
                            'alphaClr': alphaClr_CmbToIC1,
                            'axXLim': axXLim_CmbToIC1,
                            'axYLim': axYLim_CmbToIC1,
                            'plotTtl': plotTtl_CmbToIC1,
                            'sYLbl': sYLbl_CmbToIC1,
                            'wdthBar': wdthBar_CmbToIC1,
                            'lblXTck': lblXTck_CmbToIC1,
                            'stepYTck': stepYTck_CmbToIC1,
                            'plotVLines': plotVLines_CmbToIC1},
          # --- graphics parameters / combine to IC plot 2
          'plot2_CmbToIC': {'nmPlt': nmPlt_CmbToIC2,
                            'tFigSz': tFigSz_CmbToIC2,
                            'lwdBarIC': lwdBarIC_CmbToIC2,
                            'hatBarIC': hatBarIC_CmbToIC2,
                            'alphaClr': alphaClr_CmbToIC2,
                            'alphaClrLow': alphaClrLow_CmbToIC2,
                            'axXLim': axXLim_CmbToIC2,
                            'axYLim': axYLim_CmbToIC2,
                            'plotTtl': plotTtl_CmbToIC2,
                            'sYLbl': sYLbl_CmbToIC2,
                            'wdthBar': wdthBar_CmbToIC2,
                            'lblXTck': lblXTck_CmbToIC2,
                            'stepYTck': stepYTck_CmbToIC2,
                            'plotVLines': plotVLines_CmbToIC2},
          # --- names and paths of files and dirs
          'pInCSV': pInCSV,
          'pOutPDF': pOutPDF,
          'dSFIn_DevSD': dSFIn_DevSD,
          'dSFIn_ICDrv': dSFIn_ICDrv,
          'sFOut_DevSD': sFOut_DevSD,
          'sFOut_ICDrv': sFOut_ICDrv,
          'sF1Out_CmbToIC': sF1Out_CmbToIC,
          'sF2Out_CmbToIC': sF2Out_CmbToIC,
          # --- further derived values
          'dPFIn_DevSD': dPFIn_DevSD,
          'dPFIn_ICDrv': dPFIn_ICDrv}

# --- FUNCTIONS ---------------------------------------------------------------
def flattenIt(cIterable, rmvNaN=True, retArr=False):
    itFlat = list(itertools.chain.from_iterable(cIterable))
    if rmvNaN:
        itFlat = [x for x in itFlat if (x is not None and not np.isnan(x))]
    if retArr:
        itFlat = np.array(itFlat)
    return itFlat

def getSerData(dDfrIn, sGT, sMet, sPho):
    d = dDfrIn[sGT]
    return d[(d[S_MET] == sMet) & (d[S_PHO] == sPho)].squeeze()

def getL2Val(midV=0., lenLn=1., cDir='both'):
    l = [midV - lenLn, midV + lenLn]
    if cDir == 'left':
        l = [midV - lenLn, midV]
    elif cDir == 'right':
        l = [midV, midV + lenLn]
    return l

def plotWskH(dPlt, cAx, xMid=0., yB=0., yT=0., lenW=.1, lwdW=.5, cDir='both'):
    lX = getL2Val(xMid, lenW, cDir=cDir)
    for cY in [yB, yT]:
        cAx.plot(lX, [cY]*2, lw=lwdW, color=dPlt['clrDef'])

def plotWskV(dPlt, cAx, xB=0., xT=0., yMid=0., lenW=.1, lwdW=.5, cDir='both'):
    lY = getL2Val(yMid, lenW, cDir=cDir)
    for cX in [xB, xT]:
        cAx.plot([cX]*2, lY, lw=lwdW, color=dPlt['clrDef'])

def plotClrLn(dPlt, cAx, lX, lY, cLwd, cClr, dltLwd=.5):
    cAx.plot(lX, lY, lw=cLwd, color=cClr)
    cAx.plot(lX, lY, lw=max(0, cLwd - 2*dltLwd), color=dPlt['clrDef'])

def plotMean(dPlt, cAx, xMid=0., yMn=0., cClr=S_CLR0):
    lX, lY = [xMid - dPlt['lenHfMnBar'], xMid + dPlt['lenHfMnBar']], [yMn]*2
    plotClrLn(dPlt, cAx, lX, lY, cLwd=dPlt['lwdMnBar'], cClr=cClr,
              dltLwd=dPlt['dltLwdMnBar'])
    plotWskV(dPlt, cAx, xB=lX[0], xT=lX[1], yMid=yMn, lenW=dPlt['lenWMnBar'],
             lwdW=dPlt['lwdWMnBar'])

def getNumSDSteps(othYMn, yMn, ySD):
    nSDStp, sgnDiff = 0, 0
    if othYMn > yMn:
        nSDStp, sgnDiff = math.ceil((othYMn - yMn)/ySD), 1
    elif othYMn < yMn:
        nSDStp, sgnDiff = math.ceil((yMn - othYMn)/ySD), -1
    return nSDStp, math.copysign(ySD, sgnDiff)

def plotStepsSD(dPlt, cAx, tFtI, xLn=0., othYMn=0.):
    yMn, ySD, cClr = tFtI
    nSDStp, ySDSgn = getNumSDSteps(othYMn, yMn, ySD)
    for i in range(nSDStp):
        lY = [yMn + i*ySDSgn, yMn + (i + 1)*ySDSgn]
        plotClrLn(dPlt, cAx, [xLn]*2, lY, cLwd=dPlt['lwdLnSD'], cClr=cClr,
                  dltLwd=dPlt['dltLwdLnSD'])
        plotWskH(dPlt, cAx, xMid=xLn, yB=lY[0], yT=lY[1],
                 lenW=dPlt['lenWLnSD'], lwdW=dPlt['lwdWLnSD'], cDir='right')
    lX = [xLn, dPlt['dPosFt'][0] - dPlt['lenHfMnBar']]
    cAx.plot(lX, [othYMn]*2, ls=dPlt['lstCnLnSD'], lw=dPlt['lwdCnLnSD'],
             color=dPlt['clrDef'])

def modAlpha(dClr, cAlpha=1.):
    for cK, cV in dClr.items():
        if type(cV) == tuple:
            if len(cV) in [3, 4]:
                dClr[cK] = (cV[0], cV[1], cV[2], cAlpha)

def getLClr(dClr, cSer):
    serClr = pd.Series(dClr[S_IC], index=cSer.index)
    for sI, x in cSer.items():
        if x < 0:
            serClr.at[sI] = dClr[S_IC_M]
        elif x > 0:
            serClr.at[sI] = dClr[S_IC_P]
    return serClr.to_list()

def getDfrsDatClr(dClr, dAlp, cSer):
    dfrDat = pd.DataFrame(0., index=cSer.index, columns=L_S_IC_ALL)
    modAlpha(dClr, dAlp['S'])                   # 'standard' alpha
    arrClr = np.array([[dClr[sIC] for sIC in L_S_IC_ALL]
                       for _ in range(cSer.index.size)], dtype=object)
    dfrClr = pd.DataFrame(arrClr, index=cSer.index, columns=L_S_IC_ALL)
    for k, (sI, x) in enumerate(cSer.items()):
        if x < 0:
            dfrDat.at[sI, S_IC_M] = x
        elif x > 0:
            dfrDat.at[sI, S_IC_P] = x
        if k < len(L_S_IC_ALL):
            dfrClr.at[sI, S_IC] = dClr[L_S_IC_ALL[k]]
    xICm, xICp = dfrDat.loc[:, S_IC_M].sum(), dfrDat.loc[:, S_IC_P].sum()
    if cSer.index.size > 0:
        modAlpha(dClr, dAlp['L'])               # 'low' alpha
        dfrDat.at[cSer.index[0], S_IC] = xICm
        dfrClr.at[cSer.index[0], S_IC] = dClr[L_S_IC_ALL[0]]
        if cSer.index.size > 1:
            dfrDat.at[cSer.index[1], S_IC] = xICp
            dfrClr.at[cSer.index[1], S_IC] = dClr[L_S_IC_ALL[1]]
            if cSer.index.size > 2:
                dfrDat.at[cSer.index[2], S_IC] = xICm + xICp
    return dfrDat, dfrClr

def saveClosePlot(cFig, pPltF, l=None):
    if l is not None:
        cFig.savefig(pPltF, bbox_extra_artists=(l,), bbox_inches='tight')
    else:
        cFig.savefig(pPltF)
    plt.close()

def printElapsedTimeSim(stT, cT, sPre = 'Time'):
    # calculate and display elapsed time
    elT = round(cT - stT, R04)
    print(sPre, 'elapsed:', elT, 'seconds, this is', round(elT/60, R04),
          'minutes or', round(elT/3600, R04), 'hours or',
          round(elT/(3600*24), R04), 'days.')

# --- CLASSES -----------------------------------------------------------------
class BaseClass():
    def __init__(self):
        self.idO = S_BASE_CL
        self.descO = 'Base class'
        print('Initiated "Base" base object.')

    def printIDDesc(self):
        print('Object ID:', self.idO)
        print('Object description:', self.descO)

    def printAttrList(self):
        lAttr = dir(self)
        print('List of attributes:')
        for cAttr in lAttr:
            print(cAttr)

    def printAttrData(self):
        print('Attributes and attribute values:')
        d = vars(self)
        for cK, cV in d.items():
            print(cK, ':\t', cV)

# .............................................................................
class InputData(BaseClass):
    def __init__(self, dInp):
        super().__init__()
        self.idO = dInp['sInpDat']
        self.descO = 'Input data class'
        self.dI = dInp
        self.fillInp()
        print('Initiated "InputData" base object.')

    def fillInp(self):
        for sK, cV in self.dI.items():
            setattr(self, sK, cV)
        print('Set InputData attributes.')

# .............................................................................
class RootClass(BaseClass):
    def __init__(self, InpD):
        super().__init__()
        self.idO = InpD.sRoot
        self.descO = 'Root class'
        self.inpD = InpD
        self.dDfrIn = None
        print('Initiated "RootClass" base object.')

    def printObjInfo(self):
        print('-'*20, 'Object', self.descO, '(ID', self.idO, ')', '-'*20)
        print('-'*8, 'Input data:')
        self.inpD.printAttrData()
        print('-'*8, 'Attributes of', self.descO, 'class:')
        self.printAttrData()

    def printDDfrInp(self, lIdxCol=None):
        if self.dDfrIn is None:
            print('Input DataFrames dictionary does not have any content yet.')
        else:
            if lIdxCol is None:
                print('Input DataFrames dictionary:')
                print(self.dDfrIn)
            else:
                print('Input DataFrames dictionaries (columns', lIdxCol, '):')
                for cK, dfrIn in self.dDfrIn.items():
                    print('-'*8, 'Dictionary corresponding to key', cK, ':')
                    print(dfrIn.iloc[:, lIdxCol])

    def printDfrInp(self, cK, printIdx=False, printCol=False):
        if cK in self.dDfrIn:
            print('DataFrame corresponding to key ' + str(cK) + ':')
            print(self.dDfrIn[cK])
            if printIdx:
                print('DataFrame index (with length ' +
                      str(self.dDfrIn[cK].index.size) + '):')
                print(self.dDfrIn[cK].index)
            if printCol:
                print('DataFrame columns (with length ' +
                      str(self.dDfrIn[cK].columns.size) + '):')
                print('DataFrame columns:')
                print(self.dDfrIn[cK].columns)
        else:
            print('Key', cK, 'not in input DataFrame dictionary!')
            print('Keys of input DataFrame dictionary:', list(self.dDfrIn))

# .............................................................................
class Plotter(RootClass):
    def __init__(self, InpD):
        super().__init__(InpD)
        self.idO = InpD.sPltr
        self.descO = 'Class for plotting'
        self.sSp = self.inpD.sSep
        self.dPlt = self.inpD.plot_Gen
        self.dPPltF, self.dPPltF1, self.dPPltF2 = {}, {}, {}
        self.setDummyVal()
        print('Initiated "Plotter" base object.')

    def printDPPltF(self):
        for t in [(self.dPPltF, ''), (self.dPPltF1, ' 1'),
                  (self.dPPltF2, ' 2')]:
            if len(t[0]) > 0:
                print('Dictionary of plot file' + t[1] + ' paths:')
                for tK, pF in t[0].items():
                    print(tK, ':', pF)
                print('-'*64)

    def loadDDfrInp(self, idxC=None):
        # load input DataFrames, and save them in dictionary
        if hasattr(self, 'dPFIn'):
            self.dDfrIn = {}
            for sK in self.dPFIn:
                self.dDfrIn[sK] = pd.read_csv(self.dPFIn[sK], sep=self.sSp,
                                              index_col=idxC)

    def getDPPltF(self, dPlt, sFOut):
        dPPltF, mxL = {}, self.inpD.maxSLen
        for sID, t in dPlt['dTupIn'].items():
            u = S_USC.join([sFOut, sID])
            lCmpNmF = [s.replace(S_BAR, S_USC)[:mxL] for s in t]
            sPltF = S_DOT.join([S_USC.join([u] + lCmpNmF), S_PDF])
            dPPltF[sID] = (t, os.path.join(self.pDOut, sPltF))
        return dPPltF

    def setDummyVal(self):
        l = ['clrDef', 'szFontLeg', 'nCharDsp', 'posLegXY', 'dTupIn',
             'nmPlt', 'tFigSz', 'symMark', 'szMark', 'lenMnBar', 'lenHfMnBar',
             'lwdMnBar', 'dltLwdMnBar', 'lenWMnBar', 'lwdWMnBar', 'lwdLnSD',
             'dltLwdLnSD', 'lenWLnSD', 'lwdWLnSD', 'lstCnLnSD', 'lwdCnLnSD',
             'dPosFt', 'axXLim', 'axYLim', 'sXLbl', 'sYLbl', 'lblXTck',
             'stepYTck', 'adaptYLim', 'plotTtl', 'dFontTtl', 'locTtl',
             'padTtl', 'degRotXLbl', 'lwdPlt', 'lClr', 'locLegend', 'wdthGrp',
             'wdthBar', 'plotVLines', 'lwdVLine', 'clrVLine', 'dClr',
             'alphaClr', 'alphaClrLow', 'lwdBarIC', 'hatBarIC',
             'axXTck', 'axYTck']
        for s in l:
            if s not in self.dPlt:
                self.dPlt[s] = None

    def iniPlot(self, dPlt):
        cFig, cAx = plt.subplots()
        if dPlt['tFigSz'] is not None and len(dPlt['tFigSz']) == 2:
            cFig.set_size_inches(dPlt['tFigSz'])
        return cFig, cAx

    def decoratePlot(self, dPlt, cAx, sTtl=None):
        if dPlt['axXLim'] is not None:
            cAx.set_xlim(dPlt['axXLim'])
        if dPlt['axYLim'] is not None:
            cAx.set_ylim(dPlt['axYLim'])
        if dPlt['axXTck'] is not None:
            cAx.set_xticks(dPlt['axXTck'])
            if dPlt['lblXTck'] is not None:
                cAx.set_xticklabels(dPlt['lblXTck'])
            if dPlt['degRotXLbl'] is not None:
                for cXLbl in cAx.get_xticklabels():
                    cXLbl.set_rotation(dPlt['degRotXLbl'])
        if dPlt['axYTck'] is not None:
            cAx.set_yticks(dPlt['axYTck'])
        if sTtl is not None and dPlt['plotTtl']:
            cAx.set_title(sTtl, fontdict=dPlt['dFontTtl'], loc=dPlt['locTtl'],
                          pad=dPlt['padTtl'])
        if dPlt['sXLbl'] is not None:
            cAx.set_xlabel(dPlt['sXLbl'])
        if dPlt['sYLbl'] is not None:
            cAx.set_ylabel(dPlt['sYLbl'])
        plt.tight_layout()

# .............................................................................
class DevSDPlotter(Plotter):
    def __init__(self, InpD):
        super().__init__(InpD)
        self.idO = self.inpD.sPltr_DevSD
        self.descO = 'Deviations (as multiples of feature SD) plotter'
        self.dPFIn = self.inpD.dPFIn_DevSD
        self.pDOut = self.inpD.pOutPDF
        self.dPlt = self.dPlt | self.inpD.plot_DevSD
        self.dPPltF = self.getDPPltF(self.dPlt, sFOut=self.inpD.sFOut_DevSD)
        self.loadDDfrInp(idxC=0)
        print('Initiated "DevSDPlotter" base object and loaded input data.')

    def preProcData(self, sGT, sMP, sFtChg):
        d, dFt = self.dDfrIn[(S_METS, sGT)], {}
        for k, cFt in enumerate(sFtChg.split(S_BAR)):
            lSC = [s for s in d.columns if (s[-1] in L_S_0_6 and
                                            s.startswith(cFt))]
            cSer, cC = d.loc[sMP, lSC], S_C + str(k)
            dFt[k] = (cFt, lSC, cSer, np.mean(cSer), np.std(cSer, ddof=1), cC)
        assert len(dFt) == 2
        return dFt

    def setTicks(self, dFt=None):
        if self.dPlt['axXLim'] is not None and len(self.dPlt['axXLim']) >= 2:
            self.dPlt['axXTck'] = [self.dPlt['dPosFt'][k] for k in dFt]
            self.dPlt['lblXTck'] = [dFt[k][0] for k in dFt]
        if self.dPlt['adaptYLim']:
            lAllV = flattenIt([t[2] for t in dFt.values()])
            mn, mx = math.floor(min(lAllV)), math.ceil(max(lAllV))
            self.dPlt['axYLim'] = (mn, mx)
            self.dPlt['axYTck'] = np.arange(mn, mx + 1, self.dPlt['stepYTck'])

    def createPlot(self, dFtI):
        self.setTicks(dFt=dFtI)
        cFig, cAx = self.iniPlot(self.dPlt)
        for k, (cFt, lHdC, cSer, cMn, cSD, cClr) in dFtI.items():
            lX = [self.dPlt['dPosFt'][k]]*cSer.size
            cAx.scatter(lX, cSer, marker=self.dPlt['symMark'],
                        s=self.dPlt['szMark'])
            plotMean(self.dPlt, cAx, lX[0], cMn, cClr=cClr)
            if k == 1 and cSD > 0:
                xPLn = lX[0] + self.dPlt['lenHfMnBar']
                plotStepsSD(self.dPlt, cAx, tFtI=dFtI[k][3:], xLn=xPLn,
                            othYMn=dFtI[0][3])
        return cFig, cAx

    def decoratePlot(self, cAx, sGT, sMP=''):
        # sTtl = sMP + '\n(' + D_NM_GT[sGT] + ')'
        sTtl = sMP
        super().decoratePlot(self.dPlt, cAx, sTtl=sTtl)

    def plotDevSD(self):
        for sID, ((sGT, sMP, sFtChg), pPltF) in self.dPPltF.items():
            print('Plotting deviation in multiples of SD: "' + str(sID) +
                  '" for "' + sGT + '", substance ' + sMP + ', feature change',
                  sFtChg, '...')
            # if not os.path.isfile(pPltF):
            cFig, cAx = self.createPlot(self.preProcData(sGT, sMP, sFtChg))
            self.decoratePlot(cAx, sGT=sGT, sMP=sMP)
            saveClosePlot(cFig, pPltF)

# .............................................................................
class ICDerivPlotter(Plotter):
    def __init__(self, InpD):
        super().__init__(InpD)
        self.idO = self.inpD.sPltr_ICDrv
        self.descO = 'Concordance index derivation plotter'
        self.dPFIn = self.inpD.dPFIn_ICDrv
        self.pDOut = self.inpD.pOutPDF
        self.dPlt = self.dPlt | self.inpD.plot_ICDrv
        self.dPPltF = self.getDPPltF(self.dPlt, sFOut=self.inpD.sFOut_ICDrv)
        self.loadDDfrInp()
        print('Initiated "ICDerivPlotter" base object and loaded input data.')

    def setTicks(self):
        axYLim = self.dPlt['axYLim']
        axXTck, axYTck = np.arange(len(self.dPlt['lblXTck'])), None
        if axYLim is not None and len(axYLim) >= 2:
            axYTck = np.arange(-(-axYLim[0]//2*2), axYLim[1]//2*2 + 1,
                               self.dPlt['stepYTck'])
        self.dPlt['axXTck'], self.dPlt['axYTck'] = axXTck, axYTck

    def createPlot(self, cSer, lILegPlt):
        self.setTicks()
        nChD, xTck = self.dPlt['nCharDsp'], self.dPlt['axXTck']
        cFig, cAx = self.iniPlot(self.dPlt)
        for k, (sK, lSHd) in enumerate(D_S_FT_CHG.items()):
            cSerGrp = cSer.loc[lSHd]
            cSerGrp.index = L_S_FT_CHG
            xLoc = (xTck + (2*k + 1)/(2*len(D_S_FT_CHG))*self.dPlt['wdthGrp'] -
                    1/2 + self.dPlt['wdthBar']/2)
            cAx.bar(xLoc, height=cSerGrp, width=self.dPlt['wdthBar'],
                    lw=self.dPlt['lwdPlt'], color=self.dPlt['lClr'][k],
                    label=lILegPlt[k][:nChD])
        cAx.plot([-1/2, cSerGrp.index.size - 1/2], [0, 0],
                 lw=self.dPlt['lwdPlt'], color=self.dPlt['clrDef'])
        return cFig, cAx

    def decoratePlot(self, cAx, sGT):
        super().decoratePlot(self.dPlt, cAx, sTtl=D_NM_GT[sGT])
        l = cAx.legend(loc=self.dPlt['locLegend'],
                       bbox_to_anchor=self.dPlt['posLegXY'],
                       fontsize=self.dPlt['szFontLeg'])
        if self.dPlt['plotVLines']:
            yL, yU = self.dPlt['axYTck'][0], self.dPlt['axYTck'][-1]
            plt.vlines(np.arange(-1/2, len(self.dPlt['axXTck']) + 1/2), yL, yU,
                       lw=self.dPlt['lwdVLine'], colors=self.dPlt['clrVLine'])
        return l

    def plotICDeriv(self):
        for sID, ((sGT, sMet, sPho), pPltF) in self.dPPltF.items():
            print('Plotting IC derivation', sID, 'for "' + sGT +
                  '" and the pair (' + sMet + ', ' + sPho + ')...')
            cSer = getSerData(self.dDfrIn, sGT, sMet, sPho)
            lILegPlt = [sMet, sPho, S_CMB_DEV, S_IC_CMP]
            # if not os.path.isfile(pPltF):
            cFig, cAx = self.createPlot(cSer, lILegPlt)
            cLeg = self.decoratePlot(cAx, sGT=sGT)
            saveClosePlot(cFig, pPltF, cLeg)

# .............................................................................
class CmbToICPlotter(Plotter):
    def __init__(self, InpD, pltrICDrv):
        super().__init__(InpD)
        self.idO = self.inpD.sPltr_CmbToIC
        self.descO = 'Concordance index derivation plotter'
        self.pDOut = self.inpD.pOutPDF
        self.dPltB = self.dPlt | self.inpD.plotB_CmbToIC
        self.dPlt1 = self.dPltB | self.inpD.plot1_CmbToIC
        self.dPlt2 = self.dPltB | self.inpD.plot2_CmbToIC
        self.dPPltF1 = self.getDPPltF(self.dPlt1, self.inpD.sF1Out_CmbToIC)
        self.dPPltF2 = self.getDPPltF(self.dPlt2, self.inpD.sF2Out_CmbToIC)
        self.dDfrIn = pltrICDrv.dDfrIn
        print('Initiated "CmbToICPlotter" base object and obtained data.')

    def setTicks(self, dPlt, serD=None):
        if serD is not None:                # adapt y-limits
            self.dPlt['axYLim'] = (math.floor(min(serD)), math.ceil(max(serD)))
        axYLim, stpYTck = dPlt['axYLim'], dPlt['stepYTck']
        axXTck, axYTck = np.arange(len(dPlt['lblXTck'])), None
        if axYLim is not None and len(axYLim) >= 2 and stpYTck is not None:
            axYTck = np.arange(-(-axYLim[0]//2*2), axYLim[1]//2*2 + 1, stpYTck)
        dPlt['axXTck'], dPlt['axYTck'] = axXTck, axYTck

    def createPlot1(self, cSer):
        assert cSer.size == len(L_S_FT_CHG)
        cSer.index = L_S_FT_CHG
        self.setTicks(self.dPlt1)
        modAlpha(self.dPltB['dClr'], self.dPlt1['alphaClr'])
        cFig, cAx = self.iniPlot(self.dPlt1)
        cAx.bar(self.dPlt1['axXTck'], height=cSer, width=self.dPlt1['wdthBar'],
                lw=self.dPltB['lwdPlt'], ec=self.dPlt['clrDef'],
                label=S_IC_CMP, color=getLClr(self.dPltB['dClr'], cSer))
        cAx.plot([-1/2, cSer.index.size - 1/2], [0, 0],
                 lw=self.dPltB['lwdPlt'], color=self.dPlt['clrDef'])
        return cFig, cAx

    def decoratePlot1(self, cAx, sGT):
        super().decoratePlot(self.dPlt1, cAx, sTtl=D_NM_GT[sGT])
        if self.dPlt1['plotVLines']:
            yL, yU = self.dPlt1['axYTck'][0], self.dPlt1['axYTck'][-1]
            plt.vlines(np.arange(-1/2, len(self.dPlt1['axXTck']) + 1/2),
                       yL, yU, lw=self.dPlt1['lwdVLine'],
                       colors=self.dPlt1['clrVLine'])

    def createPlot2(self, cSer):
        assert cSer.size == len(L_S_FT_CHG)
        cSer.index = L_S_FT_CHG
        dAlpha = {'S': self.dPlt2['alphaClr'], 'L': self.dPlt2['alphaClrLow']}
        dfrDat, dfrClr = getDfrsDatClr(self.dPltB['dClr'], dAlpha, cSer)
        baseSerD = pd.Series([0.]*dfrDat.columns.size, index=dfrDat.columns)
        self.setTicks(self.dPlt2, dfrDat.iloc[:, 2])
        cFig, cAx = self.iniPlot(self.dPlt2)
        for cIdx, pltSerD in dfrDat.iterrows():
            cAx.bar(self.dPlt2['axXTck'], height=pltSerD,
                    width=self.dPlt2['wdthBar'], bottom=baseSerD,
                    lw=self.dPltB['lwdPlt'], ec=self.dPlt['clrDef'],
                    label=dfrDat.columns, color=dfrClr.loc[cIdx, :].to_list())
            baseSerD = baseSerD.add(pltSerD)
            baseSerD.at[S_IC] = 0.
        # plot wide black edges of the final IC bar (index of IC = 2)
        yIC, clrIC = dfrDat.iloc[2, :].at[S_IC], dfrClr.iloc[2, :].at[S_IC]
        cAx.bar(self.dPlt2['axXTck'][-1], height=yIC,
                width=self.dPlt2['wdthBar'],  lw=self.dPlt2['lwdBarIC'],
                ec=self.dPlt['clrDef'], fc=clrIC, hatch=self.dPlt2['hatBarIC'])
        cAx.plot([-1/2, baseSerD.index.size - 1/2], [0, 0],
                 lw=self.dPltB['lwdPlt'], color=self.dPlt['clrDef'])
        if self.dPlt['axYLim'] is not None and len(self.dPlt['axYLim']) == 2:
            plt.ylim(self.dPlt['axYLim'])
        return cFig, cAx

    def decoratePlot2(self, cAx, sGT):
        super().decoratePlot(self.dPlt2, cAx, sTtl=D_NM_GT[sGT])
        if self.dPlt2['plotVLines']:
            yL, yU = cAx.get_ylim()
            plt.vlines(np.arange(-1/2, len(self.dPlt2['axXTck']) + 1/2),
                       yL, yU, lw=self.dPlt2['lwdVLine'],
                       colors=self.dPlt2['clrVLine'])

    def plotCmbToIC(self):
        # plot 1 - IC component plot
        for sID, ((sGT, sMet, sPho), pPltF) in self.dPPltF1.items():
            print('Plotting combining components to IC', sID, 'for "' + sGT +
                  '" and the pair (' + sMet + ', ' + sPho + ')...')
            cSer = getSerData(self.dDfrIn, sGT, sMet, sPho)
            # if not os.path.isfile(pPltF):
            cFig, cAx = self.createPlot1(cSer.loc[D_S_FT_CHG[S_C]])
            self.decoratePlot1(cAx, sGT=sGT)
            saveClosePlot(cFig, pPltF)
        # plot 2 - combine to IC(-), IC(+), IC plot
        for sID, ((sGT, sMet, sPho), pPltF) in self.dPPltF2.items():
            cSer = getSerData(self.dDfrIn, sGT, sMet, sPho)
            cFig, cAx = self.createPlot2(cSer.loc[D_S_FT_CHG[S_C]])
            self.decoratePlot2(cAx, sGT=sGT)
            saveClosePlot(cFig, pPltF)

# --- MAIN --------------------------------------------------------------------
startTime = time.time()
print('+'*25 + ' START', time.ctime(startTime), '+'*25)

inpDat = InputData(dInput)
if inpDat.doPlot_DevSD:
    cPltrDevSD = DevSDPlotter(inpDat)
    cPltrDevSD.printIDDesc()
    # cPltrDevSD.printDPPltF()
    # cPltrDevSD.printAttrData()
    # cPltrDevSD.printObjInfo()
    cPltrDevSD.plotDevSD()
if inpDat.doPlot_ICDrv or inpDat.doPlot_CmbToIC:
    cPltrICDrv = ICDerivPlotter(inpDat)
if inpDat.doPlot_ICDrv:
    cPltrICDrv.printIDDesc()
    cPltrICDrv.plotICDeriv()
if inpDat.doPlot_CmbToIC:
    cPltrCmbToIC = CmbToICPlotter(inpDat, cPltrICDrv)
    cPltrCmbToIC.printIDDesc()
    # cPltrCmbToIC.printDPPltF()
    cPltrCmbToIC.plotCmbToIC()

print('-'*80)
printElapsedTimeSim(startTime, time.time(), 'Total time')
print('+'*37, 'DONE.', '+'*36)

###############################################################################

# -*- coding: utf-8 -*-
###############################################################################
# --- PlotNetwork.py ----------------------------------------------------------
###############################################################################
import os, itertools

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import networkx as nx

# === CONSTANTS ===============================================================
S_EXT_CSV = 'csv'
S_EXT_PDF = 'pdf'

S_SEL = 'sel'
S_THR = 'thr'
S_MIN = 'min'
S_MAX = 'max'
S_AVG = 'avg'
S_INV = 'inv'
S_RGB = 'rgb'

S_MET_D = 'MetD'
S_PHO_D = 'PhoD'

S_PROT = 'Protein'
S_BIN_C_G = 'BinCode'
S_BIN_C_3 = 'BinCode3'
S_BIN_C_2 = 'BinCode2'
S_BIN_C_1 = 'BinCode1'
L_S_PHO_CL = [S_PROT, S_BIN_C_G, S_BIN_C_3, S_BIN_C_2, S_BIN_C_1]

S_IDX_PDCSD = 'IdxPosDvClSD'

S_SRT_FLT = 'SortFilt'

MAX_CLR_VAL = 255
INDEX_DATA_EDGE_ATTR = 2
ATTR_1 = 1
ATTR_2 = 2

THR_2 = 2
THR_3 = 3
THR_4 = 4
THR_5 = 5
THR_6 = 6
THR_7 = 7
THR_8 = 8
THR_9 = 9

N_DIG_RND_04 = 4

# === COLOURS =================================================================
clrR = 'rgb(255, 0, 0)'
clrG = 'rgb(0, 255, 0)'
clrB = 'rgb(0, 0, 255)'
clrY = 'rgb(230, 179, 0)'
clrDk1P = 'rgb(204, 0, 102)'
clrDk1G ='rgb(0, 191, 0)'
clrDk2G ='rgb(0, 127, 0)'
clrGr05 = 'rgb(127, 127, 127)'
clr8 = '#888'

# === INPUT ===================================================================
# data file type and separators -----------------------------------------------
sSelBC2 = 'SelBin2sF24'                 # 'SelBin2sA14'
                                        # 'SelBin2sB22'
                                        # 'SelBin2sC21'
                                        # 'SelBin2sD22'
                                        # 'SelBin2sE23'
                                        # 'SelBin2sF24'
sTask = 'IPDC_ge2_Avg_GT5'              # 'Idx_gr8_SD2_gr8_GT0'
                                        # 'Idx_gr8_SD2_gr8_GT1'
                                        # 'Idx_gr8_SD2_gr8_GT2'
                                        # 'Idx_gr8_SD2_gr8_GT3'
                                        # 'Idx_gr8_SD2_gr8_GT4'
                                        # 'Idx_gr8_SD2_gr8_GT5'
                                        # 'SRR_AllGT_Idx_gr8'
                                        # 'IPDC_ge$X$_$YYY$_GT$Z$'
                                        # $X$ in [2, 3, 4, 5, 6, 7, 8, 9]
                                        # $YYY$ in [Avg, Max]
                                        # $Z$ in [0, 1, 5]
cSep = ';'

# derived values - tasks ------------------------------------------------------
sSelNmF = sSelBC2 + '_' + sTask
# GT0
lTaskAvgGT0 = ['IPDC_ge2_Avg_GT0', 'IPDC_ge3_Avg_GT0',
               'IPDC_ge4_Avg_GT0', 'IPDC_ge5_Avg_GT0',
               'IPDC_ge6_Avg_GT0', 'IPDC_ge7_Avg_GT0',
               'IPDC_ge8_Avg_GT0', 'IPDC_ge9_Avg_GT0']
lTaskMaxGT0 = ['IPDC_ge2_Max_GT0', 'IPDC_ge3_Max_GT0',
               'IPDC_ge4_Max_GT0', 'IPDC_ge5_Max_GT0',
               'IPDC_ge6_Max_GT0', 'IPDC_ge7_Max_GT0',
               'IPDC_ge8_Max_GT0', 'IPDC_ge9_Max_GT0']
lTaskAvgMaxGT0 = lTaskAvgGT0 + lTaskMaxGT0
lTaskGE2GT0 = ['IPDC_ge2_Avg_GT0', 'IPDC_ge2_Max_GT0']
lTaskGE3GT0 = ['IPDC_ge3_Avg_GT0', 'IPDC_ge3_Max_GT0']
lTaskGE4GT0 = ['IPDC_ge4_Avg_GT0', 'IPDC_ge4_Max_GT0']
lTaskGE5GT0 = ['IPDC_ge5_Avg_GT0', 'IPDC_ge5_Max_GT0']
lTaskGE6GT0 = ['IPDC_ge6_Avg_GT0', 'IPDC_ge6_Max_GT0']
lTaskGE7GT0 = ['IPDC_ge7_Avg_GT0', 'IPDC_ge7_Max_GT0']
lTaskGE8GT0 = ['IPDC_ge8_Avg_GT0', 'IPDC_ge8_Max_GT0']
lTaskGE9GT0 = ['IPDC_ge9_Avg_GT0', 'IPDC_ge9_Max_GT0']

# GT1
lTaskAvgGT1 = ['IPDC_ge2_Avg_GT1', 'IPDC_ge3_Avg_GT1',
               'IPDC_ge4_Avg_GT1', 'IPDC_ge5_Avg_GT1',
               'IPDC_ge6_Avg_GT1', 'IPDC_ge7_Avg_GT1',
               'IPDC_ge8_Avg_GT1', 'IPDC_ge9_Avg_GT1']
lTaskMaxGT1 = ['IPDC_ge2_Max_GT1', 'IPDC_ge3_Max_GT1',
               'IPDC_ge4_Max_GT1', 'IPDC_ge5_Max_GT1',
               'IPDC_ge6_Max_GT1', 'IPDC_ge7_Max_GT1',
               'IPDC_ge8_Max_GT1', 'IPDC_ge9_Max_GT1']
lTaskAvgMaxGT1 = lTaskAvgGT1 + lTaskMaxGT1
lTaskGE2GT1 = ['IPDC_ge2_Avg_GT1', 'IPDC_ge2_Max_GT1']
lTaskGE3GT1 = ['IPDC_ge3_Avg_GT1', 'IPDC_ge3_Max_GT1']
lTaskGE4GT1 = ['IPDC_ge4_Avg_GT1', 'IPDC_ge4_Max_GT1']
lTaskGE5GT1 = ['IPDC_ge5_Avg_GT1', 'IPDC_ge5_Max_GT1']
lTaskGE6GT1 = ['IPDC_ge6_Avg_GT1', 'IPDC_ge6_Max_GT1']
lTaskGE7GT1 = ['IPDC_ge7_Avg_GT1', 'IPDC_ge7_Max_GT1']
lTaskGE8GT1 = ['IPDC_ge8_Avg_GT1', 'IPDC_ge8_Max_GT1']
lTaskGE9GT1 = ['IPDC_ge9_Avg_GT1', 'IPDC_ge9_Max_GT1']

# GT5
lTaskAvgGT5 = ['IPDC_ge2_Avg_GT5', 'IPDC_ge3_Avg_GT5',
               'IPDC_ge4_Avg_GT5', 'IPDC_ge5_Avg_GT5',
               'IPDC_ge6_Avg_GT5', 'IPDC_ge7_Avg_GT5',
               'IPDC_ge8_Avg_GT5', 'IPDC_ge9_Avg_GT5']
lTaskMaxGT5 = ['IPDC_ge2_Max_GT5', 'IPDC_ge3_Max_GT5',
               'IPDC_ge4_Max_GT5', 'IPDC_ge5_Max_GT5',
               'IPDC_ge6_Max_GT5', 'IPDC_ge7_Max_GT5',
               'IPDC_ge8_Max_GT5', 'IPDC_ge9_Max_GT5']
lTaskAvgMaxGT5 = lTaskAvgGT5 + lTaskMaxGT5
lTaskGE2GT5 = ['IPDC_ge2_Avg_GT5', 'IPDC_ge2_Max_GT5']
lTaskGE3GT5 = ['IPDC_ge3_Avg_GT5', 'IPDC_ge3_Max_GT5']
lTaskGE4GT5 = ['IPDC_ge4_Avg_GT5', 'IPDC_ge4_Max_GT5']
lTaskGE5GT5 = ['IPDC_ge5_Avg_GT5', 'IPDC_ge5_Max_GT5']
lTaskGE6GT5 = ['IPDC_ge6_Avg_GT5', 'IPDC_ge6_Max_GT5']
lTaskGE7GT5 = ['IPDC_ge7_Avg_GT5', 'IPDC_ge7_Max_GT5']
lTaskGE8GT5 = ['IPDC_ge8_Avg_GT5', 'IPDC_ge8_Max_GT5']
lTaskGE9GT5 = ['IPDC_ge9_Avg_GT5', 'IPDC_ge9_Max_GT5']

# Bin2sA14 / All GT
lTaskAvg = lTaskAvgGT0 + lTaskAvgGT1 + lTaskAvgGT5
lTaskMax = lTaskMaxGT0 + lTaskMaxGT1 + lTaskMaxGT5
lTaskAvgMax = lTaskAvgMaxGT0 + lTaskAvgMaxGT1 + lTaskAvgMaxGT5
lTaskGE2 = lTaskGE2GT0 + lTaskGE2GT1 + lTaskGE2GT5
lTaskGE3 = lTaskGE3GT0 + lTaskGE3GT1 + lTaskGE3GT5
lTaskGE4 = lTaskGE4GT0 + lTaskGE4GT1 + lTaskGE4GT5
lTaskGE5 = lTaskGE5GT0 + lTaskGE5GT1 + lTaskGE5GT5
lTaskGE6 = lTaskGE6GT0 + lTaskGE6GT1 + lTaskGE6GT5
lTaskGE7 = lTaskGE7GT0 + lTaskGE7GT1 + lTaskGE7GT5
lTaskGE8 = lTaskGE8GT0 + lTaskGE8GT1 + lTaskGE8GT5
lTaskGE9 = lTaskGE9GT0 + lTaskGE9GT1 + lTaskGE9GT5

# names and paths -------------------------------------------------------------
pRelDatF_G = os.path.join('..', '..', '..', '12_SysBio02_DataAnalysis',
                          '11_ResultCSV_GT015', '21_R_81_BinaryOps')
pRelDatF_92 = os.path.join('..', '..', '..', '12_SysBio02_DataAnalysis',
                           '92_Networkx', '01_Data')
pRelPltF = os.path.join('..', '..', '..', '12_SysBio02_DataAnalysis',
                        '92_Networkx', '09_Plots')
nmDatF = 'Corr__BinOp_MetD_GT0_AllD_PhoD_GT0_AllD__IdxPosDvClSD_gr8_SD2_gr8'
nmPltF = 'NetworkPlot_Idx_gr8_SD2_gr8_GT0'
pRelDatF = pRelDatF_G
if sTask in ['Idx_gr8_SD2_gr8_GT0', 'Idx_gr8_SD2_gr8_GT1',
             'Idx_gr8_SD2_gr8_GT2', 'Idx_gr8_SD2_gr8_GT3',
             'Idx_gr8_SD2_gr8_GT4', 'Idx_gr8_SD2_gr8_GT5',
             'SRR_AllGT_Idx_gr8']:
    pRelDatF = pRelDatF_92
    nmPltF = 'NetworkPlot' + '_' + sSelNmF
    if sTask == 'Idx_gr8_SD2_gr8_GT0':
        nmDatF = 'Corr__BinOp_MetD_GT0_AllD_PhoD_GT0_AllD__IdxPosDvClSD_gr8_SD2_gr8'
    elif sTask == 'Idx_gr8_SD2_gr8_GT1':
        nmDatF = 'Corr__BinOp_MetD_GT1_AllD_PhoD_GT1_AllD__IdxPosDvClSD_gr8_SD2_gr8'
    elif sTask == 'Idx_gr8_SD2_gr8_GT2':
        nmDatF = 'Corr__BinOp_MetD_GT2_AllD_PhoD_GT2_AllD__IdxPosDvClSD_gr8_SD2_gr8'
    elif sTask == 'Idx_gr8_SD2_gr8_GT3':
        nmDatF = 'Corr__BinOp_MetD_GT3_AllD_PhoD_GT3_AllD__IdxPosDvClSD_gr8_SD2_gr8'
    elif sTask == 'Idx_gr8_SD2_gr8_GT4':
        nmDatF = 'Corr__BinOp_MetD_GT4_AllD_PhoD_GT4_AllD__IdxPosDvClSD_gr8_SD2_gr8'
    elif sTask == 'Idx_gr8_SD2_gr8_GT5':
        nmDatF = 'Corr__BinOp_MetD_GT5_AllD_PhoD_GT5_AllD__IdxPosDvClSD_gr8_SD2_gr8'
    elif sTask == 'SRR_AllGT_Idx_gr8':
        nmDatF = 'StronglyReducedResults_CorrDevClust_AllGT_NoMalt__SortIdxPosDvClSD'
        nmPltF = 'NetworkPlot_SRR_Top459_Idx_gr8_AllGT'
elif sTask in lTaskAvgMax:
    pRelDatF = pRelDatF_G
    nmPltF = 'NetworkPlot' + '__' + sSelNmF
    if sTask in lTaskAvgMaxGT0:
        nmDatF = 'Corr__BinOp_MetD_DvSD_GT0_AllD_PhoD_DvSD_GT0_AllD'
    elif sTask in lTaskAvgMaxGT1:
        nmDatF = 'Corr__BinOp_MetD_DvSD_GT1_AllD_PhoD_DvSD_GT1_AllD'
    elif sTask in lTaskAvgMaxGT5:
        nmDatF = 'Corr__BinOp_MetD_DvSD_GT5_AllD_PhoD_DvSD_GT5_AllD'

# strings ---------------------------------------------------------------------
sColAttr1 = S_MET_D      # attribute 1 is assigned to the metabolites
# sColAttr2 = S_PHO_D      # attribute 2 is assigned to the phosphopeptides
sColAttr2 = S_BIN_C_2    # attribute 2 is assigned to BinCode2

sNumAttr1 = S_IDX_PDCSD

thrN1 = THR_2
selOp = S_AVG

if sTask in ['Idx_gr8_SD2_gr8_GT0', 'Idx_gr8_SD2_gr8_GT1',
             'Idx_gr8_SD2_gr8_GT2', 'Idx_gr8_SD2_gr8_GT3',
             'Idx_gr8_SD2_gr8_GT4', 'Idx_gr8_SD2_gr8_GT5']:
    sNumAttr1 = 'idxPosDvClSD'
elif sTask in ['SRR_AllGT_Idx_gr8']:
    sNumAttr1 = 'idxPosDvClSD_AllGT'
elif sTask in lTaskAvgMax:
    sNumAttr1 = S_IDX_PDCSD
    if sTask in lTaskGE2:
        thrN1 = THR_2
    elif sTask in lTaskGE3:
        thrN1 = THR_3
    elif sTask in lTaskGE4:
        thrN1 = THR_4
    elif sTask in lTaskGE5:
        thrN1 = THR_5
    elif sTask in lTaskGE6:
        thrN1 = THR_6
    elif sTask in lTaskGE7:
        thrN1 = THR_7
    elif sTask in lTaskGE8:
        thrN1 = THR_8
    elif sTask in lTaskGE9:
        thrN1 = THR_9
    if sTask in lTaskAvg:
        selOp = S_AVG
    elif sTask in lTaskMax:
        selOp = S_MAX

# sorting and filter dictionaries ---------------------------------------------
dSort_G = None
dSort_MB2ND = {sColAttr1: True,   # column header (string 1): ascending (bool)
               sColAttr2: True,   # column header (string 2): ascending (bool)
               sNumAttr1: False}  # column header (num. 1): descending (bool)
dFilt_G = None
dFilt_BC2A = {S_BIN_C_2: ['1.1', '1.3', '5.3', '10.5', '12.1', '12.2',
                          '15.1', '30.11', '30.8', '31.1', '33.1', '34.19',
                          '34.5', '34.7']}
dFilt_BC2B = {S_BIN_C_2: ['1.1', '1.3', '5.3', '10.5', '12.1', '12.2',
                          '15.1', '30.11', '30.8', '31.1', '33.1', '34.19',
                          '34.5', '34.7', '30.3', '23.1', '4.1', '30.1',
                          '28.99', '2.2', '10.1', '29.6']}
dFilt_BC2C = {S_BIN_C_2: ['1.1', '1.3', '4.1', '5.3', '8.2', '10.5', '12.1',
                          '12.2', '20.1', '20.2', '29.2', '29.4', '29.5',
                          '30.1', '30.11', '30.2', '30.3', '33.99', '34.1',
                          '34.19', '34.5']}
dFilt_BC2D = {S_BIN_C_2: ['1.1', '1.3', '4.1', '5.3', '8.2', '10.5', '12.1',
                          '12.2', '20.1', '20.2', '29.2', '29.4', '29.5',
                          '30.1', '30.11', '30.2', '30.3', '31.1', '33.99',
                          '34.1', '34.19', '34.5']}
dFilt_BC2E = {S_BIN_C_2: ['1.1', '1.3', '2.2', '4.1', '5.3', '8.2', '10.5',
                          '12.1', '12.2', '20.1', '20.2', '29.2', '29.4',
                          '29.5', '30.1', '30.11', '30.2', '30.3', '31.1',
                          '33.99', '34.1', '34.19', '34.5']}
dFilt_BC2F = {S_BIN_C_2: ['1.1', '1.3', '2.2', '4.1', '5.3', '8.2', '10.5',
                          '12.1', '12.2', '20.1', '20.2', '29.2', '29.4',
                          '29.5', '30.1', '30.11', '30.2', '30.3', '31.1',
                          '33.99', '34.1', '34.19', '34.5', '34.7']}
dFilt_BC2 = dFilt_BC2A
if sSelBC2 == 'SelBin2sA14':
    dFilt_BC2 = dFilt_BC2A
elif sSelBC2 == 'SelBin2sB22':
    dFilt_BC2 = dFilt_BC2B
elif sSelBC2 == 'SelBin2sC21':
    dFilt_BC2 = dFilt_BC2C
elif sSelBC2 == 'SelBin2sD22':
    dFilt_BC2 = dFilt_BC2D
elif sSelBC2 == 'SelBin2sE23':
    dFilt_BC2 = dFilt_BC2E
elif sSelBC2 == 'SelBin2sF24':
    dFilt_BC2 = dFilt_BC2F
dThrN_IPDC = {sNumAttr1: ('>=', thrN1)}
dFilt_MBC2IPDC = {S_SEL: dFilt_BC2,
                  S_THR: dThrN_IPDC,
                  selOp: (sNumAttr1, (sColAttr1, sColAttr2))}
 
# attributes ------------------------------------------------------------------
# attrMet, attrPho, multSzAttr = 2, 1, 6

# edge input ------------------------------------------------------------------
edgeWidth = 0.5
hoverInfEdge = 'none'
modeEdge = 'lines'
dOffsClr = {S_MIN: {'R': 0.4, 'G': 0., 'B': 0.2},
            S_MAX: {'R': 0.8, 'G': 0.8, 'B': 0.8},
            S_INV: {'R': True, 'G': True, 'B': True}}
if sTask in ['Idx_gr8_SD2_gr8_GT0', 'Idx_gr8_SD2_gr8_GT1',
             'Idx_gr8_SD2_gr8_GT2', 'Idx_gr8_SD2_gr8_GT3',
             'Idx_gr8_SD2_gr8_GT4', 'Idx_gr8_SD2_gr8_GT5']:
    dOffsClr = {S_MIN: {'R': 0.4, 'G': 0., 'B': 0.2},
                S_MAX: {'R': 0.8, 'G': 0.8, 'B': 0.8},
                S_INV: {'R': True, 'G': True, 'B': True}}
elif sTask in ['SRR_AllGT_Idx_gr8']:
    dOffsClr = {S_MIN: {'R': 0., 'G': 0.6, 'B': 0.},
                S_MAX: {'R': 0.8, 'G': 0.8, 'B': 0.8},
                S_INV: {'R': True, 'G': True, 'B': True}}
elif sTask in lTaskAvgMax:
    dOffsClr = {S_MIN: {'R': 0.4, 'G': 0., 'B': 0.2},
                S_MAX: {'R': 0.8, 'G': 0.8, 'B': 0.8},
                S_INV: {'R': True, 'G': True, 'B': True}}

# node input ------------------------------------------------------------------
lwdNodeAttr1 = 1
lwdNodeAttr2 = 1.5
hoverInfNodeAttr1 = 'text'
hoverInfNodeAttr2 = 'text'
sConNodeAttr1 = '# of connections: '
sConNodeAttr2 = '# of connections: '
modeNodeAttr1 = 'markers+text'
modeNodeAttr2 = 'markers+text'
szTxtNodeAttr1 = 4
szTxtNodeAttr2 = 5
clrTxtNodeAttr1 = clrDk1P
clrTxtNodeAttr2 = clrR

# node scale input ------------------------------------------------------------
showScaleNodeAttr1 = True
showScaleNodeAttr2 = True
clrScaleNodeAttr1 = 'YlGnBu'
clrScaleNodeAttr2 = 'YlGnBu'
revScaleNodeAttr1 = True
revScaleNodeAttr2 = True
szMarkAttr1 = 12
szMarkAttr2 = 16

# colour bar input ------------------------------------------------------------
clrBarRelXAttr1 = -0.15
clrBarRelXAttr2 = 1.02
clrBarThicknAttr1 = 15
clrBarThicknAttr2 = 15
clrBarTitleAttr1 = 'Node Connections (Metabolites)'
clrBarTitleAttr2 = 'Node Connections (Phosphopeptides)'
if sColAttr2 in [S_BIN_C_2]:
    clrBarTitleAttr2 = 'Node Connections (BinCode2)'
clrBarXAnchAttr1 = 'left'
clrBarXAnchAttr2 = 'left'
clrBarTSideAttr1 = 'right'
clrBarTSideAttr2 = 'right'

# title input -----------------------------------------------------------------
figTitle = ('<br>Network graph of all genotype data for the 459 most ' +
            'positively congruent metabolite - phosphopeptide patterns')
if sTask == 'Idx_gr8_SD2_gr8_GT0':
    figTitle = ('<br>Network graph for the 266 most positively congruent me' +
                'tabolite - phosphopeptide patterns (wild type)')
elif sTask == 'Idx_gr8_SD2_gr8_GT1':
    figTitle = ('<br>Network graph for the 221 most positively congruent me' +
                'tabolite - phosphopeptide patterns (PGM)')
elif sTask == 'Idx_gr8_SD2_gr8_GT2':
    figTitle = ('<br>Network graph for the 105 most positively congruent me' +
                'tabolite - phosphopeptide patterns (SIRK1)')
elif sTask == 'Idx_gr8_SD2_gr8_GT3':
    figTitle = ('<br>Network graph for the 271 most positively congruent me' +
                'tabolite - phosphopeptide patterns (SIRK1-PGM)')
elif sTask == 'Idx_gr8_SD2_gr8_GT4':
    figTitle = ('<br>Network graph for the 73 most positively congruent me' +
                'tabolite - phosphopeptide patterns (SIRK1-SWEET)')
elif sTask == 'Idx_gr8_SD2_gr8_GT5':
    figTitle = ('<br>Network graph for the 158 most positively congruent me' +
                'tabolite - phosphopeptide patterns (SWEET)')
elif sTask == 'SRR_AllGT_Idx_gr8':
    figTitle = ('<br>Network graph for the 459 most positively congruent me' +
                'tabolite - phosphopeptide patterns (mean over all genotypes)')
elif sTask in lTaskAvgGT0:
    figTitle = ('<br>Network graph: Metabolites and phosphopeptides' +
                ' (avg. over BinCode2 groups) with similar behaviour (WT)')
elif sTask in lTaskAvgGT1:
    figTitle = ('<br>Network graph: Metabolites and phosphopeptides' +
                ' (avg. over BinCode2 groups) with similar behaviour (PGM)')
elif sTask in lTaskAvgGT5:
    figTitle = ('<br>Network graph: Metabolites and phosphopeptides' +
                ' (avg. over BinCode2 groups) with similar behaviour (SWEET)')
elif sTask in lTaskMaxGT0:
    figTitle = ('<br>Network graph: Metabolites and phosphopeptides' +
                ' (max. over BinCode2 groups) with similar behaviour (WT)')
elif sTask in lTaskMaxGT1:
    figTitle = ('<br>Network graph: Metabolites and phosphopeptides' +
                ' (max. over BinCode2 groups) with similar behaviour (PGM)')
elif sTask in lTaskMaxGT5:
    figTitle = ('<br>Network graph: Metabolites and phosphopeptides' +
                ' (max. over BinCode2 groups) with similar behaviour (SWEET)')
fontSzTitle = 11
showLegTitle = False
hoverModeTitle = 'closest'
dMarginTitle = dict(b = 20, l = 5, r = 5, t = 40)

# annotations input -----------------------------------------------------------
textAnnot = 'Dark purple'
if sTask in ['Idx_gr8_SD2_gr8_GT0', 'Idx_gr8_SD2_gr8_GT1',
             'Idx_gr8_SD2_gr8_GT2', 'Idx_gr8_SD2_gr8_GT3',
             'Idx_gr8_SD2_gr8_GT4', 'Idx_gr8_SD2_gr8_GT5']:
    textAnnot = 'Dark purple'
elif sTask in ['SRR_AllGT_Idx_gr8']:
    textAnnot = 'Green'
elif sTask in lTaskAvgMax:
    textAnnot = 'Dark purple'
textAnnot += ' edge colour corresonds to higher confidence of connection'
showArrAnnot = False
xRefAnnot = 'paper'
yRefAnnot = 'paper'
xPosAnnot = 0.005
yPosAnnot = -0.002
xAxisAnnot = dict(showgrid = False, zeroline = False, showticklabels = False)
yAxisAnnot = dict(showgrid = False, zeroline = False, showticklabels = False)

# === INPUT PROCESSING ========================================================
nmDatF += ('.' + S_EXT_CSV)
nmPltF += ('.' + S_EXT_PDF)
# lAttr = [ATTR_1, ATTR_2]
lSCol = [sColAttr1, sColAttr2]
lLwdNode = [lwdNodeAttr1, lwdNodeAttr2]
lHoverInfNode = [hoverInfNodeAttr1, hoverInfNodeAttr2]
lSConNode = [sConNodeAttr1, sConNodeAttr2]
lModeNode = [modeNodeAttr1, modeNodeAttr2]
lSzTxtNode = [szTxtNodeAttr1, szTxtNodeAttr2]
lClrTxtNode = [clrTxtNodeAttr1, clrTxtNodeAttr2]
lShowScaleNode = [showScaleNodeAttr1, showScaleNodeAttr2]
lClrScaleNode = [clrScaleNodeAttr1, clrScaleNodeAttr2]
lRevScaleNode = [revScaleNodeAttr1, revScaleNodeAttr2]
lSzMark = [szMarkAttr1, szMarkAttr2]
lClrBarRelX = [clrBarRelXAttr1, clrBarRelXAttr2]
lClrBarThickn = [clrBarThicknAttr1, clrBarThicknAttr2]
lClrBarTitle = [clrBarTitleAttr1, clrBarTitleAttr2]
lClrBarXAnch = [clrBarXAnchAttr1, clrBarXAnchAttr2]
lClrBarTSide = [clrBarTSideAttr1, clrBarTSideAttr2]

# === FUNCTIONS ===============================================================
def createDir(pF):
    if not os.path.isdir(pF):
        os.mkdir(pF)

def joinToPath(pF = '', nmF = 'Dummy.txt'):
    if len(pF) > 0:
        createDir(pF)
        return os.path.join(pF, nmF)
    else:
        return nmF

def readCSV(pF, sepD = ',', iCol = None, dDtype = None, lHStr = [],
            iSp = None, splC = True):
    if dDtype is None and lHStr is not None:
        dDtype = {s: str for s in lHStr}
    pdDfr = pd.read_csv(pF, sep = sepD, index_col = iCol, dtype = dDtype)
    if iSp is not None:
        if splC:
            addIDfr = pdDfr.iloc[:, :iSp]
            return pdDfr.iloc[:, iSp:], addIDfr
        else:
            addIDfrT = pdDfr.iloc[:iSp, :].T
            return pdDfr.iloc[iSp:, :], addIDfrT.T
    else:
        return pdDfr

def flattenIt(cIterable, retArr = False):
    itFlat = list(itertools.chain.from_iterable(cIterable))
    if retArr:
        itFlat = np.array(itFlat)
    return itFlat

def concPdDfrS(lPdDfr, concAx = 0, verInt = True, srt = False, ignIdx = False,
               dropAx = None):
    d = pd.concat(lPdDfr, axis = concAx, verify_integrity = verInt, sort = srt,
                  ignore_index = ignIdx)
    if dropAx in [0, 1, 'index', 'columns']:
        d.dropna(axis = dropAx, inplace = True)
    return d

def splitDfr(pdDfr, tHd, j = 0):
    lSubDfr, setV = [], set(pdDfr[tHd[j]])
    for cV in setV:
        lSubDfr.append(pdDfr[pdDfr[tHd[j]] == cV])
    if j == len(tHd) - 1:
        return lSubDfr
    else:
        j += 1
        return [splitDfr(cSubDfr, tHd, j) for cSubDfr in lSubDfr]

def calcFromCVal(lDfr, numK, tHd, cOp = S_AVG):
    lDfrR = []
    # flatten the possibly nested list of DataFrames first
    for k in range(len(tHd) - 1):
        lDfr = flattenIt(lDfr)
    if cOp in [S_MAX, S_AVG]:
        # get the maximum or average of all values with same tHd-def. columns
        for cDfr in lDfr:
            if cOp == S_MAX:
                cDfr = cDfr[cDfr[numK] == max(cDfr[numK])]
                if cDfr.shape[0] > 1:
                    cDfr = cDfr.iloc[0, :].to_frame().T
            elif cOp == S_AVG:
                cDfr.loc[:, numK] = round(np.mean(cDfr[numK]), N_DIG_RND_04)
                cDfr = cDfr.iloc[0, :].to_frame().T
            lDfrR.append(cDfr)
        return concPdDfrS(lDfrR, ignIdx = True)

def sortAndFilter(pdDfr, cOp = S_AVG, dSrt = None, dFlt = None, lSCN = None):
    dfrM = pdDfr.copy()
    if dFlt is not None:
        if S_SEL in dFlt:       # filter based on list of attributes
            for cK, cL in dFlt[S_SEL].items():
                dfrM = dfrM[dfrM[cK].isin(cL)].reset_index(drop = True)
        if S_THR in dFlt:       # filter based on numeric threshold
            for (sCN, (sCmp, cThr)) in dFlt[S_THR].items():
                if sCmp == '>=':
                    dfrM = dfrM[dfrM[sCN] >= cThr].reset_index(drop = True)
                elif sCmp == '<=':
                    dfrM = dfrM[dfrM[sCN] <= cThr].reset_index(drop = True)
                else:
                    print('ERROR: Operation "', sCmp, '" not implemented.')
        if cOp in dFlt:
            # max. / avg. over all entries with the same col. header key tuple
            cKNum, tHdr = dFlt[cOp]
            dfrM = calcFromCVal(splitDfr(dfrM, tHdr), cKNum, tHdr, cOp = cOp)
    if dSrt is not None:
        dfrM = dfrM.sort_values(list(dSrt), ascending = list(dSrt.values()),
                                ignore_index = True)
    return dfrM.dropna(subset = lSCN)

def getInfoFromG(G, cDfr, lSC):
    llA = []
    for k in range(len(lSC)):
        llA.append(list(set(cDfr.loc[:, lSC[k]])))
    return llA, nx.kamada_kawai_layout(G)

def getSClr(cClrC, cMult = MAX_CLR_VAL):
    return str(round(cClrC*cMult))

def getlClr(G, dOClr, iDat = 2, sEdgeAt = 'weight', sRGB = S_RGB):
    lEdgeAt = [cDat[iDat][sEdgeAt] for cDat in G.edges.data()]
    minV, maxV = min(lEdgeAt), max(lEdgeAt)
    dClr = {'R': [], 'G': [], 'B': []}
    for x in lEdgeAt:
        for sC in dClr:
            vE, vD = minV, maxV
            if dOClr[S_INV][sC]:
                vE, vD = maxV, minV
            cC = (dOClr[S_MIN][sC]*(vD - x)/(vD - vE) +
                  dOClr[S_MAX][sC]*(x - vE)/(vD - vE))
            dClr[sC].append(cC)
    return [sRGB + '(' + getSClr(dClr['R'][k]) + ', ' + getSClr(dClr['G'][k]) +
            ', ' + getSClr(dClr['B'][k]) + ')' for k in range(len(lEdgeAt))]

def getLEdgeTrace(G, dPos, sNumAt, eWdth, hovInfE, modE, dOClr = {}):
    lEdgeTr = []
    lClr = getlClr(G, dOClr, iDat = INDEX_DATA_EDGE_ATTR, sEdgeAt = sNumAt)
    assert len(lClr) == len(G.edges())
    for k, cEdge in enumerate(G.edges()):
        ((x0, y0), (x1, y1)) = (dPos[cEdge[0]], dPos[cEdge[1]])
        dLine = dict(width = eWdth, color = lClr[k])
        cEdgeTr = go.Scatter(x = [x0, x1], y = [y0, y1], line = dLine,
                             hoverinfo = hovInfE, mode = modE)
        lEdgeTr.append(cEdgeTr)
    return lEdgeTr

def getLNodeInfo(G, dPos, lSAttr, sNdX = ''):
    lAdjInfo = list(G.adjacency())
    assert len(G.nodes()) == len(lAdjInfo)
    lCNdX, lCNdY, lNAdjNd, lTxtAdjNd, lS = [], [], [], [], []
    for k, sNd in enumerate(G.nodes()):
        if sNd in lSAttr:
            lS.append(sNd)
            lCNdX.append(dPos[sNd][0])
            lCNdY.append(dPos[sNd][1])
            nAdj = len(lAdjInfo[k][1])
            lNAdjNd.append(nAdj)
            lTxtAdjNd.append(sNdX + str(nAdj))
    return lS, lCNdX, lCNdY, lNAdjNd, lTxtAdjNd

def getNodeTrace(G, dPos, lSAttr, szMk, lwdNd, hoverInfNd, sConNd, modeNd,
                 szTxtNd, clrTxtNd, showScNd, clrScNd, revScNd, clrBarRelX,
                 clrBarThickn, clrBarTitle, clrBarXAnch, clrBarTSide):
    lS, lCX, lCY, lNCon, _ = getLNodeInfo(G, dPos, lSAttr, sNdX = sConNd)
    dClrBar = dict(x = clrBarRelX, thickness = clrBarThickn,
                   title = clrBarTitle, xanchor = clrBarXAnch,
                   titleside = clrBarTSide)
    dMark = dict(showscale = showScNd, colorscale = clrScNd,
                 reversescale = revScNd, color = lNCon, size = szMk,
                 colorbar = dClrBar, line_width = lwdNd)
    nodeTr = go.Scatter(x = lCX, y = lCY, mode = modeNd,
                        hoverinfo = hoverInfNd, marker = dMark)
    nodeTr.text = lS
    nodeTr.textfont = dict(size = szTxtNd, color = clrTxtNd)
    return nodeTr

# === MAIN ====================================================================
print('Task:', sSelNmF)
dfrSF = sortAndFilter(pd.read_csv(joinToPath(pRelDatF, nmDatF), sep = cSep),
                      cOp = selOp,
                      dSrt = dSort_MB2ND,
                      dFlt = dFilt_MBC2IPDC,
                      lSCN = list(dThrN_IPDC))
dfrSF.to_csv(joinToPath(pRelDatF_92, sSelNmF + '_' + nmDatF), sep = cSep)
sCA1, sCA2, sNumA, lNodeTrace = sColAttr1, sColAttr2, sNumAttr1, []
G = nx.convert_matrix.from_pandas_edgelist(dfrSF, source = sCA1, target = sCA2,
                                           edge_attr = sNumA)
llSAttr, dPos = getInfoFromG(G, dfrSF, lSCol)
lEdgeTrace = getLEdgeTrace(G, dPos, sNumA, eWdth = edgeWidth,
                           hovInfE = hoverInfEdge, modE = modeEdge,
                           dOClr = dOffsClr)

for k in range(len(llSAttr)):
    lNodeTrace.append(getNodeTrace(G, dPos,
                                   lSAttr = llSAttr[k],
                                   szMk = lSzMark[k],
                                   lwdNd = lLwdNode[k],
                                   hoverInfNd = lHoverInfNode[k],
                                   sConNd = lSConNode[k],
                                   modeNd = lModeNode[k],
                                   szTxtNd = lSzTxtNode[k],
                                   clrTxtNd = lClrTxtNode[k],
                                   showScNd = lShowScaleNode[k],
                                   clrScNd = lClrScaleNode[k],
                                   revScNd = lRevScaleNode[k],
                                   clrBarRelX = lClrBarRelX[k],
                                   clrBarThickn = lClrBarThickn[k],
                                   clrBarTitle = lClrBarTitle[k],
                                   clrBarXAnch = lClrBarXAnch[k],
                                   clrBarTSide = lClrBarTSide[k]))

lAnnot = [dict(text = textAnnot, showarrow = showArrAnnot, xref = xRefAnnot,
               yref = yRefAnnot, x = xPosAnnot, y = yPosAnnot)]
fig = go.Figure(data = lEdgeTrace + lNodeTrace,
                layout = go.Layout(title = figTitle,
                                   titlefont_size = fontSzTitle,
                                   showlegend = showLegTitle,
                                   hovermode = hoverModeTitle,
                                   margin = dMarginTitle,
                                   annotations = lAnnot,
                                   xaxis = xAxisAnnot, yaxis = yAxisAnnot))
fig.write_image(joinToPath(pRelPltF, nmPltF))
print('-'*37, 'DONE', '-'*37)

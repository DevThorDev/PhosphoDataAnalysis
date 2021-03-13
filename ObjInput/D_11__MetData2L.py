# -*- coding: utf-8 -*-
###############################################################################
# --- D_11__MetData2L.py ------------------------------------------------------
###############################################################################
import os

import Core.C_00__GenConstants as GC
from Control.A_00__GenInput import dIDGT, dIDFt

# --- general -----------------------------------------------------------------
sOType = 'Metabolite data set (log2-transformed)'
sNmSpec = GC.S_MET_D_L + GC.S_2L
dropNaNAx_F = None          # axis to drop NaN - all values (None / 0 / 1)
dropNaNAx_Mn = 1            # axis to drop NaN - mean values (None / 0 / 1)

# --- data specific -----------------------------------------------------------
cSep = ';'
lTransf = ['log2']*len(dIDFt)
iColNm = 0                  # index of (unique) name column
iColVsS = 2                 # index of first column of data array
lIColConv2Str = []          # list of indices of columns converted to string
lHdrAIC = []                # column header list of classification info
lHdrAIR = lHdrAIC + [GC.S_MN_CONC] # column header list of retained add. info 
tIMnMx = (-1, -1)           # index = -1 --> no restriction of start/end
iFtrBase = 0                # index of feature acting as baseline for trans.
transD = 'No'              # 'No' / 'Abs' / 'Rel' / '2LQ' (trans. shift type)
stdOp = 'CMS'                # 'N'(o) / 'T'(otal) / 'R'(ow) / 'C'(ol)
                            # [+ 'M'(ean)] [+ 'S'(td.dev)] optional adjustm.
devTp = 'SD'                # 'SD' / 'SDP' (units deviation is measured in)
                            # 'SDP': used pooled deviation

# --- filtering ---------------------------------------------------------------
lSFiltOut = ['Allose_(1MEOX)_(5TMS)', 'Talose__D-_(1MEOX)_(5TMS)',
             'Maltitol_(9TMS)', 'Maltose__D-_(1MEOX)_(8TMS)']

# --- clustering --------------------------------------------------------------
minNClusters, maxNClusters = (2, 10)                # no filter
dNCl4Corr = {(GC.NM_XFEAT_C, None): [10],             # checked
             (GC.NM_XFEAT_G, None): [4],             # yet to check
             (GC.NM_XFEAT_C, GC.NM_TYPEX_S): [7],    # checked
             (GC.NM_XFEAT_C, GC.NM_TYPEX_F): [3],    # yet to check
             (GC.NM_XFEAT_G, GC.NM_TYPEX_S): [5],    # yet to check
             (GC.NM_XFEAT_G, GC.NM_TYPEX_F): [8]}    # yet to check
# minNClusters, maxNClusters = (1, 1)                # filter = 'Sucrose'
# dNCl4Corr = {(GC.NM_XFEAT_C, None): [],
#              (GC.NM_XFEAT_G, None): [],
#              (GC.NM_XFEAT_C, GC.NM_TYPEX_S): [],
#              (GC.NM_XFEAT_C, GC.NM_TYPEX_F): [],
#              (GC.NM_XFEAT_G, GC.NM_TYPEX_S): [],
#              (GC.NM_XFEAT_G, GC.NM_TYPEX_F): []}
distThrAgglo = None

# --- names and paths of files and dirs ---------------------------------------
dISPr = {'No': {0: '', 4: '_', 6: '_'},
         'Tr': {0: '', 1: '_', 4: '_', 6: '_'},
         'Dv': {0: '', 1: '_', 4: '_', 6: '_'}}
nmDFPre = 'Metabolite_log2'
dIDDF = {idGT: nmDFPre + '_' + idGT + '_' + nmGT + '_ThS' for
         (idGT, nmGT) in dIDGT.items()}

pRelDF = os.path.join('..', '..', '12_SysBio02_DataAnalysis',
                      '02_ProcessedRawData', '01_CSV')

# --- graphics parameters -----------------------------------------------------
plotHist = True         # plot a histogram of the mean data?
lISep = [2]       # list of separation indices between categories
llIOffsClr = [[0, 2], [1, 3]]   # list of colour offset indices lists
nBins = 20              # number of bins
histAlpha = 0.4         # alpha for histogram

# --- derived values and assertions -------------------------------------------
dPRelDF = {sGT: os.path.join(pRelDF, dIDDF[sGT] + '.' + GC.NM_EXT_CSV) for
           sGT in dIDDF}
lNClusters = list(range(minNClusters, maxNClusters + 1))

assert len(stdOp) >= 1

# --- create input dictionary -------------------------------------------------
dIO = {# --- general
       'sOType': sOType,
       'sNmSpec': sNmSpec,
       'dropAx_F': dropNaNAx_F,
       'dropAx_Mn': dropNaNAx_Mn,
       # --- data specific
       'cSep': cSep,
       'lTransf': lTransf,
       'iColNm': iColNm,
       'iColVsS': iColVsS,
       'lIColConv2Str': lIColConv2Str,
       'lHdrAIC': lHdrAIC,
       'lHdrAIR': lHdrAIR,
       'tIMnMx': tIMnMx,
       'iFtrBase': iFtrBase,
       'transD': transD,
       'stdOp': stdOp,
       'devTp': devTp,
       # --- filtering
       'lSFiltOut': lSFiltOut,
       # --- clustering
       'minNClusters': minNClusters,
       'maxNClusters': maxNClusters,
       'lNClusters': lNClusters,
       'dNCl4Corr': dNCl4Corr,
       'distThrAgglo': distThrAgglo,
       # --- names and paths of files and dirs
       'dISPr': dISPr,
       'nmDFPre': nmDFPre,
       'dIDDF': dIDDF,
       'dPRelDF': dPRelDF,
       # --- graphics parameters
       'plotHist': plotHist,
       'lISep': sorted(lISep),
       'llIOffsClr': llIOffsClr,
       'nBins': nBins,
       'histAlpha': histAlpha}

###############################################################################
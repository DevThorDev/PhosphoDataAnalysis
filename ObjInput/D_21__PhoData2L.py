# -*- coding: utf-8 -*-
###############################################################################
# --- D_21__PhoData2L.py ------------------------------------------------------
###############################################################################
import os

import Core.C_00__GenConstants as GC
from Control.A_00__GenInput import dIDGT, dIDFt

# --- general -----------------------------------------------------------------
sOType = 'Phosphopeptide data set (log2-transformed)'
sNmSpec = GC.S_PHO_D_L + GC.S_2L
dropNaNAx_F = None          # axis to drop NaN - all values (None / 0 / 1)
dropNaNAx_Mn = 1            # axis to drop NaN - mean values (None / 0 / 1)

# --- data specific -----------------------------------------------------------
cSep = ';'
lTransf = ['log2']*len(dIDFt)
iColNm = 0                  # index of (unique) name column
iColVsS = 14                # index of first column of data array
lIColConv2Str = [2, 3, 4, 5, 12] # list of indices of columns conv. to str
lHdrAIC = [GC.S_PROT, GC.S_BIN_C_G, GC.S_BIN_C_3, GC.S_BIN_C_2, GC.S_BIN_C_1]
lHdrAIR = lHdrAIC + [GC.S_MN_CONC] # column header list of retained add. info
tIMnMx = (-1, -1)           # index = -1 --> no restriction of start/end
iFtrBase = 0                # index of feature acting as baseline for trans.
transD = 'No'              # 'No' / 'Abs' / 'Rel' / '2LQ' (trans. shift type)
stdOp = 'CMS'                # 'N'(o) / 'T'(otal) / 'R'(ow) / 'C'(ol)
                            # [+ 'M'(ean)] [+ 'S'(td.dev)] optional adjustm.
devTp = 'SD'                # 'SD' / 'SDP' (units deviation is measured in)
                            # 'SDP': used pooled deviation

# --- filtering ---------------------------------------------------------------
lSFiltOut = []

# --- clustering --------------------------------------------------------------
minNClusters, maxNClusters = (2, 12)
dNCl4Corr = {(GC.NM_XFEAT_C, None): [11],            # checked
             (GC.NM_XFEAT_G, None): [8],             # yet to check
             (GC.NM_XFEAT_C, GC.NM_TYPEX_S): [11],   # checked
             (GC.NM_XFEAT_C, GC.NM_TYPEX_F): [9],    # yet to check
             (GC.NM_XFEAT_G, GC.NM_TYPEX_S): [10],   # yet to check
             (GC.NM_XFEAT_G, GC.NM_TYPEX_F): [12]}   # yet to check
distThrAgglo = None

# --- names and paths of files and dirs ---------------------------------------
dISPr = {'No': {0: '', 4: '_', 6: '_'},
         'Tr': {0: '', 1: '_', 4: '_', 6: '_'},
         'Dv': {0: '', 1: '_', 4: '_', 6: '_'}}
nmDFPre = 'Phosphopeptide_log2_OutlMedian5p0'
dIDDF = {idGT: nmDFPre + '_' + idGT + '_' + nmGT + '_ThS' for
         (idGT, nmGT) in dIDGT.items()}

pRelDF = os.path.join('..', '..', '12_SysBio02_DataAnalysis',
                      '05_RevisedRawData', '01_CSV')

# --- graphics parameters -----------------------------------------------------
plotHist = True         # plot a histogram of the mean data?
lISep = [2]       # list of separation indices between categories
llIOffsClr = [[5, 8], [6, 9]]   # list of colour offset indices lists
nBins = 100              # number of bins
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

# -*- coding: utf-8 -*-
###############################################################################
# --- D_12__MetData_Red1.py ---------------------------------------------------
###############################################################################
import os

import Core.C_00__GenConstants as GC

# --- general -----------------------------------------------------------------
sOType = 'Metabolite data set, reduced by removing all rows with >= one NaN'
sNmSpec = 'MetData_Red1'
sNmMut = 'WT'

# --- data specific -----------------------------------------------------------
cSep = ';'
lCat = ['DR-WT', 'DS-WT', 'NR-WT', 'NS-WT']
lDesc = ['day_root', 'day_shoot ', 'night_root', 'night_shoot']
#lNRep = [6, 6, 6, 6]
LTransf = ['log2', 'log2', 'log2', 'log2']
iColNm = 0                  # index of (unique) name column
iColVsS = 1                 # index of first column of data array
lIColConv2Str = []          # list of indices of columns converted to string
tIMnMx = (-1, -1)           # index = -1 --> no restriction of start/end
dropDupl = True             # drop duplicate lines (True/False)
dropMeth = 'nameIdent'      # 'fullIdent' / 'nameIdent'

# --- clustering --------------------------------------------------------------
minNClusters, maxNClusters = (2, 20)
lNClusters4Corr = [2, 12, 20]
distThrAgglo = None

# --- names and paths of files and dirs ---------------------------------------
nmDF_noX = 'Metabolite_log2_WT_Red1_NoNaN_ThS'
nmRF_Means = 'MetData_Means'
nmRF_SDs = 'MetData_SDs'

pRelDF = os.path.join('..', '..', '12_SysBio02_DataAnalysis',
                      '02_ProcessedRawData', '01_WT_Only', '01_CSV')
pRelResF = os.path.join('..', '..', '12_SysBio02_DataAnalysis',
                        '21_Results__11__MetData_WT')
pRelPltF = os.path.join('..', '..', '12_SysBio02_DataAnalysis',
                        '41_Figures__11__MetData_WT')

# --- graphics parameters -----------------------------------------------------
plotHist = True         # plot a histogram of the mean data?
nBins = 20              # number of bins
histAlpha = 0.4         # alpha for histogram

# --- derived values ----------------------------------------------------------
lNChrCat = [len(cCat) for cCat in lCat]
lNClusters = list(range(minNClusters, maxNClusters + 1))

# --- create input dictionary -------------------------------------------------
dIO = {# --- general
       'sOType': sOType,
       'sNmSpec': sNmSpec,
       'sNmMut': sNmMut,
       # --- data specific
       'cSep': cSep,
       'lCat': lCat,
       'lDesc': lDesc,
#       'lNRep': lNRep,
       'LTransf': LTransf,
       'iColNm': iColNm,
       'iColVsS': iColVsS,
       'lIColConv2Str': lIColConv2Str,
       'tIMnMx': tIMnMx,
       'dropDupl': dropDupl,
       'dropMeth': dropMeth,
       # --- names and paths of files and dirs
       'nmRF_Means': nmRF_Means,
       'nmRF_SDs': nmRF_SDs,
       'pRelDF': os.path.join(pRelDF, nmDF_noX + '.' + GC.N_EXT_CSV),
       'pRelResF': pRelResF,
       'pRelPltF': pRelPltF,
       # --- graphics parameters
       'plotHist': plotHist,
       'nBins': nBins,
       'histAlpha': histAlpha,
       # --- derived values
       'lNChrCat': lNChrCat}

###############################################################################

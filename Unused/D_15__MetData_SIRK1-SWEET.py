# -*- coding: utf-8 -*-
###############################################################################
# --- D_15__MetData_SIRK1-SWEET.py --------------------------------------------
###############################################################################
import os

import Core.C_00__GenConstants as GC

# --- general -----------------------------------------------------------------
sOType = 'Metabolite data set'
sNmSpec = 'MetData'
sNmMut = 'SIRK1-SWEET'

# --- data specific -----------------------------------------------------------
cSep = ';'
lCat = ['DR-SIRK1-SWEET', 'DS-SIRK1-SWEET', 'NR-SIRK1-SWEET', 'NS-SIRK1-SWEET']
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
minNClusters, maxNClusters = (2, 10)
lNClusters4Corr = [2, 5, 8]
distThrAgglo = None

# --- names and paths of files and dirs ---------------------------------------
nmDF_noX = 'Metabolite_log2_SIRK1-SWEET_ThS'
nmRF_Means = 'MetData_Means'
nmRF_SDs = 'MetData_SDs'

pRelDF = os.path.join('..', '..', '12_SysBio02_DataAnalysis',
                      '02_ProcessedRawData', '02_Mutants_All', '01_CSV')
pRelResF = os.path.join('..', '..', '12_SysBio02_DataAnalysis',
                        '21_Results__15__MetData_SIRK1-SWEET')
pRelPltF = os.path.join('..', '..', '12_SysBio02_DataAnalysis',
                        '41_Figures__15__MetData_SIRK1-SWEET')

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
       # --- clustering
       'minNClusters': minNClusters,
       'maxNClusters': maxNClusters,
       'lNClusters': lNClusters,
       'lNClusters4Corr': lNClusters4Corr,
       'distThrAgglo': distThrAgglo,
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

# -*- coding: utf-8 -*-
###############################################################################
# --- D_97__TestDataPho.py ----------------------------------------------------
###############################################################################
import os

import Core.C_00__GenConstants as GC

# --- general -----------------------------------------------------------------
strOType = 'Test Pho data set'
strNmSpec = 'TestDataPho'
strNmMut = 'WT'

# --- data specific -----------------------------------------------------------
cSep = ';'
lCat = ['DR-WT', 'NR-WT']
lDesc = ['day_root', 'night_root']
lNRep = [6, 6]
LTransf = ['log2', 'log2']
iColNm = 0
iColVsS = 7
tIMnMx = (-1, -1)

# --- names and paths of files and dirs ---------------------------------------
nDF_noX = 'TestData_Pho_log2_Intensitites_WT_ThS'

pRelDF = os.path.join('..', '..', '12_SysBio02_DataAnalysis',
                      '02_ProcessedRawData', '01_WT_Only', '01_CSV')

# --- derived values ----------------------------------------------------------
lNChrCat = [len(cCat) for cCat in lCat]

# --- create input dictionary -------------------------------------------------
dIO = {# --- general
       'strOType': strOType,
       'strNmSpec': strNmSpec,
       'strNmMut': strNmMut,
       # --- data specific
       'cSep': cSep,
       'lCat': lCat,
       'lDesc': lDesc,
       'lNRep': lNRep,
       'LTransf': LTransf,
       'iColNm': iColNm,
       'iColVsS': iColVsS,
       'tIMnMx': tIMnMx,
       # --- names and paths of files and dirs
       'pRelDF': os.path.join(pRelDF, nDF_noX + '.' + GC.N_EXT_CSV),
       # --- derived values
       'lNChrCat': lNChrCat}

###############################################################################

# -*- coding: utf-8 -*-
###############################################################################
# --- D_91__PhoData.py --------------------------------------------------------
###############################################################################
import os

import Core.C_00__GenConstants as GC

# --- general -----------------------------------------------------------------
strOType = 'Phospho data set'
strNmSpec = 'PhoData'
strNmMut = 'WT'

# --- data specific -----------------------------------------------------------
cSep = ';'
lCat = ['DR-WT', 'DS-WT', 'NR-WT', 'NS-WT']
lDesc = ['day_root', 'day_shoot ', 'night_root', 'night_shoot']
lNRep = [6, 6, 6, 6]
LTransf = ['log2', 'log2', 'log2', 'log2']
iColNm = 0
iColVsS = 7
lIColConv2Str = [6]
tIMnMx = (-1, -1)
dropDupl = True
dropMeth = 'nameIdent'      # 'fullIdent' / 'nameIdent'

# --- names and paths of files and dirs ---------------------------------------
nDF_noX = 'PhosphoData_log2_Intensitites_WT_ThS'

pRelDF = os.path.join('..', '..', '12_SysBio02_DataAnalysis',
                      '02_ProcessedRawData', '01_WT_Only', '01_CSV')
pRelResF = os.path.join('.', 'Results')

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
       'lIColConv2Str': lIColConv2Str,
       'tIMnMx': tIMnMx,
       'dropDupl': dropDupl,
       'dropMeth': dropMeth,
       # --- names and paths of files and dirs
       'pRelDF': os.path.join(pRelDF, nDF_noX + '.' + GC.N_EXT_CSV),
       'pRelResF': pRelResF,
       # --- derived values
       'lNChrCat': lNChrCat}

###############################################################################

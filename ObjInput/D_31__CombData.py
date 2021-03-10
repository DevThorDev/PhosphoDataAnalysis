# -*- coding: utf-8 -*-
###############################################################################
# --- D_31__CombData.py -------------------------------------------------------
###############################################################################
import Core.C_00__GenConstants as GC

# --- general -----------------------------------------------------------------
sOType = 'Combined data set'
sNmSpec = GC.S_COMB_D_L
dropNaNAx_F = None          # axis to drop NaN - all values (None / 0 / 1)
dropNaNAx_Mn = 1            # axis to drop NaN - mean values (None / 0 / 1)

# --- data specific -----------------------------------------------------------
cSep = ';'

# --- clustering --------------------------------------------------------------
minNClusters, maxNClusters = (2, 10)                # 5 & 10 clusters only
dNCl4Corr = {(GC.NM_XFEAT_C, None): [5],                 # yet to check
             (GC.NM_XFEAT_G, None): [6],                 # yet to check
             (GC.NM_XFEAT_C, GC.NM_TYPEX_S): [5, 10],    # yet to check
             (GC.NM_XFEAT_C, GC.NM_TYPEX_F): [7],        # yet to check
             (GC.NM_XFEAT_G, GC.NM_TYPEX_S): [8],        # yet to check
             (GC.NM_XFEAT_G, GC.NM_TYPEX_F): [9]}        # yet to check
# minNClusters, maxNClusters = (2, 14)
# dNCl4Corr = {(GC.NM_XFEAT_C, None): [2, 8],
#              (GC.NM_XFEAT_G, None): [2, 8],
#              (GC.NM_XFEAT_C, GC.NM_TYPEX_S): [3, 4, 13],
#              (GC.NM_XFEAT_C, GC.NM_TYPEX_F): [3, 4, 13],
#              (GC.NM_XFEAT_G, GC.NM_TYPEX_S): [2, 3, 4, 5],
#              (GC.NM_XFEAT_G, GC.NM_TYPEX_F): [2, 3, 4, 5]}
distThrAgglo = None

# --- names and paths of files and dirs ---------------------------------------
dISPr = {'No': {0: ''}, 'Tr': {0: ''}, 'Dv': {0: ''}}
nmCombDataF = 'CombData'

# --- derived values ----------------------------------------------------------
lNClusters = list(range(minNClusters, maxNClusters + 1))

# --- create input dictionary -------------------------------------------------
dIO = {# --- general
       'sOType': sOType,
       'sNmSpec': sNmSpec,
       'dropAx_F': dropNaNAx_F,
       'dropAx_Mn': dropNaNAx_Mn,
       # --- data specific
       'cSep': cSep,
       # --- clustering
       'minNClusters': minNClusters,
       'maxNClusters': maxNClusters,
       'lNClusters': lNClusters,
       'dNCl4Corr': dNCl4Corr,
       'distThrAgglo': distThrAgglo,
       # --- names and paths of files and dirs
       'dISPr': dISPr,
       'nmCombDataF': nmCombDataF}

###############################################################################

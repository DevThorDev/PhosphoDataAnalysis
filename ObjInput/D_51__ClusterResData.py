# -*- coding: utf-8 -*-
###############################################################################
# --- D_51__ClusterResData.py -------------------------------------------------
###############################################################################
import Core.C_00__GenConstants as GC

# --- general -----------------------------------------------------------------
sOType = 'Cluster result data'
sNmSpec = GC.S_CLR_D_L

# --- data specific -----------------------------------------------------------
cSep = ';'

# --- names and paths of files and dirs ---------------------------------------
dISPr = {'No': {0: ''}, 'Tr': {0: ''}, 'Dv': {0: ''}}

# --- create input dictionary -------------------------------------------------
dIO = {# --- general
       'sOType': sOType,
       'sNmSpec': sNmSpec,
       # --- data specific
       'cSep': cSep,
       # --- names and paths of files and dirs
       'dISPr': dISPr}

###############################################################################

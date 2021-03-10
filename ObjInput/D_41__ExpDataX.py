# -*- coding: utf-8 -*-
###############################################################################
# --- D_41__ExpDataX.py -------------------------------------------------------
###############################################################################
import Core.C_00__GenConstants as GC

# --- general -----------------------------------------------------------------
sOType = 'Extended experimental data set'
sNmSpec = GC.S_EXP_DX_L

# --- data specific -----------------------------------------------------------
cSep = ';'
lKeysSpec = ['pRelDF']      # keys in objects with special treatment

# --- comparison between different genotypes / features -----------------------
sComp = 'Comp'              # name prefix for comparison file
iObjRef = 0                 # index of reference object (for comparison)
compOp = '2_log'            # None / '2_log' / 'e_log' / '10_log'

# --- names and paths of files and dirs ---------------------------------------
dISPr = {'No': {0: '', 3: '_', 4: '_', 6: '_'},
         'Tr': {0: '', 1: '_', 3: '_', 4: '_', 6: '_'},
         'Dv': {0: '', 1: '_', 3: '_', 4: '_', 6: '_'}}

# --- create input dictionary -------------------------------------------------
dIO = {# --- general
       'sOType': sOType,
       'sNmSpec': sNmSpec,
       # --- data specific
       'cSep': cSep,
       'lKeysSpec': lKeysSpec,
       # --- comparison between different genotypes / features
       'sComp': sComp,
       'iObjRef': iObjRef,
       'compOp': compOp,
       # --- names and paths of files and dirs
       'dISPr': dISPr}

###############################################################################

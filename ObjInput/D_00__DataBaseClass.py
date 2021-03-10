# -*- coding: utf-8 -*-
###############################################################################
# --- D_00__DataBaseClass.py --------------------------------------------------
###############################################################################
import os

import Core.C_00__GenConstants as GC

# --- general -----------------------------------------------------------------
sOType = 'Data set base class'
sNmSpec = GC.S_BASE_L

# --- data specific -----------------------------------------------------------
cSep = ','
dropDupl = True             # drop duplicate lines (True/False)
dropType = 'nameIdent'      # 'fullIdent' / 'nameIdent'
transD = None               # 'No' / 'Abs' / 'Rel' / '2LQ' (trans. shift type)
stdOp = None                # 'N'(o) / 'T'(otal) / 'R'(ow) / 'C'(ol)
                            # [+ 'M'(ean)] [+ 'S'(td.dev)] optional adjustm.
devTp = None                # 'SD' (units deviation is measured in)

# --- names for plots ---------------------------------------------------------
nmHist = 'Hist'         # name of histogram plot

# --- names and paths of files and dirs ---------------------------------------
dISPr = {'No': {0: '', 4: '_', 6: '_'},
         'Tr': {0: '', 1: '_', 4: '_', 6: '_'},
         'Dv': {0: '', 1: '_', 4: '_', 6: '_'}}
#pRelResF = os.path.join('..', '..', '12_SysBio02_DataAnalysis')
#pRelPltF = os.path.join('..', '..', '12_SysBio02_DataAnalysis')
pRelResF = os.path.join('..', '..', '..', '..', '..', 'Z')
pRelPltF = os.path.join('..', '..', '..', '..', '..', 'Z')

# --- create input dictionary -------------------------------------------------
dIO = {# --- general
       'sOType': sOType,
       'sNmSpec': sNmSpec,
       # --- data specific
       'cSep': cSep,
       'dropDupl': dropDupl,
       'dropType': dropType,
       'transD': transD,
       'stdOp': stdOp,
       'devTp': devTp,
       # --- names for plots
       'nmHist': nmHist,
       # --- names and paths of files and dirs
       'dISPr': dISPr,
       'pRelResF': pRelResF,
       'pRelPltF': pRelPltF}

###############################################################################

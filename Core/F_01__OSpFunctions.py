# -*- coding: utf-8 -*-
###############################################################################
# --- F_01__OSpFunctions.py ---------------------------------------------------
###############################################################################
import Core.F_00__GenFunctions as GF

# --- Functions (general) -----------------------------------------------------

# --- Functions (O_81__BinaryOps) ---------------------------------------------
def getDMap(lODat):
    dMap, llFtHdMnOD = {}, [OD.lFtHdMn for OD in lODat]
    for cHd in llFtHdMnOD[0]:
        if cHd in llFtHdMnOD[1]:
            for cNmR1 in lODat[0].cDfr.index:
                if cNmR1.startswith(cHd) and cNmR1 in lODat[1].cDfr.index:
                    GF.addToDictL(dMap, cHd, cNmR1)
    # dMap (e.g.): {'DR_WT': ['DR_WT_1', 'DR_WT_2',...], 'DS_WT':...}
    return dMap

###############################################################################

# -*- coding: utf-8 -*-
###############################################################################
# --- F_01__OSpFunctions.py ---------------------------------------------------
###############################################################################
import Core.C_00__GenConstants as GC
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

# --- Functions (O_83__OverRep) -----------------------------------------------
def getGT(tGT, dNmGT=GC.D_NM_GT, iPos=0):
    assert type(tGT) == tuple
    if GF.elEqual(tGT):
        if tGT[0] in dNmGT:
            return tGT[0], dNmGT[tGT[0]][iPos]
        else:
            return tGT[0], None
    else:
        if tGT[0] in dNmGT and GF.elEqual([sGT in dNmGT for sGT in tGT]):
            return (GC.S_USC.join([sGT for sGT in tGT]),
                    GC.S_USC.join([dNmGT[sGT][iPos] for sGT in tGT]))
        else:
            return GC.S_USC.join([sGT for sGT in tGT]), None

###############################################################################

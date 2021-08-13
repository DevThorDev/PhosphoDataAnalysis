# -*- coding: utf-8 -*-
###############################################################################
# --- D_81__BinaryOps.py ------------------------------------------------------
###############################################################################
import Core.C_00__GenConstants as GC

# --- general and names -------------------------------------------------------
calcDictR = False           # convert the result DataFrame to a dictionary?
calcForTrD = False          # calculate the correl. for the trans. basic data?
sNmSpec = GC.S_BIN_OP_L
sCt = GC.NM_CT
sNeg = GC.S_NEG
sPos = GC.S_POS
sAv = GC.S_AV
sTop = GC.S_TOP
sCorrS = GC.S_CORR_S
sCorrL = GC.S_CORR_L
sCorrV = GC.S_CORR_V
sSpearV = GC.S_SPEAR_V
sKendV = GC.S_KEND_V
sCorrP = GC.S_CORR_P
sSpearP = GC.S_SPEAR_P
sKendP = GC.S_KEND_P
sDvSc = GC.S_DV_SC
sCI = GC.S_CI
sOccCI = GC.S_OCC_CI
# sort dictionaries: keys: columns to sort, values: ascending? (then true)
dSDvSc = {sNeg + sDvSc: True, sPos + sDvSc: False}
dSDvScX = {sNeg + sDvSc: True, sPos + sDvSc: False, sPos + sNeg + sDvSc: False}
dSCI = {sNeg + sCI: True, sPos + sCI: False}
dSCIX = {sNeg + sCI: True, sPos + sCI: False, sPos + sNeg + sCI: False}

# --- correlation -------------------------------------------------------------
dPltSCorr = {(GC.T_NM_GT0[0], GC.T_NM_GT0[0]): True,    # plot strongest corr.
             (GC.T_NM_GT1[0], GC.T_NM_GT1[0]): True,
             (GC.T_NM_GT2[0], GC.T_NM_GT2[0]): True,
             (GC.T_NM_GT3[0], GC.T_NM_GT3[0]): True,
             (GC.T_NM_GT4[0], GC.T_NM_GT4[0]): True,
             (GC.T_NM_GT5[0], GC.T_NM_GT5[0]): True}

# dSCorrBnd = {(GC.T_NM_GT0[0], GC.T_NM_GT0[0]): (-0.94, 0.93),    # (lowB, upB)
#              (GC.T_NM_GT1[0], GC.T_NM_GT1[0]): (-0.96, 0.94),
#              (GC.T_NM_GT2[0], GC.T_NM_GT2[0]): (-0.94, 0.95),
#              (GC.T_NM_GT3[0], GC.T_NM_GT3[0]): (-0.95, 0.955),
#              (GC.T_NM_GT4[0], GC.T_NM_GT4[0]): (-0.93, 0.92),
#              (GC.T_NM_GT5[0], GC.T_NM_GT5[0]): (-0.9, 0.955)}
dSCorrBnd = {(GC.T_NM_GT0[0], GC.T_NM_GT0[0]): (-0.99, 0.79),    # (lowB, upB)
             (GC.T_NM_GT1[0], GC.T_NM_GT1[0]): (-0.99, 0.8),
             (GC.T_NM_GT2[0], GC.T_NM_GT2[0]): (-0.99, 0.9),
             (GC.T_NM_GT3[0], GC.T_NM_GT3[0]): (-0.99, 0.9),
             (GC.T_NM_GT4[0], GC.T_NM_GT4[0]): (-0.99, 0.9),
             (GC.T_NM_GT5[0], GC.T_NM_GT5[0]): (-0.99, 0.8)}

# --- deviations between features ---------------------------------------------
dPosCIBnd = {'C1': (1., 0.25),    # positive concordance index bound (x SD)
             'C2': (2., 1.),      # key: string of class
             'C3': (3., 1.5),     # value: (bound, weight)
             'C4': (5., 1.75),
             'C5': (10., 2.)}
dWtDv = {'Min': 100, 'Max': 1}

# --- correlation and deviations between features -----------------------------
dNTopCrDv = {('MetD', 'MetD'): [10, 100, 500],
             ('MetDC', 'MetDC'): [10, 100, 500],
             ('MetDG', 'MetDG'): [10, 100, 500],

             ('MetD', 'PhoD'): [10, 100, 1000, 10000],
             ('MetDC', 'PhoDC'): [10, 100, 1000, 10000],
             ('MetDG', 'PhoDG'): [10, 100, 1000, 10000],

             ('PhoD', 'PhoD'): [10, 100, 1000, 10000, 100000, 1000000],
             ('PhoDC', 'PhoDC'): [10, 100, 1000, 10000, 100000, 1000000],
             ('PhoDG', 'PhoDG'): [10, 100, 1000, 10000, 100000, 1000000],

             ('CombDC', 'CombDC'): [10, 100, 1000, 10000, 100000, 1000000],
             ('CombDG', 'CombDG'): [10, 100, 1000, 10000, 100000, 1000000],

             ('MetDCX', 'MetDCX'): [10, 100, 500],
             ('MetDGX', 'MetDGX'): [10, 100, 500],

             ('MetDCX', 'PhoDCX'): [10, 100, 1000, 10000],
             ('MetDGX', 'PhoDGX'): [10, 100, 1000, 10000],

             ('PhoDCX', 'PhoDCX'): [10, 100, 1000, 10000, 100000, 1000000],
             ('PhoDGX', 'PhoDGX'): [10, 100, 1000, 10000, 100000, 1000000],

             ('CombDCX', 'CombDCX'): [10, 100, 1000, 10000, 100000, 1000000],
             ('CombDGX', 'CombDGX'): [10, 100, 1000, 10000, 100000, 1000000],

             ('ClRMetD', 'ClRMetD'): [3, 5, 10, 20],
             ('ClRMetDC', 'ClRMetDC'): [3, 5, 10, 20],
             ('ClRMetDG', 'ClRMetDG'): [3, 5, 10, 20],

             ('ClRMetD', 'ClRPhoD'): [3, 5, 10, 20],
             ('ClRMetDC', 'ClRPhoDC'): [3, 5, 10, 20],
             ('ClRMetDG', 'ClRPhoDG'): [3, 5, 10, 20],

             ('ClRPhoD', 'ClRPhoD'): [3, 5, 10, 20],
             ('ClRPhoDC', 'ClRPhoDC'): [3, 5, 10, 20],
             ('ClRPhoDG', 'ClRPhoDG'): [3, 5, 10, 20],

             ('ClRCombDC', 'ClRCombDC'): [3, 5, 10, 20],
             ('ClRCombDG', 'ClRCombDG'): [3, 5, 10, 20],

             ('ClRMetDCX', 'ClRMetDCX'): [3, 5, 10, 20],
             ('ClRMetDGX', 'ClRMetDGX'): [3, 5, 10, 20],

             ('ClRMetDCX', 'ClRPhoDCX'): [3, 5, 10, 20],
             ('ClRMetDGX', 'ClRPhoDGX'): [3, 5, 10, 20],

             ('ClRPhoDCX', 'ClRPhoDCX'): [3, 5, 10, 20],
             ('ClRPhoDGX', 'ClRPhoDGX'): [3, 5, 10, 20],

             ('ClRCombDCX', 'ClRCombDCX'): [3, 5, 10, 20],
             ('ClRCombDGX', 'ClRCombDGX'): [3, 5, 10, 20]}

# --- top correlations / deviations -------------------------------------------
lCalcTACD = [((GC.S_MET_D, GC.S_PHO_D), False),     # (columns to calc. TACD,
             ((GC.S_MET_D, GC.S_BIN_C_2), True)]    #  do weighting (bool))

# --- data specific -----------------------------------------------------------
dSrtF = {GC.S_IDX: {sAv + sPos + sCorrV: {'Asc': False}}}
dSrtR = {GC.S_IDX: {sAv + sPos + sCorrV: {'Asc': False}}}

# --- results -----------------------------------------------------------------
cSep = ';'    # separator for results files

# --- names and paths of files and dirs ---------------------------------------
dISPr = {'No': {0: ''}, 'Tr': {0: ''}, 'Dv': {0: ''}}
nmHist = 'Hist'         # name of histogram plot
nmCorrF = sCorrS        # name of correlation result file
nmTACDRF = 'TopAv' + sCorrS    # name of top / average correlation result file
nmStrong = 'Strong'     # name prefix of strong correlation plot

# --- rounding ----------------------------------------------------------------
nDigRndBO = 4
nDigRndDv = 1
nDigRndDvSc = 2
nDigRndCI = 4

# --- graphics parameters -----------------------------------------------------
nBins = 100             # number of bins
histAlpha = 0.4         # alpha for histogram
histXLim = (-1, 1)      # x-limits for histogram
corrLim = (-1, 1)       # lower and upper limit for correlation
compToNormDist = True   # compare to normal distribution?

# --- derived values ----------------------------------------------------------
lSAvCorrV = [sAv + sNeg + sCorrV, sAv + sPos + sCorrV]
lSAvSpearV = [sAv + sNeg + sSpearV, sAv + sPos + sSpearV]
lSAvKendV = [sAv + sNeg + sKendV, sAv + sPos + sKendV]
lSAvDvSc = [sAv + sNeg + sDvSc, sAv + sPos + sDvSc]
lSAvCI = [sAv + sNeg + sCI, sAv + sPos + sCI]
lSTopCorrV = [sTop + sNeg + sCorrV, sTop + sPos + sCorrV]
lSTopSpearV = [sTop + sNeg + sSpearV, sTop + sPos + sSpearV]
lSTopKendV = [sTop + sNeg + sKendV, sTop + sPos + sKendV]
lSTopDvSc = [sTop + sNeg + sDvSc, sTop + sPos + sDvSc]
lSTopCI = [sTop + sNeg + sCI, sTop + sPos + sCI]
lSCorrVals = [sCorrV, sSpearV, sKendV]
lSCorrVals2x = [sCorrV]*2 + [sSpearV]*2 + [sKendV]*2
lSCorrAll = lSCorrVals + [sCorrP, sSpearP, sKendP]
lSAvAll = lSAvCorrV + lSAvSpearV + lSAvKendV
lSTopAll = lSTopCorrV + lSTopSpearV + lSTopKendV
lPosCIBnd = [t[0] for t in dPosCIBnd.values()]
lCIWts = [t[1] for t in dPosCIBnd.values()]

# --- assertions --------------------------------------------------------------
for tGT, tSCorrBnd in dSCorrBnd.items():
    assert len(tGT) >= 2 and len(tSCorrBnd) >= 2

for t, bWt in lCalcTACD:
    assert type(bWt) == bool
    assert len(t) == 2
    assert ((t[0] == GC.S_MET_D and t[1] in [GC.S_PHO_D] + GC.L_S_PHO_CL) or
            (t[0] in [GC.S_PHO_D] + GC.L_S_PHO_CL and t[1] == GC.S_MET_D))

# --- create input dictionary -------------------------------------------------
dIO = {# --- general and names
       'calcDictR': calcDictR,
       'calcForTrD': calcForTrD,
       'sNmSpec': sNmSpec,
       'sCt': sCt,
       'sNeg': sNeg,
       'sPos': sPos,
       'sAv': sAv,
       'sTop': sTop,
       'sCorrS': sCorrS,
       'sCorrL': sCorrL,
       'sCorrV': sCorrV,
       'sSpearV': sSpearV,
       'sKendV': sKendV,
       'sCorrP': sCorrP,
       'sSpearP': sSpearP,
       'sKendP': sKendP,
       'sDvSc': sDvSc,
       'sCI': sCI,
       'sOccCI': sOccCI,
       'dSDvSc': dSDvSc,
       'dSDvScX': dSDvScX,
       'dSCI': dSCI,
       'dSCIX': dSCIX,
       # --- correlation
       'dPltSCorr': dPltSCorr,
       'dSCorrBnd': dSCorrBnd,
       # --- deviations between features
       'dPosCIBnd': dPosCIBnd,
       'dWtDv': dWtDv,
       # --- correlation and deviations between features
       'dNTopCrDv': dNTopCrDv,
       # --- top correlations / deviations
       'lCalcTACD': lCalcTACD,
       # --- data specific
       'dSrtF': dSrtF,
       'dSrtR': dSrtR,
       # --- results
       'cSep': cSep,
       # --- names and paths of files and dirs
       'dISPr': dISPr,
       'nmHist': nmHist,
       'nmCorrF': nmCorrF,
       'nmTACDRF': nmTACDRF,
       'nmStrong': nmStrong,
       # --- rounding
       'nDigRndBO': nDigRndBO,
       'nDigRndDv': nDigRndDv,
       'nDigRndDvSc': nDigRndDvSc,
       'nDigRndCI': nDigRndCI,
       # --- graphics parameters
       'nBins': nBins,
       'histAlpha': histAlpha,
       'histXLim': histXLim,
       'corrLim': corrLim,
       'cmpND': compToNormDist,
       # --- derived values
       'lSAvCorrV': lSAvCorrV,
       'lSAvSpearV': lSAvSpearV,
       'lSAvKendV': lSAvKendV,
       'lSAvDvSc': lSAvDvSc,
       'lSAvCI': lSAvCI,
       'lSTopCorrV': lSTopCorrV,
       'lSTopSpearV': lSTopSpearV,
       'lSTopKendV': lSTopKendV,
       'lSTopDvSc': lSTopDvSc,
       'lSTopCI': lSTopCI,
       'lSCorrVals': lSCorrVals,
       'lSCorrVals2x': lSCorrVals2x,
       'lSCorrAll': lSCorrAll,
       'lSAvAll': lSAvAll,
       'lSTopAll': lSTopAll,
       'lPosCIBnd': lPosCIBnd,
       'lCIWts': lCIWts}

###############################################################################

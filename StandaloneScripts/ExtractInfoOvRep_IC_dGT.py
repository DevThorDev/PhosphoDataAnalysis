# -*- coding: utf-8 -*-
###############################################################################
# --- ExtractInfoOvRep_IC_dGT.py ----------------------------------------------
###############################################################################
import os, time

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# --- CONSTANTS ---------------------------------------------------------------
S_DASH = '-'
S_USC = '_'
S_BAR = '|'
S_DOT = '.'
S_CSV = 'csv'
S_PDF = 'pdf'
S_D = 'd'
S_M = 'm'
S_P = 'p'
S_N = 'N'
S_Y = 'Y'
S_NO = 'No'
L_N = [S_N]
L_Y = [S_Y]
L_NY = [S_N, S_Y]

S_0 = '0'
S_1 = '1'
S_2 = '2'
S_3 = '3'
S_4 = '4'
S_5 = '5'

S_10 = S_1 + S_0
S_50 = S_5 + S_0
S_51 = S_5 + S_1

S_SRT = 'Sort'
S_THR = 'Threshold'
S_SEL = 'Selection'
S_ALL = 'All'
S_SGL = 'Sgl'

S_MET = 'Metabolite'
S_PHO = 'Phosphopeptide'
L_S_M_P = [S_MET, S_PHO]

S_BIN_L = 'BinCode'
S_BIN_L_2 = S_BIN_L + S_2
S_SELECTED = 'Selected'

S_GT = 'GT'
S_GT0 = S_GT + S_0
S_GT1 = S_GT + S_1
S_GT5 = S_GT + S_5
L_S_GT = [S_GT0, S_GT1, S_GT5]

S_D_GT_M = S_D + S_GT + 'M'
S_D_GT_P = S_D + S_GT + 'P'
L_S_D_GT = [S_D_GT_M, S_D_GT_P]

L_S_FT = ['DR', 'DS', 'NR', 'NS']
L_S_FT_CHG = [L_S_FT[i] + S_BAR + L_S_FT[j] for j in range(len(L_S_FT))
              for i in range(len(L_S_FT)) if i != j]

D_HD_C_PA = {sMP: {sGT: [S_USC.join([sGT, sFt, sMP[0]]) for sFt in L_S_FT]
                   for sGT in L_S_GT} for sMP in L_S_M_P}

S_MIN = 'min'
S_MAX = 'max'
L_S_MIN_MAX = [S_MIN, S_MAX]

S_IC = 'IC'
S_IC_P = S_IC + S_USC + 'P'
S_IC_N = S_IC + S_USC + 'N'
S_PEAR_CR = 'PearsonCorr'
S_SPEAR_CR = 'SpearmanCorr'
S_PEAR_P = 'Pearson_pVal'
S_SPEAR_P = 'Spearman_pVal'
S_IC_P_S = S_IC + 'P'
S_IC_N_S = S_IC + 'N'
S_PEAR_CR_S = 'PeC'
S_SPEAR_CR_S = 'SpC'
S_PEAR_P_S = 'Pep'
S_SPEAR_P_S = 'Spp'

S_D_CL = 'd' + S_USC + 'Cl'
S_D_E = 'd' + S_USC + 'E'
S_D_GT = S_D + S_USC + S_GT
S_D_CL_S = 'dCl'
S_D_E_S = 'dE'
S_D_GT_S = S_D + S_GT
L_S_D_2GT = [S_D + S_USC*2 + s1 + S_USC + s2 for s1 in L_S_GT for s2 in L_S_GT
             if int(s1[-1]) > int(s2[-1])]
S_ICC = S_IC + 'C'
S_ICC_LEG = 'Concordance index component'
D_HD_C_ICC = {sMPI: [S_USC.join([sFtChg, sMPI[0]]) for sFtChg in L_S_FT_CHG]
              for sMPI in L_S_M_P + [S_ICC]}

S_SRT_BY = S_SRT + 'edBy'
S_ORD = 'Order'
S_ASC = 'Asc'
S_DSC = 'Dsc'
S_YLBL_PAT_PLT = 'Pattern (z-score)'
S_YLBL_ICC_PLT = 'Change (xSD) / IC component'

S_IC_GT0 = S_USC.join([S_IC, S_GT0])
S_IC_GT1 = S_USC.join([S_IC, S_GT1])
S_IC_GT5 = S_USC.join([S_IC, S_GT5])
L_S_IC_GT = [S_IC_GT0, S_IC_GT1, S_IC_GT5]

S_SG = 'S5'
S_SG_M = 'M' + S_SG
S_SG_P = 'P' + S_SG
S_SG_MP = 'MP' + S_SG
S_SG_MET = 'MetSig5'
S_SG_PHO = 'PhoSig5'
L_S_SG_S = [S_SG_M, S_SG_P]
L_S_SG = [S_SG_MET, S_SG_PHO]
L_S_SG_MET_GT = [S_USC.join([S_SG_MET, sGT]) for sGT in L_S_GT]
L_S_SG_PHO_GT = [S_USC.join([S_SG_PHO, sGT]) for sGT in L_S_GT]
L_S_SG_GT = [S_USC.join([sSg, sGT]) for sGT in L_S_GT for sSg in L_S_SG]

S_SB = 'SB'
S_SB_P = 'P' + S_SB
L_S_SB_GT = [S_USC.join([sSB, sGT]) for sGT in L_S_GT for sSB in L_S_SG]
L_S_SB_S = [S_SB_P]

S_BASE_CL = 'BaseClass'
S_INP_DATA = 'InputData'
S_ROOT_CL = 'RootClass'
S_F_NM_CMP = 'FNmCmp'
S_F_NM_CMP_DGT = S_F_NM_CMP + 'DGT'
S_F_NM_CMP_IC_ALL_GT = S_F_NM_CMP + S_IC + S_ALL + S_GT
S_F_NM_CMP_IC_SGL_GT = S_F_NM_CMP + S_IC + S_SGL + S_GT
S_F_NM_CNSTR = 'FileNameConstructor'
S_EXTR_INFO = 'ExtrInfo'
S_PLTR = 'Plotter'
S_PAT_PLTR = 'PatternPlotter'
S_ICC_PLTR = 'ICCmpPlotter'

S_RMNG_COL1_IC = 'PearsonCorr'
S_NEW_IDX = 'NewIndex'
# L_S_NO_GT = L_S_M_P + ['Protein', 'BinCode', 'BinCode2', 'MapMan',
#                        'Description', 'SelBinCode']
L_S_NO_GT = L_S_M_P + []
L_S_ADD_GT = [S_RMNG_COL1_IC, 'SpearmanCorr', 'Pearson_pVal', 'Spearman_pVal',
              'IC_N', 'IC_P', 'IC', 'MetSig5', 'PhoSig5']

S_NM_PAT_PLT = 'PatternPlot'
S_NM_ICC_PLT = 'ICCmpPlot'

R04 = 4

# --- INPUT -------------------------------------------------------------------
# --- flow control ------------------------------------------------------------
doInfoExtr = False               # True / False
doPlotPat = True                # True / False
doPlotICC = True                # True / False
lSpecSel = ['S', 'F']           # list of column selections: 'S'hort / 'F'ull

# --- general input -----------------------------------------------------------
modDisp = 10000

# --- data specific input -----------------------------------------------------
sortInFNm = False               # sorting in file name? (True / False)

dUsedK = {S_IC: S_IC,           # key (col. hdr.) for the IC file
          S_D_GT_M: S_D_GT,     # key (col. hdr.) for the dGTM file
          S_D_GT_P: S_D_GT}     # key (col. hdr.) for the dGTP file
# dUsedK = {S_IC: S_IC,           # key (col. hdr.) for the IC file
#           S_D_GT_M: S_10,     # key (col. hdr.) for the dGTM file
#           S_D_GT_P: S_10}     # key (col. hdr.) for the dGTP file

dISort = {S_IC: {S_GT0: {S_SRT_BY: dUsedK[S_IC], S_ORD: S_DSC},
                 S_GT1: {S_SRT_BY: dUsedK[S_IC], S_ORD: S_DSC},
                 S_GT5: {S_SRT_BY: dUsedK[S_IC], S_ORD: S_DSC}},
          S_D_GT_M: {S_SRT_BY: dUsedK[S_D_GT_M], S_ORD: S_DSC},
          S_D_GT_P: {S_SRT_BY: dUsedK[S_D_GT_P], S_ORD: S_DSC}}

dThr = {S_IC: {S_GT0: {S_MIN: 6.0, S_MAX: None},
                S_GT1: {S_MIN: 6.0, S_MAX: None},
                S_GT5: {S_MIN: 6.0, S_MAX: None}},
        S_D_GT_M: {S_MIN: None, S_MAX: None},
        S_D_GT_P: {S_MIN: None, S_MAX: None}}

lSelSGM, lSelSGP = L_Y, L_Y
# lSelSB = L_NY                    # L_NY / L_Y (sel. bins only) / L_N
lSelSB = L_Y                    # L_NY / L_Y (sel. bins only) / L_N
dSel = {(S_SG, S_SG_MP): {S_SG_M: {S_GT0: lSelSGM, S_GT1: lSelSGM,
                                   S_GT5: lSelSGM},
                          S_SG_P: {S_GT0: lSelSGP, S_GT1: lSelSGP,
                                   S_GT5: lSelSGP}},
        (S_SB, S_SB_P): {S_SB_P: {S_GT0: lSelSB, S_GT1: lSelSB,
                                  S_GT5: lSelSB}}}

sSep = ';'

# --- graphics parameters / all plots -----------------------------------------
szFontLeg = 'small'             # font size of legend
nCharDsp = 60                   # number of chars displayed for legend item
coordAnchorBox = (0.5, 1.02)    # coordinates of the legend anchor box
lWdPlt = 0.75                   # line width in plot

# --- graphics parameters / pattern plot --------------------------------------
# dPairsPaP = {(('Leu_STTTTV'), (S_GT0, S_GT1, S_GT5)):
#              ('Leucine', 'STTTTVS(0.003)S(0.996)VHS(0.001)PTTDQDFSK')}
# dPairsPaP = {(('Tetra_S(0.001)AS(0.749)T(0.251)P'), (S_GT0, S_GT1, S_GT5)):
#               ('Tetradecanoic_acid', 'S(0.001)AS(0.749)T(0.251)PLLNSLVHVS(0.179)S(0.821)PRDS(1)PIETVESVHQIQR'),
#               (('Beta_ADKTDII'), (S_GT0, S_GT1, S_GT5)):
#               ('Beta-alanine', 'ADKTDIIS(0.607)S(0.117)S(0.12)S(0.156)DKAS(1)PPPPSAFR'),
#               (('Hexa_S(0.001)AS(0.749)T(0.251)P'), (S_GT0, S_GT1, S_GT5)):
#               ('Hexadecanoic_acid', 'S(0.001)AS(0.749)T(0.251)PLLNSLVHVS(0.179)S(0.821)PRDS(1)PIETVESVHQIQR')}
# dPairsPaP = {(('Isoleu_DLDVNE'), (S_GT0, S_GT1, S_GT5)):
#               ('Isoleucine', 'DLDVNES(1)GPPAAR'),
#               (('Val_DLDVNE'), (S_GT0, S_GT1, S_GT5)):
#               ('Valine', 'DLDVNES(1)GPPAAR'),
#               (('Aspart_SDKPLNY'), (S_GT0, S_GT1, S_GT5)):
#               ('Aspartic_acid', 'SDKPLNYS(1)PDPENESGINER'),
#               (('Aspart_DLDVNE'), (S_GT0, S_GT1, S_GT5)):
#               ('Aspartic_acid', 'DLDVNES(1)GPPAAR'),
#               (('Malic_DLDVNE'), (S_GT0, S_GT1, S_GT5)):
#               ('Malic_acid', 'DLDVNES(1)GPPAAR')}
dPairsPaP = {(('Tetradeca_KTSSST'), (S_GT0, S_GT1, S_GT5)):
              ('Tetradecanoic_acid', 'KTSSSTISTNPS(0.001)S(0.998)PIS(0.001)TASTGKPPLPR'),
              (('Isoleu_DLDVNE'), (S_GT0, S_GT1, S_GT5)):
              ('Isoleucine', 'DLDVNES(1)GPPAAR')}

nmPaP = S_NM_PAT_PLT            # name prefix of the pattern plot

# --- graphics parameters / IC component plot ---------------------------------
dPairsICP = dPairsPaP
nmICP = S_NM_ICC_PLT            # name prefix of the IC component plot
wdthBar = 0.2                   # width of single bars
wdthGrp = 0.75                   # width of bar group
degRotXLbl = 90                 # degree rotation of x-labels

# --- names and paths of files and dirs ---------------------------------------
sFIn_IC_M_P = 'IC_Met_Pho'
sFIn_dGT_M = 'DistGT_Met'
sFIn_dGT_P = 'DistGT_Pho'

sFOutS = 'S_XIOvRep'
sFOutF = 'F_XIOvRep'

sFIn_PaP = 'F_XIOvRep_IC_IC_All_7p25_No_MPS5Y_PSBY_dGTM_dGT_No_No_dGTP_dGT_No_No'
dSFIn_ICP = {sGT: ('ICCmp__BinOp_MetD_DvSD_' + sGT + '_AllD_PhoD_DvSD_'
                   + sGT + '_AllD') for sGT in L_S_GT}
sFOutPaP = nmPaP
sFOutICP = nmICP

sDirInCSV = '51_Inp_IC_DGT'
sDirOutCSV = '52_OutCSV_IC_DGT'
sDirOutPaP = '55_OutPDF_IC_DGT'
sDirOutICP = sDirOutPaP

pBaseIn = os.path.join('..', '..', '..', '25_Papers', '01_FirstAuthor',
                       '04_SysBio_DataAnalysis')
# pBaseOut = os.path.join('..', '..', '..', '25_Papers', '01_FirstAuthor',
#                         '04_SysBio_DataAnalysis')
pBaseOut = os.path.join('..', '..', '..', '..', '..', '..', 'W')
pInCSV = os.path.join(pBaseIn, sDirInCSV)
pOutCSV = os.path.join(pBaseOut, sDirOutCSV)

pInPaP = os.path.join(pBaseOut, sDirOutCSV)
pInICP = os.path.join(pBaseIn, sDirInCSV)
pOutPaP = os.path.join(pBaseOut, sDirOutPaP)
pOutICP = os.path.join(pBaseOut, sDirOutICP)

# --- derived values ----------------------------------------------------------
dMapCHdSel = {S_SG: {S_SG_M: {sGT: L_S_SG_MET_GT[k]
                              for k, sGT in enumerate(L_S_GT)},
                     S_SG_P: {sGT: L_S_SG_PHO_GT[k]
                              for k, sGT in enumerate(L_S_GT)}},
              S_SB: {S_SB_P: {sGT: S_USC.join([S_SELECTED, S_PHO[0]])
                              for sGT in L_S_GT}}}
dComprStr = {S_IC_P: S_IC_P_S,
             S_IC_N: S_IC_N_S,
             S_PEAR_CR: S_PEAR_CR_S,
             S_SPEAR_CR: S_SPEAR_CR_S,
             S_PEAR_P: S_PEAR_P_S,
             S_SPEAR_P: S_SPEAR_P_S,
             S_D_CL: S_D_CL_S,
             S_D_E: S_D_E_S,
             S_D_GT: S_D_GT_S}

# pFIGT0 = os.path.join(pInCSV, sFIn_IC_M_P + S_USC + S_GT0 + S_DOT + S_CSV)
# pFIGT1 = os.path.join(pInCSV, sFIn_IC_M_P + S_USC + S_GT1 + S_DOT + S_CSV)
# pFIGT5 = os.path.join(pInCSV, sFIn_IC_M_P + S_USC + S_GT5 + S_DOT + S_CSV)
# dPFInIC = {S_GT0: pFIGT0, S_GT1: pFIGT1, S_GT5: pFIGT5}
dPFInIC = {sGT: os.path.join(pInCSV, sFIn_IC_M_P + S_USC + sGT + S_DOT + S_CSV)
           for sGT in L_S_GT}

# --- assertions --------------------------------------------------------------
assert set(dISort) == set(dThr)
for cD in [dISort, dThr]:
    assert S_IC in cD and set(cD[S_IC]) == set(L_S_GT)

# --- INPUT DICTIONARY --------------------------------------------------------
dInput = {# --- constants
          'sBC_L': S_BIN_L,
          'sBC2_L': S_BIN_L_2,
          'lSGT': L_S_GT,
          'sBase': S_BASE_CL,
          'sInpDat': S_INP_DATA,
          'sRoot': S_ROOT_CL,
          'sFNmCmp': S_F_NM_CMP,
          'sFNmCmpDGT': S_F_NM_CMP_DGT,
          'sFNmCmpICAllGT': S_F_NM_CMP_IC_ALL_GT,
          'sFNmCmpICSglGT': S_F_NM_CMP_IC_SGL_GT,
          'sFNmCnstr': S_F_NM_CNSTR,
          'sExtrInfo': S_EXTR_INFO,
          'sPltr': S_PLTR,
          'sPatPltr': S_PAT_PLTR,
          'sICCPltr': S_ICC_PLTR,
          'sMet': S_MET,
          'sPho': S_PHO,
          'R04': R04,
          # --- flow control
          'doInfoExtr': doInfoExtr,
          'doPlotPat': doPlotPat,
          'doPlotICC': doPlotICC,
          'lSpecSel': lSpecSel,
          # --- general input
          'modDisp': modDisp,
          # --- data specific input
          'sortInFNm': sortInFNm,
          'dUsedK': dUsedK,
          'dISort': dISort,
          'dThr': dThr,
          'dSel': dSel,
          'sSep': sSep,
          # --- graphics parameters / pattern plot
          'plotOfPatterns': {'dPairsPaP': dPairsPaP,
                             'nmPaP': nmPaP,
                             'szFontLeg': szFontLeg,
                             'nCharDsp': nCharDsp,
                             'coordAnchorBox': coordAnchorBox,
                             'lWdPlt': lWdPlt},
          # --- graphics parameters / IC component plot
          'plotOfICCmp': {'dPairsICP': dPairsICP,
                          'nmICP': nmICP,
                          'szFontLeg': szFontLeg,
                          'nCharDsp': nCharDsp,
                          'coordAnchorBox': coordAnchorBox,
                          'lWdPlt': lWdPlt,
                          'wdthBar': wdthBar,
                          'wdthGrp': wdthGrp,
                          'degRotXLbl': degRotXLbl},
          # --- names and paths of files and dirs
          'pInCSV': pInCSV,
          'pOutCSV': pOutCSV,
          'pInPaP': pInPaP,
          'pInICP': pInICP,
          'pOutPaP': pOutPaP,
          'pOutICP': pOutICP,
          'dPFInIC': dPFInIC,
          'pFInM': os.path.join(pInCSV, sFIn_dGT_M + S_DOT + S_CSV),
          'pFInP': os.path.join(pInCSV, sFIn_dGT_P + S_DOT + S_CSV),
          'sFOutS': sFOutS,
          'sFOutF': sFOutF,
          'pFInPaP': os.path.join(pInPaP, sFIn_PaP + S_DOT + S_CSV),
          'dPFInICP': {sGT: os.path.join(pInICP, dSFIn_ICP[sGT] + S_DOT +
                                         S_CSV) for sGT in L_S_GT},
          'sFOutPaP': sFOutPaP,
          'sFOutICP': sFOutICP,
          # --- further derived values
          'dMapCHdSel': dMapCHdSel,
          'dComprStr': dComprStr}

# --- FUNCTIONS ---------------------------------------------------------------
def addToDictD(cD, cKMain, cKSub, cV):
    if cKMain in cD:
        assert cKSub not in cD[cKMain]
        cD[cKMain][cKSub] = cV
    else:
        cD[cKMain] = {}
        cD[cKMain][cKSub] = cV

def allElEq(iterator):
    iterator = iter(iterator)
    try:
        first = next(iterator)
    except StopIteration:
        return True
    return all(first == x for x in iterator)

def addIt(s='', sAdd='', sSep=S_USC):
    if len(s) > 0:
        if s[-1] != sSep:
            s += sSep
    return s + sAdd

def comprStr(s, dCmpr, sSep=S_USC):
    lSSplO, lSSplN = s.split(sSep), []
    for i, s1s in enumerate(lSSplO):
        s2s = s1s
        if i > 0:
            s2s = sSep.join([lSSplO[i - 1], s1s])
        if s2s in dCmpr:
            lSSplN[-1] = dCmpr[s2s]
        elif s1s in dCmpr:
            lSSplN.append(dCmpr[s1s])
        else:
            lSSplN.append(s1s)
    return sSep.join(lSSplN)

def num2StrF(cV, dRepl={S_DOT: S_P, S_DASH: S_M, 'None': S_NO}):
    cS = str(cV)
    for sK, sV in dRepl.items():
        cS = cS.replace(sK, sV)
    return cS

def getSHdCD2GT(sID, lSDef):
    for sF in lSDef:
        lSSpl = sF.split(S_USC)
        assert len(lSSpl) >= 2
        if lSSpl[-2][-1] + lSSpl[-1][-1] == sID:
            return sF
    return None

def getLVals(pdDfr, sHd, nEl=1):
    lVCC, lVRICC = pdDfr.loc[:, sHd].to_list(), list(range(1, nEl + 1))
    minCol = min(lVCC)
    iMin = lVCC.index(minCol)
    lVRICC[iMin:] = [iMin + 1]*(len(lVCC) - iMin)
    assert min(lVRICC) >= 1
    maxVRICC = max(lVRICC)
    for i in range(nEl):
        lVRICC[i] /= maxVRICC
    return lVRICC

def transcrDict2Dfr(cD, cDfr, lSHdrIni):
    cDfrRes = cDfr.loc[:, lSHdrIni]
    for cDSub in cD.values():
        for sHdr in cDSub.values():
            cDfrRes[sHdr] = cDfr.loc[:, sHdr]
    return cDfrRes

def sortDfr(pdDfr, dSrt, sSrtBy, sOrd, srtKind='stable'):
    isAsc, sHdC = (dSrt[sOrd] == S_ASC), dSrt[sSrtBy]
    if sHdC not in pdDfr.columns:
        sHdC = getSHdCD2GT(sHdC, L_S_D_2GT)
    pdDfr.sort_values(by=sHdC, ascending=isAsc, inplace=True, kind=srtKind)

def thrFilter(pdDfr, sHdC, thrMin, thrMax):
    if sHdC not in pdDfr.columns:
        sHdC = getSHdCD2GT(sHdC, L_S_D_2GT)
    if thrMin is None:
        if thrMax is None:
            return pdDfr
        else:
            return pdDfr[pdDfr[sHdC] <= thrMax]
    else:
        if thrMax is None:
            return pdDfr[pdDfr[sHdC] >= thrMin]
        else:
            return pdDfr[(pdDfr[sHdC] >= thrMin) & (pdDfr[sHdC] <= thrMax)]

def filterDfrSel(pdDfr, dGT, dMap, sKMain, sKSub, sGT):
    sHdC, lSSel = dMap[sKMain][sKSub][sGT], dGT[sGT]
    if sHdC in pdDfr.columns:
        return pdDfr[pdDfr[sHdC].isin(lSSel)]
    else:
        return pdDfr

def applySelFilter(pdDfr, dSel, dMap):
    for tK, dSub in dSel.items():
        for sKSub, dGT in dSub.items():
            if allElEq(dGT.values()) and allElEq(dMap[tK[0]][sKSub].values()):
                pdDfr = filterDfrSel(pdDfr, dGT, dMap, tK[0], sKSub, L_S_GT[0])
            else:
                for sGT in dGT:
                    pdDfr = filterDfrSel(pdDfr, dGT, dMap, tK[0], sKSub, sGT)
    return pdDfr

def getLNewIdx(pdDfr):
    assert (S_MET in pdDfr.columns) and (S_PHO in pdDfr.columns)
    return [S_USC.join([pdDfr.at[i, S_MET], pdDfr.at[i, S_PHO]])
            for i in pdDfr.index]

def modifyDfr(pdDfr, sGT):
    lHdCN = [s + S_USC + sGT for s in pdDfr.columns if s not in L_S_NO_GT]
    pdDfr.columns = L_S_NO_GT + lHdCN
    dfrMod = pd.concat([pd.Series(getLNewIdx(pdDfr), name=S_NEW_IDX),
                        pdDfr.reset_index(drop=True)],
                       axis=1, verify_integrity=True)
    return dfrMod.set_index(S_NEW_IDX, verify_integrity=True)

def concDfr(dDfr):
    dfrM = modifyDfr(dDfr[L_S_GT[0]], L_S_GT[0])
    for sGT in L_S_GT[1:]:
        dfrMR = modifyDfr(dDfr[sGT], sGT).drop(columns=L_S_NO_GT)
        dfrM = pd.concat([dfrM, dfrMR], axis=1, verify_integrity=True)
    return dfrM.drop(columns=L_S_NO_GT)

def appendToDDat(lDat, cDfrFl, sI, sHdC):
    cV = np.nan
    if sI in cDfrFl.index:
        cV = cDfrFl.at[sI, sHdC]
    lDat.append(cV)

def checkAllSetsEqSrtThr(dSTIC):
    lSets = [set(dSTIC[sGT].values()) for sGT in L_S_GT]
    return set.intersection(*lSets) == set.union(*lSets)

def checkMPSetsEqSl(dSlS, sGT):
    lSets = [set(dSlS[sKSlS][sGT]) for sKSlS in dSlS]
    return set.intersection(*lSets) == set.union(*lSets)

def checkAllSetsEqSl(dSlS):
    lSets = [set(dSlS[sKSlS][sGT]) for sKSlS in dSlS for sGT in L_S_GT]
    return set.intersection(*lSets) == set.union(*lSets)

def getSel(setSel, sIDSl):
    if setSel == {S_N, S_Y}:
        return sIDSl + S_N + S_Y
    elif  setSel == {S_Y}:
        return sIDSl + S_Y
    elif  setSel == {S_N}:
        return sIDSl + S_N
    return ''

def addSelAll(dSl, tKSl, sGT=S_GT0):
    sSubKSl, sIDSl = list(dSl[tKSl])[0], tKSl[1]
    return getSel(set(dSl[tKSl][sSubKSl][sGT]), sIDSl)

def saveDfrRes(dfrRes, dDat, pFOut, sSep, dSel=None, dMap=None):
    dfrRes = dfrRes.append(pd.DataFrame(dDat), ignore_index=True,
                           verify_integrity=True)
    if dSel is not None and dMap is not None:
        dfrRes = applySelFilter(dfrRes, dSel, dMap)
    dfrRes.dropna(axis=0, how='any', inplace=True)
    dfrRes.reset_index(drop=True).to_csv(pFOut, sep=sSep)
    return dfrRes

def decorateClosePlot(cFig, cAx, dPlt, pPltF, sYLbl=''):
    cAx.set_ylabel(sYLbl)
    l = cAx.legend(loc='lower center', bbox_to_anchor=dPlt['coordAnchorBox'],
                   fontsize=dPlt['szFontLeg'])
    if l is not None:
        cFig.savefig(pPltF, bbox_extra_artists=(l,), bbox_inches='tight')
    else:
        cFig.savefig(pPltF)
    plt.close()

def printElapsedTimeSim(stT, cT, sPre = 'Time'):
    # calculate and display elapsed time
    elT = round(cT - stT, R04)
    print(sPre, 'elapsed:', elT, 'seconds, this is', round(elT/60, R04),
          'minutes or', round(elT/3600, R04), 'hours or',
          round(elT/(3600*24), R04), 'days.')

# --- CLASSES -----------------------------------------------------------------
class BaseClass():
    def __init__(self):
        self.idO = S_BASE_CL
        self.descO = 'Base class'
        print('Initiated "Base" base object.')

    def printIDDesc(self):
        print('Object ID:', self.idO)
        print('Object description:', self.descO)

    def printAttrList(self):
        lAttr = dir(self)
        print('List of attributes:')
        for cAttr in lAttr:
            print(cAttr)

    def printAttrData(self):
        print('Attributes and attribute values:')
        d = vars(self)
        for cK, cV in d.items():
            print(cK, ':\t', cV)

class InputData(BaseClass):
    def __init__(self, dInp):
        super().__init__()
        self.idO = dInp['sInpDat']
        self.descO = 'Input data class'
        self.dI = dInp
        self.fillInp()
        print('Initiated "InputData" base object.')

    def fillInp(self):
        for sK, cV in self.dI.items():
            setattr(self, sK, cV)
        print('Set InputData attributes.')

class RootClass(BaseClass):
    def __init__(self, InpD):
        super().__init__()
        self.idO = InpD.sRoot
        self.descO = 'Root class'
        self.inpD = InpD
        self.dUsedK = self.inpD.dUsedK
        self.dSort = self.inpD.dISort
        self.dT = self.inpD.dThr
        self.dSl = self.inpD.dSel
        self.dfrIn, self.dDfrIn = None, None
        print('Initiated "RootClass" base object.')

    def printObjInfo(self):
        print('-'*20, 'Object', self.descO, '(ID', self.idO, ')', '-'*20)
        print('-'*8, 'Input data:')
        self.inpD.printAttrData()
        print('-'*8, 'Attributes of', self.descO, 'class:')
        self.printAttrData()

    def printDfrInp(self):
        if self.dfrIn is None:
            print('Input DataFrame does not have any content yet.')
        else:
            print('Input DataFrame:')
            print(self.dfrIn)

    def printDDfrInp(self):
        if self.dDfrIn is None:
            print('Input DataFrames dictionary does not have any content yet.')
        else:
            print('Input DataFrames dictionary:')
            print(self.dDfrIn)

    def printDictsSortFilt(self):
        print('_'*8, 'Sorting dictionary:', '_'*8)
        for cK, cV in self.dSort.items():
            print('* ' + str(cK) + ':\n' + str(cV))
        print('_'*8, 'Threshold dictionary:', '_'*8)
        for cK, cV in self.dT.items():
            print('* ' + str(cK) + ':\n' + str(cV))
        print('_'*8, 'Selection dictionary:', '_'*8)
        for cK, cV in self.dSl.items():
            print('* ' + str(cK) + ':\n' + str(cV))

class FNmCmp(RootClass):
    def __init__(self, InpD, dSetsEq):
        super().__init__(InpD)
        self.idO = InpD.sFNmCmp
        self.descO = 'File name component'
        self.dSetsEq = dSetsEq
        print('Initiated "FNmCmp" base object.')

class FNmCmpICSglGT(FNmCmp):
    def __init__(self, InpD, dSetsEq, sK):
        super().__init__(InpD, dSetsEq)
        self.idO = InpD.sFNmCmpICSglGT
        self.descO = 'File name component for IC valid for all GTs'
        self.buildCmpIC(sK=sK)
        print('Initiated "FNmCmpICSglGT" base object.')

    def buildCmpSrt(self, sK=S_IC):
        if self.inpD.sortInFNm:
            if not self.dSetsEq[S_SRT]:
                for sGT in self.dCmp:
                    sKSrt = S_USC.join([self.dSort[sK][sGT][S_ORD],
                                        self.dSort[sK][sGT][S_SRT_BY]])
                    self.dCmp[sGT] = addIt(self.dCmp[sGT], sKSrt)

    def buildCmpThr(self, sK=S_IC):
        if not self.dSetsEq[S_THR]:
            for sGT in self.dCmp:
                for sMM in L_S_MIN_MAX:
                    self.dCmp[sGT] = addIt(self.dCmp[sGT],
                                           num2StrF(self.dT[sK][sGT][sMM]))

    def buildCmpSel(self):
        for tSl in self.dSl:
            if not self.dSetsEq[S_SEL][tSl]:
                dDoSub = {sGT: True for sGT in self.dCmp}
                for sKSub in self.dSl[tSl]:
                    for sGT in self.dCmp:
                        if dDoSub[sGT]:
                            setSl, sID = set(self.dSl[tSl][sKSub][sGT]), tSl[1]
                            if checkMPSetsEqSl(self.dSl[tSl], sGT):
                                dDoSub[sGT] = False
                            else:
                                sID = sKSub
                            self.dCmp[sGT] = addIt(self.dCmp[sGT],
                                                   getSel(setSl, sID))

    def assembleSCmp(self):
        for sGT in self.dCmp:
            if len(self.dCmp[sGT]) > 0:
                self.dCmp[sGT] = sGT + S_USC + self.dCmp[sGT]
                self.sCmp = addIt(self.sCmp, self.dCmp[sGT])

    def buildCmpIC(self, sK=S_IC):
        self.sCmp = S_USC.join([sK, self.dUsedK[sK]])
        self.dCmp = {sGT: '' for sGT in L_S_GT}
        self.buildCmpSrt(sK=sK)
        self.buildCmpThr(sK=sK)
        self.buildCmpSel()
        self.assembleSCmp()

class FNmCmpICAllGT(FNmCmp):
    def __init__(self, InpD, dSetsEq, sK):
        super().__init__(InpD, dSetsEq)
        self.idO = InpD.sFNmCmpICAllGT
        self.descO = 'File name component for IC valid for all GTs'
        self.buildCmpIC(sK=sK)
        print('Initiated "FNmCmpICAllGT" base object.')

    def buildCmpSrt(self, sK=S_IC, sGT=S_GT0):
        if self.inpD.sortInFNm:
            if self.dSetsEq[S_SRT]:
                sKSrt = S_USC.join([self.dSort[sK][sGT][S_ORD],
                                    self.dSort[sK][sGT][S_SRT_BY]])
                self.sCmp = addIt(self.sCmp, sKSrt)

    def buildCmpThr(self, sK=S_IC, sGT=S_GT0):
        if self.dSetsEq[S_THR]:
            for sMM in L_S_MIN_MAX:
                self.sCmp = addIt(self.sCmp, num2StrF(self.dT[sK][sGT][sMM]))

    def buildCmpSel(self):
        for tSl in self.dSl:
            if self.dSetsEq[S_SEL][tSl]:
                self.sCmp = addIt(self.sCmp, addSelAll(self.dSl, tSl))
        if len(self.sCmp) > 0:
            self.sCmp = addIt(S_ALL, self.sCmp)

    def buildCmpIC(self, sK=S_IC):
        self.sCmp = ''
        self.buildCmpSrt(sK=sK)
        self.buildCmpThr(sK=sK)
        self.buildCmpSel()

class FNmCmpDGT(FNmCmp):
    def __init__(self, InpD, dSetsEq, lSK):
        super().__init__(InpD, dSetsEq)
        self.idO = InpD.sFNmCmpDGT
        self.descO = 'File name component for dGT'
        self.sCmp = ''
        for i, sK in enumerate(lSK):
            self.buildCmpDGT(sK, i < (len(lSK) - 1))
        print('Initiated "FNmCmpDGT" base object.')

    def buildCmpDGT(self, sK, addUSC):
        if self.inpD.sortInFNm:
            sKSrt = S_USC.join([self.dSort[sK][S_ORD],
                                self.dSort[sK][S_SRT_BY]])
            self.sCmp += S_USC.join([sK, self.dUsedK[sK], sKSrt])
        else:
            self.sCmp += S_USC.join([sK, self.dUsedK[sK]])
        for sMM in L_S_MIN_MAX:
            self.sCmp += S_USC + num2StrF(self.dT[sK][sMM])
        if addUSC:
            self.sCmp += S_USC

class FileNameConstructor(RootClass):
    def __init__(self, InpD, tpFNm='F'):
        super().__init__(InpD)
        self.idO = InpD.sFNmCnstr
        self.descO = 'File name constructor'
        self.constrFNm(tpFNm=tpFNm)
        print('Initiated "FileNameConstructor" base object.')

    def checkSetsEq(self):
        self.dSetsEq = {S_SRT: checkAllSetsEqSrtThr(self.dSort[S_IC]),
                        S_THR: checkAllSetsEqSrtThr(self.dT[S_IC]),
                        S_SEL: {}}
        for tSl, dSlSub in self.dSl.items():
            self.dSetsEq[S_SEL][tSl] = checkAllSetsEqSl(dSlSub)

    def joinCmpToFNm(self):
        self.sFNm = self.sFOut
        l = [self.oCmpICSglGT.sCmp, self.oCmpICAllGT.sCmp, self.oCmpDGT.sCmp]
        for s in l:
            if len(s) > 0:
                self.sFNm += S_USC + comprStr(s, self.inpD.dComprStr)

    def constrFNm(self, tpFNm='F'):
        self.sFOut = self.inpD.sFOutF
        if tpFNm == 'S':
            self.sFOut = self.inpD.sFOutS
        self.checkSetsEq()
        self.oCmpICSglGT = FNmCmpICSglGT(self.inpD, self.dSetsEq, sK=S_IC)
        self.oCmpICAllGT = FNmCmpICAllGT(self.inpD, self.dSetsEq, sK=S_IC)
        if len(self.oCmpICSglGT.sCmp) == 0:
            self.oCmpICAllGT.sCmp = S_USC.join([self.dUsedK[S_IC],
                                                self.oCmpICAllGT.sCmp])
        self.oCmpDGT = FNmCmpDGT(self.inpD, self.dSetsEq, lSK=L_S_D_GT)
        self.joinCmpToFNm()

class ExtractedInfo(RootClass):
    def __init__(self, InpD):
        super().__init__(InpD)
        self.idO = InpD.sExtrInfo
        self.descO = 'Extracted info'
        self.dMap = self.inpD.dMapCHdSel
        self.sSp = self.inpD.sSep
        self.getProcData()
        print('Initiated "ExtractedInfo" base object.')

    def getPResF(self):
        cnstrFNmS = FileNameConstructor(self.inpD, tpFNm='S')
        cnstrFNmF = FileNameConstructor(self.inpD, tpFNm='F')
        pOut = self.inpD.pOutCSV
        self.pFOutS = os.path.join(pOut, S_DOT.join([cnstrFNmS.sFNm, S_CSV]))
        self.pFOutF = os.path.join(pOut, S_DOT.join([cnstrFNmF.sFNm, S_CSV]))
        self.dSpcSel = {'S': {'dfr': None, 'pF': self.pFOutS},
                        'F': {'dfr': None, 'pF': self.pFOutF}}

    def getInf4Inp(self):
        dDatTp_IC = {sIn: str for sIn in [self.inpD.sBC_L, self.inpD.sBC2_L]}
        dDatTp_P = {sIn: str for sIn in [self.inpD.sBC2_L]}
        self.dDatTp = {S_IC: dDatTp_IC, S_D_GT_P: dDatTp_P}
        self.dPFIn = {S_IC: self.inpD.dPFInIC,
                      S_D_GT_M: self.inpD.pFInM,
                      S_D_GT_P: self.inpD.pFInP}

    def loadDDfrInp(self):
        # load input DataFrames
        dDfrIn_IC = {sGT: pd.read_csv(self.dPFIn[S_IC][sGT], sep=self.sSp,
                                      dtype=self.dDatTp[S_IC])
                     for sGT in L_S_GT}
        dfrIn_M = pd.read_csv(self.dPFIn[S_D_GT_M], sep=self.sSp)
        dfrIn_P = pd.read_csv(self.dPFIn[S_D_GT_P], sep=self.sSp,
                              dtype=self.dDatTp[S_D_GT_P])
        self.dDfrIn = {S_IC: dDfrIn_IC,
                       S_D_GT_M: dfrIn_M,
                       S_D_GT_P: dfrIn_P}

    def getDHdCol(self):
        # get dictionary of column headers
        dDfr_IC = self.dDfrIn[S_IC]
        dHdCol_IC = {sGT: [sC + S_USC + sGT for sC in
                           list(dDfr_IC[sGT].loc[:, S_RMNG_COL1_IC:].columns)]
                     for sGT in L_S_GT}
        self.dHdCol = {S_IC: dHdCol_IC,
                       S_D_GT_M: list(self.dDfrIn[S_D_GT_M].columns),
                       S_D_GT_P: list(self.dDfrIn[S_D_GT_P].columns)}

    def getDMapK(self):
        # get mapping of data dictionary keys to column headers
        self.dMapK = {}
        for i, sMP in enumerate(L_S_M_P):
            self.dMapK[sMP] = (L_S_D_GT[i], sMP)
            lHdCRed = [s for s in self.dHdCol[L_S_D_GT[i]] if s not in L_S_M_P]
            for sHdC in lHdCRed:
                self.dMapK[S_USC.join([sHdC, sMP[0]])] = (L_S_D_GT[i], sHdC)
        for sGT in L_S_GT:
            for s in self.dHdCol[S_IC][sGT]:
                self.dMapK[s] = (S_IC, s)

    def getProcData(self):
        self.getPResF()
        self.getInf4Inp()
        self.loadDDfrInp()
        self.getDHdCol()
        self.getDMapK()

    def filterDfr(self, sKey, sMin=S_MIN, sMax=S_MAX):
        thMin, thMax = self.dT[sKey][sMin], self.dT[sKey][sMax]
        return thrFilter(self.dDfrIn[sKey], self.dUsedK[sKey], thMin, thMax)

    def filterAndConc(self, sKey, sMin=S_MIN, sMax=S_MAX):
        # filter data
        self.dDfrFl = {sKey: {}}
        for sGT in L_S_GT:
            thMin, thMax = self.dT[sKey][sGT][sMin], self.dT[sKey][sGT][sMax]
            self.dDfrFl[sKey][sGT] = thrFilter(self.dDfrIn[sKey][sGT],
                                               self.dUsedK[sKey], thMin, thMax)
        # process data - concatenate IC DataFrames of the three GT
        return concDfr(self.dDfrFl[sKey])

    def sortAndFiltDfr(self):
        for sGT, cDfr in self.dDfrIn[S_IC].items():
            sortDfr(cDfr, self.dSort[S_IC][sGT], S_SRT_BY, S_ORD)
        sortDfr(self.dDfrIn[S_D_GT_M], self.dSort[S_D_GT_M], S_SRT_BY, S_ORD)
        sortDfr(self.dDfrIn[S_D_GT_P], self.dSort[S_D_GT_P], S_SRT_BY, S_ORD)
        self.dDfrFl = {S_IC: self.filterAndConc(S_IC),
                       S_D_GT_M: self.filterDfr(S_D_GT_M),
                       S_D_GT_P: self.filterDfr(S_D_GT_P)}
        self.N = self.dDfrFl[S_D_GT_M].shape[0]*self.dDfrFl[S_D_GT_P].shape[0]

    def iniDfrRes(self, specSel=None):
        lC = list(self.dMapK)
        if specSel == 'S':  # the "short" subset of data
            # select only a subset of all keys (columns of DataFrames): lC
            lC = [S_BIN_L_2, S_SELECTED] + L_S_IC_GT + L_S_D_GT
            lC = [S_USC.join([s, S_PHO[0]]) for s in [S_BIN_L_2, S_SELECTED]]
            lC = L_S_M_P + lC + L_S_IC_GT
            lC += [S_USC.join([S_D_GT, s[0]]) for s in L_S_M_P]
            self.dfrResS = pd.DataFrame(columns=lC)
            self.dSpcSel[specSel]['dfr'] = self.dfrResS
        else:
            # select all keys (columns of DataFrames) defined in self.dMapK
            self.dfrResF = pd.DataFrame(columns=lC)
            self.dSpcSel[specSel]['dfr'] = self.dfrResF
        return {cK: [] for cK in lC}

    def SUB_fillDDat(self, dDat, i, j, sIMP):
        for sKDDt, lDt in dDat.items():
            assert sKDDt in self.dMapK
            (sKDFl, sHdC) = self.dMapK[sKDDt]
            if sKDFl in L_S_D_GT:
                assert sKDFl in self.dDfrFl
                assert sHdC in self.dDfrFl[sKDFl].columns
                if sKDFl == S_D_GT_M:
                    dDat[sKDDt].append(self.dDfrFl[sKDFl].at[i, sHdC])
                else:
                    dDat[sKDDt].append(self.dDfrFl[sKDFl].at[j, sHdC])
            elif sKDFl in [S_IC]:
                appendToDDat(lDt, self.dDfrFl[S_IC], sIMP, sHdC)
            else:
                print('ERROR: Unknown key of filter dictionary:', sKDFl)
                assert False

    def fillDDat(self, dDat):
        n = 0
        for i in self.dDfrFl[S_D_GT_M].index:
            sMet = self.dDfrFl[S_D_GT_M].at[i, S_MET]
            for j in self.dDfrFl[S_D_GT_P].index:
                sPho = self.dDfrFl[S_D_GT_P].at[j, S_PHO]
                sIMP = S_USC.join([sMet, sPho])
                self.SUB_fillDDat(dDat, i, j, sIMP)
                n += 1
                if n%self.inpD.modDisp == 0:
                    print('Processed element', n, 'of', self.N, '.')

    def fillSaveDfrRes(self):
        for spcSel in self.inpD.lSpecSel:
            dDat = self.iniDfrRes(specSel=spcSel)
            self.fillDDat(dDat)
            print('Filled data dictionary for selection "' + spcSel + '".')
            dISel = self.dSpcSel[spcSel]
            dISel['dfr'] = saveDfrRes(dISel['dfr'], dDat, dISel['pF'],
                                      self.sSp, self.dSl, self.dMap)

    def extractionOfExtremes(self):
        self.sortAndFiltDfr()
        self.fillSaveDfrRes()

class Plotter(RootClass):
    def __init__(self, InpD):
        super().__init__(InpD)
        self.idO = InpD.sPltr
        self.descO = 'Class for plotting'
        self.sSp = self.inpD.sSep
        self.dPPltF = {}
        print('Initiated "Plotter" base object.')

    def printDPPltF(self):
        print('Dictionary of plot file paths:')
        for tK, pF in self.dPPltF.items():
            print(tK, ':', pF)
        print('-'*64)

    def loadDfrInp(self, iC=0):
        dDatTp = {sIn: str for sIn in L_S_SG_GT}
        # load input DataFrame
        if hasattr(self, 'pFIn'):
            self.dfrIn = pd.read_csv(self.pFIn, sep=self.sSp, index_col=iC,
                                     dtype=dDatTp)

    def loadDDfrInp(self):
        # load input DataFrames, and save them in dictionary
        if hasattr(self, 'dPFIn'):
            self.dDfrIn = {sGT: None for sGT in L_S_GT}
            for sGT in self.dDfrIn:
                self.dDfrIn[sGT] = pd.read_csv(self.dPFIn[sGT], sep=self.sSp)

class PatternPlotter(Plotter):
    def __init__(self, InpD):
        super().__init__(InpD)
        self.idO = self.inpD.sPatPltr
        self.descO = 'Pattern plotter'
        self.pFIn = self.inpD.pFInPaP
        self.pDOut = self.inpD.pOutPaP
        self.dPlt = self.inpD.plotOfPatterns
        self.getDPPltF()
        self.loadDfrInp()
        print('Initiated "PatternPlotter" base object and loaded input data.')

    def getDPPltF(self):
        self.dPPltF, sFPlt = {}, self.inpD.sFOutPaP
        for ((s1, tSGT), tMP) in self.dPlt['dPairsPaP'].items():
            for sGT in tSGT:
                sPltF = S_DOT.join([S_USC.join([sFPlt, s1, sGT]), S_PDF])
                self.dPPltF[(s1, sGT)] = (tMP, os.path.join(self.pDOut, sPltF))

    def plotPatterns(self):
        d, nChD = self.dfrIn, self.dPlt['nCharDsp']
        for ((s1, sGT), ((sM, sP), pPltF)) in self.dPPltF.items():
            print('Plotting pattern for "' + s1 + '" and', sGT, '...')
            cSer = d[(d[S_MET] == sM) & (d[S_PHO] == sP)].squeeze()
            lTDat = [(S_MET, sM), (S_PHO, sP)]
            # if not os.path.isfile(pPltF):
            cFig, cAx = plt.subplots()
            for tMP in lTDat:
                cPa = cSer.loc[D_HD_C_PA[tMP[0]][sGT]]
                cPa.index = L_S_FT
                cAx.plot(cPa, lw=self.dPlt['lWdPlt'], label=tMP[1][:nChD])
            decorateClosePlot(cFig, cAx, self.dPlt, pPltF, S_YLBL_PAT_PLT)

class ICCmpPlotter(Plotter):
    def __init__(self, InpD):
        super().__init__(InpD)
        self.idO = self.inpD.sICCPltr
        self.descO = 'Concordance index component plotter'
        self.dPFIn = self.inpD.dPFInICP
        self.pDOut = self.inpD.pOutICP
        self.dPlt = self.inpD.plotOfICCmp
        self.getDPPltF()
        self.loadDDfrInp()
        print('Initiated "ICCmpPlotter" base object and loaded input data.')

    def getDPPltF(self):
        self.dPPltF, sFPlt = {}, self.inpD.sFOutICP
        for ((s1, tSGT), tMP) in self.dPlt['dPairsICP'].items():
            for sGT in tSGT:
                sPltF = S_DOT.join([S_USC.join([sFPlt, s1, sGT]), S_PDF])
                self.dPPltF[(s1, sGT)] = (tMP, os.path.join(self.pDOut, sPltF))

    def plotICCmp(self):
        nChD, wdBar = self.dPlt['nCharDsp'], self.dPlt['wdthBar']
        for ((s1, sGT), ((sM, sP), pPltF)) in self.dPPltF.items():
            print('Plotting IC components for "' + s1 + '" and', sGT, '...')
            d = self.dDfrIn[sGT]
            cSer = d[(d[S_MET] == sM) & (d[S_PHO] == sP)].squeeze()
            lTDat = [(S_MET, sM), (S_PHO, sP), (S_ICC, S_ICC_LEG)]
            xLocG = np.arange(len(L_S_FT_CHG))
            # if not os.path.isfile(pPltF):
            cFig, cAx = plt.subplots()
            for k, tMP in enumerate(lTDat):
                cICC = cSer.loc[D_HD_C_ICC[tMP[0]]]
                cICC.index = L_S_FT_CHG
                xLoc = xLocG + (2*k + 1)/(2*len(lTDat))*self.dPlt['wdthGrp']
                cAx.bar(xLoc - 1/2, height=cICC, width=wdBar,
                        lw=self.dPlt['lWdPlt'], label=tMP[1][:nChD])
            cAx.plot([-1/2, len(L_S_FT_CHG) + 1/2], [0, 0],
                     lw=self.dPlt['lWdPlt'], color='black')
            cAx.set_xticks(xLocG)
            cAx.set_xticklabels(L_S_FT_CHG)
            for cXLbl in cAx.get_xticklabels():
              cXLbl.set_rotation(self.dPlt['degRotXLbl'])
            decorateClosePlot(cFig, cAx, self.dPlt, pPltF, S_YLBL_ICC_PLT)

# --- MAIN --------------------------------------------------------------------
startTime = time.time()
print('+'*25 + ' START', time.ctime(startTime), '+'*25)

inpDat = InputData(dInput)
if inpDat.doInfoExtr:
    cXtrInfo = ExtractedInfo(inpDat)
    cXtrInfo.printIDDesc()
    cXtrInfo.printDfrInp()
    cXtrInfo.printDDfrInp()
    cXtrInfo.printAttrData()
    cXtrInfo.printDictsSortFilt()
    # cXtrInfo.printObjInfo()
    cXtrInfo.extractionOfExtremes()
if inpDat.doPlotPat:
    cPltr = PatternPlotter(inpDat)
    # cPltr.printAttrData()
    # cPltr.printDPPltF()
    cPltr.plotPatterns()
if inpDat.doPlotICC:
    cPltr = ICCmpPlotter(inpDat)
    cPltr.plotICCmp()

print('-'*80)
printElapsedTimeSim(startTime, time.time(), 'Total time')
print('+'*37, 'DONE.', '+'*36)

###############################################################################

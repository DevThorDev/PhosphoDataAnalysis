# -*- coding: utf-8 -*-
###############################################################################
# --- CalcStats.py ------------------------------------------------------------
# Calculate some simple statistics, and the p-value of the Kruskal-Wallis test
# for all metabolites and phosphopeptides, and all genotypes
###############################################################################
import os, time, itertools

from scipy import stats
import numpy as np
import pandas as pd

# --- CONSTANTS ---------------------------------------------------------------
S_USC = '_'
S_DOT = '.'
S_CSV = 'csv'

S_ALL_TP = 'AllTp'
S_TPM = 'Met'
S_TPP = 'Pho'
S_ALL_GT = 'AllGT'
S_GT0 = 'GT0'
S_GT1 = 'GT1'
S_GT5 = 'GT5'
S_ALL_FT = 'AllFt'
S_FT1 = 'Ft1'
S_FT2 = 'Ft2'
S_FT3 = 'Ft3'
S_FT4 = 'Ft4'

S_MN_ALL = 'mean conc. overall'
S_MN_FT1 = 'mean conc. DR'
S_MN_FT2 = 'mean conc. DS'
S_MN_FT3 = 'mean conc. NR'
S_MN_FT4 = 'mean conc. NS'
D_S_MN_FT = {S_ALL_FT: S_MN_ALL, S_FT1: S_MN_FT1, S_FT2: S_MN_FT2,
             S_FT3: S_MN_FT3, S_FT4: S_MN_FT4}
S_N_OBS_ALL = 'num. obs. total'
S_N_OBS_FT1 = 'num. obs. DR'
S_N_OBS_FT2 = 'num. obs. DS'
S_N_OBS_FT3 = 'num. obs. NR'
S_N_OBS_FT4 = 'num. obs. NS'
D_S_N_OBS_FT = {S_ALL_FT: S_N_OBS_ALL, S_FT1: S_N_OBS_FT1, S_FT2: S_N_OBS_FT2,
                S_FT3: S_N_OBS_FT3, S_FT4: S_N_OBS_FT4}
S_SD_ALL = 'SD overall'
S_MAD = 'MAD'
S_SD_FT1 = 'SD DR'
S_SD_FT2 = 'SD DS'
S_SD_FT3 = 'SD NR'
S_SD_FT4 = 'SD NS'
D_S_SD_FT = {S_ALL_FT: S_SD_ALL, S_FT1: S_SD_FT1, S_FT2: S_SD_FT2,
             S_FT3: S_SD_FT3, S_FT4: S_SD_FT4}
S_SD_MNS = 'SD means'
S_SD_RAT = 'SD ratio'
S_V_BTW = 'Vbtw'
S_V_WTH = 'Vwithin'
S_V_RAT = 'Vbtw/Vwithin'
S_KW = 'Kruskal-W.'

L_S_C_RES = [S_N_OBS_ALL, S_N_OBS_FT1, S_N_OBS_FT2, S_N_OBS_FT3, S_N_OBS_FT4,
             S_MAD, S_SD_ALL, S_SD_FT1, S_SD_FT2, S_SD_FT3, S_SD_FT4, S_SD_MNS,
             S_SD_RAT, S_V_BTW, S_V_WTH, S_V_RAT, S_KW]
L_S_C_FULL = [S_MN_ALL, S_MN_FT1, S_MN_FT2, S_MN_FT3, S_MN_FT4] + L_S_C_RES

R04 = 4

# --- INPUT -------------------------------------------------------------------
cMode = 'Full'              # 'Test' / 'Full'
lSCol = L_S_C_RES

pCSV = os.path.join('..', '..', '..', '25_Papers', '01_FirstAuthor',
                    '04_SysBio_DataAnalysis', '80_ResultsCSV')
sSep = ';'
sRes = 'Res'
sFEnd = 'SortConc'
sFStMet = 'ST1A__Metabolite_log2'
sFStPho = 'ST1B__Phosphopeptide_log2'

dSGT = {S_GT0: 'WT', S_GT1: 'PGM', S_GT5: 'SWEET'}
dSFt = {S_FT1: 'DR', S_FT2: 'DS', S_FT3: 'NR', S_FT4: 'NS'}
dKFt = {(S_TPM, S_GT0, S_FT1): list(range(1, 7)),
        (S_TPM, S_GT0, S_FT2): list(range(1, 7)),
        (S_TPM, S_GT0, S_FT3): list(range(1, 7)),
        (S_TPM, S_GT0, S_FT4): list(range(1, 7)),
        (S_TPM, S_GT1, S_FT1): list(range(1, 7)),
        (S_TPM, S_GT1, S_FT2): list(range(1, 7)),
        (S_TPM, S_GT1, S_FT3): list(range(1, 6)),
        (S_TPM, S_GT1, S_FT4): list(range(1, 6)),
        (S_TPM, S_GT5, S_FT1): list(range(1, 7)),
        (S_TPM, S_GT5, S_FT2): list(range(1, 7)),
        (S_TPM, S_GT5, S_FT3): list(range(1, 7)),
        (S_TPM, S_GT5, S_FT4): list(range(1, 7)),
        (S_TPP, S_GT0, S_FT1): list(range(1, 7)),
        (S_TPP, S_GT0, S_FT2): list(range(1, 7)),
        (S_TPP, S_GT0, S_FT3): list(range(1, 7)),
        (S_TPP, S_GT0, S_FT4): list(range(1, 7)),
        (S_TPP, S_GT1, S_FT1): list(range(1, 7)),
        (S_TPP, S_GT1, S_FT2): list(range(1, 7)),
        (S_TPP, S_GT1, S_FT3): [1, 2, 3, 5],
        (S_TPP, S_GT1, S_FT4): list(range(1, 6)),
        (S_TPP, S_GT5, S_FT1): list(range(1, 7)),
        (S_TPP, S_GT5, S_FT2): list(range(1, 7)),
        (S_TPP, S_GT5, S_FT3): list(range(1, 7)),
        (S_TPP, S_GT5, S_FT4): list(range(1, 7))}
dSCFt = {(sTp, sGT, sFt): [dSFt[sFt] + S_USC + dSGT[sGT] + S_USC + str(k)
                           for k in dKFt[(sTp, sGT, sFt)]]
         for (sTp, sGT, sFt) in dKFt}

# --- DERIVED INPUT -----------------------------------------------------------
dSTp = {S_TPM: sFStMet, S_TPP: sFStPho}
dSF = {(sTp, sGT): (dSTp[sTp] + S_USC + sGT + S_USC + sFEnd + S_DOT + S_CSV)
        for sTp in dSTp for sGT in dSGT}
if cMode == 'Test':
    dSF = {(sTp, sGT): (dSTp[sTp] + S_USC + sGT + S_USC + sFEnd + S_DOT + S_CSV)
            for sTp in [S_TPM] for sGT in [S_GT0]}

# --- FUNCTIONS ---------------------------------------------------------------
def flattenIt(cIterable, retArr = False):
    itFlat = list(itertools.chain.from_iterable(cIterable))
    if retArr:
        itFlat = np.array(itFlat)
    return itFlat

def getDVal1L(pdDfr, dSCFt, dSFt, sTp, sGT, i):
    dV1L = {sFt: pdDfr.loc[i, dSCFt[(sTp, sGT, sFt)]].infer_objects().to_list()
            for sFt in dSFt}
    return dV1L, flattenIt(dV1L.values())

def calcSimpleStats(x, nMinCent=1, nMinVar=2, nanP='omit'):
    nObs, cMean, cVar, cSD = 0, np.nan, np.nan, np.nan
    try:
        t = stats.describe(x, axis=None, nan_policy=nanP)
        nObs = t[0]
        if nObs >= nMinCent:
            cMean = t[2]
        if nObs >= nMinVar:
            cVar, cSD = t[3], np.sqrt(t[3])
    except:
        nObs = np.count_nonzero(~np.isnan(x))
    return nObs, cMean, cVar, cSD

def calcMAD(x, scl=1., nanP='omit'):
    return stats.median_abs_deviation(x, scale=scl, nan_policy=nanP)

def fillDictR(d, i, nObs, cMn, cSD, sF=None, lV=None):
    sK1, sK2, sK3, sK4 = S_N_OBS_ALL, S_MN_ALL, S_SD_ALL, S_MAD
    if sF is not None:
        sK1, sK2, sK3 = D_S_N_OBS_FT[sF], D_S_MN_FT[sF], D_S_SD_FT[sF]
    if nObs >= 1:
        d[i][sK1] = nObs
        d[i][sK2] = cMn
        d[i][sK3] = cSD
        if lV is not None:
            d[i][sK4] = calcMAD(lV)

def calcSDMnsRat(xMns, xSDs):
    sdMns, sdRat = np.nan, np.nan
    if len(xMns) >= 2:
        sdMns = stats.tstd(xMns)
        if len(xSDs) >= 1:
            mnSDs = stats.tmean(xSDs)
            if mnSDs > 0:
                sdRat = sdMns/mnSDs
    return sdMns, sdRat

def calcSDVals(d, i, dSFt):
    xMeans = [d[i][D_S_MN_FT[sFt]] for sFt in dSFt]
    xSDevs = [d[i][D_S_SD_FT[sFt]] for sFt in dSFt]
    d[i][S_SD_MNS], d[i][S_SD_RAT] = calcSDMnsRat(xMeans, xSDevs)

def calcVarVals(d, i, dV1L, dSFt):
    x, y = 0., 0.
    for sFt in dSFt:
        x += d[i][D_S_N_OBS_FT[sFt]]*(d[i][D_S_MN_FT[sFt]] - d[i][S_MN_ALL])**2
        for cV in dV1L[sFt]:
            if cV != '':
                if not np.isnan(cV):
                    y += (cV - d[i][D_S_MN_FT[sFt]])**2
    x /= (len(dSFt) - 1)
    y /= (d[i][S_N_OBS_ALL] - len(dSFt))
    d[i][S_V_BTW], d[i][S_V_WTH], d[i][S_V_RAT] = x, y, x/y

def writeBackRes(pdDfr, d, i, lSC=L_S_C_RES):
    lVClc = []
    for s in lSC:
        lVClc.append(d[i][s])
    pdDfr.loc[i, lSC] = lVClc

def saveAsCSV(pdDfr, pFC, sRs, sSp, wrtNmF=True):
    pFN = pFC[:-len(S_DOT + S_CSV)] + S_USC + sRs + S_DOT + S_CSV
    pdDfr.to_csv(pFN, sep=sSp)
    if wrtNmF:
        print('Saved results to', pFN, '...')

def calcValsCF(pDat, sF, sRs, dSCFt, dSFt, lSCol, sTp, sGT, sSp, nanP='omit'):
    print('_'*8, 'Starting with', (sTp, sGT), '_'*8)
    pF = os.path.join(pDat, sF)
    if os.path.exists(pF):
        cDfr, dClc = pd.read_csv(pF, sep=sSp), {}
        for iL in cDfr.index:
            dClc[iL] = {s: np.nan for s in lSCol}
            dV1L, lVInp = getDVal1L(cDfr, dSCFt, dSFt, sTp, sGT, iL)
            # calculate simple stats for all features together
            nObsA, cMeanA, _, cSDA = calcSimpleStats(lVInp)
            fillDictR(dClc, iL, nObsA, cMeanA, cSDA, lV=lVInp)
            # calculate simple stats for single features
            for sFt in dSFt:
                nObsS, cMeanS, _, cSDS = calcSimpleStats(dV1L[sFt])
                fillDictR(dClc, iL, nObsS, cMeanS, cSDS, sF=sFt)
            # calculate SD of mean and SD ratio
            calcSDVals(dClc, iL, dSFt)
            # calculate variance ratio
            calcVarVals(dClc, iL, dV1L, dSFt)
            # calculate p-value of Kruskal-Wallis test
            dClc[iL][S_KW] = stats.kruskal(*dV1L.values(), nan_policy=nanP)[1]
            # write back calculated values to DataFrame
            writeBackRes(cDfr, dClc, iL, lSC=lSCol)
        saveAsCSV(cDfr, pF, sRs=sRs, sSp=sSp)
    else:
        print('Path', pF, 'corresponding to', (sTp, sGT), 'does not exist!')

def printElapsedTimeSim(stT, cT, sPre = 'Time'):
    # calculate and display elapsed time
    elT = round(cT - stT, R04)
    print(sPre, 'elapsed:', elT, 'seconds, this is', round(elT/60, R04),
          'minutes or', round(elT/3600, R04), 'hours or',
          round(elT/(3600*24), R04), 'days.')

# --- MAIN --------------------------------------------------------------------
startTime = time.time()
print('+'*25 + ' START', time.ctime(startTime), '+'*25)
for ((sTp, sGT), sF) in dSF.items():
    calcValsCF(pCSV, sF, sRes, dSCFt, dSFt, lSCol, sTp, sGT, sSp=sSep)

print('-'*80)
printElapsedTimeSim(startTime, time.time(), 'Total time')
print('+'*37, 'DONE.', '+'*36)

###############################################################################

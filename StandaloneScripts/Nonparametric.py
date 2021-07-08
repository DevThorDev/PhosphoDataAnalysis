# -*- coding: utf-8 -*oddsRatio, pValue -
###############################################################################
# --- Nonparametric.py --------------------------------------------------------
###############################################################################
import os, time

import numpy as np
import pandas as pd
import scipy.stats as stats


# --- CONSTANTS ---------------------------------------------------------------
NM_GT0 = 'GT0'
NM_GT1 = 'GT1'
NM_GT2 = 'GT2'
NM_GT3 = 'GT3'
NM_GT4 = 'GT4'
NM_GT5 = 'GT5'

NM_EXT_CSV = 'csv'
NM_CT = 'Ct'

S_NEG = 'Neg'
S_POS = 'Pos'
S_SUM = 'Sum'
S_IDX = 'Idx'
S_CORR_L = 'Correlation'
S_DV_SC = 'DvSc'
S_DV_CL_SD = 'DvClSD'
S_SNDS = S_SUM + S_NEG + S_DV_SC
S_SPDS = S_SUM + S_POS + S_DV_SC
S_INDC = S_IDX + S_NEG + S_DV_CL_SD
S_IPDC = S_IDX + S_POS + S_DV_CL_SD

S_OVER_REP = 'OverRep'
S_N_OCC_ABS = 'nOccurrAbs'
S_P_VAL_OV = 'pValueFOver'
S_P_VAL_UN = 'pValueFUnder'
S_P_OF = 'p'
S_M_CORR_BON = 'Bonferroni'

R04 = 4

S_MD = 'MetD'
S_PD = 'PhoD'

# --- TEMP INPUT --------------------------------------------------------------
# lSCol = ['Atlantic', 'Indian']
# lSRow = ['Whales', 'Sharks']
# arrDatF = [[100 - 9, 9], [3330*58 - 100 - 58 + 9, 58 - 9]]
# arrDatF = [[1, 9], [11, 3]]
# arrDatF = [[0, 10], [12, 2]]
# arrDatF = [[8., 2], [1, 5.]]
# arrDatF = [[8000., 2000], [1000, 5000.]]
# lAltFisher = ['less', 'greater']    # 'two-sided', 'less', 'greater'

# --- INPUT -------------------------------------------------------------------
cAttr = S_MD         # S_MD / S_PD
cGTAttr1 = NM_GT0       # NM_GT0 / NM_GT1 / NM_GT2 / NM_GT3 / NM_GT4 / NM_GT5
cGTAttr2 = NM_GT0       # NM_GT0 / NM_GT1 / NM_GT2 / NM_GT3 / NM_GT4 / NM_GT5
sNmDatF = 'Corr__BinOp_MetD_' + cGTAttr1 + '_AllD_PhoD_' + cGTAttr2 + '_AllD'
nMin, nMax = 1, 100
lSrt = [S_IPDC, S_SPDS, S_CORR_L]
lAsc = [False, False, False]
sMCorr = S_M_CORR_BON    # None / S_M_CORR_BON

pRelDatF = os.path.join('..', '..', '..', '12_SysBio02_DataAnalysis',
                        '21_R_81_BinaryOps')
pRelResF = os.path.join('..', '..', '..', 'TEMP', '12_SysBio02_DataAnalysis',
                        '80_ResultsCSV')

cSep = ';'

# --- INPUT PROCESSING --------------------------------------------------------
sIDCAttr, sIDOAttr = S_MD, S_PD
if cAttr == S_MD:
    sIDCAttr, sIDOAttr = S_MD, S_PD
elif cAttr == S_PD:
    sIDCAttr, sIDOAttr = S_PD, S_MD

sNmResF = S_OVER_REP + cAttr + '_' + sNmDatF
sNmDatF += ('.' + NM_EXT_CSV)
sNmResF += ('.' + NM_EXT_CSV)

# --- FUNCTIONS ---------------------------------------------------------------
def createDir(pF):
    if not os.path.isdir(pF):
        os.mkdir(pF)

def joinToPath(pF = '', nmF = 'Dummy.txt'):
    if len(pF) > 0:
        createDir(pF)
        return os.path.join(pF, nmF)
    else:
        return nmF

def readCSV(pF, sepD = ',', iCol = None):
    return pd.read_csv(pF, sep = sepD, index_col = iCol)

def getSrtData(pF, sNmF, sIDAttr1, sIDAttr2, lSrt, lAsc, sepD):
    pdDfr = readCSV(joinToPath(pF, sNmF), sepD = sepD)
    lAttr1 = list(set(pdDfr.loc[:, sIDAttr1]))
    lAttr2 = list(set(pdDfr.loc[:, sIDAttr2]))
    Ntot = len(lAttr1)*len(lAttr2)
    pdDfr = pdDfr.sort_values(by = lSrt, ascending = lAsc)
    return pdDfr, lAttr1, lAttr2, len(lAttr1), len(lAttr2), Ntot

def getNmFOut(sNmFBase, cSSp, sXt, sCorr = None):
    cSNm = sNmFBase[:-(len(sXt) + 1)] + '__' + cSSp
    if sCorr is not None:
        cSNm += ('_' + sCorr)
    cSNm += sNmFBase[-(len(sXt) + 1):]
    return cSNm

def iniKeyV1D(lD, cKey, lenArr):
    for cD in lD:
        cD[cKey] = np.zeros((lenArr))

def getCntTblF(dOccAbs, nCur, nMin, nOAttr, nTot, cArr, cAttr):
    nOccAbs = len(cArr[cArr == cAttr])
    dOccAbs[cAttr][nCur - nMin] = nOccAbs
    a00, a01 = nCur - nOccAbs, nOccAbs
    a10, a11 = nTot - nCur - nOAttr + nOccAbs, nOAttr - nOccAbs
    if min([a00, a01, a10, a11]) < 0:
        print('[a00, a01, a10, a11] =', [a00, a01, a10, a11])
        assert False
    return np.array([[a00, a01], [a10, a11]])

def exactFisher(arrDat, cAlt = 'two-sided'):
    arrDat = np.array(arrDat, dtype = np.uint64)
    assert arrDat.shape == (2, 2)
    # returns (oddsRatio, pVal)
    return stats.fisher_exact(arrDat, alternative = cAlt)

def calcPValsF(dPValOv, dPValUn, arrCTblF, nCur, nMin, nCAttr,
               cAttr, sCorr = None):
    pValOv = exactFisher(arrCTblF, cAlt = 'less')[1]
    pValUn = exactFisher(arrCTblF, cAlt = 'greater')[1]
    if sCorr == S_M_CORR_BON:   # Bonferroni correction
        pValOv *= nCAttr
        pValUn *= nCAttr
    dPValOv[cAttr][nCur - nMin] = pValOv
    dPValUn[cAttr][nCur - nMin] = pValUn

def calcPValProfile(dOccAbs, dPValOv, dPValUn, dfrRd, lCAttr, nMin, nMax,
                    nTot, nCAttr, nOAttr, sIDCAttr, sCorr = None):
    for k, cAttr in enumerate(lCAttr):
        iniKeyV1D([dOccAbs, dPValOv, dPValUn], cAttr, nMax - nMin + 1)
        for nCur in range(nMin, nMax + 1):
            arrS = dfrRd.iloc[:nCur, :].loc[:, sIDCAttr].values
            arrCntTblF = getCntTblF(dOccAbs, nCur, nMin, nOAttr, nTot, arrS,
                                    cAttr)
            calcPValsF(dPValOv, dPValUn, arrCntTblF, nCur, nMin, nCAttr,
                       cAttr, sCorr = sCorr)
        print('Done:', k + 1, 'of', len(lCAttr))

def insStrAtFNmEnd(sNmF, sAdd, XSep = '.'):
    lSpl = sNmF.split(XSep)
    assert len(lSpl) >= 1
    return (XSep.join([lSpl[k] for k in range(len(lSpl) - 1)]) + sAdd + XSep +
            lSpl[-1])

def saveDictAsPdDfr(pF, sNmF, cD, nMin, nMax, cSep, sP = None):
    cDfr = pd.DataFrame(cD, index = range(nMin, nMax + 1))
    cDfr.to_csv(joinToPath(pF, sNmF), sep = cSep)
    if sP is not None:
        cDfr = -np.log10(cDfr)
        cDfr.to_csv(joinToPath(pF, insStrAtFNmEnd(sNmF, sP)), sep = cSep)

def exactFisherPd(pdDfr, cAlt = 'two-sided'):
    # returns (oddsRatio, pVal)
    return exactFisher(pdDfr.values, cAlt = cAlt)

def addToCountDict(cD, cK, cInc = 1):
    if cK in cD:
        cD[cK] += cInc
    else:
        cD[cK] = cInc

def chi2Cont(arrDat, doCorr = True):
    arrDat = np.array(arrDat, dtype = np.uint64)
    # returns (chi2, pVal, dof, expFreq)
    return stats.chi2_contingency(arrDat, correction = doCorr)

def chi2ContPd(pdDfr, doCorr = True):
    # returns (chi2, pVal, dof, expFreq)
    return chi2Cont(pdDfr.values, doCorr = doCorr)

def printElapsedTimeSim(stT, cT, sPre = 'Time'):
    # calculate and display elapsed time 
    elT = round(cT - stT, R04)
    print(sPre, 'elapsed:', elT, 'seconds, this is', round(elT/60, R04),
          'minutes or', round(elT/3600, R04), 'hours or',
          round(elT/(3600*24), R04), 'days.')

# --- MAIN --------------------------------------------------------------------
# dfrT1 = pd.DataFrame(arrDatF, index = lSRow, columns = lSCol)
# oddsRatioF, pValueF = exactFisherPd(dfrT1, cAlt = lAltFisher[0])
# chiSq1, pValueC1, doF1, expFr1 = chi2ContPd(dfrT1)
# chiSq2, pValueC2, doF2, expFr2 = chi2ContPd(dfrT1, doCorr = False)
# print('oddsRatios:', oddsRatioF, chiSq1, chiSq2)
# print('pValues:', pValueF, pValueC1, pValueC2)
# print('doF:', '-', doF1, doF2)
# print('expFr:', '-', expFr1, expFr2)


startTime = time.time()
print('+'*50 + ' START', time.ctime(startTime), '+'*30)
print('Nonparametric')

t = getSrtData(pRelDatF, sNmDatF, sIDCAttr, sIDOAttr, lSrt, lAsc, cSep)
dfrRd, lCAttr, lOAttr, nCAttr, nOAttr, N = t
dOccAbs, dPValOv, dPValUn, lSNm = {}, {}, {}, []
sNmResFOccAbs = getNmFOut(sNmResF, S_N_OCC_ABS, NM_EXT_CSV)
sNmResFPValOv = getNmFOut(sNmResF, S_P_VAL_OV, NM_EXT_CSV, sCorr = sMCorr)
sNmResFPValUn = getNmFOut(sNmResF, S_P_VAL_UN, NM_EXT_CSV, sCorr = sMCorr)

calcPValProfile(dOccAbs, dPValOv, dPValUn, dfrRd, lCAttr, nMin, nMax, N,
                nCAttr, nOAttr, sIDCAttr, sCorr = sMCorr)
saveDictAsPdDfr(pRelResF, sNmResFOccAbs, dOccAbs, nMin, nMax, cSep)
saveDictAsPdDfr(pRelResF, sNmResFPValOv, dPValOv, nMin, nMax, cSep,
                sP = '_' + S_P_OF)
saveDictAsPdDfr(pRelResF, sNmResFPValUn, dPValUn, nMin, nMax, cSep,
                sP = '_' + S_P_OF)
print('-'*80)
printElapsedTimeSim(startTime, time.time(), 'Total time')
print('*'*20 + ' DONE', time.ctime(time.time()), '*'*20)

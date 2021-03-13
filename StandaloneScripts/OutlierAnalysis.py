# -*- coding: utf-8 -*-
###############################################################################
# --- OutlierAnalysis.py ------------------------------------------------------
###############################################################################
import os

import numpy as np
import pandas as pd
# import matplotlib.pyplot as plt

# --- CONSTANTS ---------------------------------------------------------------
dGT = {'GT0': 'WT', 'GT1': 'PGM', 'GT5': 'SWEET'}
lFt = ['DR', 'DS', 'NR', 'NS']
cSep = ';'
sCSV = 'csv'
sTXT = 'txt'
sP_CSV = os.path.join('..', '..', '..', '12_SysBio02_DataAnalysis',
                      '02_ProcessedRawData', '01_CSV')
S_2L = 'log2'
S_UN = 'unTransformed'

# --- INPUT -------------------------------------------------------------------
cTr = S_2L
lThrDef = [3., 3.5, 4., 4.5, 5., 6., 8.]
nDigRnd = 4
sFResultFull = 'FullOutlierResult'
sFNumOcc = 'NumOccurrences'

# --- DERIVED VALUES ----------------------------------------------------------
dDAttr = {}
for cGT, cNmGT in dGT.items():
    dDAttr[cGT] = {cFt: [cFt + '_' + dGT[cGT] + '_' + str(k) for k in
                         range(1, 6 + 1)] for cFt in lFt}
    if cGT in ['GT1']:
        # special case PGM mutant: NR and NS only 5 measurements
        for cFt in ['NR', 'NS']:
            dDAttr[cGT][cFt] = [cFt + '_' + dGT[cGT] + '_' + str(k) for k in
                                range(1, 5 + 1)]
dThrDef = {k: cThr for k, cThr in enumerate(lThrDef)}

# --- GENERAL FUNCTIONS -------------------------------------------------------
def addToDictCt(cD, cK, cCt=1):
    if cK in cD:
        cD[cK] += cCt
    else:
        cD[cK] = cCt

def addToDictD(cD, cKO, cKI=None, cEI=None):
    if cKO in cD:
        if cKI is not None:
            cD[cKO][cKI] = cEI
    else:
        if cKI is None:
            cD[cKO] = {}
        else:
            cD[cKO] = {cKI: cEI}

# --- SPECIAL FUNCTIONS -------------------------------------------------------
def getOutliersQuantile(pdSer, qntLow=.05, qntUp=.95):
    return pdSer[pdSer.between(pdSer.quantile(qntLow), pdSer.quantile(qntUp))]

def is_outlier(points, thresh=3.5):
    """
    Returns a boolean array with True if points are outliers and False
    otherwise.

    Parameters:
    -----------
        points : An numobservations by numdimensions array of observations
        thresh : The modified z-score to use as a threshold. Observations with
            a modified z-score (based on the median absolute deviation) greater
            than this value will be classified as outliers.

    Returns:
    --------
        mask : A numobservations-length boolean array.

    References:
    ----------
        Boris Iglewicz and David Hoaglin (1993), "Volume 16: How to Detect and
        Handle Outliers", The ASQC Basic References in Quality Control:
        Statistical Techniques, Edward F. Mykytka, Ph.D., Editor.
    """
    if len(points.shape) == 1:
        points = points[:, None]
    median = np.median(points, axis=0)
    diff = np.sum((points - median)**2, axis=-1)
    diff = np.sqrt(diff)
    med_abs_deviation = np.median(diff)

    modified_z_score = 0.67449 * diff/med_abs_deviation

    return modified_z_score > thresh

def getOutliersMAD(pdSer, dThr=dThrDef):
    cMedian = pdSer.median()
    cDiff = np.sqrt((pdSer - cMedian)**2)
    medAbsDev = np.median(cDiff)
    zScoreMod = 0.67449 * cDiff/medAbsDev
    return {k: (cThr, pdSer[zScoreMod > cThr]) for k, cThr in dThr.items()}

def analyseDataFrame(pdDfr, dAttr):
    dResDat, dResThr, dNOcc = {}, {}, {}
    for sI in pdDfr.index:
        for sAttr, lAttr in dAttr.items():
            cSer = pdDfr.loc[sI, lAttr]
            dResDat[(sI, sAttr)] = getOutliersMAD(cSer)
    for tK, cDOutl in dResDat.items():
        for k, tV in cDOutl.items():
            dVOutl = tV[1].to_dict()
            if len(dVOutl) > 0:
                addToDictD(dResThr, k, tK, (tV[0], dVOutl))
    for i, dOutl in dResThr.items():
        for (cKO, (cThr, cDO)) in dOutl.items():
            addToDictD(dNOcc, i)
            addToDictCt(dNOcc[i], cKO[1])
            for cKI, cVI in cDO.items():
                addToDictCt(dNOcc[i], cKI)
    return dResDat, dResThr, dNOcc

def printOutliers(dResThr, cGT):
    print('='*16, cGT, '|', 'START', '='*16)
    for dOutl in dResThr.values():
        for (i, (tSet, (cThr, dOutl))) in enumerate(dOutl.items()):
            if i == 0:
                print('-'*16, 'Threshold', cThr, '-'*16)
            print(tSet, ':', dOutl)
    print('='*16, cGT, '|', 'END', '='*16)

def saveFullOutlierResult(dResThr, cGT, sFResSt, nDgRd, sFEnd=sCSV, cSp=cSep):
    sFRes = sFResSt + '_' + cGT + '.' + sFEnd
    fRes = open(sFRes, mode='w')
    for dOutl in dResThr.values():
        for (cKO, (cThr, cDO)) in dOutl.items():
            fRes.write(str(cThr) + cSp + cKO[0] + cSp + cKO[1])
            for cKI, cVI in cDO.items():
                fRes.write(cSp + cKI + cSp + str(round(cVI, nDgRd)))
            fRes.write('\n')
    fRes.close()

def getFirstLine(dAttr, lF=lFt, cSp=cSep):
    sP = 'Threshold'
    for sF in lF:
        sP += cSp + sF
    for sF in lF:
        for sM in dAttr[sF]:
            sP += cSp + sM
    return sP + '\n'

def addToLineStr(cD, cS, nDef=0, cSp=cSep):
    if cS in cD:
        sP = cSp + str(cD[cS])
    else:
        sP = cSp + str(nDef)
    return sP

def getDataLine(dNOccThr, dAttr, k=0, lF=lFt, dThr=dThrDef, cSp=cSep):
    sL = str(dThr[k])
    for sF in lF:
        sL += addToLineStr(dNOccThr, sF, cSp=cSp)
    for sF in lF:
        for sM in dAttr[sF]:
            sL += addToLineStr(dNOccThr, sM, cSp=cSp)
    return sL

def saveNumOcc(dNOcc, dAttr, cGT, sFResSt, lF=lFt, dThr=dThrDef, sFEnd=sCSV,
               cSp=cSep):
    sFRes = sFResSt + '_' + cGT + '.' + sFEnd
    fRes = open(sFRes, mode='w')
    fRes.write(getFirstLine(dAttr, lF=lF, cSp=cSp))
    for i, dNOccThr in dNOcc.items():
        sLD = getDataLine(dNOccThr, dAttr, i, lF=lF, dThr=dThr, cSp=cSp)
        fRes.write(sLD + '\n')
    fRes.close()

# --- MAIN --------------------------------------------------------------------
print('*'*16, 'START', '*'*16)
# for cGT in ['GT0']:
for cGT in dGT:
    sF_CSV = 'Metabolite_' + cTr + '_' + cGT + '_' + dGT[cGT] + '_ThS'
    pFDat = os.path.join(sP_CSV, sF_CSV + '.' + sCSV)
    cDfr = pd.read_csv(pFDat, sep=cSep, index_col=0)
    dResData, dResThres, dNumOcc = analyseDataFrame(cDfr, dDAttr[cGT])
    # printOutliers(dResThres, cGT)
    saveFullOutlierResult(dResThres, cGT, sFResultFull, nDigRnd)
    saveNumOcc(dNumOcc, dDAttr[cGT], cGT, sFNumOcc)
    print('Saved results for genotype', cGT + '.')
print('*'*16, 'DONE', '*'*16)

# -----------------------------------------------------------------------------

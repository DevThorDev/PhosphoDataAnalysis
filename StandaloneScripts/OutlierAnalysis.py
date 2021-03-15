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
lSPR_Base = ['..', '..', '..', '12_SysBio02_DataAnalysis']
lSPR_Dat = ['02_ProcessedRawData', '01_CSV']
lSPR_Res = ['04_OutlierResults']
S_2L = 'log2'
S_UN = 'unTransformed'

# --- INPUT -------------------------------------------------------------------
cTrans = S_2L           # S_2L / S_UN
lThrDef = [3., 3.5, 4., 4.5, 5., 6., 8.]
useIQR = False          # calculate MAD with IQR if True, otherwise use median

nDigRnd = 4
sFDatSt = 'Metabolite'
sFResultFull = 'FullOutlierResult'
sFResultDfr = 'OutlierThresholds'
sFNumOcc = 'NumOccOutlier'
lSColSkip = ['MeanConc']

# --- DERIVED VALUES ----------------------------------------------------------
sMode = 'median'
if useIQR:
    sMode = 'IQR'
lSPR_Dat, lSPR_Res = lSPR_Base + lSPR_Dat, lSPR_Base + lSPR_Res
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

# --- OS FUNCTIONS ------------------------------------------------------------
def makeDirs(pDTarget):
    if not os.path.isdir(pDTarget):
        os.makedirs(pDTarget)

def joinToPath(lD4P=[], sF=None):
    if len(lD4P) > 0:
        pD = ''
        for cD in lD4P:
            pD = os.path.join(pD, cD)
        makeDirs(pD)
        if sF is None:
            return pD
        else:
            return os.path.join(pD, sF)
    else:
        if sF is None:
            return '.'
        else:
            return sF

def getPFDat(sFDatSt, lD4P=[], dGT=dGT, sTr=cTrans, cGT='GT0', sFExt=sCSV):
    sFDat = sFDatSt + '_' + sTr + '_' + cGT + '_' + dGT[cGT] + '_ThS'
    return joinToPath(lD4P=lD4P, sF=sFDat + '.' + sFExt)

def getPFRes(sFResSt, lD4P=[], sTr=cTrans, sMd=sMode, sFExt=sCSV):
    sFRes = sFResSt + '_' + sTr + '_' + sMd + '_' + cGT
    return joinToPath(lD4P=lD4P, sF=sFRes + '.' + sFExt)

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

    modified_z_score = 0.67449*diff/med_abs_deviation

    return modified_z_score > thresh

def calcMadIQR(pdSer, qtLB=0.25, qtUB=0.75):
    npArr = pdSer.to_numpy()
    IQR = np.quantile(npArr, qtUB, axis=0) - np.quantile(npArr, qtLB, axis=0)
    # as the IQR may be zero, using the variance avoids inf
    if IQR == 0:
        IQR = np.var(npArr)
    return np.median(np.abs(npArr)/IQR)

def getOutliersMAD(pdSer, dThr=dThrDef, useIQR=False, zFact=0.67449):
    cDiff = np.sqrt((pdSer - pdSer.median())**2)
    if useIQR:
        medAbsDev = calcMadIQR(pdSer)
    else:
        medAbsDev = np.median(cDiff)
    zScoreMod = zFact*cDiff/medAbsDev
    return {k: (cThr, pdSer[zScoreMod > cThr]) for k, cThr in dThr.items()}

def analyseDataFrame(pdDfr, dAttr, useIQR=False):
    dResDat, dResThr, dNOcc = {}, {}, {}
    for sI in pdDfr.index:
        for sAttr, lAttr in dAttr.items():
            cSer = pdDfr.loc[sI, lAttr]
            dResDat[(sI, sAttr)] = getOutliersMAD(cSer, useIQR=useIQR)
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

def saveFullOutlierResult(sPRes, dResThr, cGT, nDgRd, cSp=cSep):
    fRes = open(sPRes, mode='w')
    for dOutl in dResThr.values():
        for (cKO, (cThr, cDO)) in dOutl.items():
            fRes.write(str(cThr) + cSp + cKO[0] + cSp + cKO[1])
            for cKI, cVI in cDO.items():
                fRes.write(cSp + cKI + cSp + str(round(cVI, nDgRd)))
            fRes.write('\n')
    fRes.close()

def saveDataFrameThr(sPRes, dResThr, inpDfr, cGT, lSCSkip=None, cSp=cSep):
    lCOut = [s for s in inpDfr.columns if s not in lSCSkip]
    outDfr = pd.DataFrame(0., index=inpDfr.index, columns=lCOut)
    for dOutl in dResThr.values():
        for (cKO, (cThr, cDO)) in dOutl.items():
            assert cKO[0] in outDfr.index
            for cKI in cDO:
                assert cKI in outDfr.columns
                outDfr.at[cKO[0], cKI] = cThr
    outDfr.to_csv(sPRes, sep=cSp, na_rep='NA')

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

def saveNumOcc(sPRes, dNOcc, dAttr, cGT, lF=lFt, dThr=dThrDef, cSp=cSep):
    fRes = open(sPRes, mode='w')
    fRes.write(getFirstLine(dAttr, lF=lF, cSp=cSp))
    for i, dNOccThr in dNOcc.items():
        sLD = getDataLine(dNOccThr, dAttr, i, lF=lF, dThr=dThr, cSp=cSp)
        fRes.write(sLD + '\n')
    fRes.close()

# --- MAIN --------------------------------------------------------------------
print('*'*16, 'START', cTrans, 'with mode', sMode, '*'*16)
# for cGT in ['GT0']:
for cGT in dGT:
    pFDat = getPFDat(sFDatSt, lSPR_Dat, dGT, sTr=cTrans, cGT=cGT, sFExt=sCSV)
    cDfr = pd.read_csv(pFDat, sep=cSep, index_col=0)
    dResData, dResThres, dNumOcc = analyseDataFrame(cDfr, dDAttr[cGT], useIQR)
    # printOutliers(dResThres, cGT)
    pFRes = getPFRes(sFResultFull, lSPR_Res, sTr=cTrans, sMd=sMode, sFExt=sCSV)
    saveFullOutlierResult(pFRes, dResThres, cGT, nDgRd=nDigRnd, cSp=cSep)
    pFRes = getPFRes(sFResultDfr, lSPR_Res, sTr=cTrans, sMd=sMode, sFExt=sCSV)
    saveDataFrameThr(pFRes, dResThres, cDfr, cGT, lSCSkip=lSColSkip, cSp=cSep)
    pFRes = getPFRes(sFNumOcc, lSPR_Res, sTr=cTrans, sMd=sMode, sFExt=sCSV)
    saveNumOcc(pFRes, dNumOcc, dDAttr[cGT], cGT, lFt, dThr=dThrDef, cSp=cSep)
    print('Saved results for genotype', cGT + '.')
print('*'*16, 'DONE', cTrans, 'with mode', sMode, '*'*16)

# -----------------------------------------------------------------------------

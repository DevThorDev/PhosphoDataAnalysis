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
sP_CSV_GT0 = os.path.join('..', '..', '..', '12_SysBio02_DataAnalysis',
                          '02_ProcessedRawData', '01_CSV')
S_2L = 'log2'
S_UN = 'unTransformed'

# --- INPUT -------------------------------------------------------------------
cGT = 'GT5'
cTr = S_2L

sF_CSV_GT0 = 'Metabolite_' + cTr + '_' + cGT + '_' + dGT[cGT] + '_ThS'

dLAttr = {cFt: [cFt + '_' + dGT[cGT] + '_' + str(k) for k in range(1, 6 + 1)]
          for cFt in lFt}

# --- FUNCTIONS ---------------------------------------------------------------
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

# def getOutliersQuantile(npArr, qntLow=.05, qntUp=.95):
#     pdSer = pd.Series(npArr)
#     return pdSer[pdSer.between(pdSer.quantile(qntLow), pdSer.quantile(qntUp))]

def getOutliersMAD(pdSer, lThr=[2., 2.5, 3., 3.5, 4., 4.5, 5.]):
# def getOutliersMAD(pdSer, lThr=[3., 5.]):
    cMedian = pdSer.median()
    cDiff = np.sqrt((pdSer - cMedian)**2)
    medAbsDev = np.median(cDiff)
    zScoreMod = 0.67449 * cDiff/medAbsDev
    return {k: (cThr, pdSer[zScoreMod > cThr]) for k, cThr in enumerate(lThr)}

def getOutliersQuantile(pdSer, qntLow=.05, qntUp=.95):
    return pdSer[pdSer.between(pdSer.quantile(qntLow), pdSer.quantile(qntUp))]

def addToDictD(cD, cKO, cKI, cEI):
    if cKO in cD:
        cD[cKO][cKI] = cEI
    else:
        cD[cKO] = {cKI: cEI}

def analyseDataFrame(pdDfr, dAttr):
    dResDat, dResThr = {}, {}
    for sI in pdDfr.index:
        for sAttr, lAttr in dAttr.items():
            cSer = pdDfr.loc[sI, lAttr]
            dResDat[(sI, sAttr)] = getOutliersMAD(cSer)
    for tK, cDOutl in dResDat.items():
        for k, tV in cDOutl.items():
            dVOutl = tV[1].to_dict()
            if len(dVOutl) > 0:
                addToDictD(dResThr, k, tK, (tV[0], dVOutl))
    return dResDat, dResThr

def printOutliers(dResThr):
    print('='*16, 'START', '='*16)
    for dOutl in dResThr.values():
        for (i, (tSet, (cThr, lOutl))) in enumerate(dOutl.items()):
            if i == 0:
                print('-'*16, 'Threshold', cThr, '-'*16)
            print(tSet, ':', lOutl)
    print('='*16, 'END', '='*16)

# --- MAIN --------------------------------------------------------------------

pFDat = os.path.join(sP_CSV_GT0, sF_CSV_GT0 + '.' + sCSV)
cDfr = pd.read_csv(pFDat, sep=cSep, index_col=0)
dResData, dResThres = analyseDataFrame(cDfr, dLAttr)
printOutliers(dResThres)

print('-'*16, 'DONE', '-'*16)

# -----------------------------------------------------------------------------

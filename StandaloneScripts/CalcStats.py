# -*- coding: utf-8 -*-
###############################################################################
# --- CalcStats.py ------------------------------------------------------------
###############################################################################
import os

from scipy import stats
import numpy as np
import pandas as pd

# --- CONSTANTS ---------------------------------------------------------------
S_USC = '_'
S_DOT = '.'
S_CSV = 'csv'

S_KW = 'Kruskal-W.'

# --- INPUT -------------------------------------------------------------------
pCSV = os.path.join('..', '..', '..', '25_Papers', '01_FirstAuthor',
                    '04_SysBio_DataAnalysis', '80_ResultsCSV')
sSep = ';'
sRes = 'Res'
sFEnd = 'SortConc'
sFStMet = 'ST1A__Metabolite_log2'
sFStPho = 'ST1B__Phosphopeptide_log2'
dSTp = {'Met': sFStMet, 'Pho': sFStPho}
dSGT = {'GT0': 'WT', 'GT1': 'PGM', 'GT5': 'SWEET'}
dSFt = {'Ft1': 'DR', 'Ft2': 'DS', 'Ft3': 'NR', 'Ft4': 'NS'}
dSF = {(sTp, sGT): (dSTp[sTp] + S_USC + sGT + S_USC + sFEnd + S_DOT + S_CSV)
       for sTp in dSTp for sGT in dSGT}
dKFt = {('Met', 'GT0', 'Ft1'): list(range(1, 7)),
        ('Met', 'GT0', 'Ft2'): list(range(1, 7)),
        ('Met', 'GT0', 'Ft3'): list(range(1, 7)),
        ('Met', 'GT0', 'Ft4'): list(range(1, 7)),
        ('Met', 'GT1', 'Ft1'): list(range(1, 7)),
        ('Met', 'GT1', 'Ft2'): list(range(1, 7)),
        ('Met', 'GT1', 'Ft3'): list(range(1, 6)),
        ('Met', 'GT1', 'Ft4'): list(range(1, 6)),
        ('Met', 'GT5', 'Ft1'): list(range(1, 7)),
        ('Met', 'GT5', 'Ft2'): list(range(1, 7)),
        ('Met', 'GT5', 'Ft3'): list(range(1, 7)),
        ('Met', 'GT5', 'Ft4'): list(range(1, 7)),
        ('Pho', 'GT0', 'Ft1'): list(range(1, 7)),
        ('Pho', 'GT0', 'Ft2'): list(range(1, 7)),
        ('Pho', 'GT0', 'Ft3'): list(range(1, 7)),
        ('Pho', 'GT0', 'Ft4'): list(range(1, 7)),
        ('Pho', 'GT1', 'Ft1'): list(range(1, 7)),
        ('Pho', 'GT1', 'Ft2'): list(range(1, 7)),
        ('Pho', 'GT1', 'Ft3'): [1, 2, 3, 5],
        ('Pho', 'GT1', 'Ft4'): list(range(1, 6)),
        ('Pho', 'GT5', 'Ft1'): list(range(1, 7)),
        ('Pho', 'GT5', 'Ft2'): list(range(1, 7)),
        ('Pho', 'GT5', 'Ft3'): list(range(1, 7)),
        ('Pho', 'GT5', 'Ft4'): list(range(1, 7))}
dSCFt = {(sTp, sGT, sFt): [dSFt[sFt] + S_USC + dSGT[sGT] + S_USC + str(k)
                           for k in dKFt[(sTp, sGT, sFt)]]
         for (sTp, sGT, sFt) in dKFt}

# --- FUNCTIONS ---------------------------------------------------------------
def getDVal1L(pdSer, dSFt, sTp, sGT):
    return {sFt: pdSer.loc[dSCFt[(sTp, sGT, sFt)]].to_list() for sFt in dSFt}

def calcSimpleStats(x, nMin=2):
    cVar, cSD = np.nan, np.nan
    t = stats.describe(x, axis=None, nan_policy='omit')
    nObs, cMean = t[0], t[2]
    if nObs >= nMin:
        cVar, cSD = t[3], np.sqrt(t[3])
    return nObs, cMean, cVar, cSD

def calcValsCF(pCSV, sF, sRs, dSFt, sTp, sGT, sSp):
    print('_'*8, 'Starting with', (sTp, sGT), '_'*8)
    pF = os.path.join(pCSV, sF)
    if os.path.exists(pF):
        cDfr = pd.read_csv(pF, sep=sSp)
        lKWSig = []
        for iL in cDfr.index:
            lV2Add = []
            dV1L = getDVal1L(cDfr.loc[iL, :], dSFt, sTp, sGT)
            
            pKW = stats.kruskal(*dV1L.values(), nan_policy='omit')[1]
            lKWSig.append(pKW)
        cDfr[S_KW] = lKWSig
        saveAsCSV(cDfr, pF, sRs=sRs, sSp=sSp)
    else:
        print('Path', pF, 'corresponding to', (sTp, sGT), 'does not exist!')

def saveAsCSV(pdDfr, pFC, sRs, sSp, wrtNmF=True):
    pFN = pFC[:-len(S_DOT + S_CSV)] + S_USC + sRs + S_DOT + S_CSV
    pdDfr.to_csv(pFN, sep=sSp)
    if wrtNmF:
        print('Saved results to', pFN, '...')

# --- MAIN --------------------------------------------------------------------
for ((sTp, sGT), sF) in dSF.items():
    calcValsCF(pCSV, sF, sRes, dSFt, sTp, sGT, sSp=sSep)

print('+'*30, 'DONE.', '+'*30)

###############################################################################

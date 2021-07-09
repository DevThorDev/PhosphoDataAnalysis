# -*- coding: utf-8 -*-
###############################################################################
# --- CalcOverRepFromCSV.py ---------------------------------------------------
###############################################################################
import os, time

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as stats

# --- CONSTANTS ---------------------------------------------------------------
S_USC = '_'
S_DOT = '.'
S_CSV = 'csv'

S_O = '0'
S_1 = '1'
S_2 = '2'
S_3 = '3'
S_4 = '4'
S_5 = '5'

S_BIN_CODE_S = 'BC'
S_BIN_CODE_L = 'BinCode'
S_BIN_CODE_S_2 = S_BIN_CODE_S + S_2
S_BIN_CODE_L_2 = S_BIN_CODE_L + S_2

S_GT0 = 'GT' + S_O
S_GT1 = 'GT' + S_1
S_GT5 = 'GT' + S_5
L_S_GT = [S_GT0, S_GT1, S_GT5]

S_EXT_CSV = 'csv'
S_EXT_PDF = 'pdf'

S_IDX = 'Idx'
S_COL = 'Col'

S_INP_DATA = 'InputData'
S_OVER_REP = 'OverRep'
S_M_CORR_BON = 'Bonferroni'

R04 = 4

# --- INPUT -------------------------------------------------------------------
nMinI, nMaxI = 1, 832

sFInp_IDist = 'Input__OvRep_PhoD_GTX_AllD__BinCode2_1_832_IDist0'
sFOut_IDist = 'Result__OvRep_PhoD_GTX_AllD__BinCode2_1_832_IDist0'

pCSV = os.path.join('..', '..', '..', '25_Papers', '01_FirstAuthor',
                    '04_SysBio_DataAnalysis', '80_ResultsCSV')
sSep = ';'

# --- INPUT DICTIONARY --------------------------------------------------------
dInput = {'sBC_S': S_BIN_CODE_S,
          'sBC_L': S_BIN_CODE_L,
          'sBC2_S': S_BIN_CODE_S_2,
          'sBC2_L': S_BIN_CODE_L_2,
          'lSGT': L_S_GT,
          'sExtCSV': S_EXT_CSV,
          'sExtPDF': S_EXT_PDF,
          'sIdx': S_IDX,
          'sCol': S_COL,
          'sInpDat': S_INP_DATA,
          'sOvRep': S_OVER_REP,
          'sMCorr': S_M_CORR_BON,
          'R04': R04,
          'nMinI': nMinI,
          'nMaxI': nMaxI,
          'pCSV': pCSV,
          'sFInp': sFInp_IDist + S_DOT + S_CSV,
          'sFOut': sFOut_IDist + S_DOT + S_CSV,
          'pFInp': os.path.join(pCSV, sFInp_IDist + S_DOT + S_CSV),
          'pFOut': os.path.join(pCSV, sFOut_IDist + S_DOT + S_CSV),
          'sSep': sSep}

# --- FUNCTIONS ---------------------------------------------------------------
# --- General functions -------------------------------------------------------
def createDir(pF):
    if not os.path.isdir(pF):
        os.mkdir(pF)

def joinToPath(pF='', nmF='Dummy.txt'):
    if len(pF) > 0:
        createDir(pF)
        return os.path.join(pF, nmF)
    else:
        return nmF

def changePath(pFC, pN):
    return os.path.join(joinToPath(pN, os.path.basename(pFC)))

def modFName(pF, sPre='', sPost=''):
    lPF, nmF = pF.split('.'), os.path.basename(pF)
    nmF = sPre + nmF[:-(len(lPF[-1]) + 1)] + sPost + nmF[-(len(lPF[-1]) + 1):]
    return os.path.join(os.path.dirname(pF), nmF)

def changeFExt(pF, newX=S_EXT_PDF):
    lPF = pF.split('.')
    return pF[:-len(lPF[-1])] + newX

def addString(sRetY, sPre = '', sPost = '', sRetN = ''):
    if sRetY is not None:
        if len(sRetY) > 0:
            return sPre + sRetY + sPost
    return sRetN

def adaptPF4Plot(pFC, pN, sPre='', sPost='', newX=S_EXT_PDF):
    return changePath(changeFExt(modFName(pFC, sPre, sPost), newX), pN)

def readCSV(pF, sepD=',', iCol=None, dDtype=None, lHStr=[]):
    if dDtype is None and lHStr is not None:
        dDtype = {s: str for s in lHStr}
    return pd.read_csv(pF, sep=sepD, index_col=iCol, dtype=dDtype)

def reIndexSpec(d, srtMd='str', srtDg=0):
    if srtMd == 'str':
        d = d.reindex(sorted(d.columns, key=lambda x: x[srtDg:]), axis=1)
    elif srtMd == 'float':
        d = d.reindex(sorted(d.columns, key=lambda x: float(x[srtDg:])),
                      axis=1)
    elif srtMd == 'int':
        d = d.reindex(sorted(d.columns, key=lambda x: int(x[srtDg:])), axis=1)
    else:
        d = d.sort_index(axis=1)
    return d

def saveDictAsPdDfr(dI, lPF, dSrt, cD, nMin, nMax, k=0):
    cDfr = pd.DataFrame(cD, index=range(nMin, nMax + 1))
    cDfr.columns = cDfr.columns.astype(str)
    if dI['dTpY'][dI['lTpY'][k]][2]:
        cDfr = -np.log10(cDfr)
    srtMd, srtDg = dSrt[S_COL]['Srt']
    cDfr = reIndexSpec(cDfr, srtMd, srtDg)
    cDfr = cDfr.round(dI['dTpY'][dI['lTpY'][k]][1])
    cDfr.to_csv(lPF[k], sep=dI['cSep'])
    return cDfr

def printElapsedTimeSim(stT, cT, sPre='Time'):
    # calculate and display elapsed time
    elT = round(cT - stT, R04)
    print(sPre, 'elapsed:', elT, 'seconds, this is', round(elT/60, R04),
          'minutes or', round(elT/3600, R04), 'hours or',
          round(elT/(3600*24), R04), 'days.')

# --- Specific functions ------------------------------------------------------
def getPF(dI, nmF, nmFExt=S_EXT_CSV):
    return joinToPath(dI['pRelResF'], nmF + '.' + nmFExt)

def getNmF(dOIn, sF='sFNm', sSt='', sMn='', sEnd='', sPtSt='__', sPrMn='__',
           sPrEnd='__'):
    nmF = addString(sSt, sPost=sPtSt) + dOIn[sF]
    nmF += addString(sMn, sPre=sPrMn)
    nmF += addString(sEnd, sPre=sPrEnd)
    return nmF

def getPFRes(dI, dOIn, sSt='', sMn='', sEnd='', sMod='No',
             nmFExt=S_EXT_CSV):
    if 'calcForTrD' in dI:
        if dI['calcForTrD']:    # do calculations for the transformed data
            sMod = 'Tr'
    if sMod == 'No' or dOIn['isClRD']:  # ClRD is always transformed
        nmF = getNmF(dOIn, sF='sFNm', sSt=sSt, sMn=sMn, sEnd=sEnd)
    else:   # transformed (sMod == 'Tr') and deviation (sMod == 'Dv') data
        nmF = getNmF(dOIn, sF='sFNm' + sMod, sSt=sSt, sMn=sMn, sEnd=sEnd)
    return getPF(dI, nmF, nmFExt)

def iniKeyV1D(lD, cKey, lenArr, convInt=None):
    for k, cD in enumerate(lD):
        if convInt is None:
            cD[cKey] = np.zeros((lenArr))
        else:
            assert len(convInt) == len(lD)
            if convInt[k]:
                cD[cKey] = np.zeros((lenArr), dtype=np.uint64)
            else:
                cD[cKey] = np.zeros((lenArr))

def getParProf(dI, iTpO, isClR, lElCtDef):
    lElCt = None
    if iTpO == 0:
        sID = 'PD'
    elif iTpO == 1:
        sID = 'BO'
        if isClR:            # a "BinOps" from cluster result data object
            sID = 'CR'
    cDSrt, lElCt = dI['dSrt' + sID], dI['lElC' + sID]
    nMn, nMx = dI['nMin' + sID], dI['nMax' + sID]
    lSrt = list(cDSrt[S_IDX])
    lAsc = [cDSrt[S_IDX][cK]['Asc'] for cK in cDSrt[S_IDX]]
    if lElCt is None:    # no category to evaluate defined
        lElCt = lElCtDef
    return cDSrt, lSrt, lAsc, lElCt, nMn, nMx

def getCntTblF(dOccAbs, nCur, nMin, nOccAbs, nTot, cArr, cAttr):
    nOccCur = len(cArr[cArr == cAttr])
    dOccAbs[cAttr][nCur - nMin] = nOccCur
    a00, a01 = nCur - nOccCur, nOccCur
    a10, a11 = nTot - nCur - nOccAbs + nOccCur, nOccAbs - nOccCur
    if min([a00, a01, a10, a11]) < 0:
        print('ERROR: [a00, a01, a10, a11] =', [a00, a01, a10, a11])
        assert False
    return np.array([[a00, a01], [a10, a11]], dtype=np.uint64)

def exactFisher(arrDat, cAlt='two-sided'):
    arrDat = np.array(arrDat, dtype=np.uint64)
    assert arrDat.shape == (2, 2)
    # returns (oddsRatio, pVal)
    return stats.fisher_exact(arrDat, alternative=cAlt)

def calcPValsF(dPValOv, dPValUn, arrCTblF, nCur, nMin, nCAttr, cAttr,
               sCorrect=None):
    pValOv = exactFisher(arrCTblF, cAlt='less')[1]
    pValUn = exactFisher(arrCTblF, cAlt='greater')[1]
    if sCorrect == S_M_CORR_BON:    # Bonferroni correction
        pValOv *= nCAttr
        pValUn *= nCAttr
    dPValOv[cAttr][nCur - nMin] = pValOv
    dPValUn[cAttr][nCur - nMin] = pValUn

def calcPValProfiles(dI, dOccAbs, dPValOv, dPValUn, dfrRd, lSerVC, llAttr,
                     lNAttr, nMin, nMax, nTot, lElCt):
    if lElCt is not None:
        assert len(llAttr) == len(lNAttr) and len(llAttr) == len(lElCt)
    nCalc, lToInt = sum([len(l) for l in llAttr]), [True, False, False]
    for i in range(len(llAttr)):
        for k, cAt in enumerate(llAttr[i]):
            iniKeyV1D([dOccAbs, dPValOv, dPValUn], cAt, lenArr=nMax+1-nMin,
                      convInt=lToInt)
            for nCur in range(nMin, nMax + 1):
                if lElCt is not None:
                    arrS = dfrRd.iloc[:nCur, :].loc[:, lElCt[i]].values
                else:
                    arrS = dfrRd.iloc[:nCur, :].index.to_numpy()
                arrCntTblF = getCntTblF(dOccAbs, nCur, nMin, lSerVC[i].at[cAt],
                                        nTot, arrS, cAt)
                calcPValsF(dPValOv, dPValUn, arrCntTblF, nCur, nMin,
                           lNAttr[i], cAt, sCorrect=dI['sMCorrectL'])
            if (k + 1)%10 == 0 or k + 1 == nCalc:
                print('Done:', k + 1, 'of', nCalc)

def getSrtData(pF, lSrt, lAsc, lSCat=None, sepD=',', idxCol=0):
    pdDfr = readCSV(pF, sepD=sepD, iCol=idxCol)
    if lSCat == None:
        lSerVC = [cSer.value_counts(sort=False) for cSer in [pdDfr.index]]
        llAttr = [list(set(pdDfr.index))]
    else:
        lSer = [pdDfr.loc[:, s] for s in lSCat]
        lSerVC = [cSer.value_counts(sort=False) for cSer in lSer]
        llAttr = [list(set(pdDfr.loc[:, s])) for s in lSCat]
    lLenLAttr = [len(llAttr[k]) for k in range(len(llAttr))]
    pdDfr = pdDfr.sort_values(by=lSrt, ascending=lAsc)
    return pdDfr, lSerVC, llAttr, lLenLAttr, pdDfr.shape[0]

def selDataThr(dI, cDfr, dSrt, k=0):
    dfrSelD = pd.DataFrame(index=cDfr.index, columns=cDfr.columns)
    cThr = dI['thrProf']
    if dI['dTpY'][dI['lTpY'][k]][2]:
        cThr = -np.log10(cThr)
    for sC in cDfr.columns:
        cSer = cDfr.loc[:, sC]
        if dI['sComp'] == '>':
            dfrSelD.loc[:, sC] = cSer[cSer > cThr]
        elif dI['sComp'] == '>=':
            dfrSelD.loc[:, sC] = cSer[cSer >= cThr]
        elif dI['sComp'] == '==':
            dfrSelD.loc[:, sC] = cSer[cSer == cThr]
        elif dI['sComp'] == '<=':
            dfrSelD.loc[:, sC] = cSer[cSer <= cThr]
        elif dI['sComp'] == '<':
            dfrSelD.loc[:, sC] = cSer[cSer < cThr]
        else:
            dfrSelD.loc[:, sC] = cSer
    dfrSelD = dfrSelD.dropna(axis=1, how='all')
    srtMd, srtDg = dSrt[S_COL]['Srt']
    return reIndexSpec(dfrSelD, srtMd, srtDg)

# --- Plot functions ----------------------------------------------------------
def plotProfile(dI, cDfr, pF, k = 0, tpPr = 'PD'):
    i, j = 0, 0
    while j < len(cDfr.columns):
        d = cDfr.iloc[:, j:(j + dI['jIncr'])]
        i += dI['iIncr']
        j += dI['jIncr']
        pFN = adaptPF4Plot(pF, dI['pRelPltF'],
                           sPre = dI['nmPlt_Prof'] + str(i) + '_')
        if not os.path.isfile(pFN):
            cFig, cAx = plt.subplots()
            for sC in d.columns:
                if tpPr == 'PD' and sC in dI['dClrBinC']:
                    cAx.plot(d.loc[:, sC], lw = 0.75, label = sC,
                             color = dI['dClrBinC'][sC])
                else:
                    cAx.plot(d.loc[:, sC], lw = 0.75, label = sC)
            cAx.set_xlabel(dI['dTpX'][tpPr])
            cAx.set_ylabel(dI['dTpY'][dI['lTpY'][k]][0])
            l = cAx.legend(loc = 'center',
                           bbox_to_anchor = dI['coordAnchorBox'],
                           fontsize = dI['szFontLeg'])
            if l is not None:
                cFig.savefig(pFN, bbox_extra_artists = (l,),
                             bbox_inches = 'tight')
            else:
                cFig.savefig(pFN)
            plt.close()

# --- CLASSES -----------------------------------------------------------------
class InputData():
    def __init__(self, dInp):
        self.idO = dInp['sInpDat']
        self.descO = 'Input data class'
        self.dI = dInp
        self.fillInp()
        print('Initiated "InputData" base object.')
    
    def fillInp(self):
        for sK, cV in self.dI.items():
            setattr(self, sK, cV)
        print('Set InputData attributes.')

class OverRep():
    def __init__(self, InpD):
        self.idO = InpD.sOvRep
        self.descO = 'Over-representation'
        self.inpD = InpD
        self.loadDfrInp()
        # self.getPResF()
        print('Initiated "OverRep" base object.')

    def printIDDesc(self):
        print('Object ID:', self.idO)
        print('Object description:', self.descO)
    
    def printObjInfo(self):
        print('-'*20, 'Object', self.descO, '(ID', self.idO, ')', '-'*20)
        print('Input dictionary:', self.dI)
    
    def loadDfrInp(self):
        self.dfrIn = pd.read_csv(self.inpD.pFInp, sep=self.inpD.sSep,
                                 dtype={self.inpD.sBC2_L: str})
    
    def printDfrInp(self):
        print(self.dfrIn)

    def getPResF(self, sMn='', sMod='No'):
        self.lPRF, lDef = [], None
        if self.iTpOD == 1:
            lDef = [self.OD.lOD[k].idO for k in range(2)]
        t = getParProf(self.dI, self.iTpOD, self.OD.isClRD, lDef)
        self.dSort, self.lSrt, self.lAsc, self.lElCat, self.nMin, self.nMax = t
        sEPr = str(self.nMin) + '_' + str(self.nMax) + '_'
        for sHd in self.dSort[S_IDX]:
            sEPr += sHd + str(int(self.dSort[S_IDX][sHd]['Asc'])) + '_'
        for sCat in self.lElCat:
            sEPr = addString(str(sCat), sPost='_') + sEPr
        for sE in [self.dI['sNOcc'], self.dI['sPValOvPOf'],
                   self.dI['sPValUnPOf']]:
            sE = addString(sEPr + self.dI['sMCorrectS'], sPost='_') + sE
            self.lPRF.append(getPFRes(self.dI, self.dOIn, sMn=sMn, sEnd=sE,
                                      sMod=sMod))

    def plotProfiles(self, dfrR, k=0, tpPr='PD'):
        d = selDataThr(self.dI, dfrR, self.dSort, k)
        plotProfile(self.dI, d, self.lPRF[k], k, tpPr=tpPr)

    def calcProfiles(self, tpPr='PD', useMn=True):
        pdDfr = readCSV(self.dI['pToF'], sepD=self.dI['cSep'], iCol=None)
        srtDfr = pdDfr.sort_values(by=self.dI['lSrt'], ascending=self.dI['lAsc'])
        
        [dOccAbs, dPValOv, dPValUn] = [{}, {}, {}]
        t = getSrtData(self.getPFData(useMn=useMn), self.lSrt, self.lAsc,
                       self.lElCat, sepD=self.dI['cSep'])
        dfrRd, lSerVC, llAttr, lNAttr, N = t
        calcPValProfiles(self.dI, dOccAbs, dPValOv, dPValUn, srtDfr,
                            lSerVC, llAttr, lNAttr, self.nMin, self.nMax, N,
                            self.lElCat)
        for k, cD in enumerate([dOccAbs, dPValOv, dPValUn]):
            dfrR = saveDictAsPdDfr(self.dI, cD, self.nMin, self.nMax, k)
            self.plotProfiles(dfrR, k, tpPr=tpPr)

# --- MAIN --------------------------------------------------------------------
startTime = time.time()
print('+'*25 + ' START', time.ctime(startTime), '+'*25)

cOvRepAnalysis = OverRep(dInput)

print('-'*80)
printElapsedTimeSim(startTime, time.time(), 'Total time')
print('+'*37, 'DONE.', '+'*36)

###############################################################################

# -*- coding: utf-8 -*-
###############################################################################
# --- CalcOverRepFromCSV.py ---------------------------------------------------
# Calculate the over-representation profile of an entity (e.g. a metabolite
# or a bin) in a table sorted according to a numeric column
# [used for "IDist"; "Input__OvRep_PhoD_GTX_AllD__BinCode2_1_832_IDist0"]
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

S_BASE_CL = 'BaseClass'
S_INP_DATA = 'InputData'
S_OVER_REP = 'OverRep'

S_N_OCC_ABS = 'nOccurrAbs'
S_P_VAL_OV = 'pValFOver'
S_P_VAL_UN = 'pValFUnder'
S_P_OF = 'p'
S_M_CORR_BON = 'Bonferroni'
S_IDIST = 'IDist'
S_N_OCC = 'NOcc'
S_OVER_REP = 'ORp'
S_UNDER_REP = 'URp'
S_Y_N_OCC = 'Number of occurrences ($\it{k}$)'
S_Y_P_VAL = '$-\log_{10}$(p-value)'

R04 = 4
R06 = 6

# --- INPUT -------------------------------------------------------------------
# --- data specific input -----------------------------------------------------
reverseIt = False

nMinIDist, nMaxIDist = 1, 832
nMinIDist_rev, nMaxIDist_rev = 1, 289

dSrtIDist = {S_IDX: {S_IDIST: {'Asc': False}},
              S_COL: {'Srt': ('float', 0), 'Asc': True}}
dSrtIDist_rev = {S_IDX: {S_IDIST: {'Asc': True}},
                 S_COL: {'Srt': ('float', 0), 'Asc': True}}

lElCIDist = [S_BIN_CODE_L_2]

sMCorrectL = S_M_CORR_BON    # None / S_M_CORR_BON
sMCorrectS = sMCorrectL[:3]

sSep = ';'

# --- profile-type specific input ---------------------------------------------
lTpX = [S_IDIST]
lTpY = [S_N_OCC, S_OVER_REP, S_UNDER_REP]

lSXAx = ['Top $\it{n}$ of the highest distance indices']
lSXAx_rev = ['Top $\it{n}$ of the lowest distance indices']

lSYAx = [S_Y_N_OCC, S_Y_P_VAL, S_Y_P_VAL]
lNDigRndYAx = [0, R06, R06]
lDoPYAx = [False, True, True]

# --- names and paths of files and dirs ---------------------------------------
sFInp_IDist = 'Input__OvRep_PhoD_GTX_AllD__BinCode2_1_832_IDist0'
sFInp_IDist_rev = 'Input__OvRep_PhoD_GTX_AllD__BinCode2_1_289_IDist1'

sFOut_IDist = 'Result__OvRep_PhoD_GTX_AllD'

sDirCSV = '80_ResultsCSV'
sDirPDF = '82_ResultsPDF'

pBase = os.path.join('..', '..', '..', '25_Papers', '01_FirstAuthor',
                    '04_SysBio_DataAnalysis')
pCSV = os.path.join(pBase, sDirCSV)

# --- graphics parameters -----------------------------------------------------
nmPlt_Prf = 'Profile'   # name prefix of the plot
thrProf = 0.05          # plot threshold for profiles
sComp = '>='            # comparison string (value with threshold)
szFontLeg = 'small'     # font size of legend
iIncr = 1               # increase of file number
jIncr = 10              # number of entities (e.g. metabolites) per plot
posLegXY = (1.1, 0.5)   # coordinates of the legend anchor box

lWdPlt = 0.75
dClrBinC = {'2.1': (0.12, 0.47, 0.71),
            '4.1': (1.0, 0.5, 0.05),
            '17.2': (0.17, 0.63, 0.17),
            '29.2': (0.84, 0.15, 0.16),
            '31.4': (0.58, 0.4, 0.74),
            '33.99': (0.55, 0.34, 0.29)}

# --- derived values ----------------------------------------------------------
if reverseIt:
    nMinIDist, nMaxIDist = nMinIDist_rev, nMaxIDist_rev
    dSrtIDist = dSrtIDist_rev
    lSXAx = lSXAx_rev
    sFInp_IDist = sFInp_IDist_rev

# --- assertions --------------------------------------------------------------
assert len(lSXAx) == len(lTpX)
assert (len(lSYAx) == len(lTpY) and len(lNDigRndYAx) == len(lTpY) and
        len(lDoPYAx) == len(lTpY))

# --- INPUT DICTIONARY --------------------------------------------------------
dInput = {# --- constants
          'sBC_S': S_BIN_CODE_S,
          'sBC_L': S_BIN_CODE_L,
          'sBC2_S': S_BIN_CODE_S_2,
          'sBC2_L': S_BIN_CODE_L_2,
          'lSGT': L_S_GT,
          'sExtCSV': S_EXT_CSV,
          'sExtPDF': S_EXT_PDF,
          'sIdx': S_IDX,
          'sCol': S_COL,
          'sBase': S_BASE_CL,
          'sInpDat': S_INP_DATA,
          'sOvRep': S_OVER_REP,
          'sNOcc': S_N_OCC_ABS,
          'sPValOv': S_P_VAL_OV,
          'sPValUn': S_P_VAL_UN,
          'sPOf': S_P_OF,
          'sPValOvPOf': S_P_VAL_OV + '_' + S_P_OF,
          'sPValUnPOf': S_P_VAL_UN + '_' + S_P_OF,
          'R04': R04,
          'R06': R06,
          # --- data specific input
          'nMin': nMinIDist,
          'nMax': nMaxIDist,
          'dSrt': dSrtIDist,
          'lElC': lElCIDist,
          'sMCorrectL': sMCorrectL,
          'sMCorrectS': sMCorrectS,
          'sSep': sSep,
          # --- profile-type specific input
          'lTpX': lTpX,
          'lTpY': lTpY,
          'dTpX': {lTpX[k]: lSXAx[k] for k in range(len(lTpX))},
          'dTpY': {lTpY[k]: (lSYAx[k], lNDigRndYAx[k], lDoPYAx[k]) for k in
                   range(len(lTpY))},
          # --- names and paths of files and dirs
          'pCSV': pCSV,
          'pPDF': os.path.join(pBase, sDirPDF),
          'pFInp': os.path.join(pCSV, sFInp_IDist + S_DOT + S_CSV),
          'sFOut': sFOut_IDist,
          # --- graphics parameters
          'nmPlt_Prf': nmPlt_Prf,
          'thrProf': thrProf,
          'sComp': sComp,
          'szFontLeg': szFontLeg,
          'iIncr': iIncr,
          'jIncr': jIncr,
          'posLegXY': posLegXY,
          'lWdPlt': lWdPlt,
          'dClrBinC': dClrBinC}

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

def saveDictAsPdDfr(inpD, lPF, dSrt, cD, nMin, nMax, k=0):
    cDfr = pd.DataFrame(cD, index=range(nMin, nMax + 1))
    cDfr.columns = cDfr.columns.astype(str)
    if inpD.dTpY[inpD.lTpY[k]][2]:
        cDfr = -np.log10(cDfr)
    srtMd, srtDg = dSrt[S_COL]['Srt']
    cDfr = reIndexSpec(cDfr, srtMd, srtDg)
    cDfr = cDfr.round(inpD.dTpY[inpD.lTpY[k]][1])
    cDfr.to_csv(lPF[k], sep=inpD.sSep)
    return cDfr

def printElapsedTimeSim(stT, cT, sPre='Time'):
    # calculate and display elapsed time
    elT = round(cT - stT, R04)
    print(sPre, 'elapsed:', elT, 'seconds, this is', round(elT/60, R04),
          'minutes or', round(elT/3600, R04), 'hours or',
          round(elT/(3600*24), R04), 'days.')

# --- Specific functions ------------------------------------------------------
def getPF(inpD, nmF, nmFExt=S_EXT_CSV):
    return joinToPath(inpD.pCSV, nmF + '.' + nmFExt)

def getNmF(sF='sFNm', sSt='', sEnd='', sPtSt='__', sPrMn='__', sPrEnd='__'):
    nmF = addString(sSt, sPost=sPtSt) + sF
    nmF += addString(sEnd, sPre=sPrEnd)
    return nmF

def getPFRes(inpD, sSt='', sEnd='', nmFExt=S_EXT_CSV):
    return getPF(inpD, getNmF(sF=inpD.sFOut, sSt=sSt, sEnd=sEnd), nmFExt)

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

def getParProf(inpD):
    cDSrt, lElCt, nMn, nMx = inpD.dSrt, inpD.lElC, inpD.nMin, inpD.nMax
    lSrt = list(cDSrt[inpD.sIdx])
    lAsc = [cDSrt[inpD.sIdx][cK]['Asc'] for cK in cDSrt[inpD.sIdx]]
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

def calcPValProfiles(inpD, dOccAbs, dPValOv, dPValUn, dfrRd, lSerVC, llAttr,
                     lNAttr, nMin, nMax, nTot, lElCt):
    print('-'*20, 'Starting calculation of p-value profiles...')
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
                           lNAttr[i], cAt, sCorrect=inpD.sMCorrectL)
            if (k + 1)%10 == 0 or k + 1 == nCalc:
                print('Done:', k + 1, 'of', nCalc)

def getSrtData(dfrIn, lSrt, lAsc, lSCat=None):
    if lSCat == None:
        lSerVC = [cSer.value_counts(sort=False) for cSer in [dfrIn.index]]
        llAttr = [list(set(dfrIn.index))]
    else:
        lSer = [dfrIn.loc[:, s] for s in lSCat]
        lSerVC = [cSer.value_counts(sort=False) for cSer in lSer]
        llAttr = [list(set(dfrIn.loc[:, s])) for s in lSCat]
    lLenLAttr = [len(llAttr[k]) for k in range(len(llAttr))]
    pdDfr = dfrIn.sort_values(by=lSrt, ascending=lAsc)
    return pdDfr, lSerVC, llAttr, lLenLAttr, pdDfr.shape[0]

def selDataThr(inpD, cDfr, dSort, k=0):
    dfrSelD = pd.DataFrame(index=cDfr.index, columns=cDfr.columns)
    cThr = inpD.thrProf
    if inpD.dTpY[inpD.lTpY[k]][2]:
        cThr = -np.log10(cThr)
    for sC in cDfr.columns:
        cSer = cDfr.loc[:, sC]
        if inpD.sComp == '>':
            dfrSelD.loc[:, sC] = cSer[cSer > cThr]
        elif inpD.sComp == '>=':
            dfrSelD.loc[:, sC] = cSer[cSer >= cThr]
        elif inpD.sComp == '==':
            dfrSelD.loc[:, sC] = cSer[cSer == cThr]
        elif inpD.sComp == '<=':
            dfrSelD.loc[:, sC] = cSer[cSer <= cThr]
        elif inpD.sComp == '<':
            dfrSelD.loc[:, sC] = cSer[cSer < cThr]
        else:
            dfrSelD.loc[:, sC] = cSer
    dfrSelD = dfrSelD.dropna(axis=1, how='all')
    srtMd, srtDg = dSort[S_COL]['Srt']
    return reIndexSpec(dfrSelD, srtMd, srtDg)

# --- Plot functions ----------------------------------------------------------
def plotProfile(inpD, cDfr, pF, k=0, tpPr=S_IDIST):
    i, j = 0, 0
    while j < len(cDfr.columns):
        d = cDfr.iloc[:, j:(j + inpD.jIncr)]
        i += inpD.iIncr
        j += inpD.jIncr
        pFN = adaptPF4Plot(pF, inpD.pPDF, sPre=inpD.nmPlt_Prf + str(i) + S_USC)
        if not os.path.isfile(pFN):
            cFig, cAx = plt.subplots()
            for sC in d.columns:
                if tpPr == S_IDIST and sC in inpD.dClrBinC:
                    cAx.plot(d.loc[:, sC], lw=inpD.lWdPlt, label=sC,
                             color=inpD.dClrBinC[sC])
                else:
                    cAx.plot(d.loc[:, sC], lw=inpD.lWdPlt, label=sC)
            cAx.set_xlabel(inpD.dTpX[tpPr])
            cAx.set_ylabel(inpD.dTpY[inpD.lTpY[k]][0])
            l = cAx.legend(loc='center', bbox_to_anchor=inpD.posLegXY,
                           fontsize=inpD.szFontLeg)
            if l is not None:
                cFig.savefig(pFN, bbox_extra_artists=(l,),
                             bbox_inches='tight')
            else:
                cFig.savefig(pFN)
            plt.close()

# --- CLASSES -----------------------------------------------------------------
class BaseClass():
    def __init__(self):
        self.idO = S_BASE_CL
        self.descO = 'Base class'
        print('Initiated "Base" base object.')

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

class OverRep(BaseClass):
    def __init__(self, InpD):
        super().__init__()
        self.idO = InpD.sOvRep
        self.descO = 'Over-representation'
        self.inpD = InpD
        self.loadDfrInp()
        self.getPResF()
        print('Initiated "OverRep" base object.')

    def printIDDesc(self):
        print('Object ID:', self.idO)
        print('Object description:', self.descO)

    def printObjInfo(self):
        print('-'*20, 'Object', self.descO, '(ID', self.idO, ')', '-'*20)
        print('-'*8, 'Input data:')
        self.inpD.printAttrData()
        print('-'*8, 'Attributes of', self.descO, 'class:')
        self.printAttrData()

    def printDfrInp(self):
        print(self.dfrIn)

    def loadDfrInp(self):
        self.dfrIn = pd.read_csv(self.inpD.pFInp, sep=self.inpD.sSep,
                                 dtype={self.inpD.sBC2_L: str})

    def getPResF(self):
        self.lPRF, sIdx = [], self.inpD.sIdx
        t = getParProf(self.inpD)
        self.dSort, self.lSrt, self.lAsc, self.lElCat, self.nMin, self.nMax = t
        sEPr = str(self.nMin) + S_USC + str(self.nMax) + S_USC
        for sHd in self.dSort[sIdx]:
            sEPr += sHd + str(int(self.dSort[sIdx][sHd]['Asc'])) + S_USC
        for sCat in self.lElCat:
            sEPr = addString(str(sCat), sPost=S_USC) + sEPr
        for sE in [self.inpD.sNOcc, self.inpD.sPValOvPOf,
                   self.inpD.sPValUnPOf]:
            sE = addString(sEPr + self.inpD.sMCorrectS, sPost=S_USC) + sE
            self.lPRF.append(getPFRes(self.inpD, sEnd=sE))

    def plotProfiles(self, dfrR, k=0, tpPr=S_IDIST):
        d = selDataThr(self.inpD, dfrR, self.dSort, k)
        plotProfile(self.inpD, d, self.lPRF[k], k, tpPr=tpPr)

    def calcProfiles(self, tpPr=S_IDIST):
        srtDfr = self.dfrIn.sort_values(by=self.lSrt, ascending=self.lAsc)
        [dOccAbs, dPValOv, dPValUn] = [{}, {}, {}]
        t = getSrtData(self.dfrIn, self.lSrt, self.lAsc, self.lElCat)
        dfrRd, lSerVC, llAttr, lNAttr, N = t
        calcPValProfiles(self.inpD, dOccAbs, dPValOv, dPValUn, srtDfr, lSerVC,
                         llAttr, lNAttr, self.nMin, self.nMax, N, self.lElCat)
        for k, cD in enumerate([dOccAbs, dPValOv, dPValUn]):
            dfrR = saveDictAsPdDfr(self.inpD, self.lPRF, self.dSort, cD,
                                   self.nMin, self.nMax, k)
            self.plotProfiles(dfrR, k, tpPr=tpPr)

# --- MAIN --------------------------------------------------------------------
startTime = time.time()
print('+'*25 + ' START', time.ctime(startTime), '+'*25)
if reverseIt:
    print('-'*8, 'Reverse mode', '-'*8)

inpDat = InputData(dInput)
cOvRepAnalysis = OverRep(inpDat)
cOvRepAnalysis.printObjInfo()
cOvRepAnalysis.calcProfiles()

print('-'*80)
printElapsedTimeSim(startTime, time.time(), 'Total time')
print('+'*37, 'DONE.', '+'*36)

###############################################################################

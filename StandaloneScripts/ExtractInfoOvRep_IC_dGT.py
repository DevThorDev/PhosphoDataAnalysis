# -*- coding: utf-8 -*-
###############################################################################
# --- ExtractInfoOvRep_IC_dGT.py ----------------------------------------------
###############################################################################
import os, time

import pandas as pd

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
sFInp_IC_Met_Pho = 'IC_Met_Pho'
sFInp_dGT_Met = 'DistGT_Met'
sFInp_dGT_Pho = 'DistGT_Pho'

sFOut = 'ExtractedInfoOvRep_IC_dGT'

sDirCSV = '51_CSV_DistGT'
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
coordAnchorBox = (1.1, 0.5)         # coordinates of the legend anchor box

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
          'coordAnchorBox': coordAnchorBox,
          'lWdPlt': lWdPlt,
          'dClrBinC': dClrBinC}

# --- FUNCTIONS ---------------------------------------------------------------
def addToDictD(cD, cKMain, cKSub, cV):
    if cKMain in cD:
        assert cKSub not in cD[cKMain]
        cD[cKMain][cKSub] = cV
    else:
        cD[cKMain] = {}
        cD[cKMain][cKSub] = cV

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

def calcRIs(pdDfr, lSHdr, sHdrRef, sIOrig=S_I_ORIG, sUSC=S_USC):
    assert pdDfr.columns.to_list()[:(len(lSHdr) + 1)] == [sHdrRef] + lSHdr
    dRI, nL = {}, pdDfr.shape[0]
    pdDfr[sIOrig] = pdDfr.index.to_list()
    for sHdr in lSHdr:
        pdDfr.sort_values(by=[sHdr, sIOrig], ascending=[False, True],
                          inplace=True, kind='stable', ignore_index=True)
        pdDfr.loc[:, sHdr] = getLVals(pdDfr, sHdr, nEl=nL)
        sSpec, sTp, sGT = sHdr.split(sUSC)
        addToDictD(dRI, cKMain=(sSpec, sGT), cKSub=sTp, cV=sHdr)
    return transcrDict2Dfr(dRI, pdDfr, [sIOrig, sHdrRef])

def loopInpDataFrames(dInp):
    for cTpDat in dInp[S_GENERAL]['lTpDat']:
        print('Processing data of type', cTpDat, '...')
        cDfrV = pd.read_csv(dInp[cTpDat]['pFInp'], sep=dInp[S_GENERAL]['sSep'],
                            dtype={dInp[cTpDat]['sHdrRef']: str})
        cDfrR = calcRIs(cDfrV, lSHdr=dInp[cTpDat]['lSHdr'],
                        sHdrRef=dInp[cTpDat]['sHdrRef'])
        cDfrR.sort_values(by=S_I_ORIG, ascending=True, inplace=True,
                          kind='stable', ignore_index=True)
        cDfrR.to_csv(dInp[cTpDat]['pFOut'], sep=dInp[S_GENERAL]['sSep'])

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

class ExtractedInfo(BaseClass):
    def __init__(self, InpD):
        super().__init__()
        self.idO = InpD.ExtrInfo
        self.descO = 'Extracted info'
        self.inpD = InpD
        self.loadDfrInp()
        self.getPResF()
        print('Initiated "ExtractedInfo" base object.')

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

loopInpDataFrames(dInput)

print('-'*80)
printElapsedTimeSim(startTime, time.time(), 'Total time')
print('+'*37, 'DONE.', '+'*36)

###############################################################################

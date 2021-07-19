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

S_BIN_CODE_L = 'BinCode'
S_BIN_CODE_L_2 = S_BIN_CODE_L + S_2

S_GT0 = 'GT' + S_O
S_GT1 = 'GT' + S_1
S_GT5 = 'GT' + S_5
L_S_GT = [S_GT0, S_GT1, S_GT5]

S_EXT_CSV = 'csv'
S_EXT_PDF = 'pdf'

S_IDX = 'Idx'
S_COL = 'Col'

S_IC_MET_PHO = 'ICMetPho'
S_D_GT_MET = 'dGTMet'
S_D_GT_PHO = 'dGTPho'

S_SRT_BY = 'SortedBy'
S_ASC = 'Asc'
S_IC = 'IC'
S_D_GT = 'd_GT'
S_MIN = 'min'
S_MAX = 'max'

S_BASE_CL = 'BaseClass'
S_INP_DATA = 'InputData'
S_EXTR_INFO = 'ExtrInfo'

S_MET = 'Metabolite'
S_PHO = 'Phosphopeptide'
S_REMCOL1 = 'PearsonCorr'

R04 = 4
R06 = 6

# --- INPUT -------------------------------------------------------------------
# --- general input -----------------------------------------------------------
modDisp = 1000

# --- data specific input -----------------------------------------------------
dISort = {S_IC_MET_PHO: {S_SRT_BY: S_IC, S_ASC: False},
          S_D_GT_MET: {S_SRT_BY: S_D_GT, S_ASC: False},
          S_D_GT_PHO: {S_SRT_BY: S_D_GT, S_ASC: False}}
# dThr = {S_IC_MET_PHO: {S_MIN: 7.25, S_MAX: None},
#         S_D_GT_MET: {S_MIN: 0.6, S_MAX: None},
#         S_D_GT_PHO: {S_MIN: 0.6, S_MAX: None}}
dThr = {S_IC_MET_PHO: {S_MIN: 4.0, S_MAX: None},
        S_D_GT_MET: {S_MIN: 0.5, S_MAX: None},
        S_D_GT_PHO: {S_MIN: 0.5, S_MAX: None}}

sSep = ';'

# --- profile-type specific input ---------------------------------------------
lTpX = []
lTpY = []

lSXAx = []
# lSXAx = ['Top $\it{n}$ of the highest distance indices']
# lSXAx_rev = ['Top $\it{n}$ of the lowest distance indices']

# lSYAx = [S_Y_N_OCC]
# lNDigRndYAx = [R06]
# lDoPYAx = [True]
lSYAx = []
lNDigRndYAx = []
lDoPYAx = []

# --- names and paths of files and dirs ---------------------------------------
sFInp_IC_Met_Pho = 'IC_Met_Pho'
sFInp_dGT_Met = 'DistGT_Met'
sFInp_dGT_Pho = 'DistGT_Pho'

sFOutS = 'ExtractedInfoOvRepShort'
sFOutF = 'ExtractedInfoOvRepFull'

sDirInCSV = '51_CSV_In_DistGT'
sDirOutCSV = '52_CSV_Out_DistGT'
sDirOutPDF = '82_ResultsPDF'

pBase = os.path.join('..', '..', '..', '25_Papers', '01_FirstAuthor',
                    '04_SysBio_DataAnalysis')
pInCSV = os.path.join(pBase, sDirInCSV)
pOutCSV = os.path.join(pBase, sDirOutCSV)

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

# --- assertions --------------------------------------------------------------
assert len(lSXAx) == len(lTpX)
assert (len(lSYAx) == len(lTpY) and len(lNDigRndYAx) == len(lTpY) and
        len(lDoPYAx) == len(lTpY))

# --- INPUT DICTIONARY --------------------------------------------------------
dInput = {# --- constants
          'sBC_L': S_BIN_CODE_L,
          'sBC2_L': S_BIN_CODE_L_2,
          'lSGT': L_S_GT,
          'sExtCSV': S_EXT_CSV,
          'sExtPDF': S_EXT_PDF,
          'sIdx': S_IDX,
          'sCol': S_COL,
          'sBase': S_BASE_CL,
          'sInpDat': S_INP_DATA,
          'sExtrInfo': S_EXTR_INFO,
          'sMet': S_MET,
          'sPho': S_PHO,
          'R04': R04,
          'R06': R06,
          # --- general input
          'modDisp': modDisp,
          # --- data specific input
          'dISort': dISort,
          'dThr': dThr,
          'sSep': sSep,
          # --- profile-type specific input
          'lTpX': lTpX,
          'lTpY': lTpY,
          'dTpX': {lTpX[k]: lSXAx[k] for k in range(len(lTpX))},
          'dTpY': {lTpY[k]: (lSYAx[k], lNDigRndYAx[k], lDoPYAx[k]) for k in
                   range(len(lTpY))},
          # --- names and paths of files and dirs
          'pInCSV': pInCSV,
          'pOutCSV': pOutCSV,
          'pOutPDF': os.path.join(pBase, sDirOutPDF),
          'pFInp_IC': os.path.join(pInCSV, sFInp_IC_Met_Pho + S_DOT + S_CSV),
          'pFInp_Met': os.path.join(pInCSV, sFInp_dGT_Met + S_DOT + S_CSV),
          'pFInp_Pho': os.path.join(pInCSV, sFInp_dGT_Pho + S_DOT + S_CSV),
          'sFOutS': sFOutS + S_DOT + S_CSV,
          'sFOutF': sFOutF + S_DOT + S_CSV,
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

def sortDfr(pdDfr, dSrt, sSrtBy, sAsc, srtKind='stable'):
    pdDfr.sort_values(by=dSrt[sSrtBy], ascending=dSrt[sAsc], inplace=True,
                      kind=srtKind)

# def calcRIs(pdDfr, lSHdr, sHdrRef, sIOrig=S_I_ORIG, sUSC=S_USC):
#     assert pdDfr.columns.to_list()[:(len(lSHdr) + 1)] == [sHdrRef] + lSHdr
#     dRI, nL = {}, pdDfr.shape[0]
#     pdDfr[sIOrig] = pdDfr.index.to_list()
#     for sHdr in lSHdr:
#         pdDfr.sort_values(by=[sHdr, sIOrig], ascending=[False, True],
#                           inplace=True, kind='stable', ignore_index=True)
#         pdDfr.loc[:, sHdr] = getLVals(pdDfr, sHdr, nEl=nL)
#         sSpec, sTp, sGT = sHdr.split(sUSC)
#         addToDictD(dRI, cKMain=(sSpec, sGT), cKSub=sTp, cV=sHdr)
#     return transcrDict2Dfr(dRI, pdDfr, [sIOrig, sHdrRef])

# def loopInpDataFrames(dInp):
#     for cTpDat in dInp[S_GENERAL]['lTpDat']:
#         print('Processing data of type', cTpDat, '...')
#         cDfrV = pd.read_csv(dInp[cTpDat]['pFInp'], sep=dInp[S_GENERAL]['sSep'],
#                             dtype={dInp[cTpDat]['sHdrRef']: str})
#         cDfrR = calcRIs(cDfrV, lSHdr=dInp[cTpDat]['lSHdr'],
#                         sHdrRef=dInp[cTpDat]['sHdrRef'])
#         cDfrR.sort_values(by=S_I_ORIG, ascending=True, inplace=True,
#                           kind='stable', ignore_index=True)
#         cDfrR.to_csv(dInp[cTpDat]['pFOut'], sep=dInp[S_GENERAL]['sSep'])

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
        self.idO = InpD.sExtrInfo
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
        dDTp_IC = {sIn: str for sIn in [self.inpD.sBC_L, self.inpD.sBC2_L]}
        dDTp_Pho = {sIn: str for sIn in [self.inpD.sBC2_L]}
        self.dfrIn_IC = pd.read_csv(self.inpD.pFInp_IC, sep=self.inpD.sSep,
                                    dtype=dDTp_IC)
        self.dfrIn_Met = pd.read_csv(self.inpD.pFInp_Met, sep=self.inpD.sSep)
        self.dfrIn_Pho = pd.read_csv(self.inpD.pFInp_Pho, sep=self.inpD.sSep,
                                     dtype=dDTp_Pho)

    def getPResF(self):
        self.dSort, self.dThres = self.inpD.dISort, self.inpD.dThr
        self.lAsc = [cV[S_ASC] for cV in self.inpD.dISort.values()]
        sFOutS, sFOutF = self.inpD.sFOutS, self.inpD.sFOutF
        for sK in self.dSort:
            sFOutS += S_USC + sK + str(int(self.dSort[sK][S_ASC]))
            sFOutF += S_USC + sK + str(int(self.dSort[sK][S_ASC]))
        self.pFOutS = os.path.join(self.inpD.pOutCSV, sFOutS + S_DOT + S_CSV)
        self.pFOutF = os.path.join(self.inpD.pOutCSV, sFOutF + S_DOT + S_CSV)
    
    def applyFilter(self, pdDfr, sHdC, sKThr, sMin=S_MIN, sMax=S_MAX):
        thrMin, thrMax = self.dThres[sKThr][sMin], self.dThres[sKThr][sMax]
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

    def sortAndFiltDfr(self):
        sortDfr(self.dfrIn_IC, self.dSort[S_IC_MET_PHO], S_SRT_BY, S_ASC)
        sortDfr(self.dfrIn_Met, self.dSort[S_D_GT_MET], S_SRT_BY, S_ASC)
        sortDfr(self.dfrIn_Pho, self.dSort[S_D_GT_PHO], S_SRT_BY, S_ASC)
        self.dfrFilt_IC = self.applyFilter(self.dfrIn_IC, S_IC, S_IC_MET_PHO)
        self.dfrFilt_Met = self.applyFilter(self.dfrIn_Met, S_D_GT, S_D_GT_MET)
        self.dfrFilt_Pho = self.applyFilter(self.dfrIn_Pho, S_D_GT, S_D_GT_PHO)

    def fillDfrResS(self):
        dApp = {S_MET: [], S_PHO: [], S_IC: [], S_D_GT_MET: [], S_D_GT_PHO: []}
        lHdCF = (list(self.dfrIn_Met.columns) + list(self.dfrIn_Pho.columns) +
                 list(self.dfrIn_IC.loc[:, S_REMCOL1:].columns))
        self.dfrResS = pd.DataFrame(columns=list(dApp))
        self.dfrResF = pd.DataFrame(columns=lHdCF)
        n, N = 0, self.dfrFilt_Met.shape[0]*self.dfrFilt_Pho.shape[0]
        for i in self.dfrFilt_Met.index:
            sMet = self.dfrFilt_Met.at[i, S_MET]
            dstGTMet = self.dfrFilt_Met.at[i, S_D_GT]
            for j in self.dfrFilt_Pho.index:
                sPho = self.dfrFilt_Pho.at[j, S_PHO]
                dstGTPho = self.dfrFilt_Pho.at[i, S_D_GT]
                cL = self.dfrFilt_IC[(self.dfrFilt_IC[S_MET] == sMet) &
                                     (self.dfrFilt_IC[S_PHO] == sPho)]
                if cL.shape[0] == 1:
                    IC = cL.squeeze(axis=0).at[S_IC]
                    print('TEMP - sMet =', sMet, '- sPho =', sPho, '- IC =', IC)
                    dV = {S_MET: sMet, S_PHO: sPho, S_IC: IC,
                          S_D_GT_MET: dstGTMet, S_D_GT_PHO: dstGTPho}
                    # print('dV:\n', dV)
                    for sK in dApp:
                        dApp[sK].append(dV[sK])
                elif cL.shape[0] > 1:
                    print('ERROR: Shape of selected line =', cL.shape)
                    assert False
                n += 1
                if n%self.inpD.modDisp == 0:
                    print('Processed element', n, 'of', N, '.')
        self.dfrResS.append(dApp, ignore_index=True)
        # self.dfrResF.

    def extractionOfExtremes(self):
        self.sortAndFiltDfr()
        self.fillDfrResS()
        self.dfrResS.to_csv(self.pFOutS, sep=self.inpD.sSep)
        self.dfrResF.to_csv(self.pFOutF, sep=self.inpD.sSep)

# --- MAIN --------------------------------------------------------------------
startTime = time.time()
print('+'*25 + ' START', time.ctime(startTime), '+'*25)

inpDat = InputData(dInput)
cXtrInfo = ExtractedInfo(inpDat)
cXtrInfo.printObjInfo()
cXtrInfo.extractionOfExtremes()

print('-'*80)
printElapsedTimeSim(startTime, time.time(), 'Total time')
print('+'*37, 'DONE.', '+'*36)


###############################################################################

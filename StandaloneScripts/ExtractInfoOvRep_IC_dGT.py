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

S_IC_M_P = 'ICMetPho'
S_D_GT_M = 'dGTMet'
S_D_GT_P = 'dGTPho'

S_SRT_BY = 'SortedBy'
S_ASC = 'Asc'
S_IC = 'IC'
S_D_GT = 'd_GT'
S_MIN = 'min'
S_MAX = 'max'

S_IC_GT0 = S_IC + S_USC + S_GT0
S_IC_GT1 = S_IC + S_USC + S_GT1
S_IC_GT5 = S_IC + S_USC + S_GT5
L_S_IC_GT = [S_IC_GT0, S_IC_GT1, S_IC_GT5]

S_BASE_CL = 'BaseClass'
S_INP_DATA = 'InputData'
S_EXTR_INFO = 'ExtrInfo'


S_MET = 'Metabolite'
S_PHO = 'Phosphopeptide'
L_S_M_P = [S_MET, S_PHO]

S_RMNG_COL1_IC = 'PearsonCorr'
S_NEW_IDX = 'NewIndex'
L_S_NO_GT = [S_MET, S_PHO, 'Protein', 'BinCode', 'BinCode2', 'MapMan',
             'Description', 'SelBinCode']
L_S_ADD_GT = [S_RMNG_COL1_IC, 'SpearmanCorr', 'Pearson_pVal', 'Spearman_pVal',
              'IC_N', 'IC_P', 'IC', 'MetSig5', 'PhoSig5']

R04 = 4

# --- INPUT -------------------------------------------------------------------
# --- general input -----------------------------------------------------------
modDisp = 1000

# --- data specific input -----------------------------------------------------
dISort = {S_IC_M_P: {S_GT0: {S_SRT_BY: S_IC, S_ASC: False},
                     S_GT1: {S_SRT_BY: S_IC, S_ASC: False},
                     S_GT5: {S_SRT_BY: S_IC, S_ASC: False}},
          S_D_GT_M: {S_SRT_BY: S_D_GT, S_ASC: False},
          S_D_GT_P: {S_SRT_BY: S_D_GT, S_ASC: False}}
# dThr = {S_IC_M_P: {S_GT0: {S_MIN: 7.25, S_MAX: None},
#                    S_GT1: {S_MIN: 7.25, S_MAX: None},
#                    S_GT5: {S_MIN: 7.25, S_MAX: None}},
#         S_D_GT_M: {S_MIN: 0.6, S_MAX: None},
#         S_D_GT_P: {S_MIN: 0.6, S_MAX: None}}
dThr = {S_IC_M_P: {S_GT0: {S_MIN: 6.0, S_MAX: None},
                   S_GT1: {S_MIN: 6.0, S_MAX: None},
                   S_GT5: {S_MIN: 6.0, S_MAX: None}},
        S_D_GT_M: {S_MIN: 0.6, S_MAX: None},
        S_D_GT_P: {S_MIN: 0.6, S_MAX: None}}

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
sFInp_IC_M_P = 'IC_Met_Pho'
sFInp_dGT_M = 'DistGT_Met'
sFInp_dGT_P = 'DistGT_Pho'

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

# --- derived values ----------------------------------------------------------
pFIGT0 = os.path.join(pInCSV, sFInp_IC_M_P + S_USC + S_GT0 + S_DOT + S_CSV)
pFIGT1 = os.path.join(pInCSV, sFInp_IC_M_P + S_USC + S_GT1 + S_DOT + S_CSV)
pFIGT5 = os.path.join(pInCSV, sFInp_IC_M_P + S_USC + S_GT5 + S_DOT + S_CSV)

# --- assertions --------------------------------------------------------------
assert len(lSXAx) == len(lTpX)
assert (len(lSYAx) == len(lTpY) and len(lNDigRndYAx) == len(lTpY) and
        len(lDoPYAx) == len(lTpY))

# --- INPUT DICTIONARY --------------------------------------------------------
dInput = {# --- constants
          'sBC_L': S_BIN_CODE_L,
          'sBC2_L': S_BIN_CODE_L_2,
          'lSGT': L_S_GT,
          'sBase': S_BASE_CL,
          'sInpDat': S_INP_DATA,
          'sExtrInfo': S_EXTR_INFO,
          'sMet': S_MET,
          'sPho': S_PHO,
          'R04': R04,
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
          'pFInp_IC' + S_USC + S_GT0: pFIGT0,
          'pFInp_IC' + S_USC + S_GT1: pFIGT1,
          'pFInp_IC' + S_USC + S_GT5: pFIGT5,
          'pFInp_M': os.path.join(pInCSV, sFInp_dGT_M + S_DOT + S_CSV),
          'pFInp_P': os.path.join(pInCSV, sFInp_dGT_P + S_DOT + S_CSV),
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

def applyFilter(pdDfr, sHdC, thrMin, thrMax):
    # if sKThr in [S_D_GT_M, S_D_GT_P]:
    #     pdDfr = dDfrIn[sKThr]
    #     thrMin, thrMax = dThr[sKThr][sMin], dThr[sKThr][sMax]
    # else:
    #     assert sGT is not None
    #     pdDfr = dDfrIn[sKThr][sGT]
    #     thrMin, thrMax = dThr[sKThr][sGT][sMin], dThr[sKThr][sGT][sMax]
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

def getLNewIdx(pdDfr):
    assert (S_MET in pdDfr.columns) and (S_PHO in pdDfr.columns)
    return [sM + S_USC + pdDfr[S_PHO][i] for i, sM in enumerate(pdDfr[S_MET])]

def concDfr(dDfr):
    cDfr = dDfr[L_S_GT[0]]
    lHdCN = L_S_NO_GT + [s + S_USC + L_S_GT[0] for s in cDfr.columns
                         if s not in L_S_NO_GT]
    cDfr.columns = lHdCN
    dfrM = pd.concat([pd.Series(getLNewIdx(cDfr), name=S_NEW_IDX), cDfr],
                     axis=1, verify_integrity=True)
    dfrM.set_index(S_NEW_IDX, inplace=True, verify_integrity=True)
    for sGT in L_S_GT[1:]:
        cDfr = dDfr[sGT]
        lHdCN = L_S_NO_GT + [s + S_USC +sGT for s in cDfr.columns
                             if s not in L_S_NO_GT]
        cDfr.columns = lHdCN
        cDfrM = pd.concat([pd.Series(getLNewIdx(cDfr), name=S_NEW_IDX), cDfr],
                          axis=1, verify_integrity=True)
        cDfrM.set_index(S_NEW_IDX, inplace=True, verify_integrity=True)
        dfrM = pd.concat([dfrM, cDfrM.drop(columns=[S_MET, S_PHO])], axis=1,
                         verify_integrity=True)
    return dfrM

def fillDDat(dDat, dDfrFl, sMet, sPho, dstGTM, dstGTP):
    for sGT in L_S_GT:
        cDfrFl = dDfrFl[S_IC_M_P][sGT]
        cL = cDfrFl[(cDfrFl[S_MET] == sMet) & (cDfrFl[S_PHO] == sPho)]
        if cL.shape[0] == 1:
            IC = cL.squeeze(axis=0).at[S_IC]
            dDat[S_IC_M_P][sGT].append(IC)
        elif cL.shape[0] > 1:
            print('ERROR: Shape of selected line =', cL.shape)
            assert False

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
        self.getProcData()
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

    def getPResF(self):
        self.dSort, self.dT = self.inpD.dISort, self.inpD.dThr
        self.lAsc = [cV[S_ASC] for cV in self.inpD.dISort.values()]
        sFOutS, sFOutF, sEnd = self.inpD.sFOutS, self.inpD.sFOutF, ''
        for sK in self.dSort:
            if sK in [S_D_GT_M, S_D_GT_P]:
                sEnd += S_USC + sK + str(int(self.dSort[sK][S_ASC]))
            else:
                for i, sGT in enumerate(self.dSort[sK]):
                    sAscDsc = str(int(self.dSort[sK][sGT][S_ASC]))
                    if i == 0:
                        sEnd += (S_USC + sK + sAscDsc)
                    else:
                        sEnd += (S_USC + sAscDsc)
        sFOutS, sFOutF = sFOutS + sEnd, sFOutF + sEnd
        self.pFOutS = os.path.join(self.inpD.pOutCSV, sFOutS + S_DOT + S_CSV)
        self.pFOutF = os.path.join(self.inpD.pOutCSV, sFOutF + S_DOT + S_CSV)

    def getInf4Inp(self):
        dDatTp_IC = {sIn: str for sIn in [self.inpD.sBC_L, self.inpD.sBC2_L]}
        dDatTp_P = {sIn: str for sIn in [self.inpD.sBC2_L]}
        self.dDatTp = {S_IC_M_P: dDatTp_IC, S_D_GT_P: dDatTp_P}
        self.dPFInp = {S_IC_M_P: {S_GT0: self.inpD.pFInp_IC_GT0,
                                  S_GT1: self.inpD.pFInp_IC_GT1,
                                  S_GT5: self.inpD.pFInp_IC_GT5},
                       S_D_GT_M: self.inpD.pFInp_M,
                       S_D_GT_P: self.inpD.pFInp_P}

    def loadDfrInp(self):
        # load input DataFrames
        sSep = self.inpD.sSep
        dfrIn_IC_GT0 = pd.read_csv(self.dPFInp[S_IC_M_P][S_GT0], sep=sSep,
                                   dtype=self.dDatTp[S_IC_M_P])
        dfrIn_IC_GT1 = pd.read_csv(self.dPFInp[S_IC_M_P][S_GT1], sep=sSep,
                                   dtype=self.dDatTp[S_IC_M_P])
        dfrIn_IC_GT5 = pd.read_csv(self.dPFInp[S_IC_M_P][S_GT5], sep=sSep,
                                   dtype=self.dDatTp[S_IC_M_P])
        dfrIn_M = pd.read_csv(self.dPFInp[S_D_GT_M], sep=sSep)
        dfrIn_P = pd.read_csv(self.dPFInp[S_D_GT_P], sep=sSep,
                              dtype=self.dDatTp[S_D_GT_P])
        self.dDfrIn = {S_IC_M_P: {S_GT0: dfrIn_IC_GT0,
                                  S_GT1: dfrIn_IC_GT1,
                                  S_GT5: dfrIn_IC_GT5},
                       S_D_GT_M: dfrIn_M,
                       S_D_GT_P: dfrIn_P}
        # get dictionary of column headers
        dDfrIn_IC = self.dDfrIn[S_IC_M_P]
        lH_GT0 = list(dDfrIn_IC[S_GT0].loc[:, S_RMNG_COL1_IC:].columns)
        lH_GT1 = list(dDfrIn_IC[S_GT1].loc[:, S_RMNG_COL1_IC:].columns)
        lH_GT5 = list(dDfrIn_IC[S_GT5].loc[:, S_RMNG_COL1_IC:].columns)
        self.dHdCol = {S_IC_M_P: {S_GT0: [s + S_USC + S_GT0 for s in lH_GT0],
                                  S_GT1: [s + S_USC + S_GT1 for s in lH_GT1],
                                  S_GT5: [s + S_USC + S_GT5 for s in lH_GT5]},
                       S_D_GT_M: list(self.dDfrIn[S_D_GT_M].columns),
                       S_D_GT_P: list(self.dDfrIn[S_D_GT_P].columns)}

    def getProcData(self):
        self.getPResF()
        self.getInf4Inp()
        self.loadDfrInp()

    def simpleFilter(self, sHdC, sKey, sMin=S_MIN, sMax=S_MAX):
        thMin, thMax = self.dT[sKey][sMin], self.dT[sKey][sMax]
        return applyFilter(self.dDfrIn[sKey], sHdC, thMin, thMax)
    
    def procAndFilter(self, sHdC, sKey, sMin=S_MIN, sMax=S_MAX):
        # process data - concatenate IC DataFrames of the three GT
        concDfr(self.dDfrIn[sKey])
        # filter data
        self.dDfrFl = {sKey: {}}
        for sGT in L_S_GT:
            thMin, thMax = self.dT[sKey][sGT][sMin], self.dT[sKey][sGT][sMax]
            dfrFl = applyFilter(self.dDfrIn[sKey][sGT], sHdC, thMin, thMax)
            self.dDfrFl[sKey][sGT] = dfrFl
    
    def sortAndFiltDfr(self):
        for sGT, cDfr in self.dDfrIn[S_IC_M_P].items():
            sortDfr(cDfr, self.dSort[S_IC_M_P][sGT], S_SRT_BY, S_ASC)
        sortDfr(self.dDfrIn[S_D_GT_M], self.dSort[S_D_GT_M], S_SRT_BY, S_ASC)
        sortDfr(self.dDfrIn[S_D_GT_P], self.dSort[S_D_GT_P], S_SRT_BY, S_ASC)
        dfrFl_M = self.simpleFilter(S_D_GT, S_D_GT_M)
        dfrFl_P = self.simpleFilter(S_D_GT, S_D_GT_P)
        dDfrFl_IC = self.procAndFilter(S_IC, S_IC_M_P)
        dfrFl_IC_GT0 = applyFilter(self.dDfrIn, self.dT, S_IC, S_IC_M_P, S_GT0)
        dfrFl_IC_GT1 = applyFilter(self.dDfrIn, self.dT, S_IC, S_IC_M_P, S_GT1)
        dfrFl_IC_GT5 = applyFilter(self.dDfrIn, self.dT, S_IC, S_IC_M_P, S_GT5)
        self.dDfrFl = {S_IC_M_P: {S_GT0: dfrFl_IC_GT0,
                                  S_GT1: dfrFl_IC_GT1,
                                  S_GT5: dfrFl_IC_GT5},
                       S_D_GT_M: dfrFl_M,
                       S_D_GT_P: dfrFl_P}

    def iniDfrResS(self):
        lHCS = [S_MET, S_PHO] + L_S_IC_GT + [S_D_GT_M, S_D_GT_P]
        dDat = {cHCS: [] for cHCS in lHCS if cHCS not in L_S_IC_GT}
        dDat[S_IC_M_P] = {sGT: [] for sGT in L_S_GT}
        lHCF = (self.dHdCol[S_D_GT_M] + self.dHdCol[S_D_GT_P] +
                self.dHdCol[S_IC_M_P][S_GT0] + self.dHdCol[S_IC_M_P][S_GT1] +
                self.dHdCol[S_IC_M_P][S_GT5])
        self.dfrResS = pd.DataFrame(columns=lHCS)
        self.dfrResF = pd.DataFrame(columns=lHCF)
        n, N = 0, self.dDfrFl[S_D_GT_M].shape[0]*self.dDfrFl[S_D_GT_P].shape[0]
        return dDat, n, N

    def fillDfrResS(self):
        dDat, n, N = self.iniDfrResS()
        for i in self.dDfrFl[S_D_GT_M].index:
            sMet = self.dDfrFl[S_D_GT_M].at[i, S_MET]
            dstGTM = self.dDfrFl[S_D_GT_M].at[i, S_D_GT]
            for j in self.dDfrFl[S_D_GT_P].index:
                sPho = self.dDfrFl[S_D_GT_P].at[j, S_PHO]
                dstGTP = self.dDfrFl[S_D_GT_P].at[i, S_D_GT]
                dDat[S_MET].append(sMet)
                dDat[S_D_GT_M].append(dstGTM)
                dDat[S_PHO].append(sPho)
                dDat[S_D_GT_P].append(dstGTP)
                fillDDat(dDat, self.dDfrFl, sMet, sPho, dstGTM, dstGTP)
                n += 1
                if n%self.inpD.modDisp == 0:
                    print('Processed element', n, 'of', N, '.')
        # print('TEMP - dDat =', dDat)
        print('dDat keys:', list(dDat))
        print('dDat value lengths:')
        for sK, cL in dDat.items():
            print(sK, ':', len(cL))
        dAdd = {sK: cL for sK, cL in dDat.items() if sK not in [S_IC_M_P]}
        for i, sIC_GT in enumerate(L_S_IC_GT):
            dAdd[sIC_GT] = dDat[S_IC_M_P][L_S_GT[i]]
        print('dAdd keys:', list(dAdd))
        print('dAdd value lengths:')
        for sK, cL in dAdd.items():
            print(sK, ':', len(cL))
        dfrAdd = pd.DataFrame(dAdd)
        print(dfrAdd)
        self.dfrResS = self.dfrResS.append(dfrAdd, ignore_index=True,
                                           verify_integrity=True)

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

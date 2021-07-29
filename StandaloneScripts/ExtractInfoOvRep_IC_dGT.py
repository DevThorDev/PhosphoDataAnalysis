# -*- coding: utf-8 -*-
###############################################################################
# --- ExtractInfoOvRep_IC_dGT.py ----------------------------------------------
###############################################################################
import os, time

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# --- CONSTANTS ---------------------------------------------------------------
S_DASH = '-'
S_USC = '_'
S_DOT = '.'
S_CSV = 'csv'
S_PDF = 'pdf'
S_P = 'p'
S_M = 'm'
S_NO = 'No'

S_O = '0'
S_1 = '1'
S_2 = '2'
S_3 = '3'
S_4 = '4'
S_5 = '5'

S_MET = 'Metabolite'
S_PHO = 'Phosphopeptide'
L_S_M_P = [S_MET, S_PHO]

S_BIN_L = 'BinCode'
S_BIN_L_2 = S_BIN_L + S_2
S_SEL = 'Sel'
S_SELECTED = 'Selected'

S_GT0 = 'GT' + S_O
S_GT1 = 'GT' + S_1
S_GT5 = 'GT' + S_5
L_S_GT = [S_GT0, S_GT1, S_GT5]

S_IC_M_P = 'ICMetPho'
S_D_GT_M = 'dGTMet'
S_D_GT_P = 'dGTPho'
L_S_D_GT = [S_D_GT_M, S_D_GT_P]

L_S_FT = ['DR', 'DS', 'NR', 'NS']

D_HD_C_PA = {sMP: {sGT: [S_USC.join([sGT, sFt, sMP[0]]) for sFt in L_S_FT]
                   for sGT in L_S_GT} for sMP in L_S_M_P}

S_MIN = 'min'
S_MAX = 'max'
L_S_MIN_MAX = [S_MIN, S_MAX]

S_SRT_BY = 'SortedBy'
S_ASC = 'Asc'
S_Z_SCORE = 'Pattern (z-score)'
S_IC = 'IC'
S_D_GT = 'd_GT'

S_IC_GT0 = S_USC.join([S_IC, S_GT0])
S_IC_GT1 = S_USC.join([S_IC, S_GT1])
S_IC_GT5 = S_USC.join([S_IC, S_GT5])
L_S_IC_GT = [S_IC_GT0, S_IC_GT1, S_IC_GT5]

S_SG_MET = 'MetSig5'
S_SG_PHO = 'PhoSig5'
L_S_SG = [S_SG_MET, S_SG_PHO]
L_S_SG_GT = [S_USC.join([sSg, sGT]) for sGT in L_S_GT for sSg in L_S_SG]

S_BASE_CL = 'BaseClass'
S_INP_DATA = 'InputData'
S_EXTR_INFO = 'ExtrInfo'
S_PLTR = 'Plotter'
S_PAT_PLTR = 'PatternPlotter'

S_RMNG_COL1_IC = 'PearsonCorr'
S_NEW_IDX = 'NewIndex'
# L_S_NO_GT = L_S_M_P + ['Protein', 'BinCode', 'BinCode2', 'MapMan',
#                        'Description', 'SelBinCode']
L_S_NO_GT = L_S_M_P + []
# L_S_NO_GT = []
L_S_ADD_GT = [S_RMNG_COL1_IC, 'SpearmanCorr', 'Pearson_pVal', 'Spearman_pVal',
              'IC_N', 'IC_P', 'IC', 'MetSig5', 'PhoSig5']

S_NM_PAT_PLT = 'PatternPlot'

R04 = 4

# --- INPUT -------------------------------------------------------------------
# --- flow control ------------------------------------------------------------
doInfoExtr = False               # True / False
doPlotPat = True                # True / False

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
# dThr = {S_IC_M_P: {S_GT0: {S_MIN: 6.0, S_MAX: None},
#                     S_GT1: {S_MIN: 6.0, S_MAX: None},
#                     S_GT5: {S_MIN: 6.0, S_MAX: None}},
#         S_D_GT_M: {S_MIN: 0.6, S_MAX: None},
#         S_D_GT_P: {S_MIN: 0.6, S_MAX: None}}
dThr = {S_IC_M_P: {S_GT0: {S_MIN: 6.0, S_MAX: None},
                   S_GT1: {S_MIN: 6.0, S_MAX: None},
                   S_GT5: {S_MIN: 6.0, S_MAX: None}},
        S_D_GT_M: {S_MIN: None, S_MAX: None},
        S_D_GT_P: {S_MIN: None, S_MAX: None}}
dSel = {S_SG_MET: {S_GT0: ['Y', 'N'], S_GT1: ['Y', 'N'], S_GT5: ['Y', 'N']},
        S_SG_PHO: {S_GT0: ['Y', 'N'], S_GT1: ['Y', 'N'], S_GT5: ['Y', 'N']}}

sSep = ';'

# --- graphics parameters -----------------------------------------------------
# dPairsPaP = {(('Leu_STTTTV'), (S_GT0, S_GT1, S_GT5)):
#              ('Leucine', 'STTTTVS(0.003)S(0.996)VHS(0.001)PTTDQDFSK')}
# dPairsPaP = {(('Tetra_S(0.001)AS(0.749)T(0.251)P'), (S_GT0, S_GT1, S_GT5)):
#               ('Tetradecanoic_acid', 'S(0.001)AS(0.749)T(0.251)PLLNSLVHVS(0.179)S(0.821)PRDS(1)PIETVESVHQIQR'),
#               (('Beta_ADKTDII'), (S_GT0, S_GT1, S_GT5)):
#               ('Beta-alanine', 'ADKTDIIS(0.607)S(0.117)S(0.12)S(0.156)DKAS(1)PPPPSAFR'),
#               (('Hexa_S(0.001)AS(0.749)T(0.251)P'), (S_GT0, S_GT1, S_GT5)):
#               ('Hexadecanoic_acid', 'S(0.001)AS(0.749)T(0.251)PLLNSLVHVS(0.179)S(0.821)PRDS(1)PIETVESVHQIQR')}
dPairsPaP = {(('Isoleu_DLDVNE'), (S_GT0, S_GT1, S_GT5)):
              ('Isoleucine', 'DLDVNES(1)GPPAAR'),
              (('Val_DLDVNE'), (S_GT0, S_GT1, S_GT5)):
              ('Valine', 'DLDVNES(1)GPPAAR'),
              (('Aspart_SDKPLNY'), (S_GT0, S_GT1, S_GT5)):
              ('Aspartic_acid', 'SDKPLNYS(1)PDPENESGINER'),
              (('Aspart_DLDVNE'), (S_GT0, S_GT1, S_GT5)):
              ('Aspartic_acid', 'DLDVNES(1)GPPAAR'),
              (('Malic_DLDVNE'), (S_GT0, S_GT1, S_GT5)):
              ('Malic_acid', 'DLDVNES(1)GPPAAR')}

nmPaP = S_NM_PAT_PLT            # name prefix of the plot
szFontLeg = 'small'             # font size of legend
nCharDsp = 60                   # number of chars displayed for legend item
coordAnchorBox = (0.5, 1.02)    # coordinates of the legend anchor box
lWdPlt = 0.75                   # line width in plot

# --- names and paths of files and dirs ---------------------------------------
sFIn_IC_M_P = 'IC_Met_Pho'
sFIn_dGT_M = 'DistGT_Met'
sFIn_dGT_P = 'DistGT_Pho'

sFOutS = 'ExtrIOvRepS'
sFOutF = 'ExtrIOvRepF'

# sFIn_PaP = 'ExtrIOvRepF_ICMetPho_GT0_0_No_No_GT1_0_No_No_GT5_0_No_No_dGTMet_0_0p6_No_dGTPho_0_0p6_No'
sFIn_PaP = 'ExtrIOvRepF_ICMetPho_GT0_0_6p0_No_GT1_0_6p0_No_GT5_0_6p0_No_dGTMet_0_No_No_dGTPho_0_No_No'
sFOutPaP = S_NM_PAT_PLT

sDirInCSV = '51_CSV_In_DistGT'
sDirOutCSV = '52_CSV_Out_DistGT'
sDirOutPaP = '56_PatternPlots'

pBase = os.path.join('..', '..', '..', '25_Papers', '01_FirstAuthor',
                    '04_SysBio_DataAnalysis')
pInCSV = os.path.join(pBase, sDirInCSV)
pOutCSV = os.path.join(pBase, sDirOutCSV)

pInPaP = os.path.join(pBase, sDirOutCSV)
pOutPaP = os.path.join(pBase, sDirOutPaP)

# --- derived values ----------------------------------------------------------
pFIGT0 = os.path.join(pInCSV, sFIn_IC_M_P + S_USC + S_GT0 + S_DOT + S_CSV)
pFIGT1 = os.path.join(pInCSV, sFIn_IC_M_P + S_USC + S_GT1 + S_DOT + S_CSV)
pFIGT5 = os.path.join(pInCSV, sFIn_IC_M_P + S_USC + S_GT5 + S_DOT + S_CSV)
dPFInIC = {S_GT0: pFIGT0, S_GT1: pFIGT1, S_GT5: pFIGT5}

# --- assertions --------------------------------------------------------------
assert set(dISort) == set(dThr)
for cD in [dISort, dThr]:
    assert S_IC_M_P in cD and set(cD[S_IC_M_P]) == set(L_S_GT)

# --- INPUT DICTIONARY --------------------------------------------------------
dInput = {# --- constants
          'sBC_L': S_BIN_L,
          'sBC2_L': S_BIN_L_2,
          'lSGT': L_S_GT,
          'sBase': S_BASE_CL,
          'sInpDat': S_INP_DATA,
          'sExtrInfo': S_EXTR_INFO,
          'sPltr': S_PLTR,
          'sPatPltr': S_PAT_PLTR,
          'sMet': S_MET,
          'sPho': S_PHO,
          'R04': R04,
          # --- flow control
          'doInfoExtr': doInfoExtr,
          'doPlotPat': doPlotPat,
          # --- general input
          'modDisp': modDisp,
          # --- data specific input
          'dISort': dISort,
          'dThr': dThr,
          'dSel': dSel,
          'sSep': sSep,
          # --- graphics parameters
          'plotOfPatterns': {'dPairsPaP': dPairsPaP,
                             'nmPaP': nmPaP,
                             'szFontLeg': szFontLeg,
                             'nCharDsp': nCharDsp,
                             'coordAnchorBox': coordAnchorBox,
                             'lWdPlt': lWdPlt},
          # --- names and paths of files and dirs
          'pInCSV': pInCSV,
          'pOutCSV': pOutCSV,
          'pInPaP': pInPaP,
          'pOutPaP': pOutPaP,
          'dPFInIC': dPFInIC,
          'pFInM': os.path.join(pInCSV, sFIn_dGT_M + S_DOT + S_CSV),
          'pFInP': os.path.join(pInCSV, sFIn_dGT_P + S_DOT + S_CSV),
          'sFOutS': sFOutS,
          'sFOutF': sFOutF,
          'pFInPaP': os.path.join(pInPaP, sFIn_PaP + S_DOT + S_CSV),
          'sFOutPaP': sFOutPaP}

# --- FUNCTIONS ---------------------------------------------------------------
def addToDictD(cD, cKMain, cKSub, cV):
    if cKMain in cD:
        assert cKSub not in cD[cKMain]
        cD[cKMain][cKSub] = cV
    else:
        cD[cKMain] = {}
        cD[cKMain][cKSub] = cV

def procNum4NmF(cV, dRepl={S_DOT: S_P, S_DASH: S_M, 'None': S_NO}):
    cS = str(cV)
    for sK, sV in dRepl.items():
        cS = cS.replace(sK, sV)
    return cS

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

def modifyDfr(pdDfr, sGT):
    lHdCN = [s + S_USC + sGT for s in pdDfr.columns if s not in L_S_NO_GT]
    pdDfr.columns = L_S_NO_GT + lHdCN
    dfrMod = pd.concat([pd.Series(getLNewIdx(pdDfr), name=S_NEW_IDX), pdDfr],
                       axis=1, verify_integrity=True)
    return dfrMod.set_index(S_NEW_IDX, verify_integrity=True)

def concDfr(dDfr):
    dfrM = modifyDfr(dDfr[L_S_GT[0]], L_S_GT[0])
    for sGT in L_S_GT[1:]:
        dfrMR = modifyDfr(dDfr[sGT], sGT).drop(columns=L_S_NO_GT)
        dfrM = pd.concat([dfrM, dfrMR], axis=1, verify_integrity=True)
    return dfrM.drop(columns=L_S_NO_GT)

def appendToDDat(lDat, cDfrFl, sI, sHdC):
    cV = np.nan
    if sI in cDfrFl.index:
        cV = cDfrFl.at[sI, sHdC]
    lDat.append(cV)

def decorateClosePlot(cFig, cAx, dPlt, pPltF):
    cAx.set_ylabel(S_Z_SCORE)
    l = cAx.legend(loc='lower center', bbox_to_anchor=dPlt['coordAnchorBox'],
                   fontsize=dPlt['szFontLeg'])
    if l is not None:
        cFig.savefig(pPltF, bbox_extra_artists=(l,), bbox_inches='tight')
    else:
        cFig.savefig(pPltF)
    plt.close()

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
        self.sSp = self.inpD.sSep
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
        self.dSl = self.inpD.dSel
        sFOutS, sFOutF, sEnd = self.inpD.sFOutS, self.inpD.sFOutF, ''
        for sK in self.dSort:
            if sK in L_S_D_GT:
                sEnd += S_USC + sK + S_USC + str(int(self.dSort[sK][S_ASC]))
                for sMM in L_S_MIN_MAX:
                    sEnd += S_USC + procNum4NmF(self.dT[sK][sMM])
            else:
                for i, sGT in enumerate(self.dSort[sK]):
                    sAscDsc = str(int(self.dSort[sK][sGT][S_ASC]))
                    if i == 0:
                        sEnd += S_USC + sK + S_USC + sGT + S_USC + sAscDsc
                    else:
                        sEnd += S_USC + sGT + S_USC + sAscDsc
                    for sMM in L_S_MIN_MAX:
                        sEnd += S_USC + procNum4NmF(self.dT[sK][sGT][sMM])
        sFOutS, sFOutF = sFOutS + sEnd, sFOutF + sEnd
        self.pFOutS = os.path.join(self.inpD.pOutCSV, sFOutS + S_DOT + S_CSV)
        self.pFOutF = os.path.join(self.inpD.pOutCSV, sFOutF + S_DOT + S_CSV)

    def getInf4Inp(self):
        dDatTp_IC = {sIn: str for sIn in [self.inpD.sBC_L, self.inpD.sBC2_L]}
        dDatTp_P = {sIn: str for sIn in [self.inpD.sBC2_L]}
        self.dDatTp = {S_IC_M_P: dDatTp_IC, S_D_GT_P: dDatTp_P}
        self.dPFIn = {S_IC_M_P: self.inpD.dPFInIC,
                      S_D_GT_M: self.inpD.pFInM,
                      S_D_GT_P: self.inpD.pFInP}

    def loadDfrInp(self):
        # load input DataFrames
        dDfrIn_IC = {sGT: pd.read_csv(self.dPFIn[S_IC_M_P][sGT], sep=self.sSp,
                                      dtype=self.dDatTp[S_IC_M_P])
                     for sGT in L_S_GT}
        dfrIn_M = pd.read_csv(self.dPFIn[S_D_GT_M], sep=self.sSp)
        dfrIn_P = pd.read_csv(self.dPFIn[S_D_GT_P], sep=self.sSp,
                              dtype=self.dDatTp[S_D_GT_P])
        self.dDfrIn = {S_IC_M_P: dDfrIn_IC,
                       S_D_GT_M: dfrIn_M,
                       S_D_GT_P: dfrIn_P}

    def getDHdCol(self):
        # get dictionary of column headers
        dDfr_IC = self.dDfrIn[S_IC_M_P]
        dHdCol_IC = {sGT: [sC + S_USC + sGT for sC in
                           list(dDfr_IC[sGT].loc[:, S_RMNG_COL1_IC:].columns)]
                     for sGT in L_S_GT}
        self.dHdCol = {S_IC_M_P: dHdCol_IC,
                       S_D_GT_M: list(self.dDfrIn[S_D_GT_M].columns),
                       S_D_GT_P: list(self.dDfrIn[S_D_GT_P].columns)}

    def getDMapK(self):
        # get mapping of data dictionary keys to column headers
        self.dMapK = {}
        for i, sMP in enumerate(L_S_M_P):
            self.dMapK[sMP] = (L_S_D_GT[i], sMP)
            lHdCRed = [s for s in self.dHdCol[L_S_D_GT[i]] if s not in L_S_M_P]
            for sHdC in lHdCRed:
                self.dMapK[S_USC.join([sHdC, sMP[0]])] = (L_S_D_GT[i], sHdC)
        for sGT in L_S_GT:
            for s in self.dHdCol[S_IC_M_P][sGT]:
                self.dMapK[s] = (S_IC_M_P, s)

    def getProcData(self):
        self.getPResF()
        self.getInf4Inp()
        self.loadDfrInp()
        self.getDHdCol()
        self.getDMapK()

    def simpleFilter(self, sHdC, sKey, sMin=S_MIN, sMax=S_MAX):
        thMin, thMax = self.dT[sKey][sMin], self.dT[sKey][sMax]
        return applyFilter(self.dDfrIn[sKey], sHdC, thMin, thMax)

    def filterAndConc(self, sHdC, sKey, sMin=S_MIN, sMax=S_MAX):
        # filter data
        self.dDfrFl = {sKey: {}}
        for sGT in L_S_GT:
            thMin, thMax = self.dT[sKey][sGT][sMin], self.dT[sKey][sGT][sMax]
            dfrFl = applyFilter(self.dDfrIn[sKey][sGT], sHdC, thMin, thMax)
            self.dDfrFl[sKey][sGT] = dfrFl
        # process data - concatenate IC DataFrames of the three GT
        return concDfr(self.dDfrFl[sKey])

    def sortAndFiltDfr(self):
        for sGT, cDfr in self.dDfrIn[S_IC_M_P].items():
            sortDfr(cDfr, self.dSort[S_IC_M_P][sGT], S_SRT_BY, S_ASC)
        sortDfr(self.dDfrIn[S_D_GT_M], self.dSort[S_D_GT_M], S_SRT_BY, S_ASC)
        sortDfr(self.dDfrIn[S_D_GT_P], self.dSort[S_D_GT_P], S_SRT_BY, S_ASC)
        self.dDfrFl = {S_IC_M_P: self.filterAndConc(S_IC, S_IC_M_P),
                       S_D_GT_M: self.simpleFilter(S_D_GT, S_D_GT_M),
                       S_D_GT_P: self.simpleFilter(S_D_GT, S_D_GT_P)}
        self.N = self.dDfrFl[S_D_GT_M].shape[0]*self.dDfrFl[S_D_GT_P].shape[0]

    def iniDfrRes(self, specSel=None):
        lC = list(self.dMapK)
        if specSel == 'S':  # the "short" subset of data
            # select only a subset of all keys (columns of DataFrames): lC
            lC = [S_BIN_L_2, S_SELECTED] + L_S_IC_GT + L_S_D_GT
            lC = [S_USC.join([s, S_PHO[0]]) for s in [S_BIN_L_2, S_SELECTED]]
            lC = L_S_M_P + lC + L_S_IC_GT
            lC += [S_USC.join([S_D_GT, s[0]]) for s in L_S_M_P]
            self.dfrResS = pd.DataFrame(columns=lC)
        else:
            # select all keys (columns of DataFrames) defined in self.dMapK
            self.dfrResF = pd.DataFrame(columns=lC)
        return {cK: [] for cK in lC}

    def SUB_fillDDat(self, dDat, i, j, sIMP):
        for sKDDt, lDt in dDat.items():
            assert sKDDt in self.dMapK
            (sKDFl, sHdC) = self.dMapK[sKDDt]
            if sKDFl in L_S_D_GT:
                assert sKDFl in self.dDfrFl
                assert sHdC in self.dDfrFl[sKDFl].columns
                if sKDFl == S_D_GT_M:
                    dDat[sKDDt].append(self.dDfrFl[sKDFl].at[i, sHdC])
                else:
                    dDat[sKDDt].append(self.dDfrFl[sKDFl].at[j, sHdC])
            elif sKDFl in S_IC_M_P:
                appendToDDat(lDt, self.dDfrFl[S_IC_M_P], sIMP, sHdC)
            else:
                print('ERROR: Unknown key of filter dictionary:', sKDFl)
                assert False

    def fillDDat(self, dDat):
        n = 0
        for i in self.dDfrFl[S_D_GT_M].index:
            sMet = self.dDfrFl[S_D_GT_M].at[i, S_MET]
            for j in self.dDfrFl[S_D_GT_P].index:
                sPho = self.dDfrFl[S_D_GT_P].at[j, S_PHO]
                sIMP = S_USC.join([sMet, sPho])
                self.SUB_fillDDat(dDat, i, j, sIMP)
                n += 1
                if n%self.inpD.modDisp == 0:
                    print('Processed element', n, 'of', self.N, '.')

    def fillPrintDfrRes(self):
        for spcSel in ['S', 'F']:
            dDat = self.iniDfrRes(specSel=spcSel)
            self.fillDDat(dDat)
            print('Filled data dictionary for selection "' + spcSel + '".')
            if spcSel == 'S':
                self.dfrResS = self.dfrResS.append(pd.DataFrame(dDat),
                                                   ignore_index=True,
                                                   verify_integrity=True)
                self.dfrResS.to_csv(self.pFOutS, sep=self.sSp)
            else:
                self.dfrResF = self.dfrResF.append(pd.DataFrame(dDat),
                                                   ignore_index=True,
                                                   verify_integrity=True)
                self.dfrResF.to_csv(self.pFOutF, sep=self.sSp)

    def extractionOfExtremes(self):
        self.sortAndFiltDfr()
        self.fillPrintDfrRes()

class Plotter(BaseClass):
    def __init__(self, InpD):
        super().__init__()
        self.idO = InpD.sPltr
        self.descO = 'Class for plotting'
        self.inpD = InpD
        self.sSp = self.inpD.sSep
        self.dPPltF = {}
        print('Initiated "Plotter" base object.')

    def printDPPltF(self):
        print('Dictionary of plot file paths:')
        for tK, pF in self.dPPltF.items():
            print(tK, ':', pF)
        print('-'*64)

    def loadDfrInp(self, iC=0):
        self.dfrIn, dDatTp = None, {sIn: str for sIn in L_S_SG_GT}
        # load input DataFrames
        if hasattr(self, 'pFIn'):
            self.dfrIn = pd.read_csv(self.pFIn, sep=self.sSp, index_col=iC,
                                     dtype=dDatTp)

class PatternPlotter(Plotter):
    def __init__(self, InpD):
        super().__init__(InpD)
        self.idO = self.inpD.sPatPltr
        self.descO = 'Pattern plotter'
        self.pFIn = self.inpD.pFInPaP
        self.pDOut = self.inpD.pOutPaP
        self.dPlt = self.inpD.plotOfPatterns
        self.getDPPltF()
        self.loadDfrInp()
        print('Initiated "PatternPlotter" base object and loaded input data.')

    def getDPPltF(self):
        self.dPPltF, sFPlt = {}, self.inpD.sFOutPaP
        for ((s1, tSGT), tMP) in self.dPlt['dPairsPaP'].items():
            for sGT in tSGT:
                sPltF = S_DOT.join([S_USC.join([sFPlt, s1, sGT]), S_PDF])
                self.dPPltF[(s1, sGT)] = (tMP, os.path.join(self.pDOut, sPltF))

    def plotPatterns(self):
        d, nChD = self.dfrIn, self.dPlt['nCharDsp']
        for ((s1, sGT), ((sM, sP), pPltF)) in self.dPPltF.items():
            print('Plotting pattern for "' + s1 + '" and', sGT, '...')
            # if not os.path.isfile(pPltF):
            cFig, cAx = plt.subplots()
            cSer = d[(d[S_MET] == sM) & (d[S_PHO] == sP)].squeeze()
            for tMP in [(S_MET, sM), (S_PHO, sP)]:
                cPa = cSer.loc[D_HD_C_PA[tMP[0]][sGT]]
                cPa.index = L_S_FT
                cAx.plot(cPa, lw=self.dPlt['lWdPlt'], label=tMP[1][:nChD])
            decorateClosePlot(cFig, cAx, self.dPlt, pPltF)

# --- MAIN --------------------------------------------------------------------
startTime = time.time()
print('+'*25 + ' START', time.ctime(startTime), '+'*25)

inpDat = InputData(dInput)
if inpDat.doInfoExtr:
    cXtrInfo = ExtractedInfo(inpDat)
    # cXtrInfo.printObjInfo()
    cXtrInfo.extractionOfExtremes()
if inpDat.doPlotPat:
    cPltr = PatternPlotter(inpDat)
    # cPltr.printAttrData()
    # cPltr.printDPPltF()
    cPltr.plotPatterns()

print('-'*80)
printElapsedTimeSim(startTime, time.time(), 'Total time')
print('+'*37, 'DONE.', '+'*36)

###############################################################################

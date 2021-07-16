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

S_GT = 'GT'
S_GT0 = S_GT + S_O
S_GT1 = S_GT + S_1
S_GT5 = S_GT + S_5
S_ALL_GT = 'All' + S_GT

S_MC = 'MC'
S_IC = 'IC'
S_IDIST = 'IDist'

S_MC0 = S_MC + S_O
S_MC1 = S_MC + S_1
L_S_MC = [S_MC0, S_MC1]

S_IC0 = S_IC + S_O
S_IC1 = S_IC + S_1
L_S_IC = [S_IC0, S_IC1]

S_IDIST0 = S_IDIST + S_O
S_IDIST1 = S_IDIST + S_1
L_S_IDIST = [S_IDIST0, S_IDIST1]

L_S_CTRV = ['Sign', 'WtsA', 'WtsB']
L_S_GT = [S_GT0, S_GT1, S_GT5]

L_S_HDR_MC = [s1 + S_USC + s2 + S_USC + s3 for s1 in L_S_MC
              for s2 in L_S_CTRV for s3 in L_S_GT]
L_S_HDR_IC = [s1 + S_USC + s2 + S_USC + s3 for s1 in L_S_IC
              for s2 in L_S_CTRV for s3 in L_S_GT]
L_S_HDR_IDIST = [s1 + S_USC + s2 + S_USC + S_ALL_GT for s1 in L_S_IDIST
                 for s2 in L_S_CTRV]

S_GENERAL = 'General'
S_I_ORIG = 'IdxOrig'
S_RI = 'RI'
S_SEL = 'Selections'
S_MET_S = 'Met'
S_MET_L = 'Metabolite'
S_MET_D = 'MetD'
S_PHO_D = 'PhoD'
S_BIN_OP = 'BinOp'
S_BIN_CODES_L = S_BIN_CODE_L + 's'

S_IC_MET = S_IC + S_USC + S_MET_S
S_IC_BIN_CODE_S_2 = S_IC + S_USC + S_BIN_CODE_S_2

R04 = 4

# --- INPUT -------------------------------------------------------------------
lTpDat = [S_MC, S_IC_MET, S_IC_BIN_CODE_S_2, S_IDIST]

sFInp_IC_Met_Pho = 'IC_Met_Pho'
sFInp_dGT_Met = 'DistGT_Met'
sFInp_dGT_Pho = 'DistGT_Pho'

sFOut = 'ExtractedInfoOvRep_IC_dGT'

pCSV = os.path.join('..', '..', '..', '25_Papers', '01_FirstAuthor',
                    '04_SysBio_DataAnalysis', '51_CSV_DistGT')
sSep = ';'

# --- INPUT DICTIONARY --------------------------------------------------------
dInput = {S_GENERAL: {'lTpDat': lTpDat,
                      'pCSV': pCSV,
                      'sSep': sSep},
          S_MC: {'sFInp': sFInp_MC + S_DOT + S_CSV,
                 'sFOut': sFOut_MC + S_DOT + S_CSV,
                 'lSHdr': L_S_HDR_MC,
                 'sHdrRef': S_BIN_CODE_L_2},
          S_IC_MET: {'sFInp': sFInp_ICMet + S_DOT + S_CSV,
                     'sFOut': sFOut_ICMet + S_DOT + S_CSV,
                     'lSHdr': L_S_HDR_IC,
                     'sHdrRef': S_MET_L},
          S_IC_BIN_CODE_S_2: {'sFInp': sFInp_ICBC2 + S_DOT + S_CSV,
                              'sFOut': sFOut_ICBC2 + S_DOT + S_CSV,
                              'lSHdr': L_S_HDR_IC,
                              'sHdrRef': S_BIN_CODE_L_2},
          S_IDIST: {'sFInp': sFInp_IDist + S_DOT + S_CSV,
                    'sFOut': sFOut_IDist + S_DOT + S_CSV,
                    'lSHdr': L_S_HDR_IDIST,
                    'sHdrRef': S_BIN_CODE_L_2}}
for cTpDat in dInput[S_GENERAL]['lTpDat']:
    dInput[cTpDat]['pFInp'] = os.path.join(pCSV, dInput[cTpDat]['sFInp'])
    dInput[cTpDat]['pFOut'] = os.path.join(pCSV, dInput[cTpDat]['sFOut'])

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
        self.idO = InpD.sOvRep
        self.descO = 'Over-representation'
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

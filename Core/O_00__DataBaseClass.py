# -*- coding: utf-8 -*-
###############################################################################
# --- O_00__DataBaseClass.py --------------------------------------------------
###############################################################################
import pprint

import Core.C_00__GenConstants as GC
import Core.F_00__GenFunctions as GF
import Core.F_03__OTpFunctions as TF

class DataBaseClass:
    def __init__(self, inpDat, iTp = 0, cGT = '', cFt = '', nmFAdd = ''):
        self.idO = GC.S_BASE
        self.descO = 'Data base class'
        self.dIG = inpDat.dI
        self.dITp = TF.getDITp(self.dIG, 0, iTp)
        self.cM = self.dIG['Mode']
        self.sGTFt = self.setGenTpFeat(cGT, cFt)
        self.sTrDt, self.sDvDt = self.dIG['nmPreTr'], self.dIG['nmPreDv']
        self.idOC, self.idGTC = self.dIG['idOAtDfr'], self.dIG['idGTAtDfr']
        self.nmGTC, self.idFtC = self.dIG['nmGTAtDfr'], self.dIG['idFtAtDfr']
        self.nmFtC, self.nmWtC = self.dIG['nmFtAtDfr'], self.dIG['nmWtAtDfr']
        self.lNmAttr = self.dIG['lNmAtDfr']
        lAttrIni = TF.getLAttr(self.dIG, [self.idO, self.cGT, self.cFt, 1.],
                               iCGT = 1, iCFt = 2)
        self.attrDfr = GF.iniAttrDfr(lAttrIni, self.lNmAttr)
        self.attrDfrMn = GF.iniAttrDfr(lAttrIni, self.lNmAttr)
        self.setDummys()
        self.setNmFE(nmFAdd)
        self.getDOIn()
        self.updateSFName(setIt = True)
        print('Initiated "DataBaseClass" base object.')

    def __str__(self):
        sIn = ('~'*24 + ' ' + self.descO + ' with ID ' + str(self.idO) + ' ' +
               '~'*24 + '\nMode: ' + str(self.cM) + '\n' + '~'*80 + '\n')
        return sIn
    
    def printIDDesc(self):
        print('Object ID:', self.idO)
        print('Tuple of object IDs:', self.tID)
        print('Object description:', self.descO)
    
    def printFeatureHdr(self):
        print('Object ID:', self.idO)
        print('Feature header list (all data):', self.lFtHd)
        print('Feature header list (means data):', self.lFtHdMn)
    
    def printGT(self):
        print('Genotype ID:', self.cGT)
        print('Genotype name:', TF.getNmGT(self.dIG, self.cGT))
        print('Tuple of genotype IDs:', self.tGT)
    
    def printFt(self):
        print('Feature ID:', self.cFt)
        print('Feature name:', TF.getNmFt(self.dIG, self.cFt))
        print('Tuple of feature IDs:', self.tFt)
    
    def printDType(self):
        print('-'*20, 'Type dictionary:', '-'*20)
        pprint.pprint(self.dITp)
        
    def printDObjInst(self):
        print('-'*20, 'Object instance dictionary:', '-'*20)
        pprint.pprint(self.dOIn)
        
    def printObjInfo(self):
        print('-'*20, 'Object', self.descO, '(ID', self.idO, ')', '-'*20)
        if 'iTp' in self.dITp:
            print('Type index:', self.dITp['iTp'])
        print('Mode:', self.cM)
        self.printFeatureHdr()
        print('Tuple of object IDs:', self.tID)
        self.printGT()
        self.printFt()
        print('Data transformation:', self.sTrDt)
        print('List of row names of transformed data:', self.lSRowTrD)
        print('Data deviations type:', self.sDvDt)
        print('List of row names of deviations data:', self.lSRowDvD)
        print('Predefined file name end:', self.nmFE)

    def printDfr(self, printMnSD = None, doT = False):
        print('\n---', self.descO, '(', self.idO, ') - mode', self.cM, ':')
        cDfr, sDfr = self.cDfr, 'full data'
        if printMnSD == self.dIG['nmMeans']:
            cDfr, sDfr = self.cDfrMn, 'means data'
        if printMnSD == self.dIG['nmSDs']:
            cDfr, sDfr = self.cDfrSD, 'SDs data'
        print('Data frame of object', self.idO, '(' + sDfr + '):')
        GF.printDfrOrDfrT(cDfr, doT = doT)
    
    def printDfrN(self, doT = False):
        if self.dfrN.shape == self.cDfrMn.shape:
            print('Data frame containing the number of samples:')
            GF.printDfrOrDfrT(self.dfrN, doT = doT)
        else:
            print('ERROR: Shape of number of samples DataFrame (',
                  self.dfrN.shape, ') differs from shape of means DataFrame (',
                  self.cDfrMn.shape, ').')
    
    def printAddIDfr(self, doT = False):
        print('\n---', self.descO, '(', self.idO, ') - mode', self.cM, ':')
        print('Additional info DataFrame:')
        GF.printDfrOrDfrT(self.addIDfr, doT = doT)
        
    def printDfrAIC(self, doT = False):
        print('\n---', self.descO, '(', self.idO, ') - mode', self.cM, ':')
        print('DataFrame containing classification info:')
        GF.printDfrOrDfrT(GF.getDfrSelAI(self.addIDfr, self.lCHdrAIC),
                          doT = doT)
        
    def printDfrAIR(self, doT = False):
        print('\n---', self.descO, '(', self.idO, ') - mode', self.cM, ':')
        print('DataFrame containing retained additional info:')
        GF.printDfrOrDfrT(GF.getDfrSelAI(self.addIDfr, self.lCHdrAIR),
                          doT = doT)
        
    def printExtDfrAIC(self, doT = False):
        print('\n---', self.descO, '(', self.idO, ') - mode', self.cM, ':')
        print('Current DataFrame including classification info:')
        GF.printDfrOrDfrT(GF.getExtDfrSelAI(self.addIDfr, self.cDfr,
                                            self.lCHdrAIC), doT = doT)
        
    def printExtDfrAIR(self, doT = False):
        print('\n---', self.descO, '(', self.idO, ') - mode', self.cM, ':')
        print('Current DataFrame including retained additional info:')
        GF.printDfrOrDfrT(GF.getExtDfrSelAI(self.addIDfr, self.cDfr,
                                            self.lCHdrAIR), doT = doT)
        
    def printAttrDfr(self, isMn = False, doT = False):
        print('\n---', self.descO, '(', self.idO, ') - mode', self.cM, ':')
        cDfr, sDfr = self.attrDfr, 'Attribute DataFrame (full data):'
        if isMn:
            cDfr, sDfr = self.attrDfrMn, 'Attribute DataFrame (means data):'
        print(sDfr)
        GF.printDfrOrDfrT(cDfr, doT = doT)
        print('Data types:', cDfr.dtypes)
        
    def printFNmComp(self):
        print('File name components dictionary:')
        pprint.pprint(self.dFNmComp)
        print('File name components DataFrame:')
        print(self.dfrFNmComp)
        print('File name string from its components:')
        print(self.sFNm)
        print('File name string (transformed data) from its components:')
        print(self.sFNmTr)
        print('File name string (deviations data) from its components:')
        print(self.sFNmDv)
        
    def printLDType(self):
        print('-'*20, 'List of type dictionaries:', '-'*20)
        for dITp in self.lDITp:
            pprint.pprint(dITp)

    def printOData(self):
        print('---', self.descO, '(', self.idO, ')  - mode', self.cM, ':')
        print('Type index of object:', self.cOITp)
        print('Object data:')
        self.OD.printObjInfo()
        for cArg in [None, self.dIG['nmMeans'], self.dIG['nmSDs']]:
            self.OD.printDfr(cArg)
        print('-'*80)
        
    def printLOData(self):
        print('---', self.descO, '(', self.idO, ')  - mode', self.cM, ':')
        print('List of type indices of data:', self.lITp)
        print('List of object data:')
        for i, cO in enumerate(self.lOD):
            print('\n', '-'*20, 'Object', i + 1, ':', '-'*20)
            cO.printObjInfo()
            for cArg in [None, self.dIG['nmMeans'], self.dIG['nmSDs']]:
                cO.printDfr(cArg)
        print('-'*80)

    def printTrDfr(self):
        print('Transformed DataFrame of object data:')
        pprint.pprint(self.cDfrTr)
        print('Transposed transformed DataFrame of object data:')
        pprint.pprint(self.cDfrTrT)

    def printDvDfr(self):
        print('Deviations DataFrame of object data:')
        pprint.pprint(self.cDfrDv)
        print('Transposed deviations DataFrame of object data:')
        pprint.pprint(self.cDfrDvT)

    def printClDfr(self, is4Corr = False, printSpec = None):
        dClDfr = {}
        if is4Corr:
            dClDfr = self.dClAlgDfr4Corr
        else:
            dClDfr = self.dClAlgDfr
        if printSpec == 'First' and len(dClDfr) > 0:
            print('-'*8, 'First', self.descO, 'DataFrame:','-'*8)
            print(dClDfr[min(dClDfr)])
        elif printSpec == 'Last' and len(dClDfr) > 0:
            print('-'*8, 'Last', self.descO, 'DataFrame:','-'*8)
            print(dClDfr[max(dClDfr)])
        else:
            print('+'*16, 'List of', self.descO, 'DataFrames:','+'*16)
            for cDfr in dClDfr:
                print('-'*8, 'Current', self.descO, 'DataFrame:','-'*8)
                print(cDfr)

    def printDClusterResData(self, lNCl = []):
        lNClOut = lNCl
        if len(lNCl) == 0:
            lNClOut = list(self.dClResData)
        for cK in lNClOut:
            if cK in self.dClResData:
                self.dClResData[cK].printObjInfo()
                # self.dClResData[cK].printDType()
                self.dClResData[cK].printDfr()
    
    def setGenTpFeat(self, GT, Ft, sSep = '_'):
        s = ''
        self.cGT = GT
        if self.cGT == '':
            self.cGT = self.dIG['nmAllGT']
        else:
            s = self.cGT
        self.cFt = Ft
        if self.cFt == '':
            self.cFt = self.dIG['nmAllFt']
            if self.cGT == self.dIG['nmAllGT']:
                s = self.dIG['nmAllGTFt']
        else:
            if self.cGT == self.dIG['nmAllGT']:
                s = self.cFt
            else:
                s += (sSep + self.cFt)
        return s
    
    def setDummys(self):
        self.iCNm = -1
        self.iCVS = -1
        self.lIC2Str, self.lCHdrAIC, self.lCHdrAIR = [], [], []
        self.cQry = ''
        self.tIMM = (-1, -1)
        self.lRIDp = []
        self.lCIDp = []
        self.tpX = None
        self.cXFt = 'D'
        self.XtrDt = False
        self.cDfr = GF.iniPdDfr()
        self.cDfrMn = GF.iniPdDfr()
        self.cDfrSD = GF.iniPdDfr()
        self.dfrN = GF.iniPdDfr()
        self.addIDfr = GF.iniPdDfr()
        self.lFtHd, self.lFtHdMn = [], []
        self.lOD = []
        self.lDITp = []
        self.lITp = []
        self.lODfr, self.lODfrMn = [], []
        self.tID = ('',)
        self.tGT = (self.cGT,)
        self.tFt = (self.cFt,)
        self.lICNm = []
        self.lICVS = []
        self.lGTSglO = [self.cGT]
        self.lFtSglO = [self.cFt]
        self.lODfrAI = []
        self.lODOIn = []
        self.pFAllD = ''
        self.pFMn = ''
        self.pFTrD = ''
        self.pFDvD = ''
        self.pFCrDvAllD = ''
        self.pFCrDvMn = ''
        self.lPFTACDAllD = ''
        self.lPFTACDMn = ''
        self.OD = None
        self.dOITp = {}
        self.cOITp = -1
        self.dClAlgDfr, self.dClAlgDfr4Corr = {}, {}
        self.dClResData = {}
        self.lSRowTrD,self.lSRowDvD = [], []
        self.dfrRCrDv, self.dfrRCrDvMn = GF.iniPdDfr(), GF.iniPdDfr()
        self.cDfrTr, self.cDfrTrT = GF.iniPdDfr(), GF.iniPdDfr()
        self.cDfrDv, self.cDfrDvT = GF.iniPdDfr(), GF.iniPdDfr()
        self.pRF, self.pRFF = '', ''
        self.nClst = 0
        self.isClRD = False
        
    def setNmFE(self, nmFAdd):
        self.nmFE = nmFAdd
        if self.nmFE == '':
            self.nmFE = self.dIG['nmEStd']
    
    def printSumStat(self):
        print('\n--- Summary statistics of', self.descO,
              '(', self.idO,') - mode', self.cM, ':')
        tSumStat = GF.calcSumStat(self.cDfr)
        assert len(tSumStat) >= 6
        print('Means:\n', list(tSumStat[0]), '\n')
        print('Medians:\n', list(tSumStat[1]), '\n')
        print('Mins:\n', list(tSumStat[2]), '\n')
        print('Maxs:\n', list(tSumStat[3]), '\n')
        print('Standard devs:\n', list(tSumStat[4]), '\n')
        print('Variances:\n', list(tSumStat[5]), '\n')
        
    def getDOIn(self):
        self.dOIn = {}
        lA = ['idO', 'descO', 'cM', 'cGT', 'cFt', 'sGTFt', 'sTrDt', 'sDvDt',
              'idOC', 'idGTC', 'nmGTC', 'idFtC', 'nmFtC', 'nmWtC', 'lNmAttr',
              'attrDfr', 'attrDfrMn', 'nmFE', 'dfrFNmComp', 'dFNmComp', 'sFNm',
              'sFNmTr', 'sFNmDv', 'iCNm', 'iCVS', 'lIC2Str', 'lCHdrAIC',
              'lCHdrAIR', 'cQry', 'tIMM', 'lRIDp', 'lCIDp', 'tpX', 'cXFt',
              'XtrDt', 'cDfr', 'cDfrMn', 'cDfrSD', 'dfrN', 'addIDfr', 'lFtHd',
              'lFtHdMn', 'lOD', 'lDITp', 'lITp', 'lODfr', 'lODfrMn', 'tID',
              'tGT', 'tFt', 'lICNm', 'lICVS', 'lGTSglO', 'lFtSglO', 'lODfrAI',
              'lODOIn', 'pFAllD', 'pFMn', 'pFSD', 'pFTrD', 'pFDvD',
              'pFCrDvAllD', 'pFCrDvMn', 'lPFTACDAllD', 'lPFTACDMn', 'OD',
              'ODfr', 'dOITp', 'cOITp', 'dClAlgDfr', 'dClAlgDfr4Corr',
              'dClResData', 'lSRowTrD', 'lSRowDvD', 'dfrRCrDv', 'dfrRCrDvMn',
              'ODfrTr', 'ODfrTrT', 'pRF', 'pRFF', 'nClst', 'isClRD']
        for cA in lA:
            if hasattr(self, cA):
                self.dOIn[cA] = getattr(self, cA)

    def updateDOInBasicInfo(self):
        self.dOIn['idO'] = self.idO
        self.dOIn['descO'] = self.descO
    
    def extractDFNmComp(self, lD = None): 
        dT = self.dfrFNmComp.to_dict('index')
        self.dFNmComp = {[tuple(cV.values()) for cV in dT.values()][0]: lD}
    
    def updateFNmComp(self, lD = None, sMod = 'Tr'):
        self.dfrFNmComp.loc[0, self.dIG['nmID']] = self.idO
        self.dfrFNmComp.loc[0, self.dIG['nmTrDt']] = self.sTrDt
        if sMod == 'Dv':
            self.dfrFNmComp.loc[0, self.dIG['nmTrDt']] = self.sDvDt
        self.dfrFNmComp.loc[0, self.dIG['nmTypeX']] = self.tpX
        self.dfrFNmComp.loc[0, self.dIG['nmIDGTFt']] = self.sGTFt
        self.dfrFNmComp.loc[0, self.dIG['nmFilter']] = self.nmFE
        self.extractDFNmComp(lD)

    def updateDOInFNm(self):
        lA = ['dfrFNmComp', 'dFNmComp', 'sFNm', 'sFNmTr', 'sFNmDv']
        for cA in lA:
            if hasattr(self, cA):
                self.dOIn[cA] = getattr(self, cA)
    
    def updateSFName(self, lD = None, dS = {}, setIt = False):
        if setIt:
            self.dfrFNmComp = GF.iniPdDfr(self.dIG['dFNmComp'], lSNmR = [0])
        lSFNm, lKDS = [], ['No', 'Tr', 'Dv']
        GF.complDict(dS, lK = lKDS)
        for s in lKDS:
            self.updateFNmComp(lD = lD, sMod = s)
            lSFNm.append(TF.getSFNm(self.dFNmComp,
                                    dISPr = self.dITp['dISPr'][s], lS = dS[s]))
        self.sFNm, self.sFNmTr, self.sFNmDv = lSFNm
        self.updateDOInFNm()
        
    def iniTransData(self, sMid = '|'):
        # initialises the transformed data from the DataFrame of mean values
        iB, lI = self.dITp['iFtrBase'], list(range(len(self.cDfrMn.index)))
        sRef, self.lSRowTrD = self.cDfrMn.index[iB], list(self.cDfrMn.index)
        if self.dITp['transD'] in ['Abs', 'Rel', '2LQ']:
            lI.remove(iB)
            self.lSRowTrD = [self.cDfrMn.index[k] + sMid + sRef for k in lI]
            self.cDfrTr = GF.iniPdDfr(lSNmC = self.cDfrMn.columns,
                                      lSNmR = self.lSRowTrD)
        else:
            self.cDfrTr = self.cDfrMn.copy()
        return iB, lI

    def saveTransData(self, dS = {}, dfrAI = None):
        dISPrTr, lKDS = self.dITp['dISPr']['Tr'], ['No', 'Tr', 'Dv']
        GF.complDict(dS, lK = lKDS)
        self.updateFNmComp(sMod = 'Tr')
        self.sFNmTr = TF.getSFNm(self.dFNmComp, dISPr = dISPrTr, lS = dS['Tr'])
        self.updateDOInFNm()
        if dfrAI is None:
            dfrAI = self.addIDfr
        self.pFTrD = TF.saveTrData(self.dITp, self.dOIn, self.cDfrTr, dfrAI,
                                   sMn = self.dIG['nmMeans'])

    def calcTransData(self, dS = {}, dfrAI = None, sMid = '|'):
        iB, lI = self.iniTransData(sMid = sMid)
        TF.fillTransData(self.cDfrTr, self.cDfrMn, self.dITp['transD'],
                         self.lSRowTrD, iB, lI)
        stdOp, dropAx = self.dITp['stdOp'], self.dITp['dropAx_Mn']
        if dropAx is not None:
            self.cDfrTr.dropna(axis = dropAx, inplace = True)
            self.attrDfrMn.dropna(axis = 1 - dropAx, inplace = True)
            lID = list(self.attrDfrMn.loc[:, self.dIG['idOAtDfr']])
            self.attrDfrMn.loc[:, self.dIG['nmWtAtDfr']] = GF.adjustWts(lID)
        self.cDfrTr = GF.standardPdDfr(self.cDfrTr, stdOp = stdOp)
        self.cDfrTrT = self.cDfrTr.T
        self.saveTransData(dS = dS, dfrAI = dfrAI)

    def saveDevData(self, dS = {}, dfrAI = None):
        dISPrDv, lKDS = self.dITp['dISPr']['Dv'], ['No', 'Tr', 'Dv']
        GF.complDict(dS, lK = lKDS)
        self.updateFNmComp(sMod = 'Dv')
        self.sFNmDv = TF.getSFNm(self.dFNmComp, dISPr = dISPrDv, lS = dS['Dv'])
        self.updateDOInFNm()
        if dfrAI is None:
            dfrAI = self.addIDfr
        self.pFDvD = TF.saveDvData(self.dITp, self.dOIn, self.cDfrDv, dfrAI,
                                   sMn = self.dIG['nmMeans'])
    
    def calcDevData(self, dS = {}, dfrAI = None, sMid = '|'):
        self.cDfrDv = GF.getDfrDv(self.cDfrMn, self.cDfrSD, self.dfrN,
                                  self.dITp['devTp'], sSt = self.idO + '__',
                                  sMid = sMid)
        self.cDfrDvT = self.cDfrDv.T
        self.saveDevData(dS = dS, dfrAI = dfrAI)

###############################################################################

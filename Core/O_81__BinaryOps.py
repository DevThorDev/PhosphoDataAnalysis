# -*- coding: utf-8 -*-
###############################################################################
# --- O_81__BinaryOps.py ------------------------------------------------------
###############################################################################
import os

import Core.C_00__GenConstants as GC
import Core.F_00__GenFunctions as GF
import Core.F_02__PltFunctions as PF
import Core.F_03__OTpFunctions as TF

from Core.O_00__DataBaseClass import DataBaseClass

class BinaryOps(DataBaseClass):
    def __init__(self, inpDat, iTp, lODat, nmFAdd = '', sSep = '_'):
        super().__init__(inpDat, iTp, nmFAdd = nmFAdd)
        self.idO = GC.S_BIN_OP
        self.descO = 'Binary operations'
        assert len(lODat) >= 2
        self.lOD = TF.deepCopyLO(lODat)
        self.tpX = GF.get1stFromL([OD.tpX for OD in self.lOD])
        self.cXFt = GF.get1stFromL([OD.cXFt for OD in self.lOD])
        self.lDITp = [OD.dITp for OD in self.lOD]
        self.lITp = [cD['iTp'] for cD in self.lDITp]
        self.lODfr = [OD.cDfr for OD in self.lOD]
        self.lODfrMn = [OD.cDfrMn for OD in self.lOD]
        self.tID = tuple([OD.idO for OD in self.lOD])
        self.tGT = tuple([OD.cGT for OD in self.lOD])
        self.tFt = tuple([OD.cFt for OD in self.lOD])
        self.cGT = sSep.join(self.tGT)
        self.cFt = sSep.join(self.tFt)
        self.updateDITp()
        self.isClRD = GF.get1stFromL([OD.isClRD for OD in self.lOD])
        self.updateDOIn()
        self.updateSFName(lD = [OD.dFNmComp for OD in self.lOD],
                          dS = {'No': [OD.sFNm for OD in self.lOD],
                                'Tr': [OD.sFNmTr for OD in self.lOD],
                                'Dv': [OD.sFNmDv for OD in self.lOD]})
        print('Initiated "BinaryOps" base object.')
        
    def printObjInfo(self):
        super().printObjInfo()
        print('List of type indices of data:', self.lITp)
        print('Tuple of data IDs:', self.tID)
        print('Tuple of genotypes:', self.tGT)
        print('Tuple of features:', self.tFt)
    
    def extHeader(self, sSep = '_'):
        for OD in self.lOD:
            self.dITp['lSHeadRes'] += OD.cDfrDvT.columns.to_list()
        lCmp = GF.cutNSplFromL(list(self.lOD[0].cDfrDvT.columns))
        self.dITp['lSHeadRes'] += [self.dITp['sDvSc'] + sSep + s for s in lCmp]
        self.dITp['lSHeadRes'] += [self.dITp['sCI'] + sSep + s for s in lCmp]
        self.dITp['lSHeadRes'] += list(self.dITp['dSDvScX'])
        self.dITp['lSHeadRes'] += list(self.dITp['dSCIX'])
        for cThr in (list(reversed([-x for x in self.dITp['lPosCIBnd']])) +
                     self.dITp['lPosCIBnd']):
            self.dITp['lSHeadRes'] += [self.dITp['sOccCI'] + sSep + str(cThr)]
    
    def updateDITp(self, sSep = '_'):
        self.dITp['tCNmRes'] = (self.dITp['sCt'], self.tID[0], self.tID[1])
        self.dITp['nmObj1'], self.dITp['nmObj2'] = self.dITp['tCNmRes'][1:3]
        self.dITp['lSHeadRes'] = list(self.dITp['tCNmRes'])
        self.dITp['lSHeadResMn'] = list(self.dITp['tCNmRes'])
        if not self.isClRD:         # otherwise, no dev. data exists
            self.extHeader(sSep)    # add deviation data to the header

    def updateDOIn(self):
        self.updateDOInBasicInfo()
        lA = ['tpX', 'cXFt', 'lDITp', 'lITp', 'lODfr', 'lODfrMn', 'tID', 'tGT',
              'tFt', 'cGT', 'cFt', 'isClRD', 'nmFE', 'pFCrDvAllD', 'pFCrDvMn',
              'lPFTACDAllD', 'lPFTACDMn']
        for cA in lA:
            if hasattr(self, cA):
                self.dOIn[cA] = getattr(self, cA)
    
    def importBasicData(self, lSHdR, sMn = ''):
        lDfrDvT, lDfrAIC, lenHdB = [], [], len(self.dITp['tCNmRes'])
        if not self.isClRD:         # otherwise, no dev. data exists
            for i in range(2):
                pDvF = TF.getPFRes(self.lOD[i].dITp, self.lOD[i].dOIn,
                                   sMn = self.dIG['nmMeans'], sMod = 'Dv')
                lHStr, iS = GF.getLHStrISpl(self.lOD[i].dITp, 'lHdrAIC')
                tDfr = GF.readCSV(pDvF, self.dITp['cSep'], iCol = 0,
                                  lHStr = lHStr, iSp = iS)
                lDfrDvT.append(tDfr[0])
                lDfrAIC.append(tDfr[1])
            lSt = (lSHdR[:lenHdB] + lDfrAIC[0].columns.to_list() +
                   lDfrAIC[1].columns.to_list() + self.dITp['lSCorrAll'])
            if sMn == self.dIG['nmMeans']:
                lSHdR = lSt
            else:
                lSHdR = lSt + lSHdR[lenHdB:]
        else:
            lSHdR = lSHdR[:lenHdB] +  self.dITp['lSCorrAll'] + lSHdR[lenHdB:]
        return lDfrDvT, lDfrAIC, lSHdR

    def calcCrDv(self, sMn = ''):
        # calculate the correlation between columns of the first DataFrame
        # with columns of the second (all combinations)
        dR, lSHdR, nmMn = {}, self.dITp['lSHeadRes'], self.dIG['nmMeans']
        lDfr = [self.lOD[0].cDfr, self.lOD[1].cDfr]
        if sMn == nmMn:
            lSHdR = self.dITp['lSHeadResMn']
            lDfr = [self.lOD[0].cDfrMn, self.lOD[1].cDfrMn]
        pRF = TF.getPFRes(self.dITp, self.dOIn, sSt = self.dITp['sCorrS'],
                          sMn = sMn, sMod = 'Dv')
        # read the MD / PD deviation data from the respective DataFrames
        lDfrDvT, lDfrAIC, lSHdR = self.importBasicData(lSHdR, sMn = sMn)
        if not os.path.isfile(pRF):    # calculate the correlations
            # fill the result dictionary
            TF.fillDictRCrDv(self.dIG, self.dITp, dR, lDfr, lDfrDvT, lDfrAIC,
                             sMn)
            # write the resulting dictionary to a result file
            TF.writeCrDvDatToF(self.dITp, dR, pRF, lSHdR)
        # dfrR contains all correlations in a DataFrame - read from csv file
        dfrR = GF.readCSV(pRF, self.dITp['cSep'])
        if len(dR) == 0 and self.dITp['calcDictR']:
            lNmCK = [self.dITp['nmObj1'], self.dITp['nmObj2']]
            dR = GF.convPdDfrToDict(dfrR, lNmCK,  self.dITp['lSCorrAll'])
        return dfrR, dR, lDfrAIC, pRF

    def getDfrIn(self, tC, sMn = ''):
        if tC[0] in self.dIG['lSPhoCl'] or tC[1] in self.dIG['lSPhoCl']:
            if sMn == self.dIG['nmMeans']:
                dfrIn = GF.readCSV(self.pFCrDvMn, self.dITp['cSep'],
                                   lHStr = self.dIG['lSPhoCl'])
            else:
                dfrIn = GF.readCSV(self.pFCrDvAllD, self.dITp['cSep'],
                                   lHStr = self.dIG['lSPhoCl'])
        elif tC[0] == self.dIG['sPhoD'] or tC[1] == self.dIG['sPhoD']:
            if sMn == self.dIG['nmMeans']:
                dfrIn = self.dfrRCrDvMn
            else:
                dfrIn = self.dfrRCrDv
        else:
            print('ERROR: Tuple of columns', tC, 'faulty.')
            assert False
        return dfrIn

    def getTFNmInfo(self):
        tLFNC = ([OD.dFNmComp for OD in self.lOD],
                 [OD.dFNmComp for OD in reversed(self.lOD)])
        tLSNo = ([OD.sFNm for OD in self.lOD],
                 [OD.sFNm for OD in reversed(self.lOD)])
        tLSTr = ([OD.sFNmTr for OD in self.lOD],
                 [OD.sFNmTr for OD in reversed(self.lOD)])
        tLSDv = ([OD.sFNmDv for OD in self.lOD],
                 [OD.sFNmDv for OD in reversed(self.lOD)])
        return tLFNC, tLSNo, tLSTr, tLSDv

    def getFNmTACD(self, tC, nmStBase, lDfrAIC, i, sMn = ''):
        assert len(lDfrAIC) == 0 or len(lDfrAIC) >= 2
        sSt, sC = nmStBase, ''
        if len(lDfrAIC) >= 2:
            if tC[0] in lDfrAIC[i].columns:
                sSt = GF.combToStr([nmStBase, tC[0]])
                sC = tC[0]
            elif tC[1] in lDfrAIC[i].columns:
                sSt = GF.combToStr([nmStBase, tC[1]])
                sC = tC[1]
        lS = self.lOD[i].cDfr.columns.to_list()
        if len(sC) > 0:     # sC is a column header of the AIC DataFrame
            lS = sorted(list(set(lDfrAIC[i].loc[:, sC])))
        return TF.getPFRes(self.dITp, self.dOIn, sSt, sMn, sMod = 'Dv'), lS, sC

    def calcTACD4CC(self, tCWt, lDfrAIC, sMn = ''):
        lPRF, nmTA, tC, doWt = [], self.dITp['nmTACDRF'], tCWt[0], tCWt[1]
        dfrIn = self.getDfrIn(tC, sMn = sMn)
        tLFNC, tLSNo, tLSTr, tLSDv = self.getTFNmInfo()
        for i in range(2):
            self.updateSFName(lD = tLFNC[i],
                              dS = {'No': tLSNo[i], 'Tr': tLSTr[i],
                                    'Dv': tLSDv[i]})
            pRF, lNmI, sCol = self.getFNmTACD(tC, nmTA, lDfrAIC, i, sMn = sMn)
            if not os.path.isfile(pRF):
                dfrTACD = TF.calcTACD(self.dITp, self.dOIn, dfrIn, lNmI, i,
                                      sCol, doWt, sMn)
                TF.saveTACDData(self.dITp, dfrTACD, lDfrAIC, i,
                                self.lDITp[i]['lHdrAIC'], pRF, len(sCol) == 0)
            lPRF.append(pRF)
        self.updateSFName(lD = [OD.dFNmComp for OD in self.lOD],
                          dS = {'No': [OD.sFNm for OD in self.lOD],
                                'Tr': [OD.sFNmTr for OD in self.lOD],
                                'Dv': [OD.sFNmDv for OD in self.lOD]})
        return lPRF
    
    def doCorrDev(self, doAllD = True, doMeans = True):
        self.lPFTACDAllD, self.lPFTACDMn, nmMn = [], [], self.dIG['nmMeans']
        self.updateDOIn()
        if doAllD:
            cT = self.calcCrDv()
            self.dfrRCrDv, self.dRCrDv, lDfrAIC, self.pFCrDvAllD = cT
            for tColWt in self.dITp['lCalcTACD']:
                self.lPFTACDAllD += self.calcTACD4CC(tColWt, lDfrAIC)
            if (self.lOD[0].dITp['sNmSpec'] in self.dIG['lSBasDLTf'] or
                self.lOD[1].dITp['sNmSpec'] in self.dIG['lSBasDLTf']):
                PF.pltHist1C(self.dITp, self.dOIn, self.dfrRCrDv,
                             self.pFCrDvAllD)
                # PF.pltSCorr(self.dIG, self.dITp, self.dOIn, self.dfrRCrDv,
                #             self.pFCrDvAllD, self.lOD)
        if doMeans:
            cT = self.calcCrDv(nmMn)
            self.dfrRCrDvMn, self.dRCrDvMn, lDfrAIC, self.pFCrDvMn = cT
            for tColWt in self.dITp['lCalcTACD']:
                self.lPFTACDMn += self.calcTACD4CC(tColWt, lDfrAIC, nmMn)
            if (self.lOD[0].dITp['sNmSpec'] in self.dIG['lSBasDLTf'] or
                self.lOD[1].dITp['sNmSpec'] in self.dIG['lSBasDLTf']):
                PF.pltHist1C(self.dITp, self.dOIn, self.dfrRCrDvMn,
                             self.pFCrDvMn, useMns = True)
        self.updateDOIn()

###############################################################################

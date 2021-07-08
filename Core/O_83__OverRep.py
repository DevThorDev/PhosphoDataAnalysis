# -*- coding: utf-8 -*-
###############################################################################
# --- O_83__OverRep.py --------------------------------------------------------
###############################################################################
import os

import Core.C_00__GenConstants as GC
import Core.F_00__GenFunctions as GF
import Core.F_02__PltFunctions as PF
import Core.F_03__OTpFunctions as TF

from Core.O_00__DataBaseClass import DataBaseClass

class OverRep(DataBaseClass):
    def __init__(self, inpDat, iTp, ODat, nmFAdd = '', sSep = '_'):
        super().__init__(inpDat, iTp, nmFAdd = nmFAdd)
        self.idO = GC.S_OV_REP
        self.descO = 'Over-representation'
        self.OD = ODat
        self.tpX = self.OD.tpX
        self.cXFt = self.OD.cXFt
        self.dSort = {}
        self.lSrt, self.lAsc, self.lElCat = [], [], []
        self.nMin, self.nMax = 1, 1
        if self.OD.idO == 'BinOp':          # a "BinOps" object
            assert len(self.OD.lOD) >= 2
            self.iniIfBinObj(sSep = sSep)
        else:                               # a basic Pho data object
            self.iniIfSglObj()
        self.updateDOIn()
        self.updateSFName(lD = [self.OD.dFNmComp],
                          dS = {'No': [self.OD.sFNm],
                                'Tr': [self.OD.sFNmTr],
                                'Dv': [self.OD.sFNmDv]})
        self.getPResF()
        self.updateDOIn()
        print('Initiated "OverRep" base object with index', str(self.iTpOD))

    def iniIfSglObj(self):
        self.iTpOD = 0      # a basic Pho data object
        self.lDITp = [self.OD.dITp]
        self.lITp = [self.OD.dITp['iTp']]
        self.lODfr = [self.OD.cDfr]
        self.lODfrMn = [self.OD.cDfrMn]
        self.tID = tuple([self.OD.idO])
        self.tGT = tuple([self.OD.cGT])
        self.tFt = tuple([self.OD.cFt])
        self.cGT = self.OD.cGT
        self.cFt = self.OD.cFt

    def iniIfBinObj(self, sSep = '_'):
        self.iTpOD = 1      # a "BinOps" object
        self.lDITp = [OD.dITp for OD in self.OD.lOD]
        self.lITp = [cD['iTp'] for cD in self.OD.lDITp]
        self.lODfr = [OD.cDfr for OD in self.OD.lOD]
        self.lODfrMn = [OD.cDfrMn for OD in self.OD.lOD]
        self.tID = tuple([OD.idO for OD in self.OD.lOD])
        self.tGT = tuple([OD.cGT for OD in self.OD.lOD])
        self.tFt = tuple([OD.cFt for OD in self.OD.lOD])
        self.cGT = sSep.join(self.tGT)
        self.cFt = sSep.join(self.tFt)

    def printObjInfo(self):
        super().printObjInfo()
        print('Index denoting type of object data:', self.iTpOD)
        print('List of type indices of object data:', self.lITp)
        for k, pRF in enumerate(self.lPRF):
            print('Path of result file ', k + 1, ':', pRF, sep = '')

    def getPResF(self, sMn = '', sMod = 'No'):
        self.lPRF, lDef = [], None
        if self.iTpOD == 1:
            lDef = [self.OD.lOD[k].idO for k in range(2)]
        t = TF.getParProf(self.dITp, self.iTpOD, self.OD.isClRD, lDef)
        self.dSort, self.lSrt, self.lAsc, self.lElCat, self.nMin, self.nMax = t
        sEPr = str(self.nMin) + '_' + str(self.nMax) + '_'
        for sHd in self.dSort[GC.S_IDX]:
            sEPr += sHd + str(int(self.dSort[GC.S_IDX][sHd]['Asc'])) + '_'
        for sCat in self.lElCat:
            sEPr = GF.addString(str(sCat), sPost = '_') + sEPr
        for sE in [self.dITp['sNOcc'], self.dITp['sPValOvPOf'],
                   self.dITp['sPValUnPOf']]:
            sE = GF.addString(sEPr + self.dITp['sMCorrectS'], sPost = '_') + sE
            self.lPRF.append(TF.getPFRes(self.dITp, self.dOIn, sMn = sMn,
                                         sEnd = sE, sMod = sMod))

    def updateDOIn(self):
        self.updateDOInBasicInfo()
        lA = ['tpX', 'cXFt', 'dSort', 'lSrt', 'lAsc', 'lElCat', 'nMin', 'nMax',
              'iTpOD', 'lDITp', 'lITp', 'lODfr', 'lODfrMn', 'tID', 'tGT',
              'tFt', 'cGT', 'cFt', 'nmFE', 'lPRF']
        for cA in lA:
            if hasattr(self, cA):
                self.dOIn[cA] = getattr(self, cA)

    def plotProfiles(self, dfrR, k = 0, tpPr = 'PD'):
        d = TF.selDataThr(self.dITp, dfrR, self.dSort, k)
        PF.plotProfile(self.dITp, d, self.lPRF[k], k, tpPr = tpPr)

    def plotOnlyIfCalc(self, tpPr = 'PD'):
        lCalc = [False]*len(self.lPRF)
        for k, pRF in enumerate(self.lPRF):
            if os.path.isfile(pRF):
                lCalc[k] = True
        if lCalc == [True]*len(self.lPRF):
            for k in range(len(self.dITp['lTpY'])):
                dfrR = GF.readCSV(self.lPRF[k], sepD = self.dITp['cSep'],
                                  iCol = 0)
                self.plotProfiles(dfrR, k, tpPr = tpPr)
            return True
        else:
            return False

    def getPFData(self, useMn = True):
        if self.iTpOD == 0:      # a basic Pho data object
            if useMn:
                pFD = self.OD.pFMn
            else:
                pFD = self.OD.pFAllD
        elif self.iTpOD == 1:    # a "BinOps" from basic data object
            if useMn:
                pFD = self.OD.pFCrDvMn
            else:
                pFD = self.OD.pFCrDvAllD
        else:
            print('ERROR: Index of object data type (' + str(self.iTpOD) +
                  ') not implemented.')
            assert False
        return pFD

    def calcProfiles(self, tpPr = 'PD', useMn = True):
        [dOccAbs, dPValOv, dPValUn] = [{}, {}, {}]
        if self.plotOnlyIfCalc(tpPr = tpPr):
            return None
        t = GF.getSrtData(self.getPFData(useMn = useMn), self.lSrt, self.lAsc,
                          self.lElCat, sepD = self.dITp['cSep'])
        dfrRd, lSerVC, llAttr, lNAttr, N = t
        TF.calcPValProfiles(self.dITp, dOccAbs, dPValOv, dPValUn, dfrRd,
                            lSerVC, llAttr, lNAttr, self.nMin, self.nMax, N,
                            self.lElCat)
        for k, cD in enumerate([dOccAbs, dPValOv, dPValUn]):
            dfrR = TF.saveDictAsPdDfr(self.dITp, self.lPRF, self.dSort, cD,
                                      self.nMin, self.nMax, k)
            self.plotProfiles(dfrR, k, tpPr = tpPr)

###############################################################################

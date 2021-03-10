# -*- coding: utf-8 -*-
###############################################################################
# --- O_31__CombData.py -------------------------------------------------------
###############################################################################
import Core.C_00__GenConstants as GC
import Core.F_00__GenFunctions as GF
import Core.F_03__OTpFunctions as TF

from Core.O_00__DataBaseClass import DataBaseClass

class CombData(DataBaseClass):
    def __init__(self, inpDat, iTp, lODat):
        super().__init__(inpDat, iTp)
        self.idO = GC.S_COMB_D
        self.descO = 'Combined data'
        assert len(lODat) >= 2
        self.lOD = TF.deepCopyLO(lODat)
        self.tpX = GF.get1stFromL([OD.tpX for OD in self.lOD])
        self.cXFt = GF.get1stFromL([OD.cXFt for OD in self.lOD])
        self.lDITp = [OD.dITp for OD in self.lOD]
        self.lITp = [cD['iTp'] for cD in self.lDITp]
        self.lODfr = [OD.cDfr for OD in self.lOD]
        self.lODfrMn = [OD.cDfrMn for OD in self.lOD]
        self.lODfrSD = [OD.cDfrSD for OD in self.lOD]
        self.lODfrN = [OD.dfrN for OD in self.lOD]
        self.getFeatureHdr()
        self.tID = tuple([OD.idO for OD in self.lOD])
        self.tGT = tuple([OD.cGT for OD in self.lOD])
        self.tFt = tuple([OD.cFt for OD in self.lOD])
        self.complementAttr()
        TF.addToDITp(self.dITp, self.lDITp, lKOvwr = self.dIG['lKModD'])
        self.sTrDt += ''.join([self.dITp['transD'], self.dITp['stdOp']])
        self.sDvDt += ''.join([self.dITp['devTp']])
        self.lICNm = [OD.dITp['iColNm'] for OD in self.lOD]
        self.lICVS = [OD.dITp['iColVsS'] for OD in self.lOD]
        self.updateDOIn()
        self.updateSFName(lD = [OD.dFNmComp for OD in self.lOD],
                          dS = {'No': [OD.sFNm for OD in self.lOD],
                                'Tr': [OD.sFNmTr for OD in self.lOD],
                                'Dv': [OD.sFNmDv for OD in self.lOD]})
        print('Initiated "CombData" base object.')
    
    def getFeatureHdr(self):
        self.lFtHd, self.lFtHdMn, lFtHdTot, lFtHdMnTot = [], [], [], []
        for OD in self.lOD:
            lFtHdTot += OD.lFtHd
            lFtHdMnTot += OD.lFtHdMn
        GF.addToListUnique(self.lFtHd, lFtHdTot)
        GF.addToListUnique(self.lFtHdMn, lFtHdMnTot)
    
    def printObjInfo(self):
        super().printObjInfo()
        print('List of type indices of data:', self.lITp)
        print('Tuple of data IDs:', self.tID)
        print('Tuple of genotypes:', self.tGT)
        print('Tuple of features:', self.tFt)
        print('List of indices of columns containing the names:', self.lICNm)
        print('List of starting indices of values columns:', self.lICVS)
        
    def complementAttr(self, sSep = '_'):
        self.cGT = sSep.join(self.tGT)
        self.cFt = sSep.join(self.tFt)
        self.nmFE = ''
        if GF.listsEqual([OD.nmFE for OD in self.lOD]):
            self.nmFE = self.lOD[0].nmFE
        else:
            self.nmFE = sSep.join(OD.nmFE for OD in self.lOD)
        self.lGTSglO = list(self.tGT)
        self.lFtSglO = list(self.tFt)

    def updateDOIn(self):
        self.updateDOInBasicInfo()
        lA = ['sTrDt', 'sDvDt', 'tpX', 'cXFt', 'lDITp', 'lITp', 'lODfr',
              'lODfrMn', 'lODfrSD', 'tGT', 'tFt', 'cGT', 'cFt', 'nmFE',
              'lGTSglO', 'lFtSglO', 'cDfr', 'cDfrMn', 'cDfrSD', 'dfrN',
              'attrDfr', 'attrDfrMn', 'lICNm', 'lICVS', 'pFAllD', 'pFMn',
              'lFtHd', 'lFtHdMn', 'lSRowTrD', 'lSRowDvD', 'cDfrTr', 'cDfrTrT',
              'cDfrDv', 'cDfrDvT', 'pFTrD', 'pFDvD']
        for cA in lA:
            if hasattr(self, cA):
                self.dOIn[cA] = getattr(self, cA)

    def combineData(self):
        self.mergeO()
        self.updateDOIn()
        self.saveData()
        self.calcTransData(dS = {'Tr': [OD.sFNmTr for OD in self.lOD]})
        self.calcDevData(dS = {'Dv': [OD.sFNmDv for OD in self.lOD]})
        self.updateDOIn()

    def mergeO(self):
        dropAx_F, dropAx_Mn = self.dITp['dropAx_F'], self.dITp['dropAx_Mn']
        self.attrDfr = GF.concPdDfrS([OD.attrDfr for OD in self.lOD],
                                     ignIdx = True)
        self.attrDfrMn = GF.concPdDfrS([OD.attrDfrMn for OD in self.lOD],
                                       ignIdx = True)
        t = GF.concPdDfrWt(self.lODfr, concAx = 1, dropAx = dropAx_F)
        self.cDfr, self.attrDfr.loc[:, self.nmWtC] = t
        t = GF.concPdDfrWt(self.lODfrMn, concAx = 1, dropAx = dropAx_Mn)
        self.cDfrMn, self.attrDfrMn.loc[:, self.nmWtC] = t
        self.cDfrSD = GF.concPdDfrS(self.lODfrSD, concAx = 1,
                                    dropAx = dropAx_Mn)
        self.dfrN = GF.concPdDfrS(self.lODfrN, concAx = 1, dropAx = dropAx_Mn)

    def saveData(self):
        self.pFAllD = TF.saveData(self.dITp, self.dOIn, self.cDfr.T, '')
        self.pFMn = TF.saveData(self.dITp, self.dOIn, self.cDfrMn.T, '',
                                self.dIG['nmMeans'])

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class CombDataC(CombData):
    def __init__(self, inpDat, iTp, lODat):
        super().__init__(inpDat, iTp, lODat)
        self.idO = GC.S_COMB_D + 'C'
        self.descO = 'Combined data (feature type "C")'
        self.updateDOInBasicInfo()
        self.updateSFName(lD = [OD.dFNmComp for OD in self.lOD],
                          dS = {'No': [OD.sFNm for OD in self.lOD],
                                'Tr': [OD.sFNmTr for OD in self.lOD],
                                'Dv': [OD.sFNmDv for OD in self.lOD]})
        print('Initiated "CombDataC" base object.')

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class CombDataG(CombData):
    def __init__(self, inpDat, iTp, lODat):
        super().__init__(inpDat, iTp, lODat)
        self.idO = GC.S_COMB_D + 'G'
        self.descO = 'Combined data (feature type "G")'
        self.updateDOInBasicInfo()
        self.updateSFName(lD = [OD.dFNmComp for OD in self.lOD],
                          dS = {'No': [OD.sFNm for OD in self.lOD],
                                'Tr': [OD.sFNmTr for OD in self.lOD],
                                'Dv': [OD.sFNmDv for OD in self.lOD]})
        print('Initiated "CombDataG" base object.')

###############################################################################

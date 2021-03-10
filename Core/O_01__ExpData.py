# -*- coding: utf-8 -*-
###############################################################################
# --- O_01__ExpData.py --------------------------------------------------------
###############################################################################
import Core.C_00__GenConstants as GC
import Core.F_00__GenFunctions as GF
import Core.F_03__OTpFunctions as TF

from Core.O_00__DataBaseClass import DataBaseClass

class ExpData(DataBaseClass):
    def __init__(self, inpDat, iTp = 1, cXFrt = '', cGT = '', cFt = '',
                 nmFAdd = '', cQuery = ''):
        super().__init__(inpDat, iTp, cGT = cGT, cFt = cFt, nmFAdd = nmFAdd)
        self.idO = GC.S_EXP_D
        self.descO = 'Experimental data'
        self.iCNm = self.dITp['iColNm']
        self.iCVS = self.dITp['iColVsS']
        self.lIC2Str = self.dITp['lIColConv2Str']
        self.lCHdrAIC = self.dITp['lHdrAIC']
        self.lCHdrAIR = self.dITp['lHdrAIR']
        self.cQry = cQuery
        self.tIMM = self.dITp['tIMnMx']
        self.sTrDt += ''.join([self.dITp['transD'], self.dITp['stdOp']])
        self.sDvDt += ''.join([self.dITp['devTp']])
        self.tpX = None
        self.cXFt = cXFrt
        self.dITp['pRelDF'] = self.dITp['dPRelDF'][self.cGT]
        self.updateDOIn()
        self.updateSFName()
        print('Initiated "ExpData" base object.')
        
    def printObjInfo(self):
        super().printObjInfo()
        print('Index of column containing the names:', self.iCNm)
        print('Starting index of column containing the values:', self.iCVS)
        print('List of indices of columns converted to string:', self.lIC2Str)
        print('Column header list of classification info:', self.lCHdrAIC)
        print('Column header list of retained additional info:', self.lCHdrAIR)
        print('Applied query to filter the data:', self.cQry)
        print('Minimum and maximum index for column selection:', self.tIMM)
        print('List of positions of duplicates in rows:', self.lRIDp[:20])
        print('Number of duplicates in rows:', len(self.lRIDp))
        print('List of positions of duplicates in columns:', self.lCIDp[:20])
        print('Number of duplicates in columns:', len(self.lCIDp))
        print('Path to DataFrame for "all data" file:', self.pFAllD)
        print('Path to DataFrame for "Means" file:', self.pFMn)
        print('Path to DataFrame for "SDs" file:',self.pFSD)
    
    def calcMeanSD(self, tIMnMx = (-1, -1)):
        cT = TF.loadPdDfr(self.dITp, self.dOIn, doT = True)
        self.cDfr, self.addIDfr, self.lRIDp, self.lCIDp = cT
        self.lFtHd, self.lFtHdMn = GF.getFeatureHdr(self.cDfr)
        # get attribute DataFrame for all data
        lAttrIni = TF.getLAttr(self.dIG, [self.idO, self.cGT, self.cFt, 1.],
                               iCGT = 1, iCFt = 2)
        cL = self.cDfr.shape[1]
        self.attrDfr = GF.iniAttrDfr(lAttrIni, self.lNmAttr, cLen = cL)
        # calculate mean and standard deviation for all columns
        self.updateDOIn()
        cT = TF.getMeanSD(self.dITp, self.dOIn, self.cDfr, self.addIDfr,
                          self.dIG['nmMeans'], self.dIG['nmSDs'],
                          self.dIG['lClrSpec'])
        (self.cDfrMn, self.cDfrSD, self.pFAllD, self.pFMn, self.pFSD,
         self.dfrN) = cT
        # initialise the attribute DataFrame (mean data)
        self.attrDfrMn = GF.iniAttrDfr(lAttrIni, self.lNmAttr,
                                       cLen = self.cDfrMn.shape[1])
        # in case not all column should be used, restrict data to iMin:iMax
        self.confineColumns(tIMnMx)
        self.calcTransData()
        self.calcDevData()
        self.updateDOIn()

    def confineColumns(self, tIMnMx = (-1, -1)):
         # confine columns to the start and end column indices (tIMnMx)
        assert len(tIMnMx) >= 2 and len(self.tIMM) >= 2
        iMn = GF.confineC(tIMnMx[0], self.tIMM[0])
        iMx = GF.confineC(tIMnMx[1], self.tIMM[1], True)
        self.dITp['tIMnMx'] = (iMn, iMx)
        self.tIMM = self.dITp['tIMnMx']
        self.cDfr = GF.shrinkPdDfr(self.cDfr, self.tIMM)
        self.cDfrMn = GF.shrinkPdDfr(self.cDfrMn, self.tIMM)
        self.cDfrSD = GF.shrinkPdDfr(self.cDfrSD, self.tIMM)

    def updateDOIn(self):
        self.updateDOInBasicInfo()
        lA = ['iCNm', 'iCVS', 'lIC2Str', 'lCHdrAIC', 'lCHdrAIR', 'cQry',
              'tIMM', 'sTrDt', 'sDvDt', 'tpX', 'cXFt', 'cDfr', 'cDfrMn',
              'cDfrSD', 'dfrN', 'pFAllD', 'pFMn', 'pFSD', 'addIDfr', 'lRIDp',
              'lCIDp', 'attrDfr', 'attrDfrMn', 'lFtHd', 'lFtHdMn', 'lSRowTrD',
              'lSRowDvD', 'cDfrTr', 'cDfrTrT', 'cDfrDv', 'cDfrDvT', 'pFTrD',
              'pFDvD']
        for cA in lA:
            if hasattr(self, cA):
                self.dOIn[cA] = getattr(self, cA)
        
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class ExpDataC(DataBaseClass):
    def __init__(self, inpDat, iTp = 1, cGT = '', cFt = '', nmFAdd = ''):
        super().__init__(inpDat, iTp, cGT = cGT, cFt = cFt, nmFAdd = nmFAdd)
        self.idO = GC.S_EXP_D + 'C'
        self.descO = 'Experimental data (feature type "C")'
        self.cXFt = self.dIG['nmXFtC']
        self.tIMM = self.dITp['tIMnMx']
        self.sTrDt += ''.join([self.dITp['transD'], self.dITp['stdOp']])
        self.sDvDt += ''.join([self.dITp['devTp']])
        self.updateDOIn()
        self.updateSFName()
        print('Initiated "ExpDataC" base object.')

    def updateDOIn(self):
        self.updateDOInBasicInfo()
        lA = ['tIMM', 'cXFt', 'lFtHd', 'lFtHdMn']
        for cA in lA:
            if hasattr(self, cA):
                self.dOIn[cA] = getattr(self, cA)
    
    def calcMeanSD(self, pdDfr, tIMnMx = (-1, -1)):
        # set the DataFrame, obtain and set the feature headers
        self.cDfr = pdDfr
        self.lFtHd, self.lFtHdMn = GF.getFeatureHdr(self.cDfr, lISpl = [0])
        # get attribute DataFrame for all data
        lAttrIni = TF.getLAttr(self.dIG, [self.idO, self.cGT, self.cFt, 1.],
                               iCGT = 1, iCFt = 2)
        cL = self.cDfr.shape[1]
        self.attrDfr = GF.iniAttrDfr(lAttrIni, self.lNmAttr, cLen = cL)
        self.updateDOIn()
        # calculate mean and standard deviation for all columns
        cT = TF.getMeanSD(self.dITp, self.dOIn, self.cDfr, self.addIDfr,
                          self.dIG['nmMeans'], self.dIG['nmSDs'],
                          self.dIG['lClrSpec'])
        (self.cDfrMn, self.cDfrSD, self.pFAllD, self.pFMn, self.pFSD,
         self.dfrN) = cT
        # initialise the attribute DataFrame (mean data)
        self.attrDfrMn = GF.iniAttrDfr(lAttrIni, self.lNmAttr,
                                       cLen = self.cDfrMn.shape[1])
        self.calcTransData()
        self.calcDevData()
        self.updateDOIn()

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class ExpDataG(ExpDataC):
    def __init__(self, inpDat, iTp = 1, cGT = '', cFt = '', nmFAdd = ''):
        super().__init__(inpDat, iTp, cGT = cGT, cFt = cFt, nmFAdd = nmFAdd)
        self.idO = GC.S_EXP_D + 'G'
        self.descO = 'Experimental data (feature type "G")'
        self.cXFt = self.dIG['nmXFtG']
        self.updateDOIn()
        self.updateSFName()
        print('Initiated "ExpDataG" base object.')

###############################################################################

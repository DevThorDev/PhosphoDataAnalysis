# -*- coding: utf-8 -*-
###############################################################################
# --- O_51__ClusterResData.py -------------------------------------------------
###############################################################################
import Core.C_00__GenConstants as GC
import Core.F_00__GenFunctions as GF
import Core.F_03__OTpFunctions as TF

from Core.O_00__DataBaseClass import DataBaseClass

class ClusterResData(DataBaseClass):
    def __init__(self, inpDat, iTp, ODat, idCl, nCl, cDfr, lSRTD, nmFAdd = ''):
        super().__init__(inpDat, iTp, nmFAdd = nmFAdd)
        self.idO = (GC.S_CLR_D + ODat.idO + '_' + idCl + '_' + str(nCl) +
                    self.dIG['nmClust'])
        self.descO = ('Clustering results of ' + ODat.idO + ' data with ' +
                      str(nCl) + ' clusters')
        self.OD = ODat
        self.tpX = self.OD.tpX
        self.cXFt = self.OD.cXFt
        self.cGT = self.OD.cGT
        self.cFt = self.OD.cFt
        self.nClst = nCl
        self.cDfr, self.cDfrMn = cDfr, cDfr
        self.lFtHd, self.lFtHdMn = lSRTD, lSRTD
        self.addIDfr = self.OD.addIDfr
        lAttrIni = TF.getLAttr(self.dIG, [self.idO, self.cGT, self.cFt, 1.],
                               iCGT = 1, iCFt = 2)
        self.attrDfr = GF.iniAttrDfr(lAttrIni, self.lNmAttr,
                                     cLen = self.cDfr.shape[1])
        self.attrDfrMn = GF.iniAttrDfr(lAttrIni, self.lNmAttr,
                                       cLen = self.cDfrMn.shape[1])
        lKeep = [self.dIG['nmIFtrBs'], self.dIG['nmTrD'], self.dIG['nmStdOp']]
        TF.complDITp(self.dITp, self.OD.dITp, lKOvwr = lKeep)
        self.isClRD = True
        self.updateDOIn()
        self.updateSFName(lD = [OD.dFNmComp for OD in self.lOD],
                          dS = {'No': [GF.cutNSpl(self.OD.sFNm)],
                                'Tr': [GF.cutNSpl(self.OD.sFNmTr)],
                                'Dv': [GF.cutNSpl(self.OD.sFNmDv)]})
        self.pF = TF.saveTrData(self.dITp, self.dOIn, self.cDfrMn,
                                self.addIDfr, sMn = self.dIG['nmMeans'])
        self.updateDOIn()
        print('Initiated "ClusterResData" base object.')
    
    def printObjInfo(self):
        super().printObjInfo()
        print('Number of clusters:', self.nClst)
        print('Path to cluster result data file:', self.pF)

    def updateDOIn(self):
        self.updateDOInBasicInfo()
        lA = ['tpX', 'cXFt', 'cGT', 'cFt', 'nClst', 'cDfr', 'cDfrMn',
              'addIDfr', 'attrDfr', 'attrDfrMn', 'lFtHd', 'lFtHdMn', 'pF']
        for cA in lA:
            if hasattr(self, cA):
                self.dOIn[cA] = getattr(self, cA)

###############################################################################

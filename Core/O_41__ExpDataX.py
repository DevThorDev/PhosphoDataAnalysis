# -*- coding: utf-8 -*-
###############################################################################
# --- O_41__ExpDataX.py -------------------------------------------------------
###############################################################################
import Core.C_00__GenConstants as GC
import Core.F_00__GenFunctions as GF
import Core.F_03__OTpFunctions as TF

from Core.O_00__DataBaseClass import DataBaseClass
from Core.O_01__ExpData import ExpDataC, ExpDataG
from Core.O_11__MetData import MetDataC, MetDataG
from Core.O_21__PhoData import PhoDataC, PhoDataG

class ExpDataX(DataBaseClass):
    def __init__(self, inpDat, iTp, lODat, typeX, extrDt = False):
        super().__init__(inpDat, iTp)
        self.idO = GC.S_EXP_DX
        self.descO = 'Experimental data extended'
        self.inpDG = inpDat
        self.tpX = typeX
        assert len(lODat) >= 1
        self.lOD = TF.deepCopyLO(lODat)
        self.cXFt = GF.get1stFromL([OD.cXFt for OD in self.lOD])
        self.XtrDt = extrDt
        self.lDITp = [OD.dITp for OD in self.lOD]
        self.lITp = [OD.dITp['iTp'] for OD in self.lOD]
        self.lODfr = [OD.cDfr for OD in self.lOD]
        self.lODfrMn = [OD.cDfrMn for OD in self.lOD]
        self.lODfrSD = [OD.cDfrSD for OD in self.lOD]
        self.lODfrN = [OD.dfrN for OD in self.lOD]
        self.lCHdrAIC = GF.get1stFromL([OD.dITp['lHdrAIC'] for OD in self.lOD])
        self.lODfrAIC = [OD.addIDfr.loc[self.lCHdrAIC, :] for OD in self.lOD]
        self.getFeatureHdr()
        self.complementAttr()
        self.setLAttrDfr()
        TF.prepareDITp(self.dITp, self.lOD, lKOvwr = self.dIG['lKModD'])
        self.sTrDt += ''.join([self.dITp['transD'], self.dITp['stdOp']])
        self.sDvDt += ''.join([self.dITp['devTp']])
        self.updateDOIn()
        self.updateSFName(lD = [OD.dFNmComp for OD in self.lOD])
        print('Initiated "ExpDataX" base object.')
        
    def printObjInfo(self):
        super().printObjInfo()
        print('Type of extension:', self.tpX)
        print('Tuple of data IDs:', self.tID)
        print('Tuple of genotypes:', self.tGT)
        print('Tuple of features:', self.tFt)
        print('List of genotypes of single objects:', self.lGTSglO)
        print('List of features of single objects:', self.lFtSglO)

    def printLAttrDfr(self, isMn = False):
        lDfr, sDfr = self.lAttrDfr, 'Attribute DataFrame (full data):'
        if isMn:
            lDfr, sDfr = self.lAttrDfrMn, 'Attribute DataFrame (means data):'
        for cDfr in lDfr:
            print('\n---', self.descO, '(', self.idO, ') - mode', self.cM, ':')
            print(sDfr, '\n', cDfr)
            print('Data types:', cDfr.dtypes)
    
    def getFeatureHdr(self):
        self.lFtHd, self.lFtHdMn, lFtHdTot, lFtHdMnTot = [], [], [], []
        for OD in self.lOD:
            lFtHdTot += OD.lFtHd
            lFtHdMnTot += OD.lFtHdMn
        GF.addToListUnique(self.lFtHd, lFtHdTot)
        GF.addToListUnique(self.lFtHdMn, lFtHdMnTot)

    def getSGenTpFeat(self, sSep = '_'):
        self.lGTSglO = GF.flattenIt([objGT.split(sSep) for objGT in self.tGT])
        self.lFtSglO = GF.flattenIt([objFt.split(sSep) for objFt in self.tFt])
        if (GF.setsEqual([self.lGTSglO, self.dIG['lIDGT']]) or
            (self.lGTSglO[0] == self.dIG['nmAllGT'] and
             GF.elEqual(self.lGTSglO))):
            self.cGT = self.dIG['nmAllGT']
        if (GF.setsEqual([self.lFtSglO, self.dIG['lIDFt']]) or
            (self.lFtSglO[0] == self.dIG['nmAllFt'] and
             GF.elEqual(self.lFtSglO))):
            self.cFt = self.dIG['nmAllFt']
        self.sGTFt = self.cGT + sSep + self.cFt
        if self.cGT == self.dIG['nmAllGT'] and self.cFt == self.dIG['nmAllFt']:
            self.sGTFt = self.dIG['nmAllGTFt']

    def complementAttr(self, sSep = '_'):
        self.tID = tuple([OD.idO for OD in self.lOD])
        self.tGT = tuple([OD.cGT for OD in self.lOD])
        self.tFt = tuple([OD.cFt for OD in self.lOD])
        self.cGT = sSep.join(self.tGT)
        self.cFt = sSep.join(self.tFt)
        self.getSGenTpFeat()
        self.lODOIn = [OD.dOIn for OD in self.lOD]
        self.nmFE = ''
        if GF.listsEqual([OD.nmFE for OD in self.lOD]):
            self.nmFE = self.lOD[0].nmFE
        else:
            self.nmFE = sSep.join(OD.nmFE for OD in self.lOD)
    
    def setLAttrDfr(self):
        self.lAttrDfr = [OD.attrDfr for OD in self.lOD]
        self.lAttrDfrMn = [OD.attrDfrMn for OD in self.lOD]

    def updateDOIn(self):
        self.updateDOInBasicInfo()
        lA = ['sTrDt', 'sDvDt', 'tpX', 'cXFt', 'XtrDt', 'tID', 'tGT', 'tFt',
              'cGT', 'cFt', 'nmFE', 'lGTSglO', 'lFtSglO', 'lCHdrAIC', 'dfrAIC',
              'attrDfr', 'attrDfrMn', 'lAttrDfr', 'lAttrDfrMn', 'lODOIn',
              'cDfr', 'cDfrMn', 'cDfrSD', 'dfrN', 'pFAllD', 'pFMn', 'lFtHd',
              'lFtHdMn', 'lSRowTrD', 'lSRowDvD', 'cDfrTr', 'cDfrTrT', 'cDfrDv',
              'cDfrDvT', 'pFTrD', 'pFDvD']
        for cA in lA:
            if hasattr(self, cA):
                self.dOIn[cA] = getattr(self, cA)

    def compareObj(self):
        self.iniTransData()
        cDfrCmp = GF.iniPdDfr(lSNmC = self.cDfrTr.columns,
                              lSNmR = self.cDfrTr.index)
        iB, j, sMn = self.dITp['iObjRef'], 0, self.dIG['nmMeans']
        lHStr, iS = GF.getLHStrISpl(self.lOD[iB].dITp, 'lHdrAIC')
        dfrB = TF.fromCSV(self.lDITp[iB], self.lOD[iB].dOIn, sMn = sMn,
                          sMod = 'Tr', sepD = self.dITp['cSep'], iCol = 0,
                          lHStr = lHStr, iSp = iS)[0]
        # move first string bit of column names in case of DGX
        GF.prepPdDfr(dfrB, self.idO, ['MetDGX', 'PhoDGX', 'CombDGX'])
        for i, OD in enumerate(self.lOD):
            lHStr, iS = GF.getLHStrISpl(self.lOD[i].dITp, 'lHdrAIC')
            dfrC = TF.fromCSV(self.lDITp[i], OD.dOIn, sMn = sMn, sMod = 'Tr',
                              sepD = self.dITp['cSep'], iCol = 0,
                              lHStr = lHStr, iSp = iS)[0]
            GF.prepPdDfr(dfrC, self.idO, ['MetDGX', 'PhoDGX', 'CombDGX'])
            # arrR = TF.doDivCompOp(self.dITp, OD.dITp, dfrC.T, dfrB.T)
            arrR = TF.doDivCompOp(self.dITp, dfrC.T, dfrB.T)
            cDfrCmp.iloc[:, j:(j + arrR.shape[1])] = arrR
            j += arrR.shape[1]
        if j != cDfrCmp.shape[1]:
            print('ERROR:', cDfrCmp.shape[1], 'columns unequal j =', j, '.')
            assert False
        TF.saveTrData(self.dITp, self.dOIn, cDfrCmp, self.dfrAIC,
                      sSt = self.dITp['sComp'], sMn = sMn)
    
    def prepXDfrs(self, lANm, useI, lPMv, pIn, sSep = '_'):
        for lDfr in [self.lODfr, self.lODfrMn, self.lODfrSD]:
            TF.moveStrRC(lDfr, lAllNm = lANm, useIdx = useI, lPosMove = lPMv,
                         posIns = pIn, sSep = sSep)
        lStStr = TF.moveStrRC(self.lODfrN, lAllNm = lANm, useIdx = useI,
                              lPosMove = lPMv, posIns = pIn, sSep = sSep)
        GF.insStrLDfr(self.lODfrAIC, [s + sSep for s in lStStr], useC = True)
    
    def concatToXDfrs(self, tpX, cAx = 0):
        self.cDfr = GF.concPdDfrS(self.lODfr, concAx = cAx)
        self.cDfrMn = GF.concPdDfrS(self.lODfrMn, concAx = cAx)
        self.cDfrSD = GF.concPdDfrS(self.lODfrSD, concAx = cAx)
        self.dfrN = GF.concPdDfrS(self.lODfrN, concAx = cAx)
        if tpX == 'S':
            self.dfrAIC = GF.concPdDfrS(self.lODfrAIC, concAx = cAx)
        elif tpX == 'F':
            self.dfrAIC = GF.get1stFromLDfr(self.lODfrAIC)
   
    def mergeTpXS(self, sSep = '_'):
        self.attrDfr = GF.concPdDfrS(self.lAttrDfr, concAx = 0, ignIdx = True)
        self.attrDfrMn = GF.concPdDfrS(self.lAttrDfrMn, concAx = 0,
                                       ignIdx = True)
        if self.cXFt == self.dIG['nmXFtC']:
            lANm, useI, lPMv, pIn = self.dIG['lNmGT'], True, [1, 2], 0
        else:
            lANm, useI, lPMv, pIn = self.dIG['lNmFt'], None, None, None
        if useI is not None and lPMv is not None and pIn is not None:
            self.prepXDfrs(lANm, useI = useI, lPMv = lPMv, pIn = pIn,
                           sSep = sSep)
        self.concatToXDfrs(tpX = 'S', cAx = 1)
        self.lFtHd = GF.extractFeatureHdr(self.cDfr)
        self.lFtHdMn = GF.extractFeatureHdr(self.cDfrMn)
        self.compareObj()

    def mergeTpXF(self, sSep = '_'):
        if self.cXFt == self.dIG['nmXFtC']:
            lGTFt = [OD.cGT for OD in self.lOD]
        else:
            lGTFt, lANm = [OD.cFt for OD in self.lOD], self.dIG['lNmFt']
            self.prepXDfrs(lANm, useI = False, lPMv = [0, 1], pIn = 0,
                           sSep = sSep)
        self.concatToXDfrs(tpX = 'F', cAx = 0)
        self.attrDfr = TF.getAttrDfrTpXF(self.dIG, self.lODfr, self.lAttrDfr,
                                         lGTFt)
        self.attrDfrMn = TF.getAttrDfrTpXF(self.dIG, self.lODfrMn,
                                           self.lAttrDfrMn, lGTFt)
        self.lFtHd = GF.extractFeatureHdr(self.cDfr)
        self.lFtHdMn = GF.extractFeatureHdr(self.cDfrMn)

    def extractDataXFt(self):
        if self.cXFt == self.dIG['nmXFtC']:
            dID, lPMove, lPSel = self.dIG['dIDFt'], [0, 1], [0, 1]
        if self.cXFt == self.dIG['nmXFtG']:
            TF.moveStrRC(self.lODfr, lAllNm = self.dIG['lNmFt'],
                         useIdx = False, lPosMove = [0, 1], posIns = 0)
            dID, lPMove, lPSel = self.dIG['dIDGT'], [1, 2], [0, 1]
        self.cDfr = GF.concPdDfrS(self.lODfr, concAx = 0)
        cDfrXFt = TF.getDfrXFt(self.dIG, self.cDfr, self.cXFt, useIdx = True,
                               lPosMove = lPMove, posIns = 0)
        return GF.getDDfrFtSel(cDfrXFt, dID, lPosSel = lPSel)

    def mergeO(self):
        assert self.cXFt in [self.dIG['nmXFtC'], self.dIG['nmXFtG']]
        if self.tpX == self.dIG['nmTpXS'] and not self.XtrDt:
            self.mergeTpXS()
        elif self.tpX == self.dIG['nmTpXF'] and not self.XtrDt:
            self.mergeTpXF()
        if self.XtrDt:
            # rearrange data to obtain a dict. of the other feature type data
            return self.extractDataXFt()

    def saveData(self):
        lCHdrAIC, sNmMn = self.dOIn['lCHdrAIC'], self.dIG['nmMeans']
        cXDfrAIC = GF.getExtDfrSelAI(self.dfrAIC, self.cDfr, lCHdrAIC).T
        cXDfrMnAIC = GF.getExtDfrSelAI(self.dfrAIC, self.cDfrMn, lCHdrAIC).T
        self.pFAllD = TF.saveData(self.dITp, self.dOIn, cXDfrAIC)
        self.pFMn = TF.saveData(self.dITp, self.dOIn, cXDfrMnAIC, sMn = sNmMn)
    
    def selDtXFt(self, iTp, cK):
        if self.idO == 'MetDCX':
            cDG = MetDataG(self.inpDG, iTp, cFt = cK, nmFAdd = self.nmFE)
        elif self.idO == 'PhoDCX':
            cDG = PhoDataG(self.inpDG, iTp, cFt = cK, nmFAdd = self.nmFE)
        elif self.idO == 'MetDGX':
            cDG = MetDataC(self.inpDG, iTp, cGT = cK, nmFAdd = self.nmFE)
        elif self.idO == 'PhoDGX':
            cDG = PhoDataC(self.inpDG, iTp, cGT = cK, nmFAdd = self.nmFE)
        else:
            if self.cXFt == self.dIG['nmXFtC']:
                cDG = ExpDataG(self.inpDG, iTp, cFt = cK, nmFAdd = self.nmFE)
            elif self.cXFt == self.dIG['nmXFtG']:
                cDG = ExpDataC(self.inpDG, iTp, cGT = cK, nmFAdd = self.nmFE)
            else:
                print('ERROR: Case idO ==', self.idO, 'and cXFt ==', self.cXFt,
                      'not implemented.')
                assert False
        return cDG

    def iniDtXFt(self, dDfr):
        dDtXFt = None
        if self.idO in ['ExpDX', 'MetDCX', 'PhoDCX', 'MetDGX', 'PhoDGX']:
            iTp, dDtXFt = GF.get1stFromL(self.lITp), {}
            for cK, cV in dDfr.items():
                cDtXFt = self.selDtXFt(iTp, cK)
                cDtXFt.calcMeanSD(cV)
                dDtXFt[cK] = cDtXFt    # e.g.: DR: MetDataG of DR over all GT
        return dDtXFt

    def combineData(self):
        dDfr = self.mergeO()
        self.updateDOIn()
        if not self.XtrDt:
            self.saveData()
            self.calcTransData(dfrAI = self.dfrAIC)
            self.calcDevData(dfrAI = self.dfrAIC)
            self.updateDOIn()
        if dDfr is not None and self.XtrDt:
            return self.iniDtXFt(dDfr)

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class MetDataCX(ExpDataX):
    def __init__(self, inpDat, iTp, lODat, typeX, extrDt = False):
        super().__init__(inpDat, iTp, lODat, typeX = typeX, extrDt = extrDt)
        self.idO = GC.S_MET_D + GC.NM_XFEAT_C + GC.NM_X
        self.descO = 'Metabolite data (feature type "C") extended'
        self.updateDOInBasicInfo()
        self.updateSFName(lD = [OD.dFNmComp for OD in self.lOD])
        print('Initiated "MetDataCX" base object.')

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class MetDataGX(ExpDataX):
    def __init__(self, inpDat, iTp, lODat, typeX, extrDt = False):
        super().__init__(inpDat, iTp, lODat, typeX = typeX, extrDt = extrDt)
        self.idO = GC.S_MET_D + GC.NM_XFEAT_G + GC.NM_X
        self.descO = 'Metabolite data (feature type "G") extended'
        self.updateDOInBasicInfo()
        self.updateSFName(lD = [OD.dFNmComp for OD in self.lOD])
        print('Initiated "MetDataGX" base object.')

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class PhoDataCX(ExpDataX):
    def __init__(self, inpDat, iTp, lODat, typeX, extrDt = False):
        super().__init__(inpDat, iTp, lODat, typeX = typeX, extrDt = extrDt)
        self.idO = GC.S_PHO_D + GC.NM_XFEAT_C + GC.NM_X
        self.descO = 'Phosphopeptide data (feature type "C") extended'
        self.updateDOInBasicInfo()
        self.updateSFName(lD = [OD.dFNmComp for OD in self.lOD])
        print('Initiated "PhoDataCX" base object.')

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class PhoDataGX(ExpDataX):
    def __init__(self, inpDat, iTp, lODat, typeX, extrDt = False):
        super().__init__(inpDat, iTp, lODat, typeX = typeX, extrDt = extrDt)
        self.idO = GC.S_PHO_D + GC.NM_XFEAT_G + GC.NM_X
        self.descO = 'Phosphopeptide data (feature type "G") extended'
        self.updateDOInBasicInfo()
        self.updateSFName(lD = [OD.dFNmComp for OD in self.lOD])
        print('Initiated "PhoDataGX" base object.')

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class CombDataCX(ExpDataX):
    def __init__(self, inpDat, iTp, lODat, typeX, extrDt = False):
        super().__init__(inpDat, iTp, lODat, typeX = typeX, extrDt = extrDt)
        self.idO = 'CombDCX'
        self.descO = 'Combined data (feature type "C") extended'
        self.updateDOInBasicInfo()
        self.updateSFName(lD = [OD.dFNmComp for OD in self.lOD])
        print('Initiated "CombDataCX" base object.')

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class CombDataGX(ExpDataX):
    def __init__(self, inpDat, iTp, lODat, typeX, extrDt = False):
        super().__init__(inpDat, iTp, lODat, typeX = typeX, extrDt = extrDt)
        self.idO = 'CombDGX'
        self.descO = 'Combined data (feature type "G") extended'
        self.updateDOInBasicInfo()
        self.updateSFName(lD = [OD.dFNmComp for OD in self.lOD])
        print('Initiated "CombDataGX" base object.')

###############################################################################

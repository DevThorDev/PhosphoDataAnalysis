# -*- coding: utf-8 -*-
###############################################################################
# --- O_82__Clustering.py -----------------------------------------------------
###############################################################################
import os

from sklearn import cluster

import Core.C_00__GenConstants as GC
import Core.F_00__GenFunctions as GF
import Core.F_02__PltFunctions as PF
import Core.F_03__OTpFunctions as TF

from Core.O_00__DataBaseClass import DataBaseClass
from Core.O_51__ClusterResData import ClusterResData

class Clustering(DataBaseClass):
    def __init__(self, inpDat, iTp, ODat, nmFAdd = ''):
        super().__init__(inpDat, iTp, nmFAdd = nmFAdd)
        self.idO = GC.S_CLUST
        self.descO = 'Clustering'
        self.OD = ODat
        self.tpX = self.OD.tpX
        self.cXFt = self.OD.cXFt
        self.ODfr = self.OD.cDfrMn        # clustering only useful for Mn
        self.ODfrTr = self.OD.cDfrTr
        self.ODfrTrT = self.OD.cDfrTrT
        self.getODfrAddI()
        self.dOITp = self.OD.dITp
        self.cOITp = self.dOITp['iTp']
        self.cGT = self.OD.cGT
        self.cFt = self.OD.cFt
        self.iniDicts()
        self.updateDOIn()
        self.updateSFName(lD = [self.OD.dFNmComp],
                          dS = {'No': [self.OD.sFNm],
                                'Tr': [self.OD.sFNmTr],
                                'Dv': [self.OD.sFNmDv]})
        self.getPResF()
        self.updateDOIn()
        print('Initiated "Clustering" base object.')
    
    def printObjInfo(self):
        super().printObjInfo()
        print('Path of general result file:', self.pRF)
        print('Current path of cluster result file:', self.pRFF)
    
    def getODfrAddI(self):
        self.ODfrAIC = None
        if 'lHdrAIC' in self.OD.dITp:
            lSel = self.OD.dITp['lHdrAIC']
            if set(lSel).issubset(set(self.OD.addIDfr.index)):
                self.ODfrAIC = self.OD.addIDfr.loc[lSel, :]
        if hasattr(self.OD, 'dfrAIC'):  # overwrite in this case (ExpDX)
            self.ODfrAIC = self.OD.dfrAIC
    
    def iniDicts(self):
        self.dClAlgDfr = {}
        self.dClAlgDfr4Corr = {}
        self.dClResData = {}

    def getPResF(self):
        self.pRF, self.pRFF = TF.get2xPTrF(self.dITp, self.dOIn,
                                           self.dIG['nmMeans'])
    
    def updateDOIn(self):
        self.updateDOInBasicInfo()
        lA = ['tpX', 'cXFt', 'ODfr', 'ODfrTr', 'ODfrTrT', 'ODfrAIC', 'dOITp',
              'cOITp', 'cGT', 'cFt', 'nmFE', 'dClAlgDfr', 'dClAlgDfr4Corr',
              'dClResData', 'pRF', 'pRFF']
        for cA in lA:
            if hasattr(self, cA):
                self.dOIn[cA] = getattr(self, cA)

# --- Methods for special clustering objects ----------------------------------
    def selectPRF(self, isRFF = False):
        pRF = self.pRF
        if isRFF:
            pRF = self.pRFF
        return pRF
    
    def modifyPRF(self, nFAdd = '', doAppend = False):
        pRF = self.pRF
        if doAppend:
            pRF = self.pRFF
        self.pRFF = GF.modFName(pRF, sPost = nFAdd)
    
    def plotClCent(self):
        lSpc, dDfrPlot = self.dIG['dNmPltClCent'][self.cOITp], {}
        lNmCSel = GF.getLSubIdxOrCol(self.ODfrTr, lSpec = lSpc, isIdx = False)
        addDfr = self.ODfrTr.loc[self.OD.lSRowTrD, lNmCSel]
        for cK, cClDfr in self.dClAlgDfr4Corr.items():
            dDfrPlot[cK] = GF.concPdDfrS([cClDfr, addDfr], concAx = 1)
        PF.pltClCent(self.dITp, self.dITp[self.dITp['nmPlt_ClCent']],
                     self.OD.sTrDt, dDfrPlot, self.selectPRF())

    def prepClRes4Corr(self, inpDat):
        lSRTrD = self.OD.lSRowTrD
        # create the DataFrame
        for cN in self.dOITp['dNCl4Corr'][(self.cXFt, self.tpX)]:
            lNmC = ['Cl_' + str(cE) for cE in self.dClAlgDfr[cN]['ClLabel']]
            cDfr = GF.iniPdDfr(lSNmC = lNmC, lSNmR = lSRTrD)
            for cNmR in lSRTrD:
                cDfr.loc[cNmR, :] = list(self.dClAlgDfr[cN].loc[:, cNmR])
            self.dClAlgDfr4Corr[cN] = cDfr
            # complement the attributes for BinaryOps methods
            self.dClResData[cN] = ClusterResData(inpDat, self.dIG['iClRD'],
                                                 self.OD, self.idO, cN,
                                                 cDfr, lSRTrD, self.nmFE)

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class KMeansClustering(Clustering):
    def __init__(self, inpDat, iTp, ODat, nmFAdd = ''):
        super().__init__(inpDat, iTp, ODat, nmFAdd)
        self.idO = GC.S_KM_CLUST
        self.descO = 'K-Means Clustering'
        self.updateDOInBasicInfo()
        self.updateSFName(lD = [self.OD.dFNmComp],
                          dS = {'No': [self.OD.sFNm],
                                'Tr': [self.OD.sFNmTr],
                                'Dv': [self.OD.sFNmDv]})
        self.getPResF()
        self.doKMeansCl(inpDat)
        self.plotClCent()
        self.updateDOIn()
        print('Initiated "KMeansClustering" base object.')

    def storeClRes(self, cSep, lNCl, dRes):
        t = (self.dITp['initMeth'], self.dITp['nInit'],
             self.dITp['maxIter'], self.dITp['relTol'])
        arrDat = self.ODfrTrT.values
        for i, nCl in enumerate(lNCl):
            cKMns = cluster.KMeans(n_clusters = nCl, init = t[0],
                                   n_init = t[1], max_iter = t[2], tol = t[3])
            wtMn = self.OD.attrDfrMn.loc[:, self.nmWtC]
            cResCl = cKMns.fit(arrDat, sample_weight = wtMn)
            dRes['Inertia'][i] = cResCl.inertia_
            dRes['NIter'][i] = cResCl.n_iter_
            self.ODfrTrT.loc[:, 'Labels_' + str(nCl) + 'Cl'] = cResCl.labels_
            cClDfr = GF.iniPdDfr(lSNmR = range(nCl))
            cClDfr['ClLabel'] = range(nCl)
            lNmCClC = [self.OD.lSRowTrD[j] for j in range(arrDat.shape[1])]
            for j, cNm in enumerate(lNmCClC):
                cClDfr[cNm] = cResCl.cluster_centers_[:, j]
            cClDfr['Inertia'] = [cResCl.inertia_]*nCl
            cClDfr['NIter'] = [cResCl.n_iter_]*nCl
            self.dClAlgDfr[nCl] = cClDfr
            if self.dITp['saveClInfo']:
                self.modifyPRF('__' + str(nCl) + 'Cl')
                cClDfr.to_csv(self.pRFF, sep = cSep)
        cDfrSt = self.ODfrTrT
        if self.ODfrAIC is not None:
            cDfrSt = GF.getExtDfrSelAI(self.ODfrAIC.T, self.ODfrTrT,
                                       self.OD.dITp['lHdrAIC'], isT = True)
        cDfrSt.to_csv(self.pRF, sep = cSep)

    def loadClRes(self, cSep, lNCl, dRes):
        lHStr, iS = GF.getLHStrISpl(self.OD.dITp, 'lHdrAIC')
        self.ODfrTrT = GF.readCSV(self.pRF, sepD = cSep, iCol = 0,
                                  lHStr = lHStr, iSp = iS)[0]
        self.ODfrTr = self.ODfrTrT.T
        for i, nCl in enumerate(lNCl):
            self.modifyPRF('__' + str(nCl) + 'Cl')
            if os.path.isfile(self.pRFF):
                cClDfr = GF.readCSV(self.pRFF, sepD = cSep, iCol = 0)
                if cClDfr.shape[0] > 0:
                    dRes['Inertia'][i] = cClDfr.loc[0, 'Inertia']
                    dRes['NIter'][i] = cClDfr.loc[0, 'NIter']
                self.dClAlgDfr[nCl] = cClDfr

    def doKMeansCl(self, inpDat):
        lNCl = self.dOITp['lNClusters']
        dRes = {'Inertia': [0]*len(lNCl), 'NIter': [0]*len(lNCl), 'lNCl': lNCl}
        if not os.path.isfile(self.pRF):    # do K-Means clustering
            self.storeClRes(self.dITp['cSep'], lNCl, dRes)
        else:                               # load respective file
            self.loadClRes(self.dITp['cSep'], lNCl, dRes)
        # prepare dictionary of DataFrames for correlation calculation
        self.prepClRes4Corr(inpDat)
        for cT in [('1DS_1', 'Inertia'), ('1DS_2', 'NIter')]:
            PF.plt1DDatS(self.dITp, self.dITp[self.dITp['nmPlt_' + cT[0]]],
                         GF.toPdDfr(dRes), self.selectPRF(), cT, 'lNCl')

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class AgglomerativeClustering(Clustering):
    def __init__(self, inpDat, iTp, ODat, nmFAdd = ''):
        super().__init__(inpDat, iTp, ODat, nmFAdd)
        self.idO = GC.S_AG_CLUST
        self.descO = 'Agglomerative Clustering'
        self.updateDOInBasicInfo()
        self.updateSFName(lD = [self.OD.dFNmComp],
                          dS = {'No': [self.OD.sFNm],
                                'Tr': [self.OD.sFNmTr],
                                'Dv': [self.OD.sFNmDv]})
        self.getPResF()
        self.doAggloCl(inpDat)
        self.plotClCent()
        self.updateDOIn()
        print('Initiated "AgglomerativeClustering" base object.')

    def doAggloCl(self, inpDat):
        pass

###############################################################################

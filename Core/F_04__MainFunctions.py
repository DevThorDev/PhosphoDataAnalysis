# -*- coding: utf-8 -*-
###############################################################################
# --- F_04__MainFunctions.py --------------------------------------------------
###############################################################################
import Core.F_00__GenFunctions as GF

from Core.O_01__ExpData import ExpData
from Core.O_11__MetData import MetData
from Core.O_21__PhoData import PhoData
from Core.O_31__CombData import CombDataC, CombDataG
from Core.O_41__ExpDataX import (MetDataCX, PhoDataCX, CombDataCX,
                                 MetDataGX, PhoDataGX, CombDataGX)
from Core.O_81__BinaryOps import BinaryOps
from Core.O_82__Clustering import KMeansClustering, AgglomerativeClustering
from Core.O_83__OverRep import OverRep

# --- Functions (initialisation) ----------------------------------------------
def initAnalysis(inpDG):
    lI, XFtC, XFtG = ['']*6, inpDG.dI['nmXFtC'], inpDG.dI['nmXFtG']
    XFtCTf, XFtGTf = XFtC + inpDG.dI['sTransf'], XFtG + inpDG.dI['sTransf']
    lSUsDC, lSMgDC = inpDG.dI['lSUsedDLCTf'], inpDG.dI['lSMergDLC']
    lSCRDC, lSCcD = inpDG.dI['lSClRDLC'], inpDG.dI['lSCalcDL']
    lSUsDG = GF.modStrEndLStr(lSUsDC, XFtCTf, XFtGTf)
    lSMgDG = GF.modStrEndLStr(lSMgDC, XFtC, XFtG)
    lSCRDG = GF.modStrEndLStr(lSCRDC, XFtC, XFtG)
    lI[:2] = inpDG.dI['dIInp1'][lSUsDC[0]], inpDG.dI['dIInp1'][lSUsDC[1]]
    lI[2:4] = inpDG.dI['dIInp1'][lSMgDC[0]], inpDG.dI['dIInp1'][lSMgDC[1]]
    lI[4] = inpDG.dI['dIInp1'][lSCRDC[0]]
    lI[5:] = inpDG.dI['dIInp1'][lSCcD[0]], inpDG.dI['dIInp1'][lSCcD[1]]
    return lSUsDC, lSMgDC, lSCRDC, lSUsDG, lSMgDG, lSCRDG, [], [], [], {}, lI

# --- Functions (load data) ---------------------------------------------------
def loadData(inpDG, iTp, cGT, lDt, sD, printInf = True):
    XFtC, tQryNmFE = inpDG.dI['nmXFtC'], inpDG.dI['dQryNmFE'][iTp]
    print('='*8, 'Loading ' + sD + ', and calculating means and SDs...')
    if sD in inpDG.dI['lSBasDLTf']:    # only these data can be loaded
        if sD == inpDG.dI['lSUsedDLTf'][0]:
            cDt = MetData(inpDG, iTp, cXFrt = XFtC, cGT = cGT,
                          nmFAdd = tQryNmFE[1], cQuery = tQryNmFE[0])
        elif sD == inpDG.dI['lSUsedDLTf'][1]:
            cDt = PhoData(inpDG, iTp, cXFrt = XFtC, cGT = cGT,
                          nmFAdd = tQryNmFE[1], cQuery = tQryNmFE[0])
        else:
            cDt = ExpData(inpDG, iTp, cXFrt = XFtC, cGT = cGT,
                          nmFAdd = tQryNmFE[1], cQuery = tQryNmFE[0])
        cDt.calcMeanSD()
        lDt.append(cDt)
        if printInf:
            cDt.printObjInfo()
        return cDt
    return None

# --- Functions (combine data) ------------------------------------------------
def combineBasicData(inpDG, lOD, lDt, XFtr, lSD = ['', ''],
                     printInf = True):
    assert len(lSD) >= 2 and XFtr in [inpDG.dI['nmXFtC'], inpDG.dI['nmXFtG']]
    print('='*8, 'Combining ' + lSD[0] + ' and ' + lSD[1] + '...')
    if XFtr == inpDG.dI['nmXFtC']:
        cDt = CombDataC(inpDG, inpDG.dI['iCombD'], lOD)
    elif XFtr == inpDG.dI['nmXFtG']:
        cDt = CombDataG(inpDG, inpDG.dI['iCombD'], lOD)
    cDt.combineData()
    lDt.append(cDt)
    if printInf:
        cDt.printObjInfo()
    return cDt

def combineBasicDataC(inpDG, lOD, lDt, lSD = ['', ''], printInf = True):
    XFtC = inpDG.dI['nmXFtC']
    return combineBasicData(inpDG, lOD, lDt, XFtr = XFtC, lSD = lSD,
                            printInf = printInf)

def combineBasicDataG(inpDG, lOD, lDt, lSD = ['', ''], printInf = True):
    XFtG = inpDG.dI['nmXFtG']
    return combineBasicData(inpDG, lOD, lDt, XFtr = XFtG, lSD = lSD,
                            printInf = printInf)

# --- Functions (correlations) ------------------------------------------------
def calcCorrelations(inpDG, lOD, lSD = ['', ''], doAllD = True, doMeans = True,
                     printInf = True):
    print('='*8, 'Calculating correlations between ' + lSD[0] + ' and ' +
          lSD[1] + '...')
    cBO = BinaryOps(inpDG, inpDG.dI['iBinOp'], lOD)
    cBO.doCorrDev(doAllD = doAllD, doMeans = doMeans)
    cOR = OverRep(inpDG, inpDG.dI['iOvRep'], cBO)
    if doAllD:
        cOR.calcProfiles(tpPr = 'BO', useMn = False)
    if doMeans:
        cOR.calcProfiles(tpPr = 'BO', useMn = True)
    if printInf:
        cBO.printObjInfo()

# --- Functions (clustering) --------------------------------------------------
def calcKMeansCluster(inpDG, k, cDt, dClst, printInf = True):
    cClst = KMeansClustering(inpDG, inpDG.dI['iClust'], cDt)
    dClst[(cClst.OD.idO, cClst.cOITp, k)] = cClst
    if printInf:
        cClst.printObjInfo()

# --- Functions (clustering and correlations of resulting clusters) -----------
def combToDtX(inpDG, lOD, sO, XFtr, typeX, XtrDt = False):
    iTpX = inpDG.dI['iExpDX']
    if XFtr == inpDG.dI['nmXFtC'] and sO == inpDG.dI['lSUsedDLTfCombDL'][0]:
        cDtX = MetDataCX(inpDG, iTpX, lOD, typeX = typeX, extrDt = XtrDt)
    elif XFtr == inpDG.dI['nmXFtC'] and sO == inpDG.dI['lSUsedDLTfCombDL'][1]:
        cDtX = PhoDataCX(inpDG, iTpX, lOD, typeX = typeX, extrDt = XtrDt)
    elif XFtr == inpDG.dI['nmXFtG'] and sO == inpDG.dI['lSUsedDLTfCombDL'][2]:
        cDtX = MetDataGX(inpDG, iTpX, lOD, typeX = typeX, extrDt = XtrDt)
    elif XFtr == inpDG.dI['nmXFtG'] and sO == inpDG.dI['lSUsedDLTfCombDL'][3]:
        cDtX = PhoDataGX(inpDG, iTpX, lOD, typeX = typeX, extrDt = XtrDt)
    elif XFtr == inpDG.dI['nmXFtC'] and sO == inpDG.dI['lSUsedDLTfCombDL'][4]:
        cDtX = CombDataCX(inpDG, iTpX, lOD, typeX = typeX, extrDt = XtrDt)
    elif XFtr == inpDG.dI['nmXFtG'] and sO == inpDG.dI['lSUsedDLTfCombDL'][5]:
        cDtX = CombDataGX(inpDG, iTpX, lOD, typeX = typeX, extrDt = XtrDt)
    else:
        print('ERROR: Case XFtr ==', XFtr, 'and sO ==', sO, 'not implemented.')
        assert False
    return [cDtX, cDtX.combineData()]

def combToDtCX(inpDG, lOD, sO, typeX, XtrDt = False):
    XFtC = inpDG.dI['nmXFtC']
    return combToDtX(inpDG, lOD, sO, XFtC, typeX = typeX, XtrDt = XtrDt)

def combToDtGX(inpDG, lOD, sO, typeX, XtrDt = False):
    XFtG = inpDG.dI['nmXFtG']
    return combToDtX(inpDG, lOD, sO, XFtG, typeX = typeX, XtrDt = XtrDt)

def corrClRDt(inpDG, lClDt, printInf = True):
    print('Calculating correlations between the clusters...')
    cBO_ClustDt = BinaryOps(inpDG, inpDG.dI['iBinOp'], lClDt)
    cBO_ClustDt.doCorrDev(doAllD = False)
    if printInf:
        cBO_ClustDt.printObjInfo()

def clustCorrelClRDt(inpDG, lDX, printInf = True):
    print('Clustering...')
    lClst = []
    for cDX in lDX:
        lClst.append(KMeansClustering(inpDG, inpDG.dI['iClust'], cDX))
    for nCl1 in lDX[0].dITp['dNCl4Corr'][(lDX[0].cXFt, lDX[0].tpX)]:
        cClRD1 = lClst[0].dClResData[nCl1]
        for nCl2 in lDX[1].dITp['dNCl4Corr'][(lDX[1].cXFt, lDX[1].tpX)]:
            cClRD2 = lClst[1].dClResData[nCl2]
            corrClRDt(inpDG, [cClRD1, cClRD2], printInf = printInf)

def calcClustersCorrelX(inpDG, lOD, sO, XFtr, typeX, XtrDt = False,
                        printInf = True):
    assert sO in inpDG.dI['lSUsedDLTfCombDL']
    print('='*8, 'Combining data over all genotypes for ' + sO +
          ' using type ' + typeX + '...')
    cDX, dDt = combToDtX(inpDG, lOD, sO = sO, XFtr = XFtr, typeX = typeX,
                         XtrDt = XtrDt)
    if not XtrDt:
        clustCorrelClRDt(inpDG, [cDX, cDX], printInf = printInf)
    return dDt

def calcClustersCorrel2X(inpDG, llOD, lSO, XFtr, typeX, XtrDt = False,
                         printInf = True):
    assert len(lSO) == len(llOD)
    if len(lSO) >= 1:
        lDX, lDDt = [], []
        if len(lSO) == 1:
            d = calcClustersCorrelX(inpDG, llOD[0], lSO[0], typeX = typeX,
                                    XtrDt = XtrDt, printInf = printInf)
            lDDt.append(d)
        else:
            assert (lSO[0] in inpDG.dI['lSUsedDLTfCombDL'] and
                    lSO[1] in inpDG.dI['lSUsedDLTfCombDL'])
            print('='*8, 'Combining data over all genotypes for ' + lSO[0] +
                  ' and ' + lSO[1] + ' using type ' + typeX + '...')
            for i, sO in enumerate(lSO):
                cDX, dDt = combToDtX(inpDG, llOD[i], sO = sO, XFtr = XFtr,
                                     typeX = typeX, XtrDt = XtrDt)
                lDX.append(cDX)
                lDDt.append(dDt)
            if not XtrDt:
                clustCorrelClRDt(inpDG, lDX, printInf = printInf)
        return lDDt

# --- Functions (obtaining data of one feature type from data of the other) ---
def extractXFtData(inpDG, lOD, sO, printInf = True):
    assert sO in inpDG.dI['lSUsedDLTfCombDL']
    print('='*8, 'Extracting data of the other feature type for ' + sO + '...')
    cXFt = GF.get1stFromL([OD.cXFt for OD in lOD])
    assert cXFt in [inpDG.dI['nmXFtC'], inpDG.dI['nmXFtG']]
    if cXFt == inpDG.dI['nmXFtC']:
        return combToDtCX(inpDG, lOD, sO = sO, typeX = inpDG.dI['nmTpXF'],
                          XtrDt = True)[1]
    else:
        return combToDtGX(inpDG, lOD, sO = sO, typeX = inpDG.dI['nmTpXF'],
                          XtrDt = True)[1]

def getXFtData(dMD, dPD, lMD, lPD):
    lID = []
    for idFt in dMD:
        if idFt in dPD:
            lID.append(idFt)
    for cI in lID:
        lMD.append(dMD[cI])
        lPD.append(dPD[cI])

# --- Function (full run for given XFeature) ----------------------------------
def calcCurXFt(inDG, k, sXFt, doXFt, t, cMD, cPD, lMD, lPD, lCD, dClR, sTime):
    # initialise
    assert sXFt in [inDG.dI['nmXFtC'], inDG.dI['nmXFtG']]
    lSUsDC, lSMgDC, lSCRDC, lSUsDG, lSMgDG, lSCRDG = t
    if sXFt == inDG.dI['nmXFtC']:
        lSUsD, lSMgD = lSUsDC, lSMgDC
    else:
        lSUsD, lSMgD = lSUsDG, lSMgDG
    
    # do over-representation analysis with PhoData
    cOR = OverRep(inDG, inDG.dI['iOvRep'], cPD)
    cOR.calcProfiles()
    
    # combine MetData and PhoData
    if doXFt and inDG.dI['doCombDt']:
        if sXFt == inDG.dI['nmXFtC']:
            cCD = combineBasicDataC(inDG, [cMD, cPD], lCD, lSUsD)
        else:
            cCD = combineBasicDataG(inDG, [cMD, cPD], lCD, lSUsD)
        print('='*20, 'Combined MetData and PhoData to CombData:', k, '='*20)
    
    # calculate the correlations between MetData and PhoData
    if inDG.dI['calcCorrMetPho'] and doXFt:
        calcCorrelations(inDG, [cMD, cPD], lSUsD)
        GF.showElapsedTime(sTime)
        print('='*20, 'Calculated MetData-PhoData correlations:', k, '='*20)

    # calculate the correlations within MetData, PhoData and CombData
    # WARNING - takes a few hours to complete
    if inDG.dI['calcCorrWithin'] and doXFt:
        for i, cDt in enumerate([cMD, cPD]):
            calcCorrelations(inDG, [cDt, cDt], [lSUsD[i], lSUsD[i]])
        if inDG.dI['doCombDt']:
            calcCorrelations(inDG, [cCD, cCD], [lSMgD[0], lSMgD[0]])
        GF.showElapsedTime(sTime)
        print('='*20, 'Calculated correlations within Data:', k, '='*20)
    
    # calculate the K-Means cluster within MetData, PhoData and CombData
    if ((inDG.dI['calcKMeansClBasicDt'] or inDG.dI['calcCorrClDt']) and doXFt):
        for cDt in [cMD, cPD]:
            calcKMeansCluster(inDG, k, cDt, dClR)
        if inDG.dI['doCombDt']:
            calcKMeansCluster(inDG, k, cCD, dClR)
        GF.showElapsedTime(sTime)
        print('='*20, 'Calculated K-Means cluster of Data:', k, '='*20)

    # calculate the correlations between clusters of MetData and PhoData
    if inDG.dI['calcCorrClDt'] and doXFt:
        for nClMD in cMD.dITp['dNCl4Corr'][(cMD.cXFt, cMD.tpX)]:
            MDR = dClR[(cMD.idO, cMD.dITp['iTp'], k)].dClResData[nClMD]
            for nClPD in cPD.dITp['dNCl4Corr'][(cPD.cXFt, cPD.tpX)]:
                PDR = dClR[(cPD.idO, cPD.dITp['iTp'], k)].dClResData[nClPD]
                lS = [str(nClMD) + ' clusters for ' +  cMD.dITp['sNmSpec'],
                      str(nClPD) + ' clusters for ' +  cPD.dITp['sNmSpec']]
                calcCorrelations(inDG, [MDR, PDR], lS, doAllD = False)
        GF.showElapsedTime(sTime)
        print('='*20, 'Calculated MetData-PhoData cluster correl.:', k, '='*20)

    # calculate the correlations between clusters of CombData
    if inDG.dI['calcCorrClDt'] and inDG.dI['doCombDt'] and doXFt:
        for nCl1 in cCD.dITp['dNCl4Corr'][(cCD.cXFt, cCD.tpX)]:
            CD1 = dClR[(cCD.idO, cCD.dITp['iTp'], k)].dClResData[nCl1]
            for nCl2 in cCD.dITp['dNCl4Corr'][(cCD.cXFt, cCD.tpX)]:
                CD2 = dClR[(cCD.idO, cCD.dITp['iTp'], k)].dClResData[nCl2]
                lS = [str(nCl1) + ' clusters for ' +  cCD.dITp['sNmSpec'],
                      str(nCl2) + ' clusters for ' +  cCD.dITp['sNmSpec']]
                calcCorrelations(inDG, [CD1, CD2], lS, doAllD = False)
        GF.showElapsedTime(sTime)
        print('='*20, 'Calculated CombData cluster correl.:', k, '='*20)

# -----------------------------------------------------------------------------
def calcClustCorrClust(inDG, sXFt, doXFt, t, lMD, lPD, lCD, sTime):
    # calculate clusters and correlations thereof of the all-genotype data
    assert sXFt in [inDG.dI['nmXFtC'], inDG.dI['nmXFtG']]
    lSUsDC, lSMgDC, lSCRDC, lSUsDG, lSMgDG, lSCRDG = t
    if sXFt == inDG.dI['nmXFtC']:
        lSUsD, lSMgD = lSUsDC, lSMgDC
    else:
        lSUsD, lSMgD = lSUsDG, lSMgDG
    if inDG.dI['calcClustersCorrelX'] and doXFt:
        for tpX in inDG.dI['dInvXFtTpX'][sXFt]:
            # MetData / MetData
            calcClustersCorrelX(inDG, lMD, sO = lSUsD[0],
                                XFtr = sXFt, typeX = tpX, printInf = True)
            # PhoData / PhoData
            calcClustersCorrelX(inDG, lPD, sO = lSUsD[1],
                                XFtr = sXFt, typeX = tpX, printInf = True)
            # MetData / PhoData
            calcClustersCorrel2X(inDG, [lMD, lPD], lSO = lSUsD, XFtr = sXFt,
                                 typeX = tpX, printInf = True)
            # CombData / CombData
            if inDG.dI['doCombDt']:
                calcClustersCorrelX(inDG, lCD, sO = lSMgD[0], XFtr = sXFt,
                                    typeX = tpX, printInf = True)
        GF.showElapsedTime(sTime)
        print('='*20, 'Calculated clusters and correlations (DataX).', '='*20)

# -----------------------------------------------------------------------------
def extractDataXFt(inDG, doXFt, lSUsD, lMD, lPD, sTime):
    # extract data of the other feature type from the extended data
    if doXFt:
        # MetData / MetData
        dMD = extractXFtData(inDG, lMD, sO = lSUsD[0], printInf = True)
        # PhoData / PhoData
        dPD = extractXFtData(inDG, lPD, sO = lSUsD[1], printInf = True)
        GF.showElapsedTime(sTime)
        print('='*20, 'Extracted other feature type data from DataX.', '='*20)
        return dMD, dPD
    else:
        return {}, {}
###############################################################################

# -*- coding: utf-8 -*-
###############################################################################
# --- F_03__OTpFunctions.py ---------------------------------------------------
###############################################################################
import os, copy

import numpy as np
import pandas as pd
import scipy.stats as stats

import Core.C_00__GenConstants as GC
import Core.F_00__GenFunctions as GF
import Core.F_02__PltFunctions as PF

# --- Functions (general) -----------------------------------------------------
def getPF(dITp, nmF, nmFExt = GC.NM_EXT_CSV):
    return GF.joinToPath(dITp['pRelResF'], nmF + '.' + nmFExt)

def getNmF(dOIn, sF = 'sFNm', sSt = '', sMn = '', sEnd = '', sPtSt = '__',
           sPrMn = '__', sPrEnd = '__'):
    nmF = GF.addString(sSt, sPost = sPtSt) + dOIn[sF]
    nmF += GF.addString(sMn, sPre = sPrMn)
    nmF += GF.addString(sEnd, sPre = sPrEnd)
    return nmF

def getPFRes(dITp, dOIn, sSt = '', sMn = '', sEnd = '', sMod = 'No',
             nmFExt = GC.NM_EXT_CSV):
    if 'calcForTrD' in dITp:
        if dITp['calcForTrD']:    # do calculations for the transformed data
            sMod = 'Tr'
    if sMod == 'No' or dOIn['isClRD']:  # ClRD is always transformed
        nmF = getNmF(dOIn, sF = 'sFNm', sSt = sSt, sMn = sMn, sEnd = sEnd)
    else:   # transformed (sMod == 'Tr') and deviation (sMod == 'Dv') data
        nmF = getNmF(dOIn, sF = 'sFNm' + sMod, sSt = sSt, sMn = sMn,
                     sEnd = sEnd)
    return getPF(dITp, nmF, nmFExt)

def fromCSV(dITp, dOIn, sSt = '', sMn = '', sMod = 'No', sepD = ',',
            iCol = None, lHStr = [], iSp = None, splC = True,
            nmFExt = GC.NM_EXT_CSV):
    pF = getPFRes(dITp, dOIn, sSt, sMn = sMn, sMod = sMod, nmFExt = nmFExt)
    return GF.readCSV(pF, dITp['cSep'], iCol = iCol, lHStr = lHStr, iSp = iSp,
                      splC = splC)

def savePdDfr(dITp, pdDfr, nmF, nmFExt = GC.NM_EXT_CSV, saveIt = True):
    pF = getPF(dITp, nmF, nmFExt)
    if not os.path.isfile(pF) and saveIt:
        pdDfr.to_csv(pF, sep = dITp['cSep'])
    return pF

def saveData(dITp, dOIn, pdDfr, sSt = '', sMn = '', nmFExt = GC.NM_EXT_CSV):
    return savePdDfr(dITp, pdDfr, getNmF(dOIn, sSt = sSt, sMn = sMn), nmFExt)

def saveTrData(dITp, dOIn, pdDfr, pdDfrAI, sSt = '', sMn = '',
               nmFExt = GC.NM_EXT_CSV):
    nmF = getNmF(dOIn, sF = 'sFNmTr', sSt = sSt, sMn = sMn)
    pdDfrAIC = GF.getExtDfrSelAI(pdDfrAI, pdDfr, dOIn['lCHdrAIC']).T
    return savePdDfr(dITp, pdDfrAIC.dropna(how = 'all'), nmF, nmFExt)

def saveDvData(dITp, dOIn, pdDfr, pdDfrAI, sSt = '', sMn = '',
               nmFExt = GC.NM_EXT_CSV):
    nmF = getNmF(dOIn, sF = 'sFNmDv', sSt = sSt, sMn = sMn)
    pdDfrAIC = GF.getExtDfrSelAI(pdDfrAI, pdDfr, dOIn['lCHdrAIC']).T
    return savePdDfr(dITp, pdDfrAIC.dropna(how = 'all'), nmF, nmFExt)

def deepCopyLO(lO):
    return [copy.deepcopy(O) for O in lO]

def iniModDITp(dITp, lDITp, lK):
    for sK in lK:
        dITp[sK] = GF.get1stFromLDict(lDITp, sK)

# --- Functions (O_00__DataBaseClass) -----------------------------------------
def getDITp(dIG, iTp0, iTp, sSep = '_'):
    dITp = copy.deepcopy(dIG[iTp0])    # content of iTp = 0 input
    dITp.update(dIG[iTp])              # updated with iTp = iTp input
    for i, sP in enumerate(['pRelResF', 'pRelPltF']):   # update result paths
        nmDir = dIG['tNmPreDirOut'][i] + sSep + dIG['tIDResFig'][i] + sSep
        nmDir += ('0'*(2 - len(str(iTp))) + str(iTp) + sSep + dITp['sNmSpec'])
        dITp[sP] = GF.joinDirToPath(dITp[sP], nmDir)
    return dITp

def getIDGT(dIG, nmGT):
    return GF.getNm(dIG['dNmGT'], nmGT, dIG['nmAllGT'])

def getIDFt(dIG, nmFt):
    return GF.getNm(dIG['dNmFt'], nmFt, dIG['nmAllFt'])

def getNmGT(dIG, cGT):
    return GF.getNm(dIG['dIDGT'], cGT, dIG['nmAllGT'])

def getNmFt(dIG, cFt):
    return GF.getNm(dIG['dIDFt'], cFt, dIG['nmAllFt'])

def getLAttr(dIG, lVs, iCGT, iCFt):
    cGT, cFt, lVsRet = lVs[iCGT], lVs[iCFt], []
    nmGT = getNmGT(dIG, cGT)
    nmFt = getNmFt(dIG, cFt)
    for i, cV in enumerate(lVs):
        if i == iCGT:
            lVsRet += [cGT, nmGT]
        elif i == iCFt:
            lVsRet += [cFt, nmFt]
        else:
            lVsRet += [cV]
    return lVsRet

def fillTransData(pdDfrTr, pdDfr, transD, lSR, iB, lI):
    if transD in ['Abs', 'Rel', '2LQ']:     # otherwise, do nothing
        cDen = 1.
        if transD in ['Rel', '2LQ']:
            cDen = pdDfr.iloc[iB, :]
        for i, nmR in enumerate(lSR):
            if transD in ['2LQ']:
                cEn = pdDfr.iloc[lI[i], :]
                pdDfrTr.loc[nmR, :] = np.log2(cEn/cDen)
            else:
                cEn = pdDfr.iloc[lI[i], :] - pdDfr.iloc[iB, :]
                pdDfrTr.loc[nmR, :] = cEn/cDen

def getSFNm(dFNmCmp, dISPr = {}, lS = [], sSep = '_'):
    assert len(dFNmCmp) == 1
    cK, sFNm = list(dFNmCmp)[0], ''
    for i, sPr in dISPr.items():
        if i < len(cK):
            if cK[i] is not None:
                sFNm += GF.addString(cK[i], sPre = sPr)
    for cS in lS:
        sFNm += sSep + cS
    return sFNm

# --- Functions (O_01__ExpData) -----------------------------------------------
def loadPdDfr(dITp, dOIn, doT = False, formD = GC.NM_EXT_CSV):
    pdDfr, addIDfr, lRIDup, lCIDup = pd.DataFrame(), pd.DataFrame(), [], []
    if formD == GC.NM_EXT_CSV:  # currently the only implemented data format
        pdDfr = GF.readCSV(dITp['pRelDF'], dITp['cSep'])
    if pdDfr.size > 0:
        # drop duplicate rows
        if dITp['dropDupl']:
            if dITp['dropType'] == 'nameIdent':
                pdDfr = pdDfr.drop_duplicates(subset = pdDfr.columns[0])
            elif dITp['dropType'] == 'fullIdent':
                pdDfr = pdDfr.drop_duplicates()
        # convert specified columns to string data type
        for i in dOIn['lIC2Str']:
            pdDfr = pdDfr.astype({pdDfr.columns[i]: str})
        # filter the DataFrame by dropping unwanted entries
        pdDfr = pdDfr[[pdDfr.iloc[:, 0][k] not in dITp['lSFiltOut'] for k in
                       range(pdDfr.shape[0])]].reset_index(drop = True)
        # filter the DataFrame with a query string
        if len(dOIn['cQry']) > 0:
            pdDfr = pdDfr.query(dOIn['cQry'])
        # drop NaN values from the data if required
        if dITp['dropAx_F'] is not None:
            pdDfr.dropna(axis = dITp['dropAx_F'], inplace = True)
        # potentially transpose DataFrames, and drop the non-numeric info
        cDfr = GF.checkDoT(pdDfr, pdDfr.iloc[:, dOIn['iCVS']:], dOIn['iCNm'],
                           doT = doT)
        # create additional info DataFrame (non-numeric/unused numeric info)
        addIDfr = GF.checkDoT(pdDfr, pdDfr.iloc[:, dOIn['iCNm']:dOIn['iCVS']],
                              dOIn['iCNm'], doT = doT)
        # locate remaining duplicate rows and columns
        GF.locateDup(cDfr, list(cDfr.index.duplicated()), lRIDup)
        GF.locateDup(cDfr, list(cDfr.columns.duplicated()), lCIDup)
    return cDfr, addIDfr, lRIDup, lCIDup

def getMeanSD(dITp, dOIn, pdDfr, pdDfrAI, nmMs, nmSs, lClr):
    pRFDt = saveData(dITp, dOIn, GF.getExtDfrSelAI(pdDfrAI, pdDfr,
                                                   dOIn['lCHdrAIR']).T)
    lInf, lDfrMn, lDfrSD, lHd, lIS = [], [], [], dOIn['lFtHdMn'], dITp['lISep']
    dfrN = GF.iniPdDfr(np.zeros((len(lHd), pdDfr.shape[1]), dtype = int),
                       lSNmC = list(pdDfr.columns), lSNmR = lHd)
    for cFtHd in lHd:
        lISel = GF.getLISelR(pdDfr, cFtHd)
        lDfrMn.append(pdDfr.iloc[lISel, :].mean().to_frame(cFtHd).T)
        lDfrSD.append(pdDfr.iloc[lISel, :].std().to_frame(cFtHd).T)
        dfrN.loc[cFtHd, :] = len(lISel) - pdDfr.iloc[lISel, :].isna().sum()
    lCDat = GF.getLDataC(lHd, lIS)
    for (lDfr, sMS) in [(lDfrMn, nmMs), (lDfrSD, nmSs)]:
        pdDfr = GF.concPdDfrS(lDfr, concAx = 0, verInt = True, srt = False,
                              dropAx = dITp['dropAx_Mn'])
        pdDfrAIR = GF.getExtDfrSelAI(pdDfrAI, pdDfr, dOIn['lCHdrAIR']).T
        pRF = saveData(dITp, dOIn, pdDfrAIR, sMn = sMS)
        if dITp['plotHist']:
            PF.pltHistSelC(dITp, dOIn, pdDfr.T, pRF, lClr, sMS, lCDat)
        lInf.append((pdDfr, pRF))
    dfrN = dfrN.loc[:, list(lInf[0][0].columns)]
    return lInf[0][0], lInf[1][0], pRFDt, lInf[0][1], lInf[1][1], dfrN

# --- Functions (O_31__CombData) ----------------------------------------------
def addToDITp(dITp, lODITp, lKOvwr = []):
    for cK in lKOvwr:
        if cK in lODITp[0] and cK in lODITp[1]:
            if lODITp[0][cK] == lODITp[1][cK]:
                dITp[cK] = lODITp[0][cK]

# --- Functions (O_41__ExpDataX) ----------------------------------------------
def prepareDITp(dITp, lODat, lKOvwr = [], sSpec = None):
    assert GF.dictsEqual({OD.cGT: OD.dITp for OD in lODat}, dITp['lKeysSpec'])
    for cK, cV in lODat[0].dITp.items():   # to keep existing entries
        if cK not in dITp or cK in lKOvwr:
            dITp[cK] = cV
    for cK in dITp['lKeysSpec']:
        if cK in dITp:
            del dITp[cK]
    if sSpec is not None:
        dITp[sSpec] = GF.flattenIt([OD.dITp[sSpec] for OD in lODat])

def doDivCompOp(dITp, dfrE, dfrD):
    if list(dfrE.columns) != list(dfrD.columns):
        lC = list(set(dfrE.columns).intersection(set(dfrD.columns)))
        dfrE, dfrD = dfrE.loc[:, lC], dfrD.loc[:, lC]
    d = dfrE.reset_index(drop = True)/dfrD.reset_index(drop = True)
    if dITp['compOp'] == '2_log':
        d = np.log2(d.to_numpy().astype(np.float64))
    elif dITp['compOp'] == 'e_log':
        d = np.log(d.to_numpy().astype(np.float64))
    elif dITp['compOp'] == '10_log':
        d = np.log10(d.to_numpy().astype(np.float64))
    return d

def printWarningLUneq(ll):
    if not GF.listsEqual(ll):
        if not GF.setsEqual(ll):
            print('WARNING: Sets not equal in', ll)
        else:
            print('Lists not equal in', ll)

def moveStrRC(lDfr, lAllNm = None, useIdx = True, lPosMove = [1, 2],
              posIns = 0, prtW = False, sSep = '_'):
    if useIdx is not None and lPosMove is not None and posIns is not None:
        llNmA, lS = [], []
        for cDfr in lDfr:
            lNmA, cS = GF.redLNames(GF.getIdxOrCol(cDfr, useIdx),
                                    lAllNm = lAllNm, lPosMove = lPosMove)
            llNmA.append(lNmA)
            lS.append(cS)
            GF.setIdxOrCol(cDfr, lNmA, useIdx)
            lNmB = GF.expLNames(GF.getIdxOrCol(cDfr, not useIdx),
                                posIns = posIns, sIns = cS, sSep = sSep)
            GF.setIdxOrCol(cDfr, lNmB, not useIdx)
        if prtW:
            printWarningLUneq(llNmA)
    return lS

def createDictAttr(XAttrDfr, sSep = '_'):
    dAttr = {}
    for sC in XAttrDfr.columns:
        lSpl = sC.split(sSep)
        assert len(lSpl) >= 2
        cKS, cKM = sSep.join(lSpl[:-1]), lSpl[-1]
        GF.addToDictD(dAttr, cKM, cKS, list(XAttrDfr.loc[:, sC]))
    lLen = [len(list(dAttr[cKM][cKS])) for cKM in dAttr for cKS in dAttr[cKM]]
    assert GF.elEqual(lLen) and len(lLen) > 0
    return dAttr, lLen[0]

def createAttrDfr(dIG, dAttr, nR, sSep = '_'):
    sErrE, sErr2, sErr3 = 'ERROR:', 'List of values for', 'is empty!'
    sErr0, sErr1 = 'Elements in column', 'of attribute DataFrame not equal:'
    lSErrStd = [sErrE, sErr0, sErr1]
    attrDfr = GF.iniPdDfr(data = np.array([['']*len(dAttr)]*nR),
                          lSNmC = list(dAttr), shape = (nR, len(dAttr)))
    for cKM in dAttr:
        llV = [[l[k] for l in dAttr[cKM].values()] for k in range(nR)]
        for k, lV in enumerate(llV):
            if not GF.elEqual(lV):          # elements in lV not all equal
                serV = pd.Series(lV)
                if sum(serV.isna()) > 0:    # at least one NaN
                    attrDfr.at[k, cKM] = None
                else:                       # no NaN
                    if cKM in [dIG['idOAtDfr']]:    # object ID column
                        attrDfr.at[k, cKM] = None
                    elif cKM in [dIG['idGTAtDfr'], dIG['nmGTAtDfr']]:   # GT c.
                        if (GF.listsEqual([lV, dIG['lIDGT']]) or
                            GF.listsEqual([lV, dIG['lNmGT']])):
                            attrDfr.at[k, cKM] = dIG['nmAllGT']
                        else:
                            if (GF.allInList(dIG['lIDGT'], lV) or
                                GF.allInList(dIG['lNmGT'], lV)):
                                attrDfr.at[k, cKM] = sSep.join(lV)
                            else:
                                GF.printErrAndStop(lSErrStd, cKM, lV)
                    elif cKM in [dIG['idFtAtDfr'], dIG['nmFtAtDfr']]:   # Ft c.
                        if (GF.listsEqual([lV, dIG['lIDFt']]) or
                            GF.listsEqual([lV, dIG['lNmFt']])):
                            attrDfr.at[k, cKM] = dIG['nmAllFt']
                        else:
                            if (GF.allInList(dIG['lIDFt'], lV) or
                                GF.allInList(dIG['lNmFt'], lV)):
                                attrDfr.at[k, cKM] = sSep.join(lV)
                            else:
                                GF.printErrAndStop(lSErrStd, cKM, lV)
                    else:                   # weight column
                        attrDfr.at[k, cKM] = 0.    # weight adjusted later
            else:                           # elements in lV are all equal
                if len(lV) > 0:
                    attrDfr.at[k, cKM] = lV[0]
                else:
                    GF.printErrAndStop([sErrE, sErr2, sErr3], cKM, '')
    lIDObj = list(attrDfr.loc[:, dIG['idOAtDfr']])
    attrDfr.loc[:, dIG['nmWtAtDfr']] = GF.adjustWts(lIDObj) # adjust weights
    return attrDfr

def getAttrDfrTpXF(dIG, lODfr, lAttrDfr, lGTFt, sSep = '_'):
    assert GF.listsEqual([cAttrDfr.columns for cAttrDfr in lAttrDfr])
    lAttrDfrX, lCAttrDfrX = [], []
    for k, cAttrDfr in enumerate(lAttrDfr):
        cAttrDfr.index = lODfr[k].columns
        lSC = [sSep.join([lGTFt[k], sC]) for sC in cAttrDfr.columns]
        cAttrDfr.columns = lSC
        lAttrDfrX.append(cAttrDfr)
        GF.appendElFromL(lCAttrDfrX, lSC)
    lFullDfr = [GF.concPdDfrS([lODfr[k].T, lAttrDfrX[k]], concAx = 1)
                for k in range(len(lODfr))]
    XAttrDfr = GF.concPdDfrS(lFullDfr, concAx = 1).loc[:, lCAttrDfrX]
    dAttr, nRow = createDictAttr(XAttrDfr, sSep = sSep)
    return createAttrDfr(dIG, dAttr, nRow, sSep = sSep)

def getDfrXFt(dIG, pdDfr, cXFt, useIdx = True, lPosMove = [1, 2],
              posIns = 0, prtW = True, sSep = '_'):
    lDfr, I = [], pdDfr.index
    if cXFt == dIG['nmXFtC']:
        lNm, lP = dIG['lNmFt'], [0, 1]
    else:
        lNm, lP = dIG['lNmGT'], [1, 2]
    for cNm in lNm:
        lBI = [GF.checkIfStrOnPos(sI, cNm, lPos = lP, sSep = sSep) for sI in I]
        lDfr.append(pdDfr[lBI])
    if useIdx is not None and lPosMove is not None and posIns is not None:
        moveStrRC(lDfr, lAllNm = lNm, useIdx = useIdx, lPosMove = lPosMove,
                  posIns = posIns, prtW = prtW, sSep = sSep)
    return GF.concPdDfrS(lDfr, concAx = 1)

# --- Functions (O_51__ClusterResData) ----------------------------------------
def complDITp(dITp, cOITp, lKOvwr = []):
    for cK, cV in cOITp.items():   # to keep existing entries
        if cK not in dITp or cK in lKOvwr:
            dITp[cK] = cV

# --- Functions (O_81__BinaryOps) ---------------------------------------------
def compCorr(dtC1, dtC2):
    idxInters = dtC1.dropna().index.intersection(dtC2.dropna().index)
    rhoCr, pValPCr = stats.pearsonr(dtC1.loc[idxInters], dtC2.loc[idxInters])
    rhoSp, pValSp = stats.spearmanr(dtC1.loc[idxInters], dtC2.loc[idxInters])
    rhoKn, pValKn = stats.kendalltau(dtC1.loc[idxInters], dtC2.loc[idxInters])
    return [rhoCr, rhoSp, rhoKn] + [pValPCr, pValSp, pValKn]

def getDevR(dITp, dINP, lCmpR, lXSDR, arrD1, arrD2, lnBnd, k = 0):
    cMin, cMax = min(arrD1[k], arrD2[k]), max(arrD1[k], arrD2[k])
    # fill lCmpR
    x = dITp['dWtDv']['Min']*(min(abs(arrD1[k]), abs(arrD2[k])))
    x += dITp['dWtDv']['Max']*(max(abs(arrD1[k]), abs(arrD2[k])))
    lCmpR[k] = round(np.sign(arrD1[k])*np.sign(arrD2[k])*x,
                     dITp['nDigRndDvSc'])
    # fill lXSDR and calculate the indices dINP['iNeg'] and dINP['iPos']
    if cMin > 0:    # both arrD1[k] and arrD2[k] positive
        if cMin <= dITp['lPosCIBnd'][0]:
            lXSDR[k] = round(cMin, dITp['nDigRndCI'])
        else:       # both above lowest threshold
            for j, cThr in enumerate(reversed(dITp['lPosCIBnd'])):
                if cMin > cThr:
                    lXSDR[k] = cThr
                    dINP['iPos'] += dITp['lCIWts'][lnBnd - j - 1]
                    break
    elif cMax < 0:  # both arrD1[k] and arrD2[k] negative
        if cMax >= dITp['lNegCIBnd'][0]:
            lXSDR[k] = round(-cMax, dITp['nDigRndCI'])
        else:       # both below highest neg. threshold
            for j, cThr in enumerate(reversed(dITp['lNegCIBnd'])):
                if cMax < cThr:
                    lXSDR[k] = -cThr
                    dINP['iPos'] += dITp['lCIWts'][lnBnd - j - 1]
                    break
    else:           # one of arrD1[k] and arrD2[k] is <=0, one is >= 0
        if cMin >= dITp['lNegCIBnd'][0] or cMax <= dITp['lPosCIBnd'][0]:
            lXSDR[k] = round(-min(abs(cMin), abs(cMax)),
                             dITp['nDigRndCI'])
        else:       # both abs. values above lowest threshold
            for j, cThr in enumerate(reversed(dITp['lPosCIBnd'])):
                if cMax > cThr and cMin < -cThr:
                    lXSDR[k] = -cThr
                    dINP['iNeg'] -= dITp['lCIWts'][lnBnd - j - 1]
                    break

def compDv(dITp, arrD1, arrD2):
    assert len(arrD1) == len(arrD2) and len(dITp['lPosCIBnd']) > 0
    lenArr, lenBnd = len(arrD1), len(dITp['lPosCIBnd'])
    # calculate the basic score list - comparing single features: lCmpR
    # calculate which features have > x SD deviation for both data: lXSDR
    dINP, lCmpR, lXSDR = {'iNeg': 0., 'iPos': 0.}, ['']*lenArr, ['']*lenArr
    for k in range(lenArr):
        if arrD1[k] != '' and arrD2[k] != '':
            getDevR(dITp, dINP, lCmpR, lXSDR, arrD1, arrD2, lenBnd, k)
    # calculate the sums of the neg. and pos. values of the basic score list
    lCmpRS = [sum(x for x in lCmpR if x != '' and x < 0),
              sum(x for x in lCmpR if x != '' and x > 0)]
    lCmpRS.append(lCmpRS[0] + lCmpRS[1])    # sum of neg. and pos. DvSc
    # initialise lXSDRS with the indices dINP['iNeg'] and dINP['iPos']
    # and the sum of these (neg. and pos. CI)
    lXSDRS = [dINP['iNeg'], dINP['iPos'], dINP['iNeg'] + dINP['iPos']]
    # calculate the occurrences of neg. and pos. values of lXSDR
    for cThr in list(reversed(dITp['lNegCIBnd'])) + dITp['lPosCIBnd']:
        lXSDRS += [sum(1 for x in lXSDR if x == cThr)]
    return (lCmpR + lXSDR + GF.roundList(lCmpRS, dITp['nDigRndDvSc']) +
            GF.roundList(lXSDRS, dITp['nDigRndCI']))

def fillDictRCrDv(dIG, dITp, dR, lDfr, lDfrDvT, lDfrAIC, sMn = ''):
    assert len(lDfrDvT) == len(lDfrAIC)
    cSt, nSt, lDvD1, lDvD2 = 0, lDfr[0].shape[1]*lDfr[1].shape[1], [], []
    sPrnt = '\nCalculating correlation and deviation patterns for'
    if len(sMn) > 0:
        sPrnt += ' MEAN data...'
    else:
        sPrnt += ' ALL data...'
    print(sPrnt)
    for nm1, datC1 in lDfr[0].items():
        lAIC1 = []
        if len(lDfrDvT) > 0:
            lDvD1 = GF.getRowDfr(lDfrDvT[0], nm1)
            lAIC1 = GF.getRowDfr(lDfrAIC[0], nm1)
        for nm2, datC2 in lDfr[1].items():
            lAIC2 = []
            if len(lDfrDvT) > 1:
                lDvD2 = GF.getRowDfr(lDfrDvT[1], nm2)
                lAIC2 = GF.getRowDfr(lDfrAIC[1], nm2)
            dR[(nm1, nm2)] = lAIC1 + lAIC2 + compCorr(datC1, datC2)
            if sMn == '':
                if len(lDvD1) > 0 and len(lDvD2) > 0:
                    arr1, arr2 = np.array(lDvD1), np.array(lDvD2)
                    dR[(nm1, nm2)] += GF.roundList(lDvD1, dITp['nDigRndDv'])
                    dR[(nm1, nm2)] += GF.roundList(lDvD2, dITp['nDigRndDv'])
                    dR[(nm1, nm2)] += compDv(dITp, arr1, arr2)
            cSt += 1
            if cSt%dIG['nStPr'] == 0 or cSt == nSt:
                print(str(round(cSt/nSt*100., 1)) + '% calculated.')

def iniDataTACD(dITp, tID, lNmI, sMn = '', sSep = '_'):
    tID = tuple([cID.split(sSep)[0] for cID in tID])
    assert (tID in dITp['dNTopCrDv'] or (tID[1], tID[0]) in dITp['dNTopCrDv'])
    if tID in dITp['dNTopCrDv']:
        lNTop = dITp['dNTopCrDv'][tID]
    elif (tID[1], tID[0]) in dITp['dNTopCrDv']:
        lNTop = dITp['dNTopCrDv'][(tID[1], tID[0])]
    # initialise the returned DataFrame and the list of columns
    lSAv = copy.deepcopy(dITp['lSAvAll'])
    lSTop = copy.deepcopy(dITp['lSTopAll'])
    if sMn == '':
        lSAv += (dITp['lSAvDvSc'] + dITp['lSAvCI'])
        lSTop += (dITp['lSTopDvSc'] + dITp['lSTopCI'])
    lNmCol = [s for s in lSAv]
    lNmCol += [s + sSep + str(n) for s in lSTop for n in lNTop]
    pdDfrR = pd.DataFrame(index = lNmI, columns = lNmCol)
    return pdDfrR, lNTop

def calcAvSc(dITp, pdDfr, lCtHd, lNmI, i = 0, sC = ''):
    if len(sC) == 0:
        sC = dITp['tCNmRes'][i + 1]
    # calculate average scores
    arrAv = np.zeros((len(lNmI), len(lCtHd)))
    for j, nmI in enumerate(lNmI):
        dfrFlt = pdDfr[pdDfr[sC] == nmI]
        lAvNP = []
        for s in dITp['lSCorrVals']:
            lAvNP.append(dfrFlt[dfrFlt[s] < 0][s].mean())
            lAvNP.append(dfrFlt[dfrFlt[s] > 0][s].mean())
        for sHd in lCtHd[len(dITp['lSCorrVals2x']):]:
            lAvNP += [dfrFlt[sHd].mean()]
        arrAv[j, :] = lAvNP
    print('Calculated average negative/positive score per item.')
    return arrAv

def fillTop(dITp, arrT, dfrRd, dOcc, lNT, lNmI, sC, iCt, j, doWt = False):
    for k, nmI in enumerate(lNmI):
        dfrFlt = dfrRd[dfrRd[sC] == nmI]
        if iCt in range(len(dITp['lSCorrVals'])):       # a "double" header
            if dfrFlt.shape[0] > 0:
                sumNegSc = dfrFlt[dfrFlt['Score'] < 0]['Score'].sum()
                sumPosSc = dfrFlt[dfrFlt['Score'] > 0]['Score'].sum()
                if doWt:
                    sumNegSc /= dOcc[nmI]
                    sumPosSc /= dOcc[nmI]
                arrT[k, 2*iCt*len(lNT) + j] = sumNegSc
                arrT[k, (2*iCt + 1)*len(lNT) + j] = sumPosSc
        else:
            sumSc = dfrFlt['Score'].sum()
            if doWt:
                sumSc /= dOcc[nmI]
            arrT[k, (len(dITp['lSCorrVals']) + iCt)*len(lNT) + j] = sumSc

def calcTopSc(dITp, pdDfr, lNTop, lCtHd, lNmI, i = 0, sC = '', doWt = False):
    if len(sC) == 0:
        sC = dITp['tCNmRes'][i + 1]
    dOcc = pdDfr[sC].value_counts().to_dict()
    # calculate top scores
    arrTop = np.zeros((len(lNmI), len(lNTop)*len(lCtHd)))
    lSHd4Top = dITp['lSCorrVals'] + lCtHd[len(dITp['lSCorrVals2x']):]
    for iCt, sHd in enumerate(lSHd4Top):
        bAsc = GF.getBoolLDict([dITp['dSDvSc'], dITp['dSCI']], sHd)
        dfrSrt = pdDfr.sort_values(sHd, ascending = bAsc)
        for j, nTop in enumerate(lNTop):
            nTop = min(nTop, dfrSrt.shape[0])
            dfrRed = dfrSrt.iloc[:nTop, :]
            lSc = list(range(nTop, 0, -1))
            if iCt in range(len(dITp['lSCorrVals'])):   # a "double" header
                dfrRed = dfrRed.append(dfrSrt.iloc[-nTop:, :],
                                       ignore_index = True)
                lSc += list(range(-1, -(nTop + 1), -1))
            dfrRed = dfrRed.assign(Score = lSc)
            fillTop(dITp, arrTop, dfrRed, dOcc, lNTop, lNmI, sC, iCt, j, doWt)
            print('Calculated', sHd, 'top scores for TOP', str(nTop) + '.')
    return arrTop

def calcTACD(dITp, dOIn, pdDfr, lNmI, i, sC = '', doWt = False, sMn = '',
             sSep = '_'):
    pdDfrR, lNTop = iniDataTACD(dITp, dOIn['tID'], lNmI, sMn, sSep)
    lH = copy.deepcopy(dITp['lSCorrVals2x'])
    if sMn == '':
        lH += (list(dITp['dSDvSc']) + list(dITp['dSCI']))
    pdDfrR.iloc[:, :len(lH)] = np.around(calcAvSc(dITp, pdDfr, lH, lNmI, i = i,
                                                  sC = sC), dITp['nDigRndBO'])
    pdDfrR.iloc[:, len(lH):] = calcTopSc(dITp, pdDfr, lNTop, lH, lNmI, i = i,
                                         sC = sC, doWt = doWt)
    print('Calculated all scores.')
    dSrt = dITp['dSrtR']
    if set(dITp['dSrtF'][GC.S_IDX]).issubset(set(pdDfrR)):
        dSrt = dITp['dSrtF']
    lSrt = list(dSrt[GC.S_IDX])
    lAsc = [dSrt[GC.S_IDX][cK]['Asc'] for cK in dSrt[GC.S_IDX]]
    return pdDfrR.sort_values(by = lSrt, ascending = lAsc)

def writeCrDvDatToF(dITp, cD, pF, lSHdr = []):
    fCD = open(pF, 'w')
    if len(lSHdr) > 0:
        for cS in lSHdr[:-1]:
            fCD.write(cS + dITp['cSep'])
        fCD.write(lSHdr[-1] + '\n')
    for i, (cK, cV) in enumerate(cD.items()):
        fCD.write(str(i + 1))
        for s in cK:
            fCD.write(dITp['cSep'] + str(s))
        for s in cV:
            fCD.write(dITp['cSep'])
            if s != '':
                if type(s) in [int, float, np.float64, np.int64, np.uint64]:
                    fCD.write(str(round(s, dITp['nDigRndBO'])))
                else:
                    fCD.write(str(s))
        fCD.write('\n')
    fCD.close()

def saveTACDData(dITp, pdDfr, lDfrAIC, k, lHdrAIC, pF, concAI = True):
    assert len(lDfrAIC) == 0 or len(lDfrAIC) >= 2
    if not os.path.isfile(pF):
        if len(lDfrAIC) >= 2 and concAI:
            pdDfr = GF.getExtDfrSelAI(lDfrAIC[k], pdDfr, lHdrAIC, isT = True)
        pdDfr.to_csv(pF, sep = dITp['cSep'])

# --- Functions (O_82__Clustering) --------------------------------------------
def get2xPF(dITp, dOIn, sMn = '', nmFExt = GC.NM_EXT_CSV):
    pF = getPF(dITp, getNmF(dOIn, sMn = sMn), nmFExt)
    return pF, pF

def get2xPTrF(dITp, dOIn, sMn = '', nmFExt = GC.NM_EXT_CSV):
    pF = getPF(dITp, getNmF(dOIn, sF = 'sFNmTr', sMn = sMn), nmFExt)
    return pF, pF

# --- Functions (O_83__OverRep) -----------------------------------------------
def getParProf(dITp, iTpO, isClR, lElCtDef):
    lElCt = None
    if iTpO == 0:
        sID = 'PD'
    elif iTpO == 1:
        sID = 'BO'
        if isClR:            # a "BinOps" from cluster result data object
            sID = 'CR'
    cDSrt, lElCt = dITp['dSrt' + sID], dITp['lElC' + sID]
    nMn, nMx = dITp['nMin' + sID], dITp['nMax' + sID]
    lSrt = list(cDSrt[GC.S_IDX])
    lAsc = [cDSrt[GC.S_IDX][cK]['Asc'] for cK in cDSrt[GC.S_IDX]]
    if lElCt is None:    # no category to evaluate defined
        lElCt = lElCtDef
    return cDSrt, lSrt, lAsc, lElCt, nMn, nMx

def getCntTblF(dOccAbs, nCur, nMin, nOccAbs, nTot, cArr, cAttr):
    nOccCur = len(cArr[cArr == cAttr])
    dOccAbs[cAttr][nCur - nMin] = nOccCur
    a00, a01 = nCur - nOccCur, nOccCur
    a10, a11 = nTot - nCur - nOccAbs + nOccCur, nOccAbs - nOccCur
    if min([a00, a01, a10, a11]) < 0:
        print('ERROR: [a00, a01, a10, a11] =', [a00, a01, a10, a11])
        assert False
    return np.array([[a00, a01], [a10, a11]], dtype = np.uint64)

def exactFisher(arrDat, cAlt = 'two-sided'):
    arrDat = np.array(arrDat, dtype = np.uint64)
    assert arrDat.shape == (2, 2)
    # returns (oddsRatio, pVal)
    return stats.fisher_exact(arrDat, alternative = cAlt)

def calcPValsF(dPValOv, dPValUn, arrCTblF, nCur, nMin, nCAttr, cAttr,
               sCorrect = None):
    pValOv = exactFisher(arrCTblF, cAlt = 'less')[1]
    pValUn = exactFisher(arrCTblF, cAlt = 'greater')[1]
    if sCorrect == GC.S_M_CORR_BON:    # Bonferroni correction
        pValOv *= nCAttr
        pValUn *= nCAttr
    dPValOv[cAttr][nCur - nMin] = pValOv
    dPValUn[cAttr][nCur - nMin] = pValUn

def calcPValProfiles(dITp, dOccAbs, dPValOv, dPValUn, dfrRd, lSerVC, llAttr,
                     lNAttr, nMin, nMax, nTot, lElCt):
    if lElCt is not None:
        assert len(llAttr) == len(lNAttr) and len(llAttr) == len(lElCt)
    nCalc, lToInt = sum([len(l) for l in llAttr]), [True, False, False]
    for i in range(len(llAttr)):
        for k, cAt in enumerate(llAttr[i]):
            GF.iniKeyV1D([dOccAbs, dPValOv, dPValUn], cAt,
                         lenArr = nMax + 1 - nMin, convInt = lToInt)
            for nCur in range(nMin, nMax + 1):
                if lElCt is not None:
                    arrS = dfrRd.iloc[:nCur, :].loc[:, lElCt[i]].values
                else:
                    arrS = dfrRd.iloc[:nCur, :].index.to_numpy()
                arrCntTblF = getCntTblF(dOccAbs, nCur, nMin, lSerVC[i].at[cAt],
                                        nTot, arrS, cAt)
                calcPValsF(dPValOv, dPValUn, arrCntTblF, nCur, nMin,
                           lNAttr[i], cAt, sCorrect = dITp['sMCorrectL'])
            if (k + 1)%10 == 0 or k + 1 == nCalc:
                print('Done:', k + 1, 'of', nCalc)

def saveDictAsPdDfr(dITp, lPF, dSrt, cD, nMin, nMax, k = 0):
    cDfr = pd.DataFrame(cD, index = range(nMin, nMax + 1))
    cDfr.columns = cDfr.columns.astype(str)
    if dITp['dTpY'][dITp['lTpY'][k]][2]:
        cDfr = -np.log10(cDfr)
    srtMd, srtDg = dSrt[GC.S_COL]['Srt']
    cDfr = GF.reIndexSpec(cDfr, srtMd, srtDg)
    cDfr = cDfr.round(dITp['dTpY'][dITp['lTpY'][k]][1])
    cDfr.to_csv(lPF[k], sep = dITp['cSep'])
    return cDfr

def selDataThr(dITp, cDfr, dSrt, k = 0):
    dfrSelD = pd.DataFrame(index = cDfr.index, columns = cDfr.columns)
    cThr = dITp['thrProf']
    if dITp['dTpY'][dITp['lTpY'][k]][2]:
        cThr = -np.log10(cThr)
    for sC in cDfr.columns:
        cSer = cDfr.loc[:, sC]
        if dITp['sComp'] == '>':
            dfrSelD.loc[:, sC] = cSer[cSer > cThr]
        elif dITp['sComp'] == '>=':
            dfrSelD.loc[:, sC] = cSer[cSer >= cThr]
        elif dITp['sComp'] == '==':
            dfrSelD.loc[:, sC] = cSer[cSer == cThr]
        elif dITp['sComp'] == '<=':
            dfrSelD.loc[:, sC] = cSer[cSer <= cThr]
        elif dITp['sComp'] == '<':
            dfrSelD.loc[:, sC] = cSer[cSer < cThr]
        else:
            dfrSelD.loc[:, sC] = cSer
    dfrSelD = dfrSelD.dropna(axis = 1, how = 'all')
    srtMd, srtDg = dSrt[GC.S_COL]['Srt']
    return GF.reIndexSpec(dfrSelD, srtMd, srtDg), cThr

###############################################################################

# -*- coding: utf-8 -*-
###############################################################################
# --- F_00__GenFunctions.py ---------------------------------------------------
###############################################################################
import os, time, itertools, copy

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

import Core.C_00__GenConstants as GC

# --- Functions ---------------------------------------------------------------
def startSimu():
    startTime = time.time()
    print('+'*50 + ' START', time.ctime(startTime), '+'*30)
    print('PhosphoData Analysis')
    plt.close('all')
    mpl.rcParams['figure.max_open_warning'] = GC.FIGURE_MAX_OPEN_WARNING
    return startTime

def seedRNG(cMode):
    if cMode == GC.M_STOCH:
        np.random.seed()
        print('Seeded RNG.')

def createDir(pF):
    if not os.path.isdir(pF):
        os.mkdir(pF)

def joinToPath(pF = '', nmF = 'Dummy.txt'):
    if len(pF) > 0:
        createDir(pF)
        return os.path.join(pF, nmF)
    else:
        return nmF

def joinDirToPath(pF = '', nmD = 'Directory'):
    if len(pF) > 0:
        pF = os.path.join(pF, nmD)
        createDir(pF)
        return pF
    else:
        return nmD

def combToStr(lS, sSt = '', sEnd = '', sSep = '_'):
    cS = sSt
    if len(lS) > 0:
        cS += lS[0]
        for s in lS[1:]:
            if len(s) > 0 and len(cS) > len(sSt):
                cS += sSep
            cS += s
    cS += sEnd
    return cS

def insStrAtFNmEnd(sNmF, sAdd, XSep = '.'):
    lSpl = sNmF.split(XSep)
    assert len(lSpl) >= 1
    return (XSep.join([lSpl[k] for k in range(len(lSpl) - 1)]) + sAdd + XSep +
            lSpl[-1])

def modStrEnd(cS, sEndCur, sEndNew = ''):
    if len(sEndCur) <= len(cS):
        if cS.endswith(sEndCur):
            return cS[:-len(sEndCur)] + sEndNew
        return cS
    return cS

def modStrEndLStr(lS, sEndCur, sEndNew = ''):
    lSRet = []
    for cS in lS:
        lSRet.append(modStrEnd(cS, sEndCur, sEndNew))
    return lSRet

def modFName(pF, sPre = '', sPost = ''):
    lPF, nmF = pF.split('.'), os.path.basename(pF)
    nmF = sPre + nmF[:-(len(lPF[-1]) + 1)] + sPost + nmF[-(len(lPF[-1]) + 1):]
    return os.path.join(os.path.dirname(pF), nmF)

def changeFExt(pF, newX = GC.NM_EXT_PDF):
    lPF = pF.split('.')
    return pF[:-len(lPF[-1])] + newX

def changePath(pFC, pN):
    return os.path.join(joinToPath(pN, os.path.basename(pFC)))

def adaptPF4Plot(pFC, pN, sPre = '', sPost = '', newX = GC.NM_EXT_PDF):
    return changePath(changeFExt(modFName(pFC, sPre, sPost), newX), pN)

def roundList(cL, nDigRnd = None):
    if nDigRnd is None:
        return cL
    else:
        return [round(x, nDigRnd) for x in cL]

def addToListUnique(cL, lElToAdd):
    for cE in lElToAdd:
        if cE not in cL:
            cL.append(cE)

def addToDictL(cD, cK, cE):
    if cK in cD:
        cD[cK].append(cE)
    else:
        cD[cK] = [cE]

def addToDictD(cD, cKMain, cKSub, lV = []):
    if cKMain in cD:
        if cKSub not in cD[cKMain]:
            cD[cKMain][cKSub] = lV
        else:
            if not listsEqual([cD[cKMain][cKSub], lV]):
                print('ERROR: The following lists are unequal:')
                print(cD[cKMain][cKSub], '\n', lV)
                assert False
    else:
        cD[cKMain] = {}
        cD[cKMain][cKSub] = lV

def calcDictOcc(cIter, cIncr = 1):
    dOcc = {}
    for cEl in cIter:
        if cEl in dOcc:
            dOcc[cEl] += cIncr
        else:
            dOcc[cEl] = 1
    return dOcc

def complDict(cD, lK = [], defEl = []):
    for cK in lK:
        if cK not in cD:
            cD[cK] = defEl

def appendElFromL(lFull, lCur):
    for cE in lCur:
        lFull.append(cE)

def flattenIt(cIterable, retArr = False):
    itFlat = list(itertools.chain.from_iterable(cIterable))
    if retArr:
        itFlat = np.array(itFlat)
    return itFlat

def getBoolLDict(lD, cK, bDef = False):
    bRet = bDef
    for cD in lD:
        if cK in cD:
            bRet = cD[cK]
            break
    return bRet

def getNm(cD, cK, nmDef = '', useDict = True, sSep = '_'):
    lSpl, cNm = cK.split(sSep), ''
    for i, cSpl in enumerate(lSpl):
        if i > 0:
            cNm += sSep
        if cSpl == nmDef or cSpl == '':
            cNm += nmDef
        else:
            if useDict:
                cNm += cD[cSpl]
            else:
                cNm += cSpl
    return cNm

def setIdxOrCol(pdDfr, newL, isIdx = True):
    if isIdx:
        pdDfr.index = list(newL)[:len(pdDfr.index)]
    else:
        pdDfr.columns = list(newL)[:len(pdDfr.columns)]

def getIdxOrCol(pdDfr, isIdx = True):
    if isIdx:
        cL = pdDfr.index
    else:
        cL = pdDfr.columns
    return list(cL)

def getLSubIdxOrCol(pdDfr, lSpec = [], isIdx = True):
    lNmSel, lAllNm = [], pdDfr.index
    if not isIdx:
        lAllNm = pdDfr.columns
    for cNmSpec in lSpec:
        for cNmFull in lAllNm:
            if cNmSpec in cNmFull:
                lNmSel.append(cNmFull)
    return lNmSel

def getLNmColDtX(pdDfr, cD, cK, useDict = False, sSep = '_'):
    s = getNm(cD, cK, useDict = useDict)
    lNm = [s + sSep + cNmC for cNmC in pdDfr.columns]
    if elEqual([cNmC.split(sSep)[0] for cNmC in pdDfr.columns]):
        if pdDfr.columns.size > 0:
            if pdDfr.columns[0].split(sSep)[0] == s:
                lNm = [cNmC for cNmC in pdDfr.columns]
    return lNm

def confineC(i1, i2, isTopI = False):
    iR = -1
    if i1 >= 0 and i2 >= 0:
        iR = max(i1, i2)
        if isTopI:
            iR = min(i1, i2)
    elif (i1 >= 0 and i2 < 0) or (i1 < 0 and i2 >= 0):
        iR = max(i1, i2)
    return iR

def getLISel(pdDfr, sSt, lPosSel = None, isRow = True, sSep = '_'):
    lISel, itItems = [], pdDfr.iterrows()
    if not isRow:
        itItems = pdDfr.items()
    for i, (cI, _) in enumerate(itItems):
        # get the selection indices corresponding to the given start string
        if lPosSel is None:
            if cI.startswith(sSt):
                lISel.append(i)
        else:
            if checkIfStrOnPos(cI, sSt, lPos = lPosSel, sSep = sSep):
                lISel.append(i)
    return lISel

def getLISelC(pdDfr, sSt, lPosSel = None, sSep = '_'):
    return getLISel(pdDfr, sSt, lPosSel = lPosSel, isRow = False, sSep = sSep)

def getLISelR(pdDfr, sSt, lPosSel = None, sSep = '_'):
    return getLISel(pdDfr, sSt, lPosSel = lPosSel, isRow = True, sSep = sSep)

def getDDfrFtSel(cDfr, dID, lPosSel = [0, 1]):
    cD = {}
    for cID, cNm in dID.items():
        cD[cID] = cDfr.iloc[:, getLISelC(cDfr, cNm, lPosSel = lPosSel)]
    return cD

def locateDup(pdDfr, lBoolDup, lIDup):
    for k in range(len(lBoolDup)):
        if lBoolDup[k] == True:
            lIDup.append(k)     # 0-based index

def printDfrOrDfrT(pdDfr, doT = False):
    if doT:
        print('*'*8, 'Transposed', '*'*8)
        print(pdDfr.T)
    else:
        print(pdDfr)

def checkDoT(pdDfrF, pdDfrR, iCNm, doT = False):
    lNm = list(pdDfrF.iloc[:, iCNm])
    if doT:
        pdDfrF = pdDfrR.T
        pdDfrF.columns = lNm
    else:
        pdDfrF = pdDfrR
        pdDfrF.index = lNm
    return pdDfrF

def iniNpArr(shape = (0)):
    return np.zeros(shape = shape)

def iniKeyV1D(lD, cKey, lenArr, convInt = None):
    for k, cD in enumerate(lD):
        if convInt is None:
            cD[cKey] = np.zeros((lenArr))
        else:
            assert len(convInt) == len(lD)
            if convInt[k]:
                cD[cKey] = np.zeros((lenArr), dtype = np.uint64)
            else:
                cD[cKey] = np.zeros((lenArr))

def iniPdDfr(data = None, lSNmC = [], lSNmR = [], shape = (0, 0)):
    assert len(shape) == 2
    nR, nC = shape
    if len(lSNmC) == 0:
        if len(lSNmR) == 0:
            if data is None:
                return pd.DataFrame(np.zeros(shape))
            else:
                return pd.DataFrame(data)
        else:
            if data is None:
                return pd.DataFrame(np.zeros((len(lSNmR), nC)), index = lSNmR)
            else:
                return pd.DataFrame(data, index = lSNmR)
    else:
        if len(lSNmR) == 0:
            if data is None:
                return pd.DataFrame(np.zeros((nR, len(lSNmC))),
                                    columns = lSNmC)
            else:
                return pd.DataFrame(data, columns = lSNmC)
        else:   # ignore nR
            if data is None:
                return pd.DataFrame(np.zeros((len(lSNmR), len(lSNmC))),
                                    index = lSNmR, columns = lSNmC)
            else:
                return pd.DataFrame(data, index = lSNmR, columns = lSNmC)

def iniAttrDfr(lAttr, lNmAttr, cLen = 1):
    assert len(lAttr) == len(lNmAttr)
    if len(lAttr) == 0:
        return iniPdDfr()
    else:
        cData = {lNmAttr[k]: [lAttr[k]]*cLen for k in range(len(lAttr))}
        return iniPdDfr(data = cData, lSNmC = lNmAttr)

def getLDataC(lHd, lIS):
    lDat = [cHd for cHd in lHd]
    if len(lIS) > 0:
        lDat = [lHd[:min(lIS[0], len(lHd))]]
        lDat += [lHd[lIS[k]:lIS[k + 1]] for k in range(len(lIS) - 1)]
        lDat += [lHd[lIS[-1]:]]
    return lDat

def getPooledSD(tIn):
    # unpack the input
    n1, n2, SD1, SD2 = tIn
    # check and potentially convert types
    assert type(n1) == type(n2) and type(SD1) == type(SD2)
    if type(n1) in [list, np.ndarray]:
        assert type(SD1) in [list, np.ndarray]
        if type(n1) == list:
            n1 = np.array(n1)
            n2 = np.array(n2)
        if type(SD1) == list:
            SD1 = np.array(SD1)
            SD2 = np.array(SD2)
        for cArr in [n2, SD1, SD2]:
            assert cArr.shape == n1.shape
    elif type(n1) in [int, float]:
        assert type(SD1) in [int, float]
    # return the pooled standard deviation
    return np.sqrt(((n1 - 1)*SD1*SD1 + (n2 - 1)*SD2*SD2)/(n1 + n2 - 2))

def prepDfrForSDP(dfr1, dfr2, j, k):
    return (dfr1.iloc[j, :].to_numpy(), dfr1.iloc[k, :].to_numpy(),
            dfr2.iloc[j, :].to_numpy(), dfr2.iloc[k, :].to_numpy())

def getDfrDv(pdDfrMn, pdDfrSD, dfrN, devTp = None, sSt = '', sMid = '|'):
    shArr = ((pdDfrMn.shape[0] - 1)*pdDfrMn.shape[0], pdDfrMn.shape[1])
    lI, lS, arrD, k = list(range(len(pdDfrMn.index))), [], np.zeros(shArr), 0
    cDenom = 1
    for iB in lI:
        lIC = copy.deepcopy(lI)
        lIC.remove(iB)
        for i in lIC:
            lS.append(sSt + pdDfrMn.index[i] + sMid + pdDfrMn.index[iB])
            if devTp == 'SD':
                cDenom = pdDfrSD.iloc[iB, :]
            elif devTp == 'SDP':
                cDenom = getPooledSD(prepDfrForSDP(dfrN, pdDfrSD, i, iB))
            v = (pdDfrMn.iloc[i, :] - pdDfrMn.iloc[iB, :])/cDenom
            arrD[k, :] = v
            k += 1
    return iniPdDfr(arrD, lSNmC = pdDfrMn.columns, lSNmR = lS)

def evalStdOp(stdOp = 'N'):
    shiftMn, shiftSD = False, False
    if len(stdOp) > 1:
        if len(stdOp) > 2:
            if stdOp[1] == 'M' or stdOp[2] == 'M':
                shiftMn = True
            if stdOp[1] == 'S' or stdOp[2] == 'S':
                shiftSD = True
        else:
            if stdOp[1] == 'M':
                shiftMn = True
            elif stdOp[1] == 'S':
                shiftSD = True
    return stdOp[0], shiftMn, shiftSD

def standardPdDfr(pdDfr, stdOp = 'N'):
    # evaluate the stdOp string
    cRef, shiftMn, shiftSD = evalStdOp(stdOp = stdOp)
    # apply the operations coded in the stdOp string
    if cRef == 'T':
        cMn, cSD, opAx = pdDfr.stack().mean(), pdDfr.stack().std(), 1
    elif cRef == 'C':
        cMn, cSD, opAx = pdDfr.mean(axis = 0), pdDfr.std(axis = 0), 1
    elif cRef == 'R':
        cMn, cSD, opAx = pdDfr.mean(axis = 1), pdDfr.std(axis = 1), 0
    else:
        return pdDfr
    if shiftMn:
        if shiftSD:
            return (pdDfr.sub(cMn, axis = opAx)).div(cSD, axis = opAx)
        else:
            return pdDfr.sub(cMn, axis = opAx)
    else:
        if shiftSD:
            return pdDfr.div(cSD, axis = opAx)
        else:
            return pdDfr

def toPdDfr(dCoord, lSNmR = []):
    assert len(dCoord) > 0
    nR = len(list(dCoord.values())[0])
    for cK, cV in dCoord.items():
        assert len(cV) == nR
    if len(lSNmR) != nR:
        lSNmR = range(nR)
    return pd.DataFrame(dCoord, index = lSNmR)

def getLHStrISpl(cD, cK, lHStrDef = [], iDef = 0):
    lHStr, iSpl = lHStrDef, iDef
    if cK in cD:
        lHStr, iSpl = cD[cK], len(cD[cK])
    return lHStr, iSpl

def readCSV(pF, sepD = ',', iCol = None, dDtype = None, lHStr = [],
            iSp = None, splC = True):
    if dDtype is None and lHStr is not None:
        dDtype = {s: str for s in lHStr}
    pdDfr = pd.read_csv(pF, sep = sepD, index_col = iCol, dtype = dDtype)
    if iSp is not None:
        if splC:
            addIDfr = pdDfr.iloc[:, :iSp]
            return pdDfr.iloc[:, iSp:], addIDfr
        else:
            addIDfrT = pdDfr.iloc[:iSp, :].T
            return pdDfr.iloc[iSp:, :], addIDfrT.T
    else:
        return pdDfr

def shrinkPdDfr(pdDfr, tIMnMx = (-1, -1)):
    assert len(tIMnMx) >= 2
    if tIMnMx[0] < 0:
        if tIMnMx[1] < 0:
            return pdDfr    # do not shrink DataFrame
        else:
            return pdDfr.iloc[:, :tIMnMx[1]]
    else:
        if tIMnMx[1] < 0:
            return pdDfr.iloc[:, tIMnMx[0]:]
        else:
            return pdDfr.iloc[:, tIMnMx[0]:tIMnMx[1]]

def concPdDfrS(lPdDfr, concAx = 0, verInt = True, srt = False, ignIdx = False,
               dropAx = None):
    d = pd.concat(lPdDfr, axis = concAx, verify_integrity = verInt, sort = srt,
                  ignore_index = ignIdx)
    if dropAx in [0, 1, 'index', 'columns']:
        d.dropna(axis = dropAx, inplace = True)
    return d

def concPdDfrWt(lPdDfr, concAx = 0, dropAx = None, verInt = True, srt = False):
    pdDfr = concPdDfrS(lPdDfr, concAx, verInt, srt, dropAx = dropAx)
    arrNDS = np.array([cPdDfr.shape[concAx] for cPdDfr in lPdDfr])
    arrWts = np.array([[sum(arrNDS)/(arrNDS.size*cN)]*cN for cN in arrNDS])
    return pdDfr, flattenIt(arrWts, retArr = True)

def adjustWts(lIDO):
    dOcc = calcDictOcc(lIDO)
    return np.array([sum(dOcc.values())/(len(dOcc)*dOcc[cID]) for cID in lIDO])

def concNpArr(lArr, concAx = 0):
    if sum([np.array(cArr).ndim for cArr in lArr]) >= len(lArr):
        return np.concatenate(lArr, axis = concAx)

def getDfrSelAI(pdDfr, lHdr, isT = False):
    if isT:
        return pdDfr.loc[:, lHdr]
    else:
        return pdDfr.loc[lHdr, :]

def getExtDfrSelAI(pdDfrA, pdDfrB, lHdrA, isT = False):
    if isT:
        return concPdDfrS([pdDfrA.loc[:, lHdrA], pdDfrB], concAx = 1)
    else:
        return concPdDfrS([pdDfrA.loc[lHdrA, :], pdDfrB])

def convPdDfrToDict(pdDfr, lNmCKey, lNmCVal):
    print('Converting DataFrame to dictionary...')
    cD = {}
    for r in pdDfr.iterrows():
        tKey, tVal = [], []
        for nmCKey in lNmCKey:
            tKey.append(r[1].loc[nmCKey])
        for nmCVal in lNmCVal:
            tVal.append(r[1].loc[nmCVal])
        cD[tuple(tKey)] = tuple(tVal)
    print('Converted DataFrame to dictionary.')
    return cD

def extractFeatureHdr(pdDfr, fromR = True):
    lFtHdr = []
    if fromR:
        lFtHdr = list(pdDfr.index)
    else:
        lFtHdr = list(pdDfr.columns)
    return lFtHdr

def getFeatureHdr(pdDfr, lISpl = [0, 1], fromR = True, sSep = '_'):
    lFtHdr, lFtHdrMn = extractFeatureHdr(pdDfr), []
    lAllEl = [sSep.join([s.split(sSep)[i] for i in lISpl]) for s in lFtHdr]
    addToListUnique(lFtHdrMn, lAllEl)
    return lFtHdr, lFtHdrMn

def getVal1(cD, cDet, cK, sPre = '', sPost = '', retVDef = ''):
    retV = retVDef
    if cDet in cD:
        if cD[cDet]:
            if cK in cD:
                if sPre == '' and sPost == '':
                    retV = cD[cK]   # to avoid string concatenation errors
                else:
                    retV = sPre + str(cD[cK]) + sPost
    return retV

def selVal1(cD, cDet, cKF, cKT, sPreF = '', sPostF = '', sPreT = '',
            sPostT = ''):
    if cKF in cD:
        if sPreF == '' and sPostF == '':
            rVD = cD[cKF]           # to avoid string concatenation errors
        else:
            rVD = sPreF + str(cD[cKF]) + sPostF
        return getVal1(cD, cDet, cKT, sPre = sPreT, sPost = sPostT,
                       retVDef = rVD)
    else:
        return getVal1(cD, cDet, cKT, sPre = sPreT, sPost = sPostT,
                       retVDef = sPreF + sPostF)

def insStrLDfr(lDfr, lS, iP = 0, useC = False):
    assert len(lDfr) == len(lS)
    for k, cDfr in enumerate(lDfr):
        if useC:
            cDfr.columns = [sC[:iP] + lS[k] + sC[iP:] for sC in cDfr.columns]
        else:
            cDfr.index = [sC[:iP] + lS[k] + sC[iP:] for sC in cDfr.index]

def addString(sRetY, sPre = '', sPost = '', sRetN = ''):
    if sRetY is not None:
        if len(sRetY) > 0:
            return sPre + sRetY + sPost
    return sRetN

def joinStrings(lNm, k, sS = '', sM = '', sE = '', sSep = '_'):
    if len(sS) > 0 and len(sM) > 0 and len(sE) > 0:
        lNm[k] = sSep.join([sS, sM, sE])
    elif len(sS) > 0 and len(sM) > 0 and len(sE) <= 0:
        lNm[k] = sSep.join([sS, sM])
    elif len(sS) > 0 and len(sM) <= 0 and len(sE) > 0:
        lNm[k] = sSep.join([sS, sE])
    elif len(sS) > 0 and len(sM) <= 0 and len(sE) <= 0:
        lNm[k] = sSep.join([sS])
    elif len(sS) <= 0 and len(sM) > 0 and len(sE) > 0:
        lNm[k] = sSep.join([sM, sE])
    elif len(sS) <= 0 and len(sM) > 0 and len(sE) <= 0:
        lNm[k] = sSep.join([sM])
    elif len(sS) <= 0 and len(sM) <= 0 and len(sE) > 0:
        lNm[k] = sSep.join([sE])
    else:   # all <= 0
        lNm[k] = sSep.join([])

def cutNSpl(cS, sSep = '_', nParts = 1, fromStart = True):
    cSCt = sSep.join(cS.split(sSep)[nParts:])
    if not fromStart:
        cSCt = sSep.join(cS.split(sSep)[:(-nParts)])
    return cSCt

def cutNSplFromL(lS, sSep = '_', nParts = 1, fromStart = True):
    return [cutNSpl(cS, sSep, nParts, fromStart) for cS in lS]

def getStrFromPos(sFull, lPos = [0, 1], sSep = '_'):
    return sSep.join(sFull.split(sSep)[lPos[0]:lPos[1]])

def checkIfStrOnPos(sFull, sPart, lPos = [0, 1], sSep = '_'):
    return getStrFromPos(sFull, lPos = lPos, sSep = sSep) == sPart

def redLNames(lCurNm, lAllNm = None, lPosMove = [1, 2], sSep = '_'):
    assert len(lPosMove)  == 2
    lSL = [sSep.join(s.split(sSep)[:lPosMove[0]]) for s in lCurNm]
    lSM = [sSep.join(s.split(sSep)[lPosMove[0]:lPosMove[1]]) for s in lCurNm]
    lSR = [sSep.join(s.split(sSep)[lPosMove[1]:]) for s in lCurNm]
    lRedNm, lSMove, sMove = [s for s in lCurNm], lSM, ''
    for k in range(len(lRedNm)):
        joinStrings(lRedNm, k, sS = lSL[k], sE = lSR[k], sSep = sSep)
    if lAllNm is not None:
        for cNm in lSMove:
            if cNm not in lAllNm:
                print('ERROR: Row name part', cNm, 'not in full list', lAllNm)
                assert False
    if elEqual(lSMove) and len(lSMove) > 0:
        sMove = lSMove[0]
    else:
        print('WARNING: Elements in list', lSMove, 'unequal or list empty!')
    return lRedNm, sMove

def expLNames(lCurNm, posIns = 0, sIns = '', sSep = '_'):
    lNewNm = ['']*len(lCurNm)
    for k, cNm in enumerate(lCurNm):
        joinStrings(lNewNm, k, sS = sSep.join(cNm.split(sSep)[:posIns]),
                    sM = sIns, sE = sSep.join(cNm.split(sSep)[posIns:]),
                    sSep = sSep)
    return lNewNm

def calcSumStat(pdDfr):
    # mean, median, min, max, std.dev, variance
    return (pdDfr.mean(), pdDfr.median(), pdDfr.min(), pdDfr.max(),
            pdDfr.std(), pdDfr.var())

def allInList(lFull, lTest):
    for cET in lTest:
        if cET not in lFull:
            return False
    return True

def thingsEqual(ll, asSets = False, asLists = True):
    lAreEq = True
    if len(ll) > 1:
        if [type(l) == type(ll[0]) for l in ll[1:]] != [True]*(len(ll) - 1):
            lAreEq = False
        else:
            if ll[0] is not None:
                if asSets:
                    ll = [set(l) for l in ll]
                else:
                    if asLists:
                        ll = [list(l) for l in ll]
                refL = ll[0]
                for cL in ll[1:]:
                    if cL != refL:
                        lAreEq = False
                        break
    return lAreEq

def elEqual(l):
    return thingsEqual(l, asLists = False)

def setsEqual(ll):
    return thingsEqual(ll, asSets = True)

def listsEqual(ll):
    return thingsEqual(ll)

def dictsEqual(dDict, lKeysExcl = []):
    for cD1 in dDict.values():
        cD1R = {cK: cV for cK, cV in cD1.items() if cK not in lKeysExcl}
        for cD2 in dDict.values():
            cD2R = {cK: cV for cK, cV in cD2.items() if cK not in lKeysExcl}
            if cD1R.keys() != cD2R.keys():
                print('Lists of keys are not equal!', cD1R.keys(), cD2R.keys())
                return False
            if cD1R != cD2R:
                print('Dictionaries are not equal!', cD1R, cD2R)
                return False
    return True

def allTrue(cIt):
    return [cB for cB in cIt] == [True]*len(cIt)

def dfrEqual(lDfr):
    if len(lDfr) > 1:
        return allTrue([lDfr[0].equals(dfr) for dfr in lDfr[1:]])
    else:
        return True

def get1stFromL(l, retVDef = None):
    if len(l) > 0:
        if elEqual(l):
            return l[0]
        else:
            print('WARNING: Elements in list', l, 'unequal!')
            return retVDef
    else:
        print('WARNING: List', l, 'is empty!')
        return retVDef

def get1stFromLDict(lD, cK):
    lV = []
    if len(lD)> 0:
        if cK in lD[0]:
            lV = copy.deepcopy(lD[0][cK])
            if not listsEqual([cD[cK] for cD in lD]):
                print('ERROR: Value for key', cK, 'unequal between dicts!')
                print([cD[cK] for cD in lD])
                assert False
        else:
            print('WARNING:', cK, 'not in dictionary', lD[0])
    return lV

def get1stFromLDfr(lDfr, retVDef = None):
    if len(lDfr) > 0:
        if dfrEqual(lDfr):
            return lDfr[0]
        else:
            print('WARNING: Elements in list of DataFrames', lDfr, 'unequal!')
            return retVDef
    else:
        print('WARNING: List of DataFrames', lDfr, 'is empty!')
        return retVDef

def printErrAndStop(lSE, s1, s2):
    assert len(lSE) >= 3
    print(lSE[0], lSE[1], s1, lSE[2], s2)
    assert False

def storeResDatInDfr(dRes, lSHdr):
    assert len(lSHdr) == 3
    pdDfr = pd.DataFrame(np.zeros((len(dRes), len(lSHdr))), columns = lSHdr)
    for i, (cK, cV) in enumerate(sorted(dRes.items())):
        pdDfr.iloc[i, :] = [cK[0], cK[1], cV]
    pdDfr = pdDfr.astype({pdDfr.columns[2]: 'float64'})
    return pdDfr

def getRowDfr(pdDfr, nmR, retLDef = [], nDigRnd = None):
    retL = retLDef
    if nmR in pdDfr.index:
        if nDigRnd is None:
            retL = list(pdDfr.loc[nmR, :])
        else:
            retL = list(round(pdDfr.loc[nmR, :], nDigRnd))
    return retL

def prepPdDfr(pdDfr, cID, lID, lPosMove = [0, 1]):
    if cID in lID:
        pdDfr.index = redLNames(list(pdDfr.index), lPosMove = lPosMove)[0]

def getSrtData(pF, lSrt, lAsc, lSCat = None, sepD = ',', idxCol = 0):
    pdDfr = readCSV(pF, sepD = sepD, iCol = idxCol)
    if lSCat == None:
        lSerVC = [cSer.value_counts(sort = False) for cSer in [pdDfr.index]]
        llAttr = [list(set(pdDfr.index))]
    else:
        lSer = [pdDfr.loc[:, s] for s in lSCat]
        lSerVC = [cSer.value_counts(sort = False) for cSer in lSer]
        llAttr = [list(set(pdDfr.loc[:, s])) for s in lSCat]
    lLenLAttr = [len(llAttr[k]) for k in range(len(llAttr))]
    pdDfr = pdDfr.sort_values(by = lSrt, ascending = lAsc)
    return pdDfr, lSerVC, llAttr, lLenLAttr, pdDfr.shape[0]

def reIndexSpec(d, srtMd = 'str', srtDg = 0):
    if srtMd == 'str':
        d = d.reindex(sorted(d.columns, key = lambda x: x[srtDg:]), axis = 1)
    elif srtMd == 'float':
        d = d.reindex(sorted(d.columns, key = lambda x: float(x[srtDg:])),
                      axis = 1)
    elif srtMd == 'int':
        d = d.reindex(sorted(d.columns, key = lambda x: int(x[srtDg:])),
                      axis = 1)
    else:
        d = d.sort_index(axis = 1)
    return d

def printElapsedTimeSim(stT, cT, sPre = 'Time'):
    # calculate and display elapsed time
    elT = round(cT - stT, GC.R04)
    print(sPre, 'elapsed:', elT, 'seconds, this is', round(elT/60, GC.R04),
          'minutes or', round(elT/3600, GC.R04), 'hours or',
          round(elT/(3600*24), GC.R04), 'days.')

def showElapsedTime(startTime):
    print('-'*80)
    printElapsedTimeSim(startTime, time.time(), 'Time')
    print('+'*3 + ' Current time:', time.ctime(time.time()), '+'*3)
    print('-'*80)

def endSimu(startTime):
    print('-'*80)
    printElapsedTimeSim(startTime, time.time(), 'Total time')
    print('*'*20 + ' DONE', time.ctime(time.time()), '*'*20)

###############################################################################

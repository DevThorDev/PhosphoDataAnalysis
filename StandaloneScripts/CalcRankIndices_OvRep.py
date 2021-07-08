# -*- coding: utf-8 -*-
###############################################################################
# --- CalcRankIndices_OvRep.py ------------------------------------------------
###############################################################################
import os, time

import pandas as pd

# --- CONSTANTS ---------------------------------------------------------------
S_USC = '_'
S_DOT = '.'
S_CSV = 'csv'

S_O = '0'
S_1 = '1'
S_2 = '2'
S_3 = '3'
S_4 = '4'
S_5 = '5'

S_BIN_CODE_S = 'BC'
S_BIN_CODE_L = 'BinCode'
S_BIN_CODE_S_2 = S_BIN_CODE_S + S_2
S_BIN_CODE_L_2 = S_BIN_CODE_L + S_2

S_GT0 = 'GT' + S_O
S_GT1 = 'GT' + S_1
S_GT5 = 'GT' + S_5

S_MC = 'MC'
S_IC = 'IC'

S_MC0 = S_MC + S_O
S_MC1 = S_MC + S_1
L_S_MC = [S_MC0, S_MC1]

S_IC0 = S_IC + S_O
S_IC1 = S_IC + S_1
L_S_IC = [S_IC0, S_IC1]

L_S_CTRV = ['Sign', 'WtsA', 'WtsB']
L_S_GT = [S_GT0, S_GT1, S_GT5]

L_S_HDR_MC = [s1 + S_USC + s2 + S_USC + s3 for s1 in L_S_MC
              for s2 in L_S_CTRV for s3 in L_S_GT]
L_S_HDR_IC_MET = [s1 + S_USC + s2 + S_USC + s3 for s1 in L_S_IC
                  for s2 in L_S_CTRV for s3 in L_S_GT]

S_GENERAL = 'General'
S_I_ORIG = 'IdxOrig'
S_RI = 'RI'
S_SEL = 'Selections'
S_MET_S = 'Met'
S_MET_L = 'Metabolite'
S_MET_D = 'MetD'
S_PHO_D = 'PhoD'
S_BIN_OP = 'BinOp'
S_BIN_CODES_L = S_BIN_CODE_L + 's'

S_IC_MET = S_IC + S_USC + S_MET_S
S_IC_BIN_CODE_S_2 = S_IC + S_USC + S_BIN_CODE_S_2

R04 = 4

# --- INPUT -------------------------------------------------------------------
lTpDat = [S_MC, S_IC_MET]

sFInp_MC = 'OvRep_PhoD_GTX_AllD__BinCode2_1_832_MeanConc01_Bon_pValFOver_p'
sFOut_MC = S_SEL + S_PHO_D + S_USC + S_BIN_CODES_L + S_USC + S_MC

sFInp_ICMet = 'OvRep_BinOp_MetD_GTX_AllD_PhoD_GTX_AllD__MetD_1_13986_IC01_Bon_pValFOver_p'
sFOut_ICMet = S_SEL + S_BIN_OP + S_USC + S_MET_D + S_USC + S_IC

pCSV = os.path.join('..', '..', '..', '25_Papers', '01_FirstAuthor',
                    '04_SysBio_DataAnalysis', '80_ResultsCSV')
sSep = ';'

# --- INPUT DICTIONARY --------------------------------------------------------
dInput = {S_GENERAL: {'lTpDat': lTpDat,
                      'pCSV': pCSV,
                      'sSep': sSep},
          S_MC: {'sFInp': sFInp_MC + S_DOT + S_CSV,
                 'sFOut': sFOut_MC + S_DOT + S_CSV,
                 'lSHdr': L_S_HDR_MC,
                 'sHdrRef': S_BIN_CODE_L_2},
          S_IC_MET: {'sFInp': sFInp_ICMet + S_DOT + S_CSV,
                     'sFOut': sFOut_ICMet + S_DOT + S_CSV,
                     'lSHdr': L_S_HDR_IC_MET,
                     'sHdrRef': S_MET_L}}
for cTpDat in dInput[S_GENERAL]['lTpDat']:
    dInput[cTpDat]['pFInp'] = os.path.join(pCSV, dInput[cTpDat]['sFInp'])
    dInput[cTpDat]['pFOut'] = os.path.join(pCSV, dInput[cTpDat]['sFOut'])

# --- FUNCTIONS ---------------------------------------------------------------
def addToDictD(cD, cKMain, cKSub, cV):
    if cKMain in cD:
        assert cKSub not in cD[cKMain]
        cD[cKMain][cKSub] = cV
    else:
        cD[cKMain] = {}
        cD[cKMain][cKSub] = cV

def getLVals(pdDfr, sHd, nEl=1):
    lVCC, lVRICC = pdDfr.loc[:, sHd].to_list(), list(range(1, nEl + 1))
    minCol = min(lVCC)
    iMin = lVCC.index(minCol)
    lVRICC[iMin:] = [iMin + 1]*(len(lVCC) - iMin)
    assert min(lVRICC) >= 1
    maxVRICC = max(lVRICC)
    for i in range(nEl):
        lVRICC[i] /= maxVRICC
    return lVRICC

def transcrDict2Dfr(cD, cDfr, lSHdrIni):
    cDfrRes = cDfr.loc[:, lSHdrIni]
    for cDSub in cD.values():
        for sHdr in cDSub.values():
            cDfrRes[sHdr] = cDfr.loc[:, sHdr]
    return cDfrRes

def calcRIs(pdDfr, lSHdr, sHdrRef, sIOrig=S_I_ORIG, sUSC=S_USC):
    assert pdDfr.columns.to_list()[:(len(lSHdr) + 1)] == [sHdrRef] + lSHdr
    dRI, nL = {}, pdDfr.shape[0]
    pdDfr[sIOrig] = pdDfr.index.to_list()
    for sHdr in lSHdr:
        pdDfr.sort_values(by=[sHdr, sIOrig], ascending=[False, True],
                          inplace=True, kind='stable', ignore_index=True)
        pdDfr.loc[:, sHdr] = getLVals(pdDfr, sHdr, nEl=nL)
        sSpec, sTp, sGT = sHdr.split(sUSC)
        addToDictD(dRI, cKMain=(sSpec, sGT), cKSub=sTp, cV=sHdr)
    return transcrDict2Dfr(dRI, pdDfr, [sIOrig, sHdrRef])

def loopInpDataFrames(dInp):
    for cTpDat in dInp[S_GENERAL]['lTpDat']:
        cDfrV = pd.read_csv(dInp[cTpDat]['pFInp'], sep=dInp[S_GENERAL]['sSep'],
                            dtype={dInp[cTpDat]['sHdrRef']: str})
        cDfrR = calcRIs(cDfrV, lSHdr=dInp[cTpDat]['lSHdr'],
                        sHdrRef=dInp[cTpDat]['sHdrRef'])
        cDfrR.sort_values(by=S_I_ORIG, ascending=True, inplace=True,
                          kind='stable', ignore_index=True)
        cDfrR.to_csv(dInp[cTpDat]['pFOut'], sep=dInp[S_GENERAL]['sSep'])

def printElapsedTimeSim(stT, cT, sPre = 'Time'):
    # calculate and display elapsed time
    elT = round(cT - stT, R04)
    print(sPre, 'elapsed:', elT, 'seconds, this is', round(elT/60, R04),
          'minutes or', round(elT/3600, R04), 'hours or',
          round(elT/(3600*24), R04), 'days.')

# --- MAIN --------------------------------------------------------------------
startTime = time.time()
print('+'*25 + ' START', time.ctime(startTime), '+'*25)

loopInpDataFrames(dInput)

print('-'*80)
printElapsedTimeSim(startTime, time.time(), 'Total time')
print('+'*37, 'DONE.', '+'*36)

###############################################################################

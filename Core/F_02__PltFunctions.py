# -*- coding: utf-8 -*-
###############################################################################
# --- F_02__PltFunctions.py ---------------------------------------------------
###############################################################################
import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
import statsmodels.formula.api as smf
import statsmodels.graphics.regressionplots as regplt

import Core.C_00__GenConstants as GC
import Core.F_00__GenFunctions as GF
import Core.F_01__OSpFunctions as SF

# --- Functions (general) -----------------------------------------------------
def pltXYAxis(dfr, nmCX = None, nmCY = None, pltAxXY = (True, True)):
    minDfr, maxDfr = min(0, dfr.stack().min()), dfr.stack().max()
    if pltAxXY[0]:
        if nmCX is not None:
            minX, maxX = min(0, min(dfr.loc[:, nmCX])), max(dfr.loc[:, nmCX])
            plt.plot([minX, maxX], [0, 0], lw = 0.75, color = 'k')
        else:
            plt.plot([0, dfr.shape[0] - 1], [0, 0], lw = 0.75, color = 'k')
    if pltAxXY[1]:
        if nmCY is not None:
            minY, maxY = min(0, min(dfr.loc[:, nmCY])), max(dfr.loc[:, nmCY])
            plt.plot([0, 0], [minY, maxY], lw = 0.75, color = 'k')
        else:
            plt.plot([0, 0], [minDfr, maxDfr], lw = 0.75, color = 'k')

def pltHalfAx(dPlt, cDfr, yBase=None, pltAxX=True, pltAxY=True):
    if pltAxX:
        plt.plot([0, cDfr.shape[0] - 1], [yBase, yBase], lw=dPlt['lnWdthX'],
                 ls=dPlt['lnStyX'], color=dPlt['lnClrX'])
    if pltAxY:
        minY, maxY = cDfr.stack().min(), cDfr.stack().max()
        if yBase is not None:
            minY = yBase + dPlt['offsMinY']
        plt.plot([0., 0.], [minY, maxY], lw=dPlt['lnWdthY'],
                 ls=dPlt['lnStyY'], color=dPlt['lnClrY'])

def pltAnchoredTxt(dPlt, cAx, sTxt='', cGT=GC.NM_GT0):
    dProp = dPlt['dPropTxt'][cGT]
    cAnTxt = AnchoredText(sTxt, bbox_to_anchor=dPlt['posTxtXY'],
                          bbox_transform=plt.gca().transAxes,
                          loc=dPlt['locTxtXY'], pad=dPlt['padTxt'],
                          borderpad=dPlt['borderPadTxt'],
                          prop=dProp, frameon=dPlt['drawFrame'])
    cAx.add_artist(cAnTxt)

def decorateSavePlot(pF, dfr = pd.DataFrame(), sTtl = None, xLbl = None,
                     yLbl = None, xLim = None, yLim = None, nmCX = None,
                     nmCY = None, pltAxXY = (False, False)):
    assert len(pltAxXY) >= 2
    if dfr.size > 0:
        pltXYAxis(dfr, nmCX, nmCY, pltAxXY = pltAxXY)
    if sTtl is not None:
        plt.title(sTtl)
    if xLbl is not None:
        plt.xlabel(xLbl)
    if yLbl is not None:
        plt.ylabel(yLbl)
    if xLim is not None:
        assert len(xLim) == 2
        plt.xlim(xLim)
    if yLim is not None:
        assert len(yLim) == 2
        plt.ylim(yLim)
    plt.savefig(pF)

def setXLim(xLim = (None, None)):
    assert len(xLim) >= 2
    if xLim[0] is None:
        if xLim[1] is not None:
            plt.xlim(top = xLim[1])
    else:
        if xLim[1] is None:
            plt.xlim(bottom = xLim[0])
        else:
            plt.xlim(xLim)

def setYLim(yLim = (None, None)):
    assert len(yLim) >= 2
    if yLim[0] is None:
        if yLim[1] is not None:
            plt.ylim(top = yLim[1])
    else:
        if yLim[1] is None:
            plt.ylim(bottom = yLim[0])
        else:
            plt.ylim(yLim)

def decorateSaveFigLegOut(pF, cFig, dfr = None, sTtl = None, xLbl = None,
                          yLbl = None, xLim = None, yLim = None, nmCX = None,
                          nmCY = None, cLeg = None, pltAxXY = (False, False)):
    assert len(pltAxXY) >= 2
    if dfr is not None:
        pltXYAxis(dfr, nmCX, nmCY, pltAxXY = pltAxXY)
    if sTtl is not None:
        plt.title(sTtl)
    if xLbl is not None:
        plt.xlabel(xLbl)
    if yLbl is not None:
        plt.ylabel(yLbl)
    setXLim(xLim)
    setYLim(yLim)
    if cLeg is not None:
        cFig.savefig(pF, bbox_extra_artists = (cLeg,), bbox_inches = 'tight')
    else:
        cFig.savefig(pF)

def pltDfrS(cDfr, pF = 'Hugo.pdf', tpMark = 'x', szMark = 5, ewMark = 2,
            ecMark = (0.95, 0., 0.), fcMark = (0.9, 0.45, 0.), styLn = 'solid',
            wdthLn = 1, colLn = 'b', sTtl = '', xLbl = '', yLbl = ''):
    if not os.path.isfile(pF):
        cFig = plt.figure()
        plt.plot(cDfr, marker = tpMark, ms = szMark, mew = ewMark,
                 mec = ecMark, mfc = fcMark, ls = styLn, lw = wdthLn,
                 color = colLn)
        decorateSavePlot(pF, cDfr, sTtl, xLbl, yLbl, pltAxXY = (False, True))
        plt.close()

def pltDfrMultClr(cDfr, pF = 'Hugo.pdf', tpMark = 'x', szMark = 5, ewMark = 2,
                  ecMark = (0., 0., 0.95), fcMark = (0.3, 0.3, 0.9),
                  styLn = 'solid', wdthLn = 1, colLn = 'RG', sTtl = '',
                  xLbl = '', yLbl = '', lClr = None, pltAxXY = (False, True)):
    nCol = cDfr.shape[1]
    if not os.path.isfile(pF):
        cFig = plt.figure()
        for k in range(nCol):
            ecM, fcM = ecMark, fcMark
            if lClr is not None:
                cClr = lClr[k % len(lClr)]
            else:
                cClr = (0.9*(1. - k/nCol), 0.9*(k/nCol), 0.2*(1. - k/nCol))
                if colLn == 'RG':
                    ecM, fcM = ecMark, fcMark
                    cClr = (0.9*(1. - k/nCol), 0.9*(k/nCol), 0.2*(1. - k/nCol))
                elif colLn == 'RB':
                    ecM, fcM = (0., 0.95, 0.), (0.3, 0.9, 0.3)
                    cClr = (0.9*(1. - k/nCol), 0.2*(1. - k/nCol), 0.9*(k/nCol))
                elif colLn == 'GB':
                    ecM, fcM = (0.95, 0., 0.), (0.9, 0.3, 0.3)
                    cClr = (0.2*(1. - k/nCol), 0.9*(1. - k/nCol), 0.9*(k/nCol))
            plt.plot(cDfr.iloc[:, k], marker = tpMark, ms = szMark,
                     mew = ewMark, mec = ecM, mfc = fcM, ls = styLn,
                     lw = wdthLn, color = cClr, label = cDfr.columns[k])
        plt.legend(loc = 'best')
        decorateSavePlot(pF, cDfr, sTtl, xLbl, yLbl, pltAxXY = pltAxXY)
        plt.close()

# --- Functions (O_01__ExpData) -----------------------------------------------
def getTitleHistSelC(dITp, dOIn, addTxt = 'AllData'):
    sTtl = dOIn['cGT'] + ' Histogram of ' + dOIn['idO']
    if len(dOIn['nmFE']) > 0:
        sTtl += (' (' + dOIn['nmFE'] + ')')
    if len(addTxt) > 0:
        sTtl += ' - ' + addTxt
    return sTtl

def SUB_pltHistSelC(cDfr, pF, lClr, sTtl = '', xLb = '', hAlp = 0.5,
                    nBns = 10, lSelC = []):
    if not os.path.isfile(pF):
        plt.figure()
        if len(lSelC) > 0:
            cDfr = cDfr.loc[:, lSelC]
        cDfr.plot.hist(alpha = hAlp, bins = nBns, color = lClr)
        decorateSavePlot(pF, dfr = cDfr, sTtl = sTtl, xLbl = xLb)
        plt.close()

def pltHistSelC(dITp, dOIn, cDfr, pRF, lClr, addS = 'AllData', llSelC = None):
    hAlp, nBns, sHPre = dITp['histAlpha'], dITp['nBins'], dITp['nmHist'] + '__'
    sTtl, llOClr = getTitleHistSelC(dITp, dOIn, addS), dITp['llIOffsClr']
    if llSelC is None:     # plot all columns
        assert max(GF.flattenIt(llOClr)) < len(lClr)
        lClrS = [lClr[k] for k in GF.flattenIt(llOClr)]
        pF = GF.adaptPF4Plot(pRF, dITp['pRelPltF'], sPre = sHPre)
        SUB_pltHistSelC(cDfr, pF, lClrS, sTtl, addS, hAlp, nBns)
    else:                   # plot selected columns (len(llSelC) plots)
        for i, lSelC in enumerate(llSelC):
            pF = GF.adaptPF4Plot(pRF, dITp['pRelPltF'], sPre = sHPre,
                                 sPost = '__' + str(i))
            assert i < len(llOClr)
            lClrS = [lClr[llOClr[i][j]] for j in range(len(llOClr[i]))]
            SUB_pltHistSelC(cDfr, pF, lClrS, sTtl, addS, hAlp, nBns, lSelC)

# --- Functions (O_31__CombData) ----------------------------------------------

# --- Functions (O_41__ExpDataX) ----------------------------------------------

# --- Functions (O_81__BinaryOps) ---------------------------------------------
def getTitleHist1C(dITp, dOIn, useMs = False, sSep = '_'):
    sTtl = ('Genotype ' + dOIn['tGT'][0] + sSep + dOIn['tGT'][1] +
              ': Histogram of ' + dITp['sCorrL'])
    if len(dOIn['nmFE']) > 0:
        sTtl += (' (' + dOIn['nmFE'] + ')')
    if useMs:
        sTtl += ' - MEANS'
    return sTtl

def pltHist1C(dITp, dOIn, cDfr, pF, useMns = False):
    nmCr, hAlp, nBns = dITp['sCorrV'], dITp['histAlpha'], dITp['nBins']
    sTtl, xLm = getTitleHist1C(dITp, dOIn, useMns), dITp['histXLim']
    pF = GF.adaptPF4Plot(pF, dITp['pRelPltF'], sPre = dITp['nmHist'])
    if not os.path.isfile(pF):
        plt.figure()
        if nmCr in cDfr.columns:
            if dITp['cmpND']:
                s_SD, s_M_SD = 'Normal_SD', 'Normal_M_SD'
                aN_SD = np.random.normal(0., cDfr[nmCr].std(), cDfr.shape[0])
                aN_M_SD = np.random.normal(cDfr[nmCr].mean(), cDfr[nmCr].std(),
                                           cDfr.shape[0])
                dfrN_SD = pd.DataFrame(aN_SD, columns = [s_SD])
                dfrN_M_SD = pd.DataFrame(aN_M_SD, columns = [s_M_SD])
                lDfr = [cDfr[nmCr], dfrN_SD[s_SD], dfrN_M_SD[s_M_SD]]
                pltDfr = pd.concat(lDfr, axis = 1)
                pltDfr.plot.hist(alpha = hAlp, bins = nBns)
            else:
                cDfr[nmCr].plot.hist(alpha = hAlp, bins = nBns)
        else:
            print('ERROR:', nmCr, 'not in column names:', list(cDfr.columns))
            assert False
        decorateSavePlot(pF, sTtl = sTtl, xLbl = nmCr, xLim = xLm)
        plt.close()

def getTitlePltCorr(dITp, cDfr, i):
    sT = dITp['sCorrTtl'] + ' = ' + str(round(cDfr.loc[i, dITp['sCorrV']], 3))
    sT += ' / ' + dITp['sSpearTtl'] + ' = '
    return sT + str(round(cDfr.loc[i, dITp['sSpearV']], 3))

def SUB_pltCorr(dfr1, dfr2, pF, lOD, nmC1, nmC2, cFml, sTtl, dMark, lSXY):
    # create a DataFrame for fitting, and fit the regression line
    serX, serY = dfr1.loc[:, nmC1], dfr2.loc[:, nmC2]
    fitDfr = pd.DataFrame([serX, serY], index = lSXY).T
    cLinM = smf.ols(cFml, fitDfr).fit()
    # plot the regression line
    cFig = regplt.abline_plot(model_results = cLinM)
    cAx = cFig.axes[0]
    # plot the data in multiple colours, according to category
    cAx.set_title(sTtl)
    try:
        cAx.set_xlim((min(np.floor(serX.dropna())),
                      max(np.ceil(serX.dropna()))))
    except:
        print('Cannot set x-limits for series', list(serX))
    try:
        cAx.set_ylim((min(np.floor(serY.dropna())),
                      max(np.ceil(serY.dropna()))))
    except:
        print('Cannot set y-limits for series', list(serY))
    for i, (cHd, lNmR) in enumerate(SF.getDMap(lOD).items()):
        sLg = cHd.split(GC.S_USC)[0]
        serX, serY = dfr1.loc[lNmR, nmC1], dfr2.loc[lNmR, nmC2]
        cAx.plot(serX, serY, ls='', marker=dMark['Symbol'][i],
                 mew=dMark['EdgeWidth'][i], mec=dMark['EdgeClr'][i],
                 mfc=dMark['FaceClr'][i], ms=dMark['Size'][i], label=sLg)
    cAx.legend(loc = 'best')
    cAx.set_xlabel(nmC1)
    cAx.set_ylabel(nmC2)
    cFig.savefig(pF)
    plt.close()

def getSelDfr(cDfr, nmCr, lBd, uBd, outBd=False):
    if lBd is not None:
        if uBd is not None:
            if outBd:
                selDfr = cDfr[(cDfr[nmCr] <= lBd) | (cDfr[nmCr] >= uBd)]
            else:
                selDfr = cDfr[(cDfr[nmCr] >= lBd) & (cDfr[nmCr] <= uBd)]
        else:
            if outBd:
                selDfr = cDfr[(cDfr[nmCr] <= lBd)]
            else:
                selDfr = cDfr[(cDfr[nmCr] >= lBd)]
    else:
        if uBd is not None:
            if outBd:
                selDfr = cDfr[(cDfr[nmCr] >= uBd)]
            else:
                selDfr = cDfr[(cDfr[nmCr] <= uBd)]
        else:
            selDfr = cDfr
    return selDfr

def prepPltCorr(dITp, cDfr, pF, lODat, lBd, uBd, outBd=False):
    assert len(lODat) >= 2
    if lBd is not None and uBd is not None:
        assert lBd <= uBd
    nmCr, nmO1, nmO2 = dITp['sCorrV'], dITp['nmObj1'], dITp['nmObj2']
    dfrO1, dfrO2, sX, sY = lODat[0].cDfr, lODat[1].cDfr, 'X', 'Y'
    selDfr, cFml = cDfr, sY + ' ~ ' + sX
    if nmCr in cDfr.columns:
        selDfr = getSelDfr(cDfr, nmCr, lBd=lBd, uBd=uBd, outBd=outBd)
        for i, cRD in selDfr.iterrows():
            nmO1C, nmO2C = selDfr.loc[i, nmO1], selDfr.loc[i, nmO2]
            sTtl = getTitlePltCorr(dITp, selDfr, i)
            sNCr = '__' + '0'*(7 - len(str(i + 1))) + str(i + 1)
            pFN = GF.adaptPF4Plot(pF, dITp['pRelPltF'], sPost = sNCr)
            if not os.path.isfile(pFN):
                SUB_pltCorr(dfrO1, dfrO2, pFN, lODat, nmO1C, nmO2C, cFml, sTtl,
                            dITp['dMarker'], [sX, sY])

def pltCorr(dITp, dOIn, cDfr, pF, lODat, tPltI):
    dPltCorr, dCorrBnd, bBd = tPltI
    if dOIn['tGT'] in dPltCorr:
        if dPltCorr[dOIn['tGT']]:
            lowBd, upBd = dCorrBnd[dOIn['tGT']][:2]
            prepPltCorr(dITp, cDfr, pF, lODat, lBd=lowBd, uBd=upBd, outBd=bBd)

# --- Functions (O_82__Clustering) --------------------------------------------
def plt1DDatS(dITp, dIPlt, cDfr, pF, tInf, nmCX, pltAxXY = (True, True)):
    s1Pre, nmCY = dITp['nmPlt_' + tInf[0]] + '__', tInf[1]
    pF = GF.adaptPF4Plot(pF, dITp['pRelPltF'], sPre = s1Pre)
    if not os.path.isfile(pF):
        cFig = plt.figure()
        plt.plot(nmCX, nmCY, data = cDfr, marker = dIPlt['tpMark'],
                 ms = dIPlt['szMark'], mew = dIPlt['ewMark'],
                 mec = dIPlt['ecMark'], mfc = dIPlt['fcMark'],
                 ls = dIPlt['styLn'], lw = dIPlt['wdthLn'],
                 color = dIPlt['colLn'])
        decorateSavePlot(pF, cDfr, dIPlt['title'], dIPlt['xLbl'],
                         dIPlt['yLbl'], nmCX = nmCX, nmCY = nmCY,
                         pltAxXY = pltAxXY)
        plt.close()

def pltClCent(dITp, dIPlt, cTrD, dClDfr, pF, lClr = None,
              pltAxXY = (True, False)):
    sTtl, xLbl, yLbl = dIPlt['title'], dIPlt['xLbl'], cTrD + dIPlt['yLbl']
    tpMark, szMark, ewMark = dIPlt['tpMark'], dIPlt['szMark'], dIPlt['ewMark']
    ecMark, fcMark, styLn = dIPlt['ecMark'], dIPlt['fcMark'], dIPlt['styLn']
    wdthLn, lClr, pltAxXY = dIPlt['wdthLn'], dIPlt['lClr'], (True, False)
    xLim = (dIPlt['xLimB'], dIPlt['xLimT'])
    yLim = (dIPlt['yLimB'], dIPlt['yLimT'])
    for cNCl, cDfr in dClDfr.items():
        sCPre = dITp['nmPlt_ClCent'] + '__' + str(cNCl) + 'Cl__'
        pFN = GF.adaptPF4Plot(pF, dITp['pRelPltF'], sPre = sCPre)
        if not os.path.isfile(pFN):
            nCol = cDfr.shape[1]
            cFig, cAx = plt.subplots()
            for k in range(nCol):
                ecM, fcM, cClr = ecMark, fcMark, lClr[k % len(lClr)]
                cAx.plot(cDfr.iloc[:, k], marker = tpMark, ms = szMark,
                         mew = ewMark, mec = ecM, mfc = fcM, ls = styLn,
                         lw = wdthLn, color = cClr, label = cDfr.columns[k])
            l = cAx.legend(loc = 'center', bbox_to_anchor = dIPlt['posLegXY'])
            decorateSaveFigLegOut(pFN, cFig, cDfr, sTtl, xLbl, yLbl, xLim,
                                  yLim, cLeg = l, pltAxXY = pltAxXY)
            plt.close()

# --- Functions (O_83__OverRep) -----------------------------------------------
def plotProfile(dITp, cDfr, thrD, tGT, pF, k = 0, tpPr = 'PD'):
    i, j = 0, 0
    while j < len(cDfr.columns):
        d = cDfr.iloc[:, j:(j + dITp['jIncr'])]
        i += dITp['iIncr']
        j += dITp['jIncr']
        pFN = GF.adaptPF4Plot(pF, dITp['pRelPltF'],
                              sPre = dITp['nmPlt_Prof'] + str(i) + '_')
        if not os.path.isfile(pFN):
            cFig, cAx = plt.subplots()
            for sC in d.columns:
                if tpPr == 'PD' and sC in dITp['dClrBinC']:
                    cAx.plot(d.loc[:, sC], lw = 0.75, label = sC,
                             color = dITp['dClrBinC'][sC])
                else:
                    cAx.plot(d.loc[:, sC], lw = 0.75, label = sC)
            cAx.set_xlabel(dITp['dTpX'][tpPr])
            cAx.set_ylabel(dITp['dTpY'][dITp['lTpY'][k]][0])
            l = cAx.legend(loc = 'center', bbox_to_anchor = dITp['posLegXY'],
                           fontsize = dITp['szFontLeg'])
            pltHalfAx(dITp, d, yBase=thrD, pltAxX=True, pltAxY=True)
            sGT, sNmGT = SF.getGT(tGT, dITp['dNmGT'], iPos=3)
            pltAnchoredTxt(dITp, cAx, sTxt=sNmGT, cGT=sGT)
            if l is not None:
                cFig.savefig(pFN, bbox_extra_artists = (l,),
                             bbox_inches = 'tight')
            else:
                cFig.savefig(pFN)
            plt.close()

###############################################################################

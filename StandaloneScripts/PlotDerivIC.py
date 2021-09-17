# -*- coding: utf-8 -*-
###############################################################################
# --- PlotDerivIC.py ----------------------------------------------------------
###############################################################################
import os, time, math, itertools

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# --- CONSTANTS ---------------------------------------------------------------
S_SPACE = ' '
S_USC = '_'
S_BAR = '|'
S_SLASH = '/'
S_DOT = '.'
S_CSV = 'csv'
S_PDF = 'pdf'
S_SPSLSP = S_SPACE + S_SLASH + S_SPACE

S_M = 'M'
S_P = 'P'
S_I = 'I'
S_C = 'C'
L_S_SPC = [S_M, S_P, S_I, S_C]

S_0 = '0'
S_1 = '1'
S_2 = '2'
S_3 = '3'
S_4 = '4'
S_5 = '5'
S_6 = '6'
L_S_0_6 = [S_0, S_1, S_2, S_3, S_4, S_5, S_6]

S_MET = 'Metabolite'
S_METS = S_MET + 's'
S_PHO = 'Phosphopeptide'
S_PHOS = S_PHO + 's'

S_IC = 'IC'
S_IC_MEX = '$I_C$'
S_IC_M = S_IC + '(-)'
S_IC_P = S_IC + '(+)'

S_GT = 'GT'
S_GT0 = S_GT + S_0
S_GT1 = S_GT + S_1
S_GT5 = S_GT + S_5
L_S_GT = [S_GT0, S_GT1, S_GT5]

S_NM_GT0 = 'wild type'
S_NM_GT1 = '$\it{pgm}$ mutant'
S_NM_GT5 = '$\it{sweet11/12}$ mutant'
D_S_NM_GT = {S_GT0: S_NM_GT0, S_GT1: S_NM_GT1, S_GT5: S_NM_GT5}

S_FT1 = 'DR'
S_FT2 = 'DS'
S_FT3 = 'NR'
S_FT4 = 'NS'
L_S_FT = [S_FT1, S_FT2, S_FT3, S_FT4]
L_S_FT_CHG = [L_S_FT[i] + S_BAR + L_S_FT[j] for j in range(len(L_S_FT))
              for i in range(len(L_S_FT)) if i != j]
D_S_FT_CHG = {s: [S_USC.join([t, s]) for t in L_S_FT_CHG] for s in L_S_SPC}

S_DEV_SD = 'DeviationSD'
S_DERIV = 'Deriv'
S_IC_DERIV = S_IC + S_DERIV
S_CMB_TO_IC = 'CombineToIC'

S_CMB_C = 'Combined'
S_DEV_S = 'deviation'
S_DEVS_C = 'Deviations between features'
S_CMP_S = 'component'
S_CMB_DEV = S_CMB_C + S_SPACE + S_DEV_S
S_IC_CMP = S_IC_MEX + S_SPACE + S_CMP_S

S_YLBL_DEV_SD_PLT = 'Difference as multiples of standard deviation'
S_YLBL_IC_DERIV_PLT = S_DEVS_C + S_SPSLSP + S_IC_CMP
S_YLBL_CMB_TO_IC_PLT = (S_IC_CMP + S_SPSLSP + S_IC_M + S_SPSLSP + S_IC_P +
                        S_SPSLSP + S_IC)

S_BASE_CL = 'BaseClass'
S_INP_DATA = 'InputData'
S_ROOT_CL = 'RootClass'
S_PLTR = 'Plotter'
S_DEV_SD_PLTR = S_DEV_SD + S_PLTR
S_IC_DERIV_PLTR = S_IC_DERIV + S_PLTR
S_CMB_TO_IC_PLTR = S_CMB_TO_IC + S_PLTR

S_PLT = 'Plot'
S_NM_DEV_SD_PLT = S_1 + S_USC + S_DEV_SD + S_PLT
S_NM_IC_DERIV_PLT = S_2 + S_USC + S_IC_DERIV + S_PLT
S_NM_CMB_TO_IC_PLT = S_3 + S_USC + S_CMB_TO_IC + S_PLT

R04 = 4

# --- INPUT -------------------------------------------------------------------
# --- flow control ------------------------------------------------------------
doPlot_DevSD = False             # True / False
doPlot_ICDrv = False             # True / False
doPlot_CmbToIC = True           # True / False

# --- data specific input -----------------------------------------------------
maxSLen = 20
sSep = ';'

# --- graphics parameters / all plots -----------------------------------------
clrDef = 'k'                    # default colour
szFontLeg = 'small'             # font size of legend
nCharDsp = 60                   # number of chars displayed for legend item
coordAnchorBox = (.5, 1.02)     # coordinates of the legend anchor box

# --- graphics parameters / deviation SD plot ---------------------------------
dTupIn_DevSD = {'A': (S_GT0, 'Alanine', S_FT2 + S_BAR + S_FT1),
                'B': (S_GT0, 'Alanine', S_FT2 + S_BAR + S_FT4),
                'C': (S_GT1, 'Docosanoic_acid', S_FT4 + S_BAR + S_FT1),
                'D': (S_GT1, 'Docosanoic_acid', S_FT1 + S_BAR + S_FT4),
                'E': (S_GT5, 'Putrescine', S_FT2 + S_BAR + S_FT1),
                'F': (S_GT5, 'Putrescine', S_FT3 + S_BAR + S_FT2)}
nmPlt_DevSD = S_NM_DEV_SD_PLT       # name prefix of the deviation SD plot
tFigSz_DevSD = (2., 4.)             # (width, height): figure size [inches]
symMark_DevSD = 'x'                 # marker symbol
szMark_DevSD = 25                   # marker size
lenMnBar_DevSD = .3                 # length of mean bar
lwdMnBar_DevSD = 1.5                # line width of mean bar
dltLwdMnBar_DevSD = .5              # delta of width of black line to col. line
lenWMnBar_DevSD = .05               # length of mean bar whiskers
lwdWMnBar_DevSD = .6                # line width of mean bar whiskers
lwdLnSD_DevSD = .8                  # line width of mean bar
dltLwdLnSD_DevSD = .25              # delta of width of black line to col. line
lenWLnSD_DevSD = .1                 # length of mean bar whiskers
lwdWLnSD_DevSD = .4                 # line width of mean bar whiskers
lstCnLnSD_DevSD = '--'              # line style of connecting line
lwdCnLnSD_DevSD = .8                # line width of connecting line
dPosFt_DevSD = {0: .75, 1: .25}     # dictionary of positions of features 0|1
axXLim_DevSD = (0., 1.)             # limits for the x-axis, or None
axYLim_DevSD = None                 # limits for the y-axis, or None
adaptYLim_DevSD = True              # adapt y-limits using values to plot?
dFontTtl_DevSD = {'size': 10}       # font dictionary of the title
locTtl_DevSD = 'left'               # location of the title
padTtl_DevSD = 10.                  # padding of the title
degRotXLbl_DevSD = 90               # degree rotation of x-labels

# --- graphics parameters / IC derivation plot --------------------------------
dTupIn_ICDrv = {'A': (S_GT0, 'Alanine', 'ESLS(1)PGQQHVSQNTAVKPEGR'),
                'B': (S_GT0, 'Alanine', 'T(0.001)PS(0.996)QRPS(0.003)TSSSSGGFNIGK'),
                'C': (S_GT1, 'Glutamic_acid', 'TADS(1)DGES(1)GEIKFEDNNLPPLGEK'),
                'D': (S_GT5, 'Putrescine', 'S(0.999)FS(0.001)VADFPR'),
                'E': (S_GT5, 'Fructose', 'TEEDENDDEDHEHKDDKT(0.854)S(0.144)PDS(0.001)IVMVEAK')}
nmPlt_ICDrv = S_NM_IC_DERIV_PLT     # name prefix of the IC derivation plot
tFigSz_ICDrv = None                 # (width, height): figure size [inches]
lwdPlt_ICDrv = .75                  # line width in plot
axXLim_ICDrv = None                 # limits for the x-axis, or None
axYLim_ICDrv = (-11, 11)            # limits for the y-axis, or None
locTtl_ICDrv = 'left'               # location of the title
padTtl_ICDrv = 70.                  # padding of the title
locLegend_ICDrv = 'lower center'    # location of the legend
wdthGrp_ICDrv = .8                  # width of bar group
wdthBar_ICDrv = wdthGrp_ICDrv/len(D_S_FT_CHG) - .01    # width of single bars
degRotXLbl_ICDrv = 90               # degree rotation of x-labels
plotVLines_ICDrv = True             # plot vertical lines between groups?
lwdVLine_ICDrv = .2                 # width of vertical line
clrVLine_ICDrv = clrDef             # colour of vertical line

# --- graphics parameters / combine to IC plot --------------------------------
dTupIn_CmbToIC = {'A': (S_GT0, 'Alanine', 'ESLS(1)PGQQHVSQNTAVKPEGR'),
                  'B': (S_GT0, 'Alanine', 'T(0.001)PS(0.996)QRPS(0.003)TSSSSGGFNIGK'),
                  'C': (S_GT1, 'Glutamic_acid', 'TADS(1)DGES(1)GEIKFEDNNLPPLGEK'),
                  'D': (S_GT5, 'Putrescine', 'S(0.999)FS(0.001)VADFPR'),
                  'E': (S_GT5, 'Fructose', 'TEEDENDDEDHEHKDDKT(0.854)S(0.144)PDS(0.001)IVMVEAK')}
nmPlt_CmbToIC = S_NM_CMB_TO_IC_PLT  # name prefix of the IC derivation plot
tFigSz_CmbToIC = None               # (width, height): figure size [inches]
lwdPlt_CmbToIC = .75                # line width in plot
axXLim_CmbToIC = None               # limits for the x-axis, or None
axYLim_CmbToIC = (-11, 11)          # limits for the y-axis, or None
locTtl_CmbToIC = 'left'             # location of the title
padTtl_CmbToIC = 70.                # padding of the title
locLegend_CmbToIC = 'lower center'  # location of the legend
wdthGrp_CmbToIC = .8                # width of bar group
wdthBar_CmbToIC = wdthGrp_CmbToIC/len(D_S_FT_CHG) - .01    # width of single bars
degRotXLbl_CmbToIC = 90             # degree rotation of x-labels
plotVLines_CmbToIC = True           # plot vertical lines between groups?
lwdVLine_CmbToIC = .2               # width of vertical line
clrVLine_CmbToIC = clrDef           # colour of vertical line

# --- names and paths of files and dirs ---------------------------------------
sMs, sPs, s1, s2, s3 = S_METS, S_PHOS, 'Means', 'SDs', 'Deviations'
dSFIn_Ms = {(sMs, sGT): S_USC.join([sMs, s1, s2, s3, sGT]) for sGT in L_S_GT}
dSFIn_Ps = {(sPs, sGT): S_USC.join([sPs, s1, s2, s3, sGT]) for sGT in L_S_GT}
dSFIn_DevSD = {**dSFIn_Ms, **dSFIn_Ps}
s1, s2, s3, s4, s5, s6 = 'Corr', 'BinOp', 'MetD', 'PhoD', 'DvSD', 'AllD'
dSFIn_ICDrv = {sGT: (S_IC_DERIV + S_USC*2 + s1 + S_USC*2 +
                     S_USC.join([s2, s3, s5, sGT, s6, s4, s5, sGT, s6]))
               for sGT in L_S_GT}
sFOut_DevSD = nmPlt_DevSD
sFOut_ICDrv = nmPlt_ICDrv
sFOut_CmbToIC = nmPlt_CmbToIC

sDirInCSV = '41_Inp_ICDeriv'
sDirOutPDF = '45_OutPDF_ICDeriv'

pBaseIn = os.path.join('..', '..', '..', '25_Papers', '01_FirstAuthor',
                       '04_SysBio_DataAnalysis')
pBaseOut = os.path.join('..', '..', '..', '25_Papers', '01_FirstAuthor',
                        '04_SysBio_DataAnalysis')
# pBaseOut = os.path.join('..', '..', '..', '..', '..', '..', 'V')
pInCSV = os.path.join(pBaseIn, sDirInCSV)
pOutPDF = os.path.join(pBaseOut, sDirOutPDF)

# --- derived values ----------------------------------------------------------
dPFIn_DevSD = {tK: os.path.join(pInCSV, dSFIn_DevSD[tK] + S_DOT + S_CSV)
               for tK in dSFIn_DevSD}
dPFIn_ICDrv = {sGT: os.path.join(pInCSV, dSFIn_ICDrv[sGT] + S_DOT + S_CSV)
               for sGT in L_S_GT}

# --- assertions --------------------------------------------------------------
assert len(L_S_SPC) == 4 and len(D_S_FT_CHG) == 4

# --- INPUT DICTIONARY --------------------------------------------------------
dInput = {# --- constants
          'lSGT': L_S_GT,
          'sBase': S_BASE_CL,
          'sInpDat': S_INP_DATA,
          'sRoot': S_ROOT_CL,
          'sPltr': S_PLTR,
          'sPltr_DevSD': S_DEV_SD_PLTR,
          'sPltr_ICDrv': S_IC_DERIV_PLTR,
          'sPltr_CmbToIC': S_CMB_TO_IC_PLTR,
          'sMet': S_MET,
          'sPho': S_PHO,
          'R04': R04,
          # --- flow control
          'doPlot_DevSD': doPlot_DevSD,
          'doPlot_ICDrv': doPlot_ICDrv,
          'doPlot_CmbToIC': doPlot_CmbToIC,
          # --- data specific input
          'maxSLen': maxSLen,
          'sSep': sSep,
          # --- graphics parameters / all plots
          'plot_Gen': {'clrDef': clrDef,
                       'szFontLeg': szFontLeg,
                       'nCharDsp': nCharDsp,
                       'coordAnchorBox': coordAnchorBox},
          # --- graphics parameters / deviation SD plot
          'plot_DevSD': {'dTupIn': dTupIn_DevSD,
                         'nmPlt': nmPlt_DevSD,
                         'tFigSz': tFigSz_DevSD,
                         'symMark': symMark_DevSD,
                         'szMark': szMark_DevSD,
                         'lenMnBar': lenMnBar_DevSD,
                         'lenHfMnBar': lenMnBar_DevSD/2.,
                         'lwdMnBar': lwdMnBar_DevSD,
                         'dltLwdMnBar': dltLwdMnBar_DevSD,
                         'lenWMnBar': lenWMnBar_DevSD,
                         'lwdWMnBar': lwdWMnBar_DevSD,
                         'lwdLnSD': lwdLnSD_DevSD,
                         'dltLwdLnSD': dltLwdLnSD_DevSD,
                         'lenWLnSD': lenWLnSD_DevSD,
                         'lwdWLnSD': lwdWLnSD_DevSD,
                         'lstCnLnSD': lstCnLnSD_DevSD,
                         'lwdCnLnSD': lwdCnLnSD_DevSD,
                         'dPosFt': dPosFt_DevSD,
                         'axXLim': axXLim_DevSD,
                         'axYLim': axYLim_DevSD,
                         'adaptYLim': adaptYLim_DevSD,
                         'dFontTtl': dFontTtl_DevSD,
                         'locTtl': locTtl_DevSD,
                         'padTtl': padTtl_DevSD,
                         'degRotXLbl': degRotXLbl_DevSD},
          # --- graphics parameters / IC derivation plot
          'plot_ICDrv': {'dTupIn': dTupIn_ICDrv,
                         'nmPlt': nmPlt_ICDrv,
                         'tFigSz': tFigSz_ICDrv,
                         'lwdPlt': lwdPlt_ICDrv,
                         'axXLim': axXLim_ICDrv,
                         'axYLim': axYLim_ICDrv,
                         'locTtl': locTtl_ICDrv,
                         'padTtl': padTtl_ICDrv,
                         'locLegend': locLegend_ICDrv,
                         'wdthGrp': wdthGrp_ICDrv,
                         'wdthBar': wdthBar_ICDrv,
                         'degRotXLbl': degRotXLbl_ICDrv,
                         'plotVLines': plotVLines_ICDrv,
                         'lwdVLine': lwdVLine_ICDrv,
                         'clrVLine': clrVLine_ICDrv},
          # --- graphics parameters / combine to IC plot
          'plot_CmbToIC': {'dTupIn': dTupIn_ICDrv,
                           'nmPlt': nmPlt_CmbToIC,
                           'tFigSz': tFigSz_CmbToIC,
                           'lwdPlt': lwdPlt_CmbToIC,
                           'axXLim': axXLim_CmbToIC,
                           'axYLim': axYLim_CmbToIC,
                           'locTtl': locTtl_CmbToIC,
                           'padTtl': padTtl_CmbToIC,
                           'locLegend': locLegend_CmbToIC,
                           'wdthGrp': wdthGrp_CmbToIC,
                           'wdthBar': wdthBar_CmbToIC,
                           'degRotXLbl': degRotXLbl_CmbToIC,
                           'plotVLines': plotVLines_CmbToIC,
                           'lwdVLine': lwdVLine_CmbToIC,
                           'clrVLine': clrVLine_CmbToIC},
          # --- names and paths of files and dirs
          'pInCSV': pInCSV,
          'pOutPDF': pOutPDF,
          'dSFIn_DevSD': dSFIn_DevSD,
          'dSFIn_ICDrv': dSFIn_ICDrv,
          'sFOut_DevSD': sFOut_DevSD,
          'sFOut_ICDrv': sFOut_ICDrv,
          'sFOut_CmbToIC': sFOut_CmbToIC,
          # --- further derived values
          'dPFIn_DevSD': dPFIn_DevSD,
          'dPFIn_ICDrv': dPFIn_ICDrv}

# --- FUNCTIONS ---------------------------------------------------------------
def flattenIt(cIterable, rmvNaN=True, retArr=False):
    itFlat = list(itertools.chain.from_iterable(cIterable))
    if rmvNaN:
        itFlat = [x for x in itFlat if (x is not None and not np.isnan(x))]
    if retArr:
        itFlat = np.array(itFlat)
    return itFlat

def getL2Val(midV=0., lenLn=1., cDir='both'):
    l = [midV - lenLn, midV + lenLn]
    if cDir == 'left':
        l = [midV - lenLn, midV]
    elif cDir == 'right':
        l = [midV, midV + lenLn]
    return l

def plotWskH(dPlt, cAx, xMid=0., yB=0., yT=0., lenW=.1, lwdW=.5, cDir='both'):
    lX = getL2Val(xMid, lenW, cDir=cDir)
    for cY in [yB, yT]:
        cAx.plot(lX, [cY]*2, lw=lwdW, color=dPlt['clrDef'])

def plotWskV(dPlt, cAx, xB=0., xT=0., yMid=0., lenW=.1, lwdW=.5, cDir='both'):
    lY = getL2Val(yMid, lenW, cDir=cDir)
    for cX in [xB, xT]:
        cAx.plot([cX]*2, lY, lw=lwdW, color=dPlt['clrDef'])

def plotClrLn(dPlt, cAx, lX, lY, cLwd, cClr, dltLwd=.5):
    cAx.plot(lX, lY, lw=cLwd, color=cClr)
    cAx.plot(lX, lY, lw=max(0, cLwd - 2*dltLwd), color=dPlt['clrDef'])

def plotMean(dPlt, cAx, xMid=0., yMn=0., cClr='C0'):
    lX, lY = [xMid - dPlt['lenHfMnBar'], xMid + dPlt['lenHfMnBar']], [yMn]*2
    plotClrLn(dPlt, cAx, lX, lY, cLwd=dPlt['lwdMnBar'], cClr=cClr,
              dltLwd=dPlt['dltLwdMnBar'])
    plotWskV(dPlt, cAx, xB=lX[0], xT=lX[1], yMid=yMn, lenW=dPlt['lenWMnBar'],
             lwdW=dPlt['lwdWMnBar'])

def getNumSDSteps(othYMn, yMn, ySD):
    nSDStp, sgnDiff = 0, 0
    if othYMn > yMn:
        nSDStp, sgnDiff = math.ceil((othYMn - yMn)/ySD), 1
    elif othYMn < yMn:
        nSDStp, sgnDiff = math.ceil((yMn - othYMn)/ySD), -1
    return nSDStp, math.copysign(ySD, sgnDiff)

def plotStepsSD(dPlt, cAx, tFtI, xLn=0., othYMn=0.):
    yMn, ySD, cClr = tFtI
    nSDStp, ySDSgn = getNumSDSteps(othYMn, yMn, ySD)
    for i in range(nSDStp):
        lY = [yMn + i*ySDSgn, yMn + (i + 1)*ySDSgn]
        plotClrLn(dPlt, cAx, [xLn]*2, lY, cLwd=dPlt['lwdLnSD'], cClr=cClr,
                  dltLwd=dPlt['dltLwdLnSD'])
        plotWskH(dPlt, cAx, xMid=xLn, yB=lY[0], yT=lY[1],
                 lenW=dPlt['lenWLnSD'], lwdW=dPlt['lwdWLnSD'], cDir='right')
    lX = [xLn, dPlt['dPosFt'][0] - dPlt['lenHfMnBar']]
    cAx.plot(lX, [othYMn]*2, ls=dPlt['lstCnLnSD'], lw=dPlt['lwdCnLnSD'],
             color=dPlt['clrDef'])

def saveClosePlot(cFig, pPltF, l=None):
    if l is not None:
        cFig.savefig(pPltF, bbox_extra_artists=(l,), bbox_inches='tight')
    else:
        cFig.savefig(pPltF)
    plt.close()

def printElapsedTimeSim(stT, cT, sPre = 'Time'):
    # calculate and display elapsed time
    elT = round(cT - stT, R04)
    print(sPre, 'elapsed:', elT, 'seconds, this is', round(elT/60, R04),
          'minutes or', round(elT/3600, R04), 'hours or',
          round(elT/(3600*24), R04), 'days.')

# --- CLASSES -----------------------------------------------------------------
class BaseClass():
    def __init__(self):
        self.idO = S_BASE_CL
        self.descO = 'Base class'
        print('Initiated "Base" base object.')

    def printIDDesc(self):
        print('Object ID:', self.idO)
        print('Object description:', self.descO)

    def printAttrList(self):
        lAttr = dir(self)
        print('List of attributes:')
        for cAttr in lAttr:
            print(cAttr)

    def printAttrData(self):
        print('Attributes and attribute values:')
        d = vars(self)
        for cK, cV in d.items():
            print(cK, ':\t', cV)

# .............................................................................
class InputData(BaseClass):
    def __init__(self, dInp):
        super().__init__()
        self.idO = dInp['sInpDat']
        self.descO = 'Input data class'
        self.dI = dInp
        self.fillInp()
        print('Initiated "InputData" base object.')

    def fillInp(self):
        for sK, cV in self.dI.items():
            setattr(self, sK, cV)
        print('Set InputData attributes.')

# .............................................................................
class RootClass(BaseClass):
    def __init__(self, InpD):
        super().__init__()
        self.idO = InpD.sRoot
        self.descO = 'Root class'
        self.inpD = InpD
        self.dDfrIn = None
        print('Initiated "RootClass" base object.')

    def printObjInfo(self):
        print('-'*20, 'Object', self.descO, '(ID', self.idO, ')', '-'*20)
        print('-'*8, 'Input data:')
        self.inpD.printAttrData()
        print('-'*8, 'Attributes of', self.descO, 'class:')
        self.printAttrData()

    def printDDfrInp(self, lIdxCol=None):
        if self.dDfrIn is None:
            print('Input DataFrames dictionary does not have any content yet.')
        else:
            if lIdxCol is None:
                print('Input DataFrames dictionary:')
                print(self.dDfrIn)
            else:
                print('Input DataFrames dictionaries (columns', lIdxCol, '):')
                for cK, dfrIn in self.dDfrIn.items():
                    print('-'*8, 'Dictionary corresponding to key', cK, ':')
                    print(dfrIn.iloc[:, lIdxCol])

    def printDfrInp(self, cK, printIdx=False, printCol=False):
        if cK in self.dDfrIn:
            print('DataFrame corresponding to key ' + str(cK) + ':')
            print(self.dDfrIn[cK])
            if printIdx:
                print('DataFrame index (with length ' +
                      str(self.dDfrIn[cK].index.size) + '):')
                print(self.dDfrIn[cK].index)
            if printCol:
                print('DataFrame columns (with length ' +
                      str(self.dDfrIn[cK].columns.size) + '):')
                print('DataFrame columns:')
                print(self.dDfrIn[cK].columns)
        else:
            print('Key', cK, 'not in input DataFrame dictionary!')
            print('Keys of input DataFrame dictionary:', list(self.dDfrIn))

# .............................................................................
class Plotter(RootClass):
    def __init__(self, InpD, axYLim):
        super().__init__(InpD)
        self.idO = InpD.sPltr
        self.descO = 'Class for plotting'
        self.sSp = self.inpD.sSep
        self.dPlt = self.inpD.plot_Gen
        self.dPPltF = {}
        self.setDummyVal()
        print('Initiated "Plotter" base object.')

    def printDPPltF(self):
        print('Dictionary of plot file paths:')
        for tK, pF in self.dPPltF.items():
            print(tK, ':', pF)
        print('-'*64)

    def loadDDfrInp(self, idxC=None):
        # load input DataFrames, and save them in dictionary
        if hasattr(self, 'dPFIn'):
            self.dDfrIn = {}
            for sK in self.dPFIn:
                self.dDfrIn[sK] = pd.read_csv(self.dPFIn[sK], sep=self.sSp,
                                              index_col=idxC)

    def getDPPltF(self, sFOut):
        self.dPPltF, mxL = {}, self.inpD.maxSLen
        for sID, t in self.dPlt['dTupIn'].items():
            u = S_USC.join([sFOut, sID])
            lCmpNmF = [s.replace(S_BAR, S_USC)[:mxL] for s in t]
            sPltF = S_DOT.join([S_USC.join([u] + lCmpNmF), S_PDF])
            self.dPPltF[sID] = (t, os.path.join(self.pDOut, sPltF))

    def setDummyVal(self):
        l = ['clrDef', 'szFontLeg', 'nCharDsp', 'coordAnchorBox', 'dTupIn',
             'nmPlt', 'tFigSz', 'dPosFt', 'symMark', 'szMark', 'lenMnBar',
             'lenHfMnBar', 'lwdMnBar', 'dltLwdMnBar', 'lenWMnBar', 'lwdWMnBar',
             'lwdLnSD', 'dltLwdLnSD', 'lenWLnSD', 'lwdWLnSD', 'lstCnLnSD',
             'lwdCnLnSD', 'lwdPlt', 'axXLim', 'axYLim', 'adaptYLim',
             'dFontTtl', 'locTtl', 'padTtl', 'locLegend', 'wdthGrp', 'wdthBar',
             'degRotXLbl', 'plotVLines', 'lwdVLine', 'clrVLine', 'axXTck',
             'axYTck', 'lblXTck']
        for s in l:
            if s not in self.dPlt:
                self.dPlt[s] = None

    def iniPlot(self):
        cFig, cAx = plt.subplots()
        if self.dPlt['tFigSz'] is not None and len(self.dPlt['tFigSz']) == 2:
            cFig.set_size_inches(self.dPlt['tFigSz'])
        return cFig, cAx

    def decoratePlot(self, cAx, sTtl=None, sYLbl=None):
        if self.dPlt['axXLim'] is not None:
            cAx.set_xlim(self.dPlt['axXLim'])
        if self.dPlt['axYLim'] is not None:
            cAx.set_ylim(self.dPlt['axYLim'])
        if self.dPlt['axXTck'] is not None:
            cAx.set_xticks(self.dPlt['axXTck'])
            if self.dPlt['lblXTck'] is not None:
                cAx.set_xticklabels(self.dPlt['lblXTck'])
            if self.dPlt['degRotXLbl'] is not None:
                for cXLbl in cAx.get_xticklabels():
                    cXLbl.set_rotation(self.dPlt['degRotXLbl'])
        if self.dPlt['axYTck'] is not None:
            cAx.set_yticks(self.dPlt['axYTck'])
        if sTtl is not None:
            cAx.set_title(sTtl, fontdict=self.dPlt['dFontTtl'],
                          loc=self.dPlt['locTtl'], pad=self.dPlt['padTtl'])
        if sYLbl is not None:
            cAx.set_ylabel(sYLbl)
        plt.tight_layout()

# .............................................................................
class DevSDPlotter(Plotter):
    def __init__(self, InpD):
        super().__init__(InpD, InpD.plot_DevSD['axYLim'])
        self.idO = self.inpD.sPltr_DevSD
        self.descO = 'Deviations (as multiples of feature SD) plotter'
        self.dPFIn = self.inpD.dPFIn_DevSD
        self.pDOut = self.inpD.pOutPDF
        self.dPlt = {**self.dPlt, **self.inpD.plot_DevSD}
        self.getDPPltF(sFOut=self.inpD.sFOut_DevSD)
        self.loadDDfrInp(idxC=0)
        print('Initiated "DevSDPlotter" base object and loaded input data.')

    def preProcData(self, sGT, sMP, sFtChg):
        d, dFt = self.dDfrIn[(S_METS, sGT)], {}
        for k, cFt in enumerate(sFtChg.split(S_BAR)):
            lSC = [s for s in d.columns if (s[-1] in L_S_0_6 and
                                            s.startswith(cFt))]
            cSer, cC = d.loc[sMP, lSC], 'C' + str(k)
            dFt[k] = (cFt, lSC, cSer, np.mean(cSer), np.std(cSer, ddof=1), cC)
        assert len(dFt) == 2
        return dFt

    def setTicks(self, dFt=None):
        if self.dPlt['axXLim'] is not None and len(self.dPlt['axXLim']) >= 2:
            self.dPlt['axXTck'] = [self.dPlt['dPosFt'][k] for k in dFt]
            self.dPlt['lblXTck'] = [dFt[k][0] for k in dFt]
        if self.dPlt['adaptYLim']:
            lAllV = flattenIt([t[2] for t in dFt.values()])
            mnV, mxV, cStep = math.floor(min(lAllV)), math.ceil(max(lAllV)), 1
            self.dPlt['axYLim'] = (mnV, mxV)
            self.dPlt['axYTck'] = range(mnV, mxV + 1, cStep)

    def createPlot(self, dFtI):
        self.setTicks(dFt=dFtI)
        cFig, cAx = self.iniPlot()
        for k, (cFt, lHdC, cSer, cMn, cSD, cClr) in dFtI.items():
            lX = [self.dPlt['dPosFt'][k]]*cSer.size
            cAx.scatter(lX, cSer, marker=self.dPlt['symMark'],
                        s=self.dPlt['szMark'])
            plotMean(self.dPlt, cAx, lX[0], cMn, cClr=cClr)
            if k == 1 and cSD > 0:
                xPLn = lX[0] + self.dPlt['lenHfMnBar']
                plotStepsSD(self.dPlt, cAx, tFtI=dFtI[k][3:], xLn=xPLn,
                            othYMn=dFtI[0][3])
        return cFig, cAx

    def decoratePlot(self, cAx, sGT, sMP='', sYLbl=None):
        sTtl = sMP + '\n(' + D_S_NM_GT[sGT] + ')'
        super().decoratePlot(cAx, sTtl=sTtl, sYLbl=sYLbl)

    def plotDevSD(self):
        for sID, ((sGT, sMP, sFtChg), pPltF) in self.dPPltF.items():
            print('Plotting deviation in multiples of SD: "' + str(sID) +
                  '" for "' + sGT + '", substance ' + sMP + ', feature change',
                  sFtChg, '...')
            # if not os.path.isfile(pPltF):
            cFig, cAx = self.createPlot(self.preProcData(sGT, sMP, sFtChg))
            self.decoratePlot(cAx, sGT=sGT, sMP=sMP, sYLbl=S_YLBL_DEV_SD_PLT)
            saveClosePlot(cFig, pPltF)

# .............................................................................
class ICDerivPlotter(Plotter):
    def __init__(self, InpD):
        super().__init__(InpD, InpD.plot_ICDrv['axYLim'])
        self.idO = self.inpD.sPltr_ICDrv
        self.descO = 'Concordance index derivation plotter'
        self.dPFIn = self.inpD.dPFIn_ICDrv
        self.pDOut = self.inpD.pOutPDF
        self.dPlt = {**self.dPlt, **self.inpD.plot_ICDrv}
        self.getDPPltF(sFOut=self.inpD.sFOut_ICDrv)
        self.loadDDfrInp()
        print('Initiated "ICDerivPlotter" base object and loaded input data.')

    def setTicks(self):
        axXTck, self.dPlt['lblXTck'] = np.arange(len(L_S_FT_CHG)), L_S_FT_CHG
        axYLim, axYTck = self.dPlt['axYLim'], None
        if axYLim is not None and len(axYLim) >= 2:
            axYTck = range(-(-axYLim[0]//2*2), axYLim[1]//2*2 + 1, 2)
        self.dPlt['axXTck'], self.dPlt['axYTck'] = axXTck, axYTck

    def createPlot(self, cSer, lILegPlt):
        self.setTicks()
        nChD, xTck = self.dPlt['nCharDsp'], self.dPlt['axXTck']
        cFig, cAx = self.iniPlot()
        for k, (sK, lSHd) in enumerate(D_S_FT_CHG.items()):
            cSerGrp = cSer.loc[lSHd]
            cSerGrp.index = L_S_FT_CHG
            xLoc = (xTck + (2*k + 1)/(2*len(D_S_FT_CHG))*self.dPlt['wdthGrp'] -
                    1/2 + self.dPlt['wdthBar']/2)
            cAx.bar(xLoc, height=cSerGrp, width=self.dPlt['wdthBar'],
                    lw=self.dPlt['lwdPlt'], label=lILegPlt[k][:nChD])
        cAx.plot([-1/2, len(L_S_FT_CHG) - 1/2], [0, 0],
                 lw=self.dPlt['lwdPlt'], color='black')
        return cFig, cAx

    def decoratePlot(self, cAx, sGT, sYLbl=None):
        super().decoratePlot(cAx, sTtl=D_S_NM_GT[sGT], sYLbl=sYLbl)
        l = cAx.legend(loc=self.dPlt['locLegend'],
                       bbox_to_anchor=self.dPlt['coordAnchorBox'],
                       fontsize=self.dPlt['szFontLeg'])
        if self.dPlt['plotVLines']:
            yLow, yUp = self.dPlt['axYTck'][0], self.dPlt['axYTck'][-1]
            plt.vlines(np.arange(-1/2, len(L_S_FT_CHG) + 1/2), yLow, yUp,
                       lw=self.dPlt['lwdVLine'], colors=self.dPlt['clrVLine'])
        return l

    def plotICDeriv(self):
        for sID, ((sGT, sMet, sPho), pPltF) in self.dPPltF.items():
            print('Plotting IC derivation', sID, 'for "' + sGT +
                  '" and the pair (' + sMet + ', ' + sPho + ')...')
            d = self.dDfrIn[sGT]
            cSer = d[(d[S_MET] == sMet) & (d[S_PHO] == sPho)].squeeze()
            lILegPlt = [sMet, sPho, S_CMB_DEV, S_IC_CMP]
            # if not os.path.isfile(pPltF):
            cFig, cAx = self.createPlot(cSer, lILegPlt)
            cLeg = self.decoratePlot(cAx, sGT=sGT, sYLbl=S_YLBL_IC_DERIV_PLT)
            saveClosePlot(cFig, pPltF, cLeg)

# .............................................................................
class CmbToICPlotter(Plotter):
    def __init__(self, InpD, pltrICDrv):
        super().__init__(InpD, InpD.plot_CmbToIC['axYLim'])
        self.idO = self.inpD.sPltr_CmbToIC
        self.descO = 'Concordance index derivation plotter'
        self.pDOut = self.inpD.pOutPDF
        self.dPlt = {**self.dPlt, **self.inpD.plot_CmbToIC}
        self.getDPPltF(sFOut=self.inpD.sFOut_CmbToIC)
        self.dDfrIn = pltrICDrv.dDfrIn
        print('Initiated "CmbToICPlotter" base object and obtained data.')

    def setTicks(self):
        axXTck, self.dPlt['lblXTck'] = np.arange(len(L_S_FT_CHG)), L_S_FT_CHG
        axYLim, axYTck = self.dPlt['axYLim'], None
        if axYLim is not None and len(axYLim) >= 2:
            axYTck = range(-(-axYLim[0]//2*2), axYLim[1]//2*2 + 1, 2)
        self.dPlt['axXTck'], self.dPlt['axYTck'] = axXTck, axYTck

    def createPlot(self, cSer, lILegPlt):
        self.setTicks()
        nChD, xTck = self.dPlt['nCharDsp'], self.dPlt['axXTck']
        cFig, cAx = self.iniPlot()
        for k, (sK, lSHd) in enumerate(D_S_FT_CHG.items()):
            cSerGrp = cSer.loc[lSHd]
            cSerGrp.index = L_S_FT_CHG
            xLoc = (xTck + (2*k + 1)/(2*len(D_S_FT_CHG))*self.dPlt['wdthGrp'] -
                    1/2 + self.dPlt['wdthBar']/2)
            cAx.bar(xLoc, height=cSerGrp, width=self.dPlt['wdthBar'],
                    lw=self.dPlt['lwdPlt'], label=lILegPlt[k][:nChD])
        cAx.plot([-1/2, len(L_S_FT_CHG) - 1/2], [0, 0],
                 lw=self.dPlt['lwdPlt'], color='black')
        return cFig, cAx

    def decoratePlot(self, cAx, sGT, sYLbl=None):
        super().decoratePlot(cAx, sTtl=D_S_NM_GT[sGT], sYLbl=sYLbl)
        l = cAx.legend(loc=self.dPlt['locLegend'],
                       bbox_to_anchor=self.dPlt['coordAnchorBox'],
                       fontsize=self.dPlt['szFontLeg'])
        if self.dPlt['plotVLines']:
            yLow, yUp = self.dPlt['axYTck'][0], self.dPlt['axYTck'][-1]
            plt.vlines(np.arange(-1/2, len(L_S_FT_CHG) + 1/2), yLow, yUp,
                       lw=self.dPlt['lwdVLine'], colors=self.dPlt['clrVLine'])
        return l

    def plotCmbToIC(self):
        for sID, ((sGT, sMet, sPho), pPltF) in self.dPPltF.items():
            print('Plotting combining components to IC', sID, 'for "' + sGT +
                  '" and the pair (' + sMet + ', ' + sPho + ')...')
            d = self.dDfrIn[sGT]
            cSer = d[(d[S_MET] == sMet) & (d[S_PHO] == sPho)].squeeze()
            lILegPlt = [sMet, sPho, S_CMB_DEV, S_IC_CMP]
            # if not os.path.isfile(pPltF):
            cFig, cAx = self.createPlot(cSer, lILegPlt)
            cLeg = self.decoratePlot(cAx, sGT=sGT, sYLbl=S_YLBL_CMB_TO_IC_PLT)
            saveClosePlot(cFig, pPltF, cLeg)

# --- MAIN --------------------------------------------------------------------
startTime = time.time()
print('+'*25 + ' START', time.ctime(startTime), '+'*25)

inpDat = InputData(dInput)
if inpDat.doPlot_DevSD:
    cPltrDevSD = DevSDPlotter(inpDat)
    cPltrDevSD.printIDDesc()
    # cPltrDevSD.printDDfrInp(lIdxCol=[1, 2, 3, 5, 10, 20, 40])
    # cPltrDevSD.printDfrInp((S_METS, S_GT0), printCol=True)
    # cPltrDevSD.printDfrInp((S_METS, S_GT1), printCol=True)
    # cPltrDevSD.printDfrInp((S_METS, S_GT5), printCol=True)
    # cPltrDevSD.printDfrInp((S_PHOS, S_GT0), printCol=True)
    # cPltrDevSD.printDfrInp((S_PHOS, S_GT1), printCol=True)
    # cPltrDevSD.printDfrInp((S_PHOS, S_GT5), printCol=True)
    # cPltrDevSD.printAttrData()
    # cPltrDevSD.printObjInfo()
    cPltrDevSD.plotDevSD()
if inpDat.doPlot_ICDrv or inpDat.doPlot_CmbToIC:
    cPltrICDrv = ICDerivPlotter(inpDat)
if inpDat.doPlot_ICDrv:
    cPltrICDrv.printIDDesc()
    cPltrICDrv.plotICDeriv()
if inpDat.doPlot_CmbToIC:
    cPltrCmbToIC = CmbToICPlotter(inpDat, cPltrICDrv)
    cPltrCmbToIC.printIDDesc()
    cPltrCmbToIC.plotCmbToIC()

print('-'*80)
printElapsedTimeSim(startTime, time.time(), 'Total time')
print('+'*37, 'DONE.', '+'*36)

###############################################################################

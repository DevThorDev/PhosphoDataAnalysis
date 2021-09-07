# -*- coding: utf-8 -*-
###############################################################################
# --- PlotDerivIC.py ----------------------------------------------------------
###############################################################################
import os, time

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# --- CONSTANTS ---------------------------------------------------------------
S_DASH = '-'
S_USC = '_'
S_BAR = '|'
S_DOT = '.'
S_CSV = 'csv'
S_PDF = 'pdf'

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

S_10 = S_1 + S_0
S_50 = S_5 + S_0
S_51 = S_5 + S_1

S_MET = 'Metabolite'
S_PHO = 'Phosphopeptide'
L_S_M_P = [S_MET, S_PHO]

S_IC = 'IC'
S_IC_M = S_IC + '(-)'
S_IC_P = S_IC + '(+)'

S_SELECTED = 'Selected'

S_GT = 'GT'
S_GT0 = S_GT + S_0
S_GT1 = S_GT + S_1
S_GT5 = S_GT + S_5
L_S_GT = [S_GT0, S_GT1, S_GT5]

L_S_FT = ['DR', 'DS', 'NR', 'NS']
L_S_FT_CHG = [L_S_FT[i] + S_BAR + L_S_FT[j] for j in range(len(L_S_FT))
              for i in range(len(L_S_FT)) if i != j]
D_S_FT_CHG = {s: [S_USC.join([t, s]) for t in L_S_FT_CHG] for s in L_S_SPC}

S_DERIV = 'Deriv'
S_IC_DERIV = S_IC + S_DERIV

S_CMB_DEV = 'Combined deviation'
S_CMB_DEV = 'Combined deviation'

S_YLBL_IC_DERIV_PLT = 'Change (xSD) / IC derivation'

S_IC_GT0 = S_USC.join([S_IC, S_GT0])
S_IC_GT1 = S_USC.join([S_IC, S_GT1])
S_IC_GT5 = S_USC.join([S_IC, S_GT5])
L_S_IC_GT = [S_IC_GT0, S_IC_GT1, S_IC_GT5]

S_BASE_CL = 'BaseClass'
S_INP_DATA = 'InputData'
S_ROOT_CL = 'RootClass'
S_PLTR = 'Plotter'
S_IC_DERIV_PLTR = S_IC_DERIV + S_PLTR
S_PLT = 'Plot'
S_NM_IC_DERIV_PLT = S_IC_DERIV + S_PLT

R04 = 4

# --- INPUT -------------------------------------------------------------------
# --- flow control ------------------------------------------------------------
doPlotICDrv = True            # True / False
# dSpecSel = {'S': {'dropNA': True},      # key: column sel.: 'S'hort / 'F'ull
#             'F': {'dropNA': False}}     # value: data for sel., e.g. dropNA

# --- general input -----------------------------------------------------------
modDisp = 10000

# --- data specific input -----------------------------------------------------
dTupIn = {'A': (S_GT0, 'Alanine', 'ESLS(1)PGQQHVSQNTAVKPEGR'),
          'B': (S_GT0, 'Alanine', 'T(0.001)PS(0.996)QRPS(0.003)TSSSSGGFNIGK'),
          'C': (S_GT1, 'Glutamic_acid', 'TADS(1)DGES(1)GEIKFEDNNLPPLGEK'),
          'D': (S_GT5, 'Putrescine', 'S(0.999)FS(0.001)VADFPR'),
          'E': (S_GT5, 'Fructose', 'TEEDENDDEDHEHKDDKT(0.854)S(0.144)PDS(0.001)IVMVEAK')}
maxSLen = 20
sSep = ';'

# --- graphics parameters / all plots -----------------------------------------
szFontLeg = 'small'             # font size of legend
nCharDsp = 60                   # number of chars displayed for legend item
coordAnchorBox = (0.5, 1.02)    # coordinates of the legend anchor box
lWdPlt = 0.75                   # line width in plot

# --- graphics parameters / IC derivation plot --------------------------------
dTupInICDeriv = dTupIn
nmICDrvPlt = S_NM_IC_DERIV_PLT      # name prefix of the IC derivation plot
wdthBar = 0.2                       # width of single bars
wdthGrp = 0.75                      # width of bar group
degRotXLbl = 90                     # degree rotation of x-labels

# --- names and paths of files and dirs ---------------------------------------
s1, s2, s3, s4, s5, s6 = 'Corr', 'BinOp', 'MetD', 'PhoD', 'DvSD', 'AllD'
dSFInICDrv = {sGT: (S_IC_DERIV + S_USC*2 + s1 + S_USC*2 +
                    S_USC.join([s2, s3, s5, sGT, s6, s4, s5, sGT, s6]))
              for sGT in L_S_GT}
sFOutICDrv = nmICDrvPlt

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
dPFInICDrv = {sGT: os.path.join(pInCSV, dSFInICDrv[sGT] + S_DOT + S_CSV)
              for sGT in L_S_GT}

# --- assertions --------------------------------------------------------------

# --- INPUT DICTIONARY --------------------------------------------------------
dInput = {# --- constants
          'lSGT': L_S_GT,
          'sBase': S_BASE_CL,
          'sInpDat': S_INP_DATA,
          'sRoot': S_ROOT_CL,
          'sPltr': S_PLTR,
          'sICDrvPltr': S_IC_DERIV_PLTR,
          'sMet': S_MET,
          'sPho': S_PHO,
          'R04': R04,
          # --- flow control
          'doPlotICDrv': doPlotICDrv,
          # --- general input
          'modDisp': modDisp,
          # --- data specific input
          'dTupIn': dTupIn,
          'maxSLen': maxSLen,
          'sSep': sSep,
          # --- graphics parameters / IC derivation plot
          'plotOfICDeriv': {'dTupInICDeriv': dTupInICDeriv,
                            'nmICDrvPlt': nmICDrvPlt,
                            'szFontLeg': szFontLeg,
                            'nCharDsp': nCharDsp,
                            'coordAnchorBox': coordAnchorBox,
                            'lWdPlt': lWdPlt,
                            'wdthBar': wdthBar,
                            'wdthGrp': wdthGrp,
                            'degRotXLbl': degRotXLbl},
          # --- names and paths of files and dirs
          'pInCSV': pInCSV,
          'pOutPDF': pOutPDF,
          'dSFInICDrv': dSFInICDrv,
          'sFOutICDrv': sFOutICDrv,
          # --- further derived values
          'dPFInICDrv': dPFInICDrv}

# --- FUNCTIONS ---------------------------------------------------------------

def decorateClosePlot(cFig, cAx, dPlt, pPltF, sYLbl=''):
    cAx.set_ylabel(sYLbl)
    l = cAx.legend(loc='lower center', bbox_to_anchor=dPlt['coordAnchorBox'],
                   fontsize=dPlt['szFontLeg'])
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

    def printDDfrInp(self):
        if self.dDfrIn is None:
            print('Input DataFrames dictionary does not have any content yet.')
        else:
            print('Input DataFrames dictionary:')
            print(self.dDfrIn)

class Plotter(RootClass):
    def __init__(self, InpD):
        super().__init__(InpD)
        self.idO = InpD.sPltr
        self.descO = 'Class for plotting'
        self.sSp = self.inpD.sSep
        self.dPPltF = {}
        print('Initiated "Plotter" base object.')

    def printDPPltF(self):
        print('Dictionary of plot file paths:')
        for tK, pF in self.dPPltF.items():
            print(tK, ':', pF)
        print('-'*64)

    def loadDDfrInp(self):
        # load input DataFrames, and save them in dictionary
        if hasattr(self, 'dPFIn'):
            self.dDfrIn = {sGT: None for sGT in L_S_GT}
            for sGT in self.dDfrIn:
                self.dDfrIn[sGT] = pd.read_csv(self.dPFIn[sGT], sep=self.sSp)

class ICDerivPlotter(Plotter):
    def __init__(self, InpD):
        super().__init__(InpD)
        self.idO = self.inpD.sICDrvPltr
        self.descO = 'Concordance index derivation plotter'
        self.dPFIn = self.inpD.dPFInICDrv
        self.pDOut = self.inpD.pOutPDF
        self.dPlt = self.inpD.plotOfICDeriv
        self.getDPPltF()
        self.loadDDfrInp()
        print('Initiated "ICDerivPlotter" base object and loaded input data.')

    def getDPPltF(self):
        self.dPPltF, mxL = {}, self.inpD.maxSLen
        for sID, t in self.dPlt['dTupInICDeriv'].items():
            u = S_USC.join([self.inpD.sFOutICDrv, sID])
            sPltF = S_DOT.join([S_USC.join([u] + [s[:mxL] for s in t]), S_PDF])
            self.dPPltF[sID] = (t, os.path.join(self.pDOut, sPltF))

    def plotICDeriv(self):
        nChD, wdBar = self.dPlt['nCharDsp'], self.dPlt['wdthBar']
        for sID, ((sGT, sMet, sPho), pPltF) in self.dPPltF.items():
            print('Plotting IC derivation', sID, 'for "' + sGT +
                  '" and the pair (' + sMet + ', ' + sPho + ')...')
            d = self.dDfrIn[sGT]
            cSer = d[(d[S_MET] == sMet) & (d[S_PHO] == sPho)].squeeze()
            # lTDat = [(S_M, sMet), (S_P, sPho), (S_IC_DERIV, S_CMB_DEV)]
            aRg = np.arange(len(L_S_FT_CHG))
            # if not os.path.isfile(pPltF):
            cFig, cAx = plt.subplots()
            for k, (sK, lSHd) in enumerate(D_S_FT_CHG.items()):
                cSerGrp = cSer.loc[lSHd]
                cSerGrp.index = L_S_FT_CHG
                xLoc = aRg + (2*k + 1)/(2*len(D_S_FT_CHG))*self.dPlt['wdthGrp']
                cAx.bar(xLoc - 1/2, height=cSerGrp, width=wdBar,
                        lw=self.dPlt['lWdPlt'], label=tMP[1][:nChD])
            cAx.plot([-1/2, len(L_S_FT_CHG) + 1/2], [0, 0],
                     lw=self.dPlt['lWdPlt'], color='black')
            cAx.set_xticks(aRg)
            cAx.set_xticklabels(L_S_FT_CHG)
            for cXLbl in cAx.get_xticklabels():
              cXLbl.set_rotation(self.dPlt['degRotXLbl'])
            decorateClosePlot(cFig, cAx, self.dPlt, pPltF, S_YLBL_IC_DERIV_PLT)

# --- MAIN --------------------------------------------------------------------
startTime = time.time()
print('+'*25 + ' START', time.ctime(startTime), '+'*25)

inpDat = InputData(dInput)
if inpDat.doPlotICDrv:
    cPltr = ICDerivPlotter(inpDat)
    cPltr.printIDDesc()
    cPltr.printDDfrInp()
    cPltr.printAttrData()
    cPltr.printObjInfo()
    # cPltr.plotICDeriv()

print('-'*80)
printElapsedTimeSim(startTime, time.time(), 'Total time')
print('+'*37, 'DONE.', '+'*36)

###############################################################################

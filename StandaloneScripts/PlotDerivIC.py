# -*- coding: utf-8 -*-
###############################################################################
# --- PlotDerivIC.py ----------------------------------------------------------
###############################################################################
import os, time

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

S_CMB_C = 'Combined'
S_DEV_S = 'deviation'
S_DEVS_C = 'Deviations between features'
S_CMP_S = 'component'
S_CMB_DEV = S_CMB_C + S_SPACE + S_DEV_S
S_IC_CMP = S_IC_MEX + S_SPACE + S_CMP_S

S_YLBL_IC_DERIV_PLT = S_DEVS_C + S_SPACE + S_SLASH + S_SPACE + S_IC_CMP

S_BASE_CL = 'BaseClass'
S_INP_DATA = 'InputData'
S_ROOT_CL = 'RootClass'
S_PLTR = 'Plotter'
S_DEV_SD_PLTR = S_DEV_SD + S_PLTR
S_IC_DERIV_PLTR = S_IC_DERIV + S_PLTR

S_PLT = 'Plot'
S_NM_DEV_SD_PLT = S_1 + S_DEV_SD + S_PLT
S_NM_IC_DERIV_PLT = S_2 + S_USC + S_IC_DERIV + S_PLT

R04 = 4

# --- INPUT -------------------------------------------------------------------
# --- flow control ------------------------------------------------------------
doPlotDevSD = True            # True / False
doPlotICDrv = False            # True / False

# --- data specific input -----------------------------------------------------
maxSLen = 20
sSep = ';'

# --- graphics parameters / all plots -----------------------------------------
szFontLeg = 'small'             # font size of legend
nCharDsp = 60                   # number of chars displayed for legend item
coordAnchorBox = (0.5, 1.02)    # coordinates of the legend anchor box
lWdPlt = 0.75                   # line width in plot

# --- graphics parameters / deviation SD plot ---------------------------------
dTupInDevSD = {'A': (S_GT0, 'Alanine', S_FT2 + S_BAR + S_FT1),
               'B': (S_GT0, 'Alanine', S_FT2 + S_BAR + S_FT4),
               'C': (S_GT1, 'Docosanoic_acid', S_FT4 + S_BAR + S_FT1),
               'D': (S_GT1, 'Docosanoic_acid', S_FT1 + S_BAR + S_FT4),
               'E': (S_GT5, 'Putrescine', S_FT2 + S_BAR + S_FT1),
               'F': (S_GT5, 'Putrescine', S_FT3 + S_BAR + S_FT2)}
nmDevSDPlt = S_NM_DEV_SD_PLT        # name prefix of the deviation SD plot
axYLim = (-11, 11)                  # limits for the y-axis, or None
locTitle = 'left'                   # location of the title
padTitle = 70.                      # padding of the title
locLegend = 'lower center'          # location of the legend
wdthGrp = 0.8                       # width of bar group
wdthBar = wdthGrp/len(D_S_FT_CHG) - 0.01   # width of single bars
degRotXLbl = 90                     # degree rotation of x-labels
plotVLines = True                   # plot vertical lines between groups?
wdthVLine = 0.2                     # width of vertical line
clrVLine = 'k'                      # colour of vertical line

# --- graphics parameters / IC derivation plot --------------------------------
dTupInICDeriv = {'A': (S_GT0, 'Alanine', 'ESLS(1)PGQQHVSQNTAVKPEGR'),
                 'B': (S_GT0, 'Alanine', 'T(0.001)PS(0.996)QRPS(0.003)TSSSSGGFNIGK'),
                 'C': (S_GT1, 'Glutamic_acid', 'TADS(1)DGES(1)GEIKFEDNNLPPLGEK'),
                 'D': (S_GT5, 'Putrescine', 'S(0.999)FS(0.001)VADFPR'),
                 'E': (S_GT5, 'Fructose', 'TEEDENDDEDHEHKDDKT(0.854)S(0.144)PDS(0.001)IVMVEAK')}
nmICDrvPlt = S_NM_IC_DERIV_PLT      # name prefix of the IC derivation plot
axYLim = (-11, 11)                  # limits for the y-axis, or None
locTitle = 'left'                   # location of the title
padTitle = 70.                      # padding of the title
locLegend = 'lower center'          # location of the legend
wdthGrp = 0.8                       # width of bar group
wdthBar = wdthGrp/len(D_S_FT_CHG) - 0.01   # width of single bars
degRotXLbl = 90                     # degree rotation of x-labels
plotVLines = True                   # plot vertical lines between groups?
wdthVLine = 0.2                     # width of vertical line
clrVLine = 'k'                      # colour of vertical line

# --- names and paths of files and dirs ---------------------------------------
sMs, sPs, s1, s2, s3 = S_METS, S_PHOS, 'Means', 'SDs', 'Deviations'
dSFIn_Ms = {(sMs, sGT): S_USC.join([sMs, s1, s2, s3, sGT]) for sGT in L_S_GT}
dSFIn_Ps = {(sPs, sGT): S_USC.join([sPs, s1, s2, s3, sGT]) for sGT in L_S_GT}
dSFInDevSD = {**dSFIn_Ms, **dSFIn_Ps}
s1, s2, s3, s4, s5, s6 = 'Corr', 'BinOp', 'MetD', 'PhoD', 'DvSD', 'AllD'
dSFInICDrv = {sGT: (S_IC_DERIV + S_USC*2 + s1 + S_USC*2 +
                    S_USC.join([s2, s3, s5, sGT, s6, s4, s5, sGT, s6]))
              for sGT in L_S_GT}
sFOutDevSD = nmDevSDPlt
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
dPFInDevSD = {tK: os.path.join(pInCSV, dSFInDevSD[tK] + S_DOT + S_CSV)
              for tK in dSFInDevSD}
dPFInICDrv = {sGT: os.path.join(pInCSV, dSFInICDrv[sGT] + S_DOT + S_CSV)
              for sGT in L_S_GT}
axYTicks = None                     # # tick locations of the y-axis, or None
if axYLim is not None:
    assert len(axYLim) >= 2
    axYTicks = range(-(-axYLim[0]//2*2), axYLim[1]//2*2 + 1, 2)
    print('axYTicks:', list(axYTicks), '(length =', len(list(axYTicks)), ')')

# --- assertions --------------------------------------------------------------
assert len(L_S_SPC) == 4 and len(D_S_FT_CHG) == 4

# --- INPUT DICTIONARY --------------------------------------------------------
dInput = {# --- constants
          'lSGT': L_S_GT,
          'sBase': S_BASE_CL,
          'sInpDat': S_INP_DATA,
          'sRoot': S_ROOT_CL,
          'sPltr': S_PLTR,
          'sDevSDPltr': S_DEV_SD_PLTR,
          'sICDrvPltr': S_IC_DERIV_PLTR,
          'sMet': S_MET,
          'sPho': S_PHO,
          'R04': R04,
          # --- flow control
          'doPlotDevSD': doPlotDevSD,
          'doPlotICDrv': doPlotICDrv,
          # --- data specific input
          'maxSLen': maxSLen,
          'sSep': sSep,
          # --- graphics parameters / IC derivation plot
          'plotOfDevSD': {'dTupInDevSD': dTupInDevSD,
                            'nmDevSDPlt': nmDevSDPlt,
                            'axYTicks': axYTicks,
                            'axYLim': axYLim,
                            'locTitle': locTitle,
                            'padTitle': padTitle,
                            'locLegend': locLegend,
                            'szFontLeg': szFontLeg,
                            'nCharDsp': nCharDsp,
                            'coordAnchorBox': coordAnchorBox,
                            'lWdPlt': lWdPlt,
                            'wdthGrp': wdthGrp,
                            'wdthBar': wdthBar,
                            'degRotXLbl': degRotXLbl,
                            'plotVLines': plotVLines,
                            'wdthVLine': wdthVLine,
                            'clrVLine': clrVLine},
          # --- graphics parameters / IC derivation plot
          'plotOfICDeriv': {'dTupInICDeriv': dTupInICDeriv,
                            'nmICDrvPlt': nmICDrvPlt,
                            'axYTicks': axYTicks,
                            'axYLim': axYLim,
                            'locTitle': locTitle,
                            'padTitle': padTitle,
                            'locLegend': locLegend,
                            'szFontLeg': szFontLeg,
                            'nCharDsp': nCharDsp,
                            'coordAnchorBox': coordAnchorBox,
                            'lWdPlt': lWdPlt,
                            'wdthGrp': wdthGrp,
                            'wdthBar': wdthBar,
                            'degRotXLbl': degRotXLbl,
                            'plotVLines': plotVLines,
                            'wdthVLine': wdthVLine,
                            'clrVLine': clrVLine},
          # --- names and paths of files and dirs
          'pInCSV': pInCSV,
          'pOutPDF': pOutPDF,
          'dSFInDevSD': dSFInDevSD,
          'dSFInICDrv': dSFInICDrv,
          'sFOutDevSD': sFOutDevSD,
          'sFOutICDrv': sFOutICDrv,
          # --- further derived values
          'dPFInDevSD': dPFInDevSD,
          'dPFInICDrv': dPFInICDrv}

# --- FUNCTIONS ---------------------------------------------------------------
def saveClosePlot(cFig, l, pPltF):
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
            self.dDfrIn = {}
            for sK in self.dPFIn:
                self.dDfrIn[sK] = pd.read_csv(self.dPFIn[sK], sep=self.sSp)

class DevSDPlotter(Plotter):
    def __init__(self, InpD):
        super().__init__(InpD)
        self.idO = self.inpD.sDevSDPltr
        self.descO = 'Deviations (as multiples of feature SD) plotter'
        self.dPFIn = self.inpD.dPFInDevSD
        self.pDOut = self.inpD.pOutPDF
        self.dPlt = self.inpD.plotOfDevSD
        self.getDPPltF()
        self.loadDDfrInp()
        print('Initiated "DevSDPlotter" base object and loaded input data.')

    def getDPPltF(self):
        self.dPPltF, mxL = {}, self.inpD.maxSLen
        for sID, t in self.dPlt['dTupInDevSD'].items():
            u = S_USC.join([self.inpD.sFOutDevSD, sID])
            sPltF = S_DOT.join([S_USC.join([u] + [s[:mxL] for s in t]), S_PDF])
            self.dPPltF[sID] = (t, os.path.join(self.pDOut, sPltF))

    def createPlot(self, cSer, lILegPlt):
        nChD, aRg = self.dPlt['nCharDsp'], np.arange(len(L_S_FT_CHG))
        cFig, cAx = plt.subplots()
        for k, (sK, lSHd) in enumerate(D_S_FT_CHG.items()):
            cSerGrp = cSer.loc[lSHd]
            cSerGrp.index = L_S_FT_CHG
            xLoc = (aRg + (2*k + 1)/(2*len(D_S_FT_CHG))*self.dPlt['wdthGrp'] -
                    1/2 + self.dPlt['wdthBar']/2)
            cAx.bar(xLoc, height=cSerGrp, width=self.dPlt['wdthBar'],
                    lw=self.dPlt['lWdPlt'], label=lILegPlt[k][:nChD])
        cAx.plot([-1/2, len(L_S_FT_CHG) - 1/2], [0, 0],
                 lw=self.dPlt['lWdPlt'], color='black')
        return cFig, cAx, aRg

    def decoratePlot(self, cAx, sGT, aRg, sYLbl=''):
        if self.dPlt['axYLim'] is not None:
            cAx.set_ylim(self.dPlt['axYLim'])
        cAx.set_xticks(aRg)
        cAx.set_xticklabels(L_S_FT_CHG)
        if self.dPlt['axYTicks'] is not None:
            cAx.set_yticks(self.dPlt['axYTicks'])
        for cXLbl in cAx.get_xticklabels():
            cXLbl.set_rotation(self.dPlt['degRotXLbl'])
        cAx.set_ylabel(sYLbl)
        cAx.set_title(D_S_NM_GT[sGT], loc=self.dPlt['locTitle'],
                      pad=self.dPlt['padTitle'])
        l = cAx.legend(loc=self.dPlt['locLegend'],
                       bbox_to_anchor=self.dPlt['coordAnchorBox'],
                       fontsize=self.dPlt['szFontLeg'])
        if self.dPlt['plotVLines']:
            yLow, yUp = self.dPlt['axYTicks'][0], self.dPlt['axYTicks'][-1]
            plt.vlines(np.arange(-1/2, len(L_S_FT_CHG) + 1/2), yLow, yUp,
                       lw=self.dPlt['wdthVLine'], colors=self.dPlt['clrVLine'])
        return l

    def plotDevSD(self):
        for sID, ((sGT, sMet, sPho), pPltF) in self.dPPltF.items():
            print('Plotting IC derivation', sID, 'for "' + sGT +
                  '" and the pair (' + sMet + ', ' + sPho + ')...')
            d = self.dDfrIn[sGT]
            cSer = d[(d[S_MET] == sMet) & (d[S_PHO] == sPho)].squeeze()
            lILegPlt = [sMet, sPho, S_CMB_DEV, S_IC_CMP]
            # if not os.path.isfile(pPltF):
            cFig, cAx, aRg = self.createPlot(cSer, lILegPlt)
            cLeg = self.decoratePlot(cAx, sGT, aRg, sYLbl=S_YLBL_IC_DERIV_PLT)
            saveClosePlot(cFig, cLeg, pPltF)

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

    def createPlot(self, cSer, lILegPlt):
        nChD, aRg = self.dPlt['nCharDsp'], np.arange(len(L_S_FT_CHG))
        cFig, cAx = plt.subplots()
        for k, (sK, lSHd) in enumerate(D_S_FT_CHG.items()):
            cSerGrp = cSer.loc[lSHd]
            cSerGrp.index = L_S_FT_CHG
            xLoc = (aRg + (2*k + 1)/(2*len(D_S_FT_CHG))*self.dPlt['wdthGrp'] -
                    1/2 + self.dPlt['wdthBar']/2)
            cAx.bar(xLoc, height=cSerGrp, width=self.dPlt['wdthBar'],
                    lw=self.dPlt['lWdPlt'], label=lILegPlt[k][:nChD])
        cAx.plot([-1/2, len(L_S_FT_CHG) - 1/2], [0, 0],
                 lw=self.dPlt['lWdPlt'], color='black')
        return cFig, cAx, aRg

    def decoratePlot(self, cAx, sGT, aRg, sYLbl=''):
        if self.dPlt['axYLim'] is not None:
            cAx.set_ylim(self.dPlt['axYLim'])
        cAx.set_xticks(aRg)
        cAx.set_xticklabels(L_S_FT_CHG)
        if self.dPlt['axYTicks'] is not None:
            cAx.set_yticks(self.dPlt['axYTicks'])
        for cXLbl in cAx.get_xticklabels():
            cXLbl.set_rotation(self.dPlt['degRotXLbl'])
        cAx.set_ylabel(sYLbl)
        cAx.set_title(D_S_NM_GT[sGT], loc=self.dPlt['locTitle'],
                      pad=self.dPlt['padTitle'])
        l = cAx.legend(loc=self.dPlt['locLegend'],
                       bbox_to_anchor=self.dPlt['coordAnchorBox'],
                       fontsize=self.dPlt['szFontLeg'])
        if self.dPlt['plotVLines']:
            yLow, yUp = self.dPlt['axYTicks'][0], self.dPlt['axYTicks'][-1]
            plt.vlines(np.arange(-1/2, len(L_S_FT_CHG) + 1/2), yLow, yUp,
                       lw=self.dPlt['wdthVLine'], colors=self.dPlt['clrVLine'])
        return l

    def plotICDeriv(self):
        for sID, ((sGT, sMet, sPho), pPltF) in self.dPPltF.items():
            print('Plotting IC derivation', sID, 'for "' + sGT +
                  '" and the pair (' + sMet + ', ' + sPho + ')...')
            d = self.dDfrIn[sGT]
            cSer = d[(d[S_MET] == sMet) & (d[S_PHO] == sPho)].squeeze()
            lILegPlt = [sMet, sPho, S_CMB_DEV, S_IC_CMP]
            # if not os.path.isfile(pPltF):
            cFig, cAx, aRg = self.createPlot(cSer, lILegPlt)
            cLeg = self.decoratePlot(cAx, sGT, aRg, sYLbl=S_YLBL_IC_DERIV_PLT)
            saveClosePlot(cFig, cLeg, pPltF)

# --- MAIN --------------------------------------------------------------------
startTime = time.time()
print('+'*25 + ' START', time.ctime(startTime), '+'*25)

inpDat = InputData(dInput)
if inpDat.doPlotDevSD:
    cPltr = DevSDPlotter(inpDat)
    cPltr.printIDDesc()
    cPltr.printDDfrInp(lIdxCol=[1, 2, 3, 5, 10, 20, 40])
    # cPltr.printAttrData()
    # cPltr.printObjInfo()
    # cPltr.plotDevSD()
if inpDat.doPlotICDrv:
    cPltr = ICDerivPlotter(inpDat)
    cPltr.plotICDeriv()

print('-'*80)
printElapsedTimeSim(startTime, time.time(), 'Total time')
print('+'*37, 'DONE.', '+'*36)

###############################################################################

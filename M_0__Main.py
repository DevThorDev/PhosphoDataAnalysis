# -*- coding: utf-8 -*-
###############################################################################
# --- M_0__Main.py ------------------------------------------------------------
###############################################################################
import os

import Core.F_00__GenFunctions as GF
import Core.F_04__MainFunctions as MF

from Control.A_00__GenInput import dictInpG
from Core.C_00__GenConstants import NM_OBJINP
from Core.I_01__InpData import InputData

# ### MAIN ####################################################################
startTime = GF.startSimu()
# -----------------------------------------------------------------------------
print('='*80, '\n', '-'*33, 'M_0__Main.py', '-'*33, '\n')
print('Current working directory:', os.getcwd())
inDG = InputData(dictInpG)
inDG.addObjTps(NM_OBJINP)
print('Added object types.')
# -----------------------------------------------------------------------------
# loop through genotypes | WT  | PGM | SIRK1 | SIRK1-PGM | SIRK1-SWEET | SWEET
#                        | GT0 | GT1 | GT2   | GT3       | GT4         | GT5
t = MF.initAnalysis(inDG)
lSUsDC, lSMgDC, lSCRDC, lSUsDG, lSMgDG, lSCRDG, lMD, lPD, lCD, dClR, lIDt = t
for k, cGT in enumerate(inDG.dI['lInvGT']):
    print('-'*80, '\n', '*'*20, 'Current genotype:', '|', cGT, '|',
          inDG.dI['dIDGT'][cGT], '|', '*'*20)
    
    # load MetData, and calculate means and SDs
    if inDG.dI['doXFtC'] or inDG.dI['doXFtG']:
        cMD = MF.loadData(inDG, inDG.dI['iMetD'], cGT, lMD, sD = lSUsDC[0])
        print('='*20, 'Loaded MetDataC, calculated means and SDs.', '='*20)
    
    # load PhoData, and calculate means and SDs
    if inDG.dI['doXFtC'] or inDG.dI['doXFtG']:
        cPD = MF.loadData(inDG, inDG.dI['iPhoD'], cGT, lPD, sD = lSUsDC[1])
        print('='*20, 'Loaded PhoDataC, calculated means and SDs.', '='*20)
    
    # do a full run of the model for each genotype (XFeature type "C")
    MF.calcCurXFt(inDG, k, inDG.dI['nmXFtC'], inDG.dI['doXFtC'], t[:6], cMD,
                  cPD, lMD, lPD, lCD, dClR, startTime)

# -----------------------------------------------------------------------------
# calculate clusters and correlations thereof of the all-genotype data
MF.calcClustCorrClust(inDG, inDG.dI['nmXFtC'], inDG.dI['doXFtC'], t[:6], lMD,
                      lPD, lCD, startTime)

# -----------------------------------------------------------------------------
# extract the feature type "G" data from the extended data
dMDG, dPDG = MF.extractDataXFt(inDG, inDG.dI['doXFtG'], lSUsDC, lMD, lPD,
                               startTime)

# -----------------------------------------------------------------------------
# loop through original features | DR  | DS  | NR  | NS
#                                | FT1 | FT2 | FT3 | FT4
lMD, lPD, lCD, dClR = [], [], [], {}
# get feature type "G" data from the dictionaries
if inDG.dI['doXFtG']:
    MF.getXFtData(dMDG, dPDG, lMD, lPD)
for k, cMD in enumerate(lMD):
    cPD = lPD[k]
    # do a full run of the model for each data object (XFeature type "G")
    MF.calcCurXFt(inDG, k, inDG.dI['nmXFtG'], inDG.dI['doXFtG'], t[:6], cMD,
                  cPD, lMD, lPD, lCD, dClR, startTime)

# -----------------------------------------------------------------------------
# calculate clusters and correlations thereof of the all-genotype data
MF.calcClustCorrClust(inDG, inDG.dI['nmXFtG'], inDG.dI['doXFtG'], t[:6], lMD,
                      lPD, lCD, startTime)

# # -----------------------------------------------------------------------------
# # extract the feature type "C" data from the extended data
# doXFtCG = inDG.dI['doXFtC'] and inDG.dI['doXFtG']
# dMDC, dPDC = MF.extractDataXFt(inDG, doXFtCG, lSUsDG, lMD, lPD, startTime)

# # TEST for consistency ONLY ---------------------------------------------------
# # loop through genotypes | WT  | PGM | SIRK1 | SIRK1-PGM | SIRK1-SWEET | SWEET
# #                        | GT0 | GT1 | GT2   | GT3       | GT4         | GT5
# lMD, lPD, lCD, dClR = [], [], [], {}
# # get feature type "C" data from the dictionaries
# if inDG.dI['doXFtC'] and inDG.dI['doXFtG']:
#     MF.getXFtData(dMDC, dPDC, lMD, lPD)
# for k, cMDC in enumerate(lMD):
#     cPDC = lPD[k]
#     # combine MetData and PhoData
#     if inDG.dI['doXFtC'] and inDG.dI['doXFtG'] and inDG.dI['doCombDt']:
#         cCDC = MF.combineBasicDataC(inDG, [cMDC, cPDC], lCD, lSUsDC)
#         print('='*20, 'Combined MetDataC and PhoDataC to CombDataC.', '='*20)
    
#     # calculate the correlations between MetData and PhoData
#     if inDG.dI['calcCorrMetPho'] and inDG.dI['doXFtC'] and inDG.dI['doXFtG']:
#         MF.calcCorrelations(inDG, [cMDC, cPDC], lSUsDC)
#         GF.showElapsedTime(startTime)
#         print('='*20, 'Calculated MetDataC-PhoDataC correlations.', '='*20)
# -----------------------------------------------------------------------------
GF.endSimu(startTime)

# -----------------------------------------------------------------------------
#sFiltSTD = ''
#sFiltA = '2_4_30_34'
# -----------------------------------------------------------------------------
#print('-'*80)
#print(inDG)
#print('-'*80, '\n', '*'*20, 'MetData11:', '*'*20)
#MetData21 = MetData(inDG, 21, sFiltSTD, '', (-1, 54))
#MetData11 = MetData(inDG, 11, sFiltSTD)
#print(MetData11)
#MetData11.printObjInfo()
#MetData11.printAddIDfr()
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
#print('-'*80, '\n', '*'*20, 'PhoData31_All (NEW):', '*'*20)
#PhoData31 = PhoData(inDG, 31, sFiltSTD)
#sQry = ''
#sQry = '~(bin_new == 38021 | bin_new == 38721) & (bin1_new == 4) & (bin_short_new == "4.1")'
#PhoData31 = PhoData(inDG, 31, sFiltSTD, sQry)
#PhoData31 = PhoData(inDG, 31, sFiltSTD, sQry, (-1, 35))
#print(PhoData31)
#PhoData31.printObjInfo()
#PhoData31.printDfr()
#PhoData31.printAddIDfr()
###############################################################################

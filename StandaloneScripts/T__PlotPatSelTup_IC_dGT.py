# -*- coding: utf-8 -*-
###############################################################################
# --- PlotPatSelTup_IC_dGT.py -------------------------------------------------
# Select data according to specified limits for the concordance index and the
# distance between genotypes, and plot the resulting patterns
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
S_D = 'd'
S_M = 'm'
S_P = 'p'
S_N = 'N'
S_Y = 'Y'
S_NO = 'No'
L_N = [S_N]
L_Y = [S_Y]
L_NY = [S_N, S_Y]

S_0 = '0'
S_1 = '1'
S_2 = '2'
S_3 = '3'
S_4 = '4'
S_5 = '5'

S_10 = S_1 + S_0
S_50 = S_5 + S_0
S_51 = S_5 + S_1

S_SRT = 'Sort'
S_THR = 'Threshold'
S_SEL = 'Selection'
S_ALL = 'All'
S_SGL = 'Sgl'

S_MET = 'Metabolite'
S_PHO = 'Phosphopeptide'
L_S_M_P = [S_MET, S_PHO]

S_BIN_L = 'BinCode'
S_BIN_L_2 = S_BIN_L + S_2
S_SELECTED = 'Selected'

S_GT = 'GT'
S_GT0 = S_GT + S_0
S_GT1 = S_GT + S_1
S_GT5 = S_GT + S_5
L_S_GT = [S_GT0, S_GT1, S_GT5]
T_S_GT = tuple(L_S_GT)
T_S_GT0 = tuple([S_GT0])
T_S_GT1 = tuple([S_GT1])
T_S_GT5 = tuple([S_GT5])

S_D_GT_M = S_D + S_GT + 'M'
S_D_GT_P = S_D + S_GT + 'P'
L_S_D_GT = [S_D_GT_M, S_D_GT_P]

L_S_FT = ['DR', 'DS', 'NR', 'NS']
L_S_FT_CHG = [L_S_FT[i] + S_BAR + L_S_FT[j] for j in range(len(L_S_FT))
              for i in range(len(L_S_FT)) if i != j]

D_HD_C_PA = {sMP: {sGT: [S_USC.join([sGT, sFt, sMP[0]]) for sFt in L_S_FT]
                   for sGT in L_S_GT} for sMP in L_S_M_P}

S_MIN = 'min'
S_MAX = 'max'
L_S_MIN_MAX = [S_MIN, S_MAX]

S_IC = 'IC'
S_IC_P = S_IC + S_USC + 'P'
S_IC_N = S_IC + S_USC + 'N'
S_PEAR_CR = 'PearsonCorr'
S_SPEAR_CR = 'SpearmanCorr'
S_PEAR_P = 'Pearson_pVal'
S_SPEAR_P = 'Spearman_pVal'
S_IC_P_S = S_IC + 'P'
S_IC_N_S = S_IC + 'N'
S_PEAR_CR_S = 'PeC'
S_SPEAR_CR_S = 'SpC'
S_PEAR_P_S = 'Pep'
S_SPEAR_P_S = 'Spp'

S_D_CL = 'd' + S_USC + 'Cl'
S_D_E = 'd' + S_USC + 'E'
S_D_GT = S_D + S_USC + S_GT
S_D_CL_S = 'dCl'
S_D_E_S = 'dE'
S_D_GT_S = S_D + S_GT
L_S_D_2GT = [S_D + S_USC*2 + s1 + S_USC + s2 for s1 in L_S_GT for s2 in L_S_GT
             if int(s1[-1]) > int(s2[-1])]
S_ICC = S_IC + 'C'
S_ICC_LEG = 'Concordance index component'
D_HD_C_ICC = {sMPI: [S_USC.join([sFtChg, sMPI[0]]) for sFtChg in L_S_FT_CHG]
              for sMPI in L_S_M_P + [S_ICC]}

S_SRT_BY = S_SRT + 'edBy'
S_ORD = 'Order'
S_ASC = 'Asc'
S_DSC = 'Dsc'
S_YLBL_PAT_PLT = 'Pattern (z-score)'
S_YLBL_ICC_PLT = 'Change (xSD) / IC component'

S_IC_GT0 = S_USC.join([S_IC, S_GT0])
S_IC_GT1 = S_USC.join([S_IC, S_GT1])
S_IC_GT5 = S_USC.join([S_IC, S_GT5])
L_S_IC_GT = [S_IC_GT0, S_IC_GT1, S_IC_GT5]

S_SG = 'S5'
S_SG_M = 'M' + S_SG
S_SG_P = 'P' + S_SG
S_SG_MP = 'MP' + S_SG
S_SG_MET = 'MetSig5'
S_SG_PHO = 'PhoSig5'
L_S_SG_S = [S_SG_M, S_SG_P]
L_S_SG = [S_SG_MET, S_SG_PHO]
L_S_SG_GT = [S_USC.join([sSg, sGT]) for sGT in L_S_GT for sSg in L_S_SG]

S_SB = 'SB'
S_SB_P = 'P' + S_SB
L_S_SB_GT = [S_USC.join([sSB, sGT]) for sGT in L_S_GT for sSB in L_S_SG]
L_S_SB_S = [S_SB_P]

S_BASE_CL = 'BaseClass'
S_INP_DATA = 'InputData'
S_ROOT_CL = 'RootClass'
S_F_NM_CMP = 'FNmCmp'
S_F_NM_CMP_DGT = S_F_NM_CMP + 'DGT'
S_F_NM_CMP_IC_ALL_GT = S_F_NM_CMP + S_IC + S_ALL + S_GT
S_F_NM_CMP_IC_SGL_GT = S_F_NM_CMP + S_IC + S_SGL + S_GT
S_F_NM_CNSTR = 'FileNameConstructor'
S_EXTR_INFO = 'ExtrInfo'
S_PLTR = 'Plotter'
S_PAT_PLTR = 'PatternPlotter'
S_ICC_PLTR = 'ICCmpPlotter'

S_RMNG_COL1_IC = 'PearsonCorr'
S_NEW_IDX = 'NewIndex'
L_S_NO_GT = L_S_M_P + []
L_S_ADD_GT = [S_RMNG_COL1_IC, 'SpearmanCorr', 'Pearson_pVal', 'Spearman_pVal',
              'IC_N', 'IC_P', 'IC', 'MetSig5', 'PhoSig5']

S_NM_PAT_PLT = 'PatternPlot'
S_NM_ICC_PLT = 'ICCmpPlot'

R04 = 4

# --- INPUT -------------------------------------------------------------------
# --- flow control ------------------------------------------------------------
doInfoExtr = True               # True / False
doPlotPat = True                # True / False
doPlotICC = True                # True / False
dSpecSel = {'S': {'dropNA': True},      # key: column sel.: 'S'hort / 'F'ull
            'F': {'dropNA': False}}     # value: data for sel., e.g. dropNA

# --- general input -----------------------------------------------------------
modDisp = 10000

# --- data specific input -----------------------------------------------------
sortInFNm = False               # sorting in file name? (True / False)

dUsedK = {S_IC: S_IC,           # key (col. hdr.) for the IC file
          S_D_GT_M: S_D_GT,     # key (col. hdr.) for the dGTM file
          S_D_GT_P: S_D_GT}     # key (col. hdr.) for the dGTP file
# dUsedK = {S_IC: S_IC,           # key (col. hdr.) for the IC file
#           S_D_GT_M: S_50,     # key (col. hdr.) for the dGTM file
#           S_D_GT_P: S_50}     # key (col. hdr.) for the dGTP file

dISort = {S_IC: {S_GT0: {S_SRT_BY: dUsedK[S_IC], S_ORD: S_DSC},
                 S_GT1: {S_SRT_BY: dUsedK[S_IC], S_ORD: S_DSC},
                 S_GT5: {S_SRT_BY: dUsedK[S_IC], S_ORD: S_DSC}},
          S_D_GT_M: {S_SRT_BY: dUsedK[S_D_GT_M], S_ORD: S_DSC},
          S_D_GT_P: {S_SRT_BY: dUsedK[S_D_GT_P], S_ORD: S_DSC}}

dThr = {S_IC: {S_GT0: {S_MIN: None, S_MAX: None},
               S_GT1: {S_MIN: None, S_MAX: None},
               S_GT5: {S_MIN: None, S_MAX: None}},
        S_D_GT_M: {S_MIN: None, S_MAX: None},
        S_D_GT_P: {S_MIN: None, S_MAX: None}}
# dThr = {S_IC: {S_GT0: {S_MIN: 7.25, S_MAX: None},
#                S_GT1: {S_MIN: 7.25, S_MAX: None},
#                S_GT5: {S_MIN: 7.25, S_MAX: None}},
#         S_D_GT_M: {S_MIN: None, S_MAX: None},
#         S_D_GT_P: {S_MIN: None, S_MAX: None}}

lSelSGM_GT0, lSelSGM_GT1, lSelSGM_GT5 = L_NY, L_NY, L_NY
lSelSGP_GT0, lSelSGP_GT1, lSelSGP_GT5 = L_NY, L_NY, L_NY
lSelSB = L_NY                    # L_NY / L_Y (sel. bins only) / L_N
dSel = {(S_SG, S_SG_MP): {S_SG_M: {S_GT0: lSelSGM_GT0, S_GT1: lSelSGM_GT1,
                                   S_GT5: lSelSGM_GT5},
                          S_SG_P: {S_GT0: lSelSGP_GT0, S_GT1: lSelSGP_GT1,
                                   S_GT5: lSelSGP_GT5}},
        (S_SB, S_SB_P): {S_SB_P: {S_GT0: lSelSB, S_GT1: lSelSB,
                                  S_GT5: lSelSB}}}

sSep = ';'

# --- graphics parameters / all plots -----------------------------------------
szFontLeg = 'small'             # font size of legend
nCharDsp = 60                   # number of chars displayed for legend item
posLegXY = (0.5, 1.02)          # coordinates of the legend anchor box
lWdPlt = 1.                     # line width in plot
clrDef = 'k'                    # default colour

# --- graphics parameters / pattern plot --------------------------------------
# dPairsPaP = {(('Leu_STTTTV'), T_S_GT):
#              ('Leucine', 'STTTTVS(0.003)S(0.996)VHS(0.001)PTTDQDFSK')}
# dPairsPaP = {(('Tetra_S(0.001)AS(0.749)T(0.251)P'), T_S_GT):
#               ('Tetradecanoic_acid', 'S(0.001)AS(0.749)T(0.251)PLLNSLVHVS(0.179)S(0.821)PRDS(1)PIETVESVHQIQR'),
#               (('Beta_ADKTDII'), T_S_GT):
#               ('Beta-alanine', 'ADKTDIIS(0.607)S(0.117)S(0.12)S(0.156)DKAS(1)PPPPSAFR'),
#               (('Hexa_S(0.001)AS(0.749)T(0.251)P'), T_S_GT):
#               ('Hexadecanoic_acid', 'S(0.001)AS(0.749)T(0.251)PLLNSLVHVS(0.179)S(0.821)PRDS(1)PIETVESVHQIQR')}
# dPairsPaP = {(('Isoleu_DLDVNE'), T_S_GT):
#               ('Isoleucine', 'DLDVNES(1)GPPAAR'),
#               (('Val_DLDVNE'), T_S_GT):
#               ('Valine', 'DLDVNES(1)GPPAAR'),
#               (('Aspart_SDKPLNY'), T_S_GT):
#               ('Aspartic_acid', 'SDKPLNYS(1)PDPENESGINER'),
#               (('Aspart_DLDVNE'), T_S_GT):
#               ('Aspartic_acid', 'DLDVNES(1)GPPAAR'),
#               (('Malic_DLDVNE'), T_S_GT):
#               ('Malic_acid', 'DLDVNES(1)GPPAAR')}

# IC_No_No_SBNY_dGT_No_No_GT0_MPS5Y_GT1_MPS5NY_GT5_MPS5NY
dPairsPaP = {# GT0 / WT
             (('Nona_T(0.011)FDELS(0.759)D'), T_S_GT0): ('Nonanoic_acid', 'T(0.011)FDELS(0.759)DT(0.23)EVYEDS(1)D'),
             (('Citric_DNKEVTF'), T_S_GT0): ('Citric_acids', 'DNKEVTFGDLGS(1)KR'),
             (('myo-Ino_TFDELS(1)D'), T_S_GT0): ('myo-Inositol', 'TFDELS(1)DTEVYEDS(1)D'),
             (('Lys_DPS(1)PPPLS'), T_S_GT0): ('Lysine', 'DPS(1)PPPLSSLGK'),
             (('Sucr_IQEGPEGS(1)L'), T_S_GT0): ('Sucrose', 'IQEGPEGS(1)LQS(1)EMK')}

# IC_No_No_SBNY_dGT_No_No_GT0_MPS5NY_GT1_MPS5Y_GT5_MPS5NY
dPairsPaP = {# GT1 / PGM
             (('Aspartic_QGTLPTVIE'), T_S_GT1): ('Aspartic_acid', 'QGTLPTVIEEDDS(0.016)S(0.977)ET(0.007)'),
             (('Phenylala_S(0.825)RS(0.137)VDE'), T_S_GT1): ('Phenylalanine', 'S(0.825)RS(0.137)VDES(0.039)FANSFSPR'),
             (('Phenylala_S(0.001)GRT(0.004)S(0.996)E'), T_S_GT1): ('Phenylalanine', 'S(0.001)GRT(0.004)S(0.996)EPNS(1)EDEAAGVGK'),
             (('Ornith_SYS(1)GSLYR'), T_S_GT1): ('Ornithine', 'SYS(1)GSLYR')}

# IC_No_No_SBNY_dGT_No_No_GT0_MPS5NY_GT1_MPS5NY_GT5_MPS5Y
dPairsPaP = {# GT5 / SWEET
             (('Aspartic_IGS(0.999)S(0.001)EML'), T_S_GT5): ('Aspartic_acid', 'IGS(0.999)S(0.001)EMLIEGEDVR'),
             (('Aspartic_VTLVPPS(0.407)D'), T_S_GT5): ('Aspartic_acid', 'VTLVPPS(0.407)DS(0.593)PELS(0.999)PINT(0.001)PK'),
             (('Beta-ala_DIS(1)PTAAG'), T_S_GT5): ('Beta-alanine', 'DIS(1)PTAAGLGLPVTGGK'),
             (('Beta-ala_IGS(0.999)S(0.001)EML'), T_S_GT5): ('Beta-alanine', 'IGS(0.999)S(0.001)EMLIEGEDVR'),
             (('Docosan_LSRPGS(1)G'), T_S_GT5): ('Docosanoic_acid', 'LSRPGS(1)GS(1)VSGLASQR'),
             (('Lys_HPQWQSDDG'), T_S_GT5): ('Lysine', 'HPQWQSDDGGDNS(1)EPESPSDSLR'),
             (('Lys_VSS(1)FEAL'), T_S_GT5): ('Lysine', 'VSS(1)FEALQPATR'),
             (('Lys_ETLNRPAAP'), T_S_GT5): ('Lysine', 'ETLNRPAAPTNYVAISKEEAASSPVSGAADHQVPAS(1)P'),
             (('Ornith_TDSEVTSLA'), T_S_GT5): ('Ornithine', 'TDSEVTSLAAS(0.024)S(0.976)PARS(1)PR'),
             (('Phosphoric_VTLVPPS(0.006)D'), T_S_GT5): ('Phosphoric_acid', 'VTLVPPS(0.407)DS(0.593)PELS(0.999)PINT(0.001)PK')}

# IC_No_No_dGT_No_No_MPS5NY_PSBNY
dPairsPaP = {# ConcStrongPos
             (('Docosan_SLEELS(1)GEA'), T_S_GT0): ('Docosanoic_acid', 'SLEELS(1)GEAEVS(1)HDEK'),
             (('Val_SLEELS(1)GEA'), T_S_GT0): ('Valine', 'SLEELS(1)GEAEVS(1)HDEK'),
             (('Val_SDKPLNYS(1)P'), T_S_GT0): ('Valine', 'SDKPLNYS(1)PDPENESGINER'),
             (('Putres_T(1)AILERR'), T_S_GT0): ('Putrescine', 'T(1)AILERR'),
             (('Putres_TKDELT(1)EE'), T_S_GT0): ('Putrescine', 'TKDELT(1)EEES(1)LSGKDYLDPPPVK'),
             (('Nicotin_SDKPLNYS(1)P'), T_S_GT0): ('Nicotinic_acid', 'SDKPLNYS(1)PDPENESGINER'),
             (('Docosan_SDKPLNYS(1)P'), T_S_GT0): ('Docosanoic_acid', 'SDKPLNYS(1)PDPENESGINER'),
             (('Val_ALGSFGS(1)F'), T_S_GT0): ('Valine', 'ALGSFGS(1)FGS(0.999)FRS(0.001)FA'),
             (('Fumar_TKDELT(1)EE'), T_S_GT1): ('Fumaric_acid', 'TKDELT(1)EEES(1)LSGKDYLDPPPVK'),
             (('Docosan_ALGSFGS(1)F'), T_S_GT1): ('Docosanoic_acid', 'ALGSFGS(1)FGS(0.999)FRS(0.001)FA'),
             (('Nicotin_ALGSFGS(1)F'), T_S_GT1): ('Nicotinic_acid', 'ALGSFGS(1)FGS(0.999)FRS(0.001)FA'),
             (('Proline_ALGSFGS(1)F'), T_S_GT1): ('Proline', 'ALGSFGS(1)FGS(0.999)FRS(0.001)FA'),
             (('Docosan_SLEELS(1)GEA'), T_S_GT1): ('Docosanoic_acid', 'SLEELS(1)GEAEVS(1)HDEK'),
             (('Hexadecan_TKDELT(1)EE'), T_S_GT5): ('Hexadecanoic_acid', 'TKDELT(1)EEES(1)LSGKDYLDPPPVK'),
             (('Hexadecan_AYGS(1)VRS(1)Q'), T_S_GT5): ('Hexadecanoic_acid', 'AYGS(1)VRS(1)QLHELHA'),
             (('Hexadecan_T(1)AILERR'), T_S_GT5): ('Hexadecanoic_acid', 'T(1)AILERR'),
             (('Aspartic_SLEELS(1)GEA'), T_S_GT5): ('Aspartic_acid', 'SLEELS(1)GEAEVS(1)HDEK'),
             (('Docosan_YVS(1)PEGS(1)P'), T_S_GT5): ('Docosanoic_acid', 'YVS(1)PEGS(1)PFKIENPK'),
             (('Docosan_S(0.007)ES(0.993)LGH'), T_S_GT5): ('Docosanoic_acid', 'S(0.007)ES(0.993)LGHRS(1)DVS(1)S(1)PEAK'),
             (('Octadecan_AYGS(1)VRS(1)Q'), T_S_GT5): ('Octadecanoic_acid', 'AYGS(1)VRS(1)QLHELHA'),
             (('Putres_T(1)AILERR'), T_S_GT5): ('Putrescine', 'T(1)AILERR'),
             # ConcZero
             (('Ala_TFDELS(1)DG'), T_S_GT1): ('Alanine', 'TFDELS(1)DGEVYEDS(1)D'),
             # ConcStrongNeg
             (('Docosan_S(0.003)PS(0.997)YKEV'), T_S_GT1): ('Docosanoic_acid', 'S(0.003)PS(0.997)YKEVALAPPGSIAK'),
             # Scenarios with varying characteristics of corr., rank corr., IC
             # IC >= 7.25, but corr. / rank corr. low
             (('Z_GT0_1_Fructose_DHYDMYEPEANT'), T_S_GT0): ('Fructose', 'DHYDMYEPEANTMLQNS(1)APGS(1)PFHPAGSR'),
             (('Z_GT0_2_Glycine_NS(0.023)S(0.974)KDD'), T_S_GT0): ('Glycine', 'NS(0.023)S(0.974)KDDADRET(0.002)LEAEHTTLVTPHH'),
             # Other scenarios with varying characteristics of corr., rank corr., IC
             (('Z_GT0_3_Beta-ala_SLEELS(1)GEA'), T_S_GT0): ('Beta-alanine', 'SLEELS(1)GEAEVS(1)HDEK'),
             (('Z_GT0_4_Docosan_SLEELS(1)GEA'), T_S_GT0): ('Docosanoic_acid', 'SLEELS(1)GEAEVS(1)HDEK'),
             (('Z_GT0_5_Docosan_T(1)AILERR'), T_S_GT0): ('Docosanoic_acid', 'T(1)AILERR'),
             (('Z_GT0_6_Octadecan_HLEENGS(1)DGE'), T_S_GT0): ('Octadecanoic_acid', 'HLEENGS(1)DGEQGPGGSNGWITTINDVEMENQIVLPEDKK'),
             (('Z_GT0_7_Beta-ala_S(1)VLDT(1)PL'), T_S_GT0): ('Beta-alanine', 'S(1)VLDT(1)PLS(0.853)S(0.147)AR'),
             (('Z_GT0_8_Tetradecan_TVAAVAGS(1)PG'), T_S_GT0): ('Tetradecanoic_acid', 'TVAAVAGS(1)PGT(0.002)PT(0.963)T(0.035)PGSAR'),
             (('Z_GT0_9_Hexadecan_TVAAVAGS(1)PG'), T_S_GT0): ('Hexadecanoic_acid', 'TVAAVAGS(1)PGT(0.002)PT(0.963)T(0.035)PGSAR')}

# IC10_7p25_No_SBY_dGT10_0p2_No
# dPairsPaP = {(('Leu_SSFQEDHE'), T_S_GT): ('Leucine', 'SSFQEDHS(1)NIGGPGFSR'),
#              (('Beta-alan_SSFQEDHE'), T_S_GT): ('Beta-alanine', 'SSFQEDHS(1)NIGGPGFSR'),
#              (('Ala_VSS(1)AGL'), T_S_GT): ('Alanine', 'VSS(1)AGLRTESVLQR')}

# IC50_7p25_No_SBY_dGT50_0p2_No
# dPairsPaP = {(('Isoleu_ASGAGPN'), T_S_GT): ('Isoleucine', 'ASGAGPNSLVS(1)PQR'),
#              (('Isoleu_DYEDPPP'), T_S_GT): ('Isoleucine', 'DYEDPPPT(1)PFFDADELTK'),
#              (('Isoleu_DLDVNES(1)G'), T_S_GT): ('Isoleucine', 'DLDVNES(1)GPPAAR'),
#              (('Val_ASGAGPN'), T_S_GT): ('Valine', 'ASGAGPNSLVS(1)PQR'),
#              (('Val_DYEDPPP'), T_S_GT): ('Valine', 'DYEDPPPT(1)PFFDADELTK'),
#              (('Val_DLDVNES(1)G'), T_S_GT): ('Valine', 'DLDVNES(1)GPPAAR')}

nmPltPaP = S_NM_PAT_PLT         # name prefix of the pattern plot
tFigSzPaP = (3., 4.)            # (width, height): figure size [inches]
lWdPltPatPaP = 1.5              # line width of pattern in pattern plot
lWdPltXAxPaP = .5               # line width of x-axis in pattern plot

# --- graphics parameters / IC component plot ---------------------------------
dPairsICP = dPairsPaP           # dictionary of metabolite-phosphopeptide pairs
nmPltICP = S_NM_ICC_PLT         # name prefix of the IC component plot
tFigSzICP = (6., 4.)            # (width, height): figure size [inches]
wdthBar = 0.2                   # width of single bars
wdthGrp = 0.75                  # width of bar group
degRotXLbl = 90                 # degree rotation of x-labels

# --- names and paths of files and dirs ---------------------------------------
sFIn_IC_M_P = 'IC_Met_Pho'
sFIn_dGT_M = 'DistGT_Met'
sFIn_dGT_P = 'DistGT_Pho'

sFOutS = 'S_XIOvRep'
sFOutF = 'F_XIOvRep'

# IC_No_No_SBNY_dGT_No_No_GT0_MPS5Y_GT1_MPS5NY_GT5_MPS5NY
# sFIn_PaP = 'F_XIOvRep_IC_IC_GT0_MPS5Y_GT1_MPS5NY_GT5_MPS5NY_All_No_No_PSBNY_dGTM_dGT_No_No_dGTP_dGT_No_No'

# IC_No_No_SBNY_dGT_No_No_GT0_MPS5NY_GT1_MPS5Y_GT5_MPS5NY
# sFIn_PaP = 'F_XIOvRep_IC_IC_GT0_MPS5NY_GT1_MPS5Y_GT5_MPS5NY_All_No_No_PSBNY_dGTM_dGT_No_No_dGTP_dGT_No_No'

# IC_No_No_SBNY_dGT_No_No_GT0_MPS5NY_GT1_MPS5NY_GT5_MPS5Y
# sFIn_PaP = 'F_XIOvRep_IC_IC_GT0_MPS5NY_GT1_MPS5NY_GT5_MPS5Y_All_No_No_PSBNY_dGTM_dGT_No_No_dGTP_dGT_No_No'

# IC_GT0_12p0_No_GT1_11p0_No_GT5_10p0_No_dGT_No_No_MPS5Y_PSBNY
# sFIn_PaP = 'F_XIOvRep_IC_IC_GT0_12p0_No_GT1_11p0_No_GT5_10p0_No_All_MPS5Y_PSBNY_dGTM_dGT_No_No_dGTP_dGT_No_No'

# IC_GT0_12p0_No_GT1_11p0_No_GT5_10p0_No_dGT_No_No_MPS5NY_PSBNY
# sFIn_PaP = 'F_XIOvRep_IC_IC_GT0_12p0_No_GT1_11p0_No_GT5_10p0_No_All_MPS5NY_PSBNY_dGTM_dGT_No_No_dGTP_dGT_No_No'

# IC_7p25_No_dGT_No_No_MPS5NY_PSBNY
# sFIn_PaP = 'F_XIOvRep_IC_IC_All_7p25_No_MPS5NY_PSBNY_dGTM_dGT_No_No_dGTP_dGT_No_No'

# IC_No_No_SBNY_dGT_No_No
# sFIn_PaP = 'F_XIOvRep_IC_IC_All_No_No_MPS5Y_PSBNY_dGTM_dGT_No_No_dGTP_dGT_No_No'

# IC_No_No_dGT_No_No_MPS5NY_PSBNY
sFIn_PaP = 'F_XIOvRep_IC_IC_All_No_No_MPS5NY_PSBNY_dGTM_dGT_No_No_dGTP_dGT_No_No'

# IC_7p25_No_SBY_dGT_No_No
# sFIn_PaP = 'F_XIOvRep_IC_IC_All_7p25_No_MPS5Y_PSBY_dGTM_dGT_No_No_dGTP_dGT_No_No'

# IC_10_7p25_No_dGT_10_0p2_No
# sFIn_PaP = 'F_XIOvRep_IC_IC_GT0_7p25_No_GT1_7p25_No_GT5_No_No_All_MPS5Y_PSBY_dGTM_10_0p2_No_dGTP_10_0p2_No'

# IC_50_7p25_No_dGT_50_0p2_No
# sFIn_PaP = 'F_XIOvRep_IC_IC_GT0_7p25_No_GT1_No_No_GT5_7p25_No_All_MPS5Y_PSBY_dGTM_50_0p2_No_dGTP_50_0p2_No'

dSFIn_ICP = {sGT: ('ICCmp__BinOp_MetD_DvSD_' + sGT + '_AllD_PhoD_DvSD_'
                   + sGT + '_AllD') for sGT in L_S_GT}
sFOutPaP = nmPltPaP
sFOutICP = nmPltICP

sDirInCSV = 'InpDat'
sDirOutCSV = 'OutCSV'
sDirOutPaP = 'OutPDF'
sDirOutICP = sDirOutPaP

pBaseIn = os.path.join('..', '..', '..', '25_Papers', '01_FirstAuthor',
                       '04_SysBio_DataAnalysis', '15_Figures',
                       'B__Pattern_ICComp')
# pBaseOut = os.path.join('..', '..', '..', '25_Papers', '01_FirstAuthor',
#                         '04_SysBio_DataAnalysis', '15_Figures',
#                         'B__Pattern_ICComp')
pBaseOut = os.path.join('..', '..', '..', '..', '..', '..', 'W')
pInCSV = os.path.join(pBaseIn, sDirInCSV)
pOutCSV = os.path.join(pBaseOut, sDirOutCSV)

pInPaP = os.path.join(pBaseOut, sDirOutCSV)
pInICP = os.path.join(pBaseIn, sDirInCSV)
pOutPaP = os.path.join(pBaseOut, sDirOutPaP)
pOutICP = os.path.join(pBaseOut, sDirOutICP)

# --- derived values ----------------------------------------------------------
dMapCHdSel = {S_SG: {S_SG_M: {sGT: S_SG_MET for sGT in L_S_GT},
                     S_SG_P: {sGT: S_SG_PHO for sGT in L_S_GT}},
              S_SB: {S_SB_P: {sGT: S_USC.join([S_SELECTED, S_PHO[0]])
                              for sGT in L_S_GT}}}
dComprStr = {S_IC_P: S_IC_P_S,
             S_IC_N: S_IC_N_S,
             S_PEAR_CR: S_PEAR_CR_S,
             S_SPEAR_CR: S_SPEAR_CR_S,
             S_PEAR_P: S_PEAR_P_S,
             S_SPEAR_P: S_SPEAR_P_S,
             S_D_CL: S_D_CL_S,
             S_D_E: S_D_E_S,
             S_D_GT: S_D_GT_S}

dPFInIC = {sGT: os.path.join(pInCSV, sFIn_IC_M_P + S_USC + sGT + S_DOT + S_CSV)
           for sGT in L_S_GT}

# --- assertions --------------------------------------------------------------
assert set(dISort) == set(dThr)
for cD in [dISort, dThr]:
    assert S_IC in cD and set(cD[S_IC]) == set(L_S_GT)

# --- INPUT DICTIONARY --------------------------------------------------------
dInput = {# --- constants
          'sBC_L': S_BIN_L,
          'sBC2_L': S_BIN_L_2,
          'lSGT': L_S_GT,
          'sBase': S_BASE_CL,
          'sInpDat': S_INP_DATA,
          'sRoot': S_ROOT_CL,
          'sFNmCmp': S_F_NM_CMP,
          'sFNmCmpDGT': S_F_NM_CMP_DGT,
          'sFNmCmpICAllGT': S_F_NM_CMP_IC_ALL_GT,
          'sFNmCmpICSglGT': S_F_NM_CMP_IC_SGL_GT,
          'sFNmCnstr': S_F_NM_CNSTR,
          'sExtrInfo': S_EXTR_INFO,
          'sPltr': S_PLTR,
          'sPatPltr': S_PAT_PLTR,
          'sICCPltr': S_ICC_PLTR,
          'sMet': S_MET,
          'sPho': S_PHO,
          'R04': R04,
          # --- flow control
          'doInfoExtr': doInfoExtr,
          'doPlotPat': doPlotPat,
          'doPlotICC': doPlotICC,
          'dSpecSel': dSpecSel,
          # --- general input
          'modDisp': modDisp,
          # --- data specific input
          'sortInFNm': sortInFNm,
          'dUsedK': dUsedK,
          'dISort': dISort,
          'dThr': dThr,
          'dSel': dSel,
          'sSep': sSep,
          # --- graphics parameters / pattern plot
          'plotOfPatterns': {'dPairs': dPairsPaP,
                             'nmPlt': nmPltPaP,
                             'tFigSz': tFigSzPaP,
                             'szFontLeg': szFontLeg,
                             'nCharDsp': nCharDsp,
                             'posLegXY': posLegXY,
                             'lWdPltPat': lWdPltPatPaP,
                             'lWdPltXAx': lWdPltXAxPaP,
                             'clrDef': clrDef},
          # --- graphics parameters / IC component plot
          'plotOfICCmp': {'dPairs': dPairsICP,
                          'nmPlt': nmPltICP,
                          'tFigSz': tFigSzICP,
                          'szFontLeg': szFontLeg,
                          'nCharDsp': nCharDsp,
                          'posLegXY': posLegXY,
                          'lWdPlt': lWdPlt,
                          'wdthBar': wdthBar,
                          'wdthGrp': wdthGrp,
                          'degRotXLbl': degRotXLbl,
                          'clrDef': clrDef},
          # --- names and paths of files and dirs
          'pInCSV': pInCSV,
          'pOutCSV': pOutCSV,
          'pInPaP': pInPaP,
          'pInICP': pInICP,
          'pOutPaP': pOutPaP,
          'pOutICP': pOutICP,
          'dPFInIC': dPFInIC,
          'pFInM': os.path.join(pInCSV, sFIn_dGT_M + S_DOT + S_CSV),
          'pFInP': os.path.join(pInCSV, sFIn_dGT_P + S_DOT + S_CSV),
          'sFOutS': sFOutS,
          'sFOutF': sFOutF,
          'pFInPaP': os.path.join(pInPaP, sFIn_PaP + S_DOT + S_CSV),
          'dPFInICP': {sGT: os.path.join(pInICP, dSFIn_ICP[sGT] + S_DOT +
                                         S_CSV) for sGT in L_S_GT},
          'sFOutPaP': sFOutPaP,
          'sFOutICP': sFOutICP,
          # --- further derived values
          'dMapCHdSel': dMapCHdSel,
          'dComprStr': dComprStr}

# --- FUNCTIONS ---------------------------------------------------------------
def addToDictD(cD, cKMain, cKSub, cV):
    if cKMain in cD:
        assert cKSub not in cD[cKMain]
        cD[cKMain][cKSub] = cV
    else:
        cD[cKMain] = {}
        cD[cKMain][cKSub] = cV

def allElEq(iterator):
    iterator = iter(iterator)
    try:
        first = next(iterator)
    except StopIteration:
        return True
    return all(first == x for x in iterator)

def addIt(s='', sAdd='', sSep=S_USC):
    if len(s) > 0:
        if s[-1] != sSep:
            s += sSep
    return s + sAdd

def comprStr(s, dCmpr, sSep=S_USC):
    lSSplO, lSSplN = s.split(sSep), []
    for i, s1s in enumerate(lSSplO):
        s2s = s1s
        if i > 0:
            s2s = sSep.join([lSSplO[i - 1], s1s])
        if s2s in dCmpr:
            lSSplN[-1] = dCmpr[s2s]
        elif s1s in dCmpr:
            lSSplN.append(dCmpr[s1s])
        else:
            lSSplN.append(s1s)
    return sSep.join(lSSplN)

def num2StrF(cV, dRepl={S_DOT: S_P, S_DASH: S_M, 'None': S_NO}):
    cS = str(cV)
    for sK, sV in dRepl.items():
        cS = cS.replace(sK, sV)
    return cS

def getSHdCD2GT(sID, lSDef):
    for sF in lSDef:
        lSSpl = sF.split(S_USC)
        assert len(lSSpl) >= 2
        if lSSpl[-2][-1] + lSSpl[-1][-1] == sID:
            return sF
    return None

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

def sortDfr(pdDfr, dSrt, sSrtBy, sOrd, srtKind='stable'):
    isAsc, sHdC = (dSrt[sOrd] == S_ASC), dSrt[sSrtBy]
    if sHdC not in pdDfr.columns:
        sHdC = getSHdCD2GT(sHdC, L_S_D_2GT)
    pdDfr.sort_values(by=sHdC, ascending=isAsc, inplace=True, kind=srtKind)

def thrFilter(pdDfr, sHdC, thrMin, thrMax):
    if sHdC not in pdDfr.columns:
        sHdC = getSHdCD2GT(sHdC, L_S_D_2GT)
    if thrMin is None:
        if thrMax is None:
            return pdDfr
        else:
            return pdDfr[pdDfr[sHdC] <= thrMax]
    else:
        if thrMax is None:
            return pdDfr[pdDfr[sHdC] >= thrMin]
        else:
            return pdDfr[(pdDfr[sHdC] >= thrMin) & (pdDfr[sHdC] <= thrMax)]

def filterDfrSel(pdDfr, dGT, dMap, sKMain, sKSub, sGT):
    sHdC, lSSel = dMap[sKMain][sKSub][sGT], dGT[sGT]
    if sHdC in pdDfr.columns:
        return pdDfr[pdDfr[sHdC].isin(lSSel)]
    else:
        return pdDfr

def applySelFilter(pdDfr, dSel, dMap):
    for tK, dSub in dSel.items():
        for sKSub, dGT in dSub.items():
            if allElEq(dGT.values()) and allElEq(dMap[tK[0]][sKSub].values()):
                pdDfr = filterDfrSel(pdDfr, dGT, dMap, tK[0], sKSub, L_S_GT[0])
            else:
                for sGT in dGT:
                    pdDfr = filterDfrSel(pdDfr, dGT, dMap, tK[0], sKSub, sGT)
    return pdDfr

def getLNewIdx(pdDfr):
    assert (S_MET in pdDfr.columns) and (S_PHO in pdDfr.columns)
    return [S_USC.join([pdDfr.at[i, S_MET], pdDfr.at[i, S_PHO]])
            for i in pdDfr.index]

def modifyDfr(pdDfr, sGT):
    lHdCN = [s + S_USC + sGT for s in pdDfr.columns if s not in L_S_NO_GT]
    pdDfr.columns = L_S_NO_GT + lHdCN
    dfrMod = pd.concat([pd.Series(getLNewIdx(pdDfr), name=S_NEW_IDX),
                        pdDfr.reset_index(drop=True)],
                       axis=1, verify_integrity=True)
    return dfrMod.set_index(S_NEW_IDX, verify_integrity=True)

def concDfr(dDfr):
    dfrM = modifyDfr(dDfr[L_S_GT[0]], L_S_GT[0])
    for sGT in L_S_GT[1:]:
        dfrMR = modifyDfr(dDfr[sGT], sGT).drop(columns=L_S_NO_GT)
        dfrM = pd.concat([dfrM, dfrMR], axis=1, verify_integrity=True)
    return dfrM.drop(columns=L_S_NO_GT)

def appendToDDat(lDat, cDfrFl, sI, sHdC):
    cV = np.nan
    if sI in cDfrFl.index:
        cV = cDfrFl.at[sI, sHdC]
    lDat.append(cV)

def checkAllSetsEqSrtThr(dSTIC):
    lSets = [set(dSTIC[sGT].values()) for sGT in L_S_GT]
    return set.intersection(*lSets) == set.union(*lSets)

def checkMPSetsEqSl(dSlS, sGT):
    lSets = [set(dSlS[sKSlS][sGT]) for sKSlS in dSlS]
    return set.intersection(*lSets) == set.union(*lSets)

def checkAllSetsEqSl(dSlS):
    lSets = [set(dSlS[sKSlS][sGT]) for sKSlS in dSlS for sGT in L_S_GT]
    return set.intersection(*lSets) == set.union(*lSets)

def getSel(setSel, sIDSl):
    if setSel == {S_N, S_Y}:
        return sIDSl + S_N + S_Y
    elif  setSel == {S_Y}:
        return sIDSl + S_Y
    elif  setSel == {S_N}:
        return sIDSl + S_N
    return ''

def addSelAll(dSl, tKSl, sGT=S_GT0):
    sSubKSl, sIDSl = list(dSl[tKSl])[0], tKSl[1]
    return getSel(set(dSl[tKSl][sSubKSl][sGT]), sIDSl)

def saveDfrRes(dfrRes, dDat, pFOut, sSep, dSel=None, dMap=None, dropNA=True):
    dfrRes = dfrRes.append(pd.DataFrame(dDat), ignore_index=True,
                           verify_integrity=True)
    if dSel is not None and dMap is not None:
        dfrRes = applySelFilter(dfrRes, dSel, dMap)
    if dropNA:
        dfrRes.dropna(axis=0, how='any', inplace=True)
    dfrRes.reset_index(drop=True).to_csv(pFOut, sep=sSep)
    return dfrRes

def decorateClosePlot(cFig, cAx, dPlt, pPltF, sYLbl=''):
    cAx.set_ylabel(sYLbl)
    l = cAx.legend(loc='lower center', bbox_to_anchor=dPlt['posLegXY'],
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
        self.dUsedK = self.inpD.dUsedK
        self.dSort = self.inpD.dISort
        self.dT = self.inpD.dThr
        self.dSl = self.inpD.dSel
        self.dfrIn, self.dDfrIn = None, None
        print('Initiated "RootClass" base object.')

    def printObjInfo(self):
        print('-'*20, 'Object', self.descO, '(ID', self.idO, ')', '-'*20)
        print('-'*8, 'Input data:')
        self.inpD.printAttrData()
        print('-'*8, 'Attributes of', self.descO, 'class:')
        self.printAttrData()

    def printDfrInp(self):
        if self.dfrIn is None:
            print('Input DataFrame does not have any content yet.')
        else:
            print('Input DataFrame:')
            print(self.dfrIn)

    def printDDfrInp(self):
        if self.dDfrIn is None:
            print('Input DataFrames dictionary does not have any content yet.')
        else:
            print('Input DataFrames dictionary:')
            print(self.dDfrIn)

    def printDictsSortFilt(self):
        print('_'*8, 'Sorting dictionary:', '_'*8)
        for cK, cV in self.dSort.items():
            print('* ' + str(cK) + ':\n' + str(cV))
        print('_'*8, 'Threshold dictionary:', '_'*8)
        for cK, cV in self.dT.items():
            print('* ' + str(cK) + ':\n' + str(cV))
        print('_'*8, 'Selection dictionary:', '_'*8)
        for cK, cV in self.dSl.items():
            print('* ' + str(cK) + ':\n' + str(cV))

class FNmCmp(RootClass):
    def __init__(self, InpD, dSetsEq):
        super().__init__(InpD)
        self.idO = InpD.sFNmCmp
        self.descO = 'File name component'
        self.dSetsEq = dSetsEq
        print('Initiated "FNmCmp" base object.')

class FNmCmpICSglGT(FNmCmp):
    def __init__(self, InpD, dSetsEq, sK):
        super().__init__(InpD, dSetsEq)
        self.idO = InpD.sFNmCmpICSglGT
        self.descO = 'File name component for IC valid for all GTs'
        self.buildCmpIC(sK=sK)
        print('Initiated "FNmCmpICSglGT" base object.')

    def buildCmpSrt(self, sK=S_IC):
        if self.inpD.sortInFNm:
            if not self.dSetsEq[S_SRT]:
                for sGT in self.dCmp:
                    sKSrt = S_USC.join([self.dSort[sK][sGT][S_ORD],
                                        self.dSort[sK][sGT][S_SRT_BY]])
                    self.dCmp[sGT] = addIt(self.dCmp[sGT], sKSrt)

    def buildCmpThr(self, sK=S_IC):
        if not self.dSetsEq[S_THR]:
            for sGT in self.dCmp:
                for sMM in L_S_MIN_MAX:
                    self.dCmp[sGT] = addIt(self.dCmp[sGT],
                                           num2StrF(self.dT[sK][sGT][sMM]))

    def buildCmpSel(self):
        for tSl in self.dSl:
            if not self.dSetsEq[S_SEL][tSl]:
                dDoSub = {sGT: True for sGT in self.dCmp}
                for sKSub in self.dSl[tSl]:
                    for sGT in self.dCmp:
                        if dDoSub[sGT]:
                            setSl, sID = set(self.dSl[tSl][sKSub][sGT]), tSl[1]
                            if checkMPSetsEqSl(self.dSl[tSl], sGT):
                                dDoSub[sGT] = False
                            else:
                                sID = sKSub
                            self.dCmp[sGT] = addIt(self.dCmp[sGT],
                                                   getSel(setSl, sID))

    def assembleSCmp(self):
        for sGT in self.dCmp:
            if len(self.dCmp[sGT]) > 0:
                self.dCmp[sGT] = sGT + S_USC + self.dCmp[sGT]
                self.sCmp = addIt(self.sCmp, self.dCmp[sGT])

    def buildCmpIC(self, sK=S_IC):
        self.sCmp = S_USC.join([sK, self.dUsedK[sK]])
        self.dCmp = {sGT: '' for sGT in L_S_GT}
        self.buildCmpSrt(sK=sK)
        self.buildCmpThr(sK=sK)
        self.buildCmpSel()
        self.assembleSCmp()

class FNmCmpICAllGT(FNmCmp):
    def __init__(self, InpD, dSetsEq, sK):
        super().__init__(InpD, dSetsEq)
        self.idO = InpD.sFNmCmpICAllGT
        self.descO = 'File name component for IC valid for all GTs'
        self.buildCmpIC(sK=sK)
        print('Initiated "FNmCmpICAllGT" base object.')

    def buildCmpSrt(self, sK=S_IC, sGT=S_GT0):
        if self.inpD.sortInFNm:
            if self.dSetsEq[S_SRT]:
                sKSrt = S_USC.join([self.dSort[sK][sGT][S_ORD],
                                    self.dSort[sK][sGT][S_SRT_BY]])
                self.sCmp = addIt(self.sCmp, sKSrt)

    def buildCmpThr(self, sK=S_IC, sGT=S_GT0):
        if self.dSetsEq[S_THR]:
            for sMM in L_S_MIN_MAX:
                self.sCmp = addIt(self.sCmp, num2StrF(self.dT[sK][sGT][sMM]))

    def buildCmpSel(self):
        for tSl in self.dSl:
            if self.dSetsEq[S_SEL][tSl]:
                self.sCmp = addIt(self.sCmp, addSelAll(self.dSl, tSl))
        if len(self.sCmp) > 0:
            self.sCmp = addIt(S_ALL, self.sCmp)

    def buildCmpIC(self, sK=S_IC):
        self.sCmp = ''
        self.buildCmpSrt(sK=sK)
        self.buildCmpThr(sK=sK)
        self.buildCmpSel()

class FNmCmpDGT(FNmCmp):
    def __init__(self, InpD, dSetsEq, lSK):
        super().__init__(InpD, dSetsEq)
        self.idO = InpD.sFNmCmpDGT
        self.descO = 'File name component for dGT'
        self.sCmp = ''
        for i, sK in enumerate(lSK):
            self.buildCmpDGT(sK, i < (len(lSK) - 1))
        print('Initiated "FNmCmpDGT" base object.')

    def buildCmpDGT(self, sK, addUSC):
        if self.inpD.sortInFNm:
            sKSrt = S_USC.join([self.dSort[sK][S_ORD],
                                self.dSort[sK][S_SRT_BY]])
            self.sCmp += S_USC.join([sK, self.dUsedK[sK], sKSrt])
        else:
            self.sCmp += S_USC.join([sK, self.dUsedK[sK]])
        for sMM in L_S_MIN_MAX:
            self.sCmp += S_USC + num2StrF(self.dT[sK][sMM])
        if addUSC:
            self.sCmp += S_USC

class FileNameConstructor(RootClass):
    def __init__(self, InpD, tpFNm='F'):
        super().__init__(InpD)
        self.idO = InpD.sFNmCnstr
        self.descO = 'File name constructor'
        self.constrFNm(tpFNm=tpFNm)
        print('Initiated "FileNameConstructor" base object.')

    def checkSetsEq(self):
        self.dSetsEq = {S_SRT: checkAllSetsEqSrtThr(self.dSort[S_IC]),
                        S_THR: checkAllSetsEqSrtThr(self.dT[S_IC]),
                        S_SEL: {}}
        for tSl, dSlSub in self.dSl.items():
            self.dSetsEq[S_SEL][tSl] = checkAllSetsEqSl(dSlSub)

    def joinCmpToFNm(self):
        self.sFNm = self.sFOut
        l = [self.oCmpICSglGT.sCmp, self.oCmpICAllGT.sCmp, self.oCmpDGT.sCmp]
        for s in l:
            if len(s) > 0:
                self.sFNm += S_USC + comprStr(s, self.inpD.dComprStr)

    def constrFNm(self, tpFNm='F'):
        self.sFOut = self.inpD.sFOutF
        if tpFNm == 'S':
            self.sFOut = self.inpD.sFOutS
        self.checkSetsEq()
        self.oCmpICSglGT = FNmCmpICSglGT(self.inpD, self.dSetsEq, sK=S_IC)
        self.oCmpICAllGT = FNmCmpICAllGT(self.inpD, self.dSetsEq, sK=S_IC)
        if len(self.oCmpICSglGT.sCmp) == 0:
            self.oCmpICAllGT.sCmp = S_USC.join([self.dUsedK[S_IC],
                                                self.oCmpICAllGT.sCmp])
        self.oCmpDGT = FNmCmpDGT(self.inpD, self.dSetsEq, lSK=L_S_D_GT)
        self.joinCmpToFNm()

class ExtractedInfo(RootClass):
    def __init__(self, InpD):
        super().__init__(InpD)
        self.idO = InpD.sExtrInfo
        self.descO = 'Extracted info'
        self.dMap = self.inpD.dMapCHdSel
        self.sSp = self.inpD.sSep
        self.getProcData()
        print('Initiated "ExtractedInfo" base object.')

    def getPResF(self):
        cnstrFNmS = FileNameConstructor(self.inpD, tpFNm='S')
        cnstrFNmF = FileNameConstructor(self.inpD, tpFNm='F')
        pOut = self.inpD.pOutCSV
        self.pFOutS = os.path.join(pOut, S_DOT.join([cnstrFNmS.sFNm, S_CSV]))
        self.pFOutF = os.path.join(pOut, S_DOT.join([cnstrFNmF.sFNm, S_CSV]))
        self.dSpcSel = {'S': {'dfr': None, 'pF': self.pFOutS},
                        'F': {'dfr': None, 'pF': self.pFOutF}}

    def getInf4Inp(self):
        dDatTp_IC = {sIn: str for sIn in [self.inpD.sBC_L, self.inpD.sBC2_L]}
        dDatTp_P = {sIn: str for sIn in [self.inpD.sBC2_L]}
        self.dDatTp = {S_IC: dDatTp_IC, S_D_GT_P: dDatTp_P}
        self.dPFIn = {S_IC: self.inpD.dPFInIC,
                      S_D_GT_M: self.inpD.pFInM,
                      S_D_GT_P: self.inpD.pFInP}

    def loadDDfrInp(self):
        # load input DataFrames
        dDfrIn_IC = {sGT: pd.read_csv(self.dPFIn[S_IC][sGT], sep=self.sSp,
                                      dtype=self.dDatTp[S_IC])
                     for sGT in L_S_GT}
        dfrIn_M = pd.read_csv(self.dPFIn[S_D_GT_M], sep=self.sSp)
        dfrIn_P = pd.read_csv(self.dPFIn[S_D_GT_P], sep=self.sSp,
                              dtype=self.dDatTp[S_D_GT_P])
        self.dDfrIn = {S_IC: dDfrIn_IC,
                       S_D_GT_M: dfrIn_M,
                       S_D_GT_P: dfrIn_P}

    def getDHdCol(self):
        # get dictionary of column headers
        dDfr_IC = self.dDfrIn[S_IC]
        dHdCol_IC = {sGT: [sC + S_USC + sGT for sC in
                           list(dDfr_IC[sGT].loc[:, S_RMNG_COL1_IC:].columns)]
                     for sGT in L_S_GT}
        self.dHdCol = {S_IC: dHdCol_IC,
                       S_D_GT_M: list(self.dDfrIn[S_D_GT_M].columns),
                       S_D_GT_P: list(self.dDfrIn[S_D_GT_P].columns)}

    def getDMapK(self):
        # get mapping of data dictionary keys to column headers
        self.dMapK = {}
        for i, sMP in enumerate(L_S_M_P):
            self.dMapK[sMP] = (L_S_D_GT[i], sMP)
            lHdCRed = [s for s in self.dHdCol[L_S_D_GT[i]] if s not in L_S_M_P]
            for sHdC in lHdCRed:
                self.dMapK[S_USC.join([sHdC, sMP[0]])] = (L_S_D_GT[i], sHdC)
        for sGT in L_S_GT:
            for s in self.dHdCol[S_IC][sGT]:
                self.dMapK[s] = (S_IC, s)

    def getProcData(self):
        self.getPResF()
        self.getInf4Inp()
        self.loadDDfrInp()
        self.getDHdCol()
        self.getDMapK()

    def filterDfr(self, sKey, sMin=S_MIN, sMax=S_MAX):
        thMin, thMax = self.dT[sKey][sMin], self.dT[sKey][sMax]
        return thrFilter(self.dDfrIn[sKey], self.dUsedK[sKey], thMin, thMax)

    def filterAndConc(self, sKey, sMin=S_MIN, sMax=S_MAX):
        # filter data
        self.dDfrFl = {sKey: {}}
        for sGT in L_S_GT:
            thMin, thMax = self.dT[sKey][sGT][sMin], self.dT[sKey][sGT][sMax]
            self.dDfrFl[sKey][sGT] = thrFilter(self.dDfrIn[sKey][sGT],
                                               self.dUsedK[sKey], thMin, thMax)
        # process data - concatenate IC DataFrames of the three GT
        return concDfr(self.dDfrFl[sKey])

    def sortAndFiltDfr(self):
        for sGT, cDfr in self.dDfrIn[S_IC].items():
            sortDfr(cDfr, self.dSort[S_IC][sGT], S_SRT_BY, S_ORD)
        sortDfr(self.dDfrIn[S_D_GT_M], self.dSort[S_D_GT_M], S_SRT_BY, S_ORD)
        sortDfr(self.dDfrIn[S_D_GT_P], self.dSort[S_D_GT_P], S_SRT_BY, S_ORD)
        self.dDfrFl = {S_IC: self.filterAndConc(S_IC),
                       S_D_GT_M: self.filterDfr(S_D_GT_M),
                       S_D_GT_P: self.filterDfr(S_D_GT_P)}
        self.N = self.dDfrFl[S_D_GT_M].shape[0]*self.dDfrFl[S_D_GT_P].shape[0]

    def iniDfrRes(self, specSel=None):
        lC = list(self.dMapK)
        if specSel == 'S':  # the "short" subset of data
            # select only a subset of all keys (columns of DataFrames): lC
            lC = [S_BIN_L_2, S_SELECTED] + L_S_IC_GT + L_S_D_GT
            lC = [S_USC.join([s, S_PHO[0]]) for s in [S_BIN_L_2, S_SELECTED]]
            lC = L_S_M_P + lC + L_S_IC_GT
            lC += [S_USC.join([S_D_GT, s[0]]) for s in L_S_M_P]
            self.dfrResS = pd.DataFrame(columns=lC)
            self.dSpcSel[specSel]['dfr'] = self.dfrResS
        else:
            # select all keys (columns of DataFrames) defined in self.dMapK
            self.dfrResF = pd.DataFrame(columns=lC)
            self.dSpcSel[specSel]['dfr'] = self.dfrResF
        return {cK: [] for cK in lC}

    def SUB_fillDDat(self, dDat, i, j, sIMP):
        for sKDDt, lDt in dDat.items():
            assert sKDDt in self.dMapK
            (sKDFl, sHdC) = self.dMapK[sKDDt]
            if sKDFl in L_S_D_GT:
                assert sKDFl in self.dDfrFl
                assert sHdC in self.dDfrFl[sKDFl].columns
                if sKDFl == S_D_GT_M:
                    dDat[sKDDt].append(self.dDfrFl[sKDFl].at[i, sHdC])
                else:
                    dDat[sKDDt].append(self.dDfrFl[sKDFl].at[j, sHdC])
            elif sKDFl in [S_IC]:
                appendToDDat(lDt, self.dDfrFl[S_IC], sIMP, sHdC)
            else:
                print('ERROR: Unknown key of filter dictionary:', sKDFl)
                assert False

    def fillDDat(self, dDat):
        n = 0
        for i in self.dDfrFl[S_D_GT_M].index:
            sMet = self.dDfrFl[S_D_GT_M].at[i, S_MET]
            for j in self.dDfrFl[S_D_GT_P].index:
                sPho = self.dDfrFl[S_D_GT_P].at[j, S_PHO]
                sIMP = S_USC.join([sMet, sPho])
                self.SUB_fillDDat(dDat, i, j, sIMP)
                n += 1
                if n%self.inpD.modDisp == 0:
                    print('Processed element', n, 'of', self.N, '.')

    def fillSaveDfrRes(self):
        for spcSel, dISpcSel in self.inpD.dSpecSel.items():
            dDat = self.iniDfrRes(specSel=spcSel)
            self.fillDDat(dDat)
            print('Filled data dictionary for selection "' + spcSel + '".')
            dISel, dropNA = self.dSpcSel[spcSel], dISpcSel['dropNA']
            dISel['dfr'] = saveDfrRes(dISel['dfr'], dDat, dISel['pF'],
                                      self.sSp, self.dSl, self.dMap, dropNA)

    def extractionOfExtremes(self):
        self.sortAndFiltDfr()
        self.fillSaveDfrRes()

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

    def loadDfrInp(self, iC=0):
        dDatTp = {sIn: str for sIn in L_S_SG_GT}
        # load input DataFrame
        if hasattr(self, 'pFIn'):
            self.dfrIn = pd.read_csv(self.pFIn, sep=self.sSp, index_col=iC,
                                     dtype=dDatTp)

    def loadDDfrInp(self):
        # load input DataFrames, and save them in dictionary
        if hasattr(self, 'dPFIn'):
            self.dDfrIn = {sGT: None for sGT in L_S_GT}
            for sGT in self.dDfrIn:
                self.dDfrIn[sGT] = pd.read_csv(self.dPFIn[sGT], sep=self.sSp)

    def iniPlot(self, dPlt):
        cFig, cAx = plt.subplots()
        if dPlt['tFigSz'] is not None and len(dPlt['tFigSz']) == 2:
            cFig.set_size_inches(dPlt['tFigSz'])
        return cFig, cAx

class PatternPlotter(Plotter):
    def __init__(self, InpD):
        super().__init__(InpD)
        self.idO = self.inpD.sPatPltr
        self.descO = 'Pattern plotter'
        self.pFIn = self.inpD.pFInPaP
        self.pDOut = self.inpD.pOutPaP
        self.dPlt = self.inpD.plotOfPatterns
        self.getDPPltF()
        self.loadDfrInp()
        print('Initiated "PatternPlotter" base object and loaded input data.')

    def getDPPltF(self):
        self.dPPltF, sFPlt = {}, self.inpD.sFOutPaP
        for ((s1, tSGT), tMP) in self.dPlt['dPairs'].items():
            for sGT in tSGT:
                sPltF = S_DOT.join([S_USC.join([sFPlt, s1, sGT]), S_PDF])
                self.dPPltF[(s1, sGT)] = (tMP, os.path.join(self.pDOut, sPltF))

    def plotPatterns(self):
        d, nChD = self.dfrIn, self.dPlt['nCharDsp']
        serXAx = pd.Series([0.]*len(L_S_FT), index=L_S_FT)
        for ((s1, sGT), ((sM, sP), pPltF)) in self.dPPltF.items():
            print('Plotting pattern for "' + s1 + '" and', sGT, '...')
            cSer = d[(d[S_MET] == sM) & (d[S_PHO] == sP)].squeeze()
            lTDat = [(S_MET, sM), (S_PHO, sP)]
            # if not os.path.isfile(pPltF):
            cFig, cAx = self.iniPlot(self.dPlt)
            for tMP in lTDat:
                cPa = cSer.loc[D_HD_C_PA[tMP[0]][sGT]]
                cPa.index = L_S_FT
                cAx.plot(cPa, lw=self.dPlt['lWdPltPat'], label=tMP[1][:nChD])
            cAx.plot(serXAx, lw=self.dPlt['lWdPltXAx'],
                     color=self.dPlt['clrDef'])
            decorateClosePlot(cFig, cAx, self.dPlt, pPltF, S_YLBL_PAT_PLT)

class ICCmpPlotter(Plotter):
    def __init__(self, InpD):
        super().__init__(InpD)
        self.idO = self.inpD.sICCPltr
        self.descO = 'Concordance index component plotter'
        self.dPFIn = self.inpD.dPFInICP
        self.pDOut = self.inpD.pOutICP
        self.dPlt = self.inpD.plotOfICCmp
        self.getDPPltF()
        self.loadDDfrInp()
        print('Initiated "ICCmpPlotter" base object and loaded input data.')

    def getDPPltF(self):
        self.dPPltF, sFPlt = {}, self.inpD.sFOutICP
        for ((s1, tSGT), tMP) in self.dPlt['dPairs'].items():
            for sGT in tSGT:
                sPltF = S_DOT.join([S_USC.join([sFPlt, s1, sGT]), S_PDF])
                self.dPPltF[(s1, sGT)] = (tMP, os.path.join(self.pDOut, sPltF))

    def plotICCmp(self):
        nChD, wdBar = self.dPlt['nCharDsp'], self.dPlt['wdthBar']
        for ((s1, sGT), ((sM, sP), pPltF)) in self.dPPltF.items():
            print('Plotting IC components for "' + s1 + '" and', sGT, '...')
            d = self.dDfrIn[sGT]
            cSer = d[(d[S_MET] == sM) & (d[S_PHO] == sP)].squeeze()
            lTDat = [(S_MET, sM), (S_PHO, sP), (S_ICC, S_ICC_LEG)]
            xLocG = np.arange(len(L_S_FT_CHG))
            # if not os.path.isfile(pPltF):
            cFig, cAx = self.iniPlot(self.dPlt)
            for k, tMP in enumerate(lTDat):
                cICC = cSer.loc[D_HD_C_ICC[tMP[0]]]
                cICC.index = L_S_FT_CHG
                xLoc = xLocG + (2*k + 1)/(2*len(lTDat))*self.dPlt['wdthGrp']
                cAx.bar(xLoc - 1/2, height=cICC, width=wdBar,
                        lw=self.dPlt['lWdPlt'], label=tMP[1][:nChD])
            cAx.plot([-1/2, len(L_S_FT_CHG) + 1/2], [0, 0],
                     lw=self.dPlt['lWdPlt'], color=self.dPlt['clrDef'])
            cAx.set_xticks(xLocG)
            cAx.set_xticklabels(L_S_FT_CHG)
            for cXLbl in cAx.get_xticklabels():
                cXLbl.set_rotation(self.dPlt['degRotXLbl'])
            decorateClosePlot(cFig, cAx, self.dPlt, pPltF, S_YLBL_ICC_PLT)

# --- MAIN --------------------------------------------------------------------
startTime = time.time()
print('+'*25 + ' START', time.ctime(startTime), '+'*25)

inpDat = InputData(dInput)
if inpDat.doInfoExtr:
    cXtrInfo = ExtractedInfo(inpDat)
    cXtrInfo.printIDDesc()
    cXtrInfo.printDfrInp()
    cXtrInfo.printDDfrInp()
    cXtrInfo.printAttrData()
    cXtrInfo.printDictsSortFilt()
    # cXtrInfo.printObjInfo()
    cXtrInfo.extractionOfExtremes()
if inpDat.doPlotPat:
    cPltr = PatternPlotter(inpDat)
    # cPltr.printAttrData()
    # cPltr.printDPPltF()
    cPltr.plotPatterns()
if inpDat.doPlotICC:
    cPltr = ICCmpPlotter(inpDat)
    cPltr.plotICCmp()

print('-'*80)
printElapsedTimeSim(startTime, time.time(), 'Total time')
print('+'*37, 'DONE.', '+'*36)

###############################################################################

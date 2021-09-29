# -*- coding: utf-8 -*-
###############################################################################
# --- C_00__GenConstants.py ---------------------------------------------------
###############################################################################

# --- Information w.r.t. genotypes / features investigated --------------------
NM_GT0 = 'GT0'
NM_GT1 = 'GT1'
NM_GT2 = 'GT2'
NM_GT3 = 'GT3'
NM_GT4 = 'GT4'
NM_GT5 = 'GT5'
L_NM_GT = [NM_GT0, NM_GT1, NM_GT2, NM_GT3, NM_GT4, NM_GT5]

NM_FT1 = 'FT1'
NM_FT2 = 'FT2'
NM_FT3 = 'FT3'
NM_FT4 = 'FT4'
L_NM_FT = [NM_FT1, NM_FT2, NM_FT3, NM_FT4]

T_NM_GT0 = (NM_GT0, 'WT', 'wild type', 'wt')
T_NM_GT1 = (NM_GT1, 'PGM', 'PGM mutant', 'pgm')
T_NM_GT2 = (NM_GT2, 'SIRK1', 'SIRK1 mutant', 'sirk1')
T_NM_GT3 = (NM_GT3, 'SIRK1-PGM', 'PGM/SIRK1 double mutant', 'pgm|sirk1')
T_NM_GT4 = (NM_GT4, 'SIRK1-SWEET', 'SIRK1/SWEET double mutant',
            'pgm|sweet11/12')
T_NM_GT5 = (NM_GT5, 'SWEET', 'SWEET mutant', 'sweet11/12')
L_T_NM_GT = [T_NM_GT0, T_NM_GT1, T_NM_GT2, T_NM_GT3, T_NM_GT4, T_NM_GT5]
D_NM_GT = {nmGT: L_T_NM_GT[k] for k, nmGT in enumerate(L_NM_GT)}

T_NM_FT1 = (NM_FT1, 'DR', 'day root', 'dr')
T_NM_FT2 = (NM_FT2, 'DS', 'day shoot', 'ds')
T_NM_FT3 = (NM_FT3, 'NR', 'night root', 'nr')
T_NM_FT4 = (NM_FT4, 'NS', 'night shoot', 'ns')
L_T_NM_FT = [T_NM_FT1, T_NM_FT2, T_NM_FT3, T_NM_FT4]
D_NM_FT = {nmFt: L_T_NM_FT[k] for k, nmFt in enumerate(L_NM_FT)}

NM_ALL_GT = 'AllGT'
NM_ALL_FT = 'AllFt'
NM_ALL_GT_FT = 'All'

NM_X = 'X'

NM_XFEAT_C = 'C'
NM_XFEAT_G = 'G'

NM_TYPEX_S = 'S'
NM_TYPEX_F = 'F'

# --- Information w.r.t. objects ----------------------------------------------
S_BASE = 'DtBC'
S_EXP_D = 'ExpD'
S_MET_D = 'MetD'
S_PHO_D = 'PhoD'
S_COMB_D = 'CombD'
S_EXP_DX = 'ExpDX'
S_MET_DX = 'MetDX'
S_PHO_DX = 'PhoDX'
S_CLR_D = 'ClR'
S_BIN_OP = 'BinOp'
S_CLUST = 'Cl'
S_KM_CLUST = 'KMeansCl'
S_AG_CLUST = 'AggloCl'
S_OV_REP = 'OvRep'

S_BASE_L = 'DataBaseClass'
S_EXP_D_L = 'ExpData'
S_MET_D_L = 'MetData'
S_PHO_D_L = 'PhoData'
S_COMB_D_L = 'CombData'
S_EXP_DX_L = 'ExpDataX'
S_MET_DX_L = 'MetDataX'
S_PHO_DX_L = 'PhoDataX'
S_CLR_D_L = 'ClResData'
S_BIN_OP_L = 'BinaryOps'
S_CLUST_L = 'Clustering'
S_OV_REP_L = 'OverRep'

S_2L = '2L'
S_UN = 'Un'

I_EXP_D_2L = 1
I_MET_D_2L = 11
I_MET_D_UN = 12
I_PHO_D_2L = 21
I_PHO_D_UN = 22
I_COMB_D = 31
I_EXP_DX = 41
I_CLR_D = 51
I_BIN_OP = 81
I_CLUST = 82
I_OV_REP = 83

L_S_BASIC_D_L_IN = [s for s in [S_EXP_D_L, S_MET_D_L, S_PHO_D_L]]
L_S_BASIC_D_L_C = [s + NM_XFEAT_C for s in [S_EXP_D_L, S_MET_D_L, S_PHO_D_L]]
L_S_BASIC_D_L_G = [s + NM_XFEAT_G for s in [S_EXP_D_L, S_MET_D_L, S_PHO_D_L]]

L_S_COMB_D_L_C = [S_COMB_D_L + NM_XFEAT_C]
L_S_COMB_D_L_G = [S_COMB_D_L + NM_XFEAT_G]

L_S_MERGED_D_L_C = L_S_COMB_D_L_C + [S_EXP_DX_L]
L_S_MERGED_D_L_G = L_S_COMB_D_L_G + [S_EXP_DX_L]

L_S_CLR_D_L_C = [S_CLR_D_L + NM_XFEAT_C]
L_S_CLR_D_L_G = [S_CLR_D_L + NM_XFEAT_G]

L_S_CALC_D_L = [S_BIN_OP_L, S_CLUST_L, S_OV_REP_L]

# --- names of files and dirs -------------------------------------------------
NM_OBJINP = 'ObjInput'
NMF_OBJINP_PRE = 'D_'
T_NM_RES_FIG_2L = ('21', '41')
T_NM_RES_FIG_UN = ('22', '42')
T_ID_RES_FIG = ('R', 'F')

NM_EXT_CSV = 'csv'
NM_EXT_PY = 'py'
NM_EXT_PDF = 'pdf'

# --- names -------------------------------------------------------------------
NM_E_STD = 'AllD'
NM_I_FTR_BASE = 'iFtrBase'
NM_TRANS_D = 'transD'
NM_STD_OP = 'stdOp'
NM_DEV_TP = 'devTp'
NM_CL_RES = S_CLR_D
NM_CLUSTERS = S_CLUST
NM_MEANS = 'Means'
NM_SDS = 'SDs'
ID_O_ATTRDFR = 'idO'
ID_GT_ATTRDFR = 'idGT'
NM_GT_ATTRDFR = 'nmGT'
ID_FT_ATTRDFR = 'idFt'
NM_FT_ATTRDFR = 'nmFt'
NM_WT_ATTRDFR = 'wtDS'
NM_PRE_TRANS = 'Tr'
NM_PRE_DEV = 'Dv'

# --- info for result files ---------------------------------------------------
NM_ID = 'ID'
NM_CT = 'Ct'
NM_TRANS_DATA = 'TrDt'
NM_TYPEX = 'typeX'
NM_ID_GT_FT = 'GTFt'
NM_NUM_CLUST = 'numClust'
NM_FILTER = 'Filter'

# --- other strings -----------------------------------------------------------
S_USC = '_'

S_PROT = 'Protein'
S_BIN_C_G = 'BinCode'
S_BIN_C_3 = 'BinCode3'
S_BIN_C_2 = 'BinCode2'
S_BIN_C_1 = 'BinCode1'
L_S_PHO_CL = [S_PROT, S_BIN_C_G, S_BIN_C_3, S_BIN_C_2, S_BIN_C_1]

S_MN_CONC = 'MeanConc'

S_IDX = 'Idx'
S_COL = 'Col'
S_NEG = 'Neg'
S_POS = 'Pos'
S_AV = 'Av'
S_TOP = 'Top'
S_CORR_S = 'Corr'
S_CORR_L = 'Correlation'
S_CORR_V = 'CorrV'
S_CORR_P = 'CorrP'
S_SPEAR_V = 'SpearV'
S_SPEAR_P = 'SpearP'
S_KEND_V = 'KendV'
S_KEND_P = 'KendP'
S_CORR_TTL = 'Correlation'
S_SPEAR_TTL = 'Rank Correlation'
S_DV_SC = 'DvSc'
S_CI = 'CI'
S_OCC_CI = 'OccCI'
S_N_OCC_ABS = 'nOccurrAbs'
S_P_VAL_OV = 'pValFOver'
S_P_VAL_UN = 'pValFUnder'
S_P_OF = 'p'
S_M_CORR_BON = 'Bonferroni'
S_PD = 'PD'
S_BO = 'BO'
S_N_OCC = 'NOcc'
S_OVER_REP = 'ORp'
S_UNDER_REP = 'URp'
S_Y_N_OCC = 'Number of occurrences ($\it{k}$)'
S_Y_P_VAL = 'over-representation score\n' + '-' + '$\log_{10}$' + '(p-value)'

# --- constants for special plots ---------------------------------------------
L_BIN_C_2 = ['1.1', '1.2', '1.3', '2.1', '2.2', '3.1', '3.2', '3.4', '3.5',
             '3.6', '3.99', '4.1', '4.2', '4.3', '5.3', '6.4', '6.5', '7.1',
             '8.1', '8.2', '8.3', '10.1', '10.2', '10.3', '10.5', '10.8',
             '11.1', '11.3', '11.8', '11.9', '12.1', '12.2', '12.4', '13.1',
             '13.2', '15', '15.1', '15.2', '16.2', '16.5', '16.8', '17.1',
             '17.2', '17.3', '17.5', '17.6', '18.2', '18.4', '20', '20.1',
             '20.2', '21.1', '21.2', '21.4', '21.5', '21.6', '23.1', '23.2',
             '23.3', '23.4', '24.1', '25', '25.6', '26.13', '26.16', '26.17',
             '26.22', '26.3', '26.4', '26.5', '26.6', '26.7', '26.8', '26.9',
             '27.1', '27.2', '27.3', '27.4', '28.1', '28.2', '28.99', '29.1',
             '29.2', '29.3', '29.4', '29.5', '29.6', '30.1', '30.11', '30.2',
             '30.3', '30.4', '30.5', '30.6', '30.7', '30.8', '31.1', '31.2',
             '31.3', '31.4', '31.5', '32', '33.1', '33.99', '34.1', '34.11',
             '34.12', '34.13', '34.14', '34.15', '34.16', '34.17', '34.18',
             '34.19', '34.2', '34.21', '34.3', '34.4', '34.5', '34.7', '34.8',
             '34.9', '34.98', '34.99', '35.1', '35.2']
L_CLR_CYC_RGB = [(0.12, 0.47, 0.71), (1.0, 0.5, 0.05), (0.17, 0.63, 0.17),
                 (0.84, 0.15, 0.16), (0.58, 0.4, 0.74), (0.55, 0.34, 0.29),
                 (0.89, 0.47, 0.76), (0.5, 0.5, 0.5), (0.74, 0.74, 0.13),
                 (0.09, 0.75, 0.81)]

# --- other constants ---------------------------------------------------------
M_DETER = 'Deterministic'
M_STOCH = 'Stochastic'

R04 = 4
R06 = 6

FIGURE_MAX_OPEN_WARNING = 200

# --- derived values ----------------------------------------------------------
L_S_BASIC_D_L = L_S_BASIC_D_L_IN + L_S_BASIC_D_L_C + L_S_BASIC_D_L_G
L_S_USED_D_L_C = L_S_BASIC_D_L_C[1:]
L_S_USED_D_L = L_S_USED_D_L_C + L_S_BASIC_D_L_G[1:]
L_S_COMB_D_L = L_S_COMB_D_L_C + L_S_COMB_D_L_G
L_S_MERGED_D_L = L_S_MERGED_D_L_C + L_S_MERGED_D_L_G
L_S_CLR_D_L = L_S_CLR_D_L_C + L_S_CLR_D_L_G
L_S_BASIC_COMB_D_L = (L_S_BASIC_D_L_C + L_S_COMB_D_L_C +
                      L_S_BASIC_D_L_G + L_S_COMB_D_L_C)
L_S_USED_COMB_D_L = (L_S_BASIC_D_L_C[1:] + L_S_COMB_D_L_C +
                     L_S_BASIC_D_L_G[1:] + L_S_COMB_D_L_G)
S_SNDS = S_NEG + S_DV_SC
S_SPDS = S_POS + S_DV_SC
S_SPNDS = S_POS + S_NEG + S_DV_SC
S_NCI = S_NEG + S_CI
S_PCI = S_POS + S_CI
S_PNCI = S_POS + S_NEG + S_CI

###############################################################################

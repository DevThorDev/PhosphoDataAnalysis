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

NM_FT1 = 'FT1'
NM_FT2 = 'FT2'
NM_FT3 = 'FT3'
NM_FT4 = 'FT4'

T_NM_GT0 = (NM_GT0, 'WT', 'wild type')
T_NM_GT1 = (NM_GT1, 'PGM', 'PGM mutant')
T_NM_GT2 = (NM_GT2, 'SIRK1', 'SIRK1 mutant')
T_NM_GT3 = (NM_GT3, 'SIRK1-PGM', 'PGM/SIRK1 double mutant')
T_NM_GT4 = (NM_GT4, 'SIRK1-SWEET', 'SIRK1/SWEET double mutant')
T_NM_GT5 = (NM_GT5, 'SWEET', 'SWEET mutant')

T_NM_FT1 = (NM_FT1, 'DR', 'day root')
T_NM_FT2 = (NM_FT2, 'DS', 'day shoot')
T_NM_FT3 = (NM_FT3, 'NR', 'night root')
T_NM_FT4 = (NM_FT4, 'NS', 'night shoot')

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
S_PROT = 'Protein'
S_BIN_C_G = 'BinCode'
S_BIN_C_3 = 'BinCode3'
S_BIN_C_2 = 'BinCode2'
S_BIN_C_1 = 'BinCode1'
L_S_PHO_CL = [S_PROT, S_BIN_C_G, S_BIN_C_3, S_BIN_C_2, S_BIN_C_1]

S_MN_CONC = 'MeanConc'

S_NEG = 'Neg'
S_POS = 'Pos'
S_AV = 'Av'
S_TOP = 'Top'
S_SUM = 'Sum'
S_IDX = 'Idx'
S_CORR_S = 'Corr'
S_CORR_L = 'Correlation'
S_CORR_V = 'CorrV'
S_CORR_P = 'CorrP'
S_SPEAR_V = 'SpearV'
S_SPEAR_P = 'SpearP'
S_KEND_V = 'KendV'
S_KEND_P = 'KendP'
S_DV_SC = 'DvSc'
S_DV_CL_SD = 'DvClSD'
S_OCC_CL_SD = 'OccClSD'
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
S_Y_P_VAL = '$-\log_{10}$(p-value)'

# --- other constants ---------------------------------------------------------
M_DETER = 'Deterministic'
M_STOCH = 'Stochastic'

R04 = 4

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
S_SNDS = S_SUM + S_NEG + S_DV_SC
S_SPDS = S_SUM + S_POS + S_DV_SC
S_INDC = S_IDX + S_NEG + S_DV_CL_SD
S_IPDC = S_IDX + S_POS + S_DV_CL_SD

###############################################################################

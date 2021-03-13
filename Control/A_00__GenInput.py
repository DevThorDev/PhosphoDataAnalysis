# -*- coding: utf-8 -*-
###############################################################################
# --- A_00__GenInput.py -------------------------------------------------------
###############################################################################
import Core.C_00__GenConstants as GC

# --- Input: Flow control -----------------------------------------------------
calcCorrMetPho = False   # correlations between MetData and PhoData
calcCorrWithin = False  # correlations within MetData, PhoData and CombData (!)
calcKMeansClBasicDt = True      # K-Means cluster of basic data (MetData,...)
calcCorrClDt = True             # correl. btw. clusters of Met/Pho/CombData
calcClustersCorrelX = True      # clustering, correlations of clusters (DataX)

doCombDt = False                # do calculations with CombData

# --- Input: Output and debug info --------------------------------------------
nStepsProgOut = 10000   # number of steps when progress output is displayed
levelDebugOut = 2       # level of debug output (0: no debug output)

# --- Input: Transformations and settings investigated ------------------------
nmX = GC.NM_X
nmXFtC = GC.NM_XFEAT_C
nmXFtG = GC.NM_XFEAT_G
nmTpXS = GC.NM_TYPEX_S
nmTpXF = GC.NM_TYPEX_F

dInvXFtTpX = {nmXFtC: [nmTpXS]}
# dInvXFtTpX = {nmXFtC: [nmTpXS, nmTpXF]}
# dInvXFtTpX = {nmXFtG: [nmTpXS, nmTpXF]}
# dInvXFtTpX = {nmXFtC: [nmTpXS], nmXFtG: [nmTpXS]}
# dInvXFtTpX = {nmXFtC: [nmTpXS, nmTpXF], nmXFtG: [nmTpXS, nmTpXF]}
    # keys: 'C'(classic features) / 'G'(genotypes as features)
    # values: 'S': num. of samples multiplied by num. of genotypes
    #         'F': num. of features multiplied by num. of genotypes

sTransf = '2L'                              # '2L'(2log) / 'Un'(no transf.)
lIInp1 = [11, 21, 31, 41, 51, 81, 82, 83]   # redef. dep. on sTransf (below)

# --- Input: Names of genotypes and features ----------------------------------
dIDGT = {GC.T_NM_GT0[0]: GC.T_NM_GT0[1],
         GC.T_NM_GT1[0]: GC.T_NM_GT1[1],
         GC.T_NM_GT2[0]: GC.T_NM_GT2[1],
         GC.T_NM_GT3[0]: GC.T_NM_GT3[1],
         GC.T_NM_GT4[0]: GC.T_NM_GT4[1],
         GC.T_NM_GT5[0]: GC.T_NM_GT5[1]}
dIDFt = {GC.T_NM_FT1[0]: GC.T_NM_FT1[1],
         GC.T_NM_FT2[0]: GC.T_NM_FT2[1],
         GC.T_NM_FT3[0]: GC.T_NM_FT3[1],
         GC.T_NM_FT4[0]: GC.T_NM_FT4[1]}
lIDGT = list(dIDGT)
lIDFt = list(dIDFt)
lNmGT = list(dIDGT.values())
lNmFt = list(dIDFt.values())
dNmGT = {cV: cK for cK, cV in dIDGT.items()}
dNmFt = {cV: cK for cK, cV in dIDFt.items()}

# lInvGT = [lIDGT[0]]                 # first genotype only
# lInvGT = [lIDGT[0], lIDGT[-1]]      # first and last genotype
lInvGT = [lIDGT[0], lIDGT[1], lIDGT[-1]]    # first, second and last genotype
# lInvGT = lIDGT                      # all genotypes

lDescGT = [GC.T_NM_GT0[2], GC.T_NM_GT1[2], GC.T_NM_GT2[2], GC.T_NM_GT3[2],
           GC.T_NM_GT4[2], GC.T_NM_GT5[2]]
lDescFt = [GC.T_NM_FT1[2], GC.T_NM_FT2[2], GC.T_NM_FT3[2], GC.T_NM_FT4[2]]
nmAllGT = GC.NM_ALL_GT
nmAllFt = GC.NM_ALL_FT
nmAllGTFt = GC.NM_ALL_GT_FT

# --- Input: Filter dictionary ------------------------------------------------
# iInp1: (query (or '': no filter), name end (or '': standard name end))
dQryNmFE = {11: ('', ''),   # 11: ('Name_Analyte == "Sucrose__D-_(8TMS)"', 'Sucrose')
            12: ('', ''),
            21: ('', ''),
            22: ('', '')}
dNmPltClCent = {11: [],     # 11: ['Sucrose__D-_(8TMS)'],
                21: [],
                31: [],
                41: []}

# --- Input: General ----------------------------------------------------------
tIDResFig = GC.T_ID_RES_FIG                     # ('R'(esults), 'F'(igures))
nDigObj = 2         # number of digits reserved for all input objects
cMode = GC.M_STOCH  # GC.M_DETER / GC.M_STOCH

# --- Input: Names, indices and file name extensions --------------------------
sBase = GC.S_BASE
sExpD = GC.S_EXP_D
sMetD = GC.S_MET_D
sPhoD = GC.S_PHO_D
sCombD = GC.S_COMB_D
sExpDX = GC.S_EXP_DX
sMetDX = GC.S_MET_DX
sPhoDX = GC.S_PHO_DX
sClRD = GC.S_CLR_D
sBinOp = GC.S_BIN_OP
sClust = GC.S_CLUST
sKMClust = GC.S_KM_CLUST
sAgClust = GC.S_AG_CLUST
sOvRep = GC.S_OV_REP
lSIDs = [sExpD, sMetD, sPhoD, sCombD, sExpDX, sClRD, sBinOp, sClust, sOvRep]

sBaseL = GC.S_BASE_L
sExpDL = GC.S_EXP_D_L
sMetDL = GC.S_MET_D_L
sPhoDL = GC.S_PHO_D_L
sCombDL = GC.S_COMB_D_L
sExpDLX = GC.S_EXP_DX_L
sMetDLX = GC.S_MET_DX_L
sPhoDLX = GC.S_PHO_DX_L
sClRDL = GC.S_CLR_D_L
sBinOpL = GC.S_BIN_OP_L
sClustL = GC.S_CLUST_L
sOvRepL = GC.S_OV_REP_L
lSClass = [sExpDL, sMetDL, sPhoDL, sCombDL, sExpDLX, sClRDL, sBinOpL, sClustL,
           sOvRepL]

iExpD = GC.I_EXP_D_2L
iMetD = GC.I_MET_D_2L
iPhoD = GC.I_PHO_D_2L
if sTransf == 'Un':
    iMetD = GC.I_MET_D_UN
    iPhoD = GC.I_PHO_D_UN
iCombD = GC.I_COMB_D
iExpDX = GC.I_EXP_DX
iClRD = GC.I_CLR_D
iBinOp = GC.I_BIN_OP
iClust = GC.I_CLUST
iOvRep = GC.I_OV_REP
lIClass = [iExpD, iMetD, iPhoD, iCombD, iExpDX, iClRD, iBinOp, iClust, iOvRep]

lSBasicDLC = GC.L_S_BASIC_D_L_C
lSBasicDL = GC.L_S_BASIC_D_L
lSUsedDLC = GC.L_S_USED_D_L_C
lSUsedDL = GC.L_S_USED_D_L
lSCombDLC = GC.L_S_COMB_D_L_C
lSCombDL = GC.L_S_COMB_D_L
lSMergDLC = GC.L_S_MERGED_D_L_C
lSMergDL = GC.L_S_MERGED_D_L
lSClRDLC = GC.L_S_CLR_D_L_C
lSClRDL = GC.L_S_CLR_D_L
lSCalcDL = GC.L_S_CALC_D_L

nmEStd = GC.NM_E_STD
nmIFtrBs = GC.NM_I_FTR_BASE
nmTrD = GC.NM_TRANS_D
nmStdOp = GC.NM_STD_OP
nmDevTp = GC.NM_DEV_TP
nmClRes = GC.NM_CL_RES
nmClust = GC.NM_CLUSTERS
nmMeans = GC.NM_MEANS
nmSDs = GC.NM_SDS
idOAtDfr = GC.ID_O_ATTRDFR
idGTAtDfr = GC.ID_GT_ATTRDFR
nmGTAtDfr = GC.NM_GT_ATTRDFR
idFtAtDfr = GC.ID_FT_ATTRDFR
nmFtAtDfr = GC.NM_FT_ATTRDFR
nmWtAtDfr = GC.NM_WT_ATTRDFR
nmPreTr = GC.NM_PRE_TRANS
nmPreDv = GC.NM_PRE_DEV

nmExtCSV = GC.NM_EXT_CSV
nmExtPY = GC.NM_EXT_PY
nmExtPDF = GC.NM_EXT_PDF

# --- Input: Other strings ----------------------------------------------------
sProt = GC.S_PROT
sBinCG = GC.S_BIN_C_G
sBinC3 = GC.S_BIN_C_3
sBinC2 = GC.S_BIN_C_2
sBinC1 = GC.S_BIN_C_1
lSPhoCl = GC.L_S_PHO_CL
sMnConc = GC.S_MN_CONC

# --- Input: Info for result files --------------------------------------------
nmID = GC.NM_ID
nmTrDt = GC.NM_TRANS_DATA
nmClResData = GC.S_CLR_D_L
nmTypeX = GC.NM_TYPEX
nmIDGTFt = GC.NM_ID_GT_FT
nmNumClust = GC.NM_NUM_CLUST
nmFilter = GC.NM_FILTER
dFNmComp = {nmID: '',           # 0 (name of ID)
            nmTrDt: '',         # 1 (name of the data transformation)
            nmClResData: '',    # 2 (name of data resulting from clustering)
            nmTypeX: '',        # 3 (name of type describing X-combination)
            nmIDGTFt: '',       # 4 (name of ID of the genotype or feature)
            nmNumClust: '',     # 5 (name of the number of clusters (unused))
            nmFilter: ''}       # 6 (name of the data filter (or 'AllD'))

# --- Input: Colours ----------------------------------------------------------
lClr4Cat = [(0.8, 0.2, 0.2), (0.2, 0.8, 0.2), (0.3, 0., 0.), (0., 0.3, 0.)]
lClrMax = [(0.9, 0., 0.), (0.8, 0.4, 0.), (0.6, 0.6, 0.), (0., 0.8, 0.),
           (0., 0.4, 0.4), (0., 0., 0.9), (0.5, 0., 0.5), (0.25, 0.25, 0.25),
           (0.5, 0.5, 0.5), (0.75, 0.75, 0.75), (1., 0.4, 0.4), (1., 0.7, 0.),
           (1., 1., 0.), (0.3, 1., 0.3), (0., 1., 1.), (0.5, 0.5, 1.),
           (1., 0., 1.)]
lClrRGB4 = [(0.8, 0., 0.), (1., 0.1, 0.1), (1., 0.4, 0.4), (1., 0.7, 0.7),
            (0., 0.8, 0.), (0.1, 1., 0.1), (0.4, 1., 0.4), (0.7, 1., 0.7),
            (0., 0., 0.8), (0.1, 0.1, 1.), (0.4, 0.4, 1.), (0.7, 0.7, 1.),
            (0.1, 0.1, 0.1), (0.3, 0.3, 0.3), (0.5, 0.5, 0.5), (0.7, 0.7, 0.7),
            (0.9, 0.9, 0.9)]
lClrSpec = [(0.9, 0., 0.), (0.8, 0.2, 0.), (0.7, 0.4, 0.), (0.6, 0.6, 0.),
            (0.3, 0.7, 0.), (0., 0.8, 0.), (0., 0.6, 0.3), (0., 0.4, 0.6),
            (0., 0.2, 0.9), (0., 0., 1.), (0.3, 0., 0.7), (0.6, 0., 0.4)]

# --- Input: Derived values and assertions ------------------------------------
doXFtC, doXFtG = False, False
if nmXFtC in dInvXFtTpX:
    doXFtC = True
if nmXFtG in dInvXFtTpX:
    doXFtG = True
lSBasDLTf = [s + sTransf for s in lSBasicDL]
lSUsedDLCTf = [s + sTransf for s in lSUsedDLC]
lSUsedDLTf = [s + sTransf for s in lSUsedDL]
lSBasDLTfCombDL = lSBasDLTf + lSCombDL
lSUsedDLTfCombDL = lSUsedDLTf + lSCombDL
lSInpDLC1 = lSUsedDLCTf + lSMergDLC + lSClRDLC + lSCalcDL
if sTransf == '2L':
    lIInp1 = [11, 21, 31, 41, 51, 81, 82, 83]
    tNmPreDirOut = GC.T_NM_RES_FIG_2L
elif sTransf == 'Un':
    lIInp1[:2] = [12, 22]
    tNmPreDirOut = GC.T_NM_RES_FIG_UN
assert len(lSClass) == len(lIClass)
assert len(lSInpDLC1) == len(lIInp1)
dIInp0 = {lSClass[k]: lIClass[k] for k in range(len(lSClass))}
dIInp1 = {lSInpDLC1[k]: lIInp1[k] for k in range(len(lSInpDLC1))}
dGTFtC = {idGT: [nmFt + '_' + nmGT for nmFt in lNmFt] for (idGT, nmGT) in
          dIDGT.items()}

# names of attributes DataFrame
lNmAtDfr = [idOAtDfr, idGTAtDfr, nmGTAtDfr, idFtAtDfr, nmFtAtDfr, nmWtAtDfr]

# --- create input dictionary -------------------------------------------------
dictInpG = {# --- Input: Flow control
            'calcCorrMetPho': calcCorrMetPho,
            'calcCorrWithin': calcCorrWithin,
            'calcKMeansClBasicDt': calcKMeansClBasicDt,
            'calcCorrClDt': calcCorrClDt,
            'calcClustersCorrelX': calcClustersCorrelX,
            'doCombDt': doCombDt,
            # --- Input: Output and debug info
            'nStPr': nStepsProgOut,
            'lvlDbg': levelDebugOut,
            # --- Input: Transformations and settings investigated
            'nmX': nmX,
            'nmXFtC': nmXFtC,
            'nmXFtG': nmXFtG,
            'nmTpXS': nmTpXS,
            'nmTpXF': nmTpXF,
            'dInvXFtTpX': dInvXFtTpX,
            'sTransf': sTransf,
            'lIInp1': lIInp1,
            'tNmPreDirOut': tNmPreDirOut,
            # --- Input: Names of genotypes and features
            'dIDGT': dIDGT,
            'dIDFt': dIDFt,
            'lIDGT': lIDGT,
            'lIDFt': lIDFt,
            'lNmGT': lNmGT,
            'lNmFt': lNmFt,
            'dNmGT': dNmGT,
            'dNmFt': dNmFt,
            'lInvGT': lInvGT,
            'lDescGT': lDescGT,
            'lDescFt': lDescFt,
            'nmAllGT': nmAllGT,
            'nmAllFt': nmAllFt,
            'nmAllGTFt': nmAllGTFt,
            # --- Input: Filter dictionary
            'dQryNmFE': dQryNmFE,
            'dNmPltClCent': dNmPltClCent,
            # --- Input: General
            'tIDResFig': tIDResFig,
            'nDigObj': nDigObj,
            'Mode': cMode,
            # --- Input: Names, indices and file name extensions
            'sBase': sBase,
            'sExpD': sExpD,
            'sMetD': sMetD,
            'sPhoD': sPhoD,
            'sCombD': sCombD,
            'sExpDX': sExpDX,
            'sMetDX': sMetDX,
            'sPhoDX': sPhoDX,
            'sClRD': sClRD,
            'sBinOp': sBinOp,
            'sClust': sClust,
            'sKMClust': sKMClust,
            'sAgClust': sAgClust,
            'sOvRep': sOvRep,
            'lSIDs': lSIDs,
            'sBaseL': sBaseL,
            'sExpDL': sExpDL,
            'sMetDL': sMetDL,
            'sPhoDL': sPhoDL,
            'sCombDL': sCombDL,
            'sExpDLX': sExpDLX,
            'sMetDLX': sMetDLX,
            'sPhoDLX': sPhoDLX,
            'sClRDL': sClRDL,
            'sBinOpL': sBinOpL,
            'sClustL': sClustL,
            'sOvRepL': sOvRepL,
            'lSClass': lSClass,
            'iExpD': iExpD,
            'iMetD': iMetD,
            'iPhoD': iPhoD,
            'iCombD': iCombD,
            'iExpDX': iExpDX,
            'iClRD': iClRD,
            'iBinOp': iBinOp,
            'iClust': iClust,
            'iOvRep': iOvRep,
            'lIClass': lIClass,
            'lSBasicDLC': lSBasicDLC,
            'lSBasicDL': lSBasicDL,
            'lSUsedDLC': lSUsedDLC,
            'lSUsedDL': lSUsedDL,
            'lSCombDLC': lSCombDLC,
            'lSCombDL': lSCombDL,
            'lSMergDLC': lSMergDLC,
            'lSMergDL': lSMergDL,
            'lSClRDLC': lSClRDLC,
            'lSClRDL': lSClRDL,
            'lSCalcDL': lSCalcDL,
            'nmEStd': nmEStd,
            'nmIFtrBs': nmIFtrBs,
            'nmTrD': nmTrD,
            'nmStdOp': nmStdOp,
            'nmDevTp': nmDevTp,
            'nmClRes': nmClRes,
            'nmClust': nmClust,
            'nmMeans': nmMeans,
            'nmSDs': nmSDs,
            'idOAtDfr': idOAtDfr,
            'idGTAtDfr': idGTAtDfr,
            'nmGTAtDfr': nmGTAtDfr,
            'idFtAtDfr': idFtAtDfr,
            'nmFtAtDfr': nmFtAtDfr,
            'nmWtAtDfr': nmWtAtDfr,
            'nmPreTr': nmPreTr,
            'nmPreDv': nmPreDv,
            'nmExtCSV': nmExtCSV,
            'nmExtPY': nmExtPY,
            'nmExtPDF': nmExtPDF,
            # --- Input: Other strings
            'sProt': sProt,
            'sBinCG': sBinCG,
            'sBinC3': sBinC3,
            'sBinC2': sBinC2,
            'sBinC1': sBinC1,
            'lSPhoCl': lSPhoCl,
            'sMnConc': sMnConc,
            # --- Input: Info for result files
            'nmID': nmID,
            'nmTrDt': nmTrDt,
            'nmClResData': nmClResData,
            'nmTypeX': nmTypeX,
            'nmIDGTFt': nmIDGTFt,
            'nmNumClust': nmNumClust,
            'nmFilter': nmFilter,
            'dFNmComp': dFNmComp,
            # --- Input: Colours
            'lClr4Cat': lClr4Cat,
            'lClrMax': lClrMax,
            'lClrRGB4': lClrRGB4,
            'lClrSpec': lClrSpec,
            # --- Input: Derived values and assertions
            'doXFtC': doXFtC,
            'doXFtG': doXFtG,
            'lSBasDLTf': lSBasDLTf,
            'lSUsedDLCTf': lSUsedDLCTf,
            'lSUsedDLTf': lSUsedDLTf,
            'lSBasDLTfCombDL': lSBasDLTfCombDL,
            'lSUsedDLTfCombDL': lSUsedDLTfCombDL,
            'lSInpDLC1': lSInpDLC1,
            'dIInp0': dIInp0,
            'dIInp1': dIInp1,
            'dGTFtC': dGTFtC,
            'lNmAtDfr': lNmAtDfr,
            'lKModD': [nmIFtrBs, nmTrD, nmStdOp, nmDevTp]}

###############################################################################

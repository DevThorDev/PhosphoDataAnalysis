# -*- coding: utf-8 -*-
###############################################################################
# --- D_82__Clustering.py ------------------------------------------------------
###############################################################################
import Core.C_00__GenConstants as GC

# --- general -----------------------------------------------------------------
sNmSpec = GC.S_CLUST_L

# === all clustering methods ==================================================
saveClInfo = True           # save the clustering info output

# === K-Means clustering specific =============================================
initMeth = 'k-means++'      # method for initialization
nInit = 100                 # number of runs with different centroid seeds
maxIter = 300               # maximum number of iterations for a single run
relTol = 1.E-4              # rel. tolerance w.r.t. inertia for convergence
useWtsCombD = True          # use weights for combined data?

# === Agglomerative clustering specific =======================================
affinMet = 'euclidean'
linkageCrit = 'ward'

# === other input =============================================================
# --- results -----------------------------------------------------------------
cSep = ';'    # separator for results files
tCNmRes = ('Ct', 'ClusterLabel', 'ClusterCentre')

# --- names and paths of files and dirs ---------------------------------------
dISPr = {'No': {0: '', 2: '_'}, 'Tr': {0: '', 2: '_'}, 'Dv': {0: '', 2: '_'}}

# --- graphics parameters 1DS_1 -----------------------------------------------
nmPlt_1DS_1 = '1DS_1'   # name of the plot and key of the input dict
title_1DS_1 = 'Inertia'     # title of 1D-simple plot 1
xLbl_1DS_1 = 'Number of clusters'   # x-label of 1D-simple plot 1
yLbl_1DS_1 = 'Inertia'  # y-label of 1D-simple plot 1
tpMark_1DS_1 = ''       # marker type of 1D-simple plot 1
szMark_1DS_1 = 1        # marker size of 1D-simple plot 1
ewMark_1DS_1 = 1        # marker edge width of 1D-simple plot 1
ecMark_1DS_1 = (1., 0., 0.)     # marker edge colour of 1D-simple plot 1
fcMark_1DS_1 = (1., 0.5, 0.)    # marker face colour of 1D-simple plot 1
styLn_1DS_1 = 'solid'   # line style of 1D-simple plot 1
wdthLn_1DS_1 = 1        # line width of 1D-simple plot 1
colLn_1DS_1 = (1., 0.5, 0.)     # line colour of 1D-simple plot 1
# --- graphics parameters 1DS_2 -----------------------------------------------
nmPlt_1DS_2 = '1DS_2'   # name of the plot and key of the input dict
title_1DS_2 = 'NIter'       # title of 1D-simple plot 2
xLbl_1DS_2 = 'Number of clusters'   # x-label of 1D-simple plot 2
yLbl_1DS_2 = 'Number of iterations' # y-label of 1D-simple plot 2
tpMark_1DS_2 = 'x'      # marker type of 1D-simple plot 2
szMark_1DS_2 = 4        # marker size of 1D-simple plot 2
ewMark_1DS_2 = 1.5      # marker edge width of 1D-simple plot 2
ecMark_1DS_2 = (1., 0., 0.)     # marker edge colour of 1D-simple plot 2
fcMark_1DS_2 = (1., 0.5, 0.)    # marker face colour of 1D-simple plot 2
styLn_1DS_2 = ''        # line style of 1D-simple plot 2
wdthLn_1DS_2 = 1        # line width of 1D-simple plot 2
colLn_1DS_2 = (0.75, 0.1, 0.75) # line colour of 1D-simple plot 2
# --- graphics parameters ClCent ----------------------------------------------
nmPlt_ClCent = 'ClCent' # name of the plot and key of the input dict
title_ClCent = 'Cluster centres'    # title of plot
xLbl_ClCent = 'Coordinates of clusters'             # x-label of plot
yLbl_ClCent =  '. change to ref. (standardised)'    # y-label of plot
xLimB = None            # bottom x-limit of plot (None / number)
xLimT = None            # top x-limit of plot (None / number)
yLimB = None              # bottom y-limit of plot (None / number (e.g. 16))
yLimT = None            # top y-limit of plot (None / number)
tpMark_ClCent = 'o'     # marker type of plot
szMark_ClCent = 4       # marker size of plot
ewMark_ClCent = 1       # marker edge width of plot
ecMark_ClCent = (0., 0., 0.95)  # marker edge colour of plot (overwritten)
fcMark_ClCent = (0.3, 0.3, 0.9) # marker face colour of plot (overwritten)
styLn_ClCent = 'solid'  # line style of plot
wdthLn_ClCent = 1       # line width of plot
colLns_ClCent = 'RG'    # line colour scheme of plot ('RG'/'RB'/'GB')
# if the list of colours lClr_ClCent is not None, it overwrites colLns_ClCent
lClr_ClCent = [(0.9, 0., 0.), (0.8, 0.4, 0.), (0.6, 0.6, 0.),
               (0., 0.8, 0.), (0., 0.4, 0.4), (0., 0., 0.9),
               (0.5, 0., 0.5), (0.25, 0.25, 0.25), (0.5, 0.5, 0.5),
               (0.75, 0.75, 0.75), (1., 0.4, 0.4), (1., 0.7, 0.),
               (1., 1., 0.), (0.3, 1., 0.3), (0., 1., 1.),
               (0.5, 0.5, 1.), (1., 0., 1.)]
coordAnchorBox = (1.1, 0.5)         # coordinates of the legend anchor box

# --- derived values ----------------------------------------------------------
nmCt, nmLbl, nmCtr = tCNmRes[0], tCNmRes[1], tCNmRes[2]

# --- assertions --------------------------------------------------------------

# --- create input dictionary -------------------------------------------------
dIO = {# --- general
       'sNmSpec': sNmSpec,
       # === all clustering methods
       'saveClInfo': saveClInfo,
       # === K-Means clustering specific
       'initMeth': initMeth,
       'nInit': nInit,
       'maxIter': maxIter,
       'relTol': relTol,
       'useWtsCombD': useWtsCombD,
       # === Agglomerative clustering specific
       'affinMet': affinMet,
       'linkageCrit': linkageCrit,
       # === other input
       # --- results
       'cSep': cSep,
       'tCNmRes': tCNmRes,
       # --- names and paths of files and dirs
       'dISPr': dISPr,
       # --- graphics parameters 1DS_1
       'nmPlt_1DS_1': nmPlt_1DS_1,
       nmPlt_1DS_1: {'title': title_1DS_1,
                     'xLbl': xLbl_1DS_1,
                     'yLbl': yLbl_1DS_1,
                     'tpMark': tpMark_1DS_1,
                     'szMark': szMark_1DS_1,
                     'ewMark': ewMark_1DS_1,
                     'ecMark': ecMark_1DS_1,
                     'fcMark': fcMark_1DS_1,
                     'styLn': styLn_1DS_1,
                     'wdthLn': wdthLn_1DS_1,
                     'colLn': colLn_1DS_1},
       # --- graphics parameters 1DS_2
       'nmPlt_1DS_2': nmPlt_1DS_2,
       nmPlt_1DS_2: {'title': title_1DS_2,
                     'xLbl': xLbl_1DS_2,
                     'yLbl': yLbl_1DS_2,
                     'tpMark': tpMark_1DS_2,
                     'szMark': szMark_1DS_2,
                     'ewMark': ewMark_1DS_2,
                     'ecMark': ecMark_1DS_2,
                     'fcMark': fcMark_1DS_2,
                     'styLn': styLn_1DS_2,
                     'wdthLn': wdthLn_1DS_2,
                     'colLn': colLn_1DS_2},
       # --- graphics parameters ClCent
       'nmPlt_ClCent': nmPlt_ClCent,
       nmPlt_ClCent: {'title': title_ClCent,
                      'xLbl': xLbl_ClCent,
                      'yLbl': yLbl_ClCent,
                      'xLimB': xLimB,
                      'xLimT': xLimT,
                      'yLimB': yLimB,
                      'yLimT': yLimT,
                      'tpMark': tpMark_ClCent,
                      'szMark': szMark_ClCent,
                      'ewMark': ewMark_ClCent,
                      'ecMark': ecMark_ClCent,
                      'fcMark': fcMark_ClCent,
                      'styLn': styLn_ClCent,
                      'wdthLn': wdthLn_ClCent,
                      'colLns': colLns_ClCent,
                      'lClr': lClr_ClCent,
                      'coordAnchorBox': coordAnchorBox},
       # --- derived values
       'nmCt': nmCt,
       'nmLbl': nmLbl,
       'nmCtr': nmCtr}

###############################################################################

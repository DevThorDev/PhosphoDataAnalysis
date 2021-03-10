# -*- coding: utf-8 -*-
###############################################################################
# --- RenameFiles.py ----------------------------------------------------------
###############################################################################
import os

# --- INPUT -------------------------------------------------------------------
pMD = os.path.join('..', '..', '..', '12_SysBio02_DataAnalysis', '92_Networkx')
lPDSub = ['01_Data', '09_Plots']

lSID = ['SelBin2sC22', 'SelBin2sD23', 'SelBin2sE24']
lPosID = [0, 2]
lSRepl = ['C21', 'D22', 'E23']

# --- ASSERTIONS --------------------------------------------------------------
assert len(lPDSub) == len(lPosID)
assert len(lSID) == len(lSRepl)

# --- FUNCTIONS ---------------------------------------------------------------
def getNewFileNameS1(sF, sId, cPosId, sNew = '', sSep = '_', sDot = '.'):
    assert len(sF.split(sDot)) == 2    # only works for file names without dot
    [sNmF, sXtF] = sF.split(sDot)
    lSSpl, sFN = sNmF.split(sSep), ''
    if lSSpl[cPosId].startswith(sId):
        for cSSpl in lSSpl[:-1]:
            cS = cSSpl
            if cSSpl == lSSpl[cPosId]:
                if len(cSSpl) > len(sNew):
                    cS = cSSpl[:-len(sNew)] + sNew
                else:
                    cS = sNew
            sFN += cS + sSep
        sFN += lSSpl[-1] + sDot + sXtF
    return sFN

def renameFilesInDirs(pMainDir, lPSubDir, sId, lPosId, sRpl = '', sSep = '_'):
    for k, cSubDir in enumerate(lPSubDir):
        pD = os.path.join(pMainDir, cSubDir)
        for sF in os.listdir(pD):
            pF = os.path.join(pD, sF)
            sFN = getNewFileNameS1(sF, sId, lPosId[k], sRpl, sSep = sSep)
            if len(sFN) > 0 and sFN != sF:
                os.rename(pF, os.path.join(pD, sFN))
                print('Renamed', sF, 'to', sFN, 'in directory', pD + '.')

# --- MAIN --------------------------------------------------------------------
for k, sID in enumerate(lSID):
    renameFilesInDirs(pMD, lPDSub, sID, lPosID, lSRepl[k])
    print('-'*24, 'Replacement', k + 1, 'by', lSRepl[k], 'DONE.', '-'*24)
print('+'*30, 'DONE.', '+'*30)

###############################################################################

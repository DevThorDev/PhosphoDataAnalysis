# -*- coding: utf-8 -*-
###############################################################################
# --- O_21__PhoData.py --------------------------------------------------------
###############################################################################
import Core.C_00__GenConstants as GC
from Core.O_01__ExpData import ExpData, ExpDataC, ExpDataG

class PhoData(ExpData):
    def __init__(self, inpDat, iTp, cXFrt = '', cGT = '', cFt = '',
                 nmFAdd = '', cQuery = ''):
        super().__init__(inpDat, iTp, cXFrt, cGT, cFt, nmFAdd, cQuery)
        self.idO = GC.S_PHO_D
        self.descO = 'Phosphopeptide data'
        self.updateDOInBasicInfo()
        self.updateSFName()
        print('Initiated "PhoData" base object.')

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class PhoDataC(ExpDataC):
    def __init__(self, inpDat, iTp, cGT = '', cFt = '', nmFAdd = ''):
        super().__init__(inpDat, iTp, cGT = cGT, cFt = cFt, nmFAdd = nmFAdd)
        self.idO = GC.S_PHO_D + 'C'
        self.descO = 'Phosphopeptide data (feature type "C")'
        self.updateDOInBasicInfo()
        self.updateSFName()
        print('Initiated "PhoDataC" base object.')

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class PhoDataG(ExpDataG):
    def __init__(self, inpDat, iTp, cGT = '', cFt = '', nmFAdd = ''):
        super().__init__(inpDat, iTp, cGT = cGT, cFt = cFt, nmFAdd = nmFAdd)
        self.idO = GC.S_PHO_D + 'G'
        self.descO = 'Phosphopeptide data (feature type "G")'
        self.updateDOInBasicInfo()
        self.updateSFName()
        print('Initiated "PhoDataG" base object.')

###############################################################################

# -*- coding: utf-8 -*-
###############################################################################
# --- O_11__MetData.py --------------------------------------------------------
###############################################################################
import Core.C_00__GenConstants as GC
from Core.O_01__ExpData import ExpData, ExpDataC, ExpDataG

class MetData(ExpData):
    def __init__(self, inpDat, iTp, cXFrt = '', cGT = '', cFt = '',
                 nmFAdd = '', cQuery = ''):
        super().__init__(inpDat, iTp, cXFrt, cGT, cFt, nmFAdd, cQuery)
        self.idO = GC.S_MET_D
        self.descO = 'Metabolite data'
        self.updateDOInBasicInfo()
        self.updateSFName()
        print('Initiated "MetData" base object.')

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class MetDataC(ExpDataC):
    def __init__(self, inpDat, iTp, cGT = '', cFt = '', nmFAdd = ''):
        super().__init__(inpDat, iTp, cGT = cGT, cFt = cFt, nmFAdd = nmFAdd)
        self.idO = GC.S_MET_D + 'C'
        self.descO = 'Metabolite data (feature type "C")'
        self.updateDOInBasicInfo()
        self.updateSFName()
        print('Initiated "MetDataC" base object.')

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class MetDataG(ExpDataG):
    def __init__(self, inpDat, iTp, cGT = '', cFt = '', nmFAdd = ''):
        super().__init__(inpDat, iTp, cGT = cGT, cFt = cFt, nmFAdd = nmFAdd)
        self.idO = GC.S_MET_D + 'G'
        self.descO = 'Metabolite data (feature type "G")'
        self.updateDOInBasicInfo()
        self.updateSFName()
        print('Initiated "MetDataG" base object.')

###############################################################################

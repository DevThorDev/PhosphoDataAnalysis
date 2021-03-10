# -*- coding: utf-8 -*-
###############################################################################
# --- O_99__TestData.py -------------------------------------------------------
###############################################################################
from Core.O_01__ExpData import ExpData
#import Core.F_00__GenFunctions as GF

class TestData(ExpData):
    def __init__(self, inpDat, iTp, cQuery = '', tIMnMx = (-1, -1)):
        super().__init__(inpDat, iTp, cQuery, tIMnMx)
        self.idO = 'TestD'
        self.descO = 'Test data'
        print('Initiated "TestData" base object.')

###############################################################################

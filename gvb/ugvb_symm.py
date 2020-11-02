#!/usr/bin/env python 
# -*- coding: utf-8 -*- 
#
#      Copyright:   Qingchun Wang @ NJU
#      File Name:   ugvb_symm.py
#            Des:   ugvb_symm
#           Mail:   qingchun720@foxmail.com
#   Created Time:   11:25 二月-25/2019 
#


from gvb import ugvb

def func():
    '''
    func
    
    '''
    pass


class UGVB(ugvb.UGVB):
    '''

    '''
    def __init__(self, mf, np=None, method='JR', auxbasis='sto-6g',
                 max_cycle=(500, 200, 50), conv_tol=(1e-8, 1e-6, 1e-8),
                 local='Boys', pop='mulliken'):
        self._scf = mf
        self.np = np


if __name__ == "__main__":
    pass




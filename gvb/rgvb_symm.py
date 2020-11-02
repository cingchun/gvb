#!/usr/bin/env python 
# -*- coding: utf-8 -*- 
#
#      Copyright:   Qingchun Wang @ NJU
#      File Name:   rgvb_symm.py
#            Des:   rgvb_symm
#           Mail:   qingchun720@foxmail.com
#   Created Time:   11:23 二月-25/2019 
#


from gvb import rgvb


def func():
    '''
    func
    
    '''
    pass


class GVB(rgvb.GVB):
    __doc__ = rgvb.GVB.__doc__ + \
              '''
              Attributes for symmetry allowed RGVB:
                  irrep_nelec : dict
                      Specify the number of electrons for particular irrep {'ir_name':int,...}.
                      For the irreps not listed in this dict, the program will choose the
                      occupancy based on the orbital energies.
          
              Examples:
          
          
          
              '''

    def __init__(self, mf, np=None, init=None, init_level=0, method='JR', auxbasis='sto-6g',
                 max_cycle=(500, 200, 50), conv_tol=(1e-8, 1e-6, 1e-8),
                 local='Boys', pop='mulliken'):
        assert (mf.mol.symmetry)
        rgvb.GVB.__init__(self, mf, np)
        self.irrep_nelec = {}
        self._keys = self._keys.union(['irrep_nelec'])


if __name__ == "__main__":
    pass





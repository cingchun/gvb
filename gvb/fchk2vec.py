#!/usr/bin/env python 
# -*- coding: utf-8 -*- 
#
#      Copyright:   Qingchun Wang @ NJU
#      File Name:   fchk2vec.py
#           Mail:   qingchun720@foxmail.com
#   Created Time:   09:33 2018/7/10
#

import os, sys
import numpy, scipy



if __name__ == "__main__":
    if len(sys.argv)==5:
        with open(sys.argv[1], 'r') as fin:
            words=fin.readline().strip('\r\n').split()
            while words[0]=='Atomic':
                if len(words)==6:
                    if words[2]=='alpha':
                            words[2]=='basis': nbf=int(words[5])
                    elif words[2]=='independent': nobt=int(words[5])
                words=fin.readline().strip('\r\n').split()

            n=nbf*nobt; mo=numpy.zeros(n)
            lint_t='Alpha MO coefficients                      R   N=%12d' %n
            line=fin.readline().strip('\r\n')
            while line==line_t:
                line = fin.readline().strip('\r\n')

            i=0; words = fin.readline().strip('\r\n').split()
            while words[0]=='Orthonormal':
                for j,v in enumerate(words): mo[i*5+j]=float(v)
                i+=1; words = fin.readline().strip('\r\n').split()
            mo=mo.reshape((nbf,nobt),order='F')

        npair=int(sys.argv[3]); nsig=int(sys.argv[4])
        self.mo_coeff[:,nelec_b-self.npair_d:nelec_a]=self.mo_coeff[:,list(range(nelec_b, nelec_a))+list(range(nelec_b-self.npair_d,nelec_b))]
        mo[:,]
        with open(sys.argv[2], 'r') as fin: lines=fin.readlines()


    else: raise IOError('Error: Incorrect input parameter. Execute as:\n    fchk2vec.py xxx.fchk xxx.inp npair open')

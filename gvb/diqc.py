#!/usr/bin/env python 
# -*- coding: utf-8 -*- 
#
#      Authon:   Qingchun Wang @ NJU
#      E-Mail:   qingchun720@foxmail.com
#

'''
diqc: QC calcuation

    main program for gvb and bccc calculation
    
    diqc xxx.gjf
'''


import sys, time
from gvb.dump import setargv,dumpfile,cal,run,extv






if __name__ == '__main__':
    setargv(sys.argv[1:])
    from gvb import dump  # for shared variable
    
    infile=dump.infile; outfile=dump.inpre+'.out'
    fout=open(outfile, 'w'); stdout=sys.stdout; stderr=sys.stderr; sys.stdout=fout; sys.stderr=fout
    print(time.strftime('Start diqc program at %H:%M:%S %m/%d %Y', time.localtime()))
    print( ' Authon:   Qingchun Wang @ NJU    ')
    print( ' Mail :   qingchun720@foxmail.com ')
    print(F' input  file:   {infile} ')
    print(F' output file:   {outfile}')
    print( '\n')

    print( 'The input parameters are as follows:')
    print(F'basis,auxbasis,cart = {dump.mol.basis},{dump.auxbasis},{dump.mol.cart}')
    print(F'metd,np,nt = {dump.metd},{dump.np},{dump.nt}')
    print(F'sch,pop,weight = {dump.sch},{dump.pop},{dump.weight}')
    print( 'mol = ')
    print(F'{dump.mol.charge},{dump.mol.spin}')
    print(F'{dump.mol.atom}')
    print(F'\n',flush=True)

    print(' BCCC calcualtion')
    inbpy = dumpfile(dump.mol, dump.inpre, 'pyscf', dump.metd, dump.np, dump.auxbasis, dump.pop, dump.weight, dump.sch)
    cal(inbpy, '.py')
    print('\n',flush=True)

    # inbinp = dumpfile(dump.mol, dump.inpre, 'gamess', 'gvb', dump.np, dump.auxbasis, outsuf='init')
    # run(F"fchk2vec {inbinp}.fchk {inbinp}.inp -gvb {extv(inbpy+'.out', 'np= %d    npa= %d')}")
    # cal(inbinp, '.inp')
    # run(F'cp {inbinp}.fchk {inbinp}.dat.fchk')
    # run(F'dat2fchk {inbinp}.dat {inbinp}.dat.fchk')
    # print('\n',flush=True)
        
    print(time.strftime('End diqc program at %H:%M:%S %m/%d %Y', time.localtime()))
    fout.close()


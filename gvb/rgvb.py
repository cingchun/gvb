#!/usr/bin/env python 
# -*- coding: utf-8 -*- 
#
#      Author:   Qingchun Wang @ NJU
#      E-Mail:   qingchun720@foxmail.com
#

'''
rgvb: Restricted GVB calculation

    Example:
    
'''


import sys, os, subprocess
import tempfile
import numpy, math, scipy
from functools import reduce

from pyscf import lib
from pyscf import gto, scf
from pyscf.lib import logger

from gvb import dump
from gvb import stability
from gvb import newton_ah


au2debye = 2.541766; lT0 = None
in_pre2b=None; reffile=None


def sortobt(mo, order):
    i=min(order); j=max(order)
    mo[:,i:j]=mo[:,order]
def swapobt(mo, i, j):
    mo[:, [i, j]] = mo[:, [j, i]]
def sortbf(mo, order):
    i = min(order); j = max(order)
    mo[i:j,:] = mo[order,:]
def swapbf(mo, i, j):
    mo[[i, j], :] = mo[[j, i], :]

def dump_mo(_mf, _file, _reffile=None, _filetype='fchk'):
    '''
    Args :
        filetype: str
            One of com, inp, and fchk
        file: str
            Name of the output file.
        reffile: str, for filetype is fchk only.
            The reference file of the output file if filetype is fchk.

    Example :

    '''
    
    if not _reffile is None: dump.run(F'cp -f {_reffile} {_file}')
    logger.note(_mf, F'dump MO to {_file}')
    # f=f.copy()
    from gvb import py2fchk
    if _filetype == 'fchk':
        if isinstance(_mf, scf.uhf.UHF):
            mo = (_mf.mo_coeff[0].copy(), _mf.mo_coeff[1].copy())
            l = numpy.sqrt((_mf.mol.intor_symmetric('int1e_ovlp')).diagonal())
            for i in range(_mf.mo_coeff[0].shape[1]):
                mo[0][:, i] *= l
                # mo[1][:, i] *= l
            logger.debug(_mf, F'shape={_mf.mo_coeff[0].shape}')
            py2fchk.py2fchk(_file, _mf.mo_coeff[0].shape[0], _mf.mo_coeff[0].shape[1], mo[0], 'a')
            # py2fchk.py2fchk(_file, _mf.mo_coeff[1].shape[0], _mf.mo_coeff[1].shape[1], mo[1], 'b')
        else:
            mo = _mf.mo_coeff.copy()
            l = numpy.sqrt((_mf.mol.intor_symmetric('int1e_ovlp')).diagonal())
            for i in range(_mf.mo_coeff.shape[1]): mo[:, i] *= l
            logger.debug(_mf, F'shape={_mf.mo_coeff.shape}')
            py2fchk.py2fchk(_file, _mf.mo_coeff.shape[0], _mf.mo_coeff.shape[1], mo, 'a')
    elif _filetype == 'com': pass
    elif _filetype == 'inp': pass
    else: raise RuntimeError(F'The {_filetype} type doesn\'t exist')
def dump_ndarr(f, c, str, fmt='tri',stdout=sys.stdout, label_row=None, label_col=None,
               ncol=5, digits=14, start_row=1, start_col=1):
    '''
    Format print for the lower triangular part of an array

    Args:
        fmt : str, one of tri, rec, mo.
            Defualt is tri
        c : numpy.ndarray
            coefficients
        stdout : file object
            eg sys.stdout, or stdout = open('/path/to/file') or
            mol.stdout if mol is an object initialized from :class:`gto.Mole`

    Kwargs:
        label_row : list of strings
            Row labels (default is 1,2,3,4,...)
        label_col : list of strings
            Col labels (default is 1,2,3,4,...)
        ncol : int
            Number of columns in the format output (default 5)
        digits : int
            Number of digits of precision for floating point output (default 5)
        start : int
            The number to start to count the index (default 0)

    Examples:

        >>> import sys, numpy
        >>> dm = numpy.eye(3)
        >>> dump_tri(sys.stdout, dm)
                #0        #1        #2
        0       1.00000
        1       0.00000   1.00000
        2       0.00000   0.00000   1.00000
        >>> from pyscf import gto
        >>> mol = gto.M(atom='C 0 0 0')
        >>> dm = numpy.eye(mol.nao_nr())
        >>> dump_tri(sys.stdout, dm, label=mol.ao_labels(), ncol=9, digits=2)
                    #0     #1     #2     #3     #4     #5     #6     #7     #8
        0  C 1s     1.00
        0  C 2s     0.00   1.00
        0  C 3s     0.00   0.00   1.00
        0  C 2px    0.00   0.00   0.00   1.00
        0  C 2py    0.00   0.00   0.00   0.00   1.00
        0  C 2pz    0.00   0.00   0.00   0.00   0.00   1.00
        0  C 3px    0.00   0.00   0.00   0.00   0.00   0.00   1.00
        0  C 3py    0.00   0.00   0.00   0.00   0.00   0.00   0.00   1.00
        0  C 3pz    0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00   1.00
    '''
    if c.ndim == 2:
        logger.note(f, str)
        from gvb import dump_mat
        if fmt == 'tri': dump_mat.dump_tri(stdout,c,label_row,ncol,digits,start_row)
        elif fmt == 'rec': dump_mat.dump_rec(stdout,c,label_row,label_col,ncol,digits,start_row,start_col)
        elif fmt == 'mo': dump_mat.dump_mo(stdout,c,label_row,ncol,digits,start_row)
        else: raise RuntimeError(F'The {fmt} format doesn\'t exist')
    else:
        from gvb import dump_mat
        logger.note(f, str+' [1]')
        if fmt == 'tri': dump_mat.dump_tri(stdout, c[0], label_row, ncol, digits, start_row)
        elif fmt == 'rec': dump_mat.dump_rec(stdout, c[0], label_row, label_col, ncol, digits, start_row, start_col)
        elif fmt == 'mo': dump_mat.dump_mo(stdout, c[0], label_row, ncol, digits, start_row)
        else: raise RuntimeError(F'The {fmt} format doesn\'t exist')
        logger.note(f, str+' [2]')
        if fmt == 'tri': dump_mat.dump_tri(stdout, c[1], label_row, ncol, digits, start_row)
        elif fmt == 'rec': dump_mat.dump_rec(stdout, c[1], label_row, label_col, ncol, digits, start_row, start_col)
        elif fmt == 'mo': dump_mat.dump_mo(stdout, c[1], label_row, ncol, digits, start_row)
        else: raise RuntimeError(F'The {fmt} format doesn\'t exist')
def load_mo(_mf, _file, _filetype='fchk'):
    logger.note(_mf, F'load mo from {_file} file')
    if not os.path.exists(_file): raise FileNotFoundError(F'The {_file} file doesn\'t exist')

    from gvb import fchk2py
    if _filetype == 'fchk':
        if isinstance(_mf, scf.uhf.UHF):
            _mf.mo_coeff = (fchk2py.fchk2py(_file, _mf.mo_coeff[0].shape[0], _mf.mo_coeff[0].shape[1], 'a'),
                            fchk2py.fchk2py(_file, _mf.mo_coeff[0].shape[0], _mf.mo_coeff[0].shape[1], 'b'))
            l = numpy.sqrt((_mf.mol.intor_symmetric('int1e_ovlp')).diagonal())
            for i in range( _mf.mo_coeff[0].shape[1]):
                _mf.mo_coeff[0][:, i] /= l; _mf.mo_coeff[1][:, i] /= l
        else:
            _mf.mo_coeff = fchk2py.fchk2py(_file, _mf.mo_coeff.shape[0], _mf.mo_coeff.shape[1], 'a')
            l = numpy.sqrt((_mf.mol.intor_symmetric('int1e_ovlp')).diagonal())
            for i in range( _mf.mo_coeff.shape[1]): _mf.mo_coeff[:, i] /= l
    elif _filetype == 'com': raise NotImplementedError('TODO: load_mo for *.com')
    elif _filetype == 'inp': raise NotImplementedError('TODO: load_mo for *.inp')
    else: logger.error(_mf, F'The {_filetype} type doesn\'t exist')

def init_by_auto():
    pass
def init_by_file():
    pass


def lowTri0(no):
    lT0 = numpy.array([0]*no, dtype=numpy.uint32)

    for i in range(1, no):
        lT0[i] = lT0[i-1]+i

    return lT0
def dei(i, j, k, l):
    '''
    Index of Lower Triangular

    :param i:
    :param j:
    :param k:
    :param l:

    :return:
    '''

    ij = lT0[j]+i if i<=j else lT0[i]+j
    kl = lT0[l]+k if k<=l else lT0[k]+l

    return ij, kl

def delta(i,j):
    if i==j: return 1.0
    else: return 0.0
def fun(q4, q3, q2, q1, x):
    sinx = math.sin(x); cosx = math.cos(x)
    sinx2 = sinx*sinx; sinx3=sinx2*sinx; sinx4 = sinx3*sinx
    cosx2 = cosx*cosx; cosx3= cosx2*cosx
    return q4*sinx4 + q3*sinx3*cosx + q2*sinx2*cosx2 + q1*sinx*cosx3
def root(q4, q3, q2, q1):
    roots=numpy.poly1d([q3, 2*(q2-2*q4), 3*(q1-q3), -2*q2, -q1]).r
    roots = roots[abs(roots.imag) < 1.0e-12].real
    if len(roots)<1: return math.atan(0)
    else:
        roots=numpy.arctan(roots)
        funs = [fun(q4, q3, q2, q3, x) for x in roots]
        i=funs.index(min(funs))
        if funs[i]<0.: return roots[i]
        else: return math.atan(0)

def parse(mol):
    '''
    parse of mol in pyscf
    
    :param mol: gto.Mole
                mol.basis: str or dict
                mol.output: str, xxxx.out or xxx_basis_gvb.out
    
    :return: basis
    '''
    
    if mol.output is None:
        raise IOError('Error: mol.output is None\n'
                      'Set mol.output for convenience, please !')
    # words = os.path.splitext(mol.output)[0].split('_')
    # if words[-1] in dump.metds: in_pre='_'.join(words[:-2])
    # else: in_pre=os.path.splitext(mol.output)[0]
    # in_pre =  '_'.join(words[:-2]) if words[-1] in dump.metds else os.path.splitext(mol.output)[0]
    in_pre = in_pre2b[:in_pre2b.rfind('_')]
    # if isinstance(mol.basis, str): fort7='%s_%s.fort.7' %(in_pre, mol.basis)
    # else: fort7=in_pre+'_gen.fort.7'
    logger.note(mol, F'in_pre = {in_pre}')
    basiss = {}
    
    lines=dump.datfort7(mol, in_pre)[1]; nline=len(lines)
    for elem in mol._basis.keys():
        bas=''
        for i in range(nline):
            words = lines[i].split()
            if len(words)>2 and elem==words[0]:
                i=i+1
                while lines[i].strip(' \r\n'):
                    words = lines[i].split()
                    if words[0]=='L': bas+=elem+'    SP\n'
                    else: bas+=elem+'    '+words[0]+'\n'
                    i=i+1; nd=int(words[1])
                    cut=True; shell=[]
                    for j in range(nd):
                        words = lines[i+j].split()
                        if nd==1 and words[2]=='0.000000000E+00':   # dubug for g shell of Cr2
                            words[2]='0.100000000E+01'
                        shell += [words[1:]]
                        if cut and words[-1]!='0.000000000E+00': cut=False
                    for k in range(nd):
                        if cut: bas+='    '+'  '.join(shell[k][:-1])+'\n'
                        else: bas+='    '+'  '.join(shell[k])+'\n'
                    i=i+nd
                break
        else: raise RuntimeError(F' not found basis for {elem}')
        logger.note(mol, F'bas=\n{bas}\n')
        basiss[elem]=gto.basis.parse(bas)

    return basiss
def stabilize(_mf):
    logger.note(_mf, 'entry stabilize()')
    def stable(_mf, _internal=True, _external=False, _verbose=None):
        if isinstance(_mf, scf.rhf.RHF):
            if isinstance(_mf, scf.rohf.ROHF): return stability.rohf_stability(_mf, _internal, _external, _verbose)
            else: return stability.rhf_stability(_mf, _internal, _external, _verbose)
        else: return stability.uhf_stability(_mf, _internal, _external, _verbose)
    mo_i, mo_e, stable_i, stable_e = stable(_mf)
    if not stable_i:
        for i in range(10):  # stable ?
            _mf.max_cycle=50
            _mf = newton_ah.newton(_mf).run(mo_i, _mf.mo_occ)
            mo_i, mo_e, stable_i, stable_e = stable(_mf)
            if stable_i: break
        else: raise RuntimeError('Error: stabilize failure')
    return _mf
def converge(_mf):
    if not _mf.converged: _mf = newton_ah.newton(_mf).run(_mf.mo_coeff, _mf.mo_occ, max_cycle=50)
    if _mf.converged: return _mf
    else: raise RuntimeError('Error: UHF SCF not converged using the Newton technique')
# def rhf()
def hf(mol, basis, cart, metd):
    '''
    :param mol: gto.Mole
                mol.basis: dict, { x: gto.basis.parse() }
                mol.output: str, xxxx.out or xxx_basis_gvb.out
    :param basis: str or dict
    :param cart: bool
    :param metd: str, uhf or rhf
    :return:
    '''
    
    mol=gto.copy(mol).build(basis=basis, cart=cart)
    if mol.output is None:
        raise IOError('Error: mol.output is None\n'
                      'Set mol.output for convenience, please !')
    # words = os.path.splitext(mol.output)[0].split('_')
    # if words[-1] in dump.metds: in_pre='_'.join(words[:-2])
    # else: in_pre=os.path.splitext(mol.output)[0]
    # in_pre='_'.join(words[:-2]) if words[-1] in dump.metds else os.path.splitext(mol.output)[0]
    in_pre=in_pre2b[:in_pre2b.rfind('_')]
    # if isinstance(basis, str): in_pre2m='%s_%s_%s' %(in_pre, basis, metd)
    # else: in_pre2m='%s_gen_%s' %(in_pre, metd)
    in_pre2m=F'{in_pre}_{mol.basis}_{metd}' if isinstance(mol.basis, str) else F'{in_pre}_gen_{metd}'

    if not os.path.exists(in_pre2m+'.fchk'): dump.mkfchk(mol, in_pre, metd)
    mol.build(basis=parse(mol))
    if metd=='uhf': mf = scf.UHF(mol)
    else: mf = scf.RHF(mol)
    mf = scf.remove_linear_dep_(mf, threshold=1e-6)
    mf.run(max_cycle=1)
    load_mo(mf, in_pre2m+'.fchk')
    mf.kernel(mf.make_rdm1())
    mf=converge(mf)
    mf=stabilize(mf)

    return mf

def localize(mol, mo, o1, v2, pop_method='dipole'):
    v1, o2 = mol.nelec
    # from pyscf.lo.boys import dipole_integral
    # dipo = dipole_integral(mol, mo)
    # o1=0; o2=nco
    # mo=boys(mo,o1,o2,dipo)
    # o1=nco; o2=nocc
    # mo=boys(mo,o1,o2,dipo)
    # o1=mol.nelec[0]; o2=nobt_s
    # mo=boys(mo,o1,o2,dipo)
    logger.note(mol, F'Pipek-Mezey pop = {pop_method}')
    if pop_method == 'dipole':
        from pyscf.lo.boys import dipole_integral
        pop = dipole_integral(mol, mo)
    else:
        from pyscf.lo.pipek import atomic_pops
        pop=atomic_pops(mol,mo,pop_method)

    if o1>1:
        mo = pm(mo, 0, o1, pop)  # unkown bug，GVB will diverge after localization when Cr2 @ 1.6Å
    # if (o2-o1)==(v2-v1):
    #     actocc = mo[:, o1:o2].copy()
    #     mo = pm(mo, o1, o2, pop)
    #     from gvb.assoc_rot import assoc_rot
    #     mo[:, v1:v2] = assoc_rot(actocc.shape[0], o2-o1, actocc, mo[:, o1:o2], mo[:, v1:v2])
    # else:
    mo = pm(mo, o1, o2, pop)
    mo = pm(mo, v1, v2, pop)
    return mo
def boys(mo, o1, o2, dipo):
    print(F'Boys localize orbital    {o1+1:3d} -> {o2:3d}')
    sum = 0
    for i in range(o1, o2):
        for k in range(dipo.shape[0]):
            sum += dipo[k,i,i]**2
    sum *= au2debye ** 2; sum_old=sum
    for it in range(50000):
        Delta = 0
        for i in range(o1, o2-1):
            for j in range(i+1, o2):
                riix=dipo[0,i,i]; riiy=dipo[1,i,i]; riiz=dipo[2,i,i]
                rjjx=dipo[0,j,j]; rjjy=dipo[1,j,j]; rjjz=dipo[2,j,j]
                rijx=dipo[0,j,i]; rijy=dipo[1,j,i]; rijz=dipo[2,j,i]
                px=riix-rjjx; py=riiy-rjjy; pz=riiz-rjjz
                Aij=rijx*rijx+rijy*rijy+rijz*rijz-(px*px+py*py+pz*pz)/4.0
                Bij=rijx*px+rijy*py+rijz*pz
                if math.fabs(Aij)<1e-10 and math.fabs(Bij)<1e-10: continue
                p1 = math.sqrt(Aij*Aij+Bij*Bij)
                cos4theta = -Aij / p1
                sin4theta = Bij / p1
                cos2theta = math.sqrt((1.0 + cos4theta) / 2.0)
                sin2theta = math.sqrt((1.0 - cos4theta) / 2.0)
                costheta = math.sqrt((1.0 + cos2theta) / 2.0)
                sintheta = math.sqrt((1.0 - cos2theta) / 2.0)
                if sin4theta < 0.0:
                    cos2theta = -cos2theta
                    # TT = costheta
                    # costheta = sintheta
                    # sintheta = TT
                    sintheta, costheta = costheta, sintheta
                if math.fabs(costheta - 1.0) < 1.0e-13: continue
                if math.fabs(sintheta - 1.0) < 1.0e-13: continue
                p = p1 * (1.0 - cos4theta)
                Delta = Delta + p
                ## TT = eigen(i);
                ## eigen(i) = costheta*TT + sintheta*eigen(j);
                ## eigen(j) = costheta*eigen(j) - sintheta*TT;
                # for k in mo.shape[0]:
                #    TT = mo[k, i]
                #    mo[k, i] = costheta * TT + sintheta * mo[k, j]
                #    mo[k, j] = costheta * mo[k, j] - sintheta * TT
                TT = mo[:, i].copy()
                mo[:, i] = costheta * TT + sintheta * mo[:, j]
                mo[:, j] = costheta * mo[:, j] - sintheta * TT

                # Transform Dipole Integrals over MOs
                p2 = sin2theta
                p1 = p2 / 2.0
                p3 = (1.0 + cos2theta) / 2.0
                p4 = (1.0 - cos2theta) / 2.0
                p5 = cos2theta
                # X Direct
                p6 = rijx * p2
                dipo[0, i, i] = riix * p3 + rjjx * p4 + p6
                dipo[0, j, j] = riix * p4 + rjjx * p3 - p6
                dipo[0, j, i] = (rjjx - riix) * p1 + rijx * p5
                dipo[0, i, j] = dipo[0, j, i]
                for n in range(o1, o2):
                    if n == i or n == j: continue
                    p = dipo[0, i, n]
                    dipo[0, i, n] = p * costheta + dipo[0, j, n] * sintheta
                    dipo[0, j, n] = dipo[0, j, n] * costheta - p * sintheta
                    dipo[0, n, i] = dipo[0, i, n]
                    dipo[0, n, j] = dipo[0, j, n]
                # Y Direct
                p6 = rijy * p2
                dipo[1, i, i] = riiy * p3 + rjjy * p4 + p6
                dipo[1, j, j] = riiy * p4 + rjjy * p3 - p6
                dipo[1, j, i] = (rjjy - riiy) * p1 + rijy * p5
                dipo[1, i, j] = dipo[1, j, i]
                for n in range(o1, o2):
                    if n == i or n == j: continue
                    p = dipo[1, i, n]
                    dipo[1, i, n] = p * costheta + dipo[1, j, n] * sintheta
                    dipo[1, j, n] = dipo[1, j, n] * costheta - p * sintheta
                    dipo[1, n, i] = dipo[1, i, n]
                    dipo[1, n, j] = dipo[1, j, n]
                # Z Direct
                p6 = rijz * p2
                dipo[2, i, i] = riiz * p3 + rjjz * p4 + p6
                dipo[2, j, j] = riiz * p4 + rjjz * p3 - p6
                dipo[2, j, i] = (rjjz - riiz) * p1 + rijz * p5
                dipo[2, i, j] = dipo[2, j, i]
                for n in range(o1, o2):
                    if n == i or n == j: continue
                    p = dipo[2, i, n]
                    dipo[2, i, n] = p * costheta + dipo[2, j, n] * sintheta
                    dipo[2, j, n] = dipo[2, j, n] * costheta - p * sintheta
                    dipo[2, n, i] = dipo[2, i, n]
                    dipo[2, n, j] = dipo[2, j, n]
        k = o2 - o1 + 1
        sum += Delta * au2debye ** 2
        p0 = math.sqrt(2.0 * Delta / (k * (k - 1)))
        if p0 < 1.0e-8:
            deg=0
            for i in range(o1,o2):
                for k in range(dipo.shape[0]):
                    deg += dipo[k,i,i]**2
            deg=math.sqrt(deg/(o2-o1))
            print(F'Up to {it+1} cycle, the global locality has reached {deg}')
            break
    else: raise RuntimeError(F'delta = {Delta:9f}: Boys Localization Failed')
    return mo
def pm(mo, o1, o2, pop):
    print(F'localize {o1+1:3d} -> {o2:3d}')
    sum = 0
    for i in range(o1, o2):
        for k in range(pop.shape[0]):
            sum += pop[k, i, i]**2
    for it in range(500000):
        Delta = 0.0
        for i in range(o1, o2-1):
            for j in range(i + 1, o2):
                Aij = 0; Bij = 0
                for k in range(pop.shape[0]):
                    piik = pop[k, i, i]; pjjk = pop[k, j, j]; pijk = pop[k, j, i]
                    P1 = piik - pjjk
                    Aij += pijk*pijk - P1*P1/4.0
                    Bij += pijk*P1
                if math.fabs(Aij) < 1e-10 and math.fabs(Bij) < 1e-10: continue

                P1 = math.sqrt(Aij*Aij + Bij*Bij)
                cos4A = -Aij/P1; sin4A = Bij/P1
                cos2A = math.sqrt((1.0 + cos4A)/2.0); sin2A = math.sqrt((1.0 - cos4A)/2.0)
                cosA = math.sqrt((1.0 + cos2A)/2.0); sinA = math.sqrt((1.0 - cos2A)/2.0)
                if sin4A < 0.0:
                    cos2A = -cos2A
                    # TT = cosA
                    # cosA = sinA
                    # sinA = TT
                    cosA, sinA = sinA, cosA
                if math.fabs(cosA-1.0)<1e-13 or math.fabs(sinA-1.0)<1e-13: continue

                PP = P1*(1.0 - cos4A)
                Delta += PP
                # for (int k = 0; k < nbs; ++k)
                # {
                #	TT = _mo(k, i + o1);
                #	_mo(k, i + o1) = cosA * TT + sinA * _mo(k, j + o1);
                #	_mo(k, j + o1) = -sinA * TT + cosA * _mo(k, j + o1);
                # }
                TT = mo[:, i].copy()
                mo[:, i] = cosA*TT + mo[:, j]*sinA
                mo[:, j] = -sinA*TT + mo[:, j]*cosA

                P2 = sin2A
                P1 = P2/2.0
                P3 = (1.0 + cos2A)/2.0
                P4 = (1.0 - cos2A)/2.0
                P5 = cos2A
                for k in range(pop.shape[0]):
                    piik = pop[k, i, i]; pjjk = pop[k, j, j]; pijk = pop[k, j, i]
                    P6 = pijk * P2
                    pop[k, i, i] = piik*P3 + pjjk*P4 + P6
                    pop[k, j, j] = piik*P4 + pjjk*P3 - P6
                    pop[k, j, i] = (pjjk - piik)*P1 + pijk*P5
                    pop[k, i, j] = pop[k, j, i]
                    for l in range(o1,o2):
                        if l == i or l == j: continue
                        PP = pop[k, i, l]
                        pop[k, i, l] = PP*cosA + pop[k, j, l]*sinA
                        pop[k, j, l] = -PP*sinA + pop[k, j, l]*cosA
                        pop[k, l, i] = pop[k, i, l]
                        pop[k, l, j] = pop[k, j, l]
        sum += Delta
        k = o2-o1+1
        p0=math.sqrt(2.0*Delta/(k*(k - 1)))
        if p0 < 1e-8:
            deg = math.sqrt(sum/(o2 - o1))
            print(F'Up to {it+1} cycle, the global locality has reached {deg}')
            break
    else: raise RuntimeError('Pipek-Mezey Localization Failed')
    return mo

def km(weight, best=True, search='dfs'):
    '''
        w1 = numpy.array(
            [[3., 4., 6., 4., 9.],
             [6., 4., 5., 3., 8.],
             [7., 5., 3., 4., 2.],
             [6., 3., 2., 2., 5.],
             [8., 4., 5., 4., 7.]])
        w2 = numpy.array(
           [[7., 6., 4., 6., 1.],
            [4., 6., 5., 7., 2.],
            [3., 5., 7., 6., 8.],
            [4., 7., 8., 8., 5.],
            [2., 6., 5., 6., 3.]])
        km(w1)
        km(w1,False)
        km(w2)
        km(w2,False)
    '''
    inf = numpy.max(weight)+1
    if best: _weight = weight.copy()
    else: _weight = inf-weight
    nx, ny = _weight.shape; transpose = False
    if nx > ny:
        transpose=True
        _weight = _weight.T; nx, ny = ny, nx
    lx = numpy.array([0.]*nx); ly = numpy.array([0.]*ny); slack = [inf]*ny
    for i in range(nx): lx[i] = numpy.max(_weight[i,:])
    px = numpy.array([None]*nx); py = numpy.array([None]*ny)
    vx = numpy.array([False]*nx); vy = numpy.array([False]*ny)

    def dfs(x): # Depth-first search
        vx[x] = True
        for y in range(ny):
            if vy[y]: continue;
            delta = lx[x] + ly[y] - _weight[x][y]
            if abs(delta) < 1.0e-10:
                vy[y] = True
                if py[y] == None or dfs(py[y]):
                    py[y] = x; px[x] = y
                    return True
            else:
                slack[y] = min(slack[y],delta)
        else: return False
    def bfs(x):
        pass
    def anti(bl): # bl: bool list
        # print(bl)
        nbl=bl.copy(); nl=len(bl)
        for i in range(nl):
            nbl[i] = not bl[i]
        return nbl

    if search == 'dfs':
        for i in range(nx):
            slack = numpy.array([inf]*ny)
            while True:
                vx = [False]*nx; vy = [False]*ny
                if dfs(i): break
                delta = numpy.min(slack[anti(vy)])
                lx[vx] -= delta
                ly[vy] += delta
                slack[anti(vy)] -= delta
    else: pass

    if not transpose: cost = sum(weight[i, px[i]] for i in range(nx))
    else: cost = sum(weight[px[i],i] for i in range(nx))
    
    return py, px, cost

def uno(mf, threshold=1):
    '''
    uno
    
    :param mf:
    :param threshold:

    :return:
    '''
    
    logger.note(mf,'entry uno()')

    s = mf.mol.intor_symmetric('int1e_ovlp')

    # nco=dump.ncore(mf.mol); nelec_a, nelec_b = mf.mol.nelec   # nelec_a > nelec_b
    # if nco>1:
    #     occ_a = mf.mo_coeff[0][:, :nco]; occ_b = mf.mo_coeff[1][:, :nco]
    #     U, l, Vt = numpy.linalg.svd(reduce(numpy.dot, (occ_a.T, s, occ_b)))
    #     mf.mo_coeff[0][:, :nco] = numpy.dot(occ_a, U); mf.mo_coeff[1][:, :nco] = numpy.dot(occ_b, Vt.T)
    #     logger.note(mf, 'uno svd(core) s =\n%s' %l)
    nco = 0; nelec_a, nelec_b = mf.mol.nelec   # nelec_a > nelec_b

    vir_a = mf.mo_coeff[0][:,nelec_a:].copy(); vir_b = mf.mo_coeff[1][:,nelec_b:].copy()
    U, l, Vt = numpy.linalg.svd(reduce(numpy.dot,(vir_a.T, s, vir_b)))
    # mf.mo_coeff[0][:, nelec_a:] = numpy.dot(vir_a, U); mf.mo_coeff[1][:,nelec_b:] = numpy.dot(vir_b, Vt.T)
    vir_a = numpy.dot(vir_a, U); vir_b = numpy.dot(vir_b, Vt.T)
    logger.note(mf, F'uno svd(vir) s =\n{l}')

    occ_a = mf.mo_coeff[0][:,nco:nelec_a].copy(); occ_b = mf.mo_coeff[1][:,nco:nelec_b].copy()
    U, l, Vt = numpy.linalg.svd(reduce(numpy.dot,(occ_a.T, s, occ_b)))
    # mf.mo_coeff[0][:, nco:nelec_a] = numpy.dot(occ_a, U); mf.mo_coeff[1][:,nco:nelec_b] = numpy.dot(occ_b, Vt.T)
    occ_a = numpy.dot(occ_a, U); occ_b = numpy.dot(occ_b, Vt.T)
    logger.note(mf, F'uno svd(occ) s =\n{l}')
    # nactocc = nelec_b-dump.ncore(mf.mol); nactvir = dump.nrefobt(mf.mol, '_'.join(mf.mol.output.split('_')[:-2]))-nelec_a
    nactocc = nelec_b-dump.ncore(mf.mol); nactvir = dump.nrefobt(mf.mol, in_pre2b[:in_pre2b.rfind('_')])-nelec_a
    np = min(nactocc, nactvir) # max np
    if threshold is None:
        npu = numpy.count_nonzero((l-0.99999)<0)
        if npu>np: npu=np
    elif threshold>=1:
        if threshold<=np: npu = numpy.count_nonzero((l-0.99999)<0)
        else:
            raise RuntimeError(F'np({threshold}) > min(nactocc,nactvir)({np})')
    else:
        npu=numpy.count_nonzero((l-threshold)<0)
        if npu>np:
            raise RuntimeError(F'npu({npu}) > min(nactocc,nactvir)({np})')
    n=-1; a_b=nelec_a-nelec_b
    while n>-(npu+1):
        o_a = occ_a[:, n-a_b].copy(); v_a = vir_a[:,n].copy()
        o_b = occ_b[:, n].copy(); v_b = vir_b[:,n-a_b].copy()
        occ_a[:,n-a_b] = (o_a+o_b)/((2*(1+l[n]))**.5); vir_a[:,n] = (o_a-o_b)/((2*(1-l[n]))**.5)
        occ_b[:,n] = (v_a-v_b)/((2*(1-l[n]))**.5); vir_b[:,n-a_b] = (v_a+v_b)/((2*(1+l[n]))**.5)
        n-=1
    mf.mo_coeff[0][:,nco:nelec_a] = occ_a; mf.mo_coeff[0][:, nelec_a:] = vir_a[:,::-1]
    mf.mo_coeff[1][:,nco:nelec_b] = occ_b; mf.mo_coeff[1][:, nelec_b:] = vir_b[:,::-1]

    return npu

class GVB(lib.StreamObject):
    ''' GVB base class. Non-relativistic RGVB

        Attributes:

        verbose : int
            Print level.  Default value equals to :class:`Mole.verbose`
        max_memory : float or int
            Allowed memory in MB.  Default value equals to :class:`Mole.max_memory`
        chkfile : str
            checkpoint file to save MOs, orbital energies etc.

        np : None, int or str
            The number of pair in GVB. if it is None, this indicate auto detetmine napir and init guess
            and str come from file
            Default is None.
        init : None, ndarray or str
            The init guess for GVB calculation. if it is None, this indicate detetminate automatically,
            and str reading from file.
            Default is None
        init_in : int
            The levle of init from input orbitals
            The levle of init. As shown below
            0(mol) ->  1(rhf/uno) -> 2(pao) -> 3(svo) -> 4(local) -> 5(init-I) -> 6(init-II) -> 7(gvb)
            Default is 0.
        init_in : end
            The levle of init for output orbitals
            The levle of init. As shown below
            0(mol) ->  1(rhf/uno) -> 2(pao) -> 3(svo) -> 4(local) -> 5(init-I) -> 6(init-II) -> 7(gvb)
            Default is 5.
        mothod : str
            orbital optimization method. It can be one of 'JR', 'OCBSE', 'SOSCF'.
            Default is 'JR'
        auxbasis : str
            Auxiliary basis is reference basis set used for projection, e.g
            sto-3g or sto-6g (default) in general.
        max_cycle : (int, int, int), for JR method only
            method max number for (global, local, ci) optimization.
            Default is (500, 200, 50)
        conv_tol : (float, float, float), for JR method only
            converge threshold to (global, local, theta) optimization.
            Default is (1e-8, 1e-4, 1e-8)

        local : str, for np = None only
            orbital localization method. It can be one of 'Boys', 'PM', 'ER', 'NAO', 'IAO'
            Default is 'Boys'
        pop: str, for Pipek-Mezey localization only
            Population analysis method. It can be one of 'mulliken', 'lowdin', 'meta_lowdin'
            Default is 'mulliken'

    Saved results

        converged : bool
            SCF converged or not
        e_tot : float
            Total GVB energy (electronic energy plus nuclear repulsion)
        ci : ndarray
            GVB pair CI coefficients
        mo_coeff
            Orbital coefficients

    Examples:

        # start from mol object
        >>> from pyscf import gto, scf, gvb
        >>> mol = gto.M(atom='H 0 0 0; H 0 0 1')
        >>> mf = gvb.GVB(mol).run()

        or:

        # start from HF/ROHF/UHF
        >>> from pyscf import gto, scf, gvb
        >>> mol = gto.M(atom='H 0 0 0; H 0 0 1')
        >>> _mf = scf.RHF(mol).run()
        >>> mf = gvb.GVB(_mf).run()

    '''
    
    def __init__(self, mf, np=None):
        # basis parameters
        self.mol = mf
        self.verbose = mf.verbose
        self.stdout = mf.stdout
        self.max_memory = mf.max_memory

        # parameters
        self.chkfile = None
        self.max_cycle = (500, 200, 50)
        self.conv_tol = (1e-8, 1e-4, 1e-8)
        self.np = np
        self.nco = 0
        self.npa = None  # n pair all
        self.init = None
        self.init_in = 0
        self.sch = 'sch2'
        self.init_out = 5
        self.basis=None
        self.auxbasis = 'sto-6g'
        self.pop = 'meta_Lowdin'
        self.weight = 'trans_dipo'
        
        # saved results
        self.converged = False
        self.e_tot = 0.
        self.mo_coeff = None
        self.mo_occ = None
        self.ci = None

        # additional parameters under this implementation
        self.nbs = None
        self.nobt = None
        self.nocc = None
        self.nvir = None
        self.pair = None
        # active space fix
        self.nactocc = None
        self.nactvir = None
        # self.npr = None
        # self.npu = 0

        self.e_nuc = 0.
        self.e_ele = 0.
        self.h_mo = None
        self.g_mo = None
        self._keys = set(self.__dict__.keys())
 
    def get_init(self):
        if self.init is None: self.init_by_auto()
        else: self.init_by_file()
    def init_by_auto(self):
        global in_pre2b, reffile
        if isinstance(self.mol, gto.Mole):  # mf -> mol
            self._chkfile = tempfile.NamedTemporaryFile(dir=lib.param.TMPDIR)
            self.chkfile = self._chkfile.name
            self.basis=self.mol.basis
            in_pre2b = self.mol.output[:self.mol.output.rfind('_')]; reffile = in_pre2b+'_rhf.fchk'
            self.mol.build(basis=parse(self.mol))
    
            mf=hf(self.mol, self.basis, self.mol.cart, 'uhf')
            if self.sch=='sch1' or mf.spin_square()[1]-mf.mol.spin-1<0.1:  # 2S+1 # Sch-I
                mf=hf(self.mol, self.basis, self.mol.cart, 'rhf')
        else:  # mf -> scf
            mf = self.mol; self.mol = mf.mol
            self.chkfile=mf.chkfile
            self.basis=self.mol.basis
            in_pre2b = self.mol.output[:self.mol.output.rfind('_')]; reffile = in_pre2b+'_rhf.fchk'
            self.mol.build(basis=parse(self.mol))
            converge(mf)
            stabilize(mf)
        if isinstance(mf, scf.uhf.UHF):
            npu = uno(mf, threshold=self.np)
            if not os.path.exists(reffile):
                dump.run(F'fchk_uhf2rhf {in_pre2b}_uhf.fchk')
                dump.run(F'mv {in_pre2b}_uhf_r.fchk {reffile}')
            # dump_mo(mf, F'{in_pre2b}_gvb_uno.fchk', reffile)
            mf.mo_coeff=mf.mo_coeff[0]; mf.mo_occ=mf.mo_occ[0]
        else: npu = 0
        self.mo_coeff = mf.mo_coeff; self.mo_occ = mf.mo_occ
        self.e_tot = mf.e_tot; self.e_nuc = self.mol.energy_nuc(); self.e_ele = self.e_tot-self.e_nuc
        
        self.nbs,self.nobt = self.mo_coeff.shape; nelec_a,nelec_b = self.mol.nelec
        self.nocc,self.nvir = nelec_a,self.nobt-nelec_a
        nfrez,nobt_s = dump.ncore(self.mol),dump.nrefobt(self.mol, in_pre2b[:in_pre2b.rfind('_')])
        self.nactocc = nelec_b-nfrez; self.nactvir = nobt_s-nelec_a
        npr = min(self.nactocc,self.nactvir)
        if not self.np: self.np = npr
        elif self.np<1: self.np = npu
        self.nco,self.npa = nelec_b-self.np, self.nobt-self.np
        
        logger.note(self, F'nbs= {self.nbs}    nobt= {self.nobt}')
        logger.note(self, F'nactocc= {self.nactocc}    nactvir= {self.nactvir}')
        logger.note(self, F'np= {self.np}    npa= {self.npa}')
        
        self.pair = numpy.zeros((self.npa, 2), dtype=numpy.uint32)
        for i in range(self.nocc): self.pair[i, 0] = i
        for i in range(self.np): self.pair[self.nocc-1-i, 1] = self.nocc+i
        for i in range(self.nocc, self.npa): self.pair[i, 0] = i+self.np
        self.ci = numpy.zeros((self.npa, 2)); self.ci[:,0] = 1.
        
        if npu<self.np:
            self.svo()
            # dump_mo(self, F'{in_pre2b}_gvb_svo.fchk', reffile)
        
        np0,np1 = self.nocc-self.np,self.nocc+self.np
        self.mo_coeff = localize(self.mol, self.mo_coeff, np0, np1, self.pop)
        # dump_mo(self, F'{in_pre2b}_gvb_localize.fchk', reffile)
        self.pairing(np0, np1, self.weight)
        
        dump_mo(self, F'{in_pre2b}_gvb_init.fchk', reffile)
    def init_by_file(self):
        global in_pre2b, reffile
        self._chkfile = tempfile.NamedTemporaryFile(dir=lib.param.TMPDIR)
        self.chkfile = self._chkfile.name
        self.basis=self.mol.basis
        in_pre2b = self.mol.output[:self.mol.output.rfind('_')]; reffile = in_pre2b+'_rhf.fchk'
        self.mol.build(basis=parse(self.mol))
        if not os.path.exists(reffile):
            dump.run(F'fchk_uhf2rhf {in_pre2b}_uhf.fchk')
            dump.run(F'mv {in_pre2b}_uhf_r.fchk {reffile}')

        mf = scf.RHF(self.mol).run(max_cycle=1)
        load_mo(mf, self.init)
        self.mo_coeff = mf.mo_coeff; self.mo_occ = mf.mo_occ
        self.e_tot = mf.e_tot; self.e_nuc = self.mol.energy_nuc(); self.e_ele = self.e_tot-self.e_nuc
        
        self.nbs,self.nobt = self.mo_coeff.shape; nelec_a,nelec_b = self.mol.nelec
        self.nocc,self.nvir = nelec_a,self.nobt-nelec_a
        nfrez,nobt_s = dump.ncore(self.mol),dump.nrefobt(self.mol, in_pre2b[:in_pre2b.rfind('_')])
        self.nactocc,self.nactvir = nelec_b-nfrez,nobt_s-nelec_a
        if not self.np: self.np = min(self.nactocc,self.nactvir)
        if self.pair:
            vl = [o+2*(self.nocc-o)-1 for o in self.pair[::-1]]
            pair = zip(self.pair, vl)
            ol,vl = zip(*pair)
            _ol,_vl = tuple(o for o in range(self.nocc) if o not in ol),tuple(o for o in range(self.nocc, self.nobt) if o not in vl)
            self.mo_coeff[:,:] = self.mo_coeff[:,[*(_ol+ol),*(vl+_vl)]]
        self.nco,self.npa = nelec_b-self.np,self.nobt-self.np
        
        logger.note(self, F'nbs= {self.nbs}    nobt= {self.nobt}')
        logger.note(self, F'nactocc= {self.nactocc}    nactvir= {self.nactvir}')
        logger.note(self, F'np= {self.np}    npa= {self.npa}')

        self.pair = numpy.zeros((self.npa, 2), dtype=numpy.uint32)
        for i in range(self.nocc): self.pair[i, 0] = i
        for i in range(self.np): self.pair[self.nocc-1-i, 1] = self.nocc+i
        for i in range(self.nocc, self.npa): self.pair[i, 0] = i+self.np
        self.ci = numpy.zeros((self.npa, 2)); self.ci[:,0] = 1.
        
        dump_mo(self, F'{in_pre2b}_gvb_init.fchk', reffile)
    def uno(self, mf, threshold=0.98):
        logger.note(self, 'entry uno()')

        s = mf.mol.intor_symmetric('int1e_ovlp')
        nco=dump.ncore(self.mol); nelec_a, nelec_b = mf.mol.nelec  # defult nelec_a>nelec_b

        occ_a = mf.mo_coeff[0][:, :nco]; occ_b = mf.mo_coeff[1][:, :nco]
        U, l, Vt = numpy.linalg.svd(reduce(numpy.dot, (occ_a.T, s, occ_b)))
        mf.mo_coeff[0][:, :nco] = numpy.dot(occ_a, U); mf.mo_coeff[1][:, :nco] = numpy.dot(occ_b, Vt.T)
        logger.note(self, F'uno svd(core) s =\n{l}')

        vir_a = mf.mo_coeff[0][:, nelec_a:]; vir_b = mf.mo_coeff[1][:, nelec_b:]
        U, l, Vt = numpy.linalg.svd(reduce(numpy.dot, (vir_a.T, s, vir_b)))
        vir_a = numpy.dot(vir_a, U); vir_b = numpy.dot(vir_b, Vt.T)
        occ_a = mf.mo_coeff[0][:, nco:nelec_a]; occ_b = mf.mo_coeff[1][:, nco:nelec_b]
        U, l, Vt = numpy.linalg.svd(reduce(numpy.dot, (occ_a.T, s, occ_b)))
        occ_a = numpy.dot(occ_a, U); occ_b = numpy.dot(occ_b, Vt.T)
        logger.note(self, F'uno svd(occ) s =\n{l}')
        npu = numpy.count_nonzero((l - threshold) < 0)
        if npu > self.npr: npu = self.npr
        n = -1; a_b = nelec_a - nelec_b
        while n > -(npu+1):
            o_a = occ_a[:, n-a_b].copy(); v_a = vir_a[:, n].copy()
            o_b = occ_b[:, n].copy(); v_b = vir_b[:, n-a_b].copy()
            occ_a[:, n-a_b] = (o_a+o_b)/((2*(1+l[n]))**.5); vir_a[:,n] = (o_a-o_b)/((2*(1-l[n]))**.5)
            occ_b[:, n] = (v_a+v_b)/((2*(1+l[n]))**.5); vir_b[:,n-a_b] = (v_a-v_b)/((2*(1-l[n]))**.5)
            n -= 1
        mf.mo_coeff[0][:, nco:nelec_a] = occ_a; mf.mo_coeff[0][:, nelec_a:] = vir_a[:, ::-1]
        mf.mo_coeff[1][:, nco:nelec_b] = occ_b; mf.mo_coeff[1][:, nelec_b:] = vir_b[:, ::-1]
        return npu

    def pao(self):
        logger.note(self,'entry pao()')
        nobt=self.nocc
        d=numpy.dot(self.mo_coeff[:,:nobt],self.mo_coeff[:,:nobt].T)
        s=self.mol.intor_symmetric('int1e_ovlp')
        v=numpy.eye(self.nbs)-numpy.dot(d,s)

        l, u = numpy.linalg.eigh(reduce(numpy.dot,(v.T,s,v)))
        logger.note(self, F'l={l}')
        u = u[:, abs(l.real)>1.0e-6].real; l=l[abs(l.real)>1.0e-6].real
        l=1/numpy.sqrt(l)
        self.mo_coeff[:,nobt:] =reduce(numpy.dot,(v,u,numpy.diag(l)))

    def svd(self, mo1, mo2, s):
        U, l, Vt = numpy.linalg.svd(reduce(numpy.dot, (mo1.T, s, mo2)))
        mo1 = numpy.dot(mo1, U); mo2=numpy.dot(mo2, Vt.T)
        logger.note(self, F'svo SVD overlap singular:\n{l}')

    def svo(self):
        logger.note(self, 'entry svo()')
        mol = gto.copy(self.mol).build(basis=self.auxbasis, cart=False)
        mol.basis=parse(mol)
        mf=hf(mol, self.auxbasis, False, 'uhf')
        if self.sch=='sch1' or mf.spin_square()[1]-mf.mol.spin-1 < 0.1:  # 2S+1 # Sch-I
            mf=hf(mol, self.auxbasis, False, 'rhf')
        if isinstance(mf, scf.uhf.UHF):
            uno(mf)
            mf.mo_coeff = mf.mo_coeff[0]
            mf.mo_occ = mf.mo_occ[0]
        s = gto.intor_cross('int1e_ovlp', self.mol, mol)
        U, l, Vt = numpy.linalg.svd(reduce(numpy.dot, (self.mo_coeff[:, self.nocc:].T, s, mf.mo_coeff[:, self.nocc:])))
        self.mo_coeff[:, self.nocc:] = numpy.dot(self.mo_coeff[:, self.nocc:], U)
        logger.note(self, F'svo SVD overlap singular:\n{l}')

    def localize(self, o1, v2):
        o2=self.mol.nelec[1]; v1=self.mol.nelec[0]
        if o1>1: logger.note(self, F'entry localize(): |1->{o1}|,  |{o1+1}->{o2}|,  |{v1+1}->{v2}|')
        else: logger.note(self, F'entry localize(): |{o1+1}->{o2}|,  |{v1+1}->{v2}|')
        from pyscf import lo
        if self.pop in set(['meta_lowdin', 'lowdin', 'mulliken']):
            localizer = lo.PipekMezey(self.mol); localizer.pop_method = self.pop
        elif self.pop == 'dipole': localizer = lo.Boys(self.mol)
        else: localizer = lo.EdmistonRuedenberg(self.mol)
        if o1>1: self.mo_coeff[:, :o1] = localizer.kernel(self.mo_coeff[:, :o1])
        # if (o2-o1)==(v2-v1):
        #     actocc = self.mo_coeff[:, o1:o2].copy()
        #     self.mo_coeff[:, o1:o2] = localizer.kernel(self.mo_coeff[:, o1:o2])
        #     from gvb.assoc_rot import assoc_rot
        #     self.mo_coeff[:, v1:v2] = assoc_rot(actocc.shape[0], o2-o1, actocc, self.mo_coeff[:, o1:o2], self.mo_coeff[:, v1:v2])
        # else:
        self.mo_coeff[:, o1:o2] = localizer.kernel(self.mo_coeff[:, o1:o2])
        self.mo_coeff[:, v1:v2] = localizer.kernel(self.mo_coeff[:, v1:v2])
        
        from pyscf.lo.boys import dipole_integral
        dipo = dipole_integral(self.mol, self.mo_coeff)
        deg = 0
        if o2>o1:
            for i in range(o1, o2):
                deg += dipo[0, i, i]**2 + dipo[1, i, i]**2 + dipo[2, i, i]**2
            deg = math.sqrt(deg/(o2 - o1))
        logger.note(self, F'the global locality is {deg}')
    def pairing(self, o1, v2, _weight='distance', metd='KM'):
        nelec_a, nelec_b = self.mol.nelec
        logger.note(self, F'Match |{o1+1}->{nelec_b}| <-> |{nelec_a+1}->{v2}|, weight={_weight}')
        o2 = nelec_b; no = o2-o1
        v1 = nelec_a; nv = v2-v1

        def distance():
            from pyscf.lo.boys import dipole_integral
            # dip_ao = self.mol.intor_symmetric('int1e_r', comp=3)
            dip = dipole_integral(self.mol, self.mo_coeff)
            # o2 = self.nocc; no = o2 - o1
            # v1 = self.nocc; nv = v2 - v1

            # n=15; o=self.nocc-n
            # _dip = numpy.zeros((3,n,n))
            # for k in range(3):
            #     for i in range(n):
            #         for j in range(n):
            #             _dip[k][i][j] = dip[k][i+o][j+v1]
            # dump_ndarr(self, _dip[0], 'dipole moment X:', fmt='rec', start_row=o+1, start_col=v1+1)
            # dump_ndarr(self, _dip[1], 'dipole moment Y:', fmt='rec', start_row=o+1, start_col=v1+1)
            # dump_ndarr(self, _dip[2], 'dipole moment Z:', fmt='rec', start_row=o+1, start_col=v1+1)

            rov = numpy.zeros((no, nv))
            for i in range(no):
                o = o1 + i
                for j in range(nv):
                    v = v1 + j
                    rov[i, j] +=  (dip[0, o, o] - dip[0, v, v])**2 \
                                 + (dip[1, o, o] - dip[1, v, v])**2 \
                                 + (dip[2, o, o] - dip[2, v, v])**2
            # dump_ndarr(self, rov, 'orbital distance =', fmt='rec', start_row=o1+1, start_col=v1+1)
            # n=15
            # dump_ndarr(self, rov[no-n:,:n], 'orbital distance =', fmt='rec',start_row=o1+1,start_col=v1+1)
            return rov
        def distance2():
            from pyscf.lo.boys import dipole_integral
            # dip_ao = self.mol.intor_symmetric('int1e_r', comp=3)
            dip = dipole_integral(self.mol, self.mo_coeff)
            # o2 = self.nocc; no = o2 - o1
            # v1 = self.nocc; nv = v2 - v1
            rov = numpy.zeros((no, nv))
            for i in range(no):
                o = o1 + i
                for j in range(nv):
                    v = v1 + j
                    rov[i, j] += (dip[0, o, o] - dip[0, v, v])**2 \
                                 + (dip[1, o, o] - dip[1, v, v])**2 \
                                 + ( dip[2, o, o] - dip[2, v, v])**2 \
                                 - 2*(dip[0, o, v] ** 2 + dip[1, o, v]**2 + dip[2, o, v]**2)
            # dump_ndarr(self, rov, 'orbital distance2 =', fmt='rec',start_row=o1+1,start_col=v1+1)
            # n=15
            # dump_ndarr(self, rov[no-n:,:n], 'orbital distance2 =', fmt='rec',start_row=o1+1,start_col=v1+1)
            return rov - 2*trans_dipo()
        def angle():
            # o2 = self.nocc; no = o2 - o1
            # v1 = self.nocc; nv = v2 - v1
            mo_occ = self.mo_coeff[:, o1:o2].copy(); mo_vir = self.mo_coeff[:, v1:v2].copy()
            mo_occ = numpy.fabs(mo_occ); mo_vir = numpy.fabs(mo_vir)  # 取绝对
            aov = numpy.zeros((no, nv))
            # s = self.mol.intor_symmetric('int1e_ovlp')
            # # s_orth=reduce(numpy.dot, (mo_occ.T.conj(), s, mo_vir)) # 正交基下的重叠
            # # dump_ndarr(self, s_orth, 'orth basis overlap =', fmt='rec', start_row=o1 + 1, start_col=v1 + 1)
            # s_diag = s.diagonal()   # 非正交化
            # for i in range(self.npu):
            #     mo_occ[:,i] /= s_diag
            #     mo_vir[:,i] /= s_diag
            # s_nonorth=reduce(numpy.dot, (mo_occ.T.conj(), s, mo_vir)) # 非正交基下的重叠
            # dump_ndarr(self, s_nonorth, 'orth basis overlap =', fmt='rec', start_row=o1 + 1, start_col=v1 + 1)
            for i in range(no):
                for j in range(nv):
                    x = mo_occ[:, i]; y = mo_vir[:, j]
                    aov[i, j] = x.dot(y) / numpy.sqrt(x.dot(x) * y.dot(y))
            # dump_ndarr(self, aov, 'Cancel phase difference angle =', fmt='rec', start_row=o1+1, start_col=v1+1)
            # n=15
            # dump_ndarr(self, aov[no-15:,:n], 'Cancel phase difference angle =', fmt='rec', start_row=o1+1, start_col=v1+1)
            return aov
        def exch_inte():
            # o2 = self.nocc; no = o2 - o1
            # v1 = self.nocc; nv = v2 - v1
            # hcore = scf.hf.get_hcore(self.mol)
            # self.h_mo = reduce(numpy.dot, (self.mo_coeff.T, hcore, self.mo_coeff))
            # dump_ndarr(self, self.h_mo, 'h_mo =', fmt='rec')
            eri = self.mol.intor('int2e', aosym='s8')
            global lT0; lT0 = lowTri0(self.nobt)
            from pyscf import ao2mo
            # self.g_mo=ao2mo.outcore.full_iofree(self.mo_coeff[:,:self.nocc+self.npr],self.mo_coeff[:,:self.nocc+self.npr])
            self.g_mo = ao2mo.kernel(eri, self.mo_coeff)  # s4
            inteov = numpy.zeros((no, nv))
            for i in range(no):
                o = o1 + i
                for j in range(nv):
                    v = v1 + j
                    inteov[i, j] = self.g_mo[dei(o, v, o, v)]
            # dump_ndarr(self, inteov, 'exchange integral matix =', fmt='rec',start_row=o1+1,start_col=v1+1)
            return inteov
        def trans_dipo():
            from pyscf.lo.boys import dipole_integral
            # dip_ao = self.mol.intor_symmetric('int1e_r', comp=3)
            dip = dipole_integral(self.mol, self.mo_coeff)
            # o2 = self.nocc; no = o2 - o1
            # v1 = self.nocc; nv = v2 - v1
            dipov = numpy.zeros((no, nv))
            for i in range(no):
                o = o1 + i
                for j in range(nv):
                    v = v1 + j
                    dipov[i, j] = dip[0, o, v] ** 2 + dip[1, o, v] ** 2 + dip[2, o, v] ** 2
            # dump_ndarr(self, dipov, 'transition dipole strength =', fmt='rec', start_row=o1+1, start_col=v1+1)
            # n=15
            # dump_ndarr(self, dipov[no-n:,:n], 'transition dipole strength =', fmt='rec', start_row=o1+1, start_col=v1+1)
            return dipov

        from scipy.optimize import linear_sum_assignment
        if _weight=='trans_dipo':
            wov=trans_dipo(); pairocc, pairvir = linear_sum_assignment(-wov)
        elif _weight=='exch_inte':
            wov=exch_inte(); pairocc, pairvir = linear_sum_assignment(-wov)
        elif _weight=='angle':
            wov=angle(); pairocc, pairvir = linear_sum_assignment(-wov)
        elif _weight=='distance':
            wov=distance(); pairocc, pairvir = linear_sum_assignment(wov)
        elif _weight=='distance2':
            wov=distance2(); pairocc, pairvir = linear_sum_assignment(wov)
        else: raise ModuleNotFoundError(F'{_weight} module can\'t found')

        # dump_ndarr(self, wov, 'weight matrix =', fmt='rec', start_row=o1+1, start_col=v1+1)
        np=len(pairocc)
        for i in range(np): logger.note(self, F' {i+1:3d}: ({o1+pairocc[i]+1:3d},{v1+pairvir[i]+1:3d})   {wov[pairocc[i],pairvir[i]]}')
        logger.note(self, F'sum = {wov[pairocc,pairvir].sum()}')
        pairocc = numpy.array(list(set(range(no)) - set(pairocc)) + list(pairocc))
        pairvir = numpy.array(list(reversed(pairvir)) + list(set(range(nv)) - set(pairvir)))
        self.mo_coeff[:, o1:o2] = self.mo_coeff[:, o1+pairocc]
        self.mo_coeff[:, v1:v2] = self.mo_coeff[:, pairvir+v1]


    def chkfile(mol, chkfile_name, project=True):
        '''Read the HF results from checkpoint file, then project it to the
        basis defined by ``mol``

        Returns:
            Density matrix, 2D ndarray
        '''
        from pyscf.scf import uhf
        dm = uhf.init_guess_by_chkfile(mol, chkfile_name, project)
        return dm[0] + dm[1]

    def energy_nuc(self):
        return gto.mole.energy_nuc(self.mol)
    def energy_elec(self):
        energy = 0.0
        for i in range(2):
            for u in range(self.nocc):
                energy += 2.0 * self.ci[u,i]**2 * self.h_mo[self.pair[u,i],self.pair[u,i]]
                for j in range(2):
                    energy += self.ci[u,i] * self.ci[u,j] * \
                              self.g_mo[dei(self.pair[u,i],self.pair[u,j],self.pair[u,i],self.pair[u,j])]
                    for v in range(u+1, self.nocc):
                        energy += self.ci[u,i]**2 * self.ci[v,j]**2 * \
                                  (4*self.g_mo[dei(self.pair[u,i],self.pair[u,i],self.pair[v,j],self.pair[v,j])]
                                   -2*self.g_mo[dei(self.pair[u,i],self.pair[v,j],self.pair[u,i],self.pair[v,j])])
        return energy

    def opt_ci(self):
        for it in range(self.max_cycle[-1]):
            _ci=self.ci; _e_ele = self.e_ele
            for u in range(self.nocc-self.npr, self.nocc):
                hp = numpy.zeros((2,2))
                for i in range(2):
                    for j in range(2):
                        hp[i,j] += 2*self.h_mo[self.pair[u,i],self.pair[u,i]]*delta(i,j) \
                                    + self.g_mo[dei(self.pair[u,i],self.pair[u,j],self.pair[u,i],self.pair[u,j])]
                        if i==j:
                            for v in range(self.nocc):
                                if v!=u:
                                    for jj in range(2):
                                        hp[i,j] += _ci[v,jj]**2 * \
                                        (4 * self.g_mo[dei(self.pair[u, i], self.pair[u, i], self.pair[v, jj], self.pair[v, jj])]
                                        - 2 * self.g_mo[dei(self.pair[u, i], self.pair[v, jj], self.pair[u, i],self.pair[v, jj])])
                s, a = numpy.linalg.eig(hp)
                a=a[:,numpy.argsort(s.real)]
                self.ci[u,:] = a[:,0].T
            self.e_ele=self.energy_elec()
            if (abs(_e_ele-self.e_ele)<self.conv_tol[-1]): break
        else: logger.error(self, F'Error: Not converge after {self.max_cycle[-1]} cycles ci optimization')

    def b0(self, u, i, v, j):
        return 8 * self.ci[u, i] ** 2 * self.ci[v, j] ** 2 * (1 - delta(u, v))
    def b1(self, u, i, v, j):
        return self.ci[u, i] ** 2 + self.ci[v, j] ** 2 - 2 * self.ci[u, i] * self.ci[v, j] * delta(u, v) \
               - 2 * self.ci[u, i] ** 2 * self.ci[v, j] ** 2 * (1 - delta(u, v))
    def b2(self, u, i, v, j):
        return self.ci[u, i] ** 2 - self.ci[v, j] ** 2
    def xm(self, m, u, i, v, j):
        _xm=0.0
        for e in range(2):
            _xm += self.ci[m, e] * self.g_mo[dei(self.pair[m, e], self.pair[u, i], self.pair[m, e], self.pair[v, j])]
        return _xm
    def ym(self, _m, u, i, v, j):
        _ym=0
        for m in range(_m):
            for e in range(2):
                _ym += self.ci[m, e] ** 2 * (
                            2 * self.g_mo[dei(self.pair[m, e], self.pair[m, e], self.pair[u, i], self.pair[v, j])]
                            - self.g_mo[dei(self.pair[m, e], self.pair[u, i], self.pair[m, e], self.pair[v, j])])
        for m in range(_m+1, self.nocc):
            for e in range(2):
                _ym += self.ci[m, e] ** 2 * (
                            2 * self.g_mo[dei(self.pair[m, e], self.pair[m, e], self.pair[u, i], self.pair[v, j])]
                            - self.g_mo[dei(self.pair[m, e], self.pair[u, i], self.pair[m, e], self.pair[v, j])])
        return _ym
    def auv(self, u, i, v, j):
        return -self.b0(u, i, v, j) * (
                    self.g_mo[dei(self.pair[u, i], self.pair[v, j], self.pair[u, i], self.pair[v, j])]
                    - self.g_mo[dei(self.pair[u, i], self.pair[u, i], self.pair[v, j], self.pair[v, j])]) \
               + 2 * self.b2(u, i, v, j) * (
                           self.h_mo[self.pair[v, j], self.pair[v, j]] - self.h_mo[self.pair[u, i], self.pair[u, i]]) \
               + 2 * self.ci[u, i] * (self.xm(u, v, j, v, j) - self.xm(u, u, i, u, i)) \
               - 2 * self.ci[v, j] * (self.xm(v, v, j, v, j) - self.xm(v, u, i, u, i)) \
               + 2 * self.ci[u, i] ** 2 * (self.ym(u, v, j, v, j) - self.ym(u, u, i, u, i)) \
               - 2 * self.ci[v, j] ** 2 * (self.ym(v, v, j, v, j) - self.ym(v, u, i, u, i))
    def buv(self, u, i, v, j):
        temp = 4 * self.b2(u, i, v, j) * self.h_mo[self.pair[u, i], self.pair[v, j]] \
               + 4 * (self.ci[u, i] * self.xm(u, u, i, v, j) - self.ci[v, j] * self.xm(v, u, i, v, j)) \
               + 4 * (self.ci[u, i] ** 2 * self.ym(u, u, i, v, j) - self.ci[v, j] ** 2 * self.ym(v, u, i, v, j))
        # print('buv = ', temp)
        return temp
    def q1(self, u, i, v, j):
        return self.buv(u, i, v, j)
    def q2(self, u, i, v, j):
        return self.auv(u, i, v, j) + 2 * self.b1(u, i, v, j) \
               * (self.g_mo[dei(self.pair[u, i], self.pair[u, i], self.pair[v, j], self.pair[v, j])]
                  + self.g_mo[dei(self.pair[u, i], self.pair[v, j], self.pair[u, i], self.pair[v, j])])
    def q3(self, u, i, v, j):
        return self.buv(u, i, v, j) - 4 * self.b1(u, i, v, j) \
               * (self.g_mo[dei(self.pair[u, i], self.pair[u, i], self.pair[v, j], self.pair[u, i])]
                  - self.g_mo[dei(self.pair[v, j], self.pair[v, j], self.pair[u, i], self.pair[v, j])])
    def q4(self, u, i, v, j):
        return self.auv(u, i, v, j) + self.b1(u, i, v, j) \
               * (self.g_mo[dei(self.pair[u, i], self.pair[u, i], self.pair[u, i], self.pair[u, i])]
                  - 2 * self.g_mo[dei(self.pair[u, i], self.pair[v, j], self.pair[u, i], self.pair[v, j])]
                  + self.g_mo[dei(self.pair[v, j], self.pair[v, j], self.pair[v, j], self.pair[v, j])])

    def kernel(self):
        self.get_init()
        logger.note(self, '\n\n\n\n')

        logger.note(self, '********* Entry GVB Calculation ********')
        initfile,gvbfchk=F'{in_pre2b}_gvb_init.fchk',F'{in_pre2b}_gvb.fchk'
        inp = dump.dumpfile(gto.copy(self.mol).build(basis=self.basis), in_pre2b[:in_pre2b.rfind('_')], 'gamess', 'gvb', self.np, self.auxbasis, outsuf='init')
        datfile=F'{inp}.dat'; datfchk = F'{datfile}.fchk'
        dump.run(F"fchk2vec {initfile} {inp}.inp -gvb {self.np}")
        dump.cal(inp, '.inp')
        dump.run(F'cp {initfile} {datfchk}')
        dump.run(F'dat2fchk {datfile} {datfchk}')
        dump.run(F'cp {initfile} {gvbfchk}')
        dump.run(F'dat2fchk {datfile} {gvbfchk} -gvb {self.np}')
        load_mo(self, gvbfchk)
        from bcpt2.read_pair import read_dat
        cil = read_dat(datfile)
        for i in range(self.np): self.ci[self.nco+i] = cil[2*i:2*(i+1)]
        logger.note(self, 'GVB CI:')
        for P in range(self.nco,self.nocc): logger.note(self, F'{self.pair[P,0]+1:2d}({self.ci[P,0]:6f}) <-> {self.pair[P,1]+1:2d}({self.ci[P,1]:6f})')
        logger.note(self, '\n\n\n\n')
        
        return
        # load_mo(self, in_pre2b+'_gvb.fchk')
        # self.mol.intor_symmetric('int1e_kin')+mol.intor_symmetric('int1e_nuc')
        hcore = scf.hf.get_hcore(self.mol)
        self.h_mo = reduce(numpy.dot, (self.mo_coeff.T, hcore, self.mo_coeff))
        eri = self.mol.intor('int2e', aosym='s8')
        global lT0; lT0 = lowTri0(self.nobt)
        from pyscf import ao2mo
        # self.g_mo=ao2mo.outcore.full_iofree(self.mo_coeff[:,:self.nocc+self.npr],self.mo_coeff[:,:self.nocc+self.npr])
        self.g_mo = ao2mo.kernel(eri, self.mo_coeff)  # s4
        # eri = self.mol.intor('int2e_cart', aosym='s8')
        # eri = self.mol.intor('int2e', aosym='s8')
        # self.g_mo = ao2mo.incore.full(ao2mo.restore(8, eri, self.nbs), self.mo_coeff)

        # check HF energy
        self.e_ele = self.energy_elec()
        self.e_tot = self.e_ele+self.e_nuc
        for i in range(self.nocc-self.npr,self.nocc):
            logger.note(self, F'{self.pair[i,0]+1:2d}({self.ci[i,0]:6f}) <-> {self.pair[i,1]+1:2d}({self.ci[i,1]:6f})')
        logger.note(self, F'e_ele = {self.e_ele:10f}, e_tot = {self.e_tot:10f}')

        self.opt_ci()
        self.e_tot = self.e_ele + self.e_nuc
        for i in range(self.nocc-self.npr,self.nocc):
            logger.note(self, F'{self.pair[i,0]+1:2d}({self.ci[i,0]:6f}) <-> {self.pair[i,1]+1:2d}({self.ci[i,1]:6f})')
        logger.note(self, F'e_ele = {self.e_ele:10f}, e_tot = {self.e_tot:10f}')
        # dump_mo(self, '%s_init.fchk' %in_pre2b, reffile)

        logger.note(self, '  it       energy         delta')
        logger.note(self, '  it  obt1 orb2       q4           q3           q2           q1         theta       energy       delta(E)        fun        delta(f)')

        logger.note(self, 'start local optimization iteration')
        for it in range(self.max_cycle[0]):
            _e_ele = self.e_ele
            for u in reversed(range(self.nocc)):
                for i in range(2):
                    if u<self.nocc-self.npr and i==1: continue
                    else:
                        for v in reversed(range(self.nocc-self.npr, self.nocc)):
                            for j in range(2):
                                if u==v and i==j: continue
                                else:
                                    q1 = self.q1(u,i,v,j)
                                    q2 = self.q2(u,i,v,j)
                                    q3 = self.q3(u,i,v,j)
                                    q4 = self.q4(u,i,v,j)
                                    # theta=root(q4,q3,q2,q1)%(math.pi*2)
                                    theta=root(q4,q3,q2,q1)
                                    # if abs(theta)>self.conv_tol[1]:
                                    f=fun(q4, q3, q2, q1, theta); __e_ele=self.e_ele
                                    orb1 = self.mo_coeff[:, self.pair[u, i]].copy()
                                    orb2 = self.mo_coeff[:, self.pair[v, j]].copy()
                                    self.mo_coeff[:, self.pair[u, i]] = orb1*math.cos(theta)+orb2*math.sin(theta)
                                    self.mo_coeff[:, self.pair[v, j]] = -orb1*math.sin(theta)+orb2*math.cos(theta)
                                    self.h_mo = reduce(numpy.dot, (self.mo_coeff.T, hcore, self.mo_coeff))
                                    self.g_mo = ao2mo.kernel(eri, self.mo_coeff[:,:self.nocc+self.nactvir])
                                    deltaf=self.energy_elec()-__e_ele
                                    self.opt_ci()
                                    self.e_ele=self.energy_elec()
                                    self.e_tot=self.e_ele+self.e_nuc; deltaE=self.e_ele-__e_ele
                                    logger.note(self, F' {it+1:3d} {self.pair[u,i]+1:3d}<->{self.pair[v,j]+1:3d} {q4:12.8f} {q3:12.8f} {q2:12.8f} {q1:12.8f}'
                                                      F' {theta:12.8f} {self.e_tot:12.8f} {deltaE:12.8f} {f:12.8f} {deltaf:12.8f}')
                        for v in range(self.nocc,self.nocc+self.nactvir-self.npr):
                            j = 0
                            x11=self.xm(u,u,i,u,i); x12=self.xm(u,u,i,v,j); x22=self.xm(u,v,j,v,j); x22_11=x22-x11
                            z11=self.h_mo[self.pair[u,i],self.pair[u,i]]+self.ym(u,u,i,u,i)
                            z12=self.h_mo[self.pair[u,i],self.pair[v,j]]+self.ym(u,u,i,v,j)
                            z22=self.h_mo[self.pair[v,j],self.pair[v,j]]+self.ym(u,v,j,v,j)
                            z22_11=z22-z11
                            q1=4.*self.ci[u,i]*(x12+self.ci[u,i]*z12)
                            q2=2.*self.ci[u,i]*(x22_11+self.ci[u,i]*(z22_11
                                + self.g_mo[dei(self.pair[u,i], self.pair[v,j], self.pair[u,i], self.pair[v,j])]
                                + self.g_mo[dei(self.pair[u,i], self.pair[u,i], self.pair[v,j], self.pair[v,j])]))
                            q3=4.*self.ci[u,i]*(x12+self.ci[u,i]*(z12
                                + self.g_mo[dei(self.pair[u,i], self.pair[v,j], self.pair[v,j], self.pair[v,j])]
                                + self.g_mo[dei(self.pair[u,i], self.pair[u,i], self.pair[u,i], self.pair[v,j])]))
                            q4=self.ci[u,i]*(2.*x22_11+self.ci[u,i]*(2.*z22_11+
                                + self.g_mo[dei(self.pair[v,j], self.pair[v,j], self.pair[v,j], self.pair[v,j])]
                                -2.*self.g_mo[dei(self.pair[u,i], self.pair[v,j], self.pair[u,i], self.pair[v,j])]
                                + self.g_mo[dei(self.pair[u,i], self.pair[u,i], self.pair[u,i], self.pair[u,i])]))
                            theta=root(q4,q3,q2,q1)%(math.pi*2)
                            if abs(theta)>self.conv_tol[0]:
                                f = fun(q4, q3, q2, q1, theta); __e_ele = self.e_ele
                                orb1 = self.mo_coeff[:, self.pair[u, i]].copy()
                                orb2 = self.mo_coeff[:, self.pair[v, j]].copy()
                                self.mo_coeff[:, self.pair[u, i]] = orb1*math.cos(theta)+orb2*math.sin(theta)
                                self.mo_coeff[:, self.pair[v, j]] = -orb1*math.sin(theta)+orb2*math.cos(theta)
                                self.h_mo = reduce(numpy.dot, (self.mo_coeff.T, hcore, self.mo_coeff))
                                self.g_mo = ao2mo.kernel(eri, self.mo_coeff[:,:self.nocc+self.nactvir])
                                # eri = self.mol.intor('int2e_cart', aosym='s8')
                                # self.g_mo = ao2mo.incore.full(ao2mo.restore(8, eri, self.nbs), self.mo_coeff)
                                deltaf = self.energy_elec() - __e_ele
                                self.opt_ci()
                                self.e_ele=self.energy_elec()
                                self.e_tot=self.e_ele+self.e_nuc; deltaE=self.e_ele-__e_ele
                                logger.note(self, F' {it+1:3d} {self.pair[u,i]+1:3d}<->{self.pair[v,j]+1:3d} {q4:12.8f} {q3:12.8f} {q2:12.8f} {q1:12.8f}'
                                                  F' {theta:12.8f} {self.e_tot:12.8f} {deltaE:12.8f} {f:12.8f} {deltaf:12.8f}')
            logger.note(self, F' {it+1:3d}  {self.e_tot:11.8f}   {self.e_ele-_e_ele:12.9f}')
            if abs(_e_ele-self.e_ele)<self.conv_tol[1]: break
            # else: dump_mo(self, F'{in_pre2b}_gvb_init-II.fchk', reffile)
        else: logger.error(self, 'Not converge in local optimization')
        logger.note(self, 'The GVB pair:')
        for _i in range(self.nocc-self.npr, self.nocc):
            logger.note(self, F'{self.pair[_i,0]+1:2d}({self.ci[_i,0]:6f}) <-> {self.pair[_i,1]+1:2d}({self.ci[_i,1]:6f})')
        self.e_tot=self.e_nuc+self.e_ele
        logger.note(self, F'gvb local converged energy = {self.e_tot:11.8f}')
        # dump_mo(self, F'{in_pre2b}_gvb_init-II.fchk', reffile)

        self.h_mo = reduce(numpy.dot, (self.mo_coeff.T, hcore, self.mo_coeff))
        self.g_mo = ao2mo.kernel(eri, self.mo_coeff)
        logger.note(self, 'start global optimization iteration')
        for it in range(self.max_cycle[0]):
            _e_ele = self.e_ele
            for u in reversed(range(self.nocc)):
                for i in range(2):
                    if u<self.nocc-self.npr and i==1: continue
                    else:
                        for v in reversed(range(self.nocc-self.npr, self.nocc)):
                            for j in range(2):
                                if u==v and i==j: continue
                                else:
                                    q1=self.q1(u,i,v,j)
                                    q2=self.q2(u,i,v,j)
                                    q3=self.q3(u,i,v,j)
                                    q4=self.q4(u,i,v,j)
                                    # theta=root(q4,q3,q2,q1)%(math.pi*2)
                                    theta=root(q4,q3,q2,q1)
                                    # if abs(theta-0.0)>self.conv_tol[0]:
                                    __e_ele=self.e_ele
                                    orb1 = self.mo_coeff[:, self.pair[u, i]].copy()
                                    orb2 = self.mo_coeff[:, self.pair[v, j]].copy()
                                    self.mo_coeff[:, self.pair[u, i]] = orb1*math.cos(theta)+orb2*math.sin(theta)
                                    self.mo_coeff[:, self.pair[v, j]] = -orb1*math.sin(theta)+orb2*math.cos(theta)
                                    self.h_mo = reduce(numpy.dot, (self.mo_coeff.T, hcore, self.mo_coeff))
                                    self.g_mo = ao2mo.kernel(eri, self.mo_coeff)
                                    self.opt_ci()
                                    self.e_ele=self.energy_elec()
                                    self.e_tot=self.e_ele+self.e_nuc; deltaE=self.e_ele-__e_ele
                                    logger.note(self, F' {it+1:3d} {self.pair[u,i]+1:3d}<->{self.pair[v,j]+1:3d} {q4:12.8f} {q3:12.8f} {q2:12.8f} {q1:12.8f}'
                                                      F' {theta:12.8f} {self.e_tot:12.8f} {deltaE:12.8f} {f:12.8f}')
                        for v in range(self.nocc, self.npa):
                            j = 0
                            x11=self.xm(u,u,i,u,i); x12=self.xm(u,u,i,v,j); x22=self.xm(u,v,j,v,j); x22_11=x22-x11
                            z11=self.h_mo[self.pair[u,i],self.pair[u,i]]+self.ym(u,u,i,u,i)
                            z12=self.h_mo[self.pair[u,i],self.pair[v,j]]+self.ym(u,u,i,v,j)
                            z22=self.h_mo[self.pair[v,j],self.pair[v,j]]+self.ym(u,v,j,v,j)
                            z22_11=z22-z11
                            q1=4.*self.ci[u,i]*(x12+self.ci[u,i]*z12)
                            q2=2.*self.ci[u,i]*(x22_11+self.ci[u,i]*(z22_11
                                + self.g_mo[dei(self.pair[u,i], self.pair[v,j], self.pair[u,i], self.pair[v,j])]
                                + self.g_mo[dei(self.pair[u,i], self.pair[u,i], self.pair[v,j], self.pair[v,j])]))
                            q3=4.*self.ci[u,i]*(x12+self.ci[u,i]*(z12
                                + self.g_mo[dei(self.pair[u,i], self.pair[v,j], self.pair[v,j], self.pair[v,j])]
                                + self.g_mo[dei(self.pair[u,i], self.pair[u,i], self.pair[u,i], self.pair[v,j])]))
                            q4=self.ci[u,i]*(2.*x22_11+self.ci[u,i]*(2.*z22_11+
                                + self.g_mo[dei(self.pair[v,j], self.pair[v,j], self.pair[v,j], self.pair[v,j])]
                                -2.*self.g_mo[dei(self.pair[u,i], self.pair[v,j], self.pair[u,i], self.pair[v,j])]
                                + self.g_mo[dei(self.pair[u,i], self.pair[u,i], self.pair[u,i], self.pair[u,i])]))
                            # theta=root(q4,q3,q2,q1)%(math.pi*2)
                            theta=root(q4,q3,q2,q1)
                            # if abs(theta-0.0)>self.conv_tol[0]:
                            __e_ele=self.e_ele
                            orb1 = self.mo_coeff[:, self.pair[u, i]].copy()
                            orb2 = self.mo_coeff[:, self.pair[v, j]].copy()
                            self.mo_coeff[:, self.pair[u, i]] = orb1*math.cos(theta)+orb2*math.sin(theta)
                            self.mo_coeff[:, self.pair[v, j]] = -orb1*math.sin(theta)+orb2*math.cos(theta)
                            self.h_mo = reduce(numpy.dot, (self.mo_coeff.T, hcore, self.mo_coeff))
                            self.g_mo = ao2mo.kernel(eri, self.mo_coeff)
                            # eri = self.mol.intor('int2e_cart', aosym='s8')
                            # self.g_mo = ao2mo.incore.full(ao2mo.restore(8, eri, self.nbs), self.mo_coeff)
                            self.opt_ci()
                            self.e_ele=self.energy_elec()
                            self.e_tot=self.e_ele+self.e_nuc; deltaE=self.e_ele-__e_ele
                            logger.note(self, F' {it+1:3d} {self.pair[u,i]+1:3d}<->{self.pair[v,j]+1:3d} {q4:12.8f} {q3:12.8f} {q2:12.8f} {q1:12.8f}'
                                              F' {theta:12.8f} {self.e_tot:12.8f} {deltaE:12.8f} {f:12.8f}')
            logger.note(self, F' {it+1:3d}  {self.e_tot:11.8f}   {self.e_ele-_e_ele:12.9f}')
            if abs(_e_ele-self.e_ele)<self.conv_tol[0]: break
            # else: dump_mo(self, F'{in_pre2b}_gvb_init-III.fchk', reffile)
        else: raise RuntimeError('Not converge in global optimization')
        logger.note(self, 'The GVB pair:')
        for _i in range(self.nocc-self.npr, self.nocc):
            logger.note(self, F'{self.pair[_i,0]+1:2d}({self.ci[_i,0]:6f}) <-> {self.pair[_i,1]+1:2d}({self.ci[_i,1]:6f})')
        self.e_tot=self.e_nuc+self.e_ele
        logger.note(self, F'gvb global converged energy = {self.e_tot:11.8f}')
        # dump_mo(self, F'{in_pre2b}_gvb_init-III.fchk', reffile)









if __name__ == "__main__":
    mol=gto.M(atom=
              '''
                H 0.0 0.0 0.0
                F 0.0 0.0 1.4
              ''',
              basis='6-31g',
              output='HF_1.4_6-31g_gvb.out')
    print('start')
    mf=GVB(mol)
    mf.kernel()
    print('end')










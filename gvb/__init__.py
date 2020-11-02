#!/usr/bin/env python 
# -*- coding: utf-8 -*- 
#
#      Copyright:   Qingchun Wang @ NJU
#      File Name:   __init__.py.py
#            Des:   __init__.py
#           Mail:   qingchun720@foxmail.com
#   Created Time:   11:10 二月-25/2019 
#


'''
Non-relativistic and relativistic Generalized Valence Bond Theory (GVB)
    Note:
        As CASSCF, GVB calculation require an basis HF/ROHF/UHF to obtain initial guess.
        But GVB could start with mol object, which dou't carry out HF/ROHF/UHF calculation.
===========================================================================================

Simple usage::

    # start from mol object
    >>> from pyscf import gto, scf, gvb
    >>> mol = gto.M(atom='H 0 0 0; H 0 0 1')
    >>> f = gvb.GVB(mol).run()

    or:

    # start from HF/ROHF/UHF
    >>> from pyscf import gto, scf, gvb
    >>> mol = gto.M(atom='H 0 0 0; H 0 0 1')
    >>> _f = scf.RHF(mol).run()
    >>> f = gvb.GVB(_f).run()


:func:`gvb.GVB` returns an proper instance of GVB class.
There are some parameters to control the GVB method.

    verbose : int
        Print level.  Default value equals to :class:`Mole.verbose`
    max_memory : float or int
        Allowed memory in MB.  Default value equals to :class:`Mole.max_memory`
    chkfile : str
        checkpoint file to save MOs, orbital energies etc.

    npair : None, int
        The number of pair in GVB. if it is None, this indicate auto detetmine
        Default is None.
    init : None, ndarray or str
        The init guess for GVB calculation. if it is None, this indicate detetminate automatically,
        and str reading from file.
        Default is None
    init_level : int
        The levle of init. As shown below
        mo.hf(0) -> mo.svo(1) -> mo.local(2) -> mo.pair(3)
        Default is 0.
    mothod : str
        orbital optimization method. It can be one of 'JR', 'OCBSE', 'SOSCF'.
        Default is 'JR'
    auxbasis : str
        Auxiliary basis is reference basis set used for projection, e.g
        sto-3g or sto-6g (default) in general.
    max_cycle : (int, int, int), for JR method only
        method max number for (global, local, theta) optimization.
        Default is (500, 200, 50)
    conv_tol : (float, float, float), for JR method only
        converge threshold to (global, local, theta) optimization.
        Default is (1e-8, 1e-6, 1e-8)

    local : str, for npair = None only
        orbital localization method. It can be one of 'Boys', 'PM', 'ER', 'NAO', 'IAO'
        Default is 'Boys'
    pop: str, for Pipek-Mezey localization only
        Population analysis method. It can be one of 'mulliken', 'lowdin', 'meta_lowdin'
        Default is 'mulliken'

    irrep_nelec : dict, for symmetry- RGVB/UGVB class only
        to indicate the number of electrons for each irreps.
        In RGVB, give {'ir_name':int, ...} ;
        In UGVB, give {'ir_name':(int,int), ...} .
        It is effective when :attr:`Mole.symmetry` is set ``True``.

Saved results

    converged : bool
        SCF converged or not
    e_tot : float
        Total GVB energy (electronic energy plus nuclear repulsion)
    mo_energy :
        Orbital energies
    mo_occ
        Orbital occupancy
    mo_coeff
        Orbital coefficients
'''

__version__ = '2.0'

from pyscf import gto
from gvb import rgvb
from gvb import rgvb_symm
from gvb import ugvb
from gvb import ugvb_symm

def GVB(mf, np=None, **kwargs):
    # if isinstance(mf, gto.Mole): # mf <- mol
    #     f = scf.UHF(mf).run()
    #     mo_i, mo_e, stable_i, stable_e = f.stability()
    #     if not stable_i:
    #         for i in range(10): # stable ?
    #             f = scf.newton(f).run(mo_i, f.mo_occ)
    #             mo_i, mo_e, stable_i, stable_e = f.stability()
    #             if stable_i: break
    #         else: raise RuntimeError('UHF stability opt failed over 10 times')
    #     occ = (f.mo_coeff[0][:, f.mo_occ[0] > 0], f.mo_coeff[1][:, f.mo_occ[1] > 0])
    #     if mf.cart: spin2 = scf.uhf.spin_square(occ, f.mol.intor('int1e_ovlp_cart'))
    #     else: spin2 = scf.uhf.spin_square(occ, f.mol.intor('int1e_ovlp_sph'))
    #     if abs(spin2[0] - 0.) < 1e-8:
    #         f = scf.RHF(mf).kernel(scf.rhf.uhf_to_rhf(f).make_rdm1())
    #         mo_i, mo_e, stable_i, stable_e = f.stability()
    #         if not stable_e:
    #             raise RuntimeError('Error ! After UHF stability opt, RHF stability show RHF -> UHF once again')
    #         if not stable_i:
    #             for i in range(10):  # stable ?
    #                 f = scf.newton(f).run(mo_i, f.mo_occ)
    #                 mo_i, mo_e, stable_i, stable_e = f.stability()
    #                 if stable_i: break
    #             else:
    #                 raise RuntimeError('RHF stability opt failed over 10 times')
    # else:   # mf <- scf
    #     if isinstance(mf, scf.uhf.UHF):   # mf <- uhf
    #         mo_i, mo_e, stable_i, stable_e = mf.stability()
    #         if not stable_i:
    #             print('The UHF wave function provided is\'t stable, and will be optimized automatically')
    #             for i in range(10):  # stable ?
    #                 f = scf.newton(mf).run(mo_i, f.mo_occ)
    #                 mo_i, mo_e, stable_i, stable_e = f.stability()
    #                 if stable_i: break
    #             else:
    #                 raise RuntimeError('UHF stability opt failed over 10 times')
    #         occ = (f.mo_coeff[0][:, f.mo_occ[0] > 0], f.mo_coeff[1][:, f.mo_occ[1] > 0])
    #         if mf.mol.cart:  spin2 = scf.uhf.spin_square(occ, f.mol.intor('int1e_ovlp_cart'))
    #         else: spin2 = scf.uhf.spin_square(occ, f.mol.intor('int1e_ovlp_sph'))
    #         if abs(spin2[0] - 0.) < 1e-8:
    #             f = scf.RHF(mf.mol).kernel(f.make_rdm1())
    #             mo_i, mo_e, stable_i, stable_e = f.stability()
    #             if not stable_e:
    #                 raise RuntimeError('Error ! After UHF stability opt, RHF stability show RHF -> UHF once again')
    #             if not stable_i:
    #                 for i in range(10):  # stable ?
    #                     f = scf.newton(f).run(mo_i, f.mo_occ)
    #                     mo_i, mo_e, stable_i, stable_e = f.stability()
    #                     if stable_i: break
    #                 else:
    #                     raise RuntimeError('RHF stability opt failed over 10 times')
    #     else:   # mf <- rhf
    #         mo_i, mo_e, stable_i, stable_e = mf.stability()
    #         if not stable_i:
    #             print('The RHF wave function provided is\'t stable, and will be optimized automatically')
    #             for i in range(10):  # stable ?
    #                 f = scf.newton(mf).run(mo_i, f.mo_occ)
    #                 mo_i, mo_e, stable_i, stable_e = f.stability()
    #                 if stable_i: break
    #             else:
    #                 raise RuntimeError('RHF stability opt failed over 10 times')
    #         if not stable_e:
    #             f = scf.UHF(mf.mol).kernel(scf.uhf.rhf_to_uhf(mf).make_rdm1())
    #             mo_i, mo_e, stable_i, stable_e = f.stability()
    #             if not stable_i:
    #                 for i in range(10):  # stable ?
    #                     f = scf.newton(f).run(mo_i, f.mo_occ)
    #                     mo_i, mo_e, stable_i, stable_e = f.stability()
    #                     if stable_i: break
    #                 else:
    #                     raise RuntimeError('UHF stability opt failed over 10 times')
    if isinstance(mf, gto.Mole):
        if mf.symmetry:
            return rgvb_symm.GVB(mf, np=np, **kwargs)
        else:
            return rgvb.GVB(mf, np=np)
    else:
        if mf.mol.symmetry:
            return rgvb_symm.GVB(mf, np=np, **kwargs)
        else:
            return rgvb.GVB(mf, np=np)
RGVB = GVB

def UGVB(mf, np=None, **kwargs):
    if isinstance(mf, gto.Mole):
        if mf.symmetry:
            return ugvb_symm.UGVB(mf, np=np, **kwargs)
        else:
            return ugvb.UGVB(mf, np=np, **kwargs)
    else:
        if mf.mol.symmetry:
            return ugvb_symm.UGVB(mf, np=np, **kwargs)
        else:
            return ugvb.UGVB(mf, np=np, **kwargs)




#!/usr/bin/env python 
# -*- coding: utf-8 -*- 
#
#      Copyright:   Qingchun Wang @ NJU
#      File Name:   ugvb.py
#            Des:   ugvb
#           Mail:   qingchun720@foxmail.com
#   Created Time:   11:24 二月-25/2019 
#


import sys, os, subprocess
from pyscf import lib
from pyscf import gto, scf
import tempfile
from pyscf.lib import logger
import numpy, math, scipy
from functools import reduce
from gvb import dump
from gvb import rgvb





au2debye = 2.541766; i_lowTri = []
def index(i, j, k, l):
    if i >= j:
        ij = i_lowTri[i]+j
    else:
        ij = i_lowTri[j]+i
    if k >= l:
        kl = i_lowTri[k]+l
    else:
        kl = i_lowTri[l]+k
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


in_pre2b=None; reffile=None


class UGVB(rgvb.GVB):
    ''' UGVB base class. Non-relativistic UGVB

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
            0(mol) ->  1(rhf, uno) -> 2(pao) -> 3(svo) -> 4(local) -> 5(init-I) -> 6(init-II) -> 7(gvb)
            Default is 0.
        init_in : end
            The levle of init for output orbitals
            The levle of init. As shown below
            0(mol) ->  1(rhf, uno) -> 2(pao) -> 3(svo) -> 4(local) -> 5(init-I) -> 6(init-II) -> 7(gvb)
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
        >>> from pyscf import gto, scf
        >>> from gvb import rgvb
        >>> mol = gto.M(atom='H 0 0 0; H 0 0 1')
        >>> mf = ugvb.UGVB(mol).run()

        or:

        # start from HF/ROHF/UHF
        >>> from pyscf import gto, scf
        >>> from gvb import rgvb
        >>> mol = gto.M(atom='H 0 0 0; H 0 0 1')
        >>> _mf = scf.RHF(mol).run()
        >>> mf = ugvb.UGVB(_mf).run()

    '''
    def __init__(self, mf, np=None, init=None, init_in=0, sch='sch-ii', init_out=5, auxbasis='sto-6g',
                 max_cycle=(500, 200, 50), conv_tol=(1e-8, 1e-4, 1e-8),
                 pop='dipole', weight='trans_dipo'):
        rgvb.GVB.__init__(self,mf,np)
        self.nsig=0

        self._keys = set(self.__dict__.keys())

    def init_by_auto(self):
        global in_pre2b, reffile
        if isinstance(self.mol, gto.Mole):  # mf -> mol
            self._chkfile = tempfile.NamedTemporaryFile(dir=lib.param.TMPDIR)
            self.chkfile = self._chkfile.name
            self.basis=self.mol.basis
            self.mol.build(basis=rgvb.parse(self.mol))
            in_pre2b = self.mol.output[:self.mol.output.rfind('_')]; reffile = in_pre2b+'_rhf.fchk'
            mf=rgvb.hf(self.mol, self.basis, self.mol.cart, 'uhf')
            if self.sch=='sch-i' or mf.spin_square()[1]-mf.mol.spin-1 < 0.1:  # 2S+1 # Sch-I
                mf=rgvb.hf(self.mol, self.basis, self.mol.cart, 'rhf')
        else:  # mf -> scf
            mf = self.mol; self.mol = mf.mol
            self.chkfile=mf.chkfile
            self.basis=self.mol.basis
            self.mol.build(basis=rgvb.parse(self.mol))
            in_pre2b = self.mol.output[:self.mol.output.rfind('_')]; reffile = in_pre2b+'_rhf.fchk'
            rgvb.converge(mf)
            rgvb.stabilize(mf)
        if isinstance(mf, scf.uhf.UHF):
            self.npair_u=rgvb.uno(mf)
            if not os.path.exists(reffile):
                dump.execute('fchk_uhf2rhf %s_uhf.fchk' %in_pre2b)
                dump.execute('mv %s_uhf_r.fchk %s' %(in_pre2b, reffile))
            rgvb.dump_mo(mf, '%s_gvb_uno.fchk' %in_pre2b, reffile)
            mf.mo_coeff=mf.mo_coeff[0]; mf.mo_occ=mf.mo_occ[0]
        self.mo_coeff = mf.mo_coeff; self.mo_occ = mf.mo_occ
        self.e_tot = mf.e_tot; self.e_nuc = self.mol.energy_nuc(); self.e_ele = self.e_tot-self.e_nuc

        self.nbs, self.nobt = self.mo_coeff.shape; nelec_a, nelec_b= self.mol.nelec; self.nsig=self.mol.spin
        self.nocc = nelec_a; self.nvir = self.nobt - nelec_a
        nfrez = dump.ncore(self.mol); nobt_s=dump.nrefobt(self.mol, in_pre2b[:in_pre2b.rfind('_')])
        self.nactocc = nelec_b-nfrez; self.nactvir = nobt_s-nelec_a
        self.npair_d = min(self.nactocc,self.nactvir); self.npair_t = self.nobt-self.npair_d
        logger.note(self, 'np = %d', self.npair_d)

        logger.note(self,'nbs= %3d      nobt= %3d     nobt_s= %3d', self.nbs,self.nobt,nobt_s)
        logger.note(self,'nfrez= %2d    nactocc= %2d  nactvir= %2d', nfrez,self.nactocc,self.nactvir)
        logger.note(self,'npair_u= %2d  npair_d= %2d  npair_t= %2d  nsig=%d', self.npair_u,self.npair_d,self.npair_t,self.nsig)

        self.pair = numpy.zeros((self.npair_t, 2), dtype=int)
        for i in range(self.nocc): self.pair[i, 0] = i
        for i in range(self.npair_d): self.pair[self.nocc-1-i, 1] = self.nocc+i
        for i in range(self.nocc, self.npair_t): self.pair[i, 0] = i+self.npair_d
        self.ci = numpy.zeros((self.npair_t, 2)); self.ci[:,0] = 1.
        if self.npair_d>self.npair_u:
            if self.nvir>self.nactvir:
                self.svo()
                rgvb.dump_mo(self, '%s_gvb_svo.fchk' %in_pre2b, reffile)

            # self.localize(nfrez, nobt_s)
            self.mo_coeff=rgvb.localize(self.mol,self.mo_coeff, nfrez, nobt_s, self.pop)
            rgvb.dump_mo(self, '%s_gvb_localize.fchk' %in_pre2b, reffile)
            # self.pairing(nfrez, nobt_s, weight='distance')
            # self.pairing(nfrez, nobt_s, weight='distance2')
            # self.pairing(nfrez, nobt_s, weight='angle')
            # # self.pairing(nfrez, nobt_s, weight='exch_inte')
            self.pairing(nfrez, nobt_s, self.weight)
        else:
            # self.localize(nelec_b-self.npair_u, nelec_a+self.npair_u)
            self.mo_coeff=rgvb.localize(self.mol, self.mo_coeff, nelec_b-self.npair_u, nelec_a+self.npair_u, self.pop)
            rgvb.dump_mo(self, '%s_gvb_localize.fchk' %in_pre2b, reffile)
            self.pairing(nelec_b-self.npair_u, nelec_a+self.npair_u, self.weight)
            # self.pairing(self.nocc-self.npair_u, self.nocc+self.npair_u, weight='distance2')
            # self.pairing(self.nocc-self.npair_u, self.nocc+self.npair_u, weight='angle')
            # # self.pairing(self.nocc-self.npair_u, self.nocc+self.npair_u, weight='exch_inte')
            # self.pairing(self.nocc-self.npair_u, self.nocc+self.npair_u)
        # rgvb.dump_mo(self, '%s_gvb.fchk' %in_pre2b, reffile)
        self.mo_coeff[:,nelec_b-self.npair_d:nelec_a]=self.mo_coeff[:,list(range(nelec_b, nelec_a))+list(range(nelec_b-self.npair_d,nelec_b))]
        self.ci[self.nocc-self.npair_d-self.nsig:self.nocc-self.npair_d, 0]=1/math.sqrt(2.)
        # dump_ndarr(self, self.mo_coeff[:,nfrez:nobt_s], 'vir =', fmt='rec', start_col=nfrez+1)
        rgvb.dump_mo(self, '%s_gvb_init-I.fchk' %in_pre2b, reffile)

    def s(self, u):
        sig=self.nocc-self.npair_d-self.nsig; gem=self.nocc-self.npair_d
        if sig<=u<gem: return 0.
        else: return 1.
    def t(self, u, v):
        sig=self.nocc-self.npair_d-self.nsig; gem=self.nocc-self.npair_d
        if (sig<=u<gem) and (sig<=v<gem): return 2.
        else: return 1.
    def energy_nuc(self):
        return gto.mole.energy_nuc(self.mol)
    def energy_elec(self):
        energy = 0.0; sig_s=self.nocc-self.nsig-self.npair_d; sig_e=self.nocc-self.npair_d
        for u in range(self.nocc):
            for i in range(2):
                energy += 2.0 * self.ci[u,i]**2 * self.h_mo[self.pair[u,i],self.pair[u,i]]
                s=self.s(u)
                for j in range(2):
                        energy += s*self.ci[u,i] * self.ci[u,j] * \
                                  self.g_mo[index(self.pair[u,i],self.pair[u,j],self.pair[u,i],self.pair[u,j])]
                for v in range(u+1, self.nocc):
                    t=self.t(u,v)
                    for j in range(2):
                        energy += 2.*self.ci[u,i]**2 * self.ci[v,j]**2 * \
                                  (2.*self.g_mo[index(self.pair[u,i],self.pair[u,i],self.pair[v,j],self.pair[v,j])]
                                   -t*self.g_mo[index(self.pair[u,i],self.pair[v,j],self.pair[u,i],self.pair[v,j])])
        # energy = 0.0; sig=self.nocc-self.npair_d-self.nsig; gem=self.nocc-self.npair_d; nocc=self.nocc
        # for u in range(sig): # core
        #     energy += 2.0 * self.h_mo[self.pair[u, 0], self.pair[u, 0]]
        #     for v in range(sig): # - core
        #         energy += (2.*self.g_mo[index(self.pair[u,0],self.pair[u,0],self.pair[v,0],self.pair[v,0])]
        #                         -self.g_mo[index(self.pair[u,0],self.pair[v,0],self.pair[u,0],self.pair[v,0])])
        #     for v in range(sig, gem):
        #         energy += (2.*self.g_mo[index(self.pair[u,0],self.pair[u,0],self.pair[v,0],self.pair[v,0])]
        #                       -self.g_mo[index(self.pair[u,0],self.pair[v,0],self.pair[u,0],self.pair[v,0])])
        #     for v in range(gem, nocc):                        # - geminal
        #         for j in range(2):
        #             energy += 2.*self.ci[v,j]**2 *(2.*self.g_mo[index(self.pair[u,0],self.pair[u,0],self.pair[v,j],self.pair[v,j])]
        #                                              -self.g_mo[index(self.pair[u,0],self.pair[v,j],self.pair[u,0],self.pair[v,j])])
        # for u in range(sig, gem): # open
        #     energy += self.h_mo[self.pair[u, 0], self.pair[u, 0]]
        #     for v in range(u+1, gem):
        #         energy += (self.g_mo[index(self.pair[u,0],self.pair[u,0],self.pair[v,0],self.pair[v,0])]
        #                   -self.g_mo[index(self.pair[u,0],self.pair[v,0],self.pair[u,0],self.pair[v,0])])
        #     for v in range(gem, nocc):
        #         for j in range(2):
        #             energy += self.ci[v,j]**2 *(2.*self.g_mo[index(self.pair[u,0],self.pair[u,0],self.pair[v,j],self.pair[v,j])]
        #                                           -self.g_mo[index(self.pair[u,0],self.pair[v,j],self.pair[u,0],self.pair[v,j])])
        # for u in range(gem, nocc): # geminal
        #     for i in range(2):
        #         energy += 2.0*self.ci[u, i]**2 * self.h_mo[self.pair[u, i], self.pair[u, i]]
        #         for j in range(2):
        #             energy += self.ci[u,i] * self.ci[u,j] * self.g_mo[index(self.pair[u,i],self.pair[u,j],self.pair[u,i],self.pair[u,j])]
        #         for v in range(u+1, nocc):
        #             for j in range(2):
        #                 energy += 2.*self.ci[u,i]**2 * self.ci[v,j]**2 * \
        #                            (2.*self.g_mo[index(self.pair[u,i],self.pair[u,i],self.pair[v,j],self.pair[v,j])]
        #                               -self.g_mo[index(self.pair[u,i],self.pair[v,j],self.pair[u,i],self.pair[v,j])])
        return energy
    def opt_ci(self):
        for it in range(self.max_cycle[-1]):
            _ci=self.ci; _e_ele = self.e_ele
            for u in range(self.nocc-self.npair_d, self.nocc):
                hp = numpy.zeros((2,2))
                for i in range(2):
                    for j in range(2):
                        hp[i,j] += 2*self.h_mo[self.pair[u,i],self.pair[u,i]]*delta(i,j) \
                                    + self.g_mo[index(self.pair[u,i],self.pair[u,j],self.pair[u,i],self.pair[u,j])]
                        if i==j:
                            for v in range(self.nocc):
                                t=self.t(u,v)
                                if v!=u:
                                    for jj in range(2):
                                        hp[i,j] += _ci[v,jj]**2 * \
                                        (4*self.g_mo[index(self.pair[u,i], self.pair[u,i], self.pair[v,jj], self.pair[v,jj])]
                                        - 2*t*self.g_mo[index(self.pair[u,i], self.pair[v,jj], self.pair[u,i],self.pair[v,jj])])
                s, a = numpy.linalg.eig(hp)
                a=a[:,numpy.argsort(s.real)]
                self.ci[u,:] = a[:,0].T
            self.e_ele=self.energy_elec()
            if (abs(_e_ele-self.e_ele)<self.conv_tol[-1]): break
        else: logger.error(self, 'Error: Not converge after %d cycles ci optimization' %self.max_cycle[-1])
    def x(self, p, u, i, v, j):
        temp = 0.
        for k in range(2):
            temp += self.ci[p,k]*self.g_mo[index(self.pair[p,k], self.pair[u,i], self.pair[p,k], self.pair[v,j])]
        return temp
    def y(self, p, u, i, v, j):
        temp = 0.
        for q in range(p):
            t=self.t(p, q)
            for k in range(2):
                temp += self.ci[q,k]**2 * (2*self.g_mo[index(self.pair[q,k], self.pair[q,k], self.pair[u,i], self.pair[v,j])]
                                          -t*self.g_mo[index(self.pair[q,k], self.pair[u,i], self.pair[q,k], self.pair[v,j])])
        for q in range(p+1, self.nocc):
            t = self.t(p, q)
            for k in range(2):
                temp += self.ci[q,k]**2 * (2*self.g_mo[index(self.pair[q,k], self.pair[q,k], self.pair[u,i], self.pair[v,j])]
                                          -t*self.g_mo[index(self.pair[q,k], self.pair[u,i], self.pair[q,k], self.pair[v,j])])
        return temp
    def z(self, p, u, i, v, j):
        return self.h_mo[self.pair[u,i], self.pair[v,j]] + self.y(p, u, i, v, j)
    def w(self, p, k, u, i, v, j):
        if self.s(p): return self.ci[p,k]**2 * self.z(p, u, i, v, j) + self.ci[p,k]*self.x(p, u, i, v, j)
        else: return self.ci[p,k]**2 * self.z(p, u, i, v, j)

    def kernel(self):
        self.get_init()
        sys.exit()
        # rgvb.load_mo(self, in_pre2b+'_gvb.fchk')
        # self.mo_coeff[:,[6,7,8]]=self.mo_coeff[:,[7,8,6]]  # swap orbital
        # self.ci[7,0]=0.998570; self.ci[7,1]=-0.053452
        # self.ci[6,0]=0.981710; self.ci[6,1]=-0.190382
        # self.mol.intor_symmetric('int1e_kin')+mol.intor_symmetric('int1e_nuc')
        hcore = scf.hf.get_hcore(self.mol)
        self.h_mo = reduce(numpy.dot, (self.mo_coeff.T, hcore, self.mo_coeff))
        eri = self.mol.intor('int2e', aosym='s8')
        global i_lowTri; i_lowTri = [0]*self.nobt
        for i in range(1, self.nobt): i_lowTri[i] = i_lowTri[i-1]+i
        from pyscf import ao2mo
        # self.g_mo=ao2mo.outcore.full_iofree(self.mo_coeff[:,:self.nocc+self.npair_d],self.mo_coeff[:,:self.nocc+self.npair_d])
        self.g_mo = ao2mo.kernel(eri, self.mo_coeff)  # s4
        # eri = self.mol.intor('int2e_cart', aosym='s8')
        # eri = self.mol.intor('int2e', aosym='s8')
        # self.g_mo = ao2mo.incore.full(ao2mo.restore(8, eri, self.nbs), self.mo_coeff)

        # check HF energy
        for i in range(self.nocc):
            logger.note(self, '%2d(%6f) <-> %2d(%6f)' % (self.pair[i, 0] + 1, self.ci[i, 0], self.pair[i, 1] + 1, self.ci[i, 1]))
        self.e_ele = self.energy_elec()
        self.e_tot = self.e_ele+self.e_nuc
        logger.note(self, 'e_ele = %10f, e_tot = %.10f' %(self.e_ele,self.e_tot))

        self.opt_ci()
        self.e_tot = self.e_ele + self.e_nuc
        for i in range(self.nocc-self.npair_d,self.nocc):
            logger.note(self, '%2d(%6f) <-> %2d(%6f)' %(self.pair[i,0]+1, self.ci[i,0], self.pair[i,1]+1, self.ci[i,1]))
        logger.note(self, 'e_ele = %10f, e_tot = %.10f' % (self.e_ele, self.e_tot))
        # rgvb.dump_mo(self, '%s_init.fchk' %in_pre2b, reffile)



        logger.note(self, 'start local optimization iteration')
        logger.note(self, '  it       energy         delta')
        logger.note(self, '  it  obt1 orb2       q4           q3           q2           q1         theta       energy       delta(E)        fun        delta(f)')
        for it in range(self.max_cycle[0]):
            _e_ele = self.e_ele
            for u in reversed(range(self.nocc)):
                for i in range(2):
                    if u < self.nocc-self.npair_d and i==1: continue
                    for v in reversed(range(self.nocc-self.npair_d-self.nsig, self.nocc)):
                        for j in range(2):
                            if v < self.nocc - self.npair_d and j == 1: # 注意：单占据的geminal，只能j=0，不能j=1,
                                continue # 因为默认单占据s=0，若j取到1时，也被设置为s=0，而vj=1是core轨道，双占据，就不对了
                            if u==v and i==j: continue
                            b1 = 8*self.ci[u,i]**2 * self.ci[v,j]**2 * (1.-delta(u,v))
                            b2 = (self.s(u)*self.ci[u, i]**2 + self.s(v)*self.ci[v, j]**2
                                       - 2*self.s(u)*self.ci[u, i]*self.ci[v, j]*delta(u, v)
                                       - 2 * self.ci[u, i]**2 * self.ci[v, j]**2 * (2-self.t(u,v)) * (1.-delta(u,v)))
                            auv = (b1*(self.g_mo[index(self.pair[u,i], self.pair[u,i], self.pair[v,j], self.pair[v,j])]
                                     - self.g_mo[index(self.pair[u,i], self.pair[v,j], self.pair[u,i], self.pair[v,j])])
                                     +2.*(self.w(u,i,v,j,v,j) - self.w(u,i,u,i,u,i) - self.w(v,j,v,j,v,j) + self.w(v,j,u,i,u,i)))
                            buv = self.w(u,i,u,i,v,j) - self.w(v,j,u,i,v,j)
                            q1 = 4.*buv
                            q2 = auv + 2.*b2*(self.g_mo[index(self.pair[u,i], self.pair[u,i], self.pair[v,j], self.pair[v,j])]
                                            + self.g_mo[index(self.pair[u,i], self.pair[v,j], self.pair[u,i], self.pair[v,j])])
                            q3 = 4.*(buv + b2*(self.g_mo[index(self.pair[u,i], self.pair[v,j], self.pair[v,j], self.pair[v,j])]
                                             - self.g_mo[index(self.pair[u,i], self.pair[u,i], self.pair[u,i], self.pair[v,j])]))
                            q4 = auv + b2*(self.g_mo[index(self.pair[u,i], self.pair[u,i], self.pair[u,i], self.pair[u,i])]
                                         + self.g_mo[index(self.pair[v,j], self.pair[v,j], self.pair[v,j], self.pair[v,j])]
                                       -2.*self.g_mo[index(self.pair[u,i], self.pair[v,j], self.pair[u,i], self.pair[v,j])])
                            # theta=root(q4,q3,q2,q1)%(math.pi*2)
                            theta=root(q4,q3,q2,q1)
                            # if abs(theta)>self.conv_tol[1]:
                            f=fun(q4, q3, q2, q1, theta); __e_ele=self.e_ele
                            orb1 = self.mo_coeff[:, self.pair[u, i]].copy()
                            orb2 = self.mo_coeff[:, self.pair[v, j]].copy()
                            self.mo_coeff[:, self.pair[u, i]] = orb1*math.cos(theta)+orb2*math.sin(theta)
                            self.mo_coeff[:, self.pair[v, j]] = -orb1*math.sin(theta)+orb2*math.cos(theta)
                            self.h_mo = reduce(numpy.dot, (self.mo_coeff.T, hcore, self.mo_coeff))
                            # self.g_mo = ao2mo.outcore.full_iofree(
                            #     self.mo_coeff[:, :self.nocc + self.npair_d],
                            #     self.mo_coeff[:, :self.nocc + self.npair_d])
                            self.g_mo = ao2mo.kernel(eri, self.mo_coeff[:,:self.nocc+self.nactvir])
                            # eri = self.mol.intor('int2e_cart', aosym='s8')
                            # self.g_mo = ao2mo.incore.full(ao2mo.restore(8, eri, self.nbs), self.mo_coeff)
                            deltaf=self.energy_elec()-__e_ele
                            self.opt_ci()
                            self.e_ele=self.energy_elec()
                            self.e_tot=self.e_ele+self.e_nuc; deltaE=self.e_ele-__e_ele
                            logger.note(self,' %3d %3d<->%3d %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f'
                              %(it+1,self.pair[u,i]+1,self.pair[v,j]+1,q4,q3,q2,q1,theta,self.e_tot,deltaE,f,deltaf))
                    for v in range(self.nocc,self.nocc+self.nactvir-self.npair_d):
                        j = 0
                        b2 = self.s(u)*self.ci[u,i]**2
                        auv = self.w(u,i,v,j,v,j) - self.w(u,i,u,i,u,i)
                        buv = self.w(u,i,u,i,v,j)
                        q1 = 4.*buv
                        q2 = 2.*(auv + b2*(self.g_mo[index(self.pair[u,i], self.pair[u,i], self.pair[v,j], self.pair[v,j])]
                                         + self.g_mo[index(self.pair[u,i], self.pair[v,j], self.pair[u,i], self.pair[v,j])]))
                        q3 = 4.*(buv + b2*(self.g_mo[index(self.pair[u,i], self.pair[v,j], self.pair[v,j], self.pair[v,j])]
                                         - self.g_mo[index(self.pair[u,i], self.pair[u,i], self.pair[u,i], self.pair[v,j])]))
                        q4 = 2.*auv + b2*(self.g_mo[index(self.pair[u,i], self.pair[u,i], self.pair[u,i], self.pair[u,i])]
                                        + self.g_mo[index(self.pair[v,j], self.pair[v,j], self.pair[v,j], self.pair[v,j])]
                                      -2.*self.g_mo[index(self.pair[u,i], self.pair[v,j], self.pair[u,i], self.pair[v,j])])
                        # theta=root(q4,q3,q2,q1)%(math.pi*2)
                        theta=root(q4,q3,q2,q1)
                        # if abs(theta)>self.conv_tol[1]:
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
                        # logger.note(self,' %3d %3d<->%3d %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f'
                        #   %(it+1,self.pair[u,i]+1,self.pair[v,j]+1,q4,q3,q3,q2,theta,self.e_tot,deltaE, f, deltaf))
            # for _i in range(self.nocc-self.npair_d,self.nocc):
            #     print('%2d(%6f) <-> %2d(%6f)' %(self.pair[_i,0]+1, self.ci[_i,0], self.pair[_i,1]+1, self.ci[_i,1]))
            logger.note(self,' %3d  %11.8f   %12.9f' %(it+1, self.e_tot, self.e_ele-_e_ele))
            if abs(_e_ele-self.e_ele)<self.conv_tol[1]: break
            # else: rgvb.dump_mo(self, '%s_gvb_init-II.fchk' %in_pre2b, reffile)
        else: logger.error(self, 'Not converge in local optimization')
        logger.note(self, 'The GVB pair:')
        for _i in range(self.nocc-self.npair_d, self.nocc):
            logger.note(self, '%2d(%6f) <-> %2d(%6f)'
                        %(self.pair[_i,0]+1, self.ci[_i,0], self.pair[_i,1]+1, self.ci[_i,1]))
        self.e_tot=self.e_nuc+self.e_ele
        logger.note(self, 'gvb local converged energy = %11.8f' %self.e_tot)
        rgvb.dump_mo(self, '%s_gvb_init-II.fchk' %in_pre2b, reffile)

        self.h_mo = reduce(numpy.dot, (self.mo_coeff.T, hcore, self.mo_coeff))
        self.g_mo = ao2mo.kernel(eri, self.mo_coeff)
        logger.note(self, 'start global optimization iteration')
        logger.note(self, '  it       energy         delta')
        logger.note(self, '  it  obt1 orb2       q4           q3           q2           q1         theta       energy       delta(E)        fun        delta(f)')
        for it in range(self.max_cycle[0]):
            _e_ele = self.e_ele
            for u in reversed(range(self.nocc)):
                for i in range(2):
                    if u<self.nocc-self.npair_d and i==1: continue
                    for v in reversed(range(self.nocc-self.npair_d-self.nsig, self.nocc)):
                        for j in range(2):
                            if v < self.nocc - self.npair_d and j == 1: # 注意：单占据的geminal，只能j=0，不能j=1,
                                continue # 因为默认单占据s=0，若j取到1时，也被设置为s=0，而vj=1是core轨道，双占据，就不对了
                            if u==v and i==j: continue
                            b1 = 8*self.ci[u,i]**2 * self.ci[v,j]**2 * (1.-delta(u,v))
                            b2 = (self.s(u)*self.ci[u, i]**2 + self.s(v)*self.ci[v, j]**2
                                       - 2*self.s(u)*self.ci[u, i]*self.ci[v, j]*delta(u, v)
                                       - 2 * self.ci[u, i]**2 * self.ci[v, j]**2 * (2-self.t(u,v)) * (1.-delta(u,v)))
                            auv = (b1*(self.g_mo[index(self.pair[u,i], self.pair[u,i], self.pair[v,j], self.pair[v,j])]
                                     - self.g_mo[index(self.pair[u,i], self.pair[v,j], self.pair[u,i], self.pair[v,j])])
                                     +2.*(self.w(u,i,v,j,v,j) - self.w(u,i,u,i,u,i) - self.w(v,j,v,j,v,j) + self.w(v,j,u,i,u,i)))
                            buv = self.w(u,i,u,i,v,j) - self.w(v,j,u,i,v,j)
                            q1 = 4.*buv
                            q2 = auv + 2.*b2*(self.g_mo[index(self.pair[u,i], self.pair[u,i], self.pair[v,j], self.pair[v,j])]
                                            + self.g_mo[index(self.pair[u,i], self.pair[v,j], self.pair[u,i], self.pair[v,j])])
                            q3 = 4.*(buv + b2*(self.g_mo[index(self.pair[u,i], self.pair[v,j], self.pair[v,j], self.pair[v,j])]
                                             - self.g_mo[index(self.pair[u,i], self.pair[u,i], self.pair[u,i], self.pair[v,j])]))
                            q4 = auv + b2*(self.g_mo[index(self.pair[u,i], self.pair[u,i], self.pair[u,i], self.pair[u,i])]
                                         + self.g_mo[index(self.pair[v,j], self.pair[v,j], self.pair[v,j], self.pair[v,j])]
                                       -2.*self.g_mo[index(self.pair[u,i], self.pair[v,j], self.pair[u,i], self.pair[v,j])])
                            # theta=root(q4,q3,q2,q1)%(math.pi*2)
                            theta=root(q4,q3,q2,q1)
                            # if abs(theta-0.0)>self.conv_tol[0]:
                            f=fun(q4, q3, q2, q1, theta); __e_ele=self.e_ele
                            orb1 = self.mo_coeff[:, self.pair[u, i]].copy()
                            orb2 = self.mo_coeff[:, self.pair[v, j]].copy()
                            self.mo_coeff[:, self.pair[u, i]] = orb1*math.cos(theta)+orb2*math.sin(theta)
                            self.mo_coeff[:, self.pair[v, j]] = -orb1*math.sin(theta)+orb2*math.cos(theta)
                            self.h_mo = reduce(numpy.dot, (self.mo_coeff.T, hcore, self.mo_coeff))
                            self.g_mo = ao2mo.kernel(eri, self.mo_coeff)
                            # eri = self.mol.intor('int2e_cart', aosym='s8')
                            # self.g_mo = ao2mo.incore.full(ao2mo.restore(8, eri, self.nbs), self.mo_coeff)
                            deltaf=self.energy_elec()-__e_ele
                            self.opt_ci()
                            self.e_ele=self.energy_elec()
                            self.e_tot=self.e_ele+self.e_nuc; deltaE=self.e_ele-__e_ele
                            # logger.note(self,' %3d %3d<->%3d %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f'
                            #   %(it+1,self.pair[u,i]+1,self.pair[v,j]+1,q4,q3,q2,q1,theta,self.e_tot,deltaE,f,deltaf))
                    for v in range(self.nocc, self.npair_t):
                        j = 0
                        b2 = self.s(u)*self.ci[u,i]**2
                        auv = self.w(u,i,v,j,v,j) - self.w(u,i,u,i,u,i)
                        buv = self.w(u,i,u,i,v,j)
                        q1 = 4.*buv
                        q2 = 2.*(auv + b2*(self.g_mo[index(self.pair[u,i], self.pair[u,i], self.pair[v,j], self.pair[v,j])]
                                         + self.g_mo[index(self.pair[u,i], self.pair[v,j], self.pair[u,i], self.pair[v,j])]))
                        q3 = 4.*(buv + b2*(self.g_mo[index(self.pair[u,i], self.pair[v,j], self.pair[v,j], self.pair[v,j])]
                                         - self.g_mo[index(self.pair[u,i], self.pair[u,i], self.pair[u,i], self.pair[v,j])]))
                        q4 = 2.*auv + b2*(self.g_mo[index(self.pair[u,i], self.pair[u,i], self.pair[u,i], self.pair[u,i])]
                                        + self.g_mo[index(self.pair[v,j], self.pair[v,j], self.pair[v,j], self.pair[v,j])]
                                      -2.*self.g_mo[index(self.pair[u,i], self.pair[v,j], self.pair[u,i], self.pair[v,j])])
                        # theta=root(q4,q3,q2,q1)%(math.pi*2)
                        theta=root(q4,q3,q2,q1)
                        # if abs(theta-0.0)>self.conv_tol[0]:
                        f=fun(q4, q3, q2, q1, theta); __e_ele=self.e_ele
                        orb1 = self.mo_coeff[:, self.pair[u, i]].copy()
                        orb2 = self.mo_coeff[:, self.pair[v, j]].copy()
                        self.mo_coeff[:, self.pair[u, i]] = orb1*math.cos(theta)+orb2*math.sin(theta)
                        self.mo_coeff[:, self.pair[v, j]] = -orb1*math.sin(theta)+orb2*math.cos(theta)
                        self.h_mo = reduce(numpy.dot, (self.mo_coeff.T, hcore, self.mo_coeff))
                        self.g_mo = ao2mo.kernel(eri, self.mo_coeff)
                        # eri = self.mol.intor('int2e_cart', aosym='s8')
                        # self.g_mo = ao2mo.incore.full(ao2mo.restore(8, eri, self.nbs), self.mo_coeff)
                        deltaf = self.energy_elec() - __e_ele
                        self.opt_ci()
                        self.e_ele=self.energy_elec()
                        self.e_tot=self.e_ele+self.e_nuc; deltaE=self.e_ele-__e_ele
                        # logger.note(self,' %3d %3d<->%3d %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f'
                        #   %(it+1,self.pair[u,i]+1,self.pair[v,j]+1,q4,q3,q2,q1,theta,self.e_tot,deltaE,f,deltaf))
            # for _i in range(self.nocc-self.npair_d,self.nocc):
            #     print('%2d(%6f) <-> %2d(%6f)' %(self.pair[_i, 0]+1, self.ci[_i, 0], self.pair[_i, 1]+1, self.ci[_i, 1]))
            logger.note(self, ' %3d  %11.8f   %12.9f' % (it + 1, self.e_tot, self.e_ele - _e_ele))
            if abs(_e_ele-self.e_ele)<self.conv_tol[0]: break
            # else: rgvb.dump_mo(self, '%s_gvb_init-III.fchk' %in_pre2b, reffile)
        else: raise RuntimeError('Not converge in global optimization')
        logger.note(self, 'The GVB pair:')
        for _i in range(self.nocc - self.npair_d, self.nocc):
            logger.note(self, '%2d(%6f) <-> %2d(%6f)'
                        %(self.pair[_i,0]+1, self.ci[_i,0], self.pair[_i,1]+1, self.ci[_i,1]))
        self.e_tot=self.e_nuc+self.e_ele
        logger.note(self, 'gvb global converged energy = %11.8f' %self.e_tot)
        rgvb.dump_mo(self, '%s_gvb_init-III.fchk' % in_pre2b, reffile)













if __name__ == "__main__":
    '''

    '''
    pass








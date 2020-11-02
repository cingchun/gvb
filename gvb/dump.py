#!/usr/bin/env python 
# -*- coding: utf-8 -*- 
#
#      Author:   Qingchun Wang @ NJU
#      E-Mail:   qingchun720@foxmail.com
#

'''
input a file, output a or more file

    dump infile outfile1 outfile1 opt1 opt2 ...
'''

# TODO (1) match by Regex (improt re)
# TODO (2) argv by getopt


import os, time, sys, getopt, subprocess, shutil, re
from pyscf import gto, lib


PI = 3.1415926535897932384626433832795028
bohr2ang = 0.52917724924; ang2bohr = 1.8897259877218677
# option type
shortargs='p'
longargs=['np','basis','sto-6g','sto-3g','6D','5D','10F','7F','PM','Boys','meta_Lowdin','Lowdin','Mulliken']
# set type
rhfs=set(['rhf','rohf', 'rb3lyp', 'rob3lyp']); hfs=rhfs|set(['uhf', 'ghf', 'ub3lyp'])
metds=hfs|set(['bccc', 'gvb', 'cas']); ntd={'nt1':1, 'nt2':2, 'nt3':3, 'nt4':4}
tasks=set(['sp','opt'])
auxf3='sto'; basisf3=set([auxf3,'gen','3-2', '6-3', 'cc-', 'def', 'aug']) # auxbasis and basis
polarizations=set(['*', '**']); diffusions=set(['+', '++'])
carts=set(['5d', '6d', '7f', '10f'])
pops=set(['dipole', 'mulliken', 'lowdin', 'meta_lowdin']); weights=set(['trans_dipo', 'distance','angle','exch_inte']); schs=set(['sch1','sch2'])
iexts=set(['.gjf','.com','.inp','.py']); oexts=set(['.log','.txt','.out', '.fchk']); exts=iexts|oexts
softs = set(['pyscf', 'gaussian', 'gamess']); iextd,oextd = {'pyscf':'.py', 'gaussian':'.gjf', 'gamess':'.inp'},{'pyscf':'.out', 'gaussian':'.log', 'gamess':'.txt'}  # defulte in/out for soft

infile=''; osofts=[]; inpre=''; outsuf=''
mol=gto.M(); metd='bccc'; np=None; nt=2; basis='cc-pvdz'; auxbasis='sto-6g'; cart=True; pop='meta_lowdin'; weight='trans_dipo'; sch='sch2'


def fsplit(path):
    '''
    filename split to (path, name, base, ext)
    
    :param path:

    :return: path, name, base, ext
    '''
    
    path,name = os.path.split(path)
    base,ext = os.path.splitext(name)

    return path,name,base,ext
def ncpu(type='physical_core'):
    if type=='physical_core':
        return (int(os.popen('cat /proc/cpuinfo| grep "physical id"| sort| uniq| wc -l').read().strip())
                * int(os.popen('cat /proc/cpuinfo| grep "cpu cores"| uniq').read().strip().split()[-1]))
    elif type=='physical_cpu':
        return int(os.popen('cat /proc/cpuinfo| grep "physical id"| sort| uniq| wc -l').read().strip())
    else: # logical core
        from multiprocessing import ncpu
        return ncpu()
def nmem(type='free'):
    if type=='free': return int(os.popen('cat /proc/meminfo').readlines()[1].strip().split()[1])//1232896
    elif type=='avail': return int(os.popen('cat /proc/meminfo').readlines()[2].strip().split()[1])//1232896
    else: return int(os.popen('cat /proc/meminfo').readlines()[0].strip().split()[1])//1232896

def run(command):
    print(command)
    (status, output) = subprocess.getstatusoutput(command)
    if status!=0: raise RuntimeError(output)
def cal(file_base, filetype):
    '''
    :param filetype: str. one of *.gjf(*.com), *.inp and *.py
    :param metd: str. one of uhf(uno), gvb, cas
    '''
    in_b=file_base; in_e=filetype; infile=in_b+in_e
    if in_e=='.gjf' or in_e=='.com':
        if os.path.exists(in_b+'.log'): print(F'Warning: {in_b}.log file already exists, PASS......')
        else: run(F'g16 {infile}')
    elif in_e=='.inp':
        if os.path.exists(in_b+'.txt'): print(F'Warning: {in_b}.txt file already exists, PASS......')
        else: run(F'rungms {infile} 01 {ncpu()} &> {in_b}.txt ')
    elif in_e=='.py':
        if os.path.exists(in_b+'.out'): print(F'Warning: {in_b}.out file already exists, PASS......')
        else: run(F'python {infile}')
    else: raise RuntimeError(F'{in_e} type is not recognized')
def extv(filename, format, pos=0):
    from gvb.scanf import scanf
    with open(filename, 'r') as fin:
        for line in fin.readlines():
            # words=line.strip(' \r\n').split()      # traditional mode
            # if len(words)>3:
            #     if words[0]=='SCF': return float(words[3])
            val = scanf(format, line)  # scanf mode
            if val: return val[pos]
            # import re                              # regex mode
            # regex=re.compile(format)
            # val=regex.search(format,line)
            # if val: return val.group()[pos]
        else: return .0 # raise IOError('Error: Eenery value can\'t found!')


def ncore(mol):
    '''
        return: core orbital number of a molecule in integer
    '''
    freeze=[0,2,10,18,36,54,86,118]
    n=0
    for i in mol.atom_charges():
        # for j in range(len(freeze)):
        #     if i<=freeze[j]:
        #         n+=(freeze[j-1]//2)
        #         break
        if i<=2: pass
        elif i<=10: n+=1
        elif i<=18: n+=5
        else: n+=9
        # else: n+=5
    return n
def nrefobt(mol, inpre, auxbasis='sto-6g'):
    mol = gto.copy(mol).build(basis=auxbasis, cart=False)
    nobt = mol.nao_nr()
    return nobt
def mkfchk(mol, inpre, metd):
    '''
    make .fchk file
    
    :param outfile: xxx.fchk
    :param mol:
    :param inpre: xxxx of xxxx.gjf
    
    :return:
    '''
    
    in2b = F'{inpre}_{mol.basis}' if isinstance(mol.basis, str) else F'{inpre}_gen'; outb=F'{in2b}_{metd}'
    outgjf,outchk,outfile = F'{outb}.gjf',F'{outb}.chk',F'{outb}.fchk'
    print(F'make {outfile} file')

    if not os.path.exists(outgjf): dumpfile(mol, inpre, 'gaussian', metd=metd)
    if os.path.exists(outfile):
        run(F'unfchk {outfile} {outchk}')
        cal(outb, '.gjf')
    else:
        cal(outb, '.gjf')
        run(F'formchk {outchk} {outfile}')
    run(F'rm -r {outchk}')
    run(F'mv fort.7 {in2b}.fort.7')
def datfort7(mol, inpre):
    '''
    dat from .fort.7 file
    
    :param infile: xxxx.fort.7
    :param mol:
    
    :return: nobt, dat
    '''
    
    in2b = F'{inpre}_{mol.basis}' if isinstance(mol.basis, str) else F'{inpre}_gen'
    outfchk,infile = F'{in2b}_uhf.fchk',F'{in2b}.fort.7'
    
    if not (os.path.exists(infile) and os.path.exists(outfchk)): mkfchk(mol, inpre, 'uhf')
    nobt = extv(outfchk, 'Number of independent functions            I %d')
    with open(infile, 'r') as fin:
        dat = fin.readlines()[1:]
    
    return nobt, dat


def setargv(argv):
    global infile, inpre, outsuf, metd, np, nt, basis, auxbasis, cart, pop, weight, sch

    def fixargv(argv):
        argv = {w.lower():w for w in argv}
        ios=[None]; bas=[]; opts=[]
        for arg in argv:
            if os.path.splitext(arg)[1] in iexts: ios[0] = argv[arg]  # file
            elif arg in softs: ios.append(arg)
            elif len(arg)>2 and arg[:3] in basisf3: # basis
                if arg[:3]==auxf3: bas.insert(0, arg)
                else: bas.append(arg)
            else: opts.append(arg) # opts
        return ios, bas, opts
    ios, bas, opts = fixargv(argv)  # fix command arg

    if ios[0] is None: raise IOError('Error: no input file')
    else: infile = ios[0]

    def loadargv(infile):
        '''
        load argv from infile
        
        :param infile: str, xxxx.gjf
        
        :return:
        '''
        
        # print('load argv from '+infile)
        with open(infile, 'r') as input:
            oext=os.path.splitext(infile)[1]
            if oext=='.gjf' or oext=='.com':
                while input.readline().strip('\r\n'): pass  # skip link and route section
                argv = infile; line = input.readline().strip('\r\n')
                while line:
                    argv+=' '+line
                    line = input.readline().strip('\r\n')
            else: raise NotImplementedError(F'TODO: load_opt() for {oext} type file')
        return argv.split()
    _,_bas,_opts = fixargv(loadargv(infile))  # fix arg in infile, and ios only in command

    if not bas: bas=_bas
    opts = _opts+opts

    # confirm arg
    if len(bas)==1:
        basis = bas[0]
        if basis[:3] == auxf3: auxbasis=basis
    elif len(bas)==2:
        auxbasis=bas[0]; basis=bas[1]
    else: raise RuntimeError(F'Error: too many basis {bas}')

    for arg in opts:
        if arg in metds: metd=arg
        elif arg in carts:
            if arg in ['5d', '7f']: cart=False
        elif arg in pops: pop=arg
        elif arg in weights: weight=arg
        elif arg in schs: sch=arg
        elif arg.isdigit(): np=int(arg)
        elif arg in ntd: nt=ntd[arg]
        else:
            # raise NotImplementedError('TODO: set_opt() for %s keyword' %arg)
            pass
    
    inpre,_ = os.path.splitext(infile)
    if ios[1:]: osofts.extend(ios[1:])
    else: osofts.append('pyscf')

    def loadmol(infile):
        '''
        load mol from infile
        
        :param infile: str, xxxx.gjf
        
        :return:
        '''
        
        def load_gen(infile):
            print(F'load gen from {infile}')
            global gen
            with open(infile, 'r') as input:
                _,ine = os.path.splitext()
                if ine=='.gjf' or ine=='.com':
                    while input.readline().strip('\r\n'): pass  # skip link and route section
                    while input.readline().strip('\r\n'): pass  # title section
                    while input.readline().strip('\r\n'): pass  # geometry section
                    gen=''; line=input.readline().strip('\r')
                    while line!='\n':
                        gen+=line
                        line = input.readline().strip('\r')
                else: raise NotImplementedError(F'TODO: load_opt() for {ine} type file')
        def gen2dict(gen):
            '''
            gen basis in g16 -> parse in pyscf
            :param gen: str
            :return: dict
            '''
            basis={}
            for i in range(len(gen)//3):
                pbasis=3*i+1; elems=gen[pbasis-1].split()[:-1]
                for elem in elems:
                    basis[elem]=gen[pbasis]
            return basis
        with open(infile, 'r') as input:
            _,ine = os.path.splitext(infile)
            if ine=='.gjf' or ine=='.com':
                connectivitys=[]; line=input.readline().strip('\r\n') # skip link and route section
                while line:
                    connectivitys+=[line]
                    line = input.readline().strip('\r\n')
                while input.readline().strip('\r\n'): pass  # skip title section
                # import re # simple mode
                # args_=re.split('[,/()= ]', args.lower())
                # args=re.split('[,/= ]', args.strip()); basis=args[0]
                # if '6d' in args_ or '10f' in args_: cart = True
                # if 'boys' in args_: local = 'Boys'
                # if 'lowdin' in args_: pop = 'Lowdin'
                # if 'mulliken' in args_: pop = 'Mulliken'
                # if 'sto-3g' in args_: auxbasis='STO-3G'
                # # args = 'basis='; line = input.readline().strip('\r\n')    # getopt mode
                # # while line:
                # #     args += ' ' + line
                # #     line = input.readline().strip('\r\n')
                # # args=args.split(); for arg in args: arg='--'+arg
                # # try: opts, args = getopt.getopt(args,shortargs,longargs)
                # # except getopt.GetoptError as err:
                # #     raise RuntimeError('Some options is not recognized')
                # # opts_k=[]; for k, v in opts: opts_k+=[k.strip('_')]
                line=input.readline().strip('\r\n'); words=line.split()
                charge=int(words[0]); spin=int(words[1])
                geom=''; line=input.readline().strip('\r')
                while line!='\n':
                    geom+=line
                    line = input.readline().strip('\r')

                if basis!='gen': mol.build(atom=geom, charge=charge, spin=spin-1, basis=basis, cart=cart, verbose=0)
                else:
                    findconnectivity=False
                    for line in connectivitys:
                        if line.find('geom')>=0:
                            findconnectivity=True
                            break
                    if findconnectivity:
                        while input.readline().strip('\r\n'): pass  # skip connectivity
                    gen=[]; line=input.readline().strip('\r\n')
                    while line:
                        gen+=[line]
                        line = input.readline().strip('\r\n')
                    mol.build(atom=geom, charge=charge, spin=spin-1, basis=gen2dict(gen), cart=cart, verbose=0)
            else: raise NotImplementedError(F'TODO: load_opt() for {ine} type file')
    loadmol(infile)
def dumpfile(mol,  inpre='', osoft='pyscf', metd='bccc', np=None, auxbasis='sto-6g', pop='dipole', weight='trans_dipo', sch='sch2', outsuf=''):
    '''
    dump input file for calculation
    
    :param osoft:  str( xxx_bond_basis_method_step.ext )
    :param mol:
    :param metd:
    :param inpre:
    :param auxbasis:
    :param pop:
    :param weight:
    
    :return:
    '''
    
    in_pre2b = F'{inpre}_{mol.basis}' if isinstance(mol.basis, str) else inpre+'_gen'
    in_pre2m = F'{in_pre2b}_{metd}'
    outfile = F"{in_pre2m}{'_'+outsuf if outsuf else ''}{iextd[osoft]}"
    out_b, out_e=os.path.splitext(outfile)
    print(F'dump {outfile}')

    if os.path.exists(outfile): print(F'Warning: {outfile} file already exists! PASS......')
    else:
        # cas=(0, 0); nco=0; nelec=(0,0); nobt=0; dat=''
        def det_npair(mol, np, inpre, auxbasis):
            '''
            determine np, nco, nelec
            
            :param mol: gto.Mole
            :param auxbasis: sto-3g or sto-6g
            :param inpre: xxxx.gjf
            :return: np, nco, nelec
            '''
            
            nelec=mol.nelec
            nco=ncore(mol)  # 本身 core

            nobt=nrefobt(mol, inpre, auxbasis=auxbasis)
            _np = min(nelec[1]-nco, nobt-nelec[0])  # autogvb
            
            if np:
                if np>_np: raise RuntimeError(F'np({np}) > min(nactocc,nactvir)({_np})')
            else: np = _np
            nco = nelec[0]-np  # include 孤对
            cas = (2*np, 2*np)
            
            return np, cas, nco, nelec
        def det_nobt(mol, inpre):
            '''
            :param mol: gto.Mole
            :param inpre: xxxx.gjf
            :return: nobt
            '''
            # mf = scf.remove_linear_dep_(scf.UHF(mol), threshold=1e-6).run()
            # nobt = mf.mo_coeff[0].shape[1]

            if isinstance(mol.basis, str): in_pre2m='_'.join([inpre, mol.basis, 'uhf'])
            else: in_pre2m=inpre+'_gen_uhf'
            if not os.path.exists(in_pre2m+'.fchk'): mkfchk(in_pre2m+'.gjf', mol, 'uhf', inpre)
            nobt=extv(in_pre2m+'.fchk', 'Number of independent functions            I %d')

            return nobt
        def dict2gen(dict):
            '''
            parse in pyscf -> gen basis in g16
            :param dict:
            :return: str
            '''
            gen=''
            for bas in set(dict.values()):
                elms=[]
                for elm in dict.keys():
                    if dict[elm]==bas: elms += [elm]
                gen += '\n'.join([' '.join(elms)+' 0', bas, '****'])+'\n'
            return gen
        np, cas, nco, nelec = det_npair(mol, np, inpre, auxbasis=auxbasis)
        if out_e == '.inp':
            # if isinstance(mol.basis, str): nobt, dat=datfort7(mol, b inpre+'_'+mol.basis+'.fort.7', mol)
            # else: nobt, dat=datfort7(inpre+'_gen.fort.7', mol)
            # # nobt = det_nobt(mol, inpre)
            nobt,dat = datfort7(mol,inpre)

        with open(outfile, 'w') as fout:
            if out_e=='.gjf' or out_e=='.com':
                fout.write(F'%nprocshared={ncpu()}\n')
                fout.write(F'%mem={nmem()}GB\n')
                fout.write(F'%chk={out_b}.chk\n')
                fout.write('#p int=(nobasistransform) nosymm pop=full punch=gamessinput gfprint\n')

                if metd == 'uhf': _metd='uhf'
                elif metd == 'rhf':
                    if mol.spin: _metd='rohf'
                    else: _metd='rhf'
                elif metd == 'gvb': _metd=F'gvb({np})'
                elif metd == 'cas': _metd=F'cas({cas[0]}, {cas[1]})'
                else: raise ModuleNotFoundError(F"{metd} module can't found")

                if isinstance(mol.basis, str): _basis=mol.basis
                else: _basis='gen'
                if mol.cart: _basis += ' 6D 10F'
                else: _basis += ' 5D 7F'

                if metd=='uhf':
                    if os.path.exists(out_b+'.fchk'):
                        # run('unfchk %s %s' %(out_b+'.fchk',out_b+'.chk'))
                        fout.write(F'   {_metd}/{_basis} guess=read SCF=(xqc, maxcycle=500) stable=opt\n\n')
                    elif os.path.exists(out_b+'.chk'):
                        fout.write(F'   {_metd}/{_basis} guess=read SCF=(xqc, maxcycle=500) stable=opt\n\n')
                    else:
                        fout.write(F'   {_metd}/{_basis} guess=mix SCF=(xqc, maxcycle=500) stable=opt\n\n')
                elif metd=='rhf':
                    if mol.spin==0:
                        if os.path.exists(out_b+'.fchk'):
                            # run('unfchk %s %s' %(out_b+'.fchk',out_b+'.chk'))
                            fout.write(F'   {_metd}/{_basis} guess=read SCF=(xqc, maxcycle=500)\n\n')
                        elif os.path.exists(out_b+'.chk'):
                            fout.write(F'   {_metd}/{_basis} guess=read SCF=(xqc, maxcycle=500)\n\n')
                        else:
                            fout.write(F'   {_metd}/{_basis} SCF=(xqc, maxcycle=500)\n\n')
                    else:
                        if os.path.exists(out_b+'.fchk'):
                            # run('unfchk %s %s' %(out_b+'.fchk',out_b+'.chk'))
                            fout.write(F'   {_metd}/{_basis} guess=read SCF=(maxcycle=500)\n\n')
                        elif os.path.exists(out_b+'.chk'):
                            fout.write(F'   {_metd}/{_basis} guess=read SCF=(maxcycle=500)\n\n')
                        else:
                            fout.write(F'   {_metd}/{_basis} SCF=(maxcycle=500)\n\n')
                elif metd=='gvb' or metd=='cas':
                    if os.path.exists(out_b+'.fchk'):
                        # run('unfchk %s %s' %(out_b+'.fchk',out_b+'.chk'))
                        fout.write(F'   {_metd}/{_basis} guess=read SCF=(maxcycle=500)\n\n')
                    elif os.path.exists(out_b+'.chk'):
                        fout.write(F'   {_metd}/{_basis} guess=read SCF=(maxcycle=500)\n\n')
                    else:
                        fout.write(F'   {_metd}/{_basis} SCF=(maxcycle=500)\n\n')
                else: raise ModuleNotFoundError(F"{metd} module can't found")
                fout.write( 'Title section\n\n')
                fout.write(F'{mol.charge} {mol.spin+1}\n')
                fout.writelines(mol.atom.lstrip('\n'))

                if isinstance(mol.basis, dict):
                    fout.write('\n')
                    fout.write(dict2gen(mol.basis))

                if metd=='gvb':
                    fout.write('\n')
                    l2 = []; i = 0
                    for i in range(1, np//30+1): l2.append(''.join([' 2']*30)+'\n')
                    else: l2.append(''.join([' 2'] * (np-30*i))+'\n')
                    fout.writelines(l2)
                # elif metd=='uhf':
                #     fout.write( '\n')
                #     fout.write( ' --link1--\n')
                #     fout.write(F'%nprocshared={ncpu()}\n')
                #     fout.write(F'%mem={nmem()}GB\n')
                #     fout.write(F'%chk={out_b}.chk\n')
                #     fout.write( '#p int=(nobasistransform) nosymm pop=full punch=gamessinput gfprint\n')
                #     fout.write(F'   {_metd} chkbasis geom=allcheck guess=(read,NaturalOrbitals,save,only)\n')
                else: pass
                # if metd in hfs:
                #     fout.write( '\n')
                #     fout.write( ' --link1--\n')
                #     fout.write(F'%nprocshared={ncpu()}\n')
                #     fout.write(F'%mem={nmem()}GB\n')
                #     fout.write(F'%chk={out_b}.chk\n')
                #     fout.write( '#p int=(nobasistransform) nosymm pop=full punch=gamessinput gfprint\n')
                #     fout.write(F'   {_metd} chkbasis geom=allcheck guess=(read,local,save,only)\n')
                fout.write( '\n\n')
            elif out_e=='.inp':
                if metd=='gvb':
                    fout.write(F' $CONTRL SCFTYP=gvb RUNTYP=ENERGY ICHARG={mol.charge:2d} MULT={mol.spin+1} UNITS=BOHR\n')
                    fout.write(F'         maxit=200 ispher=-1 nosym=1 $END\n')
                    fout.write(F' $SYSTEM Mwords={nmem()*125} $END\n')
                    if nelec[0]==nelec[1]: fout.write(F' $SCF Nco={nco} Npair={np} $END\n')
                    else: fout.write(F" $SCF Nco={nco} Npair={np} NSETO={nelec[0]-nelec[1]} NO(1)={','.join(['1']*(nelec[0]-nelec[1]))} $END\n")
                    fout.write(F' $GUESS guess=moread Norb={nobt} $END\n')
                    fout.writelines(dat)
                else: raise NotImplementedError(F'TODO: dump .inp file for {metd}')
            elif out_e=='.py':
                fout.write( '#!/usr/bin/env python                              \n')
                fout.write( '# -*- coding: utf-8 -*-                            \n')
                fout.write( '#                                                  \n')
                fout.write( '#         Author:  Wingchun Wang @ NJU             \n')
                fout.write( '#         E-Mail:  qingchun720@foxmail.com         \n')
                fout.write( '#                                                  \n')
                fout.write( '\n\n')

                fout.write( 'import os, sys, time\n')
                fout.write(F"fout=open('{out_b}.out', 'w')\n")
                fout.write( 'stdout=sys.stdout; stderr=sys.stderr; sys.stdout=fout; sys.stderr=fout\n')
                fout.write( "print(time.strftime('Start abic program at %H:%M:%S %m/%d %Y', time.localtime()))\n")
                fout.write( "print(' Author:  Qingchun Wang @ NJU')\n")
                fout.write( "print(' E-Mail:  qingchun720@foxmail.com')\n")
                fout.write(F"print(' input  file: {out_b}.py')\n")
                fout.write(F"print(' output file: {out_b}.out')\n")
                fout.write( "print('\\n')\n")
                fout.write( '\n\n')

                fout.write( 'from pyscf import gto\n')
                fout.write( 'mol = gto.M(\n')
                fout.write( '    atom=\n')
                fout.write( "        '''\n")
                fout.writelines(mol.atom)
                fout.write( "        ''',\n")
                fout.write(F'    charge={mol.charge}, spin={mol.charge},\n')
                fout.write(F"    output='{out_b}.out',\n")
                fout.write(F'    cart={mol.cart},\n')
                fout.write( '    verbose=4,\n')
                if isinstance(mol.basis, str): fout.write(F"    basis='{mol.basis}', \n")
                else: fout.write(F'    basis={str(mol.basis)}\n')
                fout.write( '           )\n')
                if metd=='rhf':
                    fout.write('from pyscf import scf\n')
                    fout.write('from gvb.rgvb import parse, hf, dump_ndarr, load_mo\n')
                    fout.write('basis = mol.basis\n')
                    fout.write('mol.basis = parse(mol)\n')
                    fout.write("mf = hf(mol, basis, mol.cart, 'rhf')\n")
                    fout.write('#mf = scf.RHF(mol)\n')
                    fout.write('#mf = scf.remove_linear_dep_(mf, threshold=1e-6)\n')
                    fout.write('#mf.kernel()\n')
                elif metd=='uhf':
                    fout.write('from pyscf import scf\n')
                    fout.write('from gvb.rgvb import parse, hf, dump_ndarr, load_mo\n')
                    fout.write('basis = mol.basis\n')
                    fout.write('mol.basis = parse(mol)\n')
                    fout.write("mf = hf(mol, basis, mol.cart, 'uhf')\n")
                    fout.write('#mf = scf.UHF(mol)\n')
                    fout.write('#mf = scf.remove_linear_dep_(mf, threshold=1e-6)\n')
                    fout.write('#mf.kernel()\n')
                elif metd == 'gvb':
                    fout.write('import gvb\n')
                    if mol.spin+1 is 1:
                        fout.write( 'mf=gvb.RGVB(mol)\n')
                        fout.write(F"mf.pop='{pop}'\n")
                        fout.write(F"mf.weight='{weight}'\n")
                        # fout.write(F"mf.sch='{sch}'\n")
                        # fout.write(F"mf.auxbasis='{auxbasis}'\n")
                    else:
                        fout.write(F"mf=gvb.UGVB(mol, pop='{pop}', weight='{weight}', auxbasis='{auxbasis}', sch='{sch}')\n")
                    fout.write('mf.kernel()\n')
                elif metd=='bccc':
                    fout.write('import gvb\n')
                    if mol.spin+1 is 1:
                        fout.write( 'mf=gvb.RGVB(mol)\n')
                        fout.write(F"mf.pop='{pop}'\n")
                        fout.write(F"mf.weight='{weight}'\n")
                        # fout.write(F"mf.sch='{sch}'\n")
                        # fout.write(F"mf.auxbasis='{auxbasis}'\n")
                    else:
                        fout.write(F"mf=gvb.UGVB(mol, pop='{pop}', weight='{weight}', auxbasis='{auxbasis}', sch='{sch}')\n")
                    fout.write( 'mf.kernel()\n')
                    fout.write( "print('\\n')\n")
                    fout.write( '\n\n')
                    fout.write( 'import bccc\n')
                    fout.write(F'cc = bccc.BCCC(mf, nt={nt})\n')
                    # fout.write( 'cc.xE = 0.015\n')
                    # fout.write( 'cc.ncyc = 500\n')
                    fout.write( 'cc.kernel()\n')
                elif metd=='cas':
                    fout.write( 'from pyscf import mcscf\n')
                    fout.write( '# import gvb\n')
                    fout.write( '# basis = mol.basis\n')
                    fout.write( '# mol.basis = gvb.parse(mol)\n')
                    fout.write( "# gvb.hf(mol, basis, mol.cart, 'rhf')\n")
                    
                    fout.write( 'mf=scf.RHF(mol)\n')
                    fout.write( 'mf=scf.remove_linear_dep_(mf, threshold=1e-6)\n')
                    fout.write( 'mf.kernel()\n')
                    fout.write( "print('\\n')\n")
                    fout.write( '\n\n')

                    fout.write( "print(time.strftime('Entry CASSCF calculation at %H:%M:%S %m/%d %Y', time.localtime()))\n")
                    fout.write(F"gvb.rgvb.load_mo(mf, '{in_pre2b}_gvb_init.dat.fchk')\n")
                    fout.write( '#gvb.rgvb.swapobt(mf.mo_coeff, 11-1, 33-1)\n')
                    fout.write( '#gvb.rgvb.swapobt(mf.mo_coeff, 13-1, 34-1)\n')
                    fout.write(F'mc = mcscf.CASSCF(mf, {np*2}, {np*2})\n')
                    fout.write( '#mo = mc.sort_mo([17,18,45,46,55,56,57,58])\n')
                    fout.write( '#mc.fix_spin_(ss=0)\n')
                    fout.write( '#mc.fcisolver.level_shift = 0.2\n')
                    fout.write( '#mc.fcisolver.pspace_size = 1200\n')
                    fout.write( '#mc.fcisolver.max_space = 36\n')
                    fout.write( '#mc.max_cycle=100\n')
                    fout.write( 'mc.natorb=True\n')
                    fout.write( 'mc.kernel()\n')
                    fout.write(F"gvb.rgvb.dump_mo(mc, '{out_b}.fchk', '{in_pre2b}_rhf.fchk')\n")
                    fout.write( "print('\\n')\n")
                    fout.write( '\n\n')

                    fout.write( "print(time.strftime('Entry NEVPT2 calculation at %H:%M:%S %m/%d %Y', time.localtime()))\n")
                    fout.write( 'from pyscf.mrpt import nevpt2\n')
                    fout.write( 'nevpt2.NEVPT(mc).kernel()\n')
                else: raise ModuleNotFoundError(F"{metd} module can't found\n")
                fout.write("print('\\n')\n")
                fout.write('\n\n')

                fout.write("print(time.strftime('End abic program at %H:%M:%S %m/%d %Y', time.localtime()))\n")
                fout.write('fout.close()\n')
                fout.write('sys.stdout=stdout; sys.stderr=stderr\n')
            else: raise FileNotFoundError(F"{out_e} file is't recognized\n")

    return out_b






if __name__ == '__main__':
    setargv(sys.argv[1:])
    
    for osoft in osofts:
        dumpfile(mol, inpre=inpre, osoft=osoft, metd=metd, np=np, auxbasis=auxbasis, pop=pop, weight=weight, sch=sch, outsuf=outsuf)


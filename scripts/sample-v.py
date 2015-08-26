import os
import tempfile
import shutil
import subprocess
from os.path import abspath, join
import numpy as np
import argparse
from jinja2 import Template
from collections import OrderedDict

STARFLEET = '*** PSI4 exiting successfully. Buy a developer a beer!'


template1 = Template('''{{comment}}
plugin_load('{{wavekernel}}')

molecule {
  symmetry c1
{{zmatrix}}
}

set {
  basis 6-31g*
  reference uhf
  dft_grid_name sg1
}

set wavekernel {
  filename_out vectors.npy
  mode sample_v
  num_sample_descriptors 1000
  curve mu_prime
  subtract_sad true
  curve_min -0.5
  curve_max 0.5
  temp 1000.0
  num_curve 20
  print 0
}

energy('scf')
plugin('{{wavekernel}}')

''')


def main():
    p = argparse.ArgumentParser()
    p.add_argument('qcin')
    p.add_argument('-nt', default=16, type=int)
    p.add_argument('-wavekernel')
    p.add_argument('outfile')
    p.add_argument('vectors')

    args = p.parse_args()
    qc = load_qchem(args.qcin)

    run_psi4(qc, vectors=args.vectors, outfile=args.outfile,
             wavekernel=args.wavekernel, nt=args.nt)


def load_qchem(fn):
    k = None
    data = QChem()

    with open(fn) as f:
        for line in f:
            line = line.strip()
            if line == '$end':
                k = None
            elif line.startswith('$'):
                k = line[1:]
                data[k] = []
            elif k in data:
                data[k].append(line)
    return data


class QChem(OrderedDict):
    def __str__(self):
        lines = []
        for k, v in self.items():
            lines.append('$%s' % k)
            lines.extend([l for l in v if len(l) > 0])
            lines.append('$end')
            lines.append(os.linesep)
        return os.linesep.join(lines)

    def zmatrix(self):
        zmatrix = []
        for l in self['molecule']:
            zmatrix.append('  ' + ' '.join(l.split()[:7]))
        return '\n'.join(zmatrix)

    def topsi4(self, wavekernel):
        comment = '\n'.join('# %s' % l for l in self.get('comment', []))

        return template1.render(
            comment=comment, zmatrix=self.zmatrix(), wavekernel=wavekernel)


def run_psi4(qc, vectors=None, outfile=None, wavekernel=None, nt=16):
    curdir = abspath(os.curdir)

    try:
        tempdir = tempfile.mkdtemp()
        if outfile is None:
            outfile = join(tempdir, 'psi4.out')
        if wavekernel is not None:
            wavekernel = abspath(wavekernel)
        if vectors is not None:
            vectors = abspath(vectors)

        outfile = abspath(outfile)
        psi4data = qc.topsi4(wavekernel)

        os.chdir(tempdir)
        with tempfile.NamedTemporaryFile('w', suffix='.dat') as f:
            f.write(psi4data)
            f.flush()
            cmd = ['psi4', '-n', str(nt), '-i', f.name, '-o', outfile]

            p = subprocess.Popen(cmd, stderr=subprocess.PIPE,
                                 stdout=subprocess.PIPE)
            stdout, stderr = p.communicate()

            tail = subprocess.check_output(['tail', outfile]).decode('utf-8')
            if STARFLEET not in tail.splitlines():
                raise ValueError('PSI4 Exit Failure\n\n%s' % tail)

            if p.returncode:
                raise CalledProcessError(p.returncode, cmd)

            if vectors is not None:
                shutil.copy('vectors.npy', vectors)
            with open(outfile, 'r') as f:
                return f.read(), np.load('vectors.npy')

    finally:
        os.chdir(curdir)
        shutil.rmtree(tempdir)


if __name__ == '__main__':
    main()

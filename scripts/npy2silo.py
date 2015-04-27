import os
import argparse
import numpy as np
from pyvisfile import silo


def main():
    p = argparse.ArgumentParser()
    p.add_argument('input_npy')
    p.add_argument('output_silo')
    args = p.parse_args()

    inp = np.load(args.input_npy)
    if inp.ndim != 2 or inp.shape[1] != 4:
        p.error('Wrong shape: %s' % p.input_npy)

    mesh = np.asarray(inp[:, :3].T, order='c')
    assignment = np.asarray(inp[:, 3], order='c',)

    if os.path.exists(args.output_silo):
        p.error('File exists: %s' % args.output_silo)
    with silo.SiloFile(args.output_silo) as f:
        f.put_pointmesh('mesh', mesh)
        f.put_pointvar1('var', 'mesh', assignment)


if __name__ == '__main__':
    main()


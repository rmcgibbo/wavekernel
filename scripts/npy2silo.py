import os
import argparse
import numpy as np
from pyvisfile import silo


def main():
    p = argparse.ArgumentParser()
    p.add_argument('input_npz')
    p.add_argument('output_silo')
    args = p.parse_args()

    inp = np.load(args.input_npz)
    
    mesh = np.asarray(inp['coords'].T, order='C')
    v = np.asarray(inp['v'], order='C')
    v = np.clip(v, -5, 10)
    
    if os.path.exists(args.output_silo):
        p.error('File exists: %s' % args.output_silo)
    with silo.SiloFile(args.output_silo) as f:
        f.put_pointmesh('mesh', mesh)
        for i in range(20):
            f.put_pointvar1('var_%d' % i, 'mesh', np.ascontiguousarray(v[i, :]))



if __name__ == '__main__':
    main()


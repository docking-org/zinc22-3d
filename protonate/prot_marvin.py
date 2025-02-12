#!/bin/env python
"""
prot_marvin.py

Bascially just a wrapper to run cxcalc
Outputs all prot/taut species present at
more than %DIST at pH

Usage: 
prot_marvin.py <in.smi> <out.smi> <dist> <ph>

Created March 2013 Dahlia Weiss
"""

from __future__ import with_statement, print_function

from contextlib import contextmanager
from itertools import chain
from tempfile import NamedTemporaryFile

import os
import sys
import subprocess
from optparse import OptionParser

import logging
import sh

from past.builtins import basestring

def _command(envvar, default=None):
    path = os.environ.get(envvar, default)
    if path is None:
        sys.stderr.write("Please make sure {0} is set or that {1} is in your path\n".format(envvar, default))
        sys.exit(-1)
    toks = path.split()
    path = toks[0]
    args = toks[1:]
    try:
        cmd = sh.Command(path)
    except sh.CommandNotFound:
        sys.stderr.write("Please make sure {0} is set or that {1} is in your path\n".format(envvar, default))
        sys.exit(-1)
    if len(args) > 0:
        cmd = cmd.bake(*args)
    return cmd


_cxcalc = _command('CXCALCEXE', 'cxcalc')
_molconvert = _command('MOLCONVERTEXE', 'molconvert')
#basestring = " "


def as_handle(src, mode, opener=open, default=None, *args, **kwargs):
    if src is None:
        src = default
    if src is None:
        raise ValueError("Cannot use None as handle")
    if isinstance(src, basestring):
        return opener(src, mode, *args, **kwargs)
    else:
        return src


def smiles_file(src):
    return (line.split()[0:2] for line in src)


def _process_cxcalc_distribution_results(results, format='smiles'):
    if format != 'sdf':
        results = _molconvert(format, T='dist', _iter=True, _in=results)
    if format == 'smiles':
        results = (line.split() for line in results \
                    if not line.startswith('#SMILES'))
    for result in results:
        if 2 == len(result):  # make sure dist made it through
            smiles, dist = result  # unpack here if possible
            yield smiles, float(dist)
        else:  # ignore error for now, molecule won't be built
            pass


def generate_tautomers(smiles, pH=7, mol_format='smiles', cutoff=20):
    taut = _cxcalc.dominanttautomerdistribution.bake('"'+smiles+'"',
                                                     H=pH, 
                                                     C='false',
                                                     t='dist',
                                                     _piped=True)
    logging.info("Running {!s}".format(taut))
    results = taut()
    if (results == ''):
        sys.stdout.write(smiles+" failed to generate tautomer: _cxcalc.dominanttautomerdistribution.bake(smiles, H=pH, C='false', t='dist', _piped=True)\n")
        exit()

    processed = _process_cxcalc_distribution_results(results, format=mol_format)
    for smiles, dist in processed:
        if dist > cutoff:
            logging.debug("Took Tautomer")
            yield smiles, dist
        else:
            logging.debug("Skipped Tautomer")


def generate_protomers(smiles, pH, mol_format='smiles', cutoff=20):
    prot = _cxcalc.microspeciesdistribution.bake('"'+smiles+'"',
                                                 H=pH, 
                                                 t='dist',
                                                 _piped=True)
    logging.info("Running {!s}".format(prot))
    results = prot()
    processed = _process_cxcalc_distribution_results(results, format=mol_format)
    for smiles, dist in processed:
        if dist > cutoff:
            logging.debug("Picked Protomer")
            yield smiles, dist
        else:
            logging.debug("Skipped Protomer")


def cxcalc_single(smiles, cid='query', pH=7, mol_format='smiles', tautomer_cutoff=20, protomer_cutoff=20):
    seen = set()
    tautomers = generate_tautomers(smiles, pH=pH, 
                                           mol_format='smiles', 
                                           cutoff=tautomer_cutoff)
    for i, (taut_smiles, taut_dist) in enumerate(tautomers):
        protomers = generate_protomers(taut_smiles, pH=pH, 
                                                    mol_format=mol_format, 
                                                    cutoff=protomer_cutoff)
        for j, (prot_smiles, prot_dist) in enumerate(protomers):
            if prot_smiles not in seen:
                yield prot_smiles, cid, prot_dist
                seen.add(prot_smiles)


def cxcalc(compounds, **kwargs):
    compounds = (items[:2] for items in compounds if len(items) >= 2)
    for smi, cid in compounds:
        results = cxcalc_single(smi, cid, **kwargs)
        for result in results:
            yield result


def run_cxcalc(infile, outfile, pH, threshold, output_format='smiles'): 
    threshold = float(threshold)
    with as_handle(infile, 'r', default=sys.stdin) as in_f:
        compounds = smiles_file(in_f)
        calculated = cxcalc(compounds, pH=pH, mol_format=output_format)
        with as_handle(outfile, 'w', default=sys.stdout) as out_f:
            for compound in calculated:
                if compound[2] > threshold:
                    print("\t".join(map(str, compound)), file=out_f)

def main(argv) :
    """Parse arguments."""
    description = "A wrapper to run marvin cxcalc on a list of smiles"
    usage = "%prog [options]"
    version = "%prog: version 201303 - created by Dahlia Weiss"
    parser = OptionParser(usage=usage, description=description,
                          version=version)
    parser.set_defaults(infile=None, outfile=None, pH=7.4, threshold=20)
    parser.add_option("-i", "--infile", 
                      help="input file (default: stdin)")
    parser.add_option("-o", "--outfile", 
                      help="output file (default: stdout)")
    parser.add_option("-H", "--pH", type=float, 
                      help="pH (default: %default)")
    parser.add_option("-d", "--threshold", type=float, 
                      help="distribution threshold  (default: %default)")
    options, args = parser.parse_args(args=argv[1:])
    if len(args):
        parser.error("program takes no positional arguments.\n" +
                     "  Use --help for more information.")
    flag = run_cxcalc(infile=options.infile, 
                      outfile=options.outfile,
                      pH=options.pH,
                      threshold=options.threshold)
    return flag

if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG)
    sys.exit(main(sys.argv))


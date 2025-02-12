#!/usr/bin/env python

from __future__ import print_function

import argparse
import heapq
from collections import namedtuple
import itertools
try:
    from itertools import izip as zip
    from itertools import izip_longest as zip_longest
    from itertools import ifilter as filter
except ImportError: # will be 3.x series
    from itertools import zip_longest as zip_longest
except ImportError: # will be 3.x series
    pass
import sys

RawProtomer = namedtuple('RawProtomer', ['smiles', 'name', 'extra'])
DEFAULT_PROTOMER_TEMPLATE = "{0.smiles}\t{0.name}\t{1}"

DEFAULT_EXCLUDED_PATTERNS = []  # Allow everything

def default_smiles_processing(smiles):
    return smiles


def load_protomers(source):
    for line in source:
        line = line.strip()
        if line:
            tokens = line.split()
            smiles, name, extra = tokens[0], tokens[1], tokens[2:]
            protomer = RawProtomer(smiles=smiles, name=name, extra=extra)
            yield protomer
        

# Inspired by:
# http://wordaligned.org/articles/merging-sorted-streams-in-python.html
def merge(*sequences, **kwargs):
    """ Allows multiple sequences to be combined based on a key function
        :param: key: Get identity of item
        :param: wrap: Get practial value of item (useful when combined with groupby)
        :param: fill: Value to fill empty sequences with

        Example: 
            >>> xs = [[1, 1, 2, 3, 4, 4, 5], [1, 1, 2, 3, 4, 4, 5], [3, 3, 3, 4, 4], [2, 2, 4]]
            >>> groups = [itertools.groupby(x) for x in xs]
            >>> merged = merge(groups, key=lambda x: x[0], wrap=lambda x: tuple(x[1]), fill=())
            >>> list(merged)
            [((1, 1), (1, 1), (), ()), 
             ((2,), (2,), (), (2, 2)),
             ((3,), (3,), (3, 3, 3), ()),
             ((4, 4), (4, 4), (4, 4), (4,)),
             ((5,), (5,), (), ())]
    """
    _nothing = object()
    key = kwargs.get('key', lambda x: x)
    wrap = kwargs.get('wrap', lambda x: x)
    fill = kwargs.get('fill', None)

    def _buf_item(seq):
        it = iter(seq)
        el = next(it, _nothing)
        if el is _nothing:
            item = (_nothing, ())
        else:
            item = (key(el), itertools.chain([el], it))
        return item
    
    buf = [_buf_item(seq) for seq in sequences]

    while any(k is not _nothing for k, _it in buf):
        current_key = next(iter(k for k, _it in buf if k is not _nothing), _nothing)
        if current_key is _nothing:
            break
        selected = []
        for idx, item in enumerate(buf):
            el_key, it = item
            if el_key is _nothing or el_key != current_key:
                new_el, buf_item = fill, item
            else:
                new_el = wrap(next(it))
                buf_item = _buf_item(it)
            selected.append(new_el)
            buf[idx] = buf_item
        yield tuple(selected)
    
        

def filter_sources(sources, filter_patterns):
    no_forbidden_patterns = lambda protomer: all(pattern not in protomer.smiles for pattern in filter_patterns)
    for source in sources:
        yield filter(no_forbidden_patterns, source)


def limited_sources(sources, limits, default_limit=None):
    #source_limits = list(itertools.izip_longest(sources, limits, fillvalue=default_limit))
    source_limits = list(zip_longest(sources, limits, fillvalue=default_limit))
    for source, limit in source_limits:
        limited_source = itertools.islice(source, limit)
        yield limited_source


def extract_grouped_protomers(groups, names=None, process_smiles=default_smiles_processing):
    if names is None:
        names = itertools.count()

    #named_groups = itertools.izip(names, groups)
    named_groups = zip(names, groups)
    already_seen = set()
    for group_name, group in named_groups:
        for protomer in group:
            key = protomer.smiles
            if key not in already_seen:
                already_seen.add(key)
                yield protomer, group_name


def write_protomers(protomers, output, make_line=DEFAULT_PROTOMER_TEMPLATE.format):
    for protomer in protomers:
        line = make_line(*protomer)
        print(line, file=output)


def separate_compounds(groups, sort=True):
    key = lambda g: g[0]
    # Sort protomers for each key by the value of the last extra column (as a float)
    if sort:
        wrap = lambda g: tuple(sorted(g[1], reverse=True, key=lambda p: float(p.extra[-1])))
    else:
        wrap = lambda g: tuple(g[1])
    fill = ()

    loaded = [load_protomers(group) for group in groups]
    grouped = [itertools.groupby(group, key=lambda p: p.name) for group in loaded]
    separated = merge(*grouped, key=key, wrap=wrap, fill=fill)
    return separated


def group_protomers(groups, output, group_limits=None, group_names=None, filter_patterns=[]):
    filtered_groups = filter_sources(groups, filter_patterns)
    limited_groups = limited_sources(filtered_groups, group_limits)
    grouped_protomers = extract_grouped_protomers(limited_groups, names=group_names)
    write_protomers(grouped_protomers, output=output)


def main(args, stdout=sys.stdout):
    parser = argparse.ArgumentParser()
    comma_separated_list = lambda line: line.split(',')
    comma_separated_int_list = lambda line: map(int, comma_separated_list(line))
    parser.add_argument('groups', type=argparse.FileType('r'),
                                  nargs='+',
                                  help="One or more protomer files to extract protomers from")
    parser.add_argument('-l', '--limits', type=comma_separated_int_list,
                                          nargs='?',
                                          default=[1],
                                          help="Comma separated list of maximum items to pull from. "
                                               "By deafult, one will be selected from the first group "
                                               "and the remaining groups will be exhaustively read.")
    parser.add_argument('-n', '--names', type=comma_separated_list,
                                         nargs='?',
                                         default=None,
                                         help="Comma separated list of group named. "
                                              "Groups without names will be given an integer name.")
    parser.add_argument('-p', '--filter', type=argparse.FileType('r'),
                                          nargs='?',
                                          default=DEFAULT_EXCLUDED_PATTERNS,
                                          help="A file containing a list of string patterns causing a SMILES to be ignored")
    parser.add_argument('-s', '--sort', action='store_true', 
                                        default=False,
                                        help='Sort protomers for each group by the last extra column')
    params = parser.parse_args(args)
    excluded_patterns = set(l.strip() for l in params.filter)
    separated_protomers = separate_compounds(params.groups, sort=params.sort)
    for grouped in separated_protomers:
        group_protomers(groups=grouped, 
                        group_limits=params.limits,
                        group_names=params.names,
                        filter_patterns=excluded_patterns,
                        output=stdout)
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:], sys.stdout))



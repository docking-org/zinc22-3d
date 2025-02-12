#!/usr/bin/env python

import sys


def main(delimiter, source, prefix, strip_right=True, start=1):
    index = start
    current_handle = None

    try:
        current_name = "{0}{1:d}".format(prefix, index)
        current_handle = open(current_name, 'w')
        for line in source:
            if strip_right:
                comp_line = line.rstrip()
            else:
                comp_line = line

            if comp_line == delimiter:
                current_handle.close()
                index += 1
                current_name = "{0}{1:d}".format(prefix, index)
                current_handle = open(current_name, 'w')

            current_handle.write(line)

        current_handle.close()
        current_handle = None
    except Exception:
        return 1
    else:
        return 0
    finally:
        if current_handle is not None:
            current_handle.close()

if __name__ == '__main__':
    delimiter = sys.argv[1]
    if len(sys.argv) > 2:
        prefix = sys.argv[2]
    else:
        prefix = ''
    sys.exit(main(
        delimiter=delimiter,
        source=sys.stdin,
        prefix=prefix,
        strip_right=True,
        start=1
    ))
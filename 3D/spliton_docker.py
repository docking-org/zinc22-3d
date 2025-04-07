#!/usr/bin/env python

import sys


def main(delimiter, source, prefix, strip_right=True, start=1):
    index = start
    current_handle = None
    buffer = []
    try:
        current_name = "{0}{1:d}".format(prefix, index)
        current_handle = open(current_name, 'w')
        for line in source:
            if strip_right:
                comp_line = line.rstrip()
            else:
                comp_line = line

            if comp_line == delimiter:
                current_handle.write(''.join(buffer))
                buffer = []
                current_handle.close()
                index += 1
                current_name = "{0}{1:d}".format(prefix, index)
                current_handle = open(current_name, 'w')
                continue

            if comp_line.startswith('#') or not comp_line.strip():
                continue
            
            buffer.append(line)

        current_handle.close()
        current_handle = None
    except Exception as e:
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
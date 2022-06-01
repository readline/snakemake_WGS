#!/usr/bin/env python
import os,sys

def main():
    try:
        xcov = sys.argv[1]
    except:
        sys.exit(sys.argv[0]+' [xcov file]')
    data = {}
    with open(xcov) as infile:
        while 1:
            line = infile.readline()
            if not line:
                break
            c = line.strip().split('\t')
            data[c[0]] = c[1]
    if float(data['C.chrY']) / float(data['A2. Target region average depth']) > 0.25:
        return '-y'
    else:
        return ''
if __name__ == '__main__':
    print(main())

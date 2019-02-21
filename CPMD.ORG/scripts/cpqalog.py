#!/usr/bin/env python

"""
Sort cpqa.log file wrt rel err value
"""

__author__ = 'irina'

import sys
import re

REL_ERR_WIDTH = 20

def usage():
    print('USAGE: cpqalog <path to cpqa.log>')

def parse_cmdline(args):
    filename = args[1]
    return filename

def process(filename):
    file = open(filename)

    stat = []
    for line in file:
        if re.search('Problems with', line):
            test_name = line.split()[2]
            stat.append({'name':test_name, 'err': 0})

        if re.search('rel err', line):
            rel_err = line.split()[4]
            stat[len(stat)-1]['err'] = rel_err

    sorted_stat = sorted(stat, key=lambda k: abs(float(k['err'])))
    for x in sorted_stat:
        print(
            '{0:{width}}\t{1}'.
            format(str(x['err']).rjust(REL_ERR_WIDTH), x['name'], width=REL_ERR_WIDTH)
        )

if __name__ == '__main__':
    if len(sys.argv) != 2:
        usage()
        exit()

    params = parse_cmdline(sys.argv)
    process(params)
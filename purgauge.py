"""
"""

import sys

from parse_wflow import convert_block as cb

def inname(B, conf):
    root = "/Users/atlytle/Documents/wflow_out"
    return root + '/bt{0}/fort.{1}'.format(B, conf)

def main(argv):
    print inname(6.17, 1000)

    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
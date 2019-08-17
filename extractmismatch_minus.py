## this script is to find reads whit poly A did not map reference genome but the front map reference genome.
## input filetype is SAM

import os
import gzip, argparse, re
import re

parser = argparse.ArgumentParser(usage="it's usage tip.", description="help info.")
parser.add_argument("-i", "--input", help="the port number.")
parser.add_argument("-o", "--out", help="the file type")


args = parser.parse_args()
i = args.input
o = args.out
sam_i = open(i,"rt")
sam_o = open(o, "w")


for line in sam_i:
        if re.findall("\t\d+S",line):
                if int(re.findall("\t\d+S",line)[0][1:-1]) > 4:
                        sam_o.write(line)
                else:
                        pass

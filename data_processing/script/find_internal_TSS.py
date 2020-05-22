import sys
import Levenshtein



def fastaparser(handle):
    for line in handle:
        if line.startswith(">"):
            name = line.strip().split()[0].strip(">")
            seq = ""
            break
    for line in handle:
        if line.startswith(">"):
            yield name,seq
            name = line.strip().split()[0].strip(">")
            seq = ""
            continue
        seq += line.strip("\n")
    yield name,seq

def bedparser(handle):
    for line in handle:
        list_info = line.strip().split()
        chrs = list_info[0]
        start = list_info[1]
        end= list_info[2]
        strand = list_info[5]
        yield chrs, start, end, strand, line


def reverse_complement(seq):
    ntComplement = {'A':'T', 'C':'G', 'T':'A', 'G':'C', 'N':'N', 'a':'t', 'c':'g', 't':'a', 'g':'c', 'n':'n'}
    revSeqList = list(reversed(seq))
    revComSeqList = [ntComplement[k] for k in revSeqList]
    revComSeq = ''.join(revComSeqList)
    return revComSeq


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("""Usage:
        %s <reference> <bed file>""" % __file__)
        sys.exit()
    fa_file, bed_file = sys.argv[1:3]
    fa_file = open(fa_file)
    bed_file = open(bed_file)

    genome = {}
    for name,seq in fastaparser(fa_file):
        genome[name] = seq
        #print("load genome finish")

    for chrs, start, end, strand, line in bedparser(bed_file):
        if strand == "+":
            strs_3 = genome[chrs][int(start) - 3:int(start)]
            strs_5 = genome[chrs][int(start) - 5:int(start)]
        elif strand == "-":
            strs_3 = genome[chrs][int(end):int(end) + 3]
            strs_5 = genome[chrs][int(end):int(end) + 5]
            strs_3 = reverse_complement(strs_3)
            strs_5 = reverse_complement(strs_5)
        strs_3 = strs_3.upper()
        strs_5 = strs_5.upper()
        counts_3 = strs_3.count("G") + strs_3.count("g")
        counts_5 = strs_5.count("G") + strs_5.count("g")
        counts_5_consecutive_3 = strs_5.count("GGG") + strs_5.count("ggg")
        if ((counts_3)/3 >= 1):
            is_3_3 = 1
        elif ((counts_3)/3 < 1):
            is_3_3 = 0
        if ((counts_3)/3 >= 0.66666):
            is_3_2 = 1
        elif ((counts_3)/3 < 0.66666):
            is_3_2 = 0
        if ((counts_5)/5 >= 0.8):
            is_5_4 = 1
        elif ((counts_5)/5 < 0.8):
            is_5_4 = 0
        if ((counts_5_consecutive_3)/1 >= 1):
            is_5_3 = 1
        elif ((counts_5_consecutive_3)/1 == 0):
            is_5_3 = 0
        sys.stdout.write(line.strip("\n") + "\t" + str(is_3_3)  + "\t" + str(is_3_2) + "\t" + str(is_5_3) + "\t" + str(is_5_4) + "\t" + str((counts_3)/3)  + "\t" + str((counts_5)/5) + "\n")
        #sys.stdout.write(line.strip("\n") + "\t" + strs_3  + "\n")

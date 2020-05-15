import sys

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
            strs_5 = genome[chrs][int(end):int(end) + 5]
            strs_10 = genome[chrs][int(end):int(end) + 10]
            strs_15 = genome[chrs][int(end):int(end) + 15]
            strs_20 = genome[chrs][int(end):int(end) + 20]
            strs_30 = genome[chrs][int(end):int(end) + 30]
            strs_50 = genome[chrs][int(end):int(end) + 50]
        elif strand == "-":
            strs_5 = genome[chrs][int(start) - 5:int(start)]
            strs_10 = genome[chrs][int(start) - 10:int(start)]
            strs_15 = genome[chrs][int(start) - 15:int(start)]
            strs_20 = genome[chrs][int(start) - 20:int(start)]
            strs_30 = genome[chrs][int(start) - 30:int(start)]
            strs_50 = genome[chrs][int(start) - 50:int(start)]
            strs_5 = reverse_complement(strs_5)
            strs_10 = reverse_complement(strs_10)
            strs_15 = reverse_complement(strs_15)
            strs_20 = reverse_complement(strs_20)
            strs_30 = reverse_complement(strs_30)
            strs_50 = reverse_complement(strs_50)
        counts_5 = strs_5.count("A") + strs_5.count("a")
        counts_10 = strs_10.count("A") + strs_10.count("a")
        counts_10_consecutive_6 = strs_10.count("AAAAAA") + strs_10.count("aaaaaa")
        counts_15 = strs_15.count("A") + strs_15.count("a")
        counts_20 = strs_20.count("A") + strs_20.count("a")
        counts_30 = strs_30.count("A") + strs_30.count("a")
        counts_50 = strs_50.count("A") + strs_50.count("a")
        if ((counts_5)/5 >= 1):
            is_5 = 1
        elif ((counts_5)/5 < 1):
            is_5 = 0
        if (counts_10_consecutive_6 >= 1):
            is_10_consecutive = 1
        elif (counts_10_consecutive_6 < 1):
            is_10_consecutive = 0
        if ((counts_10)/10 >= 0.7):
            is_10_7 = 1
        elif ((counts_10)/10 < 0.7):
            is_10_7 = 0
        if ((counts_10)/10 >= 0.8):
            is_10_8 = 1
        elif ((counts_10)/10 < 0.8):
            is_10_8 = 0
        if ((counts_15)/15 >= 0.8):
            is_15 = 1
        elif ((counts_15)/15 < 0.8):
            is_15 = 0
        if ((counts_20)/20 >= 0.75):
            is_20 = 1
        elif ((counts_20)/20 < 0.75):
            is_20 = 0
        if ((strs_30.count("A") + strs_30.count("a") + strs_30.count("T") + strs_30.count("t"))/30 >= 0.9):
            is_30 = 1
        elif ((strs_30.count("A") + strs_30.count("a") + strs_30.count("T") + strs_30.count("t"))/30 < 0.9):
            is_30 = 0
        if ((counts_50)/50 >= 0.65):
            is_50 = 1
        elif ((counts_50)/50 < 0.65):
            is_50 = 0
        sys.stdout.write(line.strip("\n") + "\t" + str(is_5)  + "\t" + str(is_10_consecutive) + "\t" + str(is_10_7) + "\t" + str(is_10_8) + "\t" + str(is_15) + "\t" + str(is_20) + "\t" + str(is_30) + "\t" + str(is_50) +  "\t" + str((counts_5)/5) +  "\t" + str((counts_10)/10) +  "\t" + str((counts_15)/15) +  "\t" + str((counts_20)/20) +  "\t" + str((counts_30)/30) +  "\t" + str((counts_50)/50) + "\n")
        #sys.stdout.write(line.strip("\n") + "\t" + strs_5 + "\t" + strs_6 + "\t" + strs_15 + "\t" + strs_20 + "\t" + strs_50 + "\n")


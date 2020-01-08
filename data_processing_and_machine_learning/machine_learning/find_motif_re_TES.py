import sys
import pandas as pd
import re

def fastaparser(handlei):
    with open(handlei) as handle:
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
    ntComplement = {'A':'T', 'C':'G', 'T':'A', 'G':'C'}
    revSeqList = list(reversed(seq))
    revComSeqList = [ntComplement[k] for k in revSeqList]
    revComSeq = ''.join(revComSeqList)
    return revComSeq

def find_motif(row):
    motif_dict = {}
    if row["strand"] == "+":
        seq = genome[row["chr"]][row["dominant_tes_position"] - 50:row["dominant_tes_position"]].upper()

    if row["strand"] == "-":
        seq = genome[row["chr"]][row["dominant_tes_position"]:row["dominant_tes_position"] + 50].upper()
    #    print(seq)
        seq = reverse_complement(seq)
    #    print(seq)
    for motif in motif_info:
        if re.findall(motif, seq):
            #print( [ i.start() for i in re.finditer(motif, seq)] )
            index = [ i.start() - 50 for i in re.finditer(motif, seq) ]  # motif start position, -50 bp
            #print index
            index = ';'.join(map(str, index))
        else:  # no motif result
            index = 'NA'
        motif_dict[motif] = index
    motif_series = pd.Series(motif_dict)  # motif results, "AATAAA","ATTAAA","AAGAAA","AATACA","AATAGA","AATATA","AATGAA","ACTAAA","AGTAAA","CATAAA","GATAAA","TATAAA","TTTAAA"
    return pd.concat([row, motif_series])

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("""Usage:
        %s <reference> <bed file>""" % __file__)
        sys.exit()
    fa_file, bed_file, out = sys.argv[1:4]
    motif_info = ["AATAAA","ATTAAA","AAGAAA","AATACA","AATAGA","AATATA","AATGAA","ACTAAA","AGTAAA","CATAAA","GATAAA","TATAAA","TTTAAA"]   ##if want more motif ,change here
    df = pd.read_csv(bed_file)


    genome = {}
    print("started")
    for name,seq in fastaparser(fa_file):
        genome[name] = seq
    print("load genome finish")

    result = df.apply(find_motif, axis =  1)
    result.to_csv(out, sep=",", index=False)

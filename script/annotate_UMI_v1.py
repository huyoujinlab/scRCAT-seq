#!/usr/bin/env	python

import os, sys, re
import subprocess
import argparse
#import pysam

"assume the sequence is barcode-UMI-sequence"

def get_opt():
	"""Get options and output help document
	@return args:"""

	parser = argparse.ArgumentParser(prog=sys.argv[0], description='Extract Barcode and UMI information, from 5` start.',
		version='%(prog)s 1.0',
		usage='%(prog)s [options] -F <sequence file> -ID <input file> -O <output file>',
		formatter_class=argparse.RawTextHelpFormatter,
		epilog="Author: YEQI ZHANG")

	require = parser.add_argument_group('******* Option *******')
	require.add_argument('--fastx', '-F', dest='fastx', nargs='?', type=str, required=False, help='input fasta/fastq file')
	require.add_argument('--idxfile', '-ID', dest='idf', nargs='?', type=str, required=False, help='file contained index/name_id info')
	require.add_argument('--idx', '-IDX', dest='idx', nargs='?', type=int, default=3, required=False, help="which column is the 'id name' in --idxfile, 0-based	[default: 3]")
	require.add_argument('--out', '-O', dest='out', nargs='?', type=str, required=False,
							metavar='<out file>', help='output csv file')
	require.add_argument('--type', '-T', dest='type', type=str, choices=['Fastq', 'Fasta'], default="Fastq", help="The type of input file	[default: Fastq]")  # Bam Sam
	require.add_argument('-N', dest='bar', type=int, default=16, help="the bp length of barcode	[default: 16 nt]")
	require.add_argument('-n', dest='umi', type=int, default=10, help="the bp length of UMI	[default: 10 nt]")
	args = parser.parse_args()
	#print(args)
	if not args.fastx or not args.out:
		#parser.print_help()
		parser.print_version()
		parser.print_usage()
		sys.exit(1)
	return args

def build_index(file_type, fastx):
	"""build index file if it is inexistent
	@param file_type:		the type of input file
	@param fastx:			input file"""
	fastx = os.path.abspath(fastx)
	#path = os.path.dirname(fastx)
	if file_type == "Fastq":
		if not os.path.exists("{0}.fai".format(fastx)):
			subprocess.check_call("samtools fqidx {0}".format(fastx), shell=True)
	elif file_type == "Fasta":
		if not os.path.exists("{0}.fai".format(fastx)):
			subprocess.check_call("samtools faidx {0}".format(fastx), shell=True)
	elif file_type == "Bam":
		if not os.path.exists("{0}/{1}.bai".format(path, fastx)):
			subprocess.check_call("samtools index {0}".format(fastx), shell=True)
	elif file_type == "Sam":
		sys.stderr.write("Warn: Sam file cannot build index!")
		pass
	return

def get_sequence(file_type, fastx):
	"""get sequence from fastx
	@param file_type:		the type of input file
	@param fastx:			input file
	@return seq_dict:		{seq_name: sequence}"""
	seq_dict = {}
	with open(fastx, 'r') as db:
		if file_type == 'Fastq':
			for line in db:
				#seq_name = line.strip()[1:]  # skip first '@'
				seq_name = re.search("[\w:]*", line[1:]).group()
				seq_sequence = next(db).strip()
				seq_info = next(db).strip()
				seq_qual = next(db).strip()
				seq_dict.setdefault(seq_name, seq_sequence)
		elif file_type == 'Fasta':
			for line in db:
				seq_name = line.strip()
				seq_sequence = next(db).strip()
				#print seq_name, seq_sequence
				seq_dict.setdefault(seq_name, seq_sequence)

	return seq_dict

def get_bam_sequence():
	return

def extract_info(seq_dict, seq_name, len_b, len_u):
	"""extract Barcode and UMI info
	@param seq_dict:		seq dict
	@param seq_name:		a key within seq dict
	@param len_b:			the length of barcode
	@param len_u:			the length of UMI
	@return:
			############################################
			assume structure:  5`-barcode-UMI-sequence-3`
			############################################
	"""
	barcode = seq_dict[seq_name][0:len_b]
	umi = seq_dict[seq_name][len_b:(len_b + len_u)]
	#print barcode, umi
	return barcode, umi

def main(args):
	fastx = os.path.abspath(args.fastx)
	path = os.path.dirname(fastx)
	if not os.path.exists(fastx):
		sys.stderr.write("Error: File dose not exist!\n")
		sys.exit(1)
	
	#build_index(args.type, args.fastx)
	
	seq_dict = get_sequence(args.type, args.fastx)
	#print seq_dict
	with open(args.idf, 'r') as f, open(args.out, 'w') as o1, open("{0}.barcode.csv".format(args.out[:args.out.rfind('.')]), 'w') as o2:
		count = 0
		o2.write("#ID name\tBarcode\tUMI\n")
		for line in f:
			if line.startswith('#'): o1.write(line)
			line = line.strip().split()
			ID = line[args.idx]
			#print ID
			try:
				barcode, umi = extract_info(seq_dict, ID, args.bar, args.umi)
			except KeyError:
				sys.stderr.write("Sequence ID %s is not in result !\n" % ID)
				continue
			except:
				raise
			o1.write('\t'.join(line + [barcode] + [umi]) + '\n')
			o2.write("{0}\t{1}\t{2}\n".format(ID, barcode, umi))
			
			count +=1
			if count % 10000 == 0: sys.stderr.write("Readlines:	%s\n" % count)
	return

if  __name__ == "__main__":
	args = get_opt()
	main(args)

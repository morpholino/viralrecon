#!/usr/bin/env python
import os
import sys
import re
import argparse
import gzip


def parse_args(args=None):
	Description = 'Find indels positions in bed file'
	Epilog = """Example usage: python parse_mask_bed.py <BED_IN> <BED_OUT>"""
	parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
	parser.add_argument('VCF_IN', help="Input vcf file.")
	parser.add_argument('BED_IN', help="Input bed file.")
	parser.add_argument('BED_OUT', help="Name of the output new bed file.")
	return parser.parse_args(args)

def find_indels_vcf(VcfIn):
	encoding = 'utf-8'
	indels_pos_len={}
	snp = []
	with gzip.open(VcfIn,'r') as f:
		for line in f:
			if '#' not in str(line, encoding):
				line = re.split("\t", str(line, encoding))
				var_pos=line[1]
				ref=line[3]
				alt=line[4]
				if len(alt) > len(ref):
					indels_pos_len[var_pos] = len(alt)
				elif len(ref) > len(alt):
					indels_pos_len[var_pos] = len(ref)
				elif len(ref) == len(alt):
					snp.append(var_pos)
	return(indels_pos_len, snp)

def find_dels_vcf(VcfIn):
	encoding = 'utf-8'
	dels_pos_len={}
	with gzip.open(VcfIn,'r') as f:
		for line in f:
			if '#' not in str(line, encoding):
				line = re.split("\t", str(line, encoding))
				var_pos=line[1]
				ref=line[3]
				alt=line[4]
				if len(ref) > len(alt):
					dels_pos_len[var_pos] = len(ref)
	return(dels_pos_len)

def parse_mask_bed(BedIn,BedOut,indels_pos_len,snp):
	fout = open(BedOut,'w')
	indels_postions=[]
	for pos in indels_pos_len:
		indels_postions.append(pos)
	with open(BedIn) as b:
		snp_analysis = []
		for line in b:
			line = re.split("\t", line)
			ref_genome=line[0]
			init_pos=line[1]
			end_pos=line[2]
			coverage=line[3]
			oline=ref_genome+'\t'+init_pos+'\t'+end_pos+'\t'+coverage
			if len(snp) > 0: 
				output = []
				for position in snp:
					printings = 1
					if int(position) in range(int(init_pos),int(end_pos)+1): ### it is in the range
						if len(output) > 0: ### means it was split at least once before
							### the simplest solution is just to delete the whole mask
							output = []
							printings = 0
						else:
							if int(position) == int(init_pos):
								new_start = int(init_pos) + 1
								oline=ref_genome+'\t'+str(new_start)+'\t'+end_pos+'\t'+coverage
							elif int(position) == int(end_pos):
								new_end = int(end_pos) - 1
								oline=ref_genome+'\t'+init_pos+'\t'+str(new_end)+'\t'+coverage
							else:  ### in the middle somewhere
								new_end = int(position) - 1
								new_start = int(position) + 1
								oline=ref_genome+'\t'+init_pos+'\t'+str(new_end)+'\t'+coverage
								oline2=ref_genome+'\t'+str(new_start)+'\t'+end_pos+'\t'+coverage
								output.append(oline)
								output.append(oline2)
					else:
						pass
				if len(output) > 0:
					for element in output:
						init_poso=element.split('\t')[1]
						end_poso=element.split('\t')[2]
						if int(init_poso) == int(end_poso):
							pass
						else:
							snp_analysis.append(element)
				else:
					if printings == 1:
						init_poso=oline.split('\t')[1]
						end_poso=oline.split('\t')[2]
						if int(init_poso) == int(end_poso):
							pass
						else:
							snp_analysis.append(oline)
			else: ### no snp's need to be checked
				printing = 1
				if len(indels_postions) > 0:
					for position in indels_postions:
						printing = 1
						position_end = int(position)+int(indels_pos_len[position])
						if int(init_pos) < int(position) and int(end_pos) < int(position): ### outside of the indel left side
#							print("condition 1\t" + position + '\t' + str(position_end) + '\t' + oline)
							pass
						elif int(init_pos) < int(position) and int(end_pos) >= int(position) and int(end_pos) <= int(position_end): ### partial covering of the beginning of the gap
#							print("condition 2\t" + position + '\t' + str(position_end) + '\t' + oline)
							new_end = int(position) - 1 
							if new_end == int(init_pos):
								printing = 0
								break
							else:
								oline=ref_genome+'\t'+init_pos+'\t'+str(new_end)+'\t'+coverage
						elif int(init_pos) >= int(position) and int(end_pos) <= int(position_end): ## it is in the gap don't mask anything
#							print("condition 3" + position + '\t' + str(position_end) + '\t' + oline)
							printing = 0
							break
						elif int(init_pos) >= int(position) and int(init_pos) <= int(position_end) and int(end_pos) > int(position_end): ## it covers the gap at the end of it partially
#							print("condition 4\t" + position + '\t' + str(position_end) +'\t' + oline)
							new_start= int(position_end) + 1 
							oline=ref_genome+'\t'+str(new_start)+'\t'+end_pos+'\t'+coverage
						elif int(init_pos) >= int(position_end): ### the low coverage area is after the indel ignore it
#							print("condition 5\t" + position + '\t' + str(position_end) +'\t' + oline)
							pass
						elif int(init_pos) < int(position) and int(end_pos) > int(position_end): ### a super rare situation where the gap is completely covered by a low coverage which extends before and after the gap. I think this would be impossible in real life
							### you need to split the gap
#							print("condition 6\t" + position + '\t' + str(position_end) +'\t' + oline)
							new_start2= int(position_end) + 1
							new_end1 = int(position) - 1 
							oline=ref_genome+'\t'+ init_pos +'\t'+str(new_end1)+'\t'+coverage + ref_genome+'\t'+str(new_start2)+'\t'+end_pos+'\t'+coverage
						else: 
							print('impossible')
					if printing == 1:
						fout.write(oline)
				else:
					oline=ref_genome+'\t'+init_pos+'\t'+end_pos+'\t'+coverage
					fout.write(oline)
		if len(snp_analysis) > 0:
			for line in snp_analysis:
				ref_genome=line.split('\t')[0]
				init_pos=line.split('\t')[1]
				end_pos=line.split('\t')[2]
				coverage=line.split('\t')[3]
				oline=ref_genome+'\t'+init_pos+'\t'+end_pos+'\t'+coverage
				printing = 1
				if len(indels_postions) > 0:
					for position in indels_postions:
						printing = 1
						position_end = int(position)+int(indels_pos_len[position])
						if int(init_pos) < int(position) and int(end_pos) < int(position): ### outside of the indel left side
#							print("condition 1\t" + position + '\t' + str(position_end) + '\t' + oline)
							pass
						elif int(init_pos) < int(position) and int(end_pos) >= int(position) and int(end_pos) <= int(position_end): ### partial covering of the beginning of the gap
#							print("condition 2\t" + position + '\t' + str(position_end) + '\t' + oline)
							new_end = int(position) - 1 
							if new_end == int(init_pos):
								printing = 0
								break
							else:
								oline=ref_genome+'\t'+init_pos+'\t'+str(new_end)+'\t'+coverage
						elif int(init_pos) >= int(position) and int(end_pos) <= int(position_end): ## it is in the gap don't mask anything
#							print("condition 3\t" + position + '\t' + str(position_end) + '\t' + oline)
							printing = 0
							break
						elif int(init_pos) >= int(position) and int(init_pos) <= int(position_end) and int(end_pos) > int(position_end): ## it covers the gap at the end of it partially
#							print("condition 4\t" + position + '\t' + str(position_end) +'\t' + oline)
							new_start= int(position_end) + 1 
							oline=ref_genome+'\t'+str(new_start)+'\t'+end_pos+'\t'+coverage
						elif int(init_pos) >= int(position_end): ### the low coverage area is after the indel ignore it
#							print("condition 5\t" + position + '\t' + str(position_end) +'\t' + oline)
							pass
						elif int(init_pos) < int(position) and int(end_pos) > int(position_end): ### a super rare situation where the gap is completely covered by a low coverage which extends before and after the gap. I think this would be impossible in real life
							### you need to split the gap
#							print("condition 6\t" + position + '\t' + str(position_end) +'\t' + oline)
							new_start2= int(position_end) + 1
							new_end1 = int(position) - 1 
							oline=ref_genome+'\t'+ init_pos +'\t'+str(new_end1)+'\t'+coverage + ref_genome+'\t'+str(new_start2)+'\t'+end_pos+'\t'+coverage
						else: 
							print('impossible')
					if printing == 1:
						fout.write(oline)
				else:
					oline=ref_genome+'\t'+init_pos+'\t'+end_pos+'\t'+coverage
					fout.write(oline)

########More def functions
def main(args=None):
	args = parse_args(args)
	dels_pos_len,snp = find_indels_vcf(args.VCF_IN)
	parse_mask_bed(args.BED_IN, args.BED_OUT,dels_pos_len,snp)

if __name__ == '__main__':
	sys.exit(main())

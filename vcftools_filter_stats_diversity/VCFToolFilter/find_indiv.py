import sys
import os
import re

def get_field_samples_options(dataset):
	options = []
	line=os.popen("grep '#CHROM' %s"%dataset.file_name).read()[:-1].split('\t')
	index=line.index('FORMAT')
	for opt in line[index+1:] :
		options.append((opt,opt, True))
	return options

def get_field_chrs_options(dataset):
        options = []
	chrs=os.popen("grep '##contig' %s"%dataset.file_name).read()[:-1].split('\n')
	if len(chrs)>1:
		for line in chrs:
			opt=re.search('^##contig=<ID=(\w+),length=',line).group(1)
			options.append((opt,opt, True))
	else :
		ok=0
		chrs=[]
		file_in=open(dataset.file_name,'r')
		line=file_in.readline()
		while line : 
			if ok : 
				chrom=line.split('\t')[0]
				if chrom not in chrs :
					chrs.append(chrom)
			if re.search('^#CHROM',line):
				ok=1
			line=file_in.readline()
		print chrs
        	for opt in chrs:
                	options.append((opt,opt, True))
		file_in.close()
        return options


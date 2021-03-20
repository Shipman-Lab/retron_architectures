"""Import Modules"""
from Bio import SeqIO
import os, sys
import xlsxwriter
import fuzzysearch
import re
from collections import Counter

"""Line Arguments"""
condition = sys.argv[1] #file name (minus .fastq)
library = sys.argv[2] #stem or a1a2

"""Globals"""
user_profile = os.environ ['USERPROFILE']
wd = os.getcwd()
fastq_reads = '%s/%s.fastq' % (wd,condition)

ec86_1L = 'GCATTGAA'
ec86_1R = 'GTAAGGGT'
ec86_2R = 'ACTTTCAT'
BsaI_L = 'GGTCTCA'
BsaI_R = 'AGAGACC'

count = 0
count_dict = {}

if library == 'stem': 
    lig_sites = [ec86_1L, ec86_1R]
    synth_numbers = [4235,4267] 
if library == 'a1a2': 
    lig_sites = [ec86_1L, ec86_2R]
    synth_numbers = [5267,5297]

match_count = 0
no_lig_site_seqs = []
no_synth_match_seqs = []

"""Run"""
#Create Results folder
newpath = ((r'%s/%s_Results') % (wd,condition))
if not os.path.exists(newpath): os.makedirs(newpath)
#pull in synthesis_parts
part_counts_trimmed = {} #trimmed seq is key: fasta number, UMIs (NOTE: not true UMIs)
parts_trimmed_list = []
synth_dict = SeqIO.index("%s/Library_Parts.fasta" % wd, "fasta")
for part in range(synth_numbers[0],synth_numbers[1]):
    temp_part = synth_dict[str(part)]
    bsa_trimmed = re.search('%s(.*)%s' % (BsaI_L,BsaI_R), str(temp_part.seq))
    part_counts_trimmed[bsa_trimmed.group(1)[4:-4]] = (str(part), [])
    parts_trimmed_list.append(bsa_trimmed.group(1)[4:-4])
for seq_record in SeqIO.parse(fastq_reads, "fastq"):
    temp_record = seq_record.reverse_complement()
    lig_trimmed = re.search('%s(.*)%s' % (lig_sites[0],lig_sites[1]), str(temp_record.seq))
    if lig_trimmed:
        try:
            part_counts_trimmed[lig_trimmed.group(1)][1].append(str(seq_record.seq)[:5])
            match_count += 1
        except KeyError:
            no_synth_match_seqs.append(lig_trimmed.group(1))
    else: no_lig_site_seqs.append(temp_record)
c = Counter(no_synth_match_seqs)

"""Output"""
#write excel
workbook = xlsxwriter.Workbook('%s/%s_Results/Plasmid_Library.xlsx' % (wd,condition))
worksheet = workbook.add_worksheet()
bold = workbook.add_format({'bold': True})
#Add titles (row, col: zero referenced)
worksheet.write(0,0,'fasta #')
worksheet.write(0,1,'variable sequence')
worksheet.write(0,2,'total counts')
# worksheet.write(0,3,'UMIs')
#add data
row = 1
col = 0
for part in parts_trimmed_list:
    worksheet.write(row,col,part_counts_trimmed[part][0], bold)
    worksheet.write(row,col+1,part)
    worksheet.write(row,col+2,len(part_counts_trimmed[part][1]))
    row += 1
for seq in c.most_common():
    if seq[1] >= 10:
        worksheet.write(row,col+1,seq[0])
        worksheet.write(row,col+2,seq[1])
        row += 1
workbook.close()












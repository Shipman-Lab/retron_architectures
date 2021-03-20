"""Import Modules"""
from Bio import SeqIO
from Bio.Seq import Seq
import os, sys
import xlsxwriter
import re
from collections import Counter

"""Line Arguments"""
condition = sys.argv[1] #file name (minus .fastq)
library = sys.argv[2] #stem or a1a2

"""Globals"""
user_profile = os.environ ['USERPROFILE']
wd = os.getcwd()
fastq_reads = '%s/%s.fastq' % (wd,condition)

trim_parts = ('GGTCTCATGAA','GTAAAGAGACC')

count = 0
count_dict = {}

if library == 'stem': 
    synth_numbers = [4235,4267] 
if library == 'a1a2': 
    synth_numbers = [5267,5297]

total_seqs = 0
match_count = 0
no_first_T = []
no_post_flank = []
no_lib_match = []

"""Run"""
#Create Results folder
newpath = ((r'%s/%s_Results') % (wd,condition))
if not os.path.exists(newpath): os.makedirs(newpath)
#pull in synthesis_parts
part_counts_trimmed = {} #trimmed seq is key: fasta number, sequences
parts_trimmed_list = []
synth_dict = SeqIO.index("%s/Library_Parts.fasta" % wd, "fasta")
for part in range(synth_numbers[0],synth_numbers[1]):
    temp_part = synth_dict[str(part)]
    part_trimmed = re.search('%s(.*)%s' % (trim_parts[0],trim_parts[1]), str(temp_part.seq))
    try:   #if the trimmed part exists multiple times, change fasta to ambiguous
       temp_fasta = part_counts_trimmed[part_trimmed.group(1)[1:]][0]
       part_counts_trimmed[part_trimmed.group(1)[1:]] = ('%s_%s' % (str(temp_fasta),str(part)), [])
    except KeyError:
        part_counts_trimmed[part_trimmed.group(1)[1:]] = (str(part), []) #NOTE: assuming msd ends one base earlier than lit
    parts_trimmed_list.append(part_trimmed.group(1)[1:])

for seq_record in SeqIO.parse(fastq_reads, "fastq"):
    total_seqs += 1
    if str(seq_record.seq)[0] == 'T':      
        seq_trimmed = re.split('CCCCC', str(seq_record.seq)[1:])
        if seq_trimmed: 
            rev_seq = Seq(seq_trimmed[0]).reverse_complement()
            try:
                part_counts_trimmed[rev_seq][1].append(seq_record)
                match_count += 1
            except KeyError:
                no_lib_match.append(str(rev_seq)) 
        else: no_post_flank.append(seq_record)    
    else: no_first_T.append(seq_record)
c = Counter(no_lib_match)

"""Output"""
#write excel
workbook = xlsxwriter.Workbook('%s/%s_Results/RTDNA_Library.xlsx' % (wd,condition))
worksheet = workbook.add_worksheet()
bold = workbook.add_format({'bold': True})
#Add titles (row, col: zero referenced)
worksheet.write(0,0,'fasta #')
worksheet.write(0,1,'variable sequence')
worksheet.write(0,2,'total counts')
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

print '%s total sequences' % total_seqs
print '%s matches found' % match_count
print '%s do not contain an initial T' % len(no_first_T)
print '%s do not have an appropriate nucleotide extension' % len(no_post_flank)
print '%s were trimmed, but have no match in the target set' % len(no_lib_match)











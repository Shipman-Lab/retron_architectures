"""Import Modules"""
from Bio import SeqIO
import os, sys
import xlsxwriter
import fuzzysearch
import re
from collections import Counter

"""Line Arguments"""
condition = sys.argv[1] #file name (minus .fastq)

"""Globals"""
user_profile = os.environ ['USERPROFILE']
wd = os.getcwd()
fastq_reads = '%s/%s.fastq' % (wd,condition)

lig_sites = ['GCATTGAA','ACTTTCAT']
synth_numbers = [5267,5297] 
PreBC = 'TGTTGGAA'
PostBC = 'AGCCAACC'
PreBC_part = 'TGTTGGAA'
PostBC_part = 'AGCCAACC'

BC_dict = {} #loop BC is key: [0]fasta number, [1]total matches (seqs)
BC_list = []

total_seqs = 0
no_flanking_seqs = []
no_BC_match = []

"""Run"""
#Create Results folder
newpath = ((r'%s/%s_Results') % (wd,condition))
if not os.path.exists(newpath): os.makedirs(newpath)
#pull in synthesis_parts
synth_dict = SeqIO.index("%s/Library_Parts.fasta" % wd, "fasta")
for part in range(synth_numbers[0],synth_numbers[1]):
    temp_part = synth_dict[str(part)]
    loop_trimmed = re.search('%s(.*)%s' % (PreBC_part,PostBC_part), str(temp_part.seq))
    BC_dict[loop_trimmed.group(1)] = (str(part), []) #adds fasta
    BC_list.append(loop_trimmed.group(1))
#pull in sequences
for seq_record in SeqIO.parse(fastq_reads, "fastq"):
    total_seqs += 1
    temp_record = seq_record.reverse_complement()
    loop_trimmed = re.search('%s(.*)%s' % (PreBC,PostBC), str(temp_record.seq))
    if loop_trimmed:
        try:
            BC_dict[loop_trimmed.group(1)][1].append(str(temp_record))
        except KeyError:
            no_BC_match.append(loop_trimmed.group(1))
    else: no_flanking_seqs.append(temp_record)

print '%s total sequences' % total_seqs
print '%s do not contain matching loop flanking sequences' % len(no_flanking_seqs)
print '%s do not contain a matching barcode' % len(no_BC_match)

unmatched_BC_counter = Counter(no_BC_match)

"""Output"""
#write excel
workbook = xlsxwriter.Workbook('%s/%s_Results/Barcoded_Library.xlsx' % (wd,condition))
worksheet = workbook.add_worksheet()
bold = workbook.add_format({'bold': True})
#Add titles (row, col: zero referenced)
worksheet.write(0,0,'fasta #')
worksheet.write(0,1,'BC')
worksheet.write(0,2,'total counts')
#add data
row = 1
col = 0
for BC in BC_list:
    worksheet.write(row,col,BC_dict[BC][0], bold)
    worksheet.write(row,col+1,BC)
    worksheet.write(row,col+2,len(BC_dict[BC][1]))
    row += 1

for seq in unmatched_BC_counter.most_common():
    if seq[1] >= 10:
        worksheet.write(row,col+1,seq[0])
        worksheet.write(row,col+2,seq[1])
        row += 1
workbook.close()












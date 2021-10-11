import sys,os
from Bio import SeqIO
import xlsxwriter
import fuzzysearch
import re
from collections import Counter
from Bio.Seq import Seq

"""Line Arguments"""
condition = sys.argv[1] #file name (minus .fastq)
edit_site = sys.argv[2] #e.g. rpoB_A1547G

"""Globals"""
user_profile = os.environ ['USERPROFILE']
wd = os.getcwd()
fastq_reads = '%s/%s.fastq' % (wd,condition)

Region_dict = {'rpoB_A1547G': {'flanking': ['GAGTTCTTCGGTTCCAGCCA','GTTTATGGACCAGAACAACC'], 'wt_or_edited': ['GCTGTCTCA','GCTGCCTCA']},
			   'fabH_G954A': {'flanking': ['TTGCGTCATGTTTTAATCCT','AACGAACCAGCGCGGAGCCC'], 'wt_or_edited': ['TATCCTAGA','TATCTTAGA']},
			   'fliN_G414A': {'flanking': ['GCACAGTAGCGTGGTTATTC','GGCTCAGGCGGCGCATTCGC'], 'wt_or_edited': ['ATCACTAAC','ATCATTAAC']},
			   'murF_G1359A': {'flanking': ['TGACCAAATGTTCGGCCAGC','ATGTCCCATTCTCCTGTAAA'], 'wt_or_edited': ['CAAACTAAC','CAAATTAAC']},
			   'priB_G1069A': {'flanking': ['CGACGGAAATAACGTGCCAT','CTCCAGAATCTATCAATTCA'], 'wt_or_edited': ['ATGGCTAGT','ATGGTTAGT']},
			   'othersite': {'flanking': ['TK','TK'], 'wt_or_edited': ['wt','edited']}}

all_reads_str = []
outcomes_dict = {'wt':0,
				 'edited':0,
				 'undetermined_no_flanking_match':0,
				 'undetermined_no_site_match':0}

"""Defs"""
def extract_and_match(sequence):
	left_flank = fuzzysearch.find_near_matches(Region_dict[edit_site]['flanking'][0],sequence,max_l_dist=4)
	right_flank = fuzzysearch.find_near_matches(Region_dict[edit_site]['flanking'][1],sequence,max_l_dist=4)
	if len(left_flank) == 1 and len(right_flank) == 1:
		region = sequence[left_flank[0].end:right_flank[0].start]
		if region == Region_dict[edit_site]['wt_or_edited'][0]:
			return 'wt'
		elif region == Region_dict[edit_site]['wt_or_edited'][1]:
			return 'edited'
		else: return 'undetermined_no_site_match'
	else: return 'undetermined_no_flanking_match'

"""Run"""
#Create Results folder
newpath = ((r'%s/%s_Results') % (wd,condition))
if not os.path.exists(newpath): os.makedirs(newpath)
#pull in sequences and process into counts of unique reads
for seq_record in SeqIO.parse(fastq_reads, "fastq"):
	all_reads_str.append(str(seq_record.seq))
read_counter = Counter(all_reads_str)
for read in read_counter:
	outcomes_dict[extract_and_match(read)] += read_counter[read]
edited_percent = (float(outcomes_dict['edited'])/float(outcomes_dict['edited']+outcomes_dict['wt']))*100

"""Output"""
print '%s percent edited (%s edited/%s wt)' % (edited_percent,outcomes_dict['edited'],outcomes_dict['wt'])
print '%s do not contain clean flanking regions' % outcomes_dict['undetermined_no_flanking_match']
print '%s contain flanking regions, but not a clean target site match to edited or wt' % outcomes_dict['undetermined_no_site_match']

#write excel
workbook = xlsxwriter.Workbook('%s/%s_Results/site_analysis.xlsx' % (wd,condition))
worksheet = workbook.add_worksheet()
bold = workbook.add_format({'bold': True})
#Add titles (row, col: zero referenced)
worksheet.write(0,1,condition,bold)
worksheet.write(1,0,'Percent Edited',bold)
worksheet.write(2,0,'Edited Counts',bold)
worksheet.write(3,0,'Wild-Type Counts',bold)
worksheet.write(4,0,'undetermined_no_flanking_match',bold)
worksheet.write(5,0,'undetermined_no_site_match',bold)
#add data
worksheet.write(1,1,edited_percent)
worksheet.write(2,1,outcomes_dict['edited'])
worksheet.write(3,1,outcomes_dict['wt'])
worksheet.write(4,1,outcomes_dict['undetermined_no_flanking_match'])
worksheet.write(5,1,outcomes_dict['undetermined_no_site_match'])
workbook.close()




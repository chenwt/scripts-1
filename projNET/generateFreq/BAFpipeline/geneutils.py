#! /usr/bin/env python

# utility functions for parsing .bam, .pileup, .vcf files

import re
import pickle
from operator import itemgetter
import sys
import subprocess
import os
import string
import time
import shutil
import datetime
import collections
import csv
import unittest

REF_PATH = '/ifs/data/c2b2/ngs_lab/ngs/resources/bwa_samtools_gatk_DB/human_g1k_v37.fasta'
GC_WINDOW_SIZE = 213


def getMeanAndVarianceOfD(counts_d):
    num_observations = sum(counts_d.values())
    num_counts = 0
    for count in counts_d:
        num_counts += count * counts_d[count]
    mean = num_counts / float(num_observations)
    sum_squared_deviations = 0
    for count in counts_d:
        sum_squared_deviations += (count - mean) * (count - mean) *\
             counts_d[count]
    variance = sum_squared_deviations / float(num_observations)
    return mean, variance

# given a list of exons in the format
# [['22', 16084594, 16084833], ['22', 16100468, 16100587],...] and a
# window size, returns a list of roughly equally-sized windows in the format
# [['22', 16084594, 16084713], ['22', 16084714, 16084833],...]
def getWindowL(contig_exon_l, window_size):
    exon_length_l = [e[2] - e[1] + 1 for e in contig_exon_l]
    windows_per_exon_l = [max(1, l / window_size) for l in exon_length_l]
    window_l = []
    for i in range(len(contig_exon_l)):
        window_width = exon_length_l[i] / windows_per_exon_l[i]
        for j in range(windows_per_exon_l[i]):
            window_l.append([contig_exon_l[i][0],\
                            contig_exon_l[i][1] + window_width * j,\
                            contig_exon_l[i][1] + window_width * (j + 1) - 1])
        window_l[len(window_l) - 1][2] = contig_exon_l[i][2]
    return window_l

# given a list of exons in the format
# [['1', 14642, 14881], ['1', 14943, 15062],...], and a
# window size, returns a list of equally-sized windows in the format
# [['1', 14662, 14761], ['1', 14762, 14861], ['1', 14953, 15052],...]
# discards exons smaller than window size, and places the windows in the
# center of the exon
def getEquallySizedWindowsL(contig_exons_l, window_size):
    window_l = []
    for exon in contig_exons_l:
        exon_length = exon[2] - exon[1] + 1
        windows_per_exon = exon_length / window_size
        if windows_per_exon > 0:
            extra_bases = exon_length - windows_per_exon * window_size
            extra_bases_left = extra_bases / 2
            for j in range(windows_per_exon):
                left_boundary = exon[1] + extra_bases_left + j * window_size
                right_boundary = left_boundary + window_size - 1
                window_l.append([exon[0], left_boundary, right_boundary])
    return window_l

def getGCPercentL(window_l):
    GCregex = re.compile(r"G|C")
    GCPercent_l = [] 
    contig_offset_d = getContigOffsetD()    
    for window in window_l:
        contig = window[0]
        startbase = window[1] - GC_WINDOW_SIZE / 2
        endbase = window[2] + GC_WINDOW_SIZE / 2
        number_bases = endbase - startbase + 1
        sequence = getSequence(contig, startbase, endbase, contig_offset_d)
        GCPercent_l.append(window +\
            [int ( round (100 * len (GCregex.findall(sequence)) /
                    float (number_bases) ) )])
    return GCPercent_l

def getContigExonsL(exons_l, contig):
    return filter(lambda s:s[0] == contig, exons_l)

def WindowReadCountdGenerator(windows_l, bam_fn, with_chr_flag):
    """Generator that returns a dictionary for each window in windows_l
    
    with the number of reads with a midpoint
    over each base in each window in window_l
    """
    for window in windows_l:
        count_d = {}
        startbase = window[1]
        endbase = window[2]
        midpoint_d = getMidpointD(window, bam_fn, with_chr_flag)
        for base in range(startbase, endbase + 1):
            if base in midpoint_d.keys():
                count_d[base] = len(midpoint_d[base])
            else:
                count_d[base] = 0
        yield count_d

def getMidpointD(window, bam_fn, with_chr_flag):
    contig = window[0]
    if with_chr_flag:
       contig = "chr" + contig
    region = contig + ':' + str(window[1]) + "-" + str(window[2])
    samtoolsin = subprocess.Popen(["samtools", "view", bam_fn, region],\
                                      stdout = subprocess.PIPE, bufsize = 1)
    midpoint_d = {}
    for samline in samtoolsin.stdout:
        samline = SAMLineOb(samline)
        if not samline.isUnmapped():
            midpoint = int(samline.getMidpoint()[1])
            if midpoint >= window[1] and midpoint <= window[2]:
                templateName = samline.getTemplateName()
                if midpoint not in midpoint_d.keys():
                    midpoint_d[midpoint] = set()
                midpoint_d[midpoint].add(templateName)
    return midpoint_d

def getMeanAndVariance(counts_l):
    numObservations = sum(counts_l)
    if numObservations == 0:
        return (None, None)
    numCounts = 0
    for indx in range(len(counts_l)):
        numCounts += indx * counts_l[indx]
    mean = numCounts / float(numObservations)
    sumSquaredDeviations = 0
    for indx in range(len(counts_l)):
        sumSquaredDeviations += (indx - mean) * (indx - mean) *\
                                counts_l[indx]
    variance = sumSquaredDeviations / float(numObservations)
    return((mean, variance))

def getWindowLForExon(contig, startbase, endbase, window_size):
    exon_length = endbase - startbase + 1
    num_windows = max(1, exon_length / window_size)
    window_l = []
    window_width = exon_length / num_windows
    for j in range(num_windows):
        window_l.append([contig, startbase + window_width * j, startbase + window_width * (j + 1) - 1])
    window_l[len(window_l) - 1][2] = endbase
    return window_l

# given a bedfile with coordinates for exonic regions, returns a list
# of those regions.  because the chromEnd of a feature in a bed file list
# is not included in the exonic region, the chromEnd in the bedfile
# is reduced by 1.
# the list returned is in the format [['1', 14642, 14881], ['1', 14943, 15062],...]
def getExonL(bedfile):
    exons_BED_ob = BEDFileOb(fn = bedfile)
    exon_l = []
    for r in exons_BED_ob:
        exon_l.append([r.getStartPos_t()[0], r.getStartPos_t()[1],\
            r.getEndPos_t()[1] - 1])
    return exon_l

def getAverageTemplateLength():
    sam=geneutils.SAMFileOb(fn = "/ifs/scratch/c2b2/ngs_lab/db2175/July/July23CNV/temp/s20.sam", memory_flag = False)
    total=0
    count=0
    import math
    for samline in sam:
        tlen = math.fabs(samline.getTemplateLength())
        #ignore reads whose two ends map to far away parts of the chromosome
        if tlen < 1000:
            total = total + tlen
            count = count + 1
    print total/float(count)

def getSequence(contig, startbase, endbase, contig_offset_d):
    with open(REF_PATH, 'rb') as fh:       
        start = contig_offset_d[contig] + startbase + startbase / 60
        end = contig_offset_d[contig] + endbase + endbase / 60
        fh.seek(start)
        sequence = fh.read(end - start + 1)           
    return string.replace(sequence, '\n', '')

def getContigLength_d():
    FAI = REF_PATH + ".fai"
    FAI_ob = FAIFileOb(fn = FAI)
    contig_length_d = {}
    for FAIline in FAI_ob:
        contig_length_d[FAIline.getContig()] = FAIline.getLength()
    return contig_length_d

def getContigOffsetD():
    FAI = REF_PATH + ".fai"
    FAI_ob = FAIFileOb(fn = FAI)
    contig_offset_d = {}
    for FAIline in FAI_ob:
        contig_offset_d[FAIline.getContig()] = FAIline.getOffset()
    return contig_offset_d

def ContigWindowGenerator(contig, window_size, margin, contig_length_d, start):
    contig_length = contig_length_d[contig]
    windowstartpos = start + margin + 1
    while windowstartpos < contig_length - margin:
        yield [contig, windowstartpos, windowstartpos + window_size - 1]
        windowstartpos = windowstartpos + window_size

#window_size should be odd
def WindowGCcontentdGenerator(window_l, gc_window_size):
    """Generator that returns a dictionary summarizing the GC content 
       of a genomic window.
       window_l is a list of lists, each entry in the format 
       ['22', 16200000, 16200003]
       gc_window_size should be an odd number, specifying the size 
       of the window centered around each base
    >>> window_gc_d_g = WindowGCcontentdGenerator([['22', 16200000, 16200003]], 5)
    >>> window_gc_d_g.next()
    {16200000: 3, 16200001: 3, 16200002: 3, 16200003: 3}
    >>> window_gc_d_g = WindowGCcontentdGenerator([['22', 16199800, 16199993]], 5)
    >>> window_gc_d = window_gc_d_g.next()
    >>> window_gc_d[16199993]
    0
    """
    GCregex = re.compile(r"G|C")    
    FAI = REF_PATH + ".fai"
    FAI_ob = FAIFileOb(fn = FAI)
    contig_offset_d = {}
    for c in FAI_ob:
        contig_offset_d[c.getContig()] = c.getOffset()
    for window in window_l:
        contig = window[0]
        startbase = window[1]
        endbase = window[2]
        sequence = getSequence(contig, startbase - gc_window_size / 2, \
                       endbase + gc_window_size / 2, contig_offset_d)
        GC_d = {}
        GC_d[startbase] = len(GCregex.findall(sequence[0:gc_window_size]))
        for i in range(startbase + 1, endbase + 1):
            GC_d[i] = GC_d[i-1] + \
                 (sequence[i - startbase - 1 + gc_window_size] in 'GC') -\
                 (sequence[i - 1 - startbase] in 'GC')
        yield GC_d

class GCContentOb(object):
    def __init__(self, gc_window_size):
        self.GCregex = re.compile(r"G|C")    
        FAI = REF_PATH + ".fai"
        FAI_ob = FAIFileOb(fn = FAI)
        self.gc_window_size = gc_window_size
        self.contig_offset_d = {}
        for c in FAI_ob:
            self.contig_offset_d[c.getContig()] = c.getOffset()
        
    
    def getGCContentd(self, contig, startbase, endbase):
        GC_d = {}
        sequence = getSequence(contig, startbase - self.gc_window_size / 2, \
                       endbase + self.gc_window_size / 2, self.contig_offset_d)
        GC_d = {}
        GC_d[startbase] = len(self.GCregex.findall(sequence[0:self.gc_window_size]))
        for i in range(startbase + 1, endbase + 1):
            GC_d[i] = GC_d[i-1] + \
                 (sequence[i - startbase - 1 + self.gc_window_size] in 'GC') -\
                 (sequence[i - 1 - startbase] in 'GC')


    
def MidpointdGenerator(number_lines_to_get, samorbam_fn, verbose = False, stream = False):
    """Generator that returns a dictionary summarizing number_lines_to_get lines
       sam_fn at a time.  Each key is the position of a midpoint of a read.
       Each value is a set containing the template names of the templates with
       midpoints at that position.
       
       >>> midpoint_d_g = MidpointdGenerator("/ifs/scratch/c2b2/ngs_lab/db2175/July/July23CNV/temp/s22.sam", 5)
       >>> midpoint_d = midpoint_d_g.next()
       >>> len(midpoint_d)
       5
       >>> midpoint_d[16084452]
       set(['D8GSQ5P1:2:1103:14120:54510#0'])
       >>> midpoint_d = midpoint_d_g.next()
       >>> midpoint_d[16084648]
       set(['D8GSQ5P1:2:1307:14299:173976#0'])
       >>> len(midpoint_d)
       5
    """
    if stream:
        samtoolsin = subprocess.Popen(["samtools", "view", samorbam_fn],\
                                      stdout = subprocess.PIPE, bufsize = 1)
        sam = samtoolsin.stdout
    else:
        sam = SAMFileOb(samorbam_fn, memory_flag = False)
    count = 0
    midpoint_d = {}
    #import pdb; pdb.set_trace()
    for samline in sam:
        #import pdb; pdb.set_trace()
        if stream:
            samline = SAMLineOb(samline)
        if not samline.isUnmapped():
            count += 1
            midpoint = int(samline.getMidpoint()[1])
            templateName = samline.getTemplateName()
            if samline.getTemplateLength() < 1000:
                if midpoint not in midpoint_d.keys():
                    midpoint_d[midpoint] = set()
                midpoint_d[midpoint].add(templateName)
        if count == number_lines_to_get:
            yield midpoint_d
            midpoint_d = {}
            count = 0
    yield midpoint_d
    


def NaiveWindowReadCount(sam_fn, window):
    """For testing purposes, counts reads that have a midpoint 
       in the specified window.
       Window should be specified as a list in the 
       format ['22', 16200000, 16200003]
       >>> NaiveWindowReadCount("/ifs/scratch/c2b2/ngs_lab/db2175/July/July23CNV/temp/s22.sam", ['22', 16255773, 16256012])
       86
       >>> NaiveWindowReadCount("/ifs/scratch/c2b2/ngs_lab/db2175/July/July23CNV/temp/s22.sam", ['22', 50706255, 50706436])
       9
    """
    sam = SAMFileOb(sam_fn, memory_flag = False)
    startbase = window[1]
    endbase = window[2]
    count_d = {}
    count_set = set()
    for samline in sam:
        if not samline.isUnmapped():
            if int(samline.getMidpoint()[1]) >= startbase and\
               int(samline.getMidpoint()[1]) <= endbase:
                count_set.add(samline.getTemplateName())
            if samline.getMidpoint()[1] >  endbase + 4000 and\
               samline.getTemplateLength() < 1000:
                return len(count_set)

def countReadsInWindows(window_l, sam_path):
    window_counts_l = [c + [0] for c in window_l]
    curr_left_win = 0
    curr_right_win = 0
    num_open_windows = 5
    from collections import deque
    templates_in_open_windows = deque()
    templates_in_open_windows.append(set())
    sam = SAMFileOb(fn = sam_path, memory_flag = False)
    duplicates = 0
    nonduplicates = 0
    for samline in sam:
        #import pdb; pdb.set_trace()
        if not samline.isUnmapped():
            midpoint = samline.getMidpoint()
            if midpoint[1] > window_counts_l[curr_right_win][2]:
                if curr_right_win < len(window_counts_l) - 1:
                    curr_right_win += 1
                templates_in_open_windows.append(set())
                if curr_right_win > curr_left_win + num_open_windows - 1:
                    curr_left_win += 1
                    templates_in_open_windows.popleft()
            for i in range(curr_left_win, curr_right_win + 1):
                if midpoint[1] >= window_counts_l[i][1] and\
                   midpoint[1] <= window_counts_l[i][2]:
                    if samline.getTemplateName() not in\
                        templates_in_open_windows[i - curr_left_win]:
                       window_counts_l[i][3] += 1
                       templates_in_open_windows[i - curr_left_win].add(\
                                 samline.getTemplateName())
                       nonduplicates += 1
                    elif samline.getTemplateName() in\
                        templates_in_open_windows[i - curr_left_win]:
                       duplicates += 1 
                    break
    print "duplicates", duplicates
    print "nonduplicates", nonduplicates
    return window_counts_l

def generateSamFiles(bam_path, probe_bed_path, target_dir_path):
    exon_l = getExonL(probe_bed_path)
    contig_s = get_contig_s(exon_l)
    bed_dir_path = target_dir_path
    sam_dir_path = target_dir_path
    for contig in sorted(list(contig_s)):
        contig_exon_l = filter(lambda s:s[0] == contig, exon_l)
        bed_path = bed_dir_path + '/' + 'b' + contig + '.bed'
        writeBedFileFromRegionL(contig_exon_l, bed_path)
        sam_path = sam_dir_path + '/' + 's' + contig + '.sam'
        triggerSamfile(bam_path, bed_path, sam_path, contig)


# given a list of tuples representing genomic positions, returns a set of the contigs represented by the set
def get_contig_s(pos_list):
    return set(map(itemgetter(0), pos_list))

# merges all pileup files in a directory into a single file. the destination file should not have a .pileup extension (use .PILEUP)
def merge_pileups(pileup_file_dir_path, merged_pileup_path):
    with open(merged_pileup_path, 'wb') as destination:
        for filename in os.listdir(pileup_file_dir_path):
            if filename.endswith('.pileup'):
                with open(pileup_file_dir_path + '/' + filename, 'rb')\
                  as in_:
                   shutil.copyfileobj(in_, destination)

def make_dirs(path_l):
    for path in path_l:
        if not os.path.exists(path):
            os.makedirs(path)

def getTimeStamp():
    dt_ob = datetime.datetime.now()
    return "y" + str(dt_ob.year) + "m" + str(dt_ob.month) +\
         "d" + str(dt_ob.day) + "h" + str(dt_ob.hour) +\
         "m" + str(dt_ob.minute) + "s" + str(dt_ob.second) + \
         "mu" + str(dt_ob.microsecond)

def appendListToFile(l, path):
    new_list = []
    for line in l:
        new_list.append("\t".join([str(f) for f in line]))
    with open(path, 'a') as fh:
        fh.write("\n".join([line for line in new_list]))

# this function overwrites any file it writes to
def writeListToFile(l, path):
    with open(path, 'wb') as fh:
        fh.write("\n".join([line for line in l]))

def convertToStringAndWriteListToFile(l, path):
    new_list = []
    for line in l:
        new_list.append("\t".join([str(f) for f in line]))
    writeListToFile(new_list, path)

# improve this function so it can accommodate files with headers
def readListFromFile(path, field_l):
    output_l = []
    with open(path, 'rb') as fh:
        for line in csv.reader(fh, delimiter = '\t'):
            #import pdb; pdb.set_trace()
            output_l.append([line[i] for i in field_l])
    return output_l

def writeBedFileFromRegionL(exon_l, bed_path):
    bed_l = []
    for line in exon_l:
        line[2] += 1
        bed_l.append(line)
    convertToStringAndWriteListToFile(bed_l, bed_path)

def trigger_samfile_step(bed_job_info_ob, contig, sample_name, bam_file_path):
    samfile_job_info_ob = bed_job_info_ob.getBabyJobInfoOb(\
                            task_l = ["getsam"],\
                            sample_l = [sample_name],\
                            extension_s = "sam",\
                            temp_flag = True)
    script_str = '/ifs/scratch/c2b2/ngs_lab/db2175/BashScripts/getreadscontig.sh'
    command = '   '.join([script_str,
                              bed_job_info_ob.output_file_path,\
                              samfile_job_info_ob.output_file_path,\
                              bam_file_path])
    submit_qsub_job(command, samfile_job_info_ob.job_name,\
                    samfile_job_info_ob.log_error_file_path,\
                    samfile_job_info_ob.log_stdout_file_path)
    return samfile_job_info_ob

# submits a qsub job.  the s.wait() command prevents too many processes from being active at the same time
def submit_qsub_job(command, job_n, qsub_error_path,\
                    qsub_output_path, mem='2G', time=15):
    qsub_script_str = "".join(["qsub -l mem=", mem, ",time=",\
                              str(time), ":: "])
    qsub_command = '   '.join([qsub_script_str,\
                    '-e', qsub_error_path, '-o', qsub_output_path,\
                    '-N', job_n, command])
    print qsub_command
    #raw_input()
    s = subprocess.Popen(qsub_command, shell=True, stdout=subprocess.PIPE)
    s.wait()

def triggerSamfile(bam_path, bed_path, sam_path, contig):
   qsub_error_path = sam_path + '.e'
   qsub_output_path = sam_path + '.o'
   job_name = 'j' + contig + 'getsamfile'
   script_str = '/ifs/scratch/c2b2/ngs_lab/db2175/BashScripts/bed_to_sam.sh'
   command = '   '.join([script_str, bed_path, sam_path, bam_path])
   submit_qsub_job(command, job_name, qsub_error_path, qsub_output_path)

#generates a pileup for each bed file in bed_dir; places output in out_dir
def trigger_pileups(bed_dir, bam_path_model, out_dir):
    for path in os.listdir(bed_dir):
        if path.endswith('.bed'):
            contigname = re.split('\.', path)[0]
            #import pdb; pdb.set_trace()
            bam_path = string.replace(bam_path_model, 'Contig', contigname)
            output_path = out_dir + '/' + contigname + '.pileup'
            qsub_error_path = out_dir + '/' + 'e.' + contigname + '.e'
            qsub_output_path = out_dir + '/' + 'o.' + contigname + '.o'
            job_name = 'j'+ contigname + 'getpileup'
            script_str = '/ifs/scratch/c2b2/ngs_lab/db2175/BashScripts/bed_to_mpileup.sh'
    
            command = '   '.join([script_str, bam_path, REF_PATH,\
                              bed_dir + '/' + path,\
                              output_path])
            #print command
            submit_qsub_job(command, job_name, qsub_error_path, qsub_output_path)
            #raw_input()

# object for handling files in various position-indexed file formats
# files can be loaded into memory or read from disk based on the memory_flag
# accommodates files with one record per genomic position (e.g. pileup) and files with multiple records (e.g., vcf, sam)  
# for files with multiple records it is preferable to iterate through the file; altough records can be accessed with a position tuple (returns a list of records) or a position tuple plus an index (returns a single record)
# (change disk methods so they read a chunk from the file at once instead of one line at a time?)
class GeneFileOb(object):
    def __init__(self, fn, line_type, memory_flag = True,\
                 header_symbol = None, index_flag = False,\
                 one_record_per_pos_flag = True):
        self.memory_flag = memory_flag
        self.one_record_per_pos_flag = one_record_per_pos_flag 
        self.header_symbol = header_symbol
        self.fn = fn
        self.line_type = line_type
        if self.memory_flag:
            self.d = self.get_line_d()
            self.length = len(self.d)
        else:
            self.current_contig = ""
            self.current_pos = -1
            self.current_index = 0
            self.fh = open(self.fn, "rb")
            self.index_flag = index_flag
            if index_flag:
                self.getContigPosD()
            else:
                self.length = None   
    
    def __len__(self):
        if self.length == None:
            self.getContigPosD()
            index.flag = True
        return self.length
    
    def get_all_lines_l(self):
        with open(self.fn, 'rb') as fh:
            if self.header_symbol != None:
                next_line = fh.readline()
                while next_line.startswith(self.header_symbol):
                     next_line = fh.readline()
                return [next_line] + fh.readlines()
            else:
                return fh.readlines()
                
    def get_line_d(self):
        line_d = {} 
        line_l = self.get_all_lines_l()
        if self.one_record_per_pos_flag:
            for line in line_l:
                pos_t = self.line_type(line).getPos_t()
                line_d[pos_t] = self.line_type(line)
        else:
            last_pos_t = ("", 0)
            for line in line_l:
                pos_t = self.line_type(line).getPos_t() 
                if pos_t != last_pos_t:
                    counter = 0
                    line_d[(pos_t, counter)] = self.line_type(line)
                    last_pos_t = pos_t
                elif pos_t == last_pos_t:
                    counter += 1
                    line_d[(pos_t, counter)] = self.line_type(line)
        return line_d 

    def __iter__(self):
        if self.memory_flag:
            return self.memoryLineIter()
        else:
            return self.diskLineIter()

    def memoryLineIter(self):
        line_l = sorted(self.d.values(), cmp = lambda x,y: cmp(x.getPos_t(), y.getPos_t()))
        for l in line_l:
            yield l

    def diskLineIter(self):
        self.fh.seek(0)
        next_line = self.fh.readline()
        if self.header_symbol != None:
            while next_line.startswith(self.header_symbol):
                next_line = self.fh.readline()
        while len(next_line) != 0:
            yield self.line_type(next_line)
            next_line = self.fh.readline()

    def memoryGetItem(self, key):
        if self.one_record_per_pos_flag:
            if key in self.d:
                return self.d[key]
            else:
                return None
        else:
            if len(key[0]) == 1:
                l = []
                counter = 0
                while (key, counter) in self_d:
                    l.append(self.d[(key, counter)])
                    counter += 1
                if len(l) == 0:
                    return None
                else:
                    return l
            else:
                if key in self.d:
                    return self.d[key] 
                else:
                    return None
                    
    def diskGetItem(self, key):
        if self.one_record_per_pos_flag:
            pos_t = key
            if self.current_contig != pos_t[0] or\
              self.current_contig == pos_t[0] and self.current_pos > pos_t[1]:
                if not self.index_flag:
                    self.getContigPosD()
                    self.index_flag = True
                if pos_t[0] in self.contig_pos_d:
                    self.fh.seek(self.contig_pos_d[pos_t[0]])
                    self.current_pos = 0
                    self.current_contig = pos_t[0]
                else:
                    return None
            while self.current_contig  == pos_t[0] and self.current_pos < pos_t[1]:
                self.current_line = self.fh.readline()
                if len(self.current_line) == 0:
                    return None
                (self.current_contig, self.current_pos) = \
                    self.line_type(next_line).getPos_t()
            if (self.current_contig, self.current_pos) == pos:
                return self.line_type(self.current_line)
            else:
                return None
        elif len(key[0]) == 1:
            pos_t = key
            if self.current_contig != pos_t[0] or\
              self.current_contig == pos_t[0] and self.current_pos > pos_t[1] or\
              (self.current_contig, self.current_pos) == pos_t and\
               self.current_index > 0:
                if not self.index_flag:
                    self.getContigPosD()
                    self.index_flag = True
                if pos_t[0] in self.contig_pos_d:
                    self.fh.seek(self.contig_pos_d[pos_t[0]])
                    self.current_pos = 0
                    self.current_contig = pos_t[0]
                else:
                    return None
            while self.current_contig  == pos_t[0] and self.current_pos < pos_t[1]:
                self.current_line = self.fh.readline()
                if len(self.current_line) == 0:
                    return None
                (self.current_contig, self.current_pos) = \
                    self.line_type(self.current_line).getPos_t()
            l = []
            while (self.current_contig, self.current_pos) == pos_t:
                
                l.append(self.line_type(self.current_line))
                self.current_line = self.fh.readline()
                if len(self.current_line) == 0:
                    return l
                (self.current_contig, self.current_pos) = \
                  self.line_type(next_line).getPos_t()
            if len(l) == 0:
                return None
            else:
                return l
        elif len(key[0]) == 2:
            pass
            # TO DO TO DO TO DO    
                    
    def __getitem__(self, key):
        if self.memory_flag:
            return self.memoryGetItem(key)
        else:
            return self.diskGetItem(key)
    # if user tries random access on a file in disk mode, program creates an index first so it knows where the various contigs are
    def getContigPosD(self):
        self.contig_pos_d = {}
        self.fh.seek(0)
        self.current_file_pos = self.fh.tell()
        next_line = self.fh.readline()
        if self.header_symbol != None:
            while next_line.startswith(self.header_symbol):
                self.current_file_pos = self.fh.tell()
                next_line = self.fh.readline()
        (self.current_contig, self.current_pos) = \
                      self.line_type(next_line).getPos_t()
        self.contig_pos_d[self.current_contig] = self.current_file_pos
        current_contig = self.current_contig
        self.current_file_pos = self.fh.tell()
        next_line = self.fh.readline()
        counter = 1
        while len(next_line) != 0:
            counter += 1
            (self.current_contig, self.current_pos) = \
                      self.line_type(next_line).getPos_t()
            if self.current_contig != current_contig:
                self.contig_pos_d[self.current_contig] = \
                   self.current_file_pos
                current_contig = self.current_contig
            self.current_file_pos = self.fh.tell()
            next_line = self.fh.readline()
        self.length = counter   
        self.fh.seek(0)
        self.current_contig = ""
        self.current_pos = 0
# the GeneLineOb objects know how to parse the lines in the various formats
class GeneLineOb(object):
    pass

class VCFLineOb(GeneLineOb):
    def __init__(self, line):
        self.line = line
        self.fields = re.split('\t', line)
    def getPos_t(self):
        return (str(self.fields[0]), int(self.fields[1]))
    def getNormalGenotype(self):
        return re.split(':', self.fields[9])[0]
    def getTumorGenotype(self):
        return re.split(':', self.fields[10])[0]
    def getAltNucleotide(self):
        return self.fields[4]
    def isCoding(self):
        return "exonic" in self.line and\
        "=synonymous" not in self.line
    def isIndel(self):
        return 'INDEL' in self.line

class VCFFileOb(GeneFileOb):
    def __init__(self, fn = None, line_type = VCFLineOb,\
                 memory_flag = True, header_symbol = '#',\
                 index_flag = False, one_record_per_pos_flag = False):
        super(VCFFileOb, self).__init__(fn = fn, line_type = line_type,\
              memory_flag = memory_flag, header_symbol = header_symbol,\
              index_flag = index_flag,\
              one_record_per_pos_flag = one_record_per_pos_flag)
    def __getitem__(self, key):
        return super(VCFFileOb, self).__getitem__(key)

class PILEUPLineOb(GeneLineOb):
    indel_regex = re.compile(r"(-|\+)\d+")
    startread_regex = re.compile(r"\^\S")
    def __init__(self, line):
        self.line = line
        self.fields = re.split('\t', line)
    def getPos_t(self):
        return (str(self.fields[0]), int(self.fields[1]))
    def getReadCount(self):
        return int(self.fields[3])
    def getNucleotideCount(self, nucleotide_char):
        read_data = self.fields[4]
        spans_l = [(0,0)]
        for match in PILEUPLineOb.indel_regex.finditer(read_data):
            spans_l.append( ( match.span()[0], match.span()[1]\
                    + int( read_data[match.span()[0] + 1:match.span()[1]] ) ) )
        spans_l.append((len(read_data), len(read_data)))
        read_str = []
        for i in range(0, len(spans_l)-1):
            read_str.append(read_data[(spans_l[i][1]):(spans_l[i+1][0])])
        read_str = PILEUPLineOb.startread_regex.sub("", "".join(read_str))
        return read_str.lower().count(nucleotide_char.lower())
    def getAltNucleotideCount(self):
        read_data = self.fields[4]
        spans_l = [(0,0)]
        for match in PILEUPLineOb.indel_regex.finditer(read_data):
            spans_l.append( ( match.span()[0], match.span()[1]\
                    + int( read_data[match.span()[0] + 1:match.span()[1]] ) ) )
        spans_l.append((len(read_data), len(read_data)))
        read_str = []
        for i in range(0, len(spans_l)-1):
            read_str.append(read_data[(spans_l[i][1]):(spans_l[i+1][0])])
        read_str = PILEUPLineOb.startread_regex.sub("", "".join(read_str))
        return read_str.lower().count("a") + read_str.lower().count("c") +\
               read_str.lower().count("t") + read_str.lower().count("g")
        
class PILEUPFileOb(GeneFileOb):
    def __init__(self, fn = None, line_type = PILEUPLineOb,\
                 memory_flag = True, header_symbol = None,\
                 one_record_per_pos_flag = True):
        super(PILEUPFileOb, self).__init__(\
            fn = fn, line_type = line_type, memory_flag = memory_flag,\
            header_symbol = None, index_flag = False,\
            one_record_per_pos_flag = one_record_per_pos_flag)
       
class BEDLineOb(GeneLineOb):
    def __init__(self, line):
        self.line = line
        self.fields = re.split('\t', line)
    def getStartPos_t(self):
        return (str(self.fields[0]), int(self.fields[1]))
    def getEndPos_t(self):
        return (str(self.fields[0]), int(self.fields[2]))
    def getPos_t(self):
        return self.getStartPos_t()

class BEDFileOb(GeneFileOb):
    def __init__(self, fn = None, line_type = BEDLineOb,\
                 memory_flag = True, header_symbol = None,\
                 index_flag = False, one_record_per_pos_flag = False):
        super(BEDFileOb, self).__init__(fn = fn, line_type = line_type,\
              memory_flag = memory_flag, header_symbol = header_symbol,\
              index_flag = index_flag,\
              one_record_per_pos_flag = one_record_per_pos_flag)
    def __getitem__(self, key):
        return super(BEDFileOb, self).__getitem__(key)

class SAMLineOb(GeneLineOb):
    CIGAR_regex = re.compile(r"M|I|D|N|S|H|P|=|X")

    def __init__(self, line):
        self.line = line
        self.fields = re.split('\t', line)
    def getPos_t(self):
        return (str(self.fields[2]), int(self.fields[3]))
    def getTemplateLength(self):
        return int(self.fields[8])
    def getCIGARString(self):
        return self.fields[5]
    def getCIGARCount(self):
        CIGARstring = self.getCIGARString()
        counts_l = SAMLineOb.CIGAR_regex.split(CIGARstring)
        if '' in counts_l:
            counts_l.remove('')
        printcounts_l = [int(c) for c in counts_l]
        return sum(printcounts_l)
    def isUnmapped(self):
        return int(self.fields[1]) & 4 != 0
    def isFirstRead(self):
        return int(self.fields[1]) & 64 != 0
    def getTemplateName(self):
        return re.sub('/3$', '', self.fields[0])
    def getMidpoint(self):
        if self.getTemplateLength() > 0:
            return (self.getPos_t()[0],\
                    self.getPos_t()[1] + round(self.getTemplateLength() / float(2), 2))
        else:
            return (self.getPos_t()[0],\
                    self.getPos_t()[1] + self.getCIGARCount() +\
                    round(self.getTemplateLength() / float(2), 2))
    def getPosNext(self):
        return int(self.fields[7])

class SAMFileOb(GeneFileOb):
    def __init__(self, fn = None, line_type = SAMLineOb,\
                 memory_flag = False, header_symbol = None,\
                 index_flag = False, one_record_per_pos_flag = False):
        super(SAMFileOb, self).__init__(fn = fn, line_type = line_type,\
              memory_flag = memory_flag, header_symbol = header_symbol,\
              index_flag = index_flag,\
              one_record_per_pos_flag = one_record_per_pos_flag)
    def __getitem__(self, key):
        return super(SAMFileOb, self).__getitem__(key)

class FAILineOb(GeneLineOb):
    def __init__(self, line):
        self.line = line
        self.fields = re.split('\t', line)
    def getPos_t(self):
        return (str(self.fields[0]), 0)
    def getOffset(self):
        return int(self.fields[2])
    def getContig(self):
        return str(self.fields[0])
    def getLength(self):
        return int(self.fields[1])
    def isNumeric(self):
        return str(self.fields[0]).isdigit()

class FAIFileOb(GeneFileOb):
    def __init__(self, fn = None, line_type = FAILineOb,\
                 memory_flag = True, header_symbol = None,\
                 index_flag = False, one_record_per_pos_flag = False):
        super(FAIFileOb, self).__init__(fn = fn, line_type = line_type,\
              memory_flag = memory_flag, header_symbol = header_symbol,\
              index_flag = index_flag,\
              one_record_per_pos_flag = one_record_per_pos_flag)
    def __getitem__(self, key):
        return super(FAIFileOb, self).__getitem__(key)

def unpickle(pickle_fn):
    pickled_d = pickle.load(open(pickle_fn, "rb"))
    return pickled_d

class MidpointAfterStart(unittest.TestCase):
    def test(self):
        sam = SAMFileOb("/ifs/scratch/c2b2/ngs_lab/db2175/July/July23CNV/temp/s22.sam", memory_flag = False)
        passflag = True
        counter = 0
        for samline in sam:
            if counter % 100 == 6:
                print "counter %s is unmapped %s midpoint %s posnext %s line %s" % \
                (counter, samline.isUnmapped(), samline.getMidpoint(), samline.getPosNext(), samline.line)
            
            if not samline.isUnmapped():
                midpoint = samline.getMidpoint()[1]
                pnext = samline.getPosNext()
                pos = samline.getPos_t()[1]
                if midpoint < min(pnext, pos):
                    print [counter, pos, midpoint, pnext]
                    passflag = False
            counter += 1
            if counter == 1000: break
        self.failUnless(passflag)




if __name__ == "__main__":
    unittest.main()


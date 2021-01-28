#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
A pipeline designed to identify genomic locations of RNA-bound proteins.

"""

import os
import sys
import logging
import subprocess
import argparse
from datetime import datetime
from collections import namedtuple
from itertools import chain
import ruffus

CWD = ''
CORES = 4
DATASET = ''
SPECIES = 'hg19'
RANDOMER_LENGTH = 5
BARCODES_FASTA = ''
BLACKLIST_FILE = ''
GENOME_INDEX = '/storage/vannostrand/reference_data/hg19/genome_star_index'
REPEAT_INDEX = '/storage/vannostrand/reference_data/hg19/repbase_v2_star_index'
SE_READS = []
PE_READS = []


logger = logging.getLogger("CLIP")
logger.setLevel(logging.DEBUG)
console = logging.StreamHandler(sys.stderr)
console.setFormatter(logging.Formatter('%(asctime)s [%(levelname)s] %(message)s', datefmt='%Y-%m-%d %H:%M:%S'))
console.setLevel(logging.DEBUG)
logger.addHandler(console)


def run_cmd(cmd, output_mode='wt', **kwargs):
    """
    Run cmd or throw exception if run fails.
    """
    
    def cmding(cmd):
        cmd = [str(c) for c in cmd if c]
        if '<' in cmd:
            raise ValueError('Invalid cmd, standard input via "<" not supported yet.')
        message = ' '.join(cmd).replace(' -', ' \\\n  -').replace(' >', ' \\\n  >')
        if '>' in cmd:
            output = cmd[cmd.index('>') + 1]
            cmd = cmd[:cmd.index('>')]
        else:
            output = ''
        logger.debug(f'\n{message}')
        return cmd, output
    
    cmd, output = cmding(cmd)
    if output:
        text_mode = True if 't' in output_mode else None
        with open(output, output_mode) as out:
            process = subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE, text=text_mode, **kwargs)
        if process.returncode:
            message = process.stderr or open(output).read()
            raise Exception(f'Failed to run {cmd[0]} (exit code {process.returncode}):\n{message}')
    else:
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, **kwargs)
        stdout, stderr = process.communicate()
        if process.returncode:
            raise Exception(f'Failed to run {cmd[0]} (exit code {process.returncode}):{stderr or stdout}')
        
        
def fastq_list():
    if PE_READS:
        fastqs = chain.from_iterable([(read.fastq1, read.fastq1) for read in PE_READS])
    else:
        fastqs = [read.fastq for read in SE_READS]
    return fastqs


@ruffus.jobs_limit(1)
@ruffus.transform(fastq_list(),
                  ruffus.formatter(r'.+/(?P<BASENAME>.*).f[ast]*q.gz$'),
                  os.path.join(CWD, '{BASENAME[0]}.fastq.gz'))
def soft_link_fastq(fastq, link):
    if fastq == link:
        logger.warning("No symbolic link made. You are directly working on the original data files.")
    else:
        logger.debug(f'Soft link {os.path.basename(fastq)}:\n  ln -s {fastq} {link}')
        os.symlink(fastq, link)


@ruffus.follows(soft_link_fastq)
@ruffus.active_if(len(PE_READS) == 0)
@ruffus.transform([read.fastq for read in SE_READS],
                  ruffus.formatter(r'.+/(?P<BASENAME>.*).fastq.gz$'), '{BASENAME[0]}.umi.fastq.gz')
def extract_umi(fastq, umi_extracted_fastq):
    cmd = ['umi_tools',
           'extract',
           '--random-seed', 1,
           '--stdin', fastq,
           '--bc-pattern', 'NNNNNNNNNN',
           '--log', umi_extracted_fastq.replace('.fastq.gz', '.extract.metrics'),
           '--stdout', umi_extracted_fastq]
    logger.debug(f'Extracting UMIs for {fastq} ...')
    run_cmd(cmd)
    logger.debug(f'Extracting UMIs for {fastq} completed.')


def prepare_demux():
    for fastq1, fastq2, name, barcodes in PE_READS:
        inputs = [os.path.basename(fastq1), os.path.basename(fastq2)]
        outputs = [f'{DATASET}.{name}.{barcodes[0]}.r1.fastq.gz', f'{DATASET}.{name}.{barcodes[0]}.r2.fastq.gz',
                   f'{DATASET}.{name}.{barcodes[1]}.r1.fastq.gz', f'{DATASET}.{name}.{barcodes[1]}.r2.fastq.gz']
        yield inputs, outputs, [name] + barcodes
    

@ruffus.follows(soft_link_fastq)
@ruffus.active_if(len(PE_READS) != 0)
@ruffus.files(prepare_demux)
def demux_reads(fastqs, demuxed_fastqs, barcodes):
    cmd = ['eclipdemux',
           '--fastq1', fastqs[0],
           '--fastq2', fastqs[1],
           '--newname', barcodes[0],
           '--expectedbarcodeida', barcodes[1],
           '--expectedbarcodeidb', barcodes[2],
           '--dataset', DATASET,
           '--barcodesfile', BARCODES_FASTA,
           '--length', RANDOMER_LENGTH]
    logger.debug(f'Demuxing paired reads {os.path.basename(fastqs[0])} and {os.path.basename(fastqs[1])} ...')
    run_cmd(cmd)
    logger.debug(f'Demuxing paired reads {os.path.basename(fastqs[0])} and {os.path.basename(fastqs[1])} completed.')


def prepare_adapter_cut():
    if SE_READS:
        for read in SE_READS:
            yield f'{DATASET}.{read.name}.umi.fastq.gz', f'{DATASET}.{read.name}.umi.trimmed.fastq.gz'
    else:
        for read in PE_READS:
            for barcode in read.barcodes:
                rs = ['r1', 'r2']
                infile = [f'{DATASET}.{read.name}.{barcode}.{r}.fastq.gz' for r in rs]
                outfile = [f'{DATASET}.{read.name}.{barcode}.trimmed.{r}.fastq.gz' for r in rs]
                yield infile, outfile


@ruffus.jobs_limit(1)
@ruffus.files(prepare_adapter_cut)
def cut_adapter(fastq, adapter_trimmed_fastq):
    cmd = ['cutadapt',
           '-j', CORES,
           '-O', 1,
           '--match-read-wildcards',
           '--times', 1,
           '-e', 0.1,
           '--quality-cutoff', 6,
           '-m', 18,
           '-a', 'NNAGATCGGAAGAGC',
           '-a', 'NAGATCGGAAGAGCA',
           '-a', 'AGATCGGAAGAGCAC',
           '-a', 'GATCGGAAGAGCACA',
           '-a', 'ATCGGAAGAGCACAC',
           '-a', 'TCGGAAGAGCACACG',
           '-a', 'CGGAAGAGCACACGT',
           '-a', 'GGAAGAGCACACGTC',
           '-a', 'GAAGAGCACACGTCT',
           '-a', 'AAGAGCACACGTCTG',
           '-a', 'AGAGCACACGTCTGA',
           '-a', 'GAGCACACGTCTGAA',
           '-a', 'AGCACACGTCTGAAC',
           '-a', 'GCACACGTCTGAACT',
           '-a', 'CACACGTCTGAACTC',
           '-a', 'ACACGTCTGAACTCC',
           '-a', 'CACGTCTGAACTCCA',
           '-a', 'ACGTCTGAACTCCAG',
           '-a', 'CGTCTGAACTCCAGT',
           '-a', 'GTCTGAACTCCAGTC',
           '-a', 'TCTGAACTCCAGTCA',
           '-a', 'CTGAACTCCAGTCAC']
    if isinstance(fastq, str):
        metrics = adapter_trimmed_fastq.replace('.trimmed.fastq.gz', '.cutadapt.metrics')
        args = ['-o', adapter_trimmed_fastq, fastq, '>', metrics]
        message = f'single read {fastq}'
    else:
        fastq1, fastq2 = fastq
        trimmed1, trimmed2 = adapter_trimmed_fastq
        metrics = trimmed1.replace('.r1.trimmed.fastq.gz', '.cutadapt.metrics')
        args = ['-o', trimmed1, '-p', trimmed2, fastq1, fastq2, '>', metrics]
        message = f'paired reads {fastq1} and {fastq2}'
    cmd.extend(args)
    
    logger.debug(f'Trimming adapters for {message} ...')
    run_cmd(cmd)
    logger.debug(f'Trimming adapters for {message} completed.')


def prepare_adapter_cut_again():
    if SE_READS:
        for read in SE_READS:
            yield f'{DATASET}.{read.name}.umi.trimmed.fastq.gz', f'{DATASET}.{read.name}.umi.trimmed.trimmed.fastq'
    else:
        for read in PE_READS:
            for barcode in read.barcodes:
                rs = ['r1', 'r2']
                infile = [f'{DATASET}.{read.name}.{barcode}.trimmed.{r}.fastq.gz' for r in rs]
                outfile = [f'{DATASET}.{read.name}.{barcode}.trimmed.trimmed.{r}.fastq' for r in rs]
                yield infile, outfile
                
                
@ruffus.jobs_limit(1)
@ruffus.files(prepare_adapter_cut_again)
@ruffus.follows(cut_adapter)
def cut_adapter_again(fastq, adapter_trimmed_fastq):
    cmd = ['cutadapt',
           '-j', CORES,
           '-O', 1,
           '--match-read-wildcards',
           '--times', 1,
           '-e', 0.1,
           '--quality-cutoff', 6,
           '-m', 18,
           '-a', 'NNAGATCGGAAGAGC',
           '-a', 'NAGATCGGAAGAGCA',
           '-a', 'AGATCGGAAGAGCAC',
           '-a', 'GATCGGAAGAGCACA',
           '-a', 'ATCGGAAGAGCACAC',
           '-a', 'TCGGAAGAGCACACG',
           '-a', 'CGGAAGAGCACACGT',
           '-a', 'GGAAGAGCACACGTC',
           '-a', 'GAAGAGCACACGTCT',
           '-a', 'AAGAGCACACGTCTG',
           '-a', 'AGAGCACACGTCTGA',
           '-a', 'GAGCACACGTCTGAA',
           '-a', 'AGCACACGTCTGAAC',
           '-a', 'GCACACGTCTGAACT',
           '-a', 'CACACGTCTGAACTC',
           '-a', 'ACACGTCTGAACTCC',
           '-a', 'CACGTCTGAACTCCA',
           '-a', 'ACGTCTGAACTCCAG',
           '-a', 'CGTCTGAACTCCAGT',
           '-a', 'GTCTGAACTCCAGTC',
           '-a', 'TCTGAACTCCAGTCA',
           '-a', 'CTGAACTCCAGTCAC']
    if isinstance(fastq, str):
        metrics = adapter_trimmed_fastq.replace('.trimmed.fastq.gz', '.cutadapt.metrics')
        args = ['-o', adapter_trimmed_fastq, fastq, '>', metrics]
        message = f'single read {fastq}'
    else:
        fastq1, fastq2 = fastq
        trimmed1, trimmed2 = adapter_trimmed_fastq
        metrics = trimmed1.replace('.r1.trimmed.fastq.gz', '.cutadapt.metrics')
        args = ['-o', trimmed1, '-p', trimmed2, fastq1, fastq2, '>', metrics]
        message = f'paired reads {fastq1} and {fastq2}'
    cmd.extend(args)
    
    logger.debug(f'Trimming adapters second round for {message} ...')
    run_cmd(cmd)
    logger.debug(f'Trimming adapters second round for {message} completed.')


def prepare_fastq_sort():
    if SE_READS:
        for read in SE_READS:
            infile = f'{DATASET}.{read.name}.umi.trimmed.trimmed.fastq'
            outfile = f'{DATASET}.{read.name}.umi.trimmed.trimmed.sorted.fastq'
            yield infile, outfile
    else:
        for read in PE_READS:
            for barcode in read.barcodes:
                rs = ['r1', 'r2']
                infile = [f'{DATASET}.{read.name}.{barcode}.trimmed.trimmed.{r}.fastq' for r in rs]
                outfile = [f'{DATASET}.{read.name}.{barcode}.trimmed.trimmed.sorted.{r}.fastq' for r in rs]
                yield infile, outfile


@ruffus.files(prepare_fastq_sort)
@ruffus.follows(cut_adapter_again)
def sort_fastq(fastq, sorted_fastq):
    cmd = ['fastq-sort', '--id', fastq, '>', sorted_fastq]
    logger.debug(f'Sorting adapters trimmed fastq for {fastq} ...')
    run_cmd(cmd)
    logger.debug(f'Sorting adapters trimmed fastq for {fastq} completed.')


def prepare_repeat_element_map():
    if SE_READS:
        for read in SE_READS:
            infile = f'{DATASET}.{read.name}.umi.trimmed.trimmed.sorted.fastq'
            outfile = f'{DATASET}.{read.name}.repeat.element.Unmapped.out.mate1'
            yield infile, outfile
    else:
        for read in PE_READS:
            for barcode in read.barcodes:
                infile = [f'{DATASET}.{read.name}.{barcode}.trimmed.trimmed.{r}.sorted.fastq'
                          for r in ['r1', 'r2']]
                outfile = [f'{DATASET}.{read.name}.{barcode}.repeat.elements.Unmapped.out.mate1',
                           f'{DATASET}.{read.name}.{barcode}.repeat.elements.Unmapped.out.mate2']
                yield infile, outfile

                
@ruffus.jobs_limit(1)
@ruffus.follows(sort_fastq)
@ruffus.files(prepare_repeat_element_map)
def repeat_map(sorted_fastq, mate):
    cmd = ['STAR',
           '--runMode', 'alignReads',
           '--runThreadN', CORES,
           '--alignEndsType', 'EndToEnd',
           '--genomeDir', REPEAT_INDEX,
           '--genomeLoad', 'NoSharedMemory',
           '--outBAMcompression', 10,
           '--outFilterMultimapNmax', 30,
           '--outFilterMultimapScoreRange', 1,
           '--outFilterScoreMin', 10,
           '--outFilterType', 'BySJout',
           '--outReadsUnmapped', 'Fastx',
           '--outSAMattrRGline', 'ID:foo',
           '--outSAMattributes', 'All',
           '--outSAMmode', 'Full',
           '--outSAMtype', 'BAM', 'Unsorted',
           '--outSAMunmapped', 'Within',
           '--outStd', 'Log',
           '--outFileNamePrefix', mate.replace('repeat.elements.Unmapped.out.mate1', ''),
           '--readFilesIn', sorted_fastq]
    if isinstance(sorted_fastq, str):
        cmd.extend(['--outFileNamePrefix', mate.replace('repeat.elements.Unmapped.out.mate1', '')])
        cmd.extend(['--readFilesIn', sorted_fastq])
        message = f'single read {sorted_fastq}'
    else:
        fastq1, fastq2 = sorted_fastq
        cmd.extend(['--outFileNamePrefix', mate[0].replace('repeat.elements.Unmapped.out.mate1', '')])
        cmd.extend(['--readFilesIn', fastq1, fastq2])
        message = f'paired reads {fastq1} and {fastq2}'
    logger.debug(f'Mapping {message} to repeat elements ...')
    run_cmd(cmd)
    logger.debug(f'Mapping {message} to repeat elements completed.')


def prepare_reference_genome_map():
    if SE_READS:
        for read in SE_READS:
            outfile = f'{DATASET}.{read.name}.repeat.elements.Unmapped.out.mate1'
            infile = f'{DATASET}.{read.name}.genome.map.Aligned.out.bam'
            yield infile, outfile
    else:
        for read in PE_READS:
            for barcode in read.barcodes:
                infile = [f'{DATASET}.{read.name}.{barcode}.repeat.elements.Unmapped.out.mate1',
                          f'{DATASET}.{read.name}.{barcode}.repeat.elements.Unmapped.out.mate2']
                outfile = f'{DATASET}.{read.name}.{barcode}.genome.map.Aligned.out.bam'
                yield infile, outfile
                
                
@ruffus.jobs_limit(1)
@ruffus.files(prepare_reference_genome_map())
@ruffus.follows(repeat_map)
def genomic_map(mate, bam):
    cmd = ['STAR',
           '--runMode', 'alignReads',
           '--runThreadN', CORES,
           '--alignEndsType', 'EndToEnd',
           '--genomeDir', GENOME_INDEX,
           '--genomeLoad', 'NoSharedMemory',
           '--outBAMcompression', 10,
           '--outFileNamePrefix', bam.replace('genome.map.Aligned.out.bam', ''),
           '--outFilterMultimapNmax', 1,
           '--outFilterMultimapScoreRange', 1,
           '--outFilterScoreMin', 10,
           '--outFilterType', 'BySJout',
           '--outReadsUnmapped', 'Fastx',
           '--outSAMattrRGline', 'ID:foo',
           '--outSAMattributes', 'All',
           '--outSAMmode', 'Full',
           '--outSAMtype', 'BAM', 'Unsorted',
           '--outSAMunmapped', 'Within',
           '--outStd', 'Log',
           '--readFilesIn', mate]
    if isinstance(mate, str):
        cmd.extend(['--outFileNamePrefix', bam.replace('repeat.elements.Unmapped.out.mate1', '')])
        cmd.extend(['--readFilesIn', mate])
        message = f'unmapped mate {mate}'
    else:
        fastq1, fastq2 = mate
        cmd.extend(['--outFileNamePrefix', bam.replace('repeat.elements.Unmapped.out.mate1', '')])
        cmd.extend(['--readFilesIn', fastq1, fastq2])
        message = f'paired unmapped mates {fastq1} and {fastq2}'
    logger.debug(f'Mapping {message} to genome ...')
    run_cmd(cmd)
    logger.debug(f'Mapping {message} to genome completed.')
    
    
def sort_bam(bam, sorted_bam, name=False):
    option = '-n' if name else ''
    cmd = ['samtools', 'sort', '-@', CORES, option, '-m', '4G', '-o', sorted_bam, bam]
    logger.debug(f'Name sorting BAM {os.path.basename(bam)} ...')
    run_cmd(cmd)
    logger.debug(f'Name sorting BAM {os.path.basename(bam)} completed.')
    return sorted_bam


# SE: name sort -> pos sort -> umi_tools dedup -> pos sort -> index
# PE: name sort -> barcode collapse -> pos sort r2_bam -> index


def barcode_collapse(bam):
    dedup_bam = bam.replace('.genome.map.name.sorted.bam', '.genome.map.name.sorted.deduped.bam')
    cmd = ['barcodecollapsepe.py',
           '-b', bam,
           '-o', dedup_bam,
           '-m', bam.replace('.genome.map.name.sorted.bam', '.genome.map.name.sorted.dedup.metrics')]
    logger.debug(f'Collapsing barcode for {bam} ...')
    run_cmd(cmd)
    logger.debug(f'Collapsing barcode for {bam} completed.')
    return dedup_bam


def umi_dedup(bam):
    deduped_bam = bam.replace('.sorted.bam', '.sorted.deduped.bam')
    cmd = ['umi_tools', 'dedup',
           '--random-seed', 1,
           '-I', bam,
           '--method', 'unique',
           '--output-stats', deduped_bam.replace('.deduped.bam', '.dedup'),
           '--log', deduped_bam.replace('.deduped.bam', '.dedup.metrics'),
           '-S', deduped_bam]
    logger.debug(f'Deduplcating BAM {os.path.basename(bam)} ...')
    run_cmd(cmd)
    logger.debug(f'Deduplcating BAM {os.path.basename(bam)} completed.')
    return deduped_bam


def view_read2(bam):
    r2_bam = bam.replace('.deduped.bam', '.r2.deduped.bam')
    cmd = ['samtools', 'view', '-@', CORES, '-f', 128, '-b', '-o', r2_bam, bam]
    logger.debug(f'Extracting read 2 from {bam} ...')
    run_cmd(cmd)
    logger.debug(f'Extracting read 2 from {bam} completed.')
    return r2_bam


def index_bam(bam, bai):
    cmd = ['samtools', 'index', '-@', CORES, bam, bai]
    logger.debug(f'Indexing BAM {os.path.basename(bam)} ...')
    run_cmd(cmd)
    logger.debug(f'Indexing BAM {os.path.basename(bam)} completed.')


@ruffus.transform(genomic_map, ruffus.suffix('.genome.map.Aligned.out.bam'),
                  ['*.deduped.bam', '*.deduped.pos.sorted.bam'])
def dedup_bam(bam, deduped_bam):
    sorted_bam = sort_bam(bam, bam.replace('.genome.map.Aligned.out.bam', '.genome.map.name.sorted.bam'))
    if PE_READS:
        deduped_bam = barcode_collapse(sorted_bam)
        ds_bam = sort_bam(deduped_bam, deduped_bam.replace('.deduped.bam', '.pos.sorted.bam'), name=False)
        r2_bam = view_read2(ds_bam)
    else:
        sorted_bam = sort_bam(sorted_bam, sorted_bam.replace('.sorted.bam', '.sorted.pos.sorted.bam'), name=False)
        deduped_bam = umi_dedup(sorted_bam)
        r2_bam = sort_bam(deduped_bam, deduped_bam.replace('.deduped.bam', '.deduped.pos.sorted.bam'), name=False)
    index_bam(r2_bam, f'{r2_bam}.bai')


@ruffus.transform(dedup_bam, ['*.deduped.bam', '*.deduped.pos.sorted.bam'],
                  ['.name.sorted.pos.sorted.deduped.sorted.positive.bw',
                   '.name.sorted.pos.sorted.deduped.sorted.negative.bw'])
def make_bigwig_files(bam, bigwig):
    cmd = ['makebigwigfiles',
           '--bw_pos', bam.replace('.sorted.bam', 'sorted.positive.bw'),
           '--bw_neg', bam.replace('.sorted.bam', 'sorted.negative.bw'),
           '--bam', bam,
           '--genome', os.path.join(GENOME_INDEX, 'chrNameLength.txt'),
           '--direction', 'f']
    logger.debug(f'Making BigWig files using BAM {os.path.basename(bam)} ...')
    run_cmd(cmd)
    logger.debug(f'Making BigWig files using BAM {os.path.basename(bam)} completed.')


@ruffus.follows(make_bigwig_files)
@ruffus.transform(dedup_bam, ruffus.suffix('.deduped.sorted.bam'), '.peak.clusters.bed')
def clipper(bam, bed):
    cmd = ['clipper',
           '--species', SPECIES,
           '--bam', bam,
           '--outfile', bed]
    logger.debug(f'Identifying peaks using BAM {os.path.basename(bam)} ...')
    run_cmd(cmd)
    logger.debug(f'Identifying peaks using BAM {os.path.basename(bam)} completed.')


def set_se_reads(fastqs, names, adapters):
    read = namedtuple('Read', 'fastq name adapter')
    reads = [read(fastq, name, adapter) for fastq, name, adapter in zip(fastqs, names, adapters)]
    return reads


def set_pe_reads(fastq1s, fastq2s, names, barcodes):
    read = namedtuple('Read', 'fastq1 fastq2 name barcodes')
    reads = [read(fastq1, fastq2, name, barcode.split(';'))
             for fastq1, fastq2, name, barcode in zip(fastq1s, fastq2s, names, barcodes)]
    return reads


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--fastq1', type=str, dest='fastq1',
                        help='Path to a fastq file or comma separated paths to fastq files list for read 1.')
    parser.add_argument('--fastq2', type=str, dest='fastq2',
                        help='Path to a fastq file or comma separated paths to fastq files list for read 2, PE only')
    parser.add_argument('--name', type=str, dest='name',
                        help='A string or comma separated strings for customizing outfile name field.')
    parser.add_argument('--adapters', type=str, dest='adapters',
                        help='Path to the adapters file, SE only.')
    parser.add_argument('--blacklist-file', type=str, dest='blacklist_file',
                        help='Path to the blacklist file, SE only.')
    parser.add_argument('--barcodes', type=str, dest='barcodes',
                        help='Comma and semi-colon separated barcode IDs, PE only.')
    parser.add_argument('--barcodes-fasta', type=str, dest='barcodes_fasta',
                        help='Path to FASTA format file for barcode IDs and sequences, PE only')
    parser.add_argument('--randomer-length', type=int, dest='randomer_length',
                        help='Integer number for the length of randomer, PE only')
    parser.add_argument('--dataset', type=str, dest='dataset',
                        help='A string or comma separated strings for customizing outfile dataset field.')
    parser.add_argument('--species', type=str, dest='species',
                        help='Short name of species for the dataset.')
    parser.add_argument('--genome', type=str, dest='genome',
                        help="Path to the STAR reference genome index directory.")
    parser.add_argument('--repeat', type=str, dest='repeat',
                        help="Path to the STAR repeat elements index directory.")
    parser.add_argument('--outdir', type=str, dest='outdir',
                        help="Path to a directory for saving output.")
    parser.add_argument('-cores', type=int, dest='cores',
                        help='Number of cores can be used for running pipeline.')
    parser.add_argument('--dry-run', action='store_true', dest='dry_run',
                        help='Only print out the details of each step without actually run the pipeline.')
    
    args = parser.parse_args()
    global CWD
    CWD = args.outdir
    os.chdir(CWD)
    
    global logger
    handler = logging.FileHandler(filename=f'clip_{datetime.today().strftime("%Y-%m-%d_%H:%M:%S")}.log', mode='w')
    handler.setFormatter(logging.Formatter('%(asctime)s [%(levelname)s] %(message)s', datefmt='%Y-%m-%d %H:%M:%S'))
    handler.setLevel(logging.DEBUG)
    logger.addHandler(handler)
    
    global CORES
    CORES = args.cores if args.cores else CORES
    
    global GENOME_INDEX
    GENOME_INDEX = args.genome if args.genome else GENOME_INDEX
    
    global REPEAT_INDEX
    REPEAT_INDEX = args.repeat if args.repeat else REPEAT_INDEX
    
    global DATASET
    DATASET = args.dataset if args.dataset else DATASET
    
    global SPECIES
    SPECIES = args.species if args.species else SPECIES
    
    global BARCODES_FASTA
    BARCODES_FASTA = args.barcodes_fasta if args.barcodes_fasta else BARCODES_FASTA
    
    global RANDOMER_LENGTH
    RANDOMER_LENGTH = args.randomer_length if args.randomer_length else RANDOMER_LENGTH
    
    global BLACKLIST_FILE
    BLACKLIST_FILE = args.blacklist_file if args.blacklist_file else BLACKLIST_FILE
    
    names = args.name.split(',')
    adapters = args.adapters.split(',')
    barcodes = args.barcodes.split(',')
    
    global SE_READS
    fastq1s = args.fastq1.split(',') if args.fastq1 else []
    SE_READS = set_se_reads(fastq1s, names, adapters)
    
    global PE_READS
    
    if args.fastq2s:
        fastq2s = args.fastq2.split(',')
        PE_READS = set_pe_reads(fastq1s, fastq2s, names, barcodes)
        
    if args.dry_run:
        ruffus.pipeline_printout()
    else:
        processes = 1 if CORES == 1 else int(CORES / 2)
        ruffus.pipeline_run(multiprocess=processes, verbose=1)
    
    
if __name__ == '__main__':
    main()

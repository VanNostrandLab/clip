#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
A pipeline designed to identify genomic locations of RNA-bound proteins.

"""

import os
import sys
import logging
import subprocess
import glob

import ruffus


READS_MODE = 'paired-end'
CORES = 4
GENOME_INDEX = '/storage/vannostrand/reference_data/hg19/genome_star_index'
REPEAT_INDEX = '/storage/vannostrand/reference_data/hg19/repbase_v2_star_index'
SPECIES = 'hg19'
os.chdir(os.path.dirname(os.path.abspath(__file__)))

logger = logging.getLogger("CLIP")
logger.setLevel(logging.DEBUG)
handler = logging.StreamHandler(sys.stderr)
handler.setFormatter(logging.Formatter('%(asctime)s [%(levelname)s] %(message)s', datefmt='%Y-%m-%d %H:%M:%S'))
handler.setLevel(logging.DEBUG)
logger.addHandler(handler)
note = logging.FileHandler(filename='chimeras.log', mode='w')
note.setFormatter(logging.Formatter('%(asctime)s [%(levelname)s] %(message)s', datefmt='%Y-%m-%d %H:%M:%S'))
note.setLevel(logging.DEBUG)
logger.addHandler(note)

# fastqs = ['1_R1_001.fastq.gz', '1_R2_001.fastq.gz', '2_R1_001.fastq.gz', '2_R2_001.fastq.gz']
fastqs = [('1_R1_001.fastq.gz', '1_R2_001.fastq.gz'), ('2_R1_001.fastq.gz', '2_R2_001.fastq.gz')]


def run_cmd(cmd, output_mode='wt', **kwargs):
    """
    Run cmd or throw exception if run fails.
    """
    
    def cmding(cmd):
        cmd = [str(c) for c in cmd]
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


@ruffus.jobs_limit(1)
@ruffus.transform(fastqs,
                  ruffus.formatter(r'.+/(?P<BASENAME>.*).f[ast]*q.gz$',
                                   r'.+/(?P<BASENAME>.*).f[ast]*q.gz$'),
                  [os.path.join(os.getcwd(), '{BASENAME[0]}.fastq.gz'),
                   os.path.join(os.getcwd(), '{BASENAME[1]}.fastq.gz')])
def pe_soft_link_rawdata_to_data_directory(fastqs, links):
    for fastq, link in zip(fastqs, links):
        fastq = os.path.abspath(fastq)
        if fastq == link:
            logger.warning("No symbolic link made. You are directly working on the original data files.")
        else:
            logger.debug(f'Soft link {os.path.basename(fastq)}:\n  ln -s {fastq} {link}')
            os.symlink(fastq, link)


@ruffus.jobs_limit(1)
@ruffus.active_if(READS_MODE == 'single-end')
@ruffus.transform(fastqs,
                  ruffus.formatter(r'.+/(?P<BASENAME>.*).f[ast]*q.gz$'),
                  os.path.join(os.getcwd(), '{BASENAME[0]}.fastq.gz'))
def se_soft_link_rawdata_to_data_directory(fastq, link):
    if fastq == link:
        logger.warning("No symbolic link made. You are directly working on the original data files.")
    else:
        logger.debug(f'Soft link {os.path.basename(fastq)}:\n  ln -s {fastq} {link}')
        os.symlink(fastq, link)


@ruffus.transform(se_soft_link_rawdata_to_data_directory,
                  ruffus.formatter(r'.+/(?P<BASENAME>.*).fastq.gz$'),
                  os.path.join(os.getcwd(), '{BASENAME[0]}.umi.fastq.gz'))
def extract_umi(fastq, umi_extracted_fastq):
    cmd = ['umi_tools',
           'extract',
           '--random-seed', 1,
           '--stdin', fastq,
           '--bc-pattern', 'NNNNNNNNNN',
           '--log', umi_extracted_fastq.replace('.fastq.gz', '.extract.metrics'),
           '--stdout', umi_extracted_fastq]
    logger.debug(f'Extracting UMIs for {os.path.basename(fastq)} ...')
    run_cmd(cmd)
    logger.debug(f'Extracting UMIs for {os.path.basename(fastq)} completed.')


barcodes, barcodes_file, length = [('NIL', 'NIL'), ('0B1', 'D3f')], ['', ''], [5, 5]
newname, dataset = [('clip1', 'clip2')], [('dataset1', 'dataset2')]


def generate_demux_parameters():
    pass


@ruffus.transform(pe_soft_link_rawdata_to_data_directory,
                  ruffus.formatter(r'.+/(?P<BASENAME>.*).fastq.gz$',
                                   r'.+/(?P<BASENAME>.*).fastq.gz$'),
                  [os.path.join(os.getcwd(), '{BASENAME[0]}.*.fastq.gz'),
                   os.path.join(os.getcwd(), '{BASENAME[1]}.*.fastq.gz')],
                   barcodes, barcodes_file, length, newname, dataset)
def demux_reads(fastq, demuxed_fastq, barcodes, barcodes_file, length, newname, dataset):
    fastq1, fastq2 = fastq
    print(fastq, demuxed_fastq, barcodes, barcodes_file, length, newname, dataset)
    cmd = ['eclipdemux',
           '--metrics', '',
           '--expectedbarcodeida', barcodes[0],
           '--expectedbarcodeidb', barcodes[1],
           '--fastq1', fastq1,
           '--fastq2', fastq2,
           '--newname', newname,
           '--dataset', dataset,
           '--barcodesfile', barcodes_file,
           '--length', length]
    print(cmd)


@ruffus.jobs_limit(1)
@ruffus.transform(extract_umi if READS_MODE == 'single-end' else demux_reads,
                  ruffus.suffix('.fastq.gz'), '.trimmed.fastq.gz')
def cut_adapter(fastq, adapter_trimmed_fastq):
    cmd = ['cutadapt',
           '-j', CORES,
           '-O', 1,
           '--match-read-wildcards',
           '--times', 1,
           '-e', 0.1,
           '--quality-cutoff', 6,
           '-m', 18,
           '-o', adapter_trimmed_fastq,
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
           '-a', 'CTGAACTCCAGTCAC',
           fastq, '>', adapter_trimmed_fastq.replace('.trimmed.fastq.gz', '.cutadapt.metrics')]
    logger.debug(f'Trimming adapters for {os.path.basename(fastq)} ...')
    run_cmd(cmd)
    logger.debug(f'Trimming adapters {os.path.basename(fastq)} completed.')


@ruffus.jobs_limit(1)
@ruffus.transform(cut_adapter, ruffus.suffix('.trimmed.fastq.gz'), '.trimmed.trimmed.fastq')
def cut_adapter_again(fastq, adapter_trimmed_fastq):
    cmd = ['cutadapt',
           '-j', CORES,
           '-O', 1,
           '--match-read-wildcards',
           '--times', 1,
           '-e', 0.1,
           '--quality-cutoff', 6,
           '-m', 18,
           '-o', adapter_trimmed_fastq,
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
           '-a', 'CTGAACTCCAGTCAC',
           fastq, '>', adapter_trimmed_fastq.replace('.trimmed.fastq', '.cutadapt.metrics')]
    logger.debug(f'Trimming adapters for {os.path.basename(fastq)} ...')
    run_cmd(cmd)
    logger.debug(f'Trimming adapters {os.path.basename(fastq)} completed.')


@ruffus.transform(cut_adapter_again, ruffus.suffix('.trimmed.fastq'), '.trimmed.sorted.fastq')
def sort_fastq(fastq, sorted_fastq):
    cmd = ['fastq-sort', '--id', fastq, '>', sorted_fastq]
    logger.debug(f'Sorting adapters trimmed fastq for {os.path.basename(fastq)} ...')
    run_cmd(cmd)
    logger.debug(f'Sorting adapters trimmed fastq for {os.path.basename(fastq)} completed.')


@ruffus.jobs_limit(1)
@ruffus.transform(sort_fastq, ruffus.suffix('.trimmed.trimmed.sorted.fastq'), '.repeat.elements.Unmapped.out.mate1')
def repeat_map(sorted_fastq, mate):
    cmd = ['STAR',
           '--runMode', 'alignReads',
           '--runThreadN', CORES,
           '--alignEndsType', 'EndToEnd',
           '--genomeDir', REPEAT_INDEX,
           '--genomeLoad', 'NoSharedMemory',
           '--outBAMcompression', 10,
           '--outFileNamePrefix', mate.replace('repeat.elements.Unmapped.out.mate1', ''),
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
           '--readFilesIn', sorted_fastq]
    logger.debug(f'Mapping reads to repeat elements for {os.path.basename(sorted_fastq)} ...')
    run_cmd(cmd)
    logger.debug(f'Mapping reads to repeat elements for {os.path.basename(sorted_fastq)} completed.')


@ruffus.jobs_limit(1)
@ruffus.transform(repeat_map, ruffus.suffix('.repeat.elements.Unmapped.out.mate1'), '.genome.map.Aligned.out.bam')
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
    logger.debug(f'Mapping reads to reference genome for {os.path.basename(mate)} ...')
    run_cmd(cmd)
    logger.debug(f'Mapping reads to reference genome for {os.path.basename(mate)} completed.')


@ruffus.jobs_limit(1)
@ruffus.transform(genomic_map, ruffus.suffix('.Aligned.out.bam'), '.name.sorted.bam')
def name_sort_bam(bam, sorted_bam):
    cmd = ['samtools', 'sort', '-@', CORES, '-n', '-o', sorted_bam, bam]
    logger.debug(f'Name sorting BAM {os.path.basename(bam)} ...')
    run_cmd(cmd)
    logger.debug(f'Name sorting BAM {os.path.basename(bam)} completed.')


@ruffus.jobs_limit(1)
@ruffus.transform(name_sort_bam, ruffus.suffix('.name.sorted.bam'), '.name.sorted.pos.sorted.bam')
def pos_sort_bam(bam, sorted_bam):
    cmd = ['samtools', 'sort', '-@', CORES, '-o', sorted_bam, bam]
    logger.debug(f'Position sorting BAM {os.path.basename(bam)} ...')
    run_cmd(cmd)
    logger.debug(f'Position sorting BAM {os.path.basename(bam)} completed.')


@ruffus.jobs_limit(1)
@ruffus.transform(pos_sort_bam, ruffus.suffix('.name.sorted.pos.sorted.bam'), '.name.sorted.pos.sorted.bam.bai')
def index_bam(bam, bai):
    cmd = ['samtools', 'index', '-@', CORES, bam, bai]
    logger.debug(f'Indexing BAM {os.path.basename(bam)} ...')
    run_cmd(cmd)
    logger.debug(f'Indexing BAM {os.path.basename(bam)} completed.')


@ruffus.follows(index_bam)
@ruffus.transform(pos_sort_bam, ruffus.suffix('.name.sorted.pos.sorted.bam'), '.name.sorted.pos.sorted.deduped.bam')
def dedup_bam(bam, deduped_bam):
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


@ruffus.jobs_limit(1)
@ruffus.transform(dedup_bam, ruffus.suffix('.name.sorted.pos.sorted.deduped.bam'),
                  '.name.sorted.pos.sorted.deduped.sorted.bam')
def sort_deduped_bam(bam, sorted_bam):
    cmd = ['samtools', 'sort', '-@', CORES, '-o', sorted_bam, bam]
    logger.debug(f'Sorting deduped BAM {os.path.basename(bam)} ...')
    run_cmd(cmd)
    logger.debug(f'Sorting deduped BAM {os.path.basename(bam)} completed.')


@ruffus.transform(sort_deduped_bam, ruffus.suffix('.name.sorted.pos.sorted.deduped.sorted.bam'),
                  '.name.sorted.pos.sorted.deduped.sorted.bam.bai')
def index_dedup_sorted_bam(bam, bai):
    cmd = ['samtools', 'index', '-@', CORES, bam, bai]
    logger.debug(f'Indexing deduped sorted BAM {os.path.basename(bam)} ...')
    run_cmd(cmd)
    logger.debug(f'Indexing deduped sorted BAM {os.path.basename(bam)} completed.')


@ruffus.follows(index_dedup_sorted_bam)
@ruffus.transform(sort_deduped_bam, ruffus.suffix('.name.sorted.pos.sorted.deduped.sorted.bam'),
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
@ruffus.transform(sort_deduped_bam, ruffus.suffix('.deduped.sorted.bam'), '.peak.clusters.bed')
def clipper(bam, bed):
    cmd = ['clipper',
           '--species', SPECIES,
           '--bam', bam,
           '--outfile', bed]
    logger.debug(f'Identifying peaks using BAM {os.path.basename(bam)} ...')
    run_cmd(cmd)
    logger.debug(f'Identifying peaks using BAM {os.path.basename(bam)} completed.')


if __name__ == '__main__':
    pass
    ruffus.pipeline_run(multiprocess=1, verbose=1)
    # ruffus.pipeline_printout()

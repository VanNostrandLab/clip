#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
A pipeline designed to identify genomic locations of RNA-bound proteins.

"""

import os
import sys
import subprocess
import argparse
import json
import logging
from datetime import datetime

import ruffus
from ruamel.yaml import YAML


class Sample:
    def __init__(self, ip_read, input_read):
        self.ip_read = ip_read
        self.input_read = input_read
        
        
class Read:
    def __init__(self, name, r1, r2='', adapters='', barcodes=None, barcode_fasta='',
                 dataset='', blacklist_file='', randomer_length=5):
        self.name = name
        self.r1 = r1
        self.r2 = r2
        self.r1_link = f'{dataset}_{name}_R1.fastq.gz'
        self.r2_link = f'{dataset}_{name}_R2.fastq.gz'
        
        if r2:
            self.type = 'paired'
            assert barcodes, f'Paired-End sample {name} does not have barcodes provided.'
            assert len(barcodes) == 2, f'Invalid barcodes {barcodes} for sample {name}.'
            assert barcode_fasta, f'Paired-End sample {name} does not have barcode fasta file provided.'
            self.barcodes = barcodes
            self.barcode_fasta = barcode_fasta
            self.randomer_length = randomer_length
            b1, b2 = barcodes
            if b1 == b2:
                self.r1_umi = [f'{dataset}.{name}.{b1}.R1.fastq.gz']
                self.r2_umi = [f'{dataset}.{name}.{b1}.R2.fastq.gz']
                self.r1_trim = [f'{dataset}.{name}.{b1}.trim.R1.fastq.gz']
                self.r2_trim = [f'{dataset}.{name}.{b1}.trim.R2.fastq.gz']
                self.r1_sort = [f'{dataset}.{name}.{b1}.trim.sort.R1.fastq.gz']
                self.r2_sort = [f'{dataset}.{name}.{b1}.trim.sort.R2.fastq.gz']
                self.repeat_map_prefix = [f'{dataset}.{name}.{b1}.repeat.map.']
                self.mate1 = [f'{dataset}.{name}.{b1}.repeat.map.Unmapped.out.mate1']
                self.mate2 = [f'{dataset}.{name}.{b1}.repeat.map.Unmapped.out.mate2']
                self.mate1_sort = [f'{dataset}.{name}.{b1}.repeat.map.Unmapped.sort.mate1']
                self.mate2_sort = [f'{dataset}.{name}.{b1}.repeat.map.Unmapped.sort.mate2']
                self.genome_map_prefix = [f'{dataset}.{name}.{b1}.genome.map.']
                self.bam = [f'{dataset}.{name}.{b1}.genome.map.Aligned.out.bam']
                self.bam_dedup = [f'{dataset}.{name}.{b1}.dedup.bam']
                self.bigwig = [(f'{dataset}.{name}.{b1}.positive.bw',
                                f'{dataset}.{name}.{b1}.negative.bw')]
                self.bed = [f'{dataset}.{name}.{b1}.peak.cluster.bed']
            else:
                self.r1_umi = [f'{dataset}.{name}.{b1}.R1.fastq.gz',
                               f'{dataset}.{name}.{b2}.R1.fastq.gz']
                self.r2_umi = [f'{dataset}.{name}.{b1}.R2.fastq.gz',
                               f'{dataset}.{name}.{b2}.R2.fastq.gz']
                self.r1_trim = [f'{dataset}.{name}.{b1}.trim.R1.fastq.gz',
                                f'{dataset}.{name}.{b2}.trim.R1.fastq.gz']
                self.r2_trim = [f'{dataset}.{name}.{b1}.trim.R2.fastq.gz',
                                f'{dataset}.{name}.{b2}.trim.R2.fastq.gz']
                self.r1_sort = [f'{dataset}.{name}.{b1}.trim.sort.R1.fastq.gz',
                                f'{dataset}.{name}.{b2}.trim.sort.R1.fastq.gz']
                self.r2_sort = [f'{dataset}.{name}.{b1}.trim.sort.R2.fastq.gz',
                                f'{dataset}.{name}.{b2}.trim.sort.R2.fastq.gz']
                self.repeat_map_prefix = [f'{dataset}.{name}.{b1}.repeat.map.',
                                          f'{dataset}.{name}.{b2}.repeat.map.']
                self.mate1 = [f'{dataset}.{name}.{b1}.repeat.map.Unmapped.out.mate1',
                              f'{dataset}.{name}.{b2}.repeat.map.Unmapped.out.mate1']
                self.mate2 = [f'{dataset}.{name}.{b1}.repeat.map.Unmapped.out.mate2',
                              f'{dataset}.{name}.{b2}.repeat.map.Unmapped.out.mate2']
                self.mate1_sort = [f'{dataset}.{name}.{b1}.repeat.map.Unmapped.sort.mate1',
                                   f'{dataset}.{name}.{b2}.repeat.map.Unmapped.sort.mate1']
                self.mate2_sort = [f'{dataset}.{name}.{b1}.repeat.map.Unmapped.sort.mate2',
                                   f'{dataset}.{name}.{b2}.repeat.map.Unmapped.sort.mate2']
                self.genome_map_prefix = [f'{dataset}.{name}.{b1}.genome.map.',
                                          f'{dataset}.{name}.{b2}.genome.map.']
                self.bam = [f'{dataset}.{name}.{b1}.genome.map.Aligned.out.bam',
                            f'{dataset}.{name}.{b1}.genome.map.Aligned.out.bam']
                self.bam_dedup = [f'{dataset}.{name}.{b1}.dedup.bam',
                                  f'{dataset}.{name}.{b2}.dedup.bam']
                self.bigwig = [(f'{dataset}.{name}.{b1}.positive.bw',
                                f'{dataset}.{name}.{b1}.negative.bw'),
                               (f'{dataset}.{name}.{b2}.positive.bw',
                                f'{dataset}.{name}.{b2}.negative.bw')]
                self.bed = [f'{dataset}.{name}.{b1}.peak.cluster.bed',
                            f'{dataset}.{name}.{b2}.peak.cluster.bed']
        else:
            self.type = 'single'
            assert adapters, f'Single-End sample {name} does not have adapters file provided.'
            self.adapters = adapters
            assert blacklist_file, f'Single-End sample {name} does not have blacklist file provided.'
            self.blacklist_file = blacklist_file
            self.r1_umi = [f'{dataset}_{name}_umi_R1.fastq.gz']
            self.r2_umi = [f'{dataset}_{name}_umi_R2.fastq.gz']
            self.r1_trim = [f'{dataset}_{name}_umi_trim_R1.fastq.gz']
            self.r2_trim = [f'{dataset}_{name}_umi_trim_R2.fastq.gz']
            self.r1_sort = [f'{dataset}_{name}_umi_trim_sort_R1.fastq.gz']
            self.r2_sort = [f'{dataset}_{name}_umi_trim_sort_R2.fastq.gz']
            self.repeat_map_prefix = [f'{dataset}.{name}.repeat.map.']
            self.mate1 = [f'{dataset}.{name}.repeat.map.Unmapped.out.mate1']
            self.mate2 = [f'{dataset}.{name}.repeat.map.Unmapped.out.mate2']
            self.genome_map_prefix = [f'{dataset}.{name}.genome.map.']
            self.bam = [f'{dataset}.{name}.genome.map.Aligned.out.bam']
            self.repeat_map_prefix = [f'{dataset}.{name}.repeat.map.']
            self.mate1 = [f'{dataset}.{name}.repeat.map.Unmapped.out.mate1']
            self.mate2 = [f'{dataset}.{name}.repeat.map.Unmapped.out.mate2']
            self.mate1_sort = [f'{dataset}.{name}.repeat.map.Unmapped.sort.mate1']
            self.mate2_sort = [f'{dataset}.{name}.repeat.map.Unmapped.sort.mate2']
            self.genome_map_prefix = [f'{dataset}.{name}.genome.map.']
            self.bam = [f'{dataset}.{name}.genome.map.Aligned.out.bam']
            self.bam_dedup = [f'{dataset}.{name}.dedup.bam']
            self.bigwig = [(f'{dataset}.{name}.positive.bw',
                            f'{dataset}.{name}.negative.bw')]
            self.bed = [f'{dataset}.{name}.peak.cluster.bed']


def parse_manifest(manifest):
    if os.path.isfile(manifest):
        yaml = YAML(typ='safe')
        with open(manifest) as f:
            shebang = f.readline()
            try:
                data = yaml.load(f)
            except Exception:
                data = json.load(f)
    else:
        raise ValueError(f'Manifest {manifest} may not be a file or does not exist.')
    
    def _validate_entry(data, key):
        value = data.get(key, '')
        if isinstance(value, dict):
            _class = value.get('class', '')
            _path = value.get('path', '')
            if _class == 'Directory':
                if not os.path.isdir(_path):
                    raise ValueError(f'Path {_path} for {key} may not be a {_class.lower()} or does not exist.')
            elif _class == 'File':
                if not os.path.isfile(_path):
                    raise ValueError(f'Path {_path} for {key} may not be a {_class.lower()} or does not exist.')
        return value
    
    def _validate_read(read, path):
        if path:
            if not os.path.isfile(path):
                raise ValueError(f'Path {path} for {read} may not be a file or does not exist.')
        else:
            if read == 'read1':
                raise ValueError(f'Missing keyword path for {read}.')
        return path
    
    def _validate_barcodes(read):
        try:
            barcodes = read['barcodeids']
        except KeyError:
            barcodes = []
        assert isinstance(barcodes, list), f'Invalid barcode ids for read {read}.'
        return barcodes
    
    def _validate_sample(sample, barcode_fasta, dataset, blacklist_file, randomer_length):
        ip_read = [s for s in sample if 'ip_read' in s][0]
        input_read = [s for s in sample if 'input_read' in s][0]
        ip_read_name, input_read_name = ip_read.get('name', ''), input_read.get('name', '')
        if not ip_read_name:
            raise ValueError(f'The ip_read {ip_read} does not have a name.')
        if not input_read_name:
            raise ValueError(f'The input_read {input_read} does not have a name.')
        ip_read_r1 = _validate_read('read1', ip_read.get('read1', {}).get('path', ''))
        ip_read_r2 = _validate_read('read2', ip_read.get('read1', {}).get('path', ''))
        input_read_r1 = _validate_read('read1', input_read.get('read1', {}).get('path', ''))
        input_read_r2 = _validate_read('read2', input_read.get('read1', {}).get('path', ''))
        ip_read_adapters = _validate_read('adapters', ip_read.get('read1', {}).get('path', ''))
        ip_read_barcodes = _validate_barcodes(ip_read)
        input_read_adapters = _validate_read('adapters', input_read.get('adapters', {}).get('path', ''))
        input_read_barcodes = _validate_barcodes(input_read)
        ip_read = Read(ip_read_name, ip_read_r1, ip_read_r2,
                       ip_read_adapters, ip_read_barcodes, barcode_fasta, dataset)
        input_read = Read(input_read_name, input_read_r1, input_read_r2,
                          input_read_adapters, input_read_barcodes, barcode_fasta, dataset)
        return Sample(ip_read, input_read)
    
    dataset = _validate_entry(data, 'dataset')
    species = _validate_entry(data, 'species')
    genome = _validate_entry(data, 'speciesGenomeDir')
    repeat = _validate_entry(data, 'repeatElementGenomeDir')
    samples = _validate_entry(data, 'samples')
    # SE
    blacklist_file = _validate_entry(data, 'blacklist_file')
    # PE
    barcode_fasta = _validate_entry(data, 'barcodesfasta')
    randomer_length = _validate_entry(data, 'randomer_length')
    
    samples = [_validate_sample(sample, barcode_fasta, dataset, blacklist_file, randomer_length)
               for sample in samples]
    return dataset, species, genome, repeat, samples


parser = argparse.ArgumentParser(description=__doc__, prog='clip')
parser.add_argument('MANIFEST', type=str, help='YAML or JSON manifest file describing paths of your dataset.')
parser.add_argument('-o', type=str, dest='OUTDIR',
                    help="Path of the base output directory, default: the manifest file's parent directory.")
parser.add_argument('-j', type=str, dest='JOBNAME',
                    help="Name of your job, default: eCLIP", default='eCLIP')
parser.add_argument('-e', type=str, dest='EMAIL',
                    help='Email address for notify you the start and end of you job.')
parser.add_argument('-s', type=str, dest='SCHEDULER',
                    help='Name of the scheduler on your cluster, e.g., PBS (or qsub) or SBATCH (or slurm).')
parser.add_argument('-t', type=int, dest='TIME',
                    help='Time (in integer hours) needed for your job, default: 32.', default=32)
parser.add_argument('-m', type=int, dest='MEMORY',
                    help='Amount of memory (in GB) for all cores needed for your job, default: 32.', default=32)
parser.add_argument('-n', type=int, dest='CORES',
                    help='Number of cores needed for your job, default: 8.', default=8)
parser.add_argument('--dry-run', action='store_true',
                    help='Run the analysis in debug mode and keep cache files.')

args = parser.parse_args()
manifest = os.path.abspath(args.MANIFEST)
if os.path.isfile(manifest):
    DATASET, SPECIES, GENOME, REPEAT, SAMPLES = parse_manifest(manifest)
    
    outdir = os.path.abspath(args.OUTDIR) if args.OUTDIR else os.path.dirname(manifest)
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
        
    logger = logging.getLogger("CLIP")
    logger.setLevel(logging.DEBUG)
    console = logging.StreamHandler(sys.stderr)
    console.setFormatter(logging.Formatter('%(asctime)s [%(levelname)s] %(message)s',
                                           datefmt='%Y-%m-%d %H:%M:%S'))
    console.setLevel(logging.DEBUG)
    logger.addHandler(console)
    log = f'clip_{datetime.today().strftime("%Y-%m-%d_%H:%M:%S")}.log'
    handler = logging.FileHandler(filename=log, mode='w')
    handler.setFormatter(logging.Formatter('%(asctime)s [%(levelname)s] %(message)s',
                                           datefmt='%Y-%m-%d %H:%M:%S'))
    handler.setLevel(logging.DEBUG)
    logger.addHandler(handler)
else:
    raise ValueError('The manifest "{}" may not be a file or does not exist.'.format(manifest))
    

def extract_umi(read):
    pass


def demux_barcode(read):
    pass


def prepare_fastqs():
    for sample in SAMPLES:
        fastqs = [sample.ip_read.r1, sample.ip_read.r2,
                  sample.input_read.r1, sample.ip_read.r2]
        links = [sample.ip_read.r1_link, sample.ip_read.r2_link,
                 sample.input_read.r1_link, sample.ip_read.r2_link]
        reads = [sample.ip_read, sample.ip_read, sample.input_read, sample.input_read]
        for fastq, link, read in zip(fastqs, links, reads):
            yield fastq, link, link
            
            
def touch(file):
    pass
    

@ruffus.jobs_limit(1)
@ruffus.files(prepare_fastqs)
def soft_link_fastq(fastq, link, read):
    if fastq == link:
        logger.warning("No symbolic link made. You are directly working on the original data files.")
    else:
        if fastq:
            logger.debug(f'Soft link {os.path.basename(fastq)}:\n  ln -s {fastq} {link}')
            os.symlink(fastq, link)
        else:
            touch(link)


def prepare_umi_demux():
    for sample in SAMPLES:
        for read in (sample.ip_read, sample.input_read):
            for r1_umi, r2_umi in zip(read.r1_umi, read.r2_umi):
                yield (read.r1_link, read.r2_link), (r1_umi, r2_umi), read
            
            
@ruffus.files(prepare_umi_demux)
@ruffus.follows(soft_link_fastq)
def umi_demux(inputs, outputs, read):
    if read.type == 'single':
        extract_umi(read)
    else:
        demux_barcode(read)
        
        
def trim_adapters(inputs, outputs, read):
    pass


def fastq_sort(fastq, sorted_fastq):
    pass


def sort_fast(inputs, outputs, read):
    fastq_sort(inputs, outputs)


def map__to_repeat_element(inputs, outputs, read):
    pass


def map_to_reference_genome(inputs, outputs, read):
    pass


def dedup_bam():
    pass


def make_bigwig():
    pass


def call_peaks():
    pass


SBATCH = """#!/usr/bin/env bash

#SBATCH -n {cores}                  # Number of cores (-n)
#SBATCH -N 1                        # Ensure that all cores are on one Node (-N)
#SBATCH -t {runtime}                # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH --mem={memory}G             # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --job-name={jobname}        # Short name for the job
"""
SBATCH_EMAIL = """
#SBATCH --mail-user={email}
#SBATCH --mail-type=ALL
"""

PBS = """ #!/usr/bin/env bash

#PBS -l nodes=1:ppn={cores}
#PBS -l walltime={runtime}
#PBS -l vmem={memory}gb
#PBS -j oe
#PBS -N {jobname}
"""
PBS_EMAIL = """
#PBS -M {email}
#PBS -m abe
"""

CODE = r"""
LOCAL=local
export TMPDIR={tmpdir}
export TEMP={tmpdir}
export TMP={tmpdir}
source CLIP_ENVIRONMENT

echo [$(date +"%m-%d-%Y %H:%M:%S")] "eCLIP start."
source ECLIP_ENVIRONMENT
cwltool {debug} \
    --no-container \
    --tmpdir-prefix={tmpdir}/ \
    --outdir={outdir} \
    ECLIP/cwl/{eclip_script} \
    {manifest} \
    2>&1 \
    | tee {log_dir}/clip.${cluster_job_id}.${cluster_job_name}.log {stdout}
echo [$(date +"%m-%d-%Y %H:%M:%S")] "eCLIP complete."

echo [$(date +"%m-%d-%Y %H:%M:%S")] "Merge Peak (IDR) start."
source MERGE_PEAK_ENVIRONMENT
cwltool {debug} \
    --no-container \
    --tmpdir-prefix={tmpdir}/ \
    --outdir={idr_outdir} \
    MERGE_PEAK/cwl/{idr_script} \
    {idr_config} \
    2>&1 \
    | tee {log_dir}/idr.${cluster_job_id}.${cluster_job_name}.log {stdout}
echo [$(date +"%m-%d-%Y %H:%M:%S")] "Merge Peak (IDR) complete."

echo [$(date +"%m-%d-%Y %H:%M:%S")] "Cleaning up ..."
rm -rf {tmpdir}
echo [$(date +"%m-%d-%Y %H:%M:%S")] "All done."
"""

def validate_clip_manifest(manifest):
    shebang, data = parse_manifest(manifest)
    if 'singleend' in shebang:
        return validate_se_manifest(data)
    elif 'pairedend' in shebang:
        return validate_pe_manifest(data)
    else:
        raise ValueError('Manifest file does not have a valid shebang line.')


def write_ird_manifest(species, results, eclip_outdir, script_outdir):
    input_reads = []
    lines = ['', '', f'species: {species}', 'samples:', '  -']
    for result in results:
        lines.append(f'    - name: {result[0]}')
        lines.append(f'      ip_bam:')
        lines.append(f'        class: File')
        lines.append(f'        path: {os.path.join(eclip_outdir, result[1])}')
        lines.append(f'      input_bam:')
        lines.append(f'        class: File')
        lines.append(f'        path: {os.path.join(eclip_outdir, result[2])}')
        lines.append(f'      peak_clusters:')
        lines.append(f'        class: File')
        lines.append(f'        path: {os.path.join(eclip_outdir, result[3])}')
        input_reads.append(result[4:])
    input_reads = [''.join(read) for read in input_reads]
    mode = '2inputs' if len(input_reads) == len(set(input_reads)) > 1 else '1input'
    script = f'wf_full_IDR_pipeline_{mode}_scatter.cwl'
    shebang = f'#!/usr/bin/env eCLIP_full_IDR_pipeline_{mode}_scatter_singleNode'
    lines[0] = shebang
    idr_config = os.path.join(script_outdir, 'idr.yaml')
    with open(idr_config, 'w') as o:
        o.writelines(f'{line}\n' for line in lines)
    return script, idr_config


def _mkdir(base, directory):
    folder = os.path.join(base, directory)
    if not os.path.isdir(folder):
        try:
            os.mkdir(folder)
        except OSError as e:
            print(e)
            sys.exit(1)
    return folder


def schedule(scheduler, manifest, outdir, tmpdir, script_dir, log_dir, idr_outdir, eclip_script, idr_script, idr_config,
             debug, cores, memory, hours, jobname, email, cluster_job_name, cluster_job_id, stdout):
    
    if scheduler.upper() in ('PBS', 'QSUB'):
        runtime = f'{hours}:00:00'
        directive, exe, mail = PBS, 'qsub', PBS_EMAIL
    elif scheduler.upper() in ('SLURM', 'SBATCH'):
        days, hours = divmod(hours, 24)
        runtime = f'{days}-{hours:02}:00'
        directive, exe, mail = SBATCH, 'sbatch', SBATCH_EMAIL
    else:
        raise ValueError(f'Unsupported scheduler: {scheduler}, see help for supported schedulers.')
    
    data = {'cores': cores, 'runtime': runtime, 'memory': memory, 'jobname': jobname, 'manifest': manifest,
            'outdir': outdir, 'tmpdir': tmpdir, 'script_dir': script_dir, 'log_dir': log_dir, 'idr_outdir': idr_outdir,
            'eclip_script': eclip_script, 'idr_script': idr_script, 'idr_config': idr_config, 'debug': debug,
            'cluster_job_name': cluster_job_name, 'cluster_job_id': cluster_job_id, 'stdout': stdout}
    if email:
        data['email'] = email
        text = [directive, mail, CODE]
    else:
        text = [directive, CODE]
    text = ''.join(text).format(**data)
    
    submitter = os.path.join(script_dir, 'submit.sh')
    with open(submitter, 'w') as o:
        o.write(text)
    print(f'Job submit script was saved to: {submitter}')
    subprocess.run([exe, submitter], cwd=log_dir)
    print(f'Job {jobname} was successfully submitted with the following settings:')
    data = {'Job name:': jobname, 'Manifest file:': manifest, 'Output directory:': outdir,
            'Number of cores:': cores, 'Job memory:': memory, 'Job runtime:': f'{runtime} (D-HH:MM)'}
    for k, v in data.items():
        print(f'{k:>20} {v}')


def clip(manifest, outdir, tmpdir, script_dir, log_dir, idr_outdir, eclip_script, idr_script, idr_config,
         debug, cluster_job_name, cluster_job_id, stdout):
    data = {'manifest': manifest, 'outdir': outdir, 'tmpdir': tmpdir, 'script_dir': script_dir,
            'log_dir': log_dir, 'idr_outdir': idr_outdir, 'eclip_script': eclip_script, 'idr_script': idr_script,
            'idr_config': idr_config, 'debug': debug, 'stdout': stdout,
            'cluster_job_name': cluster_job_name, 'cluster_job_id': cluster_job_id}
    text = ['#!/usr/bin/env bash', CODE]
    text = '\n'.join(text).format(**data)
    script = os.path.join(script_dir, 'script.sh')
    with open(script, 'w') as o:
        o.write(text)
    subprocess.run(['bash', script], cwd=log_dir)


# def main():
#     parser = argparse.ArgumentParser(description=__doc__, prog='clip')
#     parser.add_argument('MANIFEST', type=str, help='Manifest YAML or JSON file describing paths of your dataset.')
#     parser.add_argument('-o', type=str, dest='OUTDIR',
#                         help="Path of the base output directory, default: the manifest file's parent directory.")
#     parser.add_argument('-j', type=str, dest='JOBNAME',
#                         help="Name of your job, default: eCLIP", default='eCLIP')
#     parser.add_argument('-e', type=str, dest='EMAIL',
#                         help='Email address for notify you the start and end of you job.')
#     parser.add_argument('-s', type=str, dest='SCHEDULER',
#                         help='Name of the scheduler on your cluster, e.g., PBS (or qsub) or SBATCH (or slurm).')
#     parser.add_argument('-t', type=int, dest='TIME',
#                         help='Time (in integer hours) needed for your job, default: 32.', default=32)
#     parser.add_argument('-m', type=int, dest='MEMORY',
#                         help='Amount of memory (in GB) for all cores needed for your job, default: 32.', default=32)
#     parser.add_argument('-n', type=int, dest='CORES',
#                         help='Number of cores needed for your job, default: 8.', default=8)
#     parser.add_argument('--debug', action='store_true', dest='DEBUG',
#                         help='Run the analysis in debug mode and keep cache files.')
#
#     args = parser.parse_args()
#     manifest = os.path.abspath(args.MANIFEST)
#     if os.path.isfile(manifest):
#         results = validate_clip_manifest(manifest)
#     else:
#         raise ValueError('The manifest "{}" may not be a file or does not exist.'.format(manifest))
#
#     shebang, data = parse_manifest(manifest)
#     outdir = os.path.abspath(args.OUTDIR) if args.OUTDIR else os.path.dirname(manifest)
#     if not os.path.isdir(outdir):
#         os.mkdir(outdir)
#
#     eclip_config = os.path.join(outdir, 'clip.yaml')
#     if not os.path.isfile(eclip_config):
#         shutil.copy(manifest, eclip_config)
#     basedir = os.path.abspath(args.OUTDIR) if args.OUTDIR else os.path.dirname(manifest)
#     basedir = _mkdir(os.path.dirname(basedir), os.path.basename(basedir))
#
#     dirs = ['tmp', 'results', 'logs', 'scripts']
#     tmpdir, outdir, log_dir, script_dir = [_mkdir(basedir, d) for d in dirs]
#
#     eclip_config = shutil.copy(manifest, os.path.join(script_dir, 'clip.yaml'))
#     eclip_script = 'wf_get_peaks_scatter_pe.cwl' if 'paired' in shebang else 'wf_get_peaks_scatter_se.cwl'
#     idr_outdir = _mkdir(basedir, 'idr')
#     idr_script, idr_config = write_ird_manifest(data['species'], results, outdir, script_dir)
#     debug = f'--debug --cachedir={_mkdir(basedir, "cache")}' if args.DEBUG else ''
#     stdout = '' if args.DEBUG else '> /dev/null'
#     if args.SCHEDULER:
#         names = {'PBS': '{PBS_JOBNAME}', 'SLURM': '{SLURM_JOB_NAME}', 'SBATCH': '{SLURM_JOB_NAME}'}
#         cluster_job_name = names.get(args.SCHEDULER.upper(), '')
#         if not cluster_job_name:
#             raise ValueError(f'Unsupported scheduler: {args.SCHEDULER}, accepts pbs, qsub, slurm, sbatch.')
#         ids = {'PBS': '{PBS_JOBID}', 'QSUB': '{PBS_JOBID}', 'SLURM': '{SLURM_JOB_ID}', 'SBATCH': '{SLURM_JOB_ID}'}
#         cluster_job_id = ids.get(args.SCHEDULER.upper(), '')
#         schedule(args.SCHEDULER, eclip_config, outdir, tmpdir, script_dir, log_dir, idr_outdir,
#                  eclip_script, idr_script, idr_config, debug, args.CORES, args.MEMORY, args.TIME,
#                  args.JOBNAME, args.EMAIL, cluster_job_name, cluster_job_id, stdout)
#     else:
#         clip(eclip_config, outdir, tmpdir, script_dir, log_dir, idr_outdir, eclip_script,
#              idr_script, idr_config, debug, 'LOCAL', 'USER', stdout)

if args.dry_run:
    ruffus.pipeline_printout()
else:
    if args.scheduler:
        pass
    else:
        ruffus.pipeline_run(multiprocess=args.cores, verbose=1)


if __name__ == '__main__':
    pass

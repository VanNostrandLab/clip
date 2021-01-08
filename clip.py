#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
A shortcut for using eCLIP pipeline to process RNA-Seq data on SLURM (SBATCH) Cluster.
"""

import os
import sys
import shutil
import subprocess
import argparse

from ruamel.yaml import YAML

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
#PBS -l walltime={runtime}:00:00
#PBS -l vmem={memory}GB
#PBS -J {jobname}
"""
PBS_EMAIL = """
#PBS -M {email}
#PBS -m n
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
    > {log_dir}/clip.${cluster_job_id}.${cluster_job_id}.log 2>&1
echo [$(date +"%m-%d-%Y %H:%M:%S")] "eCLIP complete."

echo [$(date +"%m-%d-%Y %H:%M:%S")] "Merge Peak (IDR) start."
source MERGE_PEAK_ENVIRONMENT
cwltool {debug} \
    --no-container \
    --tmpdir-prefix={tmpdir}/ \
    --outdir={idr_outdir} \
    MERGE_PEAK/cwl/{idr_script} \
    {idr_config} \
    > {log_dir}/idr.${cluster_job_id}.${cluster_job_id}.log 2>&1
echo [$(date +"%m-%d-%Y %H:%M:%S")] "Merge Peak (IDR) complete."
"""

DESCRIPTION = """A wrapper for using eCLIP to process CLIP data."""


def parse_manifest(manifest):
    if os.path.isfile(manifest):
        yaml = YAML(typ='safe')
        with open(manifest) as f:
            shebang = f.readline()
            data = yaml.load(f)
        return shebang, data
    else:
        raise ValueError(f'Manifest {manifest} may not be a file or does not exist.')


def _validate_entry(data, key):
    try:
        value = data[key]
        if isinstance(value, dict):
            _class = value.get('class', '')
            _path = value.get('path', '')
            if _class == 'Directory':
                if not os.path.isdir(_path):
                    raise ValueError(f'Path {_path} for {key} may not be a {_class.lower()} or does not exist.')
            elif _class == 'File':
                if not os.path.isfile(_path):
                    raise ValueError(f'Path {_path} for {key} may not be a {_class.lower()} or does not exist.')
        else:
            if not value:
                raise ValueError(f'No valid value was found for {key}')
    except ValueError:
        raise ValueError(f'Manifest missing keyword {key}.')


def _validate_read(read, path):
    if path:
        if not os.path.isfile(path):
            raise ValueError(f'Path {path} for {read} may not be a file or does not exist.')
    else:
        if read == 'read1':
            raise ValueError(f'Missing keyword path for {read}.')
    return path


def _validate_se_sample(sample, dataset):
    ip, inp = [s for s in sample if 'ip_read' in s][0], [s for s in sample if 'input_read' in s][0]
    ip_name, inp_name = ip.get('name', ''), inp.get('name', '')
    if not ip_name:
        raise ValueError(f'The ip_read {ip} does not have a name.')
    if not inp_name:
        raise ValueError(f'The input_read {inp} does not have a name.')
    _validate_read('read1', ip.get('read1', {}).get('path', ''))
    _validate_read('adapters', ip.get('adapters', {}).get('path', ''))
    input_read = _validate_read('read1', inp.get('read1', {}).get('path', ''))
    _validate_read('adapters', inp.get('adapters', {}).get('path', ''))
    ip_bam = f'{dataset}.{ip_name}.umi.r1.fq.genome-mappedSoSo.rmDupSo.bam'
    inp_bam = f'{dataset}.{inp_name}.umi.r1.fq.genome-mappedSoSo.rmDupSo.bam'
    bed = f'{dataset}.{ip_name}.umi.r1.fq.genome-mappedSoSo.rmDupSo.peakClusters.bed'
    return ip_name, ip_bam, inp_bam, bed, input_read, ''


def _validate_pe_sample(sample, dataset):
    ip, inp = [s for s in sample if 'ip_read' in s][0], [s for s in sample if 'input_read' in s][0]
    ip_name, inp_name = ip.get('name', ''), inp.get('name', '')
    if not ip_name:
        raise ValueError(f'The ip_read {ip} does not have a name.')
    if not inp_name:
        raise ValueError(f'The input_read {inp} does not have a name.')
    _validate_read('read1', ip.get('read1', {}).get('path', ''))
    _validate_read('read2', ip.get('read2', {}).get('path', ''))
    input_read1 = _validate_read('read1', inp.get('read1', {}).get('path', ''))
    input_read2 = _validate_read('read2', inp.get('read2', {}).get('path', ''))
    ip_barcodes, inp_barcodes = ip.get('barcodeids', []), inp.get('barcodeids', [])
    if not ip_barcodes:
        raise ValueError(f'The ip_read {ip} does not have barcodes.')
    if not inp_barcodes:
        raise ValueError(f'The input_read {inp} does not have barcodes.')
    ip_bam = f'{dataset}.{ip_name}.{ip_barcodes[0]}.r1.fq.genome-mappedSo.rmDupSo.bam'
    inp_bam = f'{dataset}.{inp_name}.{inp_barcodes[0]}.r1.fq.genome-mappedSo.rmDupSo.bam'
    bed = f'{dataset}.{ip_name}.{ip_barcodes[0]}.r1.fq.genome-mappedSo.rmDupSo.peakClusters.bed'
    return ip_name, ip_bam, inp_bam, bed, input_read1, input_read2


def validate_se_manifest(data):
    for key in ('dataset', 'species', 'speciesGenomeDir', 'repeatElementGenomeDir', 'chrom_sizes',
                'samples', 'blacklist_file'):
        _validate_entry(data, key)
    results = [_validate_se_sample(sample, data['dataset']) for sample in data['samples']]
    return results


def validate_pe_manifest(data):
    for key in ('dataset', 'species', 'speciesGenomeDir', 'repeatElementGenomeDir', 'chrom_sizes',
                'samples', 'randomer_length', 'barcodesfasta'):
        _validate_entry(data, key)
    results = [_validate_pe_sample(sample, data['dataset']) for sample in data['samples']]
    return results


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
             debug, cores, memory, hours, jobname, email, cluster_job_name, cluster_job_id):
    if scheduler.upper() == 'PBS':
        runtime = f'{hours}:00:00'
        mail = PBS_EMAIL
    elif scheduler.upper() in ('SLURM', 'SBATCH'):
        days, hours = divmod(hours, 24)
        runtime = f'{days}-{hours:02}:00'
        mail = SBATCH_EMAIL
    else:
        raise ValueError(f'Unsupported scheduler: {scheduler}, see help for supported schedulers.')
    
    data = {'cores': cores, 'runtime': runtime, 'memory': memory, 'jobname': jobname, 'manifest': manifest,
            'outdir': outdir, 'tmpdir': tmpdir, 'script_dir': script_dir, 'log_dir': log_dir, 'idr_outdir': idr_outdir,
            'eclip_script': eclip_script, 'idr_script': idr_script, 'idr_config': idr_config, 'debug': debug,
            'cluster_job_name': cluster_job_name, 'cluster_job_id': cluster_job_id}
    if email:
        data['email'] = email
        text = [SBATCH, mail, CODE]
    else:
        text = [SBATCH, CODE]
    text = ''.join(text).format(**data)
    
    submitter = os.path.join(script_dir, 'submit.sh')
    with open(submitter, 'w') as o:
        o.write(text)
    print(f'Job submit script was saved to: {submitter}')
    subprocess.run(['sbatch', submitter], cwd=log_dir)
    print(f'Job {jobname} was successfully submitted with the following settings:')
    data = {'Job name:': jobname, 'Manifest file:': manifest, 'Output directory:': outdir,
            'Number of cores:': cores, 'Job memory:': memory, 'Job runtime:': f'{runtime} (D-HH:MM)'}
    for k, v in data.items():
        print(f'{k:>20} {v}')


def clip(manifest, outdir, tmpdir, script_dir, log_dir, idr_outdir, eclip_script, idr_script, idr_config,
         debug, cluster_job_name, cluster_job_id):
    data = {'manifest': manifest, 'outdir': outdir, 'tmpdir': tmpdir, 'script_dir': script_dir,
            'log_dir': log_dir, 'idr_outdir': idr_outdir, 'eclip_script': eclip_script, 'idr_script': idr_script,
            'idr_config': idr_config, 'debug': debug,
            'cluster_job_name': cluster_job_name, 'cluster_job_id': cluster_job_id}
    text = ['#!/usr/bin/env bash', CODE]
    text = '\n'.join(text).format(**data)
    script = os.path.join(script_dir, 'script.sh')
    with open(script, 'w') as o:
        o.write(text)
    subprocess.run(['bash', script], cwd=log_dir)


def main():
    parser = argparse.ArgumentParser(description=__doc__, prog='clip')
    parser.add_argument('MANIFEST', type=str, help='Manifest YAML or JSON file describing paths of your dataset.')
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
    parser.add_argument('--debug', action='store_true', dest='DEBUG',
                        help='Run the analysis in debug mode and keep cache files.')
    
    args = parser.parse_args()
    manifest = os.path.abspath(args.MANIFEST)
    if os.path.isfile(manifest):
        results = validate_clip_manifest(manifest)
    else:
        raise ValueError('The manifest "{}" may not be a file or does not exist.'.format(manifest))
    
    shebang, data = parse_manifest(manifest)
    basedir = os.path.abspath(args.OUTDIR) if args.OUTDIR else os.path.dirname(manifest)
    basedir = _mkdir(os.path.dirname(basedir), os.path.basename(basedir))
    
    dirs = ['temporary', 'results', 'logs', 'scripts']
    tmpdir, outdir, log_dir, script_dir = [_mkdir(basedir, d) for d in dirs]
    
    eclip_config = shutil.copy(manifest, os.path.join(script_dir, 'clip.yaml'))
    eclip_script = 'wf_get_peaks_scatter_pe.cwl' if 'paired' in shebang else 'wf_get_peaks_scatter_se.cwl'
    idr_outdir = _mkdir(basedir, 'idr')
    idr_script, idr_config = write_ird_manifest(data['species'], results, outdir, script_dir)
    debug = f'--debug --cachedir={_mkdir(basedir, "cache")}' if args.DEBUG else ''
    if args.SCHEDULER:
        names = {'PBS': '{PBS_JOBNAME}', 'SLURM': '{SLURM_JOB_NAME}', 'SBATCH': '{SLURM_JOB_NAME}'}
        cluster_job_name = names.get(args.SCHEDULER.upper(), '')
        if not cluster_job_name:
            raise ValueError(f'Unsupported scheduler: {args.SCHEDULER}, accepts pbs, qsub, slurm, sbatch.')
        ids = {'PBS': '{PBS_JOBID}', 'QSUB': '{PBS_JOBID}', 'SLURM': '{SLURM_JOB_ID}', 'SBATCH': '{SLURM_JOB_ID}'}
        cluster_job_id = ids.get(args.SCHEDULER.upper(), '')
        schedule(args.SCHEDULER, eclip_config, outdir, tmpdir, script_dir, log_dir, idr_outdir,
                 eclip_script, idr_script, idr_config, debug, args.CORES, args.MEMORY, args.TIME,
                 args.JOBNAME, args.EMAIL, cluster_job_name, cluster_job_id)
    else:
        clip(eclip_config, outdir, tmpdir, script_dir, log_dir, idr_outdir, eclip_script,
             idr_script, idr_config, debug, 'LOCAL', 'USER')


if __name__ == '__main__':
    main()

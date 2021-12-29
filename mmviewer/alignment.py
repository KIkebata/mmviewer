#!/usr/bin/env python3
import subprocess, os, shutil, glob
import pandas as pd
import sys
from .utility import make_directory

def bwa_mem_py(ref_fa, sample_names, out_dir, in_fq1_list, in_fq2_list=None):
    """Run bwa mem to generate sorted bam.

    :param ref_fa: Path to input reference file in fasta format.
    :type ref_fa: str
    :param sample_names: Sample name list (same length with in_fq1_list).
    :type sample_names: list
    :param out_dir: Path to output directory.
    :type out_dir: str
    :param in_fq1_list: List of trimmed reads file path (fq format or fq.gz format)
    :type in_fq1_list: list
    :param in_fq2_list: List of trimmed reads file path for paired end reads, defaults to None (optional)
    :type in_fq2_list: list, None
    :return: Output sorted bam file list.
    :rtype: list
    """
    if in_fq2_list is None:
        fq_list = in_fq1_list
    else:
        fq_list = [fq1 + ' ' + fq2 for fq1, fq2 in zip(in_fq1_list, in_fq2_list)]
    output_bam_dir = os.path.join(out_dir, 'bam')
    os.mkdir(output_bam_dir)
    bam_files = [os.path.join(output_bam_dir, ll + '.bam') for ll in sample_names]
    subprocess.run(' '.join([
        'bwa', 'index', ref_fa
    ]), shell=True)
    for i in range(len(sample_names)):
        subprocess.run(' '.join([
            'bwa', 'mem', '-t', str(os.cpu_count()), ref_fa, fq_list[i], '|',
            'samtools', 'view', '-@', str(os.cpu_count()), '-bS', '|',
            'samtools', 'fixmate', '-m', '-@', str(os.cpu_count()), '-', '-', '|',
            'samtools',  'sort', '-@', str(os.cpu_count()), '|',
            'samtools', 'markdup', '-@', str(os.cpu_count()), '-r', '-', bam_files[i]
        ]), shell=True)
        subprocess.run(' '.join([
            'samtools', 'index', bam_files[i]
        ]), shell=True)
    return bam_files


def make_alignment(al_config_csv, ref_fa, out_dir):
    """Make alignment by using bwa mem
    Before runnning this function, prepare alignment_config.csv file.
    Alignment_config.csv file consists of 2 or 3 columns with a header in first row.
    1st column is sample name.
    2nd column is path to trimmed read in fastq format or its compressed format (.gz)
    3rd column is path to trimmed read in fastq format or its compressed format (.gz) of reverse reads if paired end read (option)

    :param al_config_csv: path to alignment_config.csv
    :type al_config_csv: str
    :param ref_fa: path to a reference sequence in fasta format
    :type ref_fa: str
    :param out_dir: path to output directory
    :type out_dir: str
    """
    print_sep = '#' * 40
    out_dir = make_directory(dir_name='alignment', out_dir=out_dir)
    ref_fa_dir = os.path.join(out_dir, 'reference')
    os.mkdir(ref_fa_dir)
    shutil.copy2(src=ref_fa, dst=ref_fa_dir)
    ref_fa = glob.glob(os.path.join(ref_fa_dir, '*'))[0]
    al_config_df = pd.read_csv(al_config_csv)
    if len(al_config_df.columns) == 2:
        in_fq2_list = None
    elif len(al_config_df.columns) == 3:
        in_fq2_list = al_config_df.iloc[:, 2].values.tolist()
    else:
        print()
        sys.exit('Error! Invalid alignment config file.')
    sample_names = al_config_df.iloc[:, 0].values.tolist()
    in_fq1_list = al_config_df.iloc[:, 1].values.tolist()
    print(print_sep)
    print('Running bwa mem.......')
    bam_files = bwa_mem_py(
        ref_fa=ref_fa,
        sample_names=sample_names,
        out_dir=out_dir,
        in_fq1_list=in_fq1_list,
        in_fq2_list=in_fq2_list
    )
    pd.DataFrame({'sample': sample_names, 'bam':bam_files}).to_csv(os.path.join(out_dir, 'graph_config.csv', ), index=False)
    print(print_sep)
    print('Finish.')
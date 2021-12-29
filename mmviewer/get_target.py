#!/usr/bin/env python3
import subprocess, os, shutil, glob, re
import pandas as pd
from .utility import make_directory

def prodigal_py(in_ref_fa, out_dir):
    """Predict CDS of fasta file

    :param in_ref_fa: A fasta file as input
    :type in_ref_fa: str
    :param out_dir: Path to output directory
    :type out_dir: str
    :return: output gff file path, predicted cds (nucl.) fa file path, predicted cds (prot.) fa file path 
    :rtype: (str, str, str)
    """
    gff_file = os.path.join(out_dir, 'cds.gff')
    nucl_fa_file = os.path.join(out_dir, 'cds_nucl.fa')
    prot_fa_file = os.path.join(out_dir, 'cds_prot.fa')
    subprocess.run(' '.join([
        'prodigal', '-i', in_ref_fa,
            '-o', gff_file,
            '-f', 'gff',
            '-d', nucl_fa_file,
            '-a', prot_fa_file
    ]), shell=True)
    return gff_file, nucl_fa_file, prot_fa_file

def makeblastdb_py(seq_type, ref_fa, out_dir):
    """Run makeblastdb

    :param seq_type: ['prot' | 'nucl']
    :type seq_type: str
    :param ref_fa: Path to reference.fasta
    :type ref_fa: str
    :param out_dir: Path to output directory
    :type out_dir: str
    :return: File path to blast database
    :rtype: str
    """
    blastdb_dir = os.path.join(out_dir, 'blastdb')
    os.mkdir(blastdb_dir)
    shutil.copy2(src=ref_fa, dst=blastdb_dir)
    blastdb_path = glob.glob(os.path.join(blastdb_dir, '*'))[0]
    subprocess.run(' '.join([
        'makeblastdb', '-in', blastdb_path, '-dbtype', seq_type
    ]), shell=True)
    return blastdb_path

def blast_py(query_fa, seq_type, blastdb_path, out_dir, anno_fmt):
    """Run blast
    If seq_type is 'nucl', run BLASTN. If seq_type is 'prot', run BLASTX.

    :param query_fa: Path to cds.fasta
    :type query_fa: str
    :param seq_type: ['nucl' | 'prot'] 
    :type seq_type: str
    :param blastdb_path: File path to blast database 
    :type blastdb_path: str
    :param anno_fmt: result format
    :type anno_fmt: list
    :param out_dir: Path to output directory
    :type out_dir: str
    :return: Path to blast output
    :rtype: str
    """
    if seq_type == 'prot':
        blast_out = os.path.join(out_dir, 'blastx_out.tsv')
    elif seq_type == 'nucl':
        blast_out = os.path.join(out_dir, 'blastn_out.tsv')
    else:
        return None
    out_fmt = '"6 '+' '.join(anno_fmt) + '"'
    if seq_type == 'prot':
        subprocess.run(' '.join([
            'blastx', '-query', query_fa,
            '-db', blastdb_path,
            '-out', blast_out,
            '-outfmt', out_fmt,
            '-evalue', '1e-5',
            '-num_threads', str(os.cpu_count())
        ]), shell=True)
    elif seq_type == 'nucl':
        subprocess.run(' '.join([
            'blastn', '-query', query_fa,
            '-db', blastdb_path,
            '-out', blast_out,
            '-outfmt', out_fmt,
            '-evalue', '1e-5',
            '-num_threads', str(os.cpu_count())
        ]), shell=True)
    return blast_out

def get_bed_from_blast(gff_file, nucl_fa_file, blast_out, anno_fmt, out_dir, f_interval=0, r_interval=0):
    """Make bed file from gff file from prodigal, cds_nucl.fa from prodigal, and blast result

    :param gff_file: gff file path made by prodigal
    :type gff_file: str
    :param nucl_fa_file: cds_nucl.fa file path made by prodigal
    :type nucl_fa_file: str
    :param blast_out: blast output file path
    :type blast_out: str
    :param anno_fmt: blast annotation list
    :type anno_fmt: list
    :param out_dir: Path to output directory
    :type out_dir: str
    :param f_interval: upper region you want to analyze than objective cds , defaults to 0
    :type f_interval: int
    :param r_interval: lower region you want to analyze than objective cds, defaults to 0
    :type r_interval: int
    """
    out_bed = os.path.join(out_dir, 'target.bed')
    gff_fmt = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
    gff_df = pd.read_table(gff_file, names=gff_fmt, comment='#')
    tmp_list = [re.sub(r'\;$', '', re.search(r'ID\=.+?\;', ll).group()) for ll in gff_df['attribute']]
    gff_df2 = pd.concat([gff_df, pd.DataFrame({'ID':tmp_list})], axis=1)
    with open(nucl_fa_file, 'r') as f:
        lines = f.read().splitlines()
    tmp_list = [re.sub(r'^\>', '', ll) for ll in lines if ll.startswith('>')]
    seqid_list = [ll.split(' # ')[0] for ll in tmp_list]
    id_list = [re.sub(r'\;$', '', re.search(r'ID\=.+?\;', ll).group()) for ll in tmp_list]
    nucl_head_df = pd.DataFrame({'qseqid':seqid_list, 'ID':id_list})
    gff_df3 = pd.merge(left=gff_df2, right=nucl_head_df, on='ID', how='left')
    blast_out_df = pd.read_table(blast_out, names=anno_fmt)
    blast_out_df2 = pd.merge(left=blast_out_df, right=gff_df3, on='qseqid', how='left')
    bed_score = (10 * blast_out_df2['pident'] * blast_out_df2['length'] / \
        (blast_out_df2['slen'] + blast_out_df2['length'] - (blast_out_df2['send'] - blast_out_df2['sstart'] + 1))).astype(int)
    c_start, c_end = [], []
    for j in range(len(blast_out_df2)):
        if blast_out_df2['strand'][j] == '+':
            c_start += [blast_out_df2['start'][j] - 1 - f_interval]
            c_end += [blast_out_df2['end'][j] + r_interval]
        elif blast_out_df2['strand'][j] == '-':
            c_start += [blast_out_df2['start'][j] - 1 - r_interval]
            c_end += [blast_out_df2['end'][j] + f_interval]
        if c_start[j] < 0:
            c_start[j] = 0
        # consider end if possible
    bed_df = pd.DataFrame({
        'chrom':blast_out_df2['seqname'],
        'chromStart':c_start,
        'chromEnd':c_end,
        'name':blast_out_df2['stitle'],
        'score':bed_score,
        'strand':blast_out_df2['strand'],
        'CDS_Start':blast_out_df2['start'] - 1,
        'CDS_End':blast_out_df2['end']
    })
    bed_df.to_csv(out_bed, sep='\t', index=False, header=False)
    with open(out_bed) as f:
        s = f.read()
    s = '#' + '\t'.join(bed_df.columns.values.tolist()) + '\n' + s 
    with open(out_bed, 'w') as f:
        f.write(s)
    return out_bed

def add_fasta_entry(ref_fa, gff_file):
    """add fasta entry at the last of gff_file

    :param ref_fa: path to reference.fasta
    :type ref_fa: str
    :param gff_file: gff file path made by prodigal
    :type gff_file: str
    """
    with open(gff_file, 'a') as f:
        with open(ref_fa) as ref_f:
            f.write('###\n##FASTA\n')
            f.write(ref_f.read())

def get_target_position(ref_fa, out_dir, gene_seq_type, gene_seq, f_interval=0, r_interval=0):
    """Generate target.bed file to indicate target region in genome

    :param ref_fa: Path to a file of complete genome or contig in fasta format as a reference
    :type ref_fa: str
    :param out_dir: Path to output directory
    :type out_dir: str
    :param gene_seq_type: ['nucl' | 'prot'] : a type of gene_seq file (nucleotide or animo acid)
    :type gene_seq_type: str
    :param gene_seq: Path to a file of target gene in fasta format
    :type gene_seq: str
    :param f_interval: upper region you want to analyze from cds (corresponding to gene_seq), defaults to 0
    :type f_interval: int
    :param r_interval: lower region you want to analyze from cds (corresponding to gene_seq), defaults to 0
    :type r_interval: int
    """
    print_sep = '#' * 40
    out_dir = make_directory(dir_name='target', out_dir=out_dir)
    anno_fmt=['qseqid', 'sseqid', 'stitle', 'evalue', 'slen', 'sstart', 'send', 'qlen', 'qstart', 'qend', 'pident', 'nident', 'mismatch', 'length', 'qseq'] 
    print(print_sep)
    print('Running prodigal.......')
    gff_file, ncds_file, pcds_file = prodigal_py(
        in_ref_fa=ref_fa,
        out_dir=out_dir
    )
    print(print_sep)
    print('Running makeblastdb......')
    blastdb_path = makeblastdb_py(
        seq_type=gene_seq_type,
        ref_fa=gene_seq,
        out_dir=out_dir
    )
    print(print_sep)
    print('Running BLAST program......')
    blast_out = blast_py(
        query_fa=ncds_file,
        seq_type=gene_seq_type,
        blastdb_path=blastdb_path,
        out_dir=out_dir,
        anno_fmt=anno_fmt
    )
    print(print_sep)
    print('Converting files into bed format......')
    out_bed = get_bed_from_blast(
        gff_file=gff_file,
        nucl_fa_file=ncds_file,
        blast_out=blast_out,
        anno_fmt=anno_fmt,
        out_dir=out_dir,
        f_interval=f_interval,
        r_interval=r_interval
    )
    add_fasta_entry(
        ref_fa=ref_fa,
        gff_file=gff_file
    )
    print(print_sep)
    print('Finish.')


#!/usr/bin/env python3
import re, subprocess, os, shutil
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.font_manager import FontProperties
from joblib import Parallel, delayed
from .utility import make_directory, make_dir_files
import vcf
from Bio import SeqIO

def mpileup(bam_files, sample_list, ref_fa, out_dir, bq_threshold=0):
    """Run mpileup to generate vcf files by using bcftools

    :param bam_files: list of bam files which are fixmate, sort, and markdup processed.
    :type bam_files: list
    :param sample_list: sample name list (same length with bam_files).
    :type sample_list: list
    :param ref_fa: path to input reference file in fasta format.
    :type ref_fa: str
    :param out_dir: path to output directory
    :type out_dir: str
    :return: output raw_vcf files list.
    :rtype: list
    """
    raw_vcf_out_dir = os.path.join(out_dir, 'raw_vcf')
    os.mkdir(raw_vcf_out_dir)
    raw_vcf_files = [os.path.join(raw_vcf_out_dir, ll + '.vcf') for ll in sample_list]
    def process(bam_file, vcf_file):
        subprocess.run(' '.join([
            'bcftools', 'mpileup', '--min-BQ', str(bq_threshold), '-Ou', '-a', 'FORMAT/AD,INFO/AD','-f', ref_fa, bam_file, '|',
            'bcftools', 'call', '-mv', '-Ov', '-o', vcf_file
        ]), shell=True)
    Parallel(n_jobs=-1, verbose=3)([delayed(process)(bam_file, vcf_file) for bam_file, vcf_file in zip(bam_files, raw_vcf_files)])        
    return raw_vcf_files

def snpEff_py(out_dir, cds_gff, raw_vcf_files, sample_list):
    """annotate mutation using snpEff

    :param out_dir: path to output directory
    :type out_dir: str
    :param cds_gff: path to gff file to indicate translated region of reference strain (e.g. calculated by prodigal) and must include fasta sequence. 
    :type cds_gff: str
    :param raw_vcf_files: list of path to raw_vcf_files calculated by mpileup()
    :type raw_vcf_files: list
    :param sample_list: sample name list
    :type sample_list: list
    :return: list of vcf files
    :rtype: list
    """
    snps_ref_dir = os.path.join(out_dir, 'ref')
    os.mkdir(snps_ref_dir)
    shutil.copy2(
        src=cds_gff,
        dst=os.path.join(snps_ref_dir, 'genes.gff')
    )
    conf_file = os.path.join(out_dir, 'snpEff.config')
    with open(conf_file, 'w') as f:
        f.write('ref.genome : ref')
    subprocess.run(' '.join([
        'snpEff', 'build', '-gff3', '-c', conf_file,
        '-dataDir', out_dir, 'ref'
    ]), shell=True)
    vcf_dir, vcf_files = make_dir_files(dir_name='annotated_vcf', out_dir=out_dir, prefix_list=sample_list, ext='.vcf')
    def process(raw_vcf_file, vcf_file):
        subprocess.run(' '.join([
            'snpEff', 'ann', '-c', conf_file, '-dataDir', out_dir,
            '-no-downstream', '-no-upstream', '-no-utr',
            '-i', 'vcf', 'ref', raw_vcf_file, '>', vcf_file
        ]), shell=True)       
    Parallel(n_jobs=-1, verbose=3)([delayed(process)(raw_vcf_file, vcf_file) for raw_vcf_file, vcf_file in zip(raw_vcf_files, vcf_files)])        
    return vcf_files

def judge_mut_type(ref, alt):
    """type mutation

    Type 	Name 	Example;
    snp 	Single Nucleotide Polymorphism 	A => T
    mnp 	Multiple Nuclotide Polymorphism 	GC => AT
    ins 	Insertion 	ATT => AGTT
    del 	Deletion 	ACGG => ACG
    complex 	Combination of snp/mnp 	ATTC => GTTA

    :param ref: reference base
    :type ref: str
    :param alt: altenated base
    :type alt: str
    :return: snp: mutation type
    :rtype: str
    """
    if len(ref) == len(alt):
        if len(ref) == 1:
            ans = 'snp'
        else:
            diff_all = sum([aa != bb for aa, bb in zip(ref, alt)])
            if diff_all:
                ans = 'complex'
            else:
                ans = 'mnp'
    else:
        if len(ref) - len(alt) > 0:
            ans = 'del'
        else:
            ans = 'ins'
    return ans

def combine_snps(vcf_df, ref_fa):
    """joint snps to mnp if they are in next amino acide in the same protein, and classify mutation type

    :param vcf_df: dataframe for vcf file (for csv file)
    :type vcf_df: pd.core.frame.DataFrame
    :param ref_fa: path to a file of complete genome or contig in fasta format as a reference
    :type ref_fa: str
    :return: newly annotated dataframe for vcf file (for csv file)
    :rtype: pd.core.frame.DataFrame
    """
    seq_list = [record.seq for record in SeqIO.parse(ref_fa, 'fasta')]
    prot_pos_list = [[int(jj) for jj in re.findall(r'\d+', ll)] for ll in vcf_df.PROT_EFFECT]
    n_chrom_list, n_pos_list, n_ftype_list, n_ref_list, n_alt_list = [vcf_df.CHROM[0]], [vcf_df.POS[0]], [[vcf_df.FTYPE[0]]], [str(vcf_df.REF[0])], [str(vcf_df.ALT[0])]
    n_nucl_effect_list, n_prot_effect_list, n_depth_list, n_frequency_list = [str(vcf_df.NUCL_EFFECT[0])], [str(vcf_df.PROT_EFFECT[0])], [vcf_df.DEPTH[0]], [vcf_df.FREQUENCY[0]]
    n_type_list = [mut_type_select(REF=n_ref_list[0], ALT=n_alt_list[0])]
    for i in range(1, len(vcf_df)):
        if (vcf_df.CHROM[i - 1] == vcf_df.CHROM[i]) & (vcf_df.POS[i] - vcf_df.POS[i - 1] <30) &\
            (len(prot_pos_list[i - 1]) > 0) & (len(prot_pos_list[i]) > 0):
            n_type = mut_type_select(REF=str(vcf_df.REF[i]), ALT=str(vcf_df.ALT[i]))
            if (difference_pos(prot_pos_list[i - 1], prot_pos_list[i]) == 1) &\
                (n_type_list[len(n_type_list) - 1] not in ['ins', 'del']) &\
                (n_type not in ['ins', 'del']):
                gap_seq = str(seq_list[0][(vcf_df.POS[i] - 1):(vcf_df.POS[i] - 1 + len(n_ref_list[len(n_ref_list) - 1]))])
                n_ref_list[len(n_ref_list) - 1] += gap_seq + str(vcf_df.REF[i])
                n_alt_list[len(n_alt_list) - 1] += gap_seq + str(vcf_df.ALT[i])
                n_nucl_effect_list[len(n_nucl_effect_list) - 1] += ';' + str(vcf_df.NUCL_EFFECT[i])
                n_prot_effect_list[len(n_prot_effect_list) - 1] += ';' + str(vcf_df.PROT_EFFECT[i])
                n_depth_list[len(n_depth_list) - 1] = min(n_depth_list[len(n_depth_list) - 1], vcf_df.DEPTH[i])
                n_frequency_list[len(n_frequency_list) - 1] = min(n_frequency_list[len(n_frequency_list) - 1], vcf_df.FREQUENCY[i])
                n_ftype_list[len(n_ftype_list) - 1] = n_ftype_list[len(n_ftype_list) - 1] + [vcf_df.FTYPE[i]]
                n_type_list[len(n_type_list) - 1] = mut_type_select(REF=n_ref_list[len(n_ref_list) - 1], ALT=n_alt_list[len(n_alt_list) - 1])
            else:
                n_chrom_list += [vcf_df.CHROM[i]]
                n_pos_list += [vcf_df.POS[i]]
                n_ftype_list += [[vcf_df.FTYPE[i]]]
                n_ref_list += [str(vcf_df.REF[i])]
                n_alt_list += [str(vcf_df.ALT[i])]
                n_nucl_effect_list += [vcf_df.NUCL_EFFECT[i]]
                n_prot_effect_list += [vcf_df.PROT_EFFECT[i]]
                n_depth_list += [vcf_df.DEPTH[i]]
                n_frequency_list += [vcf_df.FREQUENCY[i]] 
                n_type_list += [mut_type_select(REF=str(vcf_df.REF[i]), ALT=str(vcf_df.ALT[i]))]
        else:
            n_chrom_list += [vcf_df.CHROM[i]]
            n_pos_list += [vcf_df.POS[i]]
            n_ftype_list += [[vcf_df.FTYPE[i]]]         
            n_ref_list += [str(vcf_df.REF[i])]
            n_alt_list += [str(vcf_df.ALT[i])]
            n_nucl_effect_list += [vcf_df.NUCL_EFFECT[i]]
            n_prot_effect_list += [vcf_df.PROT_EFFECT[i]]
            n_depth_list += [vcf_df.DEPTH[i]]
            n_frequency_list += [vcf_df.FREQUENCY[i]]       
            n_type_list += [mut_type_select(REF=str(vcf_df.REF[i]), ALT=str(vcf_df.ALT[i]))]
    for i in range(len(n_ftype_list)):
        if n_prot_effect_list[i] == '':
            n_ftype_list[i] = 'intergenic_region'
        elif not sum([ll != 'synonymous_variant' for ll in n_ftype_list[i]]):
            n_ftype_list[i] = 'synonymous_variant'
        else:
            n_ftype_list[i] = 'missense_variant'
    new_vcf_df = pd.DataFrame({
            'CHROM':n_chrom_list,
            'POS':n_pos_list,
            'TYPE':n_type_list,
            'FTYPE':n_ftype_list,
            'REF':n_ref_list,
            'ALT':n_alt_list,
            'NUCL_EFFECT':n_nucl_effect_list,
            'PROT_EFFECT':n_prot_effect_list,
            'DEPTH':n_depth_list,
            'FREQUENCY':n_frequency_list       
    })
    return new_vcf_df

def mut_type_select(REF, ALT):
    """Classify mutation type

    snp: one base variant. (e.g. C -> G)
    mnp: multi bases variant. (e.g. CA -> GT)
    complex: combination of snps and mnp with a short gap. (e.g. ATGC -> CAGC)
    ins: insertion. (e.g. CG -> CAG)
    del: deletion. (e.g. AT -> A)

    :param REF: reference sequence
    :type REF: str
    :param ALT: alternated sequece
    :type ALT: str
    :return: ['snp', 'mnp', 'complex', 'ins', 'del']
    :rtype: [type]
    """
    if (len(REF) == 1) & (len(ALT) == 1) & (REF != ALT):
        mut_type = 'snp'
    elif len(REF) > len(ALT):
        mut_type = 'del'
    elif len(REF) < len(ALT):
        mut_type = 'ins'
    else:
        if sum([aa == bb for aa, bb in zip(REF, ALT)]) == 0:
            mut_type = 'mnp'
        else:
            mut_type = 'complex'
    return mut_type

def difference_pos(a_list, b_list):
    return min(abs(max(a_list) - min(b_list)), abs(max(b_list) - min(a_list)))

def gen_snps_csv(vcf_files, sample_list, out_dir, min_depth, min_freq, ref_fa):
    """Convert vcf files into csv files and filter them.

    :param vcf_files: list of vcf files which are annotated by snpEff
    :type vcf_files: list
    :param sample_list: sample name list (same length with vcf_files)
    :type sample_list: list
    :param out_dir: path to output directory
    :type out_dir: str
    :param min_depth: SNPs filtering parameter: minimum depth
    :type min_depth: int
    :param min_freq: SNPs filtering parameter: minimum frequency of the SNPs (0.0 - 1.0)
    :type min_freq: float
    :return: raw SNPs.csv file list, filtered SNPs.csv file list
    :rtype: (list, list)
    """
    raw_snps_csv_dir, raw_snps_csv_files = \
        make_dir_files(dir_name='raw_snps', out_dir=out_dir, prefix_list=sample_list, ext='.csv')
    filt_snps_csv_dir, filt_snps_csv_files = \
        make_dir_files(dir_name='filtered_snps', out_dir=out_dir, prefix_list=sample_list, ext='.csv')
    def process(vcf_file, raw_snps_csv_file, filt_snps_csv_file):
        vcf_gen = vcf.Reader(filename=vcf_file)
        vcf_list = list(vcf_gen)
        chrom_list = [ll.CHROM for ll in vcf_list]
        pos_list = [ll.POS for ll in vcf_list]
        dp_list = [ll.INFO['DP'] for ll in vcf_list]
        ad_list = [ll.INFO['AD'] for ll in vcf_list]
        ref_list = [ll.REF for ll in vcf_list]
        alt_list = [ll.ALT for ll in vcf_list]
        ann_list = [ll.INFO['ANN'] if 'ANN' in ll.INFO.keys() else '' for ll in vcf_list]
        alt_list2 = [alt[ad[1:].index(max(ad[1:]))] for ad, alt in zip(ad_list, alt_list)]
        freq_list = [max(ad[1:]) / dp for ad, dp in zip(ad_list, dp_list)]
        type_list = [judge_mut_type(ref=ref, alt=alt) for ref, alt in zip(ref_list, alt_list)]
        ftype_list = [ann[ad[1:].index(max(ad[1:]))].split('|')[1] if ann != '' else '' for ann, ad in zip(ann_list, ad_list)]
        effect_nt_list = [ann[ad[1:].index(max(ad[1:]))].split('|')[9] if ann != '' else '' for ann, ad in zip(ann_list, ad_list)]
        effect_pr_list = [ann[ad[1:].index(max(ad[1:]))].split('|')[10] if ann != '' else '' for ann, ad in zip(ann_list, ad_list)]
        vcf_df = pd.DataFrame({
            'CHROM':chrom_list,
            'POS':pos_list,
            'TYPE':type_list,
            'FTYPE':ftype_list,
            'REF':ref_list,
            'ALT':alt_list2,
            'NUCL_EFFECT':effect_nt_list,
            'PROT_EFFECT':effect_pr_list,
            'DEPTH':dp_list,
            'FREQUENCY':freq_list
        })
        vcf_df.to_csv(raw_snps_csv_file, index=False)
        filt_vcf_df = vcf_df[(vcf_df['DEPTH'] >= min_depth) & (vcf_df['FREQUENCY'] >= min_freq)].reset_index(drop=True)
        filt_vcf_df = combine_snps(vcf_df=filt_vcf_df, ref_fa=ref_fa)
        filt_vcf_df.to_csv(filt_snps_csv_file, index=False)
    Parallel(n_jobs=-1, verbose=3)([delayed(process)(vcf_file, raw_snps_csv_file, filt_snps_csv_file)\
        for vcf_file, raw_snps_csv_file, filt_snps_csv_file\
        in zip(vcf_files, raw_snps_csv_files, filt_snps_csv_files)])
    return raw_snps_csv_files, filt_snps_csv_files

def calc_depth(out_dir, sample_list, bam_files):
    """calculate mapping depth from bam files

    :param out_dir: path to output directory
    :type out_dir: str
    :param sample_list: sample name list (same length with bam_files)
    :type sample_list: list
    :param bam_files: list of bam files path
    :type bam_files: list
    :return: list of depth files (tab-separated)
    :rtype: list
    """
    depth_dir, depth_files = make_dir_files(dir_name='depth', out_dir=out_dir, prefix_list=sample_list, ext='txt')
    def process(bam_file, depth_file):
        subprocess.run(' '.join([
            'samtools', 'depth', bam_file, '>', depth_file
        ]), shell=True)
    Parallel(n_jobs=-1, verbose=3)([delayed(process)(bam_file, depth_file)\
        for bam_file, depth_file in zip(bam_files, depth_files)])
    return depth_files

def get_not_aligned_area(depth_files, left_pos, right_pos, min_depth=1):
    """get list of not aligned region

    :param depth_files: list of depth files calculated by calc_depth
    :type depth_files: list
    :param left_pos: [description]
    :type left_pos: [type]
    :param right_pos: [description]
    :type right_pos: [type]
    :param min_depth: [description], defaults to 1
    :type min_depth: int, optional
    """
    def process(depth_file):
        depth_df = pd.read_table(depth_file, header=None, names=['CHROM', 'POS', 'DEPTH'])
        depth_df[(depth_df['POS'] >= left_pos) & depth_df['POS'] <= right_pos].reset_index(drop=True)
        no_depth_pos_list = []
        tmp_left, tmp_right = 0, 0
        for j in range(left_pos, right_pos+1):
            if (len(depth_df[depth_df['POS'] == j]) == 0)\
                | len(depth_df[(depth_df['POS'] == j) & (depth_df['DEPTH'] < min_depth)]) == 1:
                if j - tmp_right > 1:
                    tmp_left = j
                    tmp_right = j
                    no_depth_pos_list += [[tmp_left, tmp_right]]
                else:
                    tmp_right = j
                    no_depth_pos_list[len(no_depth_pos_list) - 1][1] = tmp_right
        no_depth_pos_list2 = []
        for j in range(len(no_depth_pos_list)):
            no_depth_pos_list2 += [tuple(no_depth_pos_list[j])] 
        return no_depth_pos_list2
    print('getting not aligned region.......')
    no_depth_list = Parallel(n_jobs=-1, verbose=3)([delayed(process)(depth_file) for depth_file in depth_files])
    return no_depth_list

def get_target_snps(filt_snps_csv_file, cds_regions, start_pos, end_pos):
    """filter snps

    In CDS region; rows with 'synonymous_variant' in 'FTYPE' column are removed  
    Out of CDS region; values in 'FTYPE' column converted to 'intergenic_region'

    :param filt_snps_csv_file: path to filtered_snps file generated from gen_snps_csv()
    :type filt_snps_csv_file: str
    :param cds_regions: cds region [(start - 1, end), (start - 1, end), ......]
    :type cds_regions: list
    :param start_pos: start - 1 (target region)
    :type start_pos: int
    :param end_pos: end (target region)
    :type end_pos: int
    :return: processed filtered_snps data frame
    :rtype: pandas.core.frame.DataFrame
    """
    filt_vcf_df = pd.read_csv(filt_snps_csv_file)
    filt_vcf_df = filt_vcf_df[(filt_vcf_df['POS'] > start_pos) & (filt_vcf_df['POS'] <= end_pos)].reset_index(drop=True)
    in_cds = [False] * len(filt_vcf_df)
    pos_list = filt_vcf_df['POS'].values.tolist()
    for i in range(len(in_cds)):
        for j in range(len(cds_regions)):
            if (pos_list[i] > cds_regions[j][0]) & (pos_list[i] <= cds_regions[j][1]):
                in_cds[i] = True
            else:
                pass
    # convert ftype to 'nan' out of cds
    filt_vcf_df['FTYPE'] = [ftype if in_cds_v == True else 'nan' for ftype, in_cds_v in zip(filt_vcf_df['FTYPE'].values.tolist(), in_cds)]
    # remove synonymouse_variant in cds
    filt_vcf_df = filt_vcf_df[filt_vcf_df['FTYPE'] != 'synonymous_variant'].reset_index(drop=True)
    return filt_vcf_df

def ref_dict(input_list, input_dict):
    """input_listをinput_dictに対応した値に変換して出力

    :param input_list: 変換したい値の入ったリスト
    :type input_list: list
    :param input_dict: 変換用辞書。{'org':[], 'new':[]}で記述する。
    :type input_dict: dict
    :return: 変換後のリスト
    :rtype: list
    """
    value_list = []
    for i in range(len(input_list)):
        value_list += [input_dict[input_dict['org']==input_list[i]]['new'].values.tolist()[0]]
    return value_list

def gen_legends(out_dir):
    """generate legends

    :param out_dir: path to output directory
    :type out_dir: str
    :return: output file path
    :rtype: str
    """
    out_file = os.path.join(out_dir, 'legends.svg')
    fig = plt.figure(figsize=(1, 2))
    ax = fig.add_subplot(111)

    # set xlim and ylim
    ax.set_xlim(left=20, right=110)
    ax.set_ylim(bottom=0, top=14)

    # set font properties
    font = FontProperties(family='serif')
    font.set_name('Times New Roman')

    # add ticks
    ax.yaxis.tick_right()
    ax.set_yticks(np.linspace(1, 13, 7))
    ax.set_yticklabels(['del', 'ins', 'complex', 'mnp', 'snp', 'Not alighned region', 'CDS'], fontproperties=font)
    ax.set_xticks([30, 65, 100])
    ax.set_xticklabels(['30', '65', '100'], fontproperties=font)

    # show labels
    ax.set_xlabel('Detection frequency [%]', fontproperties=font)
    
    # show cds arrow
    arrow = patches.FancyArrowPatch((30, 13), (100, 13), mutation_scale=10)
    ax.add_patch(arrow)

    # mark mutations
    shape_list = ['X', 'P', 'd', 's', '.']
    yvalue_list = [1, 3, 5, 7, 9]
    xvalue_list = [30, 65, 100]
    m_size_list = [(6.0 * 0.3) ** 2, (6.0 * 0.65) ** 2, (6.0 * 1.0) ** 2]
    m_col = 'red'
    for i in range(len(yvalue_list)):
        ax.scatter(
            x=xvalue_list,
            y=[yvalue_list[i]] * 3,
            s=m_size_list,
            marker=shape_list[i],
            c=m_col,
            edgecolor = 'k',
            linewidth=0.5,
            alpha=0.5
        )
    
    # show not aligned region
    rect = patches.Rectangle(xy=(30, 10), width=70, height=2, fc='black', alpha=0.5, linewidth=0.)
    ax.add_patch(rect)

    fig.savefig(out_file, bbox_inches='tight')
    return out_file

def make_graph(sample_list, target_bed, filt_snps_csv_files, depth_files, out_dir, min_depth=5):
    """generate graph

    target_bed
        If the values in name column are same, target regions are combined.

    :param sample_list: sample name list
    :type sample_list: list
    :param target_bed: path to target.bed file
    :type target_bed: str
    :param filt_snps_csv_files: list of path to filtered_snps files generated from gen_snps_csv()
    :type filt_snps_csv_files: list
    :param depth_files: list of depth files calculated by calc_depth()
    :type depth_files: list
    :param out_dir: output directory
    :type out_dir: str
    :return: [description]
    :rtype: [type]
    """
    y_value_dict = pd.DataFrame({'org':sample_list, 'new':[2 * (len(sample_list) - ll) - 1 for ll in range(len(sample_list))]})
    m_shape_dict = pd.DataFrame({'org':['snp', 'mnp', 'complex', 'del', 'ins'], 'new':['.', 's', 'd', 'X', 'P']})
    m_col_dict = pd.DataFrame({'org':['nan', 'missense_variant'], 'new':['white', 'red']})
    target_df = pd.read_table(
        target_bed, comment='#',
        names=['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'CDS_Start', 'CDS_End']
    )
    unique_name = list(set(target_df['name'].values.tolist()))
    # グラフへの描画
    strain_count = len(sample_list)
    svg_dir, svg_files = make_dir_files(dir='svg', out_dir=out_dir, prefix_list=unique_name)
    for  i in range(len(unique_name)):
        filt_target_df = target_df[target_df['name'] == unique_name[i]].reset_index(drop=True)
        filt_target_df = filt_target_df.sort_values(by='CDS_Start').reset_index(drop=True)
        cds_regions = [(cds_start, cds_end) for cds_start, cds_end in zip(filt_target_df['CDS_Start'].values.tolist(), filt_target_df['CDS_End'].values.tolist())]
        start_pos = min(filt_target_df['chromStart'].values.tolist())
        end_pos = max(filt_target_df['chromEnd'].values.tolist())
        fig = plt.figure(figsize=(10, 4))
        ax = fig.add_subplot(111)

        # グラフ領域の指定
        ax.set_xlim(
            left=min(filt_target_df['chromStart'].values.tolist()),
            right=max(filt_target_df['chromEnd'].values.tolist())
        )
        ax.set_ylim(bottom=0, top=strain_count*2+2)

        # font設定
        font = FontProperties(family='serif')
        font.set_name('Times New Roman')
        title_font = FontProperties(family='serif')
        title_font.set_name('Times New Roman')
        title_font.set_style('italic')
        title_font.set_size(11)

        # タイトルの追加
        ax.set_title(unique_name[i], fontproperties=title_font)

        # 目盛りの設定
        ax.set_yticks(np.linspace(1, strain_count * 2 + 1, int(strain_count + 1)))
        ax.set_yticklabels(sample_list[::-1] + ['CDS'], fontproperties=font)
        ax.set_yticks(np.linspace(0, strain_count * 2 + 2, int(strain_count * 2 + 3)), minor=True)
        ax.get_xaxis().get_major_formatter().set_scientific(False)
        ax.get_xaxis().get_major_formatter().set_useOffset(False)    
        for x_tick_labels in ax.get_xticklabels():
            x_tick_labels.set_rotation(30)
            x_tick_labels.set_fontproperties(font)
            x_tick_labels.set_horizontalalignment('right')

        # グリッドの表示
        ax.grid(b=True, which='major', axis='x', linestyle='--')
        ax.grid(b=True, which='minor', axis='y')

        # ラベルの表示
        ax.set_xlabel('Position in genome', fontproperties=font)
        ax.set_ylabel('Sample', fontproperties=font)
            
        # cds領域の描画
        arrow_y = strain_count * 2 + 1
        for j in range(len(filt_target_df)):
            if filt_target_df['strand'][j] == '+':
                arrow = patches.FancyArrowPatch(
                    (filt_target_df['CDS_Start'][j], arrow_y),
                    (filt_target_df['CDS_End'][j], arrow_y), 
                    mutation_scale=10
                )
            elif filt_target_df['strand'][j] == '-':
                arrow = patches.FancyArrowPatch(
                    (filt_target_df['CDS_End'][j], arrow_y), 
                    (filt_target_df['CDS_Start'][j], arrow_y), 
                    mutation_scale=10
                )
            ax.add_patch(arrow)
        # mutationの描画
        for j in range(strain_count):
            filt_vcf_df = get_target_snps(
                filt_snps_csv_file=filt_snps_csv_files[j],
                cds_regions=cds_regions,
                start_pos=start_pos,
                end_pos=end_pos
            )
            for k in range(len(filt_vcf_df)):
                y_value = ref_dict([sample_list[j]], y_value_dict)
                x_value = [filt_vcf_df['POS'][k]]
                m_size = [(6.0 * filt_vcf_df['FREQUENCY'][k]) ** 2]
                m_shape = ref_dict([filt_vcf_df['TYPE'][k]], m_shape_dict)[0]
                m_col = ref_dict([filt_vcf_df['FTYPE'][k]], m_col_dict)[0]
                ax.scatter(x=x_value, y=y_value, s=m_size, marker=m_shape, c=m_col, edgecolor = 'k', linewidth=0.5, alpha=0.5)
        # アラインメントされなかった部分の描画
        print('Processing not aligned area...')
        no_depth_list = get_not_aligned_area(
            depth_files=depth_files,
            left_pos=start_pos,
            right_pos=end_pos,
            min_depth=min_depth
        )
        rect_height = 2
        for j in range(len(no_depth_list)):
            for k in range(len(no_depth_list[j])):
                x_value = no_depth_list[j][k][0]
                y_value = ref_dict([sample_list[j]], y_value_dict)[0] - 1
                rect_width = no_depth_list[j][k][1] - x_value
                rect = patches.Rectangle(xy=(x_value, y_value), width=rect_width, height=rect_height, fc='black', ec='black', alpha=0.5, linewidth=0.)
                ax.add_patch(rect)
        # fig.show()
        fig.savefig(svg_files[i], bbox_inches='tight')
        print('Completed saving graph to'+ svg_files[i])
    return svg_files

def gen_graph(out_dir, graph_config_csv, target_bed, ref_fa, cds_gff, min_depth=5):
    """gen_graph main function.

    graph_config.csv file is generated from alignment.make_alignment.
    target.bed is generated from get_target.get_target_position.
        If the values in name column are same, target regions are combined.
    cds.gff file is generated from get_target.prodigal_py.

    :param out_dir: output directory
    :type out_dir: str
    :param graph_config_csv: path to graph_config.csv
    :type graph_config_csv: str
    :param target_bed: path to target.bed
    :type target_bed: str
    :param ref_fa: path to a file of complete genome or contig in fasta format as a reference
    :type ref_fa: str
    :param cds_gff: path to cds.gff
    :type cds_gff: str
    :param min_depth: minimum depth of sequence
    :type min_depth: int
    """
    out_dir = make_directory(dir_name='missense_mut_graph', out_dir=out_dir)
    graph_conf_df = pd.read_csv(graph_config_csv)
    sample_list = graph_conf_df.iloc[:, 0].values.tolist()
    bam_files = graph_conf_df.iloc[:, 1].values.tolist()
    print('generating vcf files.......')
    row_vcf_files = mpileup(
        bam_files=bam_files,
        sample_list=sample_list,
        ref_fa=ref_fa,
        out_dir=out_dir,
        bq_threshold=0
    )
    print('annotating vcf files.......')
    vcf_files = snpEff_py(
        out_dir=out_dir,
        cds_gff=cds_gff,
        raw_vcf_files=row_vcf_files,
        sample_list=sample_list
    )
    print('converting vcf files into csv files.......')
    raw_snps_csv_files, filt_snps_csv_files = gen_snps_csv(
        vcf_files=vcf_files,
        sample_list=sample_list,
        out_dir=out_dir,
        min_depth=min_depth,
        min_freq=0.3,
        ref_fa=ref_fa
    )
    print('calculating depth.......')
    depth_files = calc_depth(
        out_dir=out_dir,
        sample_list=sample_list,
        bam_files=bam_files
    )
    print('making graph.......')
    svg_files = make_graph(
        sample_list=sample_list,
        target_bed=target_bed,
        filt_snps_csv_files=filt_snps_csv_files,
        depth_files=depth_files,
        out_dir=out_dir,
        min_depth=min_depth
    )
    legend_file = gen_legends(
        out_dir=out_dir
    )
    print('Finish.')
    return None

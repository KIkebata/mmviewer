#!/usr/bin/env python3
from .get_target import get_target_position
from .alignment import make_alignment
from .utility import get_option
from .gen_graph import gen_graph
from .__version__ import __version__

def main():
    arg_option = get_option()
    if arg_option.program_name=='get_target':
        get_target_position(
            ref_fa=arg_option.complete_seq,
            out_dir=arg_option.output,
            gene_seq_type=arg_option.gene_seq_type,
            gene_seq=arg_option.gene_sequence,
            f_interval=arg_option.upper_interval,
            r_interval=arg_option.lower_interval
        )
    elif arg_option.program_name=='alignment':
        make_alignment(
            al_config_csv=arg_option.alignment_config_file,
            ref_fa=arg_option.complete_seq,
            out_dir=arg_option.output
        )
    elif arg_option.program_name=='gen_graph':
        gen_graph(
            out_dir=arg_option.output,
            graph_config_csv=arg_option.graph_config,
            target_bed=arg_option.target_bed,
            ref_fa=arg_option.complete_seq,
            cds_gff=arg_option.cds_gff,
            min_depth=arg_option.min_depth
        )

if __name__ == '__main__':
    main()

# MMViewer
Missense Mutation Viewer (Beta-version)

## Author
Kengo Ikebata

## Synopsis

MMViewer can make a graph to grab the overview of mutations, especially missense mutations in CDS.
This program consists of 3 parts. First part is finding target region, 
in which it runs blast to find target region on refrence.fasta.
Second part is mapping, in which it runs bwa mem program to map reads to reference.
Third part is generating graph, in which it classifies mutations' type and depict graph.

## Quick Start
```
% mmviewer get_target -c reference_seq.fasta -g amino_acid.fasta -o Output/directory -t prot

% mmviewer alignment -a Alignment_config.csv -c reference_seq.fasta -o Output/directory

% mmviewer gen_graph -c reference_seq.fasta -o Output/directory -a graph_config.csv -b target.bed -d cds.gff
```


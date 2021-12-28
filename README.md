# MMViewer
Missense Mutation Viewer

## Author
K. Ikebata

## Synopsis

MMViewer make a graph to grab the overview of mutations, especially missense mutations.

## Quick Start
```
% mmviewer get_target -c reference_seq.fasta -g amino_acid.fasta -o Output/directory -t prot

% mmviewer alignment -a Alignment_config.csv -c reference_seq.fasta -o Output/directory

% mmviewer gen_graph -c reference_seq.fasta -o Output/directory -a graph_config.csv -b target.bed -d cds.gff
```


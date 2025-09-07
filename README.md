# ONT long-read methylation demo (NanoMethViz and dmrseq)

In this analysis, methylation differences in the promoter regions of the Peg3, Kcnq1ot1, and Cdkn1c genes on mouse genome chr7 were examined using ONT modBAM files downloaded from Zenodo. 

First, the modBAM files were read using NanoMethViz::ModBamResult() and converted to tabix-indexed TSV format using the modbam_to_tabix() function to enable fast access.

Then, a NanoMethResult object was created and the methylation pattern of the Peg3 gene was visualized using plot_gene(). 

The data was converted to a BSseq object, retaining only CpG regions on chr7 with sufficient coverage. 

The promoter regions within the TSS ±10 kb range of the three relevant genes were selected to narrow down the dataset. 

The dmrseq() function (with 1 permutation) was run to calculate significant DMRs in these regions; p-values, q-values, and average methylation differences were determined and the results were recorded. 

The obtained DMRs were matched with the relevant genes, reported in TSV/CSV format, and biologically meaningful differences were found in Cdkn1c and Kcnq1ot1 in regions where q<0.1. 

Finally, graphs were generated highlighting different methylation regions in green for all genes.

## Run
```bash
Rscript run_nmv_workshop.R

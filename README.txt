This is the ongoing work to put together a complete suite of functions for CoDa analysis of microbiome, transcriptome and metagenome data

Some of the code here is in the form of a R markdown documents, but the stuff you want is in CodaSeq_xxx.tar.gz. This is an R package that can be installed by downloading it and running this command in your R console:

download the tar.gz file, then inside the R console type:
install.packages("path/to/CodaSeq_xxx.tar.gz", type="source", repos=NULL)

Alternatively:

install.packages('devtools')
devtools::install_github('ggloor/CoDaSeq/CoDaSeq')

tongue_saliva.txt, and tongue_vs_cheek.txt are subsets of the HMP example dataset

Other datasets will be added. There will be datasets for transcriptome, metagenome, selex and 16S rRNA gene tag sequencing

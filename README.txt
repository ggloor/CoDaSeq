This is the ongoing work to put together a complete suite of functions for CoDa analysis of microbiome, transcriptome and metagenome data

This is now an R package that can be installed by downloading it and running this command in your R console:

install.packages('devtools')
devtools::install_github('ggloor/CoDaSeq/CoDaSeq')

Alternatively (and you will need to install a lot of dependencies):
download the tar.gz file, then inside the R console type:
install.packages("path/to/CodaSeq_xxx.tar.gz", type="source", repos=NULL)

ak_op.txt are subsets of the HMP example dataset with attached keratinized gingiva and  supragingival plaque, there is an associated genera file as well . 

Other datasets will be added. There will be datasets for transcriptome, metagenome, selex and 16S rRNA gene tag sequencing

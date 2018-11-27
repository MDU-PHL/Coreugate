require(readr)
require(tidyverse)
require(readxl)

source(file.path("{{script_dir}}",'pad.R'))
# get allele profile
alleles <- read_delim(file.path("overall_alleles.tsv"), "\t", escape_double =  F, trim_ws = T)
alleles$FILE <- gsub(".fa", "", alleles$FILE)
# calculate pairwise allelic distance
pad <- distance_kh_C(alleles = alleles, prop = F, remove.NA = T )
write.table(as.data.frame(pad) , file = file.path("distances.tab"), sep = "\t", row.names = T, quote = F)
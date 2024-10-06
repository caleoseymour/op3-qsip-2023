## Read a biom file. Convert to long OTU matrix.
args = commandArgs(trailingOnly = TRUE)
infile = args[1]
outfile = args[2]

OTUs = data.table::fread(infile)
colnames(OTUs)[1] = 'rn'

mcounts = data.table::melt(OTUs, id.var = 'rn')[value != 0,]
data.table::fwrite(mcounts, outfile, sep='\t', col.names=FALSE)
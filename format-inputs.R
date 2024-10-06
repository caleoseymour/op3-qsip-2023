## A lot of studies with different setups. merge them.
runOP3 = TRUE

soilmap = c('Andesite' = 'AN', 'Basalt' = 'BS', 'Granite' = 'GR')
soilmap = c('Andesite' = 'Andesite', 'Basalt' = 'Basalt', 'Granite' = 'Granite', 'GR' = 'Granite', 'AN' = 'Andesite', 'BS' = 'Basalt')

formatTaxonomy = function(taxtable, ranks = c('Domain','Phylum','Class','Order','Family','Genus','Species'))
## Format a taxonomy table from Qiime2 for use in R.
{
    splittaxa = strsplit(taxtable[,Taxon], ';')
    reranktaxa = lapply(splittaxa, function(x) c(gsub('.*__','', x), rep('', length(ranks)-length(x))))
    taxmat = data.table::as.data.table(do.call('rbind', reranktaxa)[,1:length(ranks)])
    colnames(taxmat) = c('Domain','Phylum','Class','Order','Family','Genus','Species')
    out = cbind(taxtable[,.(id = get('Feature ID'))],taxmat)
    data.table::setkey(out, id)
    return(out)
}

## 18O samples
oDirs = c('2014-07-02_Bri_16S', '/2014-08-07_Bri_16S_plate3', '2014-08-27_Bri_Mau_16S', '2014-07-10_Bri_16S_run2_rerun', '2014-08-14_Bri_16S_plate4')
oMapFs = list.files(oDirs, pattern = '.*map.*', full.name=TRUE)
oMaps = lapply(oMapFs, data.table::fread, header=TRUE, sep='\t')
oCols = unique(unlist(lapply(oMaps, colnames)))

## Filter metadata.  Add in the blocking factor (week).
oMD = do.call('rbind', lapply(oMaps, function(x) x[,mget(oCols, ifnotfound=NA)]))
oMDf = oMD[,.(sampleid = get('#SampleID'), soil = Soil, substrate = Substrate,
    O = IsotopeTreat, C = '12C', f = X16S_copies, x = Density, j = Sample, t = Week,
    block = paste(soilmap[Soil], Substrate, Week, sep='-'))]


tax = formatTaxonomy(data.table::fread('concatenated-table-18O-taxonomy.tsv', sep='\t', header=TRUE))

if (runOP3 == TRUE)
{
    op3 = formatTaxonomy(data.table::fread('concatenated-table-18O-op3-taxonomy-2.tsv', sep='\t', header=TRUE))
    oTax = rbind(tax[!Class == 'Omnitrophia',], op3[tax[Class == 'Omnitrophia',id],])[!duplicated(id),]
} else {
    oTax = tax
}

oSVs = data.table::fread('concatenated-table-18O-long.tsv', sep='\t',header=FALSE)[!V1 %in% oTax[Domain == 'Unassigned', id], .(SV = V1, SAM = V2, n = V3)]
oSVs[, p := n/sum(n), by=SAM] 

## 13C samples
cDirs = c('SEQ_plate_G', 'SEQ_plates_C-F', 'SEQ_plates_A-B', 'whole_DNA')
cMapFs = list.files(cDirs, pattern = '.*map.*', full.name=TRUE)
cMaps = lapply(cMapFs, data.table::fread, header=TRUE, sep='\t')
cCols = unique(unlist(lapply(cMaps, colnames)))

## Filter metadata. Add in the blocking factor (replicate).
cMD = do.call('rbind', lapply(cMaps, function(x) x[,mget(cCols, ifnotfound=NA)]))
#cMD[Substrate == 'control', Substrate := 'NoC']
cMDf = cMD[,.(sampleid = get('#SampleID'), soil = soilmap[Soil], substrate = Substrate,
    O = O.iso, C = C.iso, f = avg_16S_copies, x = Density_g_ml, j = Tube+(10000)*Replicate, t = 6,
    block = paste(soilmap[Soil], Substrate, 6, sep='-'))]

tax = formatTaxonomy(data.table::fread('concatenated-table-13C-taxonomy.tsv', sep='\t', header=TRUE))

if (runOP3 == TRUE)
{
    op3 = formatTaxonomy(data.table::fread('concatenated-table-13C-op3-taxonomy-2.tsv', sep='\t', header=TRUE))
    cTax = rbind(tax[!Class == 'Omnitrophia',], op3[tax[Class == 'Omnitrophia',id],])[!duplicated(id),]
} else {
    cTax = tax
}

cSVs = data.table::fread('concatenated-table-13C-long.tsv', sep='\t',header=FALSE)[!V1 %in% cTax[Domain == 'Unassigned', id], .(SV = V1, SAM = V2, n = V3)]
cSVs[, p := n/sum(n), by=SAM] 


## Write to file.
data.table::fwrite(cMDf, 'merged-13C-metadata.tsv',sep='\t')
data.table::fwrite(oMDf, 'merged-18O-metadata.tsv',sep='\t')
data.table::fwrite(oSVs, 'merged-table-18O-long.tsv', sep='\t')

data.table::fwrite(oTax, 'merged-18O-taxonomy.tsv',sep='\t')
data.table::fwrite(cTax, 'merged-13C-taxonomy.tsv',sep='\t')
data.table::fwrite(cSVs, 'merged-table-13C-long.tsv', sep='\t')
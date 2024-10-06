#install.packages('devtools')
#install.packages('BiocManager')

# install phyloseq and qsip using the utilities on devtools and BiocManagerK
# BiocManager::install('phyloseq')
# devtools::install_github('bramstone/qsip')

## Cale Seymour
## 2021
## Process Bri Finley's qSIP data.

## This script runs qSIP on all taxonomic ranks of the data.

library('phyloseq')
library('qsip')

## Enriched = AFE conf int is above 0.
cifun = function(x)
{
    error = qnorm(0.975)*sd(x)/sqrt(length(x))
    a = mean(x)
    l = a - error
    r = a + error
    
    return(l > 0 & r > 0)
}


## Read in OTU and taxonomy tables.
cSV = data.table::fread('merged-table-13C-long.tsv',sep='\t', header=TRUE)
cTax = data.table::fread('merged-13C-taxonomy.tsv', sep='\t', header=TRUE)
data.table::setkey(cSV, SV)
data.table::setkey(cTax, id)
cDT = cTax[cSV,]

oSV = data.table::fread('merged-table-18O-long.tsv',sep='\t', header=TRUE)
oSV[,SAM := paste0(SAM, '-2014')]
oTax = data.table::fread('merged-18O-taxonomy.tsv', sep='\t', header=TRUE)
data.table::setkey(oSV, SV)
data.table::setkey(oTax, id)
oDT = oTax[oSV,]

DT = unique(rbind(oDT, cDT))
tDT = DT[,.(id, SAM, n, p)]
data.table::setkey(tDT, id)
Tax = unique(rbind(oTax, cTax))

## import and merge metadata
cMD = data.table::fread('merged-13C-metadata.tsv', sep='\t', header=TRUE)[!is.na(x),]
cMD[C == '', C := '12C']
cSM = as.data.frame(cMD)
rownames(cSM) = cMD[,sampleid]

oMD = data.table::fread('merged-18O-metadata.tsv', sep='\t', header=TRUE)[!is.na(x),]
oMD[,sampleid := paste0(sampleid, '-2014')]
oSM = as.data.frame(oMD)
rownames(oSM) = oMD[,sampleid]

MD = rbind(oSM, cSM)
mergedata = data.table::as.data.table(MD, keep.rownames=TRUE)
data.table::setkey(mergedata, rn)


results = do.call('rbind',lapply(2:(ncol(Tax)),
function(RN)
{
    RANK = colnames(Tax)[RN]
    PATH = colnames(Tax)[2:RN]
    
    FT = Tax[,.(taxon=paste0(mget(PATH),collapse=';')), by=id]
    data.table::setkey(FT, id)
    FTDT = tDT[FT,][,.(p = sum(p)),by=.(taxon,SAM)]
    sFT = data.table::dcast(FTDT, formula = taxon ~ SAM, value.var = 'p')
    mFT = as.matrix(sFT[,-1])
    rownames(mFT) = sFT[,taxon]
    mFT[is.na(mFT)] = 0
    
    TT = matrix(rownames(mFT), ncol=1)
    rownames(TT) = TT[,1]
    colnames(TT)[1] = 'Taxon'

    
    PS = phyloseq(tax_table(TT), otu_table(mFT, taxa_are_rows = TRUE), sample_data(MD))
    
    ## Separate into O and C based on experimental blocking factor
    oPS =  specify_qsip(PS, density='x', abund='f', rep_id='sampleid',
    timepoint='t', rep_group='block', iso='18O', iso_trt='O')
    oPSp = calc_excess(oPS, separate_label=T, filter=T, correction=T)
    
    cPS =  specify_qsip(PS, density='x', abund='f', rep_id='sampleid',
    timepoint='t', rep_group='block', iso='13C', iso_trt='C')
    cPSp = calc_excess(cPS, separate_label=T, filter=T, correction=T)
    
    cqm = data.table::as.data.table(qsmelt(cPSp, include='excess|label', regex=T))
    data.table::setkey(cqm, sampleid)
    cqm[mergedata, isotope := C]
    cqm[,Rank := RANK]
    
    oqm = data.table::as.data.table(qsmelt(oPSp, include='excess|label', regex=T))
    data.table::setkey(oqm, sampleid)
    oqm[mergedata, isotope := O]
    oqm[,Rank := RANK]
    
    ## Combine.
    qm = rbind(cqm[isotope == '13C',], oqm[isotope == '18O',])[!is.na(atom_excess),.(tax_id, Rank, sampleid, atom_excess, block)]
    data.table::setkey(qm, tax_id)
    
    ## Filter results on confidence interval for the taxon.
    qmt = qm[!is.na(atom_excess),][,cifun(atom_excess),by=.(tax_id, block)][V1 == TRUE,.(tax_id, block)]
    data.table::setkeyv(qmt, c('tax_id','block'))
    data.table::setkeyv(qm, c('tax_id','block'))
    
    qmf = qm[qmt,]
    
    qmf
}))


data.table::fwrite(results, 'ranked-qsip-results-op3.tsv', sep='\t')
res_op3trimmed = results[!(Rank != 'Phylum' & grepl('^Bacteria;Omnitrophota;', tax_id)),]
data.table::fwrite(res_op3trimmed, 'ranked-qsip-results.tsv', sep='\t')

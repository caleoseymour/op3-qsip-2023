#install.packages('devtools')
#install.packages('BiocManager')

# install phyloseq and qsip using the utilities on devtools and BiocManagerK
# BiocManager::install('phyloseq')
# devtools::install_github('bramstone/qsip')


## Cale Seymour
## 2021
## Process Bri Finley's qSIP data.

## This script runs qSIP on nonpredatory vs predatory bacterial taxa.
## Plot the file afterward.

library('phyloseq')
library('qsip')
library('ggplot2')
library('dplyr')
library('tidyverse')
library('grid')
library('data.table')
library(showtext)
## add the Arial font
font_add("Arial", regular = "arial.ttf",
    bold = "arialbd.ttf", italic = "ariali.ttf", bolditalic = "arialbi.ttf")

library('perm')                    # for permTS()

source("rnc_ggplot2_border_themes_2013_01.r")

source('tax_glom_fast.R')
source('/mnt/d/Ubuntu/science/OP3/try-5/taxonomy/read-taxonomy.R')

fixed = data.table::fread('qsip-groups.txt', sep='\t')

## Enriched = AFE conf int is above 0.
cifun = function(x)
{
    error = qnorm(0.975)*sd(x)/sqrt(length(x))
    a = mean(x)
    l = a - error
    r = a + error
    
    return(l > 0 & r > 0)
}

add_group_taxa = function(taxa, groups)
{
    invisible(apply(groups, 1,
        function(x)
        {
            newids = taxa[get(x[2]) == x[1], id]
            taxa[id %in% newids, Group := x[3]]
        }))
    return(0)
}

## Function for filtering boxplot outliers.
filter_lims <- function(x){
  l <- boxplot.stats(x)$stats[1]
  u <- boxplot.stats(x)$stats[5]

  for (i in 1:length(x)){
    x[i] <- ifelse(x[i]>l & x[i]<u, x[i], NA)
  }
  return(x)
}

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
Tax = unique(rbind(oTax, cTax))

## Agglomerate fixed, predatory taxa at a certain level.
fixed_table = do.call('rbind',
    apply(fixed, 1,
        function(x)
        {
            rl = which(colnames(Tax) == x[2])
            xx = DT[get(x[2]) == x[1],]
            rx = xx[,.(id = x[1], Taxon = x[1], n = sum(n), p = sum(p)), by = c('SAM')]
            return(rx)
        }))
        
fixed_table_ids = unlist(apply(fixed, 1, function(x) DT[get(x[2]) == x[1], unique(id)]))
op3_table = DT[Phylum == 'Omnitrophota',]
op3_ids = op3_table[,unique(id)]

np_table = DT[!(id %in% union(fixed_table_ids, op3_ids)),][,.(id = 'Nonpredator', n = sum(n), p = sum(p), Taxon = 'Nonpredator'), by = SAM]

FT = rbind(fixed_table, np_table)

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

rot = 0
invisible(lapply(3:(ncol(Tax)-1),
function(RN)
{
    RANK = colnames(Tax)[RN]
    rop3_table = op3_table[,.(id = get(RANK), Taxon = get(RANK), p = sum(p), n = sum(n), facet = 'Omnitrophota'), by=c(RANK, 'SAM')]
    rop3_table[,c(RANK) := NULL]
    
    op3taxa = DT[Phylum == 'Omnitrophota', unique(get(RANK))]
    op3taxa = op3taxa[op3taxa != '']
    
    rFT = rbind(rop3_table, FT[,.(id, p, n, Taxon, facet = 'Other taxa'), by='SAM'])[Taxon != '',]
    sFT = data.table::dcast(rFT, formula = id ~ SAM, value.var = 'p')
    
    mFT = as.matrix(sFT[,-1])
    rownames(mFT) = sFT[,id]
    mFT[is.na(mFT)] = 0
    
    TT = as.matrix(rFT[!duplicated(id),Taxon])
    rownames(TT) = rFT[!duplicated(id),id]
    colnames(TT) = 'Taxon'
    
    PS = phyloseq(tax_table(TT), otu_table(mFT, taxa_are_rows = TRUE), sample_data(MD))
    
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
    
    qm = rbind(cqm[isotope == '13C',], oqm[isotope == '18O',])[!is.na(atom_excess),.(tax_id, Rank, sampleid, atom_excess, block)]
    data.table::setkey(qm, tax_id)
    
    qmt = qm[!is.na(atom_excess),][,cifun(atom_excess),by=.(tax_id, block)][V1 == TRUE,.(tax_id, block)]
    data.table::setkeyv(qmt, c('tax_id','block'))
    data.table::setkeyv(qm, c('tax_id','block'))
    
    qmf = qm[qmt,]
    
    nps = qmf[tax_id == 'Nonpredator',]
    colnames(nps)[colnames(nps) == 'atom_excess'] = 'natom_excess'
    data.table::setkey(nps, sampleid)
    
    ps = qmf[tax_id != 'Nonpredator',]
    data.table::setkey(ps, sampleid)
    ps[nps, np := natom_excess]
    ps[,ratio := ((atom_excess-np)/np)]
    
    
    data.table::setkey(ps, tax_id)
    data.table::setkey(fixed, lineage)
    ps[fixed, group := group]
    ps[tax_id %in% op3taxa, group := 'Omnitrophota']
    
    data.table::setkey(ps, sampleid)
    ps[mergedata, c(colnames(mergedata)) := mget(colnames(mergedata))]    

    ps[,isotope := factor(ifelse(O == '18O' & C == '12C', 'O', 'C'))]
    ps[substrate == 'control', substrate := 'NoC']
    
    ## Filter outliers to clean up the data a little bit.
    pd = ps %>% group_by(substrate, isotope, tax_id) %>% mutate(ratio2 = filter_lims(ratio)) %>% data.table::as.data.table()
    
    substrate_transform = c('Exu' = 'Exudate', 'glucose' = 'Glucose', 'Lit' = 'Litter', 'NoC' = 'No Carbon', 'oxalic' = 'Oxalate', 'control' = 'No Carbon')
    pd[,substrate := factor(substrate_transform[as.character(substrate)], levels = c('No Carbon', 'Glucose', 'Oxalate', 'Exudate', 'Litter'))]
    
    
    pd[,f_taxid := factor(tax_id)]
    pd[,f_substrate := factor(substrate)]
    pd[,f_isotope := factor(isotope)]
    
    po = unique(pd[,.(group, f_taxid, f_isotope, f_substrate, size=.N),by=.(group, f_taxid, f_isotope, f_substrate)][order(group, f_taxid, f_isotope, f_substrate),])[size > 7,]
    new_taxid = po[,as.character(f_taxid)] != c('', po[-nrow(po), as.character(f_taxid)])
    po$taxid_offset = 2 * cumsum(new_taxid)
    new_isotope = po[,as.character(f_isotope)] != c('', po[-nrow(po), as.character(f_isotope)])
    po$isotope_offset = 2 * cumsum(new_isotope | new_taxid)
   
    po[,xv := 1:.N + taxid_offset + isotope_offset]
    data.table::setkeyv(po, c('group','f_taxid','f_isotope','f_substrate'))
    data.table::setkeyv(pd, c('group','f_taxid','f_isotope','f_substrate'))
    pd[po, xval2 := xv]
    pd = pd[order(xval2, ratio2),]
    pd[,xp := (xval2 - 0.4) + (1:length(x)/length(x))*0.8, by=xval2]
    outliers = pd[is.na(ratio2) & !is.na(ratio), .(hi=sum(ratio > 0), lo = sum(ratio < 0)),by=.(xval2, group)]
    
    isotope_labels = po[,.(xmin = min(xv), xmax = max(xv)), by=.(group, f_taxid, f_isotope)]
    isotope_labels[,xc := xmin + (xmax - xmin) / 2]
    isotope_labels[f_isotope == 'C', label := '13C']
    isotope_labels[f_isotope == 'O', label := '18O']
    
    taxon_labels = po[,.(xmin = min(xv), xmax = max(xv)), by=.(group, as.character(f_taxid))]
    taxon_labels[,color := ifelse(.I %% 2 == 0, 'white', '#DDE3FF')]
    taxon_labels[,xc := xmin + (xmax - xmin) / 2]
    message('calc')
    tp = ggplot(pd) +
        geom_boxplot(aes(x = xval2, y = ratio2, fill=substrate, group=xval2), color='black', outlier.shape = NA,
                na.rm = TRUE, coef = 1.58)
    ymin = ggplot_build(tp)$layout$panel_params[[1]]$y.range[1]
    ymax = ggplot_build(tp)$layout$panel_params[[1]]$y.range[2]
    
    plotsize = (ymax - ymin)
    ybreaks = ggplot_build(tp)$layout$panel_params[[1]]$y.sec$breaks
    
    
    
    pd[,aovs := paste0(group, '|', tax_id)]
    # ## Run the ANOVA
    # model = lm(data = pd, ratio ~ aov1)
    # ano = aov(model)
    # thsd = TukeyHSD(ano, conf.level = 0.95)$aovl %>% round(3) %>% data.table::as.data.table(keep.rownames = TRUE)
    # thsd[,signif := '']
    # thsd[get('p adj') < 0.05, signif := '*']
    # fwrite(thsd, paste0('qsip.tukey-',RANK, '.txt'), sep = '\t')
    
    # Get the medians of each group
    medians = pd[,median(na.omit(ratio)),by=aovs]
    fwrite(medians, paste0('qsip.medians-',RANK, '.txt'), sep = '\t')

    # ## Significantly different than 0?
    pts = pd[,.(p.value = permTS(as.numeric(na.omit(ratio)), 0, alternative="greater", method="exact.mc", control=permControl(nmc=10^3))$p.value), by=aovs]
    pts[,signif := '']
    pts[p.value < 0.05, signif := '*']
    fwrite(pts, paste0('qsip.permtest-',RANK, '.txt'), sep = '\t')
    
    message('tp')
    pg = ggplot(pd) +
    
        ## Draw background and major grid
        geom_rect(data=taxon_labels, aes(xmin = xmin-2, xmax = xmax+2), ymin = ymin, ymax = ymax + plotsize * 0.025, fill = taxon_labels[,color]) +
        geom_hline(yintercept = ybreaks[ybreaks != 0], linetype = 'longdash', size = 0.5, color= 'grey') +
        geom_hline(yintercept = 0, linetype = 'longdash', color= 'grey') +
        
        ## Draw data points
        geom_boxplot(aes(x = xval2, y = ratio, fill=substrate, group=xval2),
            shape=21, color='black', outlier.shape = NA, width = 0.95,
                na.rm = TRUE, coef = 1.58, show.legend=FALSE) +
        geom_jitter(aes(x = xval2, y = ratio, fill=substrate), shape=21, color='black', width = 0.4, height = 0, size=1.5, alpha=0.5) +
        #geom_point(aes(x = xp, y = ratio, fill=substrate), shape=21, color='black', size=1.5, alpha=0.5) +
        
        ## Draw x axis
        #geom_rect(data=taxon_labels, aes(xmin = xmin-2, xmax = xmax+2), ymin = ymin - (ymax - ymin)*0.05, ymax = ymin, fill = taxon_labels[,color]) +
        geom_rect(data=taxon_labels, aes(xmin = xmin-2, xmax = xmax+2), ymin = ymin - (ymax - ymin)*0.05, ymax = ymin, fill = 'grey') +
        #geom_hline(yintercept = ymin, color='grey', size=0.5) +
        geom_text(data=taxon_labels, aes(x=xc, label=as.character), y = ymin - plotsize*0.05, hjust=1, vjust=0.75, size=2.5, angle = rot + 45, family = 'Arial') +
        geom_text(data=isotope_labels, aes(x=xc, label=label), y = ymin - plotsize*0.025, hjust=0.5, vjust=0, size=2.5, angle = rot, family = 'Arial') +
        #geom_text(data=outliers, aes(x=xval2, label=lo), y = ymin, hjust=0.5, vjust=0, size=2, angle = 90) +
        #geom_text(data=outliers, aes(x=xval2, label= hi), y = ymax, hjust=0.5, vjust=0, size=2, angle = 90) +
        geom_segment(data = isotope_labels, aes(x = xmin - 1, xend = xmax + 1), y = ymin - plotsize*0.005, yend = ymin - plotsize*0.005) +
        #scale_y_continuous(breaks = c(-10, 0,10,50,100,200)/100, labels = paste0(c(-10, 0,10,50,100,200), '%'), limits = c(-100,200)/100) +
        scale_y_continuous(labels = scales::percent, limits=c(ymin - plotsize * 0.05, ymax + plotsize * 0.025), expand = c(0,0)) +
        coord_cartesian(clip = "off") +
        
        scale_fill_manual(name = 'Substrate',
            values = c('ghostwhite', 'darkorchid1', 'palegreen1', 'olivedrab3', 'seagreen3'),
            labels = c('No Carbon', 'Glucose', 'Oxalate', 'Exudate', 'Litter')) +
        guides(fill = guide_legend(override.aes = list(alpha = 1, size= 2), label.position = "top")) +

        #geom_hline(yintercept = -100/100, color='gray', width=0.5) +
        #scale_shape_manual(values=c(21, 211S), breaks=c('C','O')) +
        facet_grid(. ~ group, scales='free_x',space='free_x') +
        theme_bw() +
        ## Relative difference
        labs(y = 'Difference in Isotope Assimilation\n(P-NP)/NP * 100%') +
        theme(legend.position=c(0,-.025),
            legend.justification = c(0,0.9),
            plot.margin = unit(c(0.02,0.02,4.1275/2,1/4), 'cm'),
            panel.grid.minor.x = element_blank(),
            panel.grid.minor.y = element_blank(),
            panel.grid.major.x = element_blank(), 
            panel.grid.major.y = element_blank(), 
            panel.spacing = unit(0.1, "cm"),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_text(angle=90, vjust = 0, hjust = 0.5, size = 10, face = 'bold', family = 'Arial'),
            axis.text.y = element_text(angle=90, vjust = 0, hjust = 0.5, size = 10, family = 'Arial'),
            panel.border = theme_border(type = c("right","left"), size=0.5),
            strip.background.x = element_rect(fill = 'grey', color='grey'),
            strip.background.y = element_rect(fill = 'grey', color='grey'),
            strip.text.x = element_text(angle = rot, size = 10, face='bold', family = 'Arial'),
            legend.text = element_text(angle = 90, vjust = 0.5, hjust = 0, size = 8, family = 'Arial'),
            legend.title = element_text(angle = 90, size = 8, vjust = 1, hjust = 0.15, face='bold', family = 'Arial'),
            legend.direction = "horizontal",
            legend.spacing.y = unit(0, 'mm'),
            legend.spacing.x = unit(0, 'mm'),
            legend.key.height = unit(5, 'mm'),
            legend.key.width = unit(3, 'mm'),
            )
        
    
    # pdf(paste0('qsip-',RANK, '.pdf'), width = 11, height = 6.5)
        # plot(pg)
    # dev.off()
    fname = paste0('qsip-',RANK, '.svg')
    svglite::svglite(fname, width = 11, height = 6.5)
        plot(pg)
    dev.off()
    
    fname2 = paste0('qsip-',RANK, '-mod.svg')

    cmd = paste0('sed -e "s/Liberation Sans/Arial/g" ', fname ,' | sed -e "s/textLength/NOTEXTLENGTH/g" > ', fname2)
    message(cmd)
    system(cmd)
}))

#res = do.call('rbind', results)
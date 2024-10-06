library('phyloseq')
library('qsip')
library('ggplot2')
library('dplyr')
library('tidyverse')
library('grid')
source("rnc_ggplot2_border_themes_2013_01.r")
# The more SROs there are, the less available organic carbon pools are.

# Run the SROs

# Granite = -;
# Basalt = 50;
# Andesite = 78;

# More SROs = more oligotrophic. Plot OP3 AFE relative to SROs when utilization
# occurs.

# Granite = -;
# Basalt = 50;
# Andesite = 78;

filter_lims <- function(x){
  l <- boxplot.stats(x)$stats[1]
  u <- boxplot.stats(x)$stats[5]

  for (i in 1:length(x)){
    x[i] <- ifelse(x[i]>l & x[i]<u, x[i], NA)
  }
  return(x)
}

fixed = data.table::fread('qsip-groups.txt', sep='\t')

qsip = data.table::fread('ranked-qsip-results-op3.tsv', sep='\t')
qsip[,taxon := gsub('.*;','',tax_id)]
qsip[,c('soil', 'substrate', 'timepoint') := data.table::as.data.table(do.call('rbind', strsplit(block, '-')))]

substrate_transform = c('Exu' = 'Exudate', 'glucose' = 'Glucose', 'Lit' = 'Litter', 'NoC' = 'No Carbon', 'oxalic' = 'Oxalate', 'control' = 'No Carbon')
qsip[,substrate := factor(substrate_transform[as.character(substrate)], levels = c('No Carbon', 'Glucose', 'Oxalate', 'Exudate', 'Litter'))]

data.table::setkeyv(qsip, c('taxon', 'Rank'))
data.table::setkeyv(fixed, c('lineage', 'rank'))

Pqsip = qsip[fixed,]
Pqsip[,rank := NULL]

Oqsip = qsip[grepl('^Bacteria;Omnitrophota', tax_id),]
Oqsip[, group := 'Omnitrophota']


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

## Import metadata
cMD = data.table::fread('merged-13C-metadata.tsv', sep='\t', header=TRUE)[!is.na(x),]
cMD[C == '', C := '12C']
cSM = as.data.frame(cMD)
rownames(cSM) = cMD[,sampleid]

oMD = data.table::fread('merged-18O-metadata.tsv', sep='\t', header=TRUE)[!is.na(x),]
oMD[,sampleid := paste0(sampleid, '-2014')]
oSM = as.data.frame(oMD)
rownames(oSM) = oMD[,sampleid]

MD = rbind(oSM, cSM)
MD = data.table::as.data.table(MD, keep.rownames=TRUE)
data.table::setkey(MD, rn)

Fqsip = rbind(Oqsip, Pqsip)

data.table::setkey(Fqsip, sampleid)

Fqsip[MD, C := C]
Fqsip[MD, O := O]

invisible(lapply(3:(ncol(Tax)-1),
function(RN)
{
    RANK = colnames(Tax)[RN]
    ps = Fqsip[group != 'Omnitrophota' | Rank == RANK & taxon != '',][!is.na(sampleid),]
    ps[,isotope := factor(ifelse(O == '18O' & C == '12C', 'O', 'C'))]
    
    pd = ps %>% group_by(substrate, soil, isotope, taxon) %>% mutate(atom_excess2 = filter_lims(atom_excess)) %>% data.table::as.data.table()
    
    
    pd[,f_taxid := factor(taxon)]
    pd[,f_substrate := substrate]
    pd[,f_isotope := factor(isotope)]
    pd[,f_soil := factor(soil, levels = c('Granite', 'Basalt', 'Andesite'))]
    
    po = unique(pd[,.(group, f_taxid, f_soil, f_isotope, f_soil, size=.N), by=.(group, f_taxid, f_soil, f_isotope, f_soil)][order(group, f_taxid, f_isotope, f_soil),])[size > 7,]
    new_taxid = po[,as.character(f_taxid)] != c('', po[-nrow(po), as.character(f_taxid)])
    po$taxid_offset = 2 * cumsum(new_taxid)
    new_isotope = po[,as.character(f_isotope)] != c('', po[-nrow(po), as.character(f_isotope)])
    po$isotope_offset = 2 * cumsum(new_isotope | new_taxid)
   
    po[,xv := 1:.N + taxid_offset + isotope_offset]
    data.table::setkeyv(po, c('group','f_taxid','f_isotope','f_soil'))
    data.table::setkeyv(pd, c('group','f_taxid','f_isotope','f_soil'))
    pd[po, xval2 := xv]
    pd = pd[order(xval2, atom_excess2),]
    pd[,xp := (xval2 - 0.4) + (1:length(xval2)/length(xval2))*0.8, by=xval2]
    outliers = pd[is.na(atom_excess2) & !is.na(atom_excess), .(hi=sum(atom_excess > 0), lo = sum(atom_excess < 0)),by=.(xval2, group)]
    
    isotope_labels = po[,.(xmin = min(xv), xmax = max(xv)), by=.(group, f_taxid, f_isotope)]
    isotope_labels[,xc := xmin + (xmax - xmin) / 2]
    isotope_labels[f_isotope == 'C', label := '13C']
    isotope_labels[f_isotope == 'O', label := '18O']
    
    isotope_labels2 = data.table::copy(isotope_labels)
    
    isotope_labels[,timepoint := 1]
    isotope_labels2[,timepoint := 6]
    
    isotope_labels = rbind(isotope_labels, isotope_labels2)
    
    taxon_labels = po[,.(xmin = min(xv), xmax = max(xv)), by=.(group, as.character(f_taxid))]
    taxon_labels[,color := ifelse(.I %% 2 == 0, 'white', '#DDE3FF')]
    taxon_labels[,xc := xmin + (xmax - xmin) / 2]
    
    taxon_labels2 = data.table::copy(taxon_labels)
    
    taxon_labels[,timepoint := 1]
    taxon_labels2[,timepoint := 6]
    
    taxon_labels = rbind(taxon_labels, taxon_labels2)
    
    
    
    message('calc')
    tp = ggplot(pd) +
        geom_boxplot(aes(x = xval2, y = atom_excess2, fill=soil, group=xval2), color='black', outlier.shape = NA,
                na.rm = TRUE, coef = 1.58)
    ymin = ggplot_build(tp)$layout$panel_params[[1]]$y.range[1]
    ymax = ggplot_build(tp)$layout$panel_params[[1]]$y.range[2]
    
    plotsize = (ymax - ymin)
    ybreaks = ggplot_build(tp)$layout$panel_params[[1]]$y.sec$breaks
    
    
    rot = 0
    
    pg = ggplot(pd) +
    
        ## Draw background and major grid
        geom_rect(data=taxon_labels, aes(xmin = xmin-2, xmax = xmax+2), ymin = ymin, ymax = ymax + plotsize * 0.025, fill = taxon_labels[,color]) +
        geom_hline(yintercept = ybreaks[ybreaks != 0], linetype = 'longdash', size = 0.5, color= 'grey') +
        geom_hline(yintercept = 0, linetype = 'longdash', color= 'grey') +
        
        ## Draw data points
        geom_boxplot(aes(x = xval2, y = atom_excess, fill=soil, group=xval2),
            shape=21, color='black', outlier.shape = NA, width = 0.95,
                na.rm = TRUE, coef = 1.58, show.legend=FALSE) +
        geom_jitter(aes(x = xval2, y = atom_excess, fill=soil), shape=21, color='black', width = 0.4, height = 0, size=1.5, alpha=0.5) +
        #geom_point(aes(x = xp, y = atom_excess, fill=substrate), shape=21, color='black', size=1.5, alpha=0.5) +
        
        ## Draw x axis
        #geom_rect(data=taxon_labels, aes(xmin = xmin-2, xmax = xmax+2), ymin = ymin - (ymax - ymin)*0.05, ymax = ymin, fill = taxon_labels[,color]) +
        geom_rect(data=taxon_labels, aes(xmin = xmin-2, xmax = xmax+2), ymin = ymin - (ymax - ymin)*0.05, ymax = ymin, fill = 'grey') +
        #geom_hline(yintercept = ymin, color='grey', size=0.5) +
        geom_text(data=taxon_labels, aes(x=xc, label=as.character), y = ymin - plotsize*0.03, hjust=0.5, vjust=1, size=2.5, angle = rot) +
        geom_text(data=isotope_labels, aes(x=xc, label=label), y = ymin - plotsize*0.025, hjust=0.5, vjust=0, size=2.5, angle = rot) +
        #geom_text(data=outliers, aes(x=xval2, label=lo), y = ymin, hjust=0.5, vjust=0, size=2, angle = 90) +
        #geom_text(data=outliers, aes(x=xval2, label= hi), y = ymax, hjust=0.5, vjust=0, size=2, angle = 90) +
        geom_segment(data = isotope_labels, aes(x = xmin - 1, xend = xmax + 1), y = ymin - plotsize*0.005, yend = ymin - plotsize*0.005) +
        #scale_y_continuous(breaks = c(-10, 0,10,50,100,200)/100, labels = paste0(c(-10, 0,10,50,100,200), '%'), limits = c(-100,200)/100) +
        scale_y_continuous(labels = scales::percent, limits=c(ymin - plotsize * 0.05, ymax + plotsize * 0.025), expand = c(0,0)) +
        
        # scale_fill_manual(name = 'substrate',
            # values = c('ghostwhite', 'darkorchid1', 'palegreen1', 'olivedrab3', 'seagreen3'),
            # labels = c('No Carbon', 'Glucose', 'Oxalate', 'Exudate', 'Litter')) +
        scale_fill_manual(name = 'soil',
            values = c('burlywood3', 'darkslategrey', 'azure4'),
            labels = c('Granite', 'Basalt', 'Andesite')) +
        guides(fill = guide_legend(override.aes = list(alpha = 1, size= 2), label.position = "top")) +

        #geom_hline(yintercept = -100/100, color='gray', width=0.5) +
        #scale_shape_manual(values=c(21, 211S), breaks=c('C','O')) +
        facet_grid(timepoint ~ group, scales='free_x',space='free_x') +
        theme_bw() +
        labs(y = 'Atom Fraction Excess (AFE)') +
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
            axis.title.y = element_text(angle=90, vjust = 0, hjust = 0.5, size = 8, face = 'bold'),
            axis.text.y = element_text(angle=90, vjust = 0, hjust = 0.5, size = 8),
            panel.border = theme_border(type = c("right","left"), size=0.5),
            strip.background.x = element_rect(fill = 'grey', color='grey'),
            strip.background.y = element_rect(fill = 'grey', color='grey'),
            strip.text.x = element_text(angle = rot, size = 8),
            legend.text = element_text(angle = 90, vjust = 0.5, hjust = 0, size = 8),
            legend.title = element_text(angle = 90, size = 8, vjust = 1, hjust = 0.15, face='bold'),
            legend.direction = "horizontal",
            legend.spacing.y = unit(0, 'mm'),
            legend.spacing.x = unit(0, 'mm'),
            legend.key.height = unit(5, 'mm'),
            legend.key.width = unit(3, 'mm'),
            )
        
    
    pdf(paste0('qsip-soil-',RANK, '.pdf'), width = 20, height = 12)
        plot(pg)
    dev.off()
}))
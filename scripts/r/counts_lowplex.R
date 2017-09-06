library(ggplot2)
library(Biostrings)
library(dplyr)
library(tidyr)
library(viridis)
library(data.table)

plot.counts = function(counts, barcodePath, grouped = FALSE) {
  ds = data.frame()
  i = 0
  for (l in counts) {
    d = as.data.frame(fread(l,stringsAsFactors=FALSE))
    d = cbind(d,l)
    i = i + 1
    ds = rbind(ds,d)
  }
  colnames(ds) = c("IdxFirst", "IdxCombined", "Counts", "Run")
  barcodes = readDNAStringSet(barcodePath)
  bc_names = names(barcodes)
  ds$NameLeading = bc_names[ds$IdxFirst+1]
  ds$NameTrailing = bc_names[ds$IdxCombined+1]
  ds$BarcodePair = paste(bc_names[ds$IdxFirst+1],bc_names[ds$IdxCombined+1],sep="--")

  g = ggplot(data=ds, aes(x=BarcodePair, y=Counts, fill=Run)) +facet_wrap(~BarcodePair,scales = "free_x")+
    geom_bar(stat="identity", position=position_dodge(width=1))+
    geom_text(aes(label=Counts,y=mean(range(Counts))), color="black",
              position = position_dodge(1), size=3.5,angle = 90)+
    scale_color_brewer(palette = "Set1")+
    theme(legend.position="top", legend.direction = "vertical") +
    coord_cartesian(ylim = c(0, max(ds$Counts)*1.1))
  ggsave("counts_group.png",g,width=36,height=24,dpi=100,units="cm")

  g = ggplot(data=ds, aes(x=BarcodePair, y=Counts, fill=Run)) +
    geom_bar(stat="identity", position=position_dodge(width=1))+
    geom_text(aes(label=Counts), vjust=.4, hjust=-.1, color="black",
              position = position_dodge(0.9), size=3.5,angle = 90)+
    scale_color_brewer(palette = "Set1")+
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 0))+
    theme(legend.position="top", legend.direction = "vertical") +
    coord_cartesian(ylim = c(0, max(ds$Counts)*1.1))
  ggsave("counts_nogroup.png",g,width=36,height=24,dpi=100,units="cm")
}

plot.counts(c("m54007_170701_183412.subreadset.demux.counts",
              "m54007_170702_064558.subreadset.demux.counts",
              "m54200_170625_190247.subreadset.demux.counts",
              "m54200_170626_051342.subreadset.demux.counts"),
            "Sequel_RSII_16_barcodes_v1.fasta")

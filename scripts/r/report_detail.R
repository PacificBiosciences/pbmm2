#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

if (length(args) == 0)
  stop("Provide (1) *.lima.report file, (2) optional output type pdf or png [def: png], and (3+) optional individual barcode (pair) names.\n", call.=FALSE)

library(ggplot2,quietly = TRUE, warn.conflicts = FALSE)
library(dplyr,quietly = TRUE, warn.conflicts = FALSE)
library(tidyr,quietly = TRUE, warn.conflicts = FALSE)
library(viridis,quietly = TRUE, warn.conflicts = FALSE)
library(scales,quietly = TRUE, warn.conflicts = FALSE)
library(data.table,quietly = TRUE, warn.conflicts = FALSE)

reportPath = args[1]
output = "png"
if (length(args) >= 2)
  output = args[2]

report = as.data.frame(fread(reportPath,stringsAsFactors=FALSE))
report$IdxFirstNamed[report$IdxFirstNamed == "-1"] = "X"
report$IdxCombinedNamed[report$IdxCombinedNamed == "-1"] = "X"
report$IdxSecondCandidateNamed[report$IdxSecondCandidateNamed == "-1"] = "X"
report$IdxLowestNamed[report$IdxLowest == -1] = "X"
report$IdxHighestNamed[report$IdxHighest == -1] = "X"
report$BarcodePair = paste(report$IdxLowestNamed,report$IdxHighestNamed,sep="--")
report$Barcoded = report$IdxLowestNamed!="X" & report$IdxHighestNamed!="X"
if (length(args) >= 3) {
  barcodeNamesOfInterest = args[3:length(args)]
  combineBC = function(x) {
    s=unlist(strsplit(x,"--"))
    c(paste0(s[1],"--",s[2]),paste0(s[2],"--",s[1]))
  }
  barcodeNamesOfInterest = unlist(sapply(barcodeNamesOfInterest,function(x) {if(grepl("--",x)) { combineBC(x) } else {x }}))
  names(barcodeNamesOfInterest) = c()
  report = report %>% filter(IdxFirstNamed %in% barcodeNamesOfInterest | IdxCombinedNamed %in% barcodeNamesOfInterest | BarcodePair %in% barcodeNamesOfInterest)
}

reportCounts = count(report,BarcodePair)
count_labeller <- function(value){
  sapply(value,function(x) { paste(x,  " / ZMWs ", reportCounts[reportCounts$BarcodePair==x,]$n,sep="")})
}

tryCatch({report$ReadLengths = sapply(report$ReadLengths,function(x) list(as.numeric(unlist(strsplit(x,",")))))},error=function(e){})
names(report$ReadLengths) = c()
report$HQLength = sapply(report$ReadLengths,sum)

unique_bps = report %>% filter(Barcoded) %>% filter(PassedFilters == 1) %>% distinct(BarcodePair)
reportFilteredADP = report %>% filter(Barcoded) %>% filter(BarcodePair %in% unique_bps$BarcodePair)
reportFilteredADP$NumAdapters = as.numeric(reportFilteredADP$NumAdapters)
if (any(reportFilteredADP$NumAdapters >= 2)) reportFilteredADP[reportFilteredADP$NumAdapters >= 2,]$NumAdapters = 2
reportFilteredADP = reportFilteredADP%>% mutate(Filter = PassedFilters) %>% mutate(Filter=ifelse(Filter==0,"NONE","PASS"))

reportFilteredADP_pass = report %>% filter(Barcoded, PassedFilters == 1) %>% filter(BarcodePair %in% unique_bps$BarcodePair)
reportFilteredADP_pass$NumAdapters = as.numeric(reportFilteredADP_pass$NumAdapters)
if (any(reportFilteredADP_pass$NumAdapters >= 2)) reportFilteredADP_pass[reportFilteredADP_pass$NumAdapters >= 2,]$NumAdapters = 2

reportFiltered = report %>% filter(Barcoded) %>% filter(BarcodePair %in% unique_bps$BarcodePair) %>% mutate(Filter = "NONE", ScoreLead = ifelse(ScoreLead==-1,NA,ScoreLead))
reportFiltered_pass = reportFiltered %>% filter(PassedFilters == 1) %>% mutate(Filter = "PASS")

reportFiltered$ScoreCombinedAll = reportFiltered$ScoreCombined
if (any(reportFiltered$PassedFilters == 0)) reportFiltered[reportFiltered$PassedFilters == 0,]$ScoreCombined = NA

numadapters = report %>% filter(BarcodePair %in% unique_bps$BarcodePair) %>% select(BarcodePair, NumAdapters,Barcoded) %>% group_by(BarcodePair, NumAdapters, Barcoded)

baseYield = report %>% filter(Barcoded) %>% filter(BarcodePair %in% unique_bps$BarcodePair) %>% select(BarcodePair, ReadLengths, PassedFilters) %>% unnest(ReadLengths) %>% group_by(BarcodePair,PassedFilters) %>% mutate(MegaBases = sum(ReadLengths) / 1000000) %>% select(BarcodePair, MegaBases, PassedFilters) %>% distinct(BarcodePair,MegaBases, PassedFilters) %>% mutate(Filter = PassedFilters) %>% mutate(Filter=ifelse(Filter==0,"NONE","PASS"))
readYield = report %>% filter(Barcoded) %>% filter(BarcodePair %in% unique_bps$BarcodePair) %>% select(BarcodePair, ReadLengths, PassedFilters) %>% unnest(ReadLengths) %>% group_by(BarcodePair,PassedFilters) %>% count(BarcodePair, PassedFilters) %>% rename(NumReads = n) %>% mutate(Filter = PassedFilters) %>% mutate(Filter=ifelse(Filter==0,"NONE","PASS"))
zmwYield = report %>% filter(Barcoded) %>% filter(BarcodePair %in% unique_bps$BarcodePair) %>% select(BarcodePair, ZMW, PassedFilters) %>% count(BarcodePair, PassedFilters) %>% rename(NumZMWs = n, Filter = PassedFilters) %>% mutate(Filter=ifelse(Filter==0,"NONE","PASS"))

readLengthsUnnestedByBC = report %>% filter(Barcoded) %>% filter(BarcodePair %in% unique_bps$BarcodePair) %>% select(BarcodePair, ReadLengths, ScoreCombined) %>% unnest(ReadLengths) %>% mutate(ReadLengths = ReadLengths / 1000)
readLengthsUnnestedByBCZmw = report %>% filter(Barcoded) %>% filter(BarcodePair %in% unique_bps$BarcodePair) %>% select(BarcodePair, ZMW, ReadLengths, ScoreCombined) %>% unnest(ReadLengths) %>% group_by(BarcodePair,ZMW) %>% mutate(KiloBases = sum(ReadLengths) / 1000) %>% select(BarcodePair, KiloBases,ScoreCombined, ZMW) %>% distinct(BarcodePair,KiloBases,ScoreCombined, ZMW)

barcodeCounts = report %>% filter(Barcoded) %>% filter(BarcodePair %in% unique_bps$BarcodePair) %>% count(BarcodePair)
titration = report %>% filter(Barcoded) %>% filter(BarcodePair %in% unique_bps$BarcodePair) %>% select(BarcodePair, ScoreCombined) %>% group_by(BarcodePair) %>% arrange(BarcodePair,desc(ScoreCombined)) %>% count(BarcodePair,ScoreCombined) %>% mutate(cs = cumsum(n))
titration$Filter = "NONE"
titration_pass = report %>% filter(Barcoded, PassedFilters) %>% filter(BarcodePair %in% unique_bps$BarcodePair) %>% select(BarcodePair, ScoreCombined) %>% group_by(BarcodePair) %>% arrange(BarcodePair,desc(ScoreCombined)) %>% count(BarcodePair,ScoreCombined) %>% mutate(cs = cumsum(n))
titration_pass$Filter = "PASS"

readLengthsUnnested = report %>% select(ReadLengths, Barcoded) %>% unnest(ReadLengths)

dpi = 150
facetHeight = max(nrow(unique_bps)+1,4)/4*5+1
yieldHeight = max(nrow(unique_bps)+1,4)*0.5+3
facetWidth = 5 + min(nrow(unique_bps)+1,4)*5

g = ggplot(bind_rows(titration,titration_pass)) +
  facet_wrap(~BarcodePair, scales = "free_y", ncol = 4)+
  geom_line(aes(x=ScoreCombined,y=cs,color=Filter))+
  ylab("ZMW yield")+xlab("Mean Score")+
  coord_cartesian(xlim=c(0,100))+
  theme_light()
ggsave(paste0("detail_score_vs_yield.",output),g,width=facetWidth,height=facetHeight,dpi=dpi,units="cm",limitsize = FALSE)

g = ggplot(bind_rows(reportFiltered, reportFiltered_pass), aes(group = BarcodePair)) +
  facet_wrap(~BarcodePair, scales = "free_y", ncol = 4)+
  geom_freqpoly(binwidth=5, aes(x = ScoreLead, group=Filter, color=Filter),alpha=.5)+
  theme_light()+
  ylab("Number of ZMWs")+xlab("Score Leads")
ggsave(paste0("detail_score_lead.",output),g,width=facetWidth,height=facetHeight,dpi=dpi,units="cm",limitsize = FALSE)

g = ggplot(bind_rows(reportFiltered, reportFiltered_pass), aes(group = BarcodePair)) +
  facet_wrap(~BarcodePair, scales = "free_y", ncol = 4)+
  geom_freqpoly(binwidth=1, aes(x = SignalIncrease, group=Filter, color=Filter),alpha=.75)+
  theme_light()+
  ylab("Number of ZMWs")+xlab("Score Increase")
ggsave(paste0("detail_signal_increase.",output),g,width=facetWidth,height=facetHeight,dpi=dpi,units="cm",limitsize = FALSE)

g = ggplot(readLengthsUnnestedByBCZmw)+
  facet_wrap(~BarcodePair, labeller=as_labeller(count_labeller),ncol = 4)+
  geom_hex(aes(KiloBases, ScoreCombined, color = ..count..))+
  scale_fill_viridis()+
  scale_color_viridis()+
  coord_cartesian(xlim = c(0, quantile(readLengthsUnnestedByBCZmw$KiloBases,0.999)))+
  theme_light()+xlab("HQ Length in Kilo Bases")+ylab("Mean Score")
ggsave(paste0("detail_hq_length_vs_score.",output),g,width=facetWidth,height=facetHeight,dpi=dpi,units="cm",limitsize = FALSE)

g = ggplot(readLengthsUnnestedByBC)+
  facet_wrap(~BarcodePair, labeller=as_labeller(count_labeller),ncol = 4)+
  geom_hex(aes(ReadLengths, ScoreCombined, color = ..count..))+
  scale_fill_viridis()+
  scale_color_viridis()+
  coord_cartesian(xlim = c(0, quantile(readLengthsUnnestedByBC$ReadLengths,0.999)))+
  theme_light()+xlab("Read Length in Kilo Bases")+ylab("Mean Score")
ggsave(paste0("detail_read_length_vs_score.",output),g,width=facetWidth,height=facetHeight,dpi=dpi,units="cm",limitsize = FALSE)

g = ggplot(zmwYield) +
  geom_bar(aes(BarcodePair, NumZMWs, fill=Filter), stat='identity', width = .5)+
  scale_y_continuous(labels=comma)+ coord_flip()+
  theme_minimal() + theme(axis.text.x = element_text(hjust=1)) + ylab("Number of ZMWs")+
  ggtitle(paste("CV with filter: NONE=",(zmwYield %>% group_by(BarcodePair) %>% mutate(n=sum(NumZMWs)) %>% distinct(BarcodePair,n) %>% ungroup() %>% summarize(cv=round(100*sd(n)/mean(n))))$cv," PASS=",(zmwYield %>% filter(Filter=="PASS") %>% ungroup() %>% summarize(cv=round(100*sd(NumZMWs)/mean(NumZMWs))))$cv,sep=""))+
  theme(plot.title = element_text(hjust = 0.5))
ggsave(paste0("detail_yield_zmw.",output),g,width=25,height=yieldHeight,dpi=dpi,units="cm",limitsize = FALSE)

g = ggplot(readYield) +
  geom_bar(aes(BarcodePair, NumReads, fill=Filter), stat='identity', width = .5)+
  scale_y_continuous(labels=comma)+ coord_flip()+
  theme_minimal() + theme(axis.text.x = element_text(hjust=1)) + ylab("Number of Reads")+
  ggtitle(paste("CV with filter: NONE=",(readYield %>% group_by(BarcodePair) %>% mutate(n=sum(NumReads)) %>% distinct(BarcodePair,n) %>% ungroup() %>% summarize(cv=round(100*sd(n)/mean(n))))$cv," PASS=",(readYield %>% filter(Filter=="PASS") %>% ungroup() %>% summarize(cv=round(100*sd(NumReads)/mean(NumReads))))$cv,sep=""))+
  theme(plot.title = element_text(hjust = 0.5))
ggsave(paste0("detail_yield_read.",output),g,width=25,height=yieldHeight,dpi=dpi,units="cm",limitsize = FALSE)

g = ggplot(baseYield) +
  geom_bar(aes(BarcodePair, MegaBases, fill=Filter), stat='identity', width = .5)+
  scale_y_continuous(labels=comma)+
  theme_minimal() +coord_flip() + theme(axis.text.x = element_text(hjust=1)) + ylab("Yield in Mega Bases")+
  ggtitle(paste("CV with filter: NONE=",(baseYield %>% group_by(BarcodePair) %>% mutate(n=sum(MegaBases)) %>% distinct(BarcodePair,n) %>% ungroup() %>% summarize(cv=round(100*sd(n)/mean(n))))$cv," PASS=",(baseYield %>% filter(Filter=="PASS") %>% ungroup() %>% summarize(cv=round(100*sd(MegaBases)/mean(MegaBases))))$cv,sep=""))+
  theme(plot.title = element_text(hjust = 0.5))
ggsave(paste0("detail_yield_base.",output),g,width=25,height=yieldHeight,dpi=dpi,units="cm",limitsize = FALSE)

g = ggplot(report, aes(group = Barcoded, color = Barcoded, fill = Barcoded)) +
  geom_histogram(binwidth=2000,aes(HQLength),position = "identity", alpha=0.3)+
  coord_cartesian(xlim = c(0, quantile(report$HQLength,0.999))) +
  theme_minimal()+scale_x_continuous(labels=comma) + xlab("HQ Length") + ylab("Number of ZMWs")
ggsave(paste0("detail_hq_length_hist_barcoded_or_not.",output),g,width=25,height=15,dpi=dpi,units="cm",limitsize = FALSE)

g = ggplot(readLengthsUnnested, aes(group = Barcoded, color = Barcoded, fill = Barcoded)) +
  geom_histogram(binwidth=1000,aes(ReadLengths),position = "identity", alpha=0.3)+
  coord_cartesian(xlim = c(0, quantile(readLengthsUnnested$ReadLengths,0.999))) +
  theme_minimal() + xlab("Read Length") + ylab("Number of ZMWs")
ggsave(paste0("detail_read_length_hist_barcoded_or_not.",output),g,width=25,height=15,dpi=dpi,units="cm",limitsize = FALSE)

g = ggplot(reportFilteredADP, geom='none', aes(group = BarcodePair)) +
  facet_wrap(~BarcodePair,scales = "free_y",ncol=4)+
  coord_cartesian(xlim=c(0,100))+
  geom_histogram(data=reportFilteredADP,binwidth=5,aes(ScoreCombined,group=Filter,fill=Filter),color="grey",alpha=.15)+
  geom_freqpoly(binwidth=5,aes(ScoreCombined,color=as.factor(NumAdapters),group=as.factor(NumAdapters)))+
  theme_light() +scale_color_discrete(breaks = names(table(reportFilteredADP$NumAdapters)), name="#Adapters",
                                      labels=paste(c(replicate(length(names(table(reportFilteredADP$NumAdapters)))-1, "=="),">="),
                                                   c(as.numeric(names(table(reportFilteredADP$NumAdapters))[1:length(names(table(reportFilteredADP$NumAdapters)))]))))+
  ylab("Number of ZMWs") + xlab("Mean Score")
ggsave(paste0("detail_scores_per_adapter.",output),g,width=facetWidth,height=facetHeight,dpi=dpi,units="cm",limitsize = FALSE)

q = quantile(numadapters$NumAdapters,0.999)
g = ggplot(numadapters, aes(group = BarcodePair, x = NumAdapters)) +
  facet_wrap(~BarcodePair, scales = "free_y",ncol=4)+
  coord_cartesian(xlim = c(0, q)) +
  geom_histogram(binwidth=.5, aes(group=Barcoded))+
  theme_light()+
  scale_x_continuous(breaks=seq(0,q,round(q/min(10,q))))+
  ylab("Number of ZMWs")+xlab("Number of Adapters")
ggsave(paste0("detail_num_adapters.",output),g,width=facetWidth,height=facetHeight,dpi=dpi,units="cm",limitsize = FALSE)

x = reportFiltered %>% group_by(BarcodePair) %>% select(ReadLengths,BarcodePair) %>% unnest(ReadLengths)
g = ggplot(x, aes(group = BarcodePair, x = ReadLengths)) +
  facet_wrap(~BarcodePair,ncol=4)+
  coord_cartesian(xlim = c(0, quantile(x$ReadLengths,0.999))) +
  geom_histogram(binwidth = 1000, position="dodge", color="black", fill="gray")+
  theme_light()+ylab("Number of ZMWs")+xlab("Read Length")
ggsave(paste0("detail_read_length_hist_group_same_y.",output),g,width=facetWidth,height=facetHeight,dpi=dpi,units="cm",limitsize = FALSE)

g = ggplot(x, aes(group = BarcodePair, x = ReadLengths)) +
  facet_wrap(~BarcodePair, scales = "free_y",ncol=4)+
  coord_cartesian(xlim = c(0, quantile(x$ReadLengths,0.999))) +
  geom_histogram(binwidth = 1000, position="dodge", color="black", fill="gray")+
  theme_light()+ylab("Number of ZMWs")+xlab("Read Length")
ggsave(paste0("detail_read_length_hist_group_free_y.",output),g,width=facetWidth,height=facetHeight,dpi=dpi,units="cm",limitsize = FALSE)

g = ggplot(x, aes(group = BarcodePair, ReadLengths, color = BarcodePair)) +
  coord_cartesian(xlim = c(0, quantile(x$ReadLengths,0.999))) +
  geom_freqpoly(binwidth = 1000)+
  theme_minimal()+ylab("Number of ZMWs")+xlab("Read Length")
ggsave(paste0("detail_read_length_linehist_nogroup.",output),g,width=25,height=15,dpi=dpi,units="cm",limitsize = FALSE)

y = reportFiltered %>% group_by(BarcodePair) %>% select(HQLength,BarcodePair) %>% unnest(HQLength)
g = ggplot(y, aes(group = BarcodePair, x = HQLength)) +
  facet_wrap(~BarcodePair,ncol=4)+
  coord_cartesian(xlim = c(0, quantile(y$HQLength,0.999))) +
  geom_histogram(binwidth = 2000, position="dodge", color="black", fill="gray")+
  theme_light()+ylab("Number of ZMWs")+xlab("HQ Length")
ggsave(paste0("detail_hq_length_hist_group_same_y.",output),g,width=facetWidth,height=facetHeight,dpi=dpi,units="cm",limitsize = FALSE)

g = ggplot(y, aes(group = BarcodePair, x = HQLength)) +
  facet_wrap(~BarcodePair, scales = "free_y",ncol=4)+
  coord_cartesian(xlim = c(0, quantile(y$HQLength,0.999))) +
  geom_histogram(binwidth = 2000, position="dodge", color="black", fill="gray")+
  theme_light()+ylab("Number of ZMWs")+xlab("HQ Length")
ggsave(paste0("detail_hq_length_hist_group_free_y.",output),g,width=facetWidth,height=facetHeight,dpi=dpi,units="cm",limitsize = FALSE)

g = ggplot(y, aes(group = BarcodePair, HQLength, color = BarcodePair)) +
  coord_cartesian(xlim = c(0, quantile(y$HQLength,0.999))) +
  geom_freqpoly(binwidth = 2000)+
  theme_minimal()+ylab("Number of ZMWs")+xlab("HQ Length")
ggsave(paste0("detail_hq_length_linehist_nogroup.",output),g,width=25,height=15,dpi=dpi,units="cm",limitsize = FALSE)
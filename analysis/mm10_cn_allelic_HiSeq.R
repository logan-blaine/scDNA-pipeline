library(tidyverse)
library(rtracklayer)
library(furrr)
library(readxl)


plan(multisession(workers = 8))
# chromosomes=paste0("chr",c(1:19,"X"))
chromosomes=c(1:19,"X")


small_bin_size = 1e6
# bin_size=250e3
# bin_size_name="250kb"

bin_size=5e6
bin_size_name="5Mb"

# bin_size=1e6
# bin_size_name="1Mb"

project_dir = "/pellmanlab/logan/data/mouse_HiSeq"
plot_dir = sprintf("%s/plots/allelic_CN_%s", project_dir, bin_size_name)

dir.create(plot_dir, showWarnings=F, recursive=T)
# dir.create(file.path(plot_dir, "zoom"), showWarnings=F, recursive=T)


fam_files = list.files(file.path(project_dir, "read_depth"), full.names = T) %>%
  `names<-`(basename(.) %>% gsub("\\..*","",.)) 

bin_counts = function(x) {
  read_tsv(x,comment="@") %>% 
    dplyr::rename_all(tolower) %>% 
    mutate(bin=floor(start/small_bin_size)*small_bin_size) %>%
    group_by(contig, bin) %>%
    summarize(n=n(), dp=sum(count),.groups="drop") %>%
    mutate(contig = gsub("chr","",contig)) %>% 
    filter(contig %in% chromosomes) %>% 
    mutate(contig=factor(contig, levels=chromosomes))
}

variants=read_tsv("/pellmanlab/logan/data/references/mm10/129S_phasing/129S1.variants.filtered.table")

snp_bin_size=1e4
snp_counts=variants %>% mutate(bin=floor(position/snp_bin_size)*snp_bin_size) %>% 
  group_by(contig,bin) %>% summarize(n=n())
qplot(snp_counts$n)

binned_cts=future_map_dfr(fam_files, bin_counts, .id="sample", .progress=T) %>% 
   mutate(sample=as.factor(sample))# may take a few minutes

stam_tbl=read_xls("/pellmanlab/logan/data/mouse_embryo_MN/all_known_relationship_table.xls")
  

summary_tbl=binned_cts %>% group_by(sample, contig) %>%
  summarize(bases=max(bin), dp=sum(dp)) %>%
  summarize(fold_coverage=sum(dp)*150/sum(bases)) %>%
  mutate(sample=gsub("_L.+$","",sample)) %>% 
  # mutate(group=gsub("[A-Z]$","",sample)) %>% 
  inner_join(stam_tbl)
  # mutate(group=gsub("_[0-9]","",embryo_id)) %>% 
  # mutate(group_index = gsub("CE..","",group) %>% as.numeric) %>% 
  # mutate(group=paste0(toupper(substring(group,3,3)),group_index))
guides=unique(summary_tbl$guide)

fractional_coverage=binned_cts %>%
  mutate(sample=gsub("_L.+$","",sample)) %>% 
  inner_join(summary_tbl, by="sample") %>% 
  mutate(sample=embryo_id) %>% 
  mutate(contig=factor(contig, levels=chromosomes)) %>% 
  group_by(sample) %>%
  filter(median(dp)>0) %>% 
  mutate(cn=dp/median(dp)) %>% 
  group_by(contig, bin) %>%
  filter(mean(cn>0.25)>0.2) %>% # i.e. at least 20% of samples have appreciable coverage here
  ungroup()
  
raw_cn_estimates = fractional_coverage %>% 
  mutate(bin=floor(bin/bin_size)*bin_size) %>% 
  group_by(group, sample, contig, bin) %>%
  summarize(cn=mean(cn,na.rm=T)) %>% 
  group_by(sample) %>%
  mutate(cn=cn/median(cn)) 
scale_colors=c("#1768AF", "white", "#A51E1E", "orange", "#FFCC33")

g=raw_cn_estimates %>% 
  ggplot()+
  geom_tile(aes(x=bin,y=sample,fill=(2*cn))) +
  scale_x_continuous(expand=c(0,0))+
  # facet_grid(allele~chrom,space="free_x",scale="free_x",switch="y")+
  # facet_grid(group~contig,space="free",scale="free")+
  theme_void()+
  theme(panel.spacing = unit(0, "lines"),
        # theme(strip.text.x = element_text(size=3),  panel.spacing.x = unit(0, "lines"),
        panel.border = element_rect(size=0.1), legend.position = "bottom") +
  # scale_fill_discrete(colors=scale_colors, values=seq(0,8),name="CN")
  scale_fill_gradientn(colors=scale_colors, limits = c(0,8), name="CN")+
  facet_grid(group~contig,space="free",scale="free")

ggsave(file.path(plot_dir,"unnormalized.pdf"),g, width=7.5,height=10)

  # filter(mean(cn)>0) %>%
  # group_by(sample) %>%
  # mutate(cn=dp/median(dp))
  # group_by(sample, contig) %>%  # TEST
  # mutate(cn = ksmooth(bin, cn, "box", bin_size, x.points=bin)$y) %>% # TEST
  # ungroup()

# fractional_coverage %>% 
#   anti_join(bad_samples) %>% 
#   # filter(contig=="chrX") %>% 
#   group_by(contig, bin) %>% 
#   summarize(m=mean(cn>0.1)) %>% 
#   pull(m) %>% qplot()

# fractional_coverage %>% filter(contig=="chrX") %>% sample_n(1000) %>% 
#   ggplot() + geom_point(aes(bin, log2(cn)),alpha=0.1)

outlier_chromosomes = fractional_coverage %>% 
  group_by(sample, contig) %>% 
  summarize(m=median(cn)) %>% 
  filter(!between(m, 0.66, 1.33))
outlier_chromosomes %>% ggplot(aes(m)) + geom_histogram()+ facet_wrap(~contig)

# outlier_chromosomes = fractional_coverage %>%
#   group_by(contig, sample) %>%
#   summarize(m=median(cn)) %>%
#   mutate(m=m/median(m)) %>%
#   filter(!between(m, 0.75, 1.25))

# ggplot(aes(2*m)) + geom_histogram() + xlim(0,4) + facet_wrap(~contig)
# outlier_chromosomes
# ggplot() + geom_col(aes(contig, m)) + facet_wrap(~sample) +
# theme_void() + theme(strip.background = element_blank(), strip.text = element_blank())

# max_fold_adjust=2.5
median_cts=fractional_coverage %>%
  # anti_join(bad_samples, by="sample") %>%
  anti_join(outlier_chromosomes, by=c("contig", "sample")) %>%
  group_by(contig, bin) %>%
  summarize(med_cn=median(cn)) %>% 
  ungroup() 
  # filter(med_cn %>% between(1/max_fold_adjust, max_fold_adjust))
median_cts$med_cn %>% qplot()

median_cts %>% ungroup() %>% ggplot() + geom_histogram(aes(med_cn))

normalized_cts=inner_join(fractional_coverage, median_cts) %>% 
  mutate(cn=cn/med_cn)

# bin_size=10e6
binned_norm_cts=normalized_cts %>%
  mutate(bin=floor(bin/bin_size)*bin_size) %>% 
  group_by(group, sample, contig, bin) %>%
  summarize(cn=mean(cn,na.rm=T)) %>% 
  group_by(sample) %>%
  mutate(cn=cn/median(cn)) %>% 
  # mutate(ploidy=ifelse(grepl("EM_SC_190809_M1",sample),3,2),cn=cn*ploidy) %>%
  ungroup()

# final_frame_T=binned_norm_cts %>% 
#   transmute(sample,chrom_bin=paste(contig,bin/bin_size,sep="_"),cn) %>% 
#   spread(chrom_bin,cn) %>%
#   column_to_rownames("sample")
# 
# scale_colors=c("#1768AF", "white", "#A51E1E", "orange", "#FFCC33")
# 
# set.seed(12345)
# hc = final_frame_T %>% 
#   dist() %>% hclust("ward.D2")

# final_df=binned_norm_cts %>%
#   mutate(sample=factor(sample, levels=hc$labels[hc$order])) %>% 
#   mutate(contig=gsub("chr","",contig) %>% factor(levels=c(1:22,"X")))
  # mutate(group=gsub("_[^_]+$","",sample))

g=binned_norm_cts %>%
  # mutate(ploidy=ifelse(grepl("9_M1",sample),3,2),cn=cn*ploidy) %>% 
  ggplot()+
  geom_tile(aes(x=bin,y=sample,fill=(2*cn))) +
  scale_x_continuous(expand=c(0,0))+
  # facet_grid(allele~chrom,space="free_x",scale="free_x",switch="y")+
  facet_grid(group~contig,space="free",scale="free")+
  theme_void()+
  theme(panel.spacing = unit(0, "lines"),
        # theme(strip.text.x = element_text(size=3),  panel.spacing.x = unit(0, "lines"),
        panel.border = element_rect(size=0.1), legend.position = "bottom") +
  # scale_fill_discrete(colors=scale_colors, values=seq(0,8),name="CN")
  scale_fill_gradientn(colors=scale_colors, limits = c(0,8), name="CN")
ggsave(file.path(plot_dir, "median_normalized.pdf"),g, width=7.5,height=10)

early_rep=import.bedGraph("/pellmanlab/common/mm10_repliseq/mESC.mm10.early.bedGraph.gz")
early_df=as_tibble(early_rep) %>% mutate(label="early")
late_rep=import.bedGraph("/pellmanlab/common/mm10_repliseq/mESC.mm10.late.bedGraph.gz")
late_df=as_tibble(late_rep) %>% mutate(label="late")

combined=rbind(early_df,late_df) %>% 
  mutate(bin=floor(start/bin_size)*bin_size) %>%
  group_by(seqnames, bin,label) %>%
  summarize(score=mean(score)) %>% 
  pivot_wider(names_from="label",values_from="score") %>%
  ungroup() %>% 
  mutate(contig=gsub("chr","",seqnames)) %>% 
  mutate(contig=factor(contig,levels=chromosomes),timing=rank(late-early)/n()) %>% 
  drop_na()

loess_cn = combined  %>% inner_join(binned_norm_cts) %>% 
  group_by(sample, contig) %>%
  mutate(cn_norm_chrom=cn/median(cn)) %>% 
  group_by(sample) %>% 
  group_split() %>% 
  map_dfr(function(x) mutate(x, cn_pred=loess(cn_norm_chrom ~ early+late, span=0.7) %>% predict)) %>%
  # map_dfr(function(x) mutate(x, cn_pred=loess(cn_norm_chrom ~ timing, span=0.7) %>% predict)) %>%
  mutate(cn=cn/cn_pred * 2)

g=loess_cn %>% 
  ggplot()+
  geom_tile(aes(x=bin,y=sample,fill=cn)) +
  scale_x_continuous(expand=c(0,0))+
  # facet_grid(allele~chrom,space="free_x",scale="free_x",switch="y")+w
  facet_grid(group~contig,space="free",scale="free")+
  theme_void()+
  theme(panel.spacing = unit(0, "lines"),
        # theme(strip.text.x = element_text(size=3),  panel.spacing.x = unit(0, "lines"),
        panel.border = element_rect(size=0.1), legend.position = "bottom") +
  # scale_fill_discrete(colors=scale_colors, values=seq(0,8),name="CN")
  scale_fill_gradientn(colors=scale_colors, limits = c(0,8), name="CN")
ggsave(file.path(plot_dir,"loess.pdf"),g, width=7.5,height=10)


allele_files = list.files(file.path(project_dir, "allelic_depth"), full.names = T) %>%
  `names<-`(basename(.) %>% gsub("\\..*","",.) %>% gsub("EM_SC_","",.) %>% gsub("M[0-9]+_","",.))

read_het_coverage=function(x) {
  snp_bin_size=5e4
  a=read_tsv(x, comment="@") %>%
    dplyr::rename_all(tolower) %>%
    dplyr::rename(ref=refallele, alt=altallele) %>%
    semi_join(variants, by=c("contig","position","ref", "alt")) %>%
    mutate(bin = floor(position/snp_bin_size)*snp_bin_size) %>%
    group_by(contig, bin) %>%
    # summarize(refcount=mean(refcount), altcount=mean(altcount)) %>%
    summarize(refcount=mean(refcount), altcount=mean(altcount)) %>%
    mutate(bin = floor(bin/bin_size)*bin_size) %>%
    group_by(contig, bin) %>%
    # summarize(A=mean(refcount), B=mean(altcount), f=A/(A+B)) %>% 
    # summarize(f=refcount/(refcount+altcount))
    summarize(A=mean(refcount),B=mean(altcount), f=A/(A+B),n=n())
    # summarize(A=mean(refcount),B=mean(altcount), f=A/(A+B),n=n())
}

binned_allelic_counts=future_map_dfr(allele_files, read_het_coverage, .id="sample", .progress=T)
binned_fracs=binned_allelic_counts %>%
  filter(contig!="chrX") %>% 
  mutate(sample=gsub("_L.+$","",sample)) %>%
  mutate(group=gsub("[A-Z]$","",sample))


# summarize(A=median(A),B=median(B)) %>%
# summarize(f=median(A/(A+B), na.rm=T))
ploidy = binned_fracs %>%
  group_by(sample) %>%
  summarize(ref_allele_frac=mean(f,na.rm = T)) %>% 
  mutate(phased_maf=pmin(1-ref_allele_frac,ref_allele_frac)) %>%
  mutate(ploidy=ifelse(1/phased_maf<3.5, round(1/phased_maf), 2))


cn_deviation=loess_cn %>% 
  filter(contig!="X") %>% 
  group_by(sample) %>% 
  summarize(deviation=mad(cn/2))
# qplot(cn_deviation$deviation)

cutoff=0.15
final_summary=summary_tbl %>% 
  full_join(ploidy,by="sample") %>% 
  replace_na(list(ploidy=2)) %>% 
  mutate(sample=embryo_id) %>% 
  full_join(cn_deviation) %>%
  drop_na(sample) %>% 
  group_by(group) %>% 
  mutate(pass_deviation_filter=!is.na(deviation) & deviation<cutoff) %>%  
  mutate(pass_ploidy_filter = ploidy==median(ploidy,na.rm=T)) %>% 
  mutate(pass_all_filters=pass_deviation_filter & pass_ploidy_filter) %>% 
  ungroup()
  # select(-sample) %>% 
  # inner_join(stam_tbl) %>%
  # arrange(group, sample) %>% 
  # mutate(sample=factor(sample, levels = unique(sample)))


final_summary %>% filter(deviation<cutoff) %>% 
  ggplot()+geom_col(aes(embryo_id, phased_maf, fill=group)) +
  geom_hline(aes(yintercept=(1/2.5))) +
  geom_hline(aes(yintercept=(1/3.5))) + theme_minimal()
                            

# ggplot(final_summary)+geom_point(aes(x=deviation, y=f, color=pass_filter))

# final_summary %>% filter(pass_filter & ploidy!=2)
# binned_fracs %>% filter(sample=="190809_9H") %>% ggplot()+geom_point(aes(x=bin,y=f))+facet_grid(sample~contig)

summary_export_tbl=final_summary %>% 
  select(-sample) %>% 
  inner_join(stam_tbl)
write_excel_csv(summary_export_tbl, file.path(project_dir,"summary.csv"))

plot_cn_final = loess_cn %>%
  inner_join(final_summary) %>% 
  filter(pass_all_filters) %>% 
  mutate(cn=cn*ploidy/2) %>% 
  mutate(group_index = gsub("^.","",group) %>% as.numeric) %>% 
  arrange(guide, group_index, desc(sample)) %>% 
  mutate(group=factor(group, levels=unique(group))) %>% 
  mutate(sample=factor(sample, levels=unique(sample)))

for(condition in guides) {
  filtered_plot_tbl=plot_cn_final %>%
    filter(guide==condition)
  n_samples = length(unique(filtered_plot_tbl$sample))
  g=ggplot(filtered_plot_tbl)+
      theme_void()+
      ggtitle(condition)+
      theme(panel.spacing = unit(0, "lines"),
            axis.text.y = element_text(size=4),
            # theme(strip.text.x = element_text(size=3),  panel.spacing.x = unit(0, "lines"),
            panel.border = element_rect(size=0.1), legend.position = "bottom") +
      geom_tile(aes(x=bin,y=sample,fill=cn)) +
      scale_x_continuous(expand=c(0,0))+
      # scale_y_discrete(position="right")+
      # facet_grid(.~contig,space="free",scale="free")+
      facet_grid(group~contig,space="free",scale="free",switch="x")+
      # scale_fill_discrete(colors=scale_colors, values=seq(0,8),name="CN")
      scale_fill_gradientn(colors=scale_colors, limits = c(0,8), name="CN")
  ggsave(file.path(plot_dir,paste0(condition,"_final.pdf")),g, width=7.5,height=n_samples/20+0.8)
  ggsave(file.path(plot_dir,paste0(condition,"_final.eps")),g, width=7.5,height=n_samples/20+0.8)
}



#######
# EXPERIMENTAL: SEGMENTATION ANALYSIS
#######

library(PSCBS)
chrX_ID="99"
cn_data=plot_cn_final %>%
  group_by(contig, bin) %>% 
  mutate(w=1/mad(cn)) %>% 
  group_by(sample, contig) %>% 
  transmute(chromosome=gsub("chr","",contig), chromosome=ifelse(chromosome=="X",chrX_ID,chromosome), x=bin, y=cn,w=w/max(w)) %>%
  # transmute(chromosome=gsub("chr","",contig), x=bin, y=cn) %>% 
  arrange(sample, chromosome, x)

cn_data %>% ggplot() + geom_point(aes(x, w)) + facet_wrap(~chromosome)

split_df = cn_data %>% group_by(sample) %>% group_split(.,.keep=F) %>% `names<-`(unique(cn_data$sample))
segmented_data=future_map_dfr(split_df,
                              ~segmentByCBS(y=.,seed="1", avg="mean",
                                            undo = 1, p.method="perm") %>% 
                                getSegments, .id="sample", .progress=T)
# knownSegments=findLargeGaps(cn_data, minLength = 500e3) %>% gapsToSegments()
plot_segments=segmented_data %>% filter(!is.na(chromosome)) %>% 
  mutate(chromosome=ifelse(chromosome==chrX_ID,"X",chromosome) %>% factor(levels=chromosomes)) %>% 
  inner_join(distinct(plot_cn_final,sample, group))
  # mutate(sample=factor(sample, levels=hc$labels[hc$order]))  %>% 2
  # mutate(mean=ifelse(nbrOfLoci<=10,ifelse(is.na(lead(nbrOfLoci)) | lag(nbrOfLoci)>lead(nbrOfLoci),lag(mean),lead(mean)),mean))

col_min=0
col_max=8
g=plot_segments %>% 
  ggplot()+
  geom_tile(aes(x=(start+end)/2,width=end-start,y=sample,fill=round(mean))) +
  scale_x_continuous(expand=c(0,0))+
  # facet_grid(allele~chrom,space="free_x",scale="free_x",switch="y")+
  facet_grid(group~chromosome,space="free",scale="free")+
  theme_void()+
  theme(panel.spacing = unit(0, "lines"),
        # theme(strip.text.x = element_text(size=3),  panel.spacing.x = unit(0, "lines"),
        panel.border = element_rect(size=0.1), legend.position = "bottom") +
  scale_fill_gradientn(colors=scale_colors, limits = c(col_min,col_max), breaks=col_min:col_max, name="CN")
ggsave(file.path(plot_dir,"segmented.pdf"),g, width=7.5,height=10)

# 
#   chr_max=normalized_allelic_cn %>% group_by(contig) %>% summarize(xmax=max(bin))
#   rects=rbind(mutate(chr_max, start=0.5, end=1.5), mutate(chr_max, start=2.5, end=3.5), mutate(chr_max, start=4.5, end=5.5))
#   ymax_plot=6
#   
#   group_svs = filtered_svs %>% filter(grepl(sample, group)) 
#     
#   intrachromosomal=group_svs %>% filter(chr1==chr2) %>%
#     mutate(contig=sub("chr","",chr1) %>% factor(levels=c(1:22,"X","Y")))
#   interchromosomal=group_svs %>% filter(chr1!=chr2)
#   
#   breakends=rbind(mutate(interchromosomal, contig=chr1, pos=pos1),
#                   mutate(interchromosomal, contig=chr2, pos=pos2)) %>% 
#     mutate(contig=sub("chr","",contig) %>% factor(levels=c(1:22,"X","Y")))
#   
#   min_n = ceiling(bin_size/small_bin_size  * .15)
#   
#   g=normalized_allelic_cn %>%   
#     # filter(n>=min_n) %>%
#     mutate(copy_number=pmin(ymax_plot, copy_number)) %>%
#     ggplot()+ 
#     ggtitle(group)+
#     geom_rect(data=rects,aes(xmin=0,xmax=xmax, ymin=start, ymax=end), fill="#e6e7e8")+
#     geom_point(aes(bin, copy_number, color=allele, alpha=n),size=bin_size/25e6) +
#     scale_color_manual(values=list("A"="red","B"="blue","total"="black"))+
#     facet_grid(allele~contig, scales="free_x",space="free_x")+
#     scale_x_continuous(name=NULL,breaks=NULL,limits=c(0,NA), expand=c(0,0))+
#     scale_y_continuous(name="Haplotype Copy Number", limits=c(0,ymax_plot), breaks=0:5, minor_breaks=NULL)+
#     theme_bw() +
#     theme(panel.spacing = unit(0,"lines"), panel.grid=element_blank(), legend.position = "none")
#   if(nrow(interchromosomal)) {
#     g=g+geom_segment(data=breakends, aes(x=pos,xend=pos,y=5,yend=6), color="goldenrod1")
#   }
#   
#   if(nrow(intrachromosomal)) {
#     g=g+geom_curve(data=intrachromosomal, aes(x=pos1,xend=pos2,y=6,yend=6),color="limegreen")
#   }
#   # g
#   # g
#   
#   plot_file = sprintf("%s.pdf", group)
#   ggsave(sprintf(file.path(plot_dir, plot_file)),g,width=10,height=4,units="in")
#   
#   if (bin_size >= 1e6) {
#     return()
#   }
#   
#   for (plot_contig in unique(normalized_allelic_cn$contig)) {
#     plot_cn=normalized_allelic_cn %>% filter(contig==plot_contig)
#     plot_rects=rects %>% filter(contig==plot_contig)
#     plot_intrachromosomal = intrachromosomal %>% filter(contig==plot_contig)
#     plot_breakends = breakends %>% filter(contig==plot_contig)
#     
#     g=plot_cn %>%   
#       mutate(copy_number=pmin(ymax_plot, copy_number)) %>%
#       ggplot()+ 
#       ggtitle(group)+
#       geom_rect(data=plot_rects,aes(xmin=0,xmax=xmax, ymin=start, ymax=end), fill="#e6e7e8")+
#       geom_point(aes(bin, copy_number, color=allele, alpha=n),size=bin_size/25e6) +
#       scale_color_manual(values=list("A"="red","B"="blue","total"="black"))+
#       facet_grid(allele~contig, scales="free_x",space="free_x")+
#       scale_x_continuous(name="Genomic Location (Mb)",breaks=seq(0,1e9,1e7), labels=seq(0,1000,10), limits=c(0,NA), expand=c(0,0))+
#       scale_y_continuous(name="Haplotype Copy Number", limits=c(0,ymax_plot), breaks=0:5, minor_breaks=NULL)+
#       theme_bw() +
#       theme(panel.spacing = unit(0,"lines"), panel.grid=element_blank(), legend.position = "none")
#     if(nrow(plot_breakends)) {
#       g=g+geom_segment(data=plot_breakends, aes(x=pos,xend=pos,y=5,yend=6), color="goldenrod1")
#     }
#     if(nrow(plot_intrachromosomal)) {
#       g=g+geom_curve(data=plot_intrachromosomal, aes(x=pos1,xend=pos2,y=6,yend=6),color="limegreen")
#     }
#     plot_file = sprintf("%s.chr%s.pdf", group, plot_contig)
#     ggsave(file.path(plot_dir, "zoom", plot_file),g,width=10,height=4,units="in")
#   
#   }
#   
#   for (plot_contig in c(1:19,"X")) {
#     plot_cn=normalized_allelic_cn %>% filter(contig==plot_contig)
#     plot_rects=rects %>% filter(contig==plot_contig)
#     plot_intrachromosomal = intrachromosomal %>% filter(contig==plot_contig)
#     plot_breakends = breakends %>% filter(contig==plot_contig)
# 
#     g=plot_cn %>%   
#       mutate(copy_number=pmin(ymax_plot, copy_number)) %>%
#       ggplot()+ 
#       geom_rect(data=plot_rects,aes(xmin=0,xmax=xmax, ymin=start, ymax=end), fill="#e6e7e8")+
#       geom_point(aes(bin, copy_number, color=allele, alpha=n),size=bin_size/25e6) +
#       scale_color_manual(values=list("A"="red","B"="blue","total"="black"))+
#       facet_grid(allele~contig, scales="free_x",space="free_x")+
#       scale_x_continuous(name="Genomic Location (Mb)",breaks=seq(0,1e9,1e7), labels=seq(0,1000,10), limits=c(0,NA), expand=c(0,0))+
#       scale_y_continuous(name="Haplotype Copy Number", limits=c(0,ymax_plot), breaks=0:5, minor_breaks=NULL)+
#       theme_bw() +
#       theme(panel.spacing = unit(0,"lines"), panel.grid=element_blank(), legend.position = "none")
#     if(nrow(plot_intrachromosomal)) {
#       g=g+geom_curve(data=plot_intrachromosomal, aes(x=pos1,xend=pos2,y=6,yend=6),color="limegreen")
#     }
#     if(nrow(plot_breakends)) {
#       g=g+geom_segment(data=plot_breakends, aes(x=pos,xend=pos,y=5,yend=6), color="goldenrod1")
#     }
#     plot_file = sprintf("%s.chr%s.pdf", group, plot_contig)
#     ggsave(file.path(plot_dir, "zoom", plot_file),g,width=10,height=4,units="in")
#   }
# 
#   nrow(phased_ads)
# }, mc.cores=12)
# 
# # draw_arc = function(row) {
# #   num_points=30
# #   y1=4
# #   y2=5
# #   x1=row$pos1
# #   x2=row$pos2
# #   r = abs(x2-x1)/2
# #   m = (x1+x2)/2
# #   theta = seq(0, 1, length.out=num_points)
# #   x=m-r*cospi(theta)
# #   y=y1+ r*sinpi(theta)*(y2-y1)/r
# #   df=tibble(x,y,sample=row$sample, contig=row$chr1)
# #   geom_line(data=df, aes(x, y), lwd=0.5, color="#42AB5D")
# # }


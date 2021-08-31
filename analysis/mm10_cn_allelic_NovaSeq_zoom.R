library(tidyverse)
library(rtracklayer)
library(furrr)

plan(multisession(workers = 8))
# chromosomes=paste0("chr",c(1:19,"X"))
chromosomes=c(1:19,"X")

small_bin_size = 50e3
# bin_size=250e3
# bin_size_name="250kb"

bin_size=250e3
bin_size_name="250kb"

project_dir = "/pellmanlab/logan/data/mouse_embryo_MN"
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

snp_bin_size=5e4
snp_counts=variants %>% mutate(bin=floor(position/snp_bin_size)*snp_bin_size) %>% 
  group_by(contig,bin) %>% summarize(n=n())
# qplot(snp_counts$n)

binned_cts=future_map_dfr(fam_files, bin_counts, .id="sample", .progress=T) %>% 
   mutate(sample=gsub("EM_SC_","",sample)) %>% 
   mutate(sample=gsub("M[0-9]+_","",sample)) %>% 
   mutate(sample=as.factor(sample))# may take a few minutes

hiseq_summary=read_csv("/pellmanlab/logan/data/mouse_HiSeq/summary.csv") %>% 
  select(sample, embryo_id, guide, group, pass_all_filters)
  

summary_tbl=binned_cts %>% group_by(sample, contig) %>%
  summarize(bases=max(bin), dp=sum(dp)) %>%
  summarize(fold_coverage=sum(dp)*150/sum(bases)) %>%
  inner_join(hiseq_summary) 
guides=unique(summary_tbl$guide)

fractional_coverage=binned_cts %>%
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

outlier_chromosomes = fractional_coverage %>% 
  group_by(sample, contig) %>% 
  summarize(m=median(cn)) %>% 
  filter(!between(m, 0.66, 1.33))
outlier_chromosomes %>% ggplot(aes(m)) + geom_histogram()+ facet_wrap(~contig)

median_cts=fractional_coverage %>%
  # anti_join(bad_samples, by="sample") %>%
  anti_join(outlier_chromosomes, by=c("contig", "sample")) %>%
  group_by(contig, bin) %>%
  summarize(med_cn=median(cn)) %>% 
  ungroup() 
  # filter(med_cn %>% between(1/max_fold_adjust, max_fold_adjust))
median_cts$med_cn %>% qplot()

# median_cts %>% ungroup() %>% ggplot() + geom_histogram(aes(med_cn))

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
  future_map_dfr(function(x) mutate(x, cn_pred=loess(cn_norm_chrom ~ early+late, span=0.7) %>% predict)) %>%
  # map_dfr(function(x) mutate(x, cn_pred=loess(cn_norm_chrom ~ timing, span=0.7) %>% predict)) %>%
  mutate(cn=cn/cn_pred * 2)

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
  mutate(sample=gsub("EM_SC_","",sample)) %>% 
  mutate(sample=gsub("M[0-9]+_","",sample))
  # mutate(group=gsub("[A-Z]$","",sample))


# summarize(A=median(A),B=median(B)) %>%
# summarize(f=median(A/(A+B), na.rm=T))
ploidy = binned_fracs %>%
  filter(contig!="chrX") %>%
  group_by(sample) %>%
  summarize(f=mean(f,na.rm = T)) %>% 
  mutate(c=pmin(1-f,f)) %>%
  mutate(ploidy=ifelse(1/c<3.5, round(1/c), 2)) 


cn_deviation=loess_cn %>% 
  filter(contig!="X") %>% 
  group_by(sample) %>% 
  summarize(deviation=mad(cn/2))
# qplot(cn_deviation$deviation)

cutoff=0.3
final_summary=summary_tbl %>% 
  mutate(pass_all_filters=pass_all_filters & fold_coverage>=0.5) %>% 
  filter(!is.na(pass_all_filters))

combined_summary=read_csv("/pellmanlab/logan/data/mouse_HiSeq/summary.csv") %>% 
  left_join(final_summary, by=c("embryo_id","sample","group","guide"),suffix=c("_hiseq","_novaseq"))

# write_csv(combined_summary, file.path(project_dir, "combined_summary.csv"))

# ggplot(final_summary)+geom_point(aes(x=deviation, y=f, color=pass_all_filters))

# final_summary %>% filter(pass_all_filters & ploidy!=2)
# binned_fracs %>% filter(sample=="190809_9H") %>% ggplot()+geom_point(aes(x=bin,y=f))+facet_grid(sample~contig)

# write_excel_csv(final_summary, file.path(plot_dir,"summary.csv"))

plot_cn_final = loess_cn %>%
  dplyr::rename(embryo_id=sample) %>% 
  inner_join(combined_summary) %>% 
  filter(pass_all_filters_novaseq) %>% 
  mutate(cn=cn*ploidy/2) %>% 
  mutate(group_index = gsub("^.","",group) %>% as.numeric) %>% 
  arrange(guide, group_index, desc(sample)) %>% 
  mutate(group=factor(group, levels=unique(group))) %>% 
  mutate(sample=factor(embryo_id, levels=unique(embryo_id)))

# for(condition in guides) {
#   filtered_plot_tbl=plot_cn_final %>%
#     filter(guide==condition)
#   n_samples = length(unique(filtered_plot_tbl$sample))
#   g=ggplot(filtered_plot_tbl)+
#       ggtitle(condition)+
#       theme_void()+
#       theme(panel.spacing = unit(0, "lines"),
#             axis.text.y = element_text(size=4),
#             # theme(strip.text.x = element_text(size=3),  panel.spacing.x = unit(0, "lines"),
#             panel.border = element_rect(size=0.1), legend.position = "bottom") +
#       geom_tile(aes(x=bin,y=sample,fill=cn)) +
#       scale_x_continuous(expand=c(0,0))+
#       # facet_grid(.~contig,space="free",scale="free")+
#       facet_grid(group~contig,space="free",scale="free",switch="x")+
#       # scale_fill_discrete(colors=scale_colors, values=seq(0,8),name="CN")
#       scale_fill_gradientn(colors=scale_colors, limits = c(0,8), name="CN")
#   ggsave(file.path(plot_dir,paste0(condition,"_final.pdf")),g, width=7.5,height=n_samples/20+0.75)
# }
# 
# library(PSCBS)
# chrX_ID="99"
# cn_data=plot_cn_final %>%
#   group_by(contig, bin) %>% 
#   mutate(w=1/mad(cn)) %>% 
#   group_by(sample, contig) %>% 
#   transmute(chromosome=gsub("chr","",contig), chromosome=ifelse(chromosome=="X",chrX_ID,chromosome), x=bin, y=cn,w=w/max(w)) %>%
#   # transmute(chromosome=gsub("chr","",contig), x=bin, y=cn) %>% 
#   arrange(sample, chromosome, x)
# 
# cn_data %>% ggplot() + geom_point(aes(x, w)) + facet_wrap(~chromosome)
# 
# split_df = cn_data %>% group_by(sample) %>% group_split(.,.keep=F) %>% `names<-`(unique(cn_data$sample))
# segmented_data=future_map_dfr(split_df,
#                               ~segmentByCBS(y=.,seed="1", avg="mean",
#                                             undo = 1.75, min.width=5, p.method="perm") %>% 
#                                 getSegments, .id="sample", .progress=T)
# 
# # knownSegments=findLargeGaps(cn_data, minLength = 500e3) %>% gapsToSegments()
# plot_segments=segmented_data %>% filter(!is.na(chromosome)) %>% 
#   mutate(chromosome=ifelse(chromosome==chrX_ID,"X",chromosome) %>% factor(levels=chromosomes)) %>% 
#   inner_join(distinct(plot_cn_final,sample, group))
# 
# cnv_calls=plot_segments %>% 
#   select(-sampleName) %>% 
#   group_by(sample) %>% 
#   mutate(genome_cn=round(sum(nbrOfLoci*mean)/sum(nbrOfLoci))) %>% 
#   filter(chromosome=="X" | round(mean)!=genome_cn) %>% 
#   filter(nbrOfLoci>10) %>% 
#   group_by(group) %>% 
#   filter(chromosome!="X" | round(mean)!=median(round(mean))) %>% 
#   transmute(sample, group, contig=chromosome, start, end, copy_number=round(mean)) 
# write_excel_csv(cnv_calls,file.path(plot_dir,"cnv_calls_novaseq.csv"))
#   # mutate(sample=factor(sample, levels=hc$labels[hc$order]))  %>% 2
#   # mutate(mean=ifelse(nbrOfLoci<10,mean(mean,na.rm=T),mean))
#                      
#                      # ifelse(is.na(lag(nbrOfLoci)) | lag(nbrOfLoci)>lead(nbrOfLoci),
#                      #        ifelse(is.na(lag(mean)),lead(mean), lag(mean)),
#                      #        ifelse(is.na(lead(mean)),lag(mean), lead(mean))),
#                      # mean))
# plot_segments %>% 
#   mutate(end=ifelse(round(mean)==lead(round(mean)),lead(end),end)) %>% 
#   filter(round(mean)!=lag(round(mean)))
# 
# 
# 
# col_min=0
# col_max=8
# g=plot_segments %>% 
#   ggplot()+
#   geom_tile(aes(x=(start+end)/2,width=end-start,y=sample,fill=round(mean))) +
#   scale_x_continuous(expand=c(0,0))+
#   # facet_grid(allele~chrom,space="free_x",scale="free_x",switch="y")+
#   facet_grid(group~chromosome,space="free",scale="free")+
#   theme_void()+
#   theme(panel.spacing = unit(0, "lines"),
#         axis.text.y = element_text(size=4),
#         # theme(strip.text.x = element_text(size=3),  panel.spacing.x = unit(0, "lines"),
#         panel.border = element_rect(size=0.1), legend.position = "bottom") +
#       facet_grid(group~chromosome,space="free",scale="free",switch="x")+
#       # scale_fill_discrete(colors=scale_colors, values=seq(0,8),name="CN")
#     scale_fill_gradientn(colors=scale_colors, limits = c(0,8), name="CN")
# ggsave(file.path(plot_dir,"segmented.pdf"),g, width=7.5,height=5.75)
# ggsave(file.path(plot_dir,"segmented.eps"),g, width=7.5,height=5.75)

plot_cn_trim = plot_cn_final %>% select(sample, contig, bin, cn, group)

haplotype_cn=binned_fracs %>% 
  inner_join(combined_summary) %>% 
  filter(pass_all_filters_hiseq) %>% 
  mutate(sample=embryo_id) %>% 
  mutate(contig=gsub("chr","",contig) %>% factor(levels=c(1:19,"X"))) %>%
  inner_join(plot_cn_trim, by=c("sample","contig","bin","group")) %>%
  mutate(sample=embryo_id) %>% 
  mutate(A=cn*f, B=cn*(1-f)) %>%
  # inner_join(final_df) %>% 
  # mutate(ploidy=ifelse(group=="190809_1",3,2),cn=cn*ploidy) %>% 
  # mutate(A=cn*f, B=cn*(1-f)) %>% 
  pivot_longer(cols=c("A","B"),names_to="allele",values_to="hap_cn") %>% 
  # anti_join(uneven_samples) %>% 
  mutate(allele=factor(allele, levels=c("B","A"))) %>% 
  # mutate(sample_index = gsub("^[^_]+","",sample)) %>% 
  # mutate(embryo_id=paste0(group, sample_index)) %>% 
  # inner_join(fig2_tbl) %>% 
  arrange(sample, contig, desc(allele))
# 
# future_map(unique(haplotype_cn$group), function(plot_group) {
#   group_cn=haplotype_cn %>%
#     filter(allele=="total" | n>=3) %>%
#     filter(group==plot_group)
#   
#   ymax=ceiling(quantile(group_cn$hap_cn, .999, na.rm=T)+1)
#   ymax=max(group_cn$hap_cn, na.rm=T)+1
#   
#   g=ggplot(group_cn)+
#     geom_point(aes(x=bin,y=pmin(ymax,hap_cn),color=allele),size=0.3, shape=19)+
#     coord_cartesian(clip="off",default=T,expand = T) +
#     scale_x_continuous(breaks=NULL,name=NULL)+
#     scale_y_continuous(limits=c(0,ymax), breaks=seq(0,ymax,2), name="Haplotype Copy Number")+
#     facet_grid(embryo_id~contig,space="free_x",scale="free_x",switch="x")+
#     theme_bw()+
#     theme(axis.line.x = element_blank(),
#           # legend.position = "none",
#           panel.spacing.x = unit(0,"in"),
#           # panel.spacing.y =unit(0.15,"in"),
#           panel.border = element_rect(size=0.1, color="black", fill=NA),
#           strip.background = element_blank(), strip.placement = "outside")  +
#   # scale_color_manual(name="Haplotype", values=c("red","blue","black"), breaks=c("A","B","total"), labels=list(A="C57BL/6J",B="129Sv/Jae",total="Total"))
#     scale_color_manual(values=list(B="grey20",A="grey60"),
#                        name="Haplotype",
#                        breaks=c("B","A"), 
#                        labels=list(A="C57BL/6J",B="129Sv/Jae"))
#   ggsave(file.path(plot_dir,paste0(plot_group,".pdf")),g, width=10,height=1+0.5*length(unique(group_cn$sample)))
#   ggsave(file.path(plot_dir,paste0(plot_group,".eps")),g, width=10,height=1+0.5*length(unique(group_cn$sample)))
#   })

### OLD STUFF
cnv_calls = read_csv("/pellmanlab/logan/data/mouse_embryo_MN/plots/allelic_CN_1Mb/cnv_calls_novaseq.csv")
shared_cnvs = cnv_calls %>%
  group_by(group,contig,sample) %>% 
  summarize() %>% 
  summarize(n=n()) %>% 
  filter(n>1) %>% ungroup()

future_map(1:nrow(shared_cnvs), function(i) {
  plot_group=shared_cnvs$group[i] %>% as.character()
  plot_contig=shared_cnvs$contig[i] %>% as.character()
  group_cn=haplotype_cn %>%
    filter(n>=3) %>%
    filter(group==plot_group & contig==plot_contig)
  
  if (nrow(group_cn)==0) {
    print(paste(plot_group,plot_contig))
  }
  ymax=max(group_cn$hap_cn, na.rm = T)+1
  chr_max=group_cn %>% group_by(contig) %>% summarize(xmax=max(bin)/1e6)
  rects=rbind(mutate(chr_max, start=0.5, end=1.5), mutate(chr_max, start=2.5, end=3.5), mutate(chr_max, start=4.5, end=5.5))
  
  g=ggplot(group_cn)+
    # geom_rect(data=rects,aes(xmin=0,xmax=xmax, ymin=start, ymax=end), fill="#e6e7e8")+
    geom_point(aes(x=bin/1e6,y=pmin(ymax,hap_cn),color=allele),size=0.8, shape=19)+
    coord_cartesian(clip="off",default=T,expand = T) +
    scale_x_continuous(breaks=seq(0,1e4,10), minor_breaks=seq(0,1e4,5),name=sprintf("Position on chromosome %s (Mb)", plot_contig))+
    scale_y_continuous(limits=c(0,ymax), breaks=seq(0,ymax,2), minor_breaks=1:(ymax+1), name="Haplotype Copy Number")+
    facet_grid(embryo_id~.,space="free_x",scale="free_x",switch="x")+
    theme_bw()+
    theme(axis.line.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          # legend.position = "none",
          panel.spacing.x = unit(0,"in"), panel.spacing.y =unit(0.1,"in"),
          panel.border = element_rect(size=0.2, color="black", fill=NA),
          strip.background = element_blank(), strip.placement = "outside", legend.position = "bottom")  +
    # scale_color_manual(name="Haplotype", values=c("red","blue","black"), breaks=c("A","B","total"), labels=list(A="C57BL/6J",B="129Sv/Jae",total="Total"))
    scale_color_manual(values=list(B="grey20",A="grey60"),
                       name="Haplotype",
                       breaks=c("B","A"), 
                       labels=list(A="C57BL/6J",B="129Sv/Jae"))
  ggsave(file.path(plot_dir,sprintf("%s_chr%s.pdf",plot_group, plot_contig)),g,
         width=0.5+max(group_cn$bin/3e7),height=1+0.15*ymax*length(unique(group_cn$sample)))
  ggsave(file.path(plot_dir,sprintf("%s_chr%s.eps",plot_group, plot_contig)),g,
         width=0.5+max(group_cn$bin/3e7),height=1+0.15*ymax*length(unique(group_cn$sample)))
})

# normalized_cts %>% 
#   group_by(sample, contig) %>% 
#   summarize(cn=median(cn)) %>% 
#   ggplot(aes(contig, cn)) + geom_point()
# variants=read_tsv("/pellmanlab/logan/data/references/mm10/129S_phasing/129S1.variants.filtered.table")
# variants=read_tsv("/pellmanlab/logan/data/mouse_bulk_seq/129S1.variants.filtered.table")
groups = names(fam_files)

plots=mclapply(groups, function(group) {
  binned_fracs = list.files(file.path(project_dir, "allelic_depth"), pattern=group, full.names=T) %>%
    read_tsv(comment="@") %>%
    dplyr::rename_all(tolower) %>% 
    dplyr::rename(ref=refallele, alt=altallele) %>% 
    inner_join(variants, by=c("contig","position","ref", "alt")) %>%
    mutate(bin = floor(position/small_bin_size)*small_bin_size) %>%
    semi_join(median_cts, by=c("contig","bin")) %>% 
    group_by(contig, bin) %>%
    summarize(ref_count=mean(refcount),alt_count=mean(altcount)) %>% 
    mutate(ref_count = ksmooth(bin, ref_count, "box", bin_size,  x.points=bin)$y,
           alt_count = ksmooth(bin, alt_count, "box", bin_size,  x.points=bin)$y) %>%
    mutate(f = ref_count/(ref_count+alt_count)) 
  
  group_cts = normalized_cts %>% filter(grepl(group,sample)) %>% 
    # group_by(contig) %>% 
    # mutate(cn = ksmooth(bin, cn, "box", bin_size, x.points=bin)$y) %>%
    # ungroup() %>% 
    mutate(cn = cn/median(cn)) 

  ploidy = binned_fracs %>% drop_na(f) %>%  group_by(contig) %>%
    summarize(f=mean(f)) %>%
    summarize(c=median(pmin(1-f,f))) %>% 
    mutate(ploidy=ifelse(1/c<4.5, round(1/c), 1)) %>% pull(ploidy)
  
  normalized_allelic_cn = left_join(group_cts,binned_fracs, by=c("contig","bin")) %>% 
    select(sample, contig, bin, cn, f) %>% 
    mutate(total=cn*ploidy, A=f*total, B=(1-f)*total)  %>% 
    mutate(contig=sub("chr","",contig) %>% factor(levels=c(1:22,"X","Y"))) %>% 
    pivot_longer(cols=c("A","B","total"), names_to = "allele", values_to="copy_number", values_drop_na=T) %>% 
    mutate(bin=floor(bin/bin_size)*bin_size)  %>%
    group_by(contig, bin, allele) %>%
    summarize(copy_number=mean(copy_number, na.rm=T), n=n()) 
  
  chr_max=normalized_allelic_cn %>% group_by(contig) %>% summarize(xmax=max(bin))
  rects=rbind(mutate(chr_max, start=0.5, end=1.5), mutate(chr_max, start=2.5, end=3.5), mutate(chr_max, start=4.5, end=5.5))
  ymax_plot=6
  
  group_svs = filtered_svs %>% filter(grepl(sample, group)) 
    
  intrachromosomal=group_svs %>% filter(chr1==chr2) %>%
    mutate(contig=sub("chr","",chr1) %>% factor(levels=c(1:22,"X","Y")))
  interchromosomal=group_svs %>% filter(chr1!=chr2)
  
  breakends=rbind(mutate(interchromosomal, contig=chr1, pos=pos1),
                  mutate(interchromosomal, contig=chr2, pos=pos2)) %>% 
    mutate(contig=sub("chr","",contig) %>% factor(levels=c(1:22,"X","Y")))
  
  min_n = ceiling(bin_size/small_bin_size  * .15)
  
  g=normalized_allelic_cn %>%   
    # filter(n>=min_n) %>%
    mutate(copy_number=pmin(ymax_plot, copy_number)) %>%
    ggplot()+ 
    ggtitle(group)+
    geom_rect(data=rects,aes(xmin=0,xmax=xmax, ymin=start, ymax=end), fill="#e6e7e8")+
    geom_point(aes(bin, copy_number, color=allele, alpha=n),size=bin_size/25e6) +
    scale_color_manual(values=list("A"="red","B"="blue","total"="black"))+
    facet_grid(allele~contig, scales="free_x",space="free_x")+
    scale_x_continuous(name=NULL,breaks=NULL,limits=c(0,NA), expand=c(0,0))+
    scale_y_continuous(name="Haplotype Copy Number", limits=c(0,ymax_plot), breaks=0:5, minor_breaks=NULL)+
    theme_bw() +
    theme(panel.spacing = unit(0,"lines"), panel.grid=element_blank(), legend.position = "none")
  if(nrow(interchromosomal)) {
    g=g+geom_segment(data=breakends, aes(x=pos,xend=pos,y=5,yend=6), color="goldenrod1")
  }
  
  if(nrow(intrachromosomal)) {
    g=g+geom_curve(data=intrachromosomal, aes(x=pos1,xend=pos2,y=6,yend=6),color="limegreen")
  }
  # g
  # g
  
  plot_file = sprintf("%s.pdf", group)
  ggsave(sprintf(file.path(plot_dir, plot_file)),g,width=10,height=4,units="in")
  
  if (bin_size >= 1e6) {
    return()
  }
  
  for (plot_contig in unique(normalized_allelic_cn$contig)) {
    plot_cn=normalized_allelic_cn %>% filter(contig==plot_contig)
    plot_rects=rects %>% filter(contig==plot_contig)
    plot_intrachromosomal = intrachromosomal %>% filter(contig==plot_contig)
    plot_breakends = breakends %>% filter(contig==plot_contig)
    
    g=plot_cn %>%   
      mutate(copy_number=pmin(ymax_plot, copy_number)) %>%
      ggplot()+ 
      ggtitle(group)+
      geom_rect(data=plot_rects,aes(xmin=0,xmax=xmax, ymin=start, ymax=end), fill="#e6e7e8")+
      geom_point(aes(bin, copy_number, color=allele, alpha=n),size=bin_size/25e6) +
      scale_color_manual(values=list("A"="red","B"="blue","total"="black"))+
      facet_grid(allele~contig, scales="free_x",space="free_x")+
      scale_x_continuous(name="Genomic Location (Mb)",breaks=seq(0,1e9,1e7), labels=seq(0,1000,10), limits=c(0,NA), expand=c(0,0))+
      scale_y_continuous(name="Haplotype Copy Number", limits=c(0,ymax_plot), breaks=0:5, minor_breaks=NULL)+
      theme_bw() +
      theme(panel.spacing = unit(0,"lines"), panel.grid=element_blank(), legend.position = "none")
    if(nrow(plot_breakends)) {
      g=g+geom_segment(data=plot_breakends, aes(x=pos,xend=pos,y=5,yend=6), color="goldenrod1")
    }
    if(nrow(plot_intrachromosomal)) {
      g=g+geom_curve(data=plot_intrachromosomal, aes(x=pos1,xend=pos2,y=6,yend=6),color="limegreen")
    }
    plot_file = sprintf("%s.chr%s.pdf", group, plot_contig)
    ggsave(file.path(plot_dir, "zoom", plot_file),g,width=10,height=4,units="in")
  
  }
  
  for (plot_contig in c(1:19,"X")) {
    plot_cn=normalized_allelic_cn %>% filter(contig==plot_contig)
    plot_rects=rects %>% filter(contig==plot_contig)
    plot_intrachromosomal = intrachromosomal %>% filter(contig==plot_contig)
    plot_breakends = breakends %>% filter(contig==plot_contig)

    g=plot_cn %>%   
      mutate(copy_number=pmin(ymax_plot, copy_number)) %>%
      ggplot()+ 
      geom_rect(data=plot_rects,aes(xmin=0,xmax=xmax, ymin=start, ymax=end), fill="#e6e7e8")+
      geom_point(aes(bin, copy_number, color=allele, alpha=n),size=bin_size/25e6) +
      scale_color_manual(values=list("A"="red","B"="blue","total"="black"))+
      facet_grid(allele~contig, scales="free_x",space="free_x")+
      scale_x_continuous(name="Genomic Location (Mb)",breaks=seq(0,1e9,1e7), labels=seq(0,1000,10), limits=c(0,NA), expand=c(0,0))+
      scale_y_continuous(name="Haplotype Copy Number", limits=c(0,ymax_plot), breaks=0:5, minor_breaks=NULL)+
      theme_bw() +
      theme(panel.spacing = unit(0,"lines"), panel.grid=element_blank(), legend.position = "none")
    if(nrow(plot_intrachromosomal)) {
      g=g+geom_curve(data=plot_intrachromosomal, aes(x=pos1,xend=pos2,y=6,yend=6),color="limegreen")
    }
    if(nrow(plot_breakends)) {
      g=g+geom_segment(data=plot_breakends, aes(x=pos,xend=pos,y=5,yend=6), color="goldenrod1")
    }
    plot_file = sprintf("%s.chr%s.pdf", group, plot_contig)
    ggsave(file.path(plot_dir, "zoom", plot_file),g,width=10,height=4,units="in")
  }

  nrow(phased_ads)
}, mc.cores=12)

# draw_arc = function(row) {
#   num_points=30
#   y1=4
#   y2=5
#   x1=row$pos1
#   x2=row$pos2
#   r = abs(x2-x1)/2
#   m = (x1+x2)/2
#   theta = seq(0, 1, length.out=num_points)
#   x=m-r*cospi(theta)
#   y=y1+ r*sinpi(theta)*(y2-y1)/r
#   df=tibble(x,y,sample=row$sample, contig=row$chr1)
#   geom_line(data=df, aes(x, y), lwd=0.5, color="#42AB5D")
# }


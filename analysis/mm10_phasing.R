library(tidyverse)
library(parallel)
library(furrr)
plan(multisession(workers=8))

project_dir ="/pellmanlab/logan/data/mouse_embryo_MN/"

ad_files = list.files(paste0(project_dir,"/allelic_depth"), full.names=T) %>% 
  `names<-`(basename(.) %>% gsub("\\..*","",.))

fam_files = list.files(paste0(project_dir,"/read_depth"), full.names = T) %>%
  `names<-`(basename(.) %>% gsub("\\..*","",.)) 
cts=fam_files %>% future_map_dfr(~read_tsv(.,comment="@"), .id="sample") %>%
  dplyr::rename_all(tolower)


bad_samples = cts %>% group_by(sample) %>%
  summarize(total=sum(count), outliers=sum(count/median(count)>10),n=n()) %>%
  filter(total/median(total)<0.5 | outliers/n>0.01) %>% select(sample)

# cts %>% group_by(sample) %>%
#   summarize(total=sum(count), outliers=sum(count/median(count)>10),n=n()) %>%
#   ggplot(aes(total, log10(outliers/n), label=sample)) + geom_text(size=3,hjust="left")

# ad_files=ad_files[!(names(ad_files) %in% bad_samples$sample)]

phasing_dir='/pellmanlab/logan/data/references/mm10/129S_phasing'

chrs=paste0("chr",c(1:19,"X"))

all_snps=future_map_dfr(ad_files,read_tsv,.id = "sample", .progress = T)

all_snps = future_map_dfr(chrs, function(chr) {
  # chr_file=sprintf("%s/129S.%s.snps.tsv", phasing_dir, chr)
  map_dfr(ad_files, ~read_tsv(.,comment="@") %>% filter(contig==chr), .id="sample", .progress=T) %>%
    `names<-`(tolower(names(.))) %>%
    mutate(depth=refcount+altcount) %>%
    group_by(sample) %>%
    filter(between(mean(refcount)/mean(depth), 0.4,0.6)) %>%
    group_by(contig, position) %>% 
    summarize(n_samples=sum(depth>0), 
              ref_count=mean(refcount), alt_count=mean(altcount))
    # write_tsv(chr_file)
  # return(chr_file)
},.progress=T)
# chrs=unique(cts$contig) %>% grep("Y",., invert=T,value=T)


phasing_files=sapply(chrs, function(chr) {
  chr_file=sprintf("%s/129S.%s.snps.tsv", phasing_dir, chr)
  future_map_dfr(ad_files, ~read_tsv(.,comment="@") %>% filter(contig==chr), .id="sample", .progress=T) %>%
    `names<-`(tolower(names(.))) %>%
    mutate(depth=refcount+altcount) %>%
    group_by(sample) %>%
    # filter(between(mean(refcount)/mean(depth), 0.4,0.6)) %>% 
    mutate()
    group_by(contig, position) %>% 
    summarize(n_samples=sum(depth>0), 
              ref_count=mean(refcount), alt_count=mean(altcount)) %>% 
    write_tsv(chr_file)
  return(chr_file)
})

# map_dfr(ad_files, ~read_tsv(.,comment="@") %>% filter(contig=="chrX"), .id="sample") %>%
#   `names<-`(tolower(names(.))) %>%
#   mutate(depth=refcount+altcount) %>%
#   group_by(sample) %>%
#   mutate((mean(refcount)/mean(depth)) %>% between(0.4,0.6)) %>% 
#   group_by(contig, position) %>% 
#   summarize(f_samples=mean(depth > 0),
#             ref_count=mean(refcount), alt_count=mean(altcount))

all_snps = phasing_files %>% 
  future_map_dfr(read_tsv) %>% 
  mutate(imbalance=((ref_count-alt_count)/(ref_count+alt_count)),depth=ref_count+alt_count) %>% 
  mutate(f = ref_count/depth)
all_snps %>% sample_n(10000) %>% ggplot() + geom_point(aes(depth/median(depth), f), alpha=0.1) 
all_snps %>% summarize(mean(ref_count+alt_count<=3))

small_bin_size=5e4
bin_size=1e6
max_scale_factor=2.5
filtered_snps = all_snps %>% 
  filter(between(f, 1/max_scale_factor, 1 - 1/max_scale_factor)) %>% 
  filter(between(depth/median(depth), 1/max_scale_factor, max_scale_factor))
nrow(filtered_snps)/nrow(all_snps)
filtered_snps %>% summarize(mean(f>0.5))
filtered_snps %>%
  mutate(bin=floor(position/small_bin_size)*small_bin_size) %>%
  group_by(contig, bin) %>% 
  summarize() %>%
  mutate(bin=floor(bin/bin_size)*bin_size) %>%
  group_by(contig, bin) %>% 
  filter(n()>=2) %>% summarize() %>% count()

tbl = read_tsv("/pellmanlab/logan/data/mouse_bulk_seq/129S1.variants.table",
              col_names = c("contig", "position", "ref", "alt", "GTA", "GTB")) %>% 
  filter(GTB=="1/1") %>% semi_join(filtered_snps)

tbl %>% write_tsv(file.path(phasing_dir,"129S1.variants.filtered.table"))

# binned_ads=ads %>% filter(sample=="EM_SC_190808_M1_6A") %>%
#   semi_join(allelic_imbalance, by=c("contig","position")) %>% 
#   mutate(bin=floor(position/1e6)) %>% group_by(contig, bin) %>% 
#   mutate(n=n(), dp=alt_count+ref_count, baf=alt_count/dp) %>% 
#   summarize_all(~mean(.,na.rm=T)) %>% 
#   filter(n>10)
# binned_ads %>% ggplot(aes(bin, baf))+geom_point()
# 
# allelic_imbalance %>% ungroup %>% summarize(median((ref_count - alt_count)/(ref_count+alt_count),na.rm=T))
# 
# mean_var_cov %>% 
#   filter(abs(ref_count - alt_count)/(ref_count+alt_count) <= 0.1)
# 
# binned_ad = ads %>%
#   mutate(bin= floor(position/1e6)) %>% 
#   group_by(contig,bin) %>% 
#   mutate(count=ref_count+alt_count) %>% 
#   semi_join(bins_w_depth, by=c("contig","bin")) %>% 
#   summarize(mean_frac = mean(alt_count/count, na.rm=T), n=n()) 
# 
# binned_ad %>%
#   filter(n>100) %>% 
#   ggplot() + geom_point(aes(bin, mean_frac, color=n)) + facet_wrap(~contig)

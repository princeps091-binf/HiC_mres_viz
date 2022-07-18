library(tidyverse)
library(vroom)
library(parallel)
library(GenomicRanges)
library(igraph)
options(scipen = 999999999)
res_set <- c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num <- c(1e6,5e5,1e5,5e4,1e4,5e3)
names(res_num)<-res_set
########################################################################################################
## Utils. Fn.
cage_tbl_coord_build_fn<-function(cage_tbl,ID_col){
  
  cage_coord<-cage_tbl %>% dplyr::select(ID_col) %>% unlist
  cage_start<-as.numeric(unlist(lapply(strsplit(unlist(lapply(strsplit(unlist(lapply(strsplit(cage_coord,split = ':'),'[',2)),split=','),'[',1)),split = '\\.'),'[',1)))
  cage_end<-as.numeric(unlist(lapply(strsplit(unlist(lapply(strsplit(unlist(lapply(strsplit(cage_coord,split = ':'),'[',2)),split=','),'[',1)),split = '\\.'),'[',3)))
  cage_chr<-unlist(lapply(strsplit(cage_coord,split = ':'),'[',1))
  return(cage_tbl%>%mutate(chr=cage_chr,start=cage_start,end=cage_end))
  
}

cage_mean_compute_fn<-function(cage_tbl){
  print("compute m")
  cl<-makeCluster(5)
  tmp_m<-parallel::parApply(cl,X = as.matrix(cage_tbl[,-1]),MARGIN = 1,function(x){
    mean(x)
  })
  stopCluster(cl)
  rm(cl)
  cage_tbl<-cage_tbl%>%
    mutate(score=tmp_m)%>%
    filter(score>0)
  return(cage_tbl)
}

########################################################################################################
## Cell-line sample ID
H1<-c("CNhs14067","CNhs14068","CNhs13964")
GM12878<-c('CNhs12331','CNhs12332','CNhs12333')
HMEC<-c('CNhs11077','CNhs11382','CNhs12032')
########################################################################################################
#TSS
cage_tss_tbl<-vroom("~/Documents/multires_bhicect/data/epi_data/hg19.cage_peak_phase1and2combined_tpm_ann.osc.txt",comment = '#',col_select = contains(c("Annotation",H1,HMEC,GM12878)))
#Enhancers
cage_enh_tbl<-vroom("~/Documents/multires_bhicect/data/epi_data/human_permissive_enhancers_phase_1_and_2_expression_tpm_matrix.txt",comment = '#',col_select = contains(c("Id",H1,HMEC,GM12878)))
########################################################################################################
out_file<-"./data/H1_CAGE_peak_track_tbl.Rda"
#Enhancers
tmp_enh_tbl<-cage_enh_tbl %>%
  dplyr::select(c("Id",all_of(H1))) %>% 
  filter(if_any(where(is.numeric), ~ .x > 0)) %>% 
  cage_mean_compute_fn(.) %>% 
  mutate(chrom=str_split_fixed(Id,":|-",3)[,1],
         start=str_split_fixed(Id,":|-",3)[,2],
         end=str_split_fixed(Id,":|-",3)[,3]) %>% 
  filter(!(is.na(chrom))) %>% 
  dplyr::select(chrom,start,end,score) %>% 
  mutate(start=as.integer(start),
         end=as.integer(end))
#Tss
tmp_tss_tbl<-cage_tss_tbl %>%
  dplyr::select(contains(c("Annotation",all_of(H1)))) %>% 
  filter(if_any(where(is.numeric), ~ .x > 0)) %>% 
  cage_mean_compute_fn(.) %>% 
  mutate(chrom=str_split_fixed(`00Annotation`,":|,|\\.\\.",4)[,1],
         start=str_split_fixed(`00Annotation`,":|,|\\.\\.",4)[,2],
         end=str_split_fixed(`00Annotation`,":|,|\\.\\.",4)[,3]) %>% 
  filter(!(grepl("STAT",`00Annotation`))) %>% 
  dplyr::select(chrom,start,end,score) %>% 
  mutate(start=as.integer(start),
         end=as.integer(end))

tmp_tbl<-tmp_enh_tbl %>% 
  mutate(set="enhancer") %>% 
  bind_rows(.,tmp_tss_tbl %>% 
              mutate(set="TSS"))
## Dealing with overlapping peaks
cage_Grange<-GRanges(seqnames=tmp_tbl$chrom,
                     ranges = IRanges(start=tmp_tbl$start,
                                      end=tmp_tbl$end
                     ))
mcols(cage_Grange)<-tibble(score=tmp_tbl$score,set=tmp_tbl$set)

peak_clusters<-findOverlaps(cage_Grange,cage_Grange) %>% 
  as_tibble %>% 
  filter(queryHits != subjectHits)
peak_cmpnt<-components(graph_from_data_frame(peak_clusters))
big_peak_cluster_idx<-which.max(peak_cmpnt$csize)
big_peak_cluster_GRange<-cage_Grange[as.numeric(names(which(peak_cmpnt$membership == big_peak_cluster_idx)))]
### Illustrate configuration of these overlapping peaks
big_peak_cluster_GRange %>% as_tibble %>% 
  mutate(ID=1:n()) %>% 
  ggplot(.,aes(group=ID,y=score,yend=score,x=start,xend=end,color=set))+
  geom_segment()

findOverlaps(disjoin(cage_Grange),cage_Grange) %>% 
  as_tibble() %>% 
  mutate(score=mcols(cage_Grange)$score[subjectHits]) %>% 
  group_by(queryHits) %>% 
  summarise(m=mean(score))

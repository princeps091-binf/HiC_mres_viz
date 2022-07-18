library(tidyverse)
library(plotgardener)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)

options(scipen = 999999999)
res_set <- c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num <- c(1e6,5e5,1e5,5e4,1e4,5e3)
names(res_num)<-res_set
########################################################################################################
## Util. fn. to load R-objects
get_tbl_in_fn<-function(tmp_file){
  out_tbl<-get(base::load(tmp_file))
  tmp_obj<-names(mget(base::load(tmp_file)))
  rm(list=tmp_obj)
  rm(tmp_obj)
  return(out_tbl)
}
########################################################################################################
cl_heat_dat_file<-"./data/chr22_50kb_34_561_38100000_39750000.Rda"
CAGE_tbl_file<-"./data/H1_enh_peak_track_tbl.Rda"
tmp_cl<-"50kb_34_561_38100000_39750000"
tmp_cl_res<-str_split_fixed(tmp_cl,pattern = "_",2)[,1]
tmp_res_set<-names(which(res_num<=res_num[tmp_cl_res]))

tmp_f_dat<-get_tbl_in_fn(cl_heat_dat_file)
cage_tbl<-get_tbl_in_fn(CAGE_tbl_file) %>% 
  dplyr::rename(chromstart=start,chromend=end)
cage_tbl<-get_tbl_in_fn(CAGE_tbl_file)
  
p_breaks<-c(unlist(lapply(seq_along(tmp_res_set),function(x){
  seq((x-1)*100,(x-1)*100 +100,length.out = 101)[-101]
})),(length(res_set)-1)*100 +100)

p_color<-rev(RColorBrewer::brewer.pal(n=length(tmp_res_set),name = "Set1"))
p_col<-unlist(lapply(seq_along(tmp_res_set),function(x){
  colorRampPalette(c("black",p_color[x]))(100)
}))



pageCreate(width = 3.25, height = 3.25, default.units = "inches",showGuides = T)

test<-plotHicSquare(
  resolution = 5000,
  data = tmp_f_dat,
  chrom = "chr22", chromstart = 38100000, chromend = 39795000,
  assembly = "hg19",
  x = 0.25, y = 0.25, width = 2.5, height = 2.5, default.units = "inches",
  draw=F)
heat_grob_id<-names(test[["grobs"]][["children"]])[2]
tmp_ego<-test[["grobs"]][["children"]][[heat_grob_id]]$x
tmp_alter<-test[["grobs"]][["children"]][[heat_grob_id]]$y
tmp_col<-test[["grobs"]][["children"]][[heat_grob_id]]$gp$fill

tmp_hic_col<-tibble(ego=as.numeric(str_split_fixed(tmp_ego,"native",n=2)[,1]),
                    alter=as.numeric(str_split_fixed(tmp_alter,"native",n=2)[,1]),
                    fill=tmp_col)

tmp_col<-p_col

names(tmp_col)<-as.character(1:300)

tmp_conv<-tmp_f_dat %>% 
  bind_rows(tibble(ego=tmp_f_dat$alter,alter=tmp_f_dat$ego,color=tmp_f_dat$color)) %>% 
  mutate(code=tmp_col[as.character(round(color))])

tmp_hic_col<-tmp_hic_col %>% 
  left_join(.,tmp_conv)

test[["grobs"]][["children"]][[heat_grob_id]]$gp$fill<-tmp_hic_col$code

pagePlotPlace(
  plot=test,
  x = 0.25, y = 0.25, width = 2.5, height = 2.5,default.units = "inches"
)

annoGenomeLabel(
  plot = test,
  axis = "y",
  x = 0.075, y = 0.25, 
  width = 2.5, height = 0.25, 
  default.units = "inches",scale = "Mb"
)

plotSignal(
  data = cage_tbl,
  chrom = "chr22", chromstart = 38100000, chromend = 39795000,
  assembly = "hg19",
  x = 0.25, y = 2.8, width = 2.5, height = 0.2,
  just = c("left", "top"), default.units = "inches")

plotGenes(
  chrom = "chr22", chromstart = 38100000, chromend = 39795000,
  assembly = "hg19",
  x = 0.25, y = 3.05, width = 2.5, height = 0.2,
  just = c("left", "top"), default.units = "inches"
)

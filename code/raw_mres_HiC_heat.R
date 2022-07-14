library(Matrix)
library(data.tree)
library(caret)
library(tidyverse)
library(parallel)
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
## Fn. mapping low-res bins to Hi-res bins
#convert bins to higher resolution
bin_map_res_fn<-function(chr_spec_res_99,res_num){
  #examine for every 5kb bin the path convergence and divergence
  #given multi resolution clustering we need to map the correspondence of bins across multiple resolutions
  #all resolutions considered for the chromosome
  chr_res_set<-unique(unlist(lapply(strsplit(names(chr_spec_res_99$cl_member),split='_'),'[',1)))
  chr_res_set<-rev(names(sort(res_num[chr_res_set])))
  #Starting with the coarser resolution
  bin_top<-unique(unlist(chr_spec_res_99$cl_member[grep(chr_res_set[1],names(chr_spec_res_99$cl_member))]))
  #build the temporary list which will contain the mapping between pairs of bin sets at consecutive resolutions
  bin_res_map<-vector('list',length(chr_res_set))
  names(bin_res_map)<-chr_res_set
  #initiate the coarser resolution
  bin_res_map[[chr_res_set[1]]]<-as.data.frame(bin_top)
  colnames(bin_res_map[[chr_res_set[1]]])<-paste0('bin_',chr_res_set[1])
  #loop through the remaining resolutions
  for (r in chr_res_set[-1]){
    #derive the anticipated bins at the higher resolution fron the coarser resolution binning 
    r_bin<-lapply(as.character(bin_res_map[[which(chr_res_set==r)-1]][[grep(chr_res_set[which(chr_res_set==r)-1],names(bin_res_map[[which(chr_res_set==r)-1]]))]]),function(x){
      tmp<-seq(as.numeric(x),as.numeric(x)+res_num[chr_res_set[which(chr_res_set==r)-1]],by=res_num[r])
      return(tmp[-length(tmp)])
    })
    #build dataframe containing one column with the derived higher resolution binning and
    #one column with the mapped coarser resolution binning
    tmp_df<-data.frame(rep(bin_res_map[[chr_res_set[which(chr_res_set==r)-1]]][[paste0('bin_',chr_res_set[which(chr_res_set==r)-1])]],unlist(lapply(r_bin,length))))
    colnames(tmp_df)<-paste0('bin_',chr_res_set[which(chr_res_set==r)-1])
    tmp_df[[paste0('bin_',r)]]<-unlist(r_bin)
    bin_res_map[[r]]<-tmp_df
  }
  #build the full mapping dataframe by sequential joining of the previously generated mapping tables.
  bin_res_map_df<-data.frame(bin_top)
  colnames(bin_res_map_df)[1]<-paste0('bin_',chr_res_set[1])
  for(i in seq_along(bin_res_map)[-1]){
    
    bin_res_map_df<-left_join(bin_res_map_df,bin_res_map[[i]])
    
  }
  rm(i,r,tmp_df,r_bin,bin_res_map,bin_top)
  return(bin_res_map_df)
}
#convert low-resolution interactions to highest resolution interactions
convert_to_hires_fn<-function(tmp_dat,tmp_res,res_num,hi_res){
  fn_env<-environment()
  cl<-makeCluster(5)
  clusterEvalQ(cl, {
    library(dplyr)
    library(tidyr)
    print("node ready")
  })
  clusterExport(cl,c("tmp_dat","tmp_res","res_num","hi_res"),envir = fn_env)
  
  out_tbl<-do.call(bind_rows,parLapply(cl,1:nrow(tmp_dat),function(x){
    r_bin<-lapply(tmp_dat[x,c(1,2)],function(i){
      tmp<-seq(i,i+res_num[tmp_res],by=hi_res)
      return(tmp[-length(tmp)])
    })
    tmp_df<-expand_grid(ego=r_bin$X1,alter=r_bin$X2)
    tmp_df<-tmp_df%>%mutate(raw=unlist(tmp_dat[x,3]))
    tmp_df<-tmp_df%>%mutate(pow=unlist(tmp_dat[x,4]))
    tmp_df<-tmp_df%>%mutate(res=tmp_res)
    tmp_df<-tmp_df%>%mutate(color=unlist(tmp_dat[x,5]))
    
    return(tmp_df)  
  }))
  stopCluster(cl)
  rm(cl)
  return(out_tbl)
}
#Produce interaction matrix
full_f_mat<-function(cl_mat,res){
  
  range_5kb<-range(unique(c(cl_mat$ego,cl_mat$alter)))
  bin_5kb<-seq(range_5kb[1],range_5kb[2],by=res)
  #add the bins not present in original Hi-C dataset
  #miss_bin<-bin_5kb[which(!(bin_5kb %in% unique(c(mat_df$X1,mat_df$X2))))]
  
  id_conv<-seq_along(bin_5kb)
  names(id_conv)<-bin_5kb
  
  cl_mat$ego_id<-id_conv[as.character(cl_mat$ego)]
  cl_mat$alter_id<-id_conv[as.character(cl_mat$alter)]
  
  #chr_mat<-sparseMatrix(i=cl_mat$ego_id,cl_mat$alter_id,x=sqrt(-log10(cl_mat$pois.pval)),symmetric = T)
  chr_mat<-sparseMatrix(i=cl_mat$ego_id,cl_mat$alter_id,x=cl_mat$color,symmetric = T)
  
}

########################################################################################################
hub_file<-"~/Documents/multires_bhicect/Bootstrapp_fn/data/DAGGER_tbl/trans_res/H1_union_trans_res_dagger_tbl.Rda"
spec_res_folder<-"~/Documents/multires_bhicect/data/H1/Dekker/spec_res/"
dat_folder<-"~/Documents/multires_bhicect/data/H1/Dekker/"
########################################################################################################
# Load considered chromosome
chromo<-"chr22"
chr_spec_res<-get_tbl_in_fn(paste0(spec_res_folder,chromo,"_spec_res.Rda"))
chr_bpt<-FromListSimple(chr_spec_res$part_tree)
node_ancestor<-chr_bpt$Get(function(x){x$Get('name',traversal='ancestor')})
node_ancestor<-lapply(node_ancestor,'[',-1)
node_lvl<-sort(chr_bpt$Get('level'))[-c(1)]
# Build resolution conversion table
bin_res_map_df<-bin_map_res_fn(chr_spec_res,res_num)
bin_res_map_df[,1]<-as.numeric(as.character(bin_res_map_df[,1]))

#-----------------------
chr_hub_tbl<-get_tbl_in_fn(hub_file) %>%
  filter(chr==chromo)
tmp_cl<-"50kb_34_561_38100000_39750000"
tmp_cl_res<-str_split_fixed(tmp_cl,pattern = "_",2)[,1]
tmp_res_set<-names(which(res_num<=res_num[tmp_cl_res]))

# Subset conversion table to only consider bins contained in cluster of interest
cl_bin<-as.numeric(chr_spec_res$cl_member[[tmp_cl]])
tmp_col<-grep(tmp_cl_res[1],colnames(bin_res_map_df),value = T)
cl_bin_convert_tbl<-bin_res_map_df%>%
  filter(!!rlang::sym(tmp_col) >= min(cl_bin) & !!rlang::sym(tmp_col) <= max(cl_bin))
#----------------------------------------------------------------------------------------------------
#Extract the interactions of every constitutive child-cluster
dat_file<-'~/Documents/multires_bhicect/data/H1/Dekker/'

chr_dat_l<-lapply(tmp_res_set,function(x)read_delim(file = paste0(dat_folder,x,'/',chromo,'.txt'),delim = '\t',col_names = F))
names(chr_dat_l)<-tmp_res_set
chr_dat_l<-lapply(chr_dat_l,function(x){
  x%>%filter(!(is.nan(X3)))%>%filter(X1!=X2)
})
chr_dat_l<-lapply(chr_dat_l,function(x){
  preprocessParams <- BoxCoxTrans(x$X3,na.rm = T)
  x <- data.frame(x,weight=predict(preprocessParams, x$X3))
  x$weight<-x$weight+(1-min(x$weight,na.rm = T))
  return(x)
})
# Produce color-scale separating each resolution into separate color-channel
chr_dat_l<-lapply(seq_along(tmp_res_set),function(x){
  tmp_dat<-chr_dat_l[[tmp_res_set[x]]]
  toMin<-(x-1)*100 +1
  toMax<-(x-1)*100 +99
  tmp_dat$color<-toMin+(tmp_dat$weight-min(tmp_dat$weight))/(max(tmp_dat$weight)-min(tmp_dat$weight))*(toMax-toMin)
  return(tmp_dat)
})
names(chr_dat_l)<-tmp_res_set
#----------------------------------------------------------------------------------------------------
#Filter the Hi-C dataset to only contain parent cluster bins

cl_dat_l<-lapply(tmp_res_set,function(x){
  tmp_col<-grep(x,colnames(bin_res_map_df),value = T)
  tmp_res_max<-cl_bin_convert_tbl %>% 
    summarise(max=max(!!rlang::sym(tmp_col)),
              min=min(!!rlang::sym(tmp_col)))
  
  chr_dat_l[[x]]%>%
    filter(X1 <= tmp_res_max$max & X1 >= tmp_res_max$min & X2 <= tmp_res_max$max & X2 >= tmp_res_max$min)
})
names(cl_dat_l)<-tmp_res_set
#----------------------------------------------------------------------------------------------------
#Convert all data to 5kb resolution


cl_dat_hires_l<-lapply(names(cl_dat_l)[-length(cl_dat_l)],function(r){
  message(r)
  return(convert_to_hires_fn(cl_dat_l[[r]],r,res_num,5000))

})
cl_dat_hires_l[[length(cl_dat_l)]]<-tibble(cl_dat_l[[length(cl_dat_l)]]) %>% 
  dplyr::rename(ego=X1,alter=X2,raw=X3,pow=weight)
names(cl_dat_hires_l)<-names(cl_dat_l)

tmp_f_dat<-cl_dat_hires_l[[length(cl_dat_hires_l)]]
for(i in (length(cl_dat_hires_l)):2){
  message(i)
  tmp_f_dat<- tmp_f_dat%>% 
    dplyr::select(ego,alter,color) %>% 
    full_join(.,cl_dat_hires_l[[i-1]] %>% 
                dplyr::select(ego,alter,color) %>% 
                rename(color.b=color)) %>% 
    mutate(color=ifelse(is.na(color),color.b,color)) %>% 
    dplyr::select(-color.b)
}


p_breaks<-c(unlist(lapply(seq_along(tmp_res_set),function(x){
  seq((x-1)*100,(x-1)*100 +100,length.out = 101)[-101]
})),(length(res_set)-1)*100 +100)

p_color<-rev(RColorBrewer::brewer.pal(n=length(tmp_res_set),name = "Set1"))
p_col<-unlist(lapply(seq_along(tmp_res_set),function(x){
  colorRampPalette(c("black",p_color[x]))(100)
}))

cl_hires_bin<-unique(c(tmp_f_dat$ego,tmp_f_dat$alter))


cl_f_mat<-full_f_mat(tmp_f_dat,5000)

tick_pos<-which(cl_hires_bin %% 5e5 == 0)/nrow(cl_f_mat)
tick_label<-paste0(cl_hires_bin[which(cl_hires_bin %% 5e5 == 0)]/1e6,"Mb")

image(as.matrix(cl_f_mat),col=p_col,breaks=p_breaks,axes = FALSE)
axis(1, at = tick_pos,
     labels = tick_label,
     
)
legend(x = "top",
       inset = c(0, -0.12), # You will need to fine-tune the first
       # value depending on the windows size
       ncol=length(tmp_res_set),
       legend = tmp_res_set, 
       fill = p_color,
       xpd = TRUE)

png(filename = "./img/raw_mres_heat.png", width =40,height = 50,units = 'mm',type='cairo',res=1000)
par(mar=c(1.25,0,1.25,0),cex.axis = 0.4,mgp=c(3, 0.25, 0))

image(as.matrix(cl_f_mat),col=p_col,breaks=p_breaks,axes = FALSE)
axis(1, at = tick_pos,
     labels = tick_label,
     
)
legend(x = "top",
       inset = c(0, -0.12), # You will need to fine-tune the first
       # value depending on the windows size
       ncol=length(tmp_res_set),
       legend = tmp_res_set, 
       fill = p_color,
       cex = 0.4, # Change legend size
       xpd = TRUE)

dev.off()

save(tmp_f_dat,file = paste0("./data/",chromo,"_",tmp_cl,".Rda"))
#==============================================
library(plotgardener)
pageCreate(width = 3.25, height = 3.25, default.units = "inches",showGuides = F)

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

## Plot gene track
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
plotGenes(
  chrom = "chr22", chromstart = 38100000, chromend = 39795000,
  assembly = "hg19",
  x = 0.25, y = 2.75, width = 2.5, height = 0.3,
  just = c("left", "top"), default.units = "inches"
)

png(width = 3.25, height = 3.25,units = "in",filename = "./img/test_gardener.png",res=2000)
pageCreate(width = 3.25, height = 3.25, default.units = "inches",showGuides = F)
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
plotGenes(
  chrom = "chr22", chromstart = 38100000, chromend = 39795000,
  assembly = "hg19",
  x = 0.25, y = 2.75, width = 2.5, height = 0.3,
  just = c("left", "top"), default.units = "inches"
)
plotSignal(
  data = IMR90_ChIP_H3K27ac_signal,
  chrom = "chr21", chromstart = 28000000, chromend = 30300000,
  assembly = "hg19",
  x = 0.5, y = 2.75, width = 2, height = 0.5,
  just = c("left", "top"), default.units = "inches"
)
dev.off()
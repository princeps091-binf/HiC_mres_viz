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

#Extract children HiC data
extract_cl_dat_fn<-function(child_cl_lvl_tbl,cl_dat_l){
  fn_env<-environment()
  cl<-makeCluster(5)
  clusterEvalQ(cl, {
    library(dplyr)
  })
  clusterExport(cl,c("child_cl_lvl_tbl",'cl_dat_l'),envir = fn_env)
  cl_child_inter<-parLapply(cl,1:nrow(child_cl_lvl_tbl),function(r){
    tmp_res<-child_cl_lvl_tbl$res[r]
    tmp_bin<-as.numeric(child_cl_lvl_tbl$bins[[r]])
    tmp_ref<-cl_dat_l[[tmp_res]]
    tmp_dat<-tmp_ref%>%filter(X1 %in% tmp_bin & X2 %in% tmp_bin)
    
    return(as_tibble(tmp_dat))
  })
  stopCluster(cl)
  rm(cl)
  
  return(child_cl_lvl_tbl %>% 
           mutate(HiC=cl_child_inter))
  
}

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

tmp_cl_child<-c(tmp_cl,names(which(unlist(lapply(node_ancestor, function(x){
  tmp_cl %in% x
})))))
child_cl_lvl_tbl<-tibble(node=tmp_cl_child,lvl=node_lvl[tmp_cl_child],bins=chr_spec_res$cl_member[tmp_cl_child]) %>% 
  mutate(res=str_split_fixed(node,"_",2)[,1])

#----------------------------------------------------------------------------------------------------
#Extract the interactions of every constitutive child-cluster
dat_file<-'~/Documents/multires_bhicect/data/H1/'

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
#Filter the Hi-C dataset to only contain parent and children cluster bins

cl_dat_l<-lapply(tmp_res_set,function(x){
  tmp_res_max<-child_cl_lvl_tbl %>% 
    filter(res==x) 
  bin_range<-range(as.numeric(unique(unlist(tmp_res_max$bins))))
  chr_dat_l[[x]]%>%
    filter(X1 <= bin_range[2] & X1 >= bin_range[1] & X2 <= bin_range[2] & X2 >= bin_range[1])
})
names(cl_dat_l)<-tmp_res_set
#----------------------------------------------------------------------------------------------------
#Extract original HiC data for every cluster at the resolution they are found
child_cl_lvl_tbl<-extract_cl_dat_fn(child_cl_lvl_tbl,cl_dat_l)
child_cl_lvl_tbl<-child_cl_lvl_tbl %>% 
  mutate(i=str_split_fixed(node,"_",4)[,3]) %>% 
  filter(i>0) %>% 
  dplyr::select(-i)
#Convert all interaction to 5kb resolution 
cl_hires_dat_l<-lapply(1:nrow(child_cl_lvl_tbl),function(i){
  message(round(i/nrow(child_cl_lvl_tbl),digits = 3))
  tmp_dat<-child_cl_lvl_tbl$HiC[[i]]
  tmp_res<-child_cl_lvl_tbl$res[i]
  hi_res<-5000
  convert_to_hires_fn(tmp_dat,tmp_res,res_num,hi_res)
})
child_cl_lvl_tbl<-child_cl_lvl_tbl %>% 
  mutate(HiC.hires=cl_hires_dat_l)

lvl_seq<-sort(unique(child_cl_lvl_tbl$lvl),decreasing = T)

tmp_lvl_dat<-child_cl_lvl_tbl %>% 
  filter(lvl==lvl_seq[1])
tmp_seed<-do.call(bind_rows,tmp_lvl_dat$HiC.hires)

for(i in lvl_seq[-1]){
  message(i)
  tmp_up_lvl_dat_tbl<-child_cl_lvl_tbl %>% 
    filter(lvl==i)
  tmp_up_lvl_dat<-do.call(bind_rows,tmp_up_lvl_dat_tbl$HiC.hires)
  
  tmp_seed<- tmp_seed%>% 
    dplyr::select(ego,alter,color) %>% 
    full_join(.,tmp_up_lvl_dat %>% 
                dplyr::select(ego,alter,color) %>% 
                rename(color.b=color)) %>% 
    mutate(color=ifelse(is.na(color),color.b,color)) %>% 
    dplyr::select(-color.b)
}


p_breaks<-c(unlist(lapply(seq_along(res_set),function(x){
  seq((x-1)*100,(x-1)*100 +100,length.out = 101)[-101]
})),(length(res_set)-1)*100 +100)

p_color<-rev(RColorBrewer::brewer.pal(n=length(res_set),name = "Set1"))
p_col<-unlist(lapply(seq_along(res_set),function(x){
  colorRampPalette(c("white",p_color[x]))(100)
}))
cl_f_mat<-full_f_mat(tmp_seed,5000)
image(as.matrix(cl_f_mat),col=p_col,breaks=p_breaks)

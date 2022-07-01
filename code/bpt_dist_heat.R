library(Matrix)
library(data.tree)
library(GenomicRanges)
library(caret)
library(tidyverse)
library(parallel)
library(igraph)
library(viridis)
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
  chr_mat<-sparseMatrix(i=cl_mat$ego_id,cl_mat$alter_id,x=cl_mat$bpt_dist)
  return(chr_mat)
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
cl_bin<-as.numeric(chr_spec_res$cl_member[[tmp_cl]])

tmp_col<-grep(tmp_cl_res[1],colnames(bin_res_map_df),value = T)
cl_bin_convert_tbl<-bin_res_map_df%>%
  filter(!!rlang::sym(tmp_col) >= min(cl_bin) & !!rlang::sym(tmp_col) <= max(cl_bin))

tmp_cl_child<-c(tmp_cl,names(which(unlist(lapply(node_ancestor, function(x){
  tmp_cl %in% x
})))))
cl_leaves<-tmp_cl_child[tmp_cl_child %in% chr_bpt$Get('name',filterFun = isLeaf)]
cl_node_set<-unique(c(cl_leaves,unique(unlist(node_ancestor[cl_leaves])))) 

chr_bpt2<-FromListSimple(chr_spec_res$part_tree)

Prune(chr_bpt2, function(x) x$name %in% cl_node_set)

g_bpt<-as.igraph.Node(chr_bpt2,directed = T,direction = 'climb')

bpt_dist_mat<-distances(g_bpt,cl_leaves,cl_leaves)

cl_leaf_tbl<-tibble(leaf=cl_leaves,bin=chr_spec_res$cl_member[cl_leaves]) %>% 
  mutate(res=str_split_fixed(leaf,"_",2)[,1]) %>% 
  mutate(GRange=pmap(list(chromo,bin,res),function(chromo,bin,res){
    
    IRanges::reduce(GRanges(seqnames=chromo,
                            ranges = IRanges(start=as.numeric(bin),
                                             end=as.numeric(bin) + res_num[res] - 1
                            )))
    
    
  }))
hires_bin_tbl<-cl_bin_convert_tbl %>% 
  dplyr::select(contains("_5kb")) %>% 
  as_tibble %>% 
  dplyr::rename(start=bin_5kb) %>% 
  mutate(end=start + 4999) %>% 
  mutate(GRange=pmap(list(chromo,start,end),function(chromo,start,end){
    
    IRanges::reduce(GRanges(seqnames=chromo,
                            ranges = IRanges(start=start,
                                             end=end
                            )))
    
    
  }))
hires_bin_GrangeL<-GRangesList(hires_bin_tbl$GRange)
cl_leaf_GrangeL<-GRangesList(cl_leaf_tbl$GRange)
bin_to_leaf_tbl<-findOverlaps(hires_bin_GrangeL,cl_leaf_GrangeL) %>% 
  as_tibble %>% 
  mutate(bin=hires_bin_tbl$start[queryHits],
         leaf=cl_leaf_tbl$leaf[subjectHits]) %>% 
  dplyr::select(bin,leaf)

hi_res_inter_tbl<-expand_grid(ego=bin_to_leaf_tbl$bin,alter=bin_to_leaf_tbl$bin) %>% 
  left_join(.,bin_to_leaf_tbl,by=c("ego"="bin")) %>%
  dplyr::rename(ego.leaf=leaf) %>% 
  left_join(.,bin_to_leaf_tbl,by=c("alter"="bin")) %>% 
  dplyr::rename(alter.leaf=leaf)
hi_res_inter_tbl<-hi_res_inter_tbl %>% 
  mutate(bpt_dist=log(1/(1+bpt_dist_mat[as.matrix(hi_res_inter_tbl[,c(3,4)])])))
cl_f_mat<-full_f_mat(hi_res_inter_tbl,5000)
image(as.matrix(cl_f_mat),col=viridis(100))

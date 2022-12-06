
subc_ini <- function(){
  suppressMessages(require(Seurat))
  suppressMessages(require(SeuratDisk))
  suppressMessages(require(raster))
  suppressMessages(require(stringr))
  suppressMessages(require(dplyr))
  suppressMessages(require(parallel))
  suppressMessages(require(pbapply))
  suppressMessages(require(ggplot2))
  suppressMessages(require(readr))
  suppressMessages(require(ggpubr))
  suppressMessages(require(monocle3))
  suppressMessages(library(nleqslv))
  suppressMessages(library(magrittr))
  suppressMessages(library(anndata))
  suppressMessages(library(factoextra))
  suppressMessages(require(circlize))
  suppressMessages(require(ComplexHeatmap))
  message("Initialize Done")
}



sc_fitting <- function(gene_list,psd_bulk_list, v_n, v_c, label_vec = c(2,2,2,2), wash_vec = c(0,2,4,6)){
  for(i in names(psd_bulk_list)){
    for(j in colnames(psd_bulk_list[[i]])){
      if(j != "gene"){
        psd_bulk_list[[i]] <- psd_bulk_list[[i]][psd_bulk_list[[i]][[j]] > 0,]
      } 
    }
    gene_list %<>% intersect(psd_bulk_list[[i]]$gene)
  }
  
  for(i in names(psd_bulk_list)){
    psd_bulk_list[[i]] %<>% subset(subset = gene %in% gene_list)
  }
  #check if labeling times are same
  if(length(label_vec) != length(wash_vec) | ncol(psd_bulk_list[[i]]) != length(label_vec)+1 | ncol(psd_bulk_list[[i]]) != length(label_vec)+1) stop("incorrect vector size")
  if(sum(mean(label_vec) != label_vec) > 0 ) stop("Labeling time inconsitency")
  #Check if only one label-only sample
  if(sum(wash_vec == 0) != 1 ) stop("No labeling-only sample(or too many)")
  labeling_time = label_vec[1]
  
  
  beta_list <- apply(data.frame(gene = gene_list),1,
                     FUN = function(x){
                       time_vec2 <- wash_vec
                       expr_vec <- rep(0,4)
                       for(i in 1:length(wash_vec)){
                         if(wash_vec[i] == 0) expr_vec[i] <- 0
                         else{
                           expr_vec[i] <- log(psd_bulk_list$RNA[[paste0(labeling_time,"h_labeling")]][psd_bulk_list$RNA$gene == x[["gene"]]]) -
                             log(psd_bulk_list$RNA[[paste(paste0(labeling_time,"h_labeling_"),
                                                          wash_vec[i],"h_wash",
                                                          sep = "")]][psd_bulk_list$RNA$gene == x[["gene"]]])
                         }
                       }
                       lm_res <- lm(expr~time,data.frame(time = time_vec2, expr = expr_vec))
                       return(list(gene = x[["gene"]],
                                   beta = lm_res[["coefficients"]][["time"]],
                                   beta_r_sq = summary(lm_res)[["r.squared"]], 
                                   beta_adj_r_sq = summary(lm_res)[["adj.r.squared"]],
                                   beta_expr_vec = log(psd_bulk_list$RNA[[paste0(labeling_time,"h_labeling")]][psd_bulk_list$RNA$gene == x[["gene"]]]) - expr_vec))
                     })
  label_para_df <- do.call(rbind,beta_list) %>% as.data.frame()
  for(i in colnames(label_para_df)[1:4]) label_para_df[[i]] <- unlist(label_para_df[[i]])
  
  label_para_df[["alpha"]] <- apply(label_para_df,1,
                                    FUN = function(x){
                                      if(is.na(x[["beta"]]) | as.numeric(x[["beta"]]) < 0){
                                        NA
                                      }else{
                                        as.numeric(x[["beta"]]) * psd_bulk_list$RNA[[paste0(labeling_time,"h_labeling")]][psd_bulk_list$RNA$gene == x[["gene"]]]/
                                          (1 - exp(-1 * as.numeric(x[["beta"]])))
                                      }
                                    })
  
  ked_list <- apply(data.frame(gene = gene_list),1,
                    FUN = function(x){
                      time_vec2 <- wash_vec
                      expr_vec <- rep(0,4)
                      for(i in 1:length(wash_vec)){
                        if(wash_vec[i] == 0) expr_vec[i] <- 0
                        else{
                          expr_vec[i] <- log(psd_bulk_list$Nucleus[[paste0(labeling_time,"h_labeling")]][psd_bulk_list$Nucleus$gene == x[["gene"]]]) -
                            log(psd_bulk_list$Nucleus[[paste(paste0(labeling_time,"h_labeling_"),
                                                             wash_vec[i],"h_wash",
                                                             sep = "")]][psd_bulk_list$Nucleus$gene == x[["gene"]]])
                        }
                        
                      }
                      lm_res <- lm(expr~time,data.frame(time = time_vec2, expr = expr_vec))
                      return(list(gene = x[["gene"]],
                                  ked = lm_res[["coefficients"]][["time"]],
                                  ked_r_sq = summary(lm_res)[["r.squared"]], 
                                  ked_adj_r_sq = summary(lm_res)[["adj.r.squared"]],
                                  ked_expr_vec = log(psd_bulk_list$Nucleus[[paste0(labeling_time,"h_labeling")]][psd_bulk_list$Nucleus$gene == x[["gene"]]]) - expr_vec))
                    })
  
  wash_para_df <- do.call(rbind,ked_list) %>% as.data.frame()
  for(i in colnames(wash_para_df)[1:4]) wash_para_df[[i]] <- unlist(wash_para_df[[i]])
  
  tmp_df <- data.frame(label_para_df[,1:6],wash_para_df[,2:5])
  #Use 20h get beta_c/lambda
  #  tmp_df[["beta_c_lambda_ratio"]] <- apply(tmp_df,1,
  #                                           FUN = function(x){
  #                                             if(psd_bulk_list$Nucleus[["20h_labeling"]][psd_bulk_list$Nucleus$gene == x[["gene"]]] == 0 | psd_bulk_list$Cytoplasm[["20h_labeling"]][psd_bulk_list$Cytoplasm$gene == x[["gene"]]] == 0) NA
  #                                             else{
  #                                               psd_bulk_list$Nucleus[["20h_labeling"]][psd_bulk_list$Nucleus$gene == x[["gene"]]] * v_n[["20h_labeling"]] / (v_c[["20h_labeling"]] * psd_bulk_list$Cytoplasm[["20h_labeling"]][psd_bulk_list$Cytoplasm$gene == x[["gene"]]])
  #                                             }
  #                                             
  #                                           })
  ##Using LSS
  tmp_list <- apply(tmp_df, 1, 
                    FUN = function(x){
                      if(is.na(x[["beta"]]) | as.numeric(x[["ked"]]) < 0 | is.infinite(as.numeric(x[["ked"]])) ){# | is.na(x[["beta_c_lambda_ratio"]])
                        NA
                      }else{
                        n_index = psd_bulk_list$Nucleus$gene == x[["gene"]]
                        c_index = psd_bulk_list$Cytoplasm$gene == x[["gene"]]
                        r_index = psd_bulk_list$RNA$gene == x[["gene"]]
                        
                        n0 <- rep(psd_bulk_list$Nucleus[[paste0(labeling_time,"h_labeling")]][n_index],length(wash_vec))
                        ked <- rep(as.numeric(x[["ked"]]),  length(wash_vec))
                        b <- rep(as.numeric(x[["beta"]]),   length(wash_vec))
                        k <- rep(x[["beta_c_lambda_ratio"]],length(wash_vec))
                        c_vec = rep(0,length(wash_vec))
                        x_vec = rep(0,length(wash_vec))
                        n_vec = rep(0,length(wash_vec))
                        
                        c_r_vec = rep(0,length(wash_vec))
                        x_r_vec = rep(0,length(wash_vec))
                        n_r_vec = rep(0,length(wash_vec))
                        
                        t_vec = wash_vec
                        vc = rep(0,length(wash_vec))
                        vn = rep(0,length(wash_vec))
                        
                        for(i in 1:length(wash_vec)){
                          if(wash_vec[i] == 0){
                            c_vec[1] =  psd_bulk_list$Cytoplasm[[paste0(labeling_time,"h_labeling")]][c_index]
                            n_vec[1] =  psd_bulk_list$Nucleus[[paste0(labeling_time,"h_labeling")]][n_index]
                            x_vec[1] =  psd_bulk_list$RNA[[paste0(labeling_time,"h_labeling")]][r_index]
                            
                            c_r_vec[1] =  psd_bulk_list$Cytoplasm_r[[paste0(labeling_time,"h_labeling")]][c_index]
                            n_r_vec[1] =  psd_bulk_list$Nucleus_r[[paste0(labeling_time,"h_labeling")]][n_index]
                            x_r_vec[1] =  psd_bulk_list$RNA_r[[paste0(labeling_time,"h_labeling")]][r_index]
                            vc[1] = v_c[[paste0(labeling_time,"h_labeling")]]
                            vn[1] = v_n[[paste0(labeling_time,"h_labeling")]]
                          }else{
                            c_vec[i] = psd_bulk_list$Cytoplasm[[paste0(paste0(labeling_time,"h_labeling_"),wash_vec[i],"h_wash")]][c_index]
                            n_vec[i] = psd_bulk_list$Nucleus[[paste0(paste0(labeling_time,"h_labeling_"),wash_vec[i],"h_wash")]][n_index]
                            x_vec[i] = psd_bulk_list$RNA[[paste0(paste0(labeling_time,"h_labeling_"),wash_vec[i],"h_wash")]][r_index]
                            c_r_vec[i] = psd_bulk_list$Cytoplasm_r[[paste0(paste0(labeling_time,"h_labeling_"),wash_vec[i],"h_wash")]][c_index]
                            n_r_vec[i] = psd_bulk_list$Nucleus_r[[paste0(paste0(labeling_time,"h_labeling_"),wash_vec[i],"h_wash")]][n_index]
                            x_r_vec[i] = psd_bulk_list$RNA_r[[paste0(paste0(labeling_time,"h_labeling_"),wash_vec[i],"h_wash")]][r_index]
                            vc[i] = v_c[[paste0(paste0(labeling_time,"h_labeling_"),wash_vec[i],"h_wash")]]
                            vn[i] = v_n[[paste0(paste0(labeling_time,"h_labeling_"),wash_vec[i],"h_wash")]]
                          }
                          
                        }
                        
                        tmp_obj <- suppressWarnings(tryCatch({
                          nls(c_vec~x*exp(- beta_c * t_vec) - lambda * n0 * vn/(vc * (ked - beta_c)) * exp(-ked * t_vec),
                              start = list(x = 0.1, beta_c = 0.1, lambda = 0.1), lower = list(x = 0, beta_c = 0, lambda = 0), algorithm = "port")
                          #nls(c_vec~x*exp(- beta_c * t_vec) - (ked + beta_c * c_r_vec/n_r_vec - x_r_vec/n_r_vec * b) * n0 * vn/(vc * (ked - beta_c)) * exp(-ked * t_vec),
                          #start = list(x = 0.1, beta_c = k[1]*ked[1]/2))
                        },
                        error = function(x){return(NA) }))
                        
                        tmp_list <- list()
                        if(!sum(is.na(tmp_obj))){
                          tmp_list[["lambda"]] = coef(tmp_obj)[["lambda"]]
                          tmp_list[["beta_c"]] = coef(tmp_obj)[["beta_c"]]
                        }else {tmp_list[["lambda"]]=NA;tmp_list[["beta_c"]] =NA}
                        
                        return(tmp_list)
                      }
                    })
  
  tmp_df2 <- do.call(rbind,tmp_list) %>% as.data.frame()
  tmp_df[["lambda"]] <- as.numeric(tmp_df2$lambda)
  #tmp_df[["lambda2"]] <- as.numeric(tmp_df2$lambda2)
  tmp_df[["beta_c"]] <- as.numeric(tmp_df2$beta_c)
  #tmp_df[["beta_c2"]] <- as.numeric(tmp_df2$beta_c2)
  tmp_df[["beta_n"]] <- apply(tmp_df,1,
                              FUN = function(x){
                                ked <- as.numeric(x[["ked"]]); ke <- as.numeric(x[["lambda"]])
                                ked - ke
                              })
  #tmp_df[["lambda_vec"]] <- tmp_df2$lambda_vec
  #res_df <- subset(tmp_df, !is.na(kdn) & kdn > 0 & kdc > 0 & ke > 0 )
  res_df <- subset(tmp_df, !is.na(alpha) & alpha > 0)
}



get_psd_bulk <- function(cc_stage = NULL, kd_label = NULL, seurat.obj){
  psd_bulk_list <- list(RNA = data.frame(gene = row.names(starmap_normalized)),
                        Cytoplasm = data.frame(gene = row.names(starmap_normalized)),
                        Nucleus = data.frame(gene = row.names(starmap_normalized)))
  
  for(i in levels(starmap_normalized@meta.data$sample)){
    tmp_vec <- starmap_normalized@meta.data$sample == i 
    if(!is.null(cc_stage)) tmp_vec & starmap_normalized@meta.data$Phase == cc_stage
    if(!is.null(kd_label)) tmp_vec & starmap_normalized@meta.data$KD_label_combined == kd_label
    psd_bulk_list$RNA[[i]] <- rowSums(starmap_normalized@assays$RNA@counts[,tmp_vec])/sum(tmp_vec)
    psd_bulk_list$Cytoplasm[[i]] <- rowSums(starmap_normalized@assays$cytoplasm@data[,tmp_vec])/sum(tmp_vec)
    psd_bulk_list$Nucleus[[i]] <- rowSums(starmap_normalized@assays$nucleus@data[,tmp_vec])/sum(tmp_vec)
  }
}

get_density <- function(x, y, ...) {
  suppressMessages(library(MASS))
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

samplewise_norm <- function(starmap, #with volume n_volume info
                            raw_norm = F, #just norm by anchor genes
                            anchor_gene = NULL, #anchor gene names in assay
                            anchor_vec = NULL, #norm vec from meta.data etc.
                            ref_sample = "2h_labeling"
                            ){
  #Normalize by ctrl genes 
  starmap_normalized <- starmap
  if(sum(starmap_normalized@meta.data$volume < starmap_normalized@meta.data$n_volume) > 0){
    warning("Found some cells with larger nucleus volume than cell volume")
  }
  starmap_normalized <- starmap_normalized[,starmap_normalized@meta.data$volume > starmap_normalized@meta.data$n_volume]
  normalize_list <- list(RNA = rep(1,length(levels(as.factor(starmap@meta.data$sample)))))
  names(normalize_list[["RNA"]]) <- levels(as.factor(starmap@meta.data$sample))
  if(!is.null(anchor_vec)){
    for(i in levels(as.factor(starmap@meta.data$sample))){
      normalize_list[["RNA"]][i] <- mean(anchor_vec[starmap@meta.data$sample == i])
      normalize_list[["RNA"]][i] <- normalize_list[["RNA"]][i]/mean(anchor_vec[starmap@meta.data$sample == ref_sample])
    }
  }else if(!is.null(anchor_gene)){
    for(i in levels(as.factor(starmap@meta.data$sample))){
      normalize_list[["RNA"]][i] <- mean(colSums(starmap@assays$RNA@counts[anchor_gene,starmap@meta.data$sample == i]))
      normalize_list[["RNA"]][i] <- normalize_list[["RNA"]][i]/
        mean(colSums(starmap@assays$RNA@counts[anchor_gene,starmap@meta.data$sample == ref_sample]))
    }
  }else warning("No anchor genes specified")
  
  ##Original/1000gene Data
  #for(i in levels(starmap@meta.data$sample)){
  #  normalize_list[["RNA"]][i] <- mean(colSums(starmap@assays$RNA@counts[c(1:2,4:7),starmap@meta.data$sample == i]))
  #  normalize_list[["RNA"]][i] <- normalize_list[["RNA"]][i]/mean(colSums(starmap@assays$RNA@counts[c(1:2,4:7),starmap@meta.data$sample == "1h_labeling"]))
  #}
  ##16 gene
  #for(i in levels(starmap@meta.data$sample)){
  #  normalize_list[["RNA"]][i] <- mean(starmap@assays$RNA@counts[16,starmap@meta.data$sample == i])
  #  normalize_list[["RNA"]][i] <- normalize_list[["RNA"]][i]/mean(starmap@assays$RNA@counts[16,starmap@meta.data$sample == "1h_labeling"])
  #}
  ##cardiac/skin
  #for(i in levels(as.factor(starmap@meta.data$sample))){
  #  normalize_list[["RNA"]][i] <- mean(starmap@meta.data$AF546[starmap@meta.data$sample == i])
  #  normalize_list[["RNA"]][i] <- normalize_list[["RNA"]][i]/mean(starmap@meta.data$AF546[starmap@meta.data$sample == "2h_labeling"])
  #}
  
  for(i in levels(as.factor(starmap@meta.data$sample))){
    starmap_normalized@assays$RNA@counts[,starmap_normalized@meta.data$sample == i] <- 
      starmap_normalized@assays$RNA@counts[,starmap_normalized@meta.data$sample == i]/normalize_list[["RNA"]][i]
    starmap_normalized@assays$nucleus@data[,starmap_normalized@meta.data$sample == i] <- 
      starmap_normalized@assays$nucleus@data[,starmap_normalized@meta.data$sample == i]/normalize_list[["RNA"]][i]
    starmap_normalized@assays$cytoplasm@data[,starmap_normalized@meta.data$sample == i] <- 
      starmap_normalized@assays$cytoplasm@data[,starmap_normalized@meta.data$sample == i]/normalize_list[["RNA"]][i]
  }
  
  #Normalize by volume
  if(raw_norm == T){
    starmap_normalized
  }else{
    tmp_mat <- matrix(rep(starmap_normalized@meta.data$volume/100000,each = dim(starmap_normalized@assays$RNA@counts)[1]), byrow = F,
                      nrow = dim(starmap_normalized@assays$RNA@counts)[1])
    starmap_normalized@assays$RNA@counts <- starmap_normalized@assays$RNA@counts / tmp_mat
    
    tmp_mat <- matrix(rep(starmap_normalized@meta.data$n_volume/100000,each = dim(starmap_normalized@assays$nucleus@data)[1]), byrow = F,
                      nrow = dim(starmap_normalized@assays$nucleus@data)[1])
    starmap_normalized@assays$nucleus@data <- starmap_normalized@assays$nucleus@data / tmp_mat
    
    tmp_mat <- matrix(rep((starmap_normalized@meta.data$volume-starmap_normalized@meta.data$n_volume)/100000,
                          each = dim(starmap_normalized@assays$cytoplasm@data)[1]), byrow = F,
                      nrow = dim(starmap_normalized@assays$cytoplasm@data)[1])
    starmap_normalized@assays$cytoplasm@data <- starmap_normalized@assays$cytoplasm@data / tmp_mat
    
    rm(tmp_mat)
    starmap_normalized
  }

}

###help function
sc <- function(seq1,seq2){
  strsplit(seq1,"")[[1]] == strsplit(seq2,"")[[1]]
}
###指标
AAR <- function(gen_seq, origin_seq){
  sum(sc(gen_seq, origin_seq))/nchar(origin_seq)
}
two_AAR <- function(seq1,seq2,cores){
  seq1 <- seq1
  seq2 <- seq2
  cl <- parallel::makeCluster(cores)
  parallel::clusterExport(cl, c("seq1","AAR","seq2","sc"),envir = environment())
  res <- pbapply::pblapply(1:length(seq1),
                           function(x){
                             AAR(seq1[x],seq2[x])
                           },cl = cl)
  res <- unlist(res)
  return(res)
}
pair_AAR <- function(seq_vec){
  split_seq <- lapply(seq_vec, strsplit, "")
  # Convert to matrix
  seq_matrix <- do.call(rbind, lapply(split_seq, unlist))
  n <- nrow(seq_matrix)
  mat <- matrix(0, nrow = n, ncol = n)
  mat <- outer(1:n, 1:n, Vectorize(function(i, j) mean(seq_matrix[i, ] == seq_matrix[j, ])))
  res <- mat[lower.tri(mat,diag=FALSE)]
  return(res)
}
pair_seq_sim <- function(seq_list, cores){
  res <- protr::parSeqSim(seq_list, cores = cores)
  return(res)
}
pair_seq_dis <- function(seq_list, cores, method){
  ##method 主要四个 lv, lcs, hamming, dl
  res <- stringdist::stringdistmatrix(unlist(seq_list),unlist(seq_list),
                                      nthread = cores, method = method)
  colnames(res) <- names(seq_list)
  rownames(res) <- names(seq_list)
  return(res)
}
cal_AAR_dis <- function(gen_dt,gen_cdr3_col,origin_cdr3_col,seq_id_col="par_seqid",gen_id_col="id",need_AAR=TRUE,need_dis=TRUE){
  if (need_AAR){
    ###生成序列和真实序列 AAR 的中位数分布
    message("=======Calculate AAR between ", gen_cdr3_col, " and ",origin_cdr3_col,"========")
    gen_dt$AAR <- two_AAR(seq1 = gen_dt[,gen_cdr3_col],
                          seq2 = gen_dt[,origin_cdr3_col],cores=60)
    OAAR <- gen_dt %>% 
      group_by(get(seq_id_col)) %>% 
      summarise(median_AAR = median(AAR, na.rm = T),
                IQR_AAR = IQR(AAR, na.rm = T),
                mean_AAR = mean(AAR, na.rm = T),
                sd_AAR = sd(AAR, na.rm = T),
                AAR_F = mean(AAR < 0.1, na.rm = T)) %>% ungroup()
    ###生成序列之间的 AAR 中位数分布
    message("=======Calculate AAR between pairwise sequence of ",gen_cdr3_col,"=========")
    all_ids <- unique(gen_dt[,seq_id_col])
    IAAR <- vector("list",length = length(all_ids))
    for (i in seq_along(IAAR)){
      tt <- gen_dt %>% filter(get(seq_id_col) == all_ids[i])
      res <- pair_AAR(tt[,gen_cdr3_col])
      IAAR[[i]] <- data.frame(seqid = all_ids[i], 
                              median_AAR = median(res, na.rm = T),
                              IQR_AAR = IQR(res, na.rm = T),
                              mean_AAR = mean(res, na.rm = T),
                              sd_AAR = sd(res, na.rm = T),
                              AAR_F = mean(res < 0.1, na.rm = T))
      message("Complete ",i)
    }
    IAAR <- bind_rows(IAAR)
  }
  if (need_dis){
    ##编辑距离
    message("=========Calculate Edit Distance========")
    all_ids <- unique(gen_dt[,seq_id_col])
    gen_dt$dis_lv <- NA
    dis_lv <- vector("list",length = length(all_ids))
    for (i in seq_along(dis_lv)){
      tt <- gen_dt %>% filter(get(seq_id_col) == all_ids[i])
      seq_list <- as.list(c(tt[,gen_cdr3_col],unique(tt[,origin_cdr3_col])))
      names(seq_list) <- c(tt[,gen_id_col],"origin_seq")
      res_lv <- pair_seq_dis(seq_list,cores = 30,method="lv")
      gen_dt$dis_lv[which(gen_dt[,seq_id_col] == all_ids[i])] <- res_lv[1:(nrow(res_lv)-1),nrow(res_lv)]
      
      res_lv <- res_lv[1:(nrow(res_lv)-1),1:(nrow(res_lv)-1)]##内部相似性
      dis_lv[[i]] <- data.frame(seqid = all_ids[i], 
                                median_dis = median(res_lv[lower.tri(res_lv,diag=FALSE)], na.rm = T),
                                IQR_dis = IQR(res_lv[lower.tri(res_lv,diag=FALSE)], na.rm = T),
                                mean_dis = mean(res_lv[lower.tri(res_lv,diag=FALSE)], na.rm = T),
                                sd_dis = sd(res_lv[lower.tri(res_lv,diag=FALSE)], na.rm = T))
      message("Complete ",i)
    }
    
    odis_lv <- gen_dt %>% 
      group_by(get(seq_id_col)) %>% 
      summarise(median_dis = median(dis_lv, na.rm = T),
                IQR_dis = IQR(dis_lv, na.rm = T),
                mean_dis = mean(dis_lv, na.rm = T),
                sd_dis = sd(dis_lv, na.rm = T)) %>% ungroup()
    idis_lv <- bind_rows(dis_lv)
  }
  message("======Done=======")
  if (need_AAR & need_dis){
    return(list(IAAR = IAAR, 
                OAAR = OAAR, 
                ODIS_lv = odis_lv,
                IDIS_lv = idis_lv))
  }
  if (need_AAR & !need_dis){
    return(list(IAAR = IAAR, OAAR = OAAR))
  }
  if (!need_AAR & need_dis){
    return(list(ODIS_lv = odis_lv,
                IDIS_lv = idis_lv))
  }
}

res <- cal_AAR_dis(as.data.frame(ReprogBERT_gen), 
                   gen_cdr3_col = "gen_CDR3", origin_cdr3_col = "HCDR3",
                   seq_id_col = "seqid", gen_id_col= "id")
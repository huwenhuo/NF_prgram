color_list = list(up = '#D6604D', down = '#4393C3', neutral = '#F7F7F7')

aml_genes = c('TET2',  'DNMT3A',  'SF3B1',  'TP53',  'PPM1D',  'ASXL1',  'NRAS', 
  'JAK2',  'GNB1',  'KRAS',  'CBL',  'IDH1',  'IDH2') 
aml_fuse = c('MLL_re', 'FLT3ITD', 'AML1ETO', 'PMLRARA', 'CBFBMYH11', 'GATA2MECOM')
aml_mut = c(paste0(aml_genes, '_mut'), aml_fuse)

scale_fun = function(x, scale = 'row', na.rm = T){
    if (scale == "row") {
        x <- sweep(x, 1L, rowMeans(x, na.rm = na.rm), check.margin = FALSE)
        sx <- apply(x, 1L, sd, na.rm = na.rm)
        x <- sweep(x, 1L, sx, "/", check.margin = FALSE)
    }
    else if (scale == "column") {
        x <- sweep(x, 2L, colMeans(x, na.rm = na.rm), check.margin = FALSE)
        sx <- apply(x, 2L, sd, na.rm = na.rm)
        x <- sweep(x, 2L, sx, "/", check.margin = FALSE)
    }
}

pathway_fun = function(){ 
    setwd('/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Collaboration_Projects/for_Min/RNA-Seq_Zheng')
    anno=read.table('/research_jude/rgs01_jude/groups/jxugrp/home/common/Data_Genomics/RNA-seq/Erythroid_Metabolome/CPDB_gene2path.txt',header=T,sep="\t")
    uni_path=as.vector(unique(anno[,2]))
    # use CPDB database to do pathway enrichment.
    ## do path enrichment using CPDB database.
    norm_readcount=read.table('readcount_TPM_FPKM/normalized_readcount_CI_CD_HI_HD_NI_ND.txt',header=T,sep="\t",row.names=1)
    metab_gene=read.table('/research_jude/rgs01_jude/groups/jxugrp/home/common/Data_Genomics/RNA-seq/Erythroid_Metabolome/Metabolic_Gene_List/Sabatini_Metabolic_Gene_List_Human.txt',header=T,sep="\t")[,1,drop=FALSE]
    name=c('maSigPro_CD_CI_NI','maSigPro_HD_HI_NI')
    for (i in 1:length(name)){
        sig_gene=read.table(paste0(name[i],'/maSigPro_sig_gene_2cluster.txt'),header=T,sep="\t")
        sig_gene=merge(sig_gene,metab_gene,by.x=0,by.y='Gene.Symbol')
        all_gene=rownames(norm_readcount)
        N=length(all_gene)
        enrich_path=matrix(NA,length(uni_path),7)
        enrich_path[,1]=uni_path
        for (j in 1:2){
            k=NULL; n=NULL; M=NULL;
            for (w in 1:dim(enrich_path)[1]){
                gene_in_path=intersect(unique(as.vector(anno[anno[,2]==enrich_path[w,1],1])),all_gene)
                M[w]=length(gene_in_path)
                gene_in_cluster=sig_gene[sig_gene[,2]==j,1]
                n[w]=length(gene_in_cluster)
                k[w]=length(intersect(gene_in_path,gene_in_cluster))
            }
            p=1-phyper(k-1,M,N-M,n)
            ratio=k/M
            enrich_path[,2:7]=cbind(N,M,n,k,p,ratio)
            enrich_path=enrich_path[order(p),]
            colnames(enrich_path)=c('Path_class','#Gene_NormRead_over10','#Gene_in_Path','#Gene_in_Cluster','#Overlap','P_value','Percentage')
            write.table(enrich_path,paste0(name[i],'/maSigPro_group',j,'_enrich_CPDB_path.txt'),row.names=F,col.names=T,sep="\t",quote=F)
        }
    }
}

add_ann = function(dat, n = NULL){
    if(!is.data.table(dat)){dat = as.data.table(dat, keep.rownames = T)}
    if(is.null(n)){n = length(unlist(strsplit(unname(unlist(dat[1,1])), ':')))}
    dat[, class_id := unlist(strsplit(rn, ':'))[n], by = 1:nrow(dat) ]
    dat[, family_id := unlist(strsplit(rn, ':'))[n-1], by = 1:nrow(dat) ]
    dat[, gene_id := unlist(strsplit(rn, ':'))[n-2], by = 1:nrow(dat) ]
    if(n == 4) { dat[, tx_id := unlist(strsplit(rn, ':'))[n-3], by = 1:nrow(dat) ] }
    dat
}

deg_local_fun = function(ctype1 = 'MEP', ctype2 = 'EB', save = F){
    group1_name = paste(ctype1, collapse = '_')
    group2_name = paste(ctype2, collapse = '_')
    dsn_sel = dsn_telocal[dsn_telocal$cell_type %in% c(ctype1, ctype2), , drop = F]
    dsn_sel$group[dsn_sel$cell_type %in% ctype1] = group1_name
    dsn_sel$group[dsn_sel$cell_type %in% ctype2] = group2_name
    dsn_sel$group = factor(dsn_sel$group, levels = c(group1_name, group2_name))

    mtx_sel = telocal_dt_w_sel_mtx[, rownames(dsn_sel) ]
    row_names = rownames(mtx_sel)
    class_ids = unname(unlist(sapply(row_names, function(x) unlist(strsplit(x, ':'))[4])) )
    mtx_sel = mtx_sel[class_ids %in% TE_class_IDs, ]

    dge = DGEList(counts=mtx_sel, group = dsn_sel$group)
    dge <- calcNormFactors(dge)
    #keep = filterByExpr(dge)
    #dge = dge[keep, , keep.lib.sizes=FALSE]
    dge = normLibSizes(dge)
    design = model.matrix( ~ 0 + dsn_sel$group)
    dge = estimateDisp(dge, design)
    et <- exactTest(dge)
    deg = topTags(et, n = nrow(mtx_sel))
    deg = as.data.table(deg$table, keep.rownames = T)
    deg

}


deg_telocal_fun = function(ctype1 = 'MEP', ctype2 = 'EB', mtx, dsn, save = F){
    group1_name = paste(ctype1, collapse = '_')
    group2_name = paste(ctype2, collapse = '_')
    dsn_sel = dsn_telocal[dsn_telocal$cell_type %in% c(ctype1, ctype2), , drop = F]
    dsn_sel$group[dsn_sel$cell_type %in% ctype1] = group1_name
    dsn_sel$group[dsn_sel$cell_type %in% ctype2] = group2_name
    dsn_sel$group = factor(dsn_sel$group, levels = c(group1_name, group2_name))

    mtx_sel = telocal_dt_w_sel_mtx[, rownames(dsn_sel) ]
    row_names = rownames(mtx_sel)
    class_ids = unname(unlist(sapply(row_names, function(x) unlist(strsplit(x, ':'))[4])) )
    mtx_sel = mtx_sel[class_ids %in% TE_class_IDs, ]
    
    dge = DGEList(counts=mtx_sel, group = dsn_sel$group)
    dge <- calcNormFactors(dge)
    #keep = filterByExpr(dge)
    #dge = dge[keep, , keep.lib.sizes=FALSE]
    dge = normLibSizes(dge)
    design = model.matrix( ~ 0 + dsn_sel$group)
    dge = estimateDisp(dge, design)
    et <- exactTest(dge)
    deg = topTags(et, n = nrow(mtx_sel))
    deg = as.data.table(deg$table, keep.rownames = T)
    deg

}

deg_fun = function(group1, group2, dsn, mtx, if_cpm = F){
    group1_name = paste(group1, collapse = '_')
    group2_name = paste(group2, collapse = '_')
    dsn_sel = dsn[dsn$group %in% c(group1, group2), , drop = F]
    dsn_sel$grp[dsn_sel$group %in% group1] = group1_name
    dsn_sel$grp[dsn_sel$group %in% group2] = group2_name
    dsn_sel$grp = factor(dsn_sel$grp, levels = c(group1_name, group2_name))

    if(is.data.table(dsn_sel) == T){
        mtx_sel = mtx[, dsn_sel$sample_name]
    }else{
        mtx_sel = mtx[, rownames(dsn_sel) ]
    }
    
    mtx_sel = mtx_sel[rowSums(mtx_sel) > ncol(mtx_sel), ]

    dge = DGEList(counts=mtx_sel, group = dsn_sel$grp)
    dge = calcNormFactors(dge)
    dge = normLibSizes(dge)
    design = model.matrix( ~ 0 + dsn_sel$grp)
    dge = estimateDisp(dge, design)
    et  = exactTest(dge)
    deg = topTags(et, n = nrow(mtx_sel))
    deg = as.data.table(deg$table, keep.rownames = T)
    if(if_cpm == T){
        fname = paste0('logcpm_', paste(group1, collapse = '_'), '__', paste(group2, collapse = '__'), '.tsv')
        write.table(cpm(dge, normalize = T, lo), file = fname, log = T)
    }
    deg
}

cut_quantile = function(x, n = NULL, probs = NULL, n_names = NULL){
    if(is.null(probs) & !is.null(n)){ probs = c(0, 1:(n-1) / n, 1 ) }
    if(is.null(n)){n = length(probs) - 1}
    if(is.null(n_names)){
        if(n == 2) { n_names = c('low', 'high') }
        if(n == 3) { n_names = c('low', 'intermediate', 'high') }
        if(n == 4) { n_names = c('Q1', 'Q2', 'Q3', 'Q4') }
    }
    grps = cut(x, quantile(x, prob = probs), include = TRUE)
    lvls = sort(unique(grps))
    names(n_names) = lvls
    grps = unname(n_names[grps])
    grps = factor(grps, levels = n_names)
    grps
}

head2 = function(x, ii_row = 3, ii_col = 5){
    head(x[1:ii_row, 1:ii_col])
}
tail2 = function(x, ii_row = 3, ii_col = 5){
    tail(x[(nrow(x)-ii_row):nrow(x), (ncol(x)-ii_col):ncol(x)])
}

cpm_fun = function(mtx, dsn = NULL){
    if(is.null(dsn)){
        dge = DGEList(counts=mtx)
    }else{
        dge = DGEList(counts=mtx, group = dsn$group)
    }
    dge = calcNormFactors(dge)
    dge = normLibSizes(dge) 
    logcpm = cpm(dge, log = T)
    logcpm
}

heatmap_fun = function(mtx, dsn, deg, if_return = F, ...){
    plotdat = mtx[deg[FDR < 0.05, rn], rownames(dsn)]
    plotdat = scale_fun(plotdat)
    plotdat[plotdat > 2] = 2
    plotdat[plotdat < -2] = -2

    column_ha = HeatmapAnnotation(df = dsn[, 'group', drop = F])
    args = list(matrix = plotdat, top_annotation = column_ha, column_split = factor(dsn$group), show_column_names = F, show_row_names = T, show_column_dend = F,
       cluster_rows = F,  row_names_gp = grid::gpar(fontsize = 7),  show_row_dend = F)
    args_other = list(...)
    for(ii in 1:length(args_other)){
        args[[names(args_other)[ii]]] = args_other[[ii]]
    }
    hp = do.call(Heatmap, args)
    draw(hp)
    if(if_return == T){
        hp
    }
}

deg_fun2 = function (dsn, mtx, if_cpm = F, if_filter = F){
    groups = dsn$group
    dge = DGEList(counts = mtx, group = groups) 
    if(if_filter == T) {
        keep = filterByExpr(dge)
        dge = dge[keep,] }
    dge = calcNormFactors(dge)
    dge = normLibSizes(dge)
    design = model.matrix(~0 + groups)    
    dge = estimateDisp(dge, design)
    et = exactTest(dge) 
    deg = topTags(et, n = nrow(mtx))
    deg = as.data.table(deg$table, keep.rownames = T)
    deg = deg[order(PValue), ]
    if (if_cpm == T) { 
        logcpm = cpm(dge, log = T)      
        deg = list(deg = deg, logcpm = logcpm)
    }           
    deg         
}               

fgsea_fun = function(mtx = NULL, ranking = NULL, pathways = NULL, if_all_pathways = F, species = 'human'){
    if(is.null(ranking) & !is.null(mtx)){ ranking = ranking_fun(mtx) }

    if(is.null(pathways)){
        if(if_all_pathways == T){
            pathwaysDF = msigdbr(species)
            pathways <- split(as.character(pathwaysDF$gene_symbol), pathwaysDF$gs_name)
        }else{
            pathwaysDF <- msigdbr(species, category="H")
            pathways <- split(as.character(pathwaysDF$gene_symbol), pathwaysDF$gs_name)
        }
    }

    res = fgsea(pathways = pathways,  stats = ranking, scoreType = 'std',  minSize = 10, maxSize = 500, nproc = 10) 
    res = res[order(padj), ]
    res

}

ranking_fun = function(mtx){
    mtx[, rank := - sign(logFC) * log(FDR, 10)] 
    mtx = mtx[order(rank), ]
    mtx[rank == -Inf, rank := min(mtx[is.finite(rank), rank]) - 1:nrow(mtx[rank == -Inf,]) ]
    mtx[rank == Inf,  rank := max(mtx[is.finite(rank), rank]) + 1:nrow(mtx[rank == Inf,]) ]
    rankings = mtx$rank
    names(rankings) = mtx$rn
    rankings
}

pca_fun = function (mtx, ntop = 1000, if_plot = F, dsn = NULL, col = 'blue') {
    row_sd = rowSds(mtx)
    row_sd = sort(row_sd, decreasing = T)
    sel = names(row_sd[1:ntop]) 
    mtx = mtx[sel, ]
    pc = prcomp(t(mtx), center = TRUE, scale = TRUE)    
    lblx = paste0("PC1 ", round(100 * pc$sdev[1]^2/sum(pc$sdev^2), 1), "% explained")
    lbly = paste0("PC2 ", round(100 * pc$sdev[2]^2/sum(pc$sdev^2), 1), "% explained")
    plotdat = as.data.table(pc$x, keep.rownames = T)
    if(!is.null(dsn)){ plotdat = merge(plotdat, dsn, by.x = 'rn', by.y = 'sample_name', all.x = T) }
    if(if_plot == T){
        gg = ggscatter(plotdat, x = 'PC1', y = 'PC2', col = col) + xlab(lblx) + ylab(lbly)
        list(plotdat = plotdat, labx = lblx, laby = lbly, gg = gg)
    }else{
        list(plotdat = plotdat, labx = lblx, laby = lbly)
    }
}

logcpm_fun = function(mtx, if_plus1 = F){
    dge = DGEList(counts = mtx)
    dge = calcNormFactors(dge)
    dge = normLibSizes(dge)
    if(if_plus1 == T){
        logcpm = cpm(dge, log=F, normalized.lib.sizes = T)
        logcpm = log2(logcpm + 1)
    }else{
        logcpm = cpm(dge, log=T, normalized.lib.sizes = T)
    }
    logcpm
}

share_dir = '/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/WenhuoHu/share_data/'

create_named_list <- function(...) {
  args <- list(...)
  names(args) <- sapply(as.list(substitute(list(...)))[-1L], deparse)
  return(args)
}

genome_fun = function(genome = NULL){
    hg38 = list(bowtie2_index = '/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Data/UCSC/hg38/Bowtie2Index')
    GRCh38 = list( name = 'GRCh38',
                  genome_fasta   = '/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/WenhuoHu/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa',
                  bismark_index  = '/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/WenhuoHu/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/Bismark/', 
                  bowtie2_index  = '/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/WenhuoHu/genomes/GRCh38/bowtie2_index/bowtie2_index',
                  bwa5_index     = '/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/WenhuoHu/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/BWAIndex/version0.5.x/', 
                  bwa6_index     = '/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/WenhuoHu/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/BWAIndex/version0.6.0/', 
                  gtf_file       = '/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/WenhuoHu/igenomes/Homo_sapiens/NCBI/GRCh38/Annotation/Genes/genes.gtf',
                  chrom_size     = '/home/whu78/WenhuoHu/genomes/GRCh38/grch38_chrosomesizes', 
                  blacklist_file   = '/home/whu78/repos/atacseq/assets/blacklists/v3.0/GRCh38-blacklist.v3.bed') 
    hg19 = list( ref_dir = '/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/YuannyuZhang/Data/For_RNA-Seq/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index') 
    mm10 =  list( name = 'mm10',
                 genome_fastq  = '/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/WenhuoHu/igenomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa',
                 star_index    = '/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/WenhuoHu/genomes/mm10/star_index/',
                 gtf_file      = '/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/WenhuoHu/igenomes/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf',
                 bowtie2_index = '/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/WenhuoHu/igenomes/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/', 
                 bwa5_index    = '/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/WenhuoHu/igenomes/Mus_musculus/NCBI/mm10/Sequence/BWAIndex/version0.5.x/',
                 bwa6_index    = '/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/WenhuoHu/igenomes/Mus_musculus/NCBI/mm10/Sequence/BWAIndex/version0.6.0/',
                 blacklist_file   = '/home/whu78/WenhuoHu/repos/Blacklist/lists/mm10-blacklist.v2.bed.gz') 
    t2t =  list( name = 't2t',
                rsem_index    = '/home/whu78/WenhuoHu/genomes/t2t/rsem/',
                genome_fasta  = '/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/WenhuoHu/genomes/t2t/hs1.fa', 
                gtf_file      = '/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/WenhuoHu/genomes/t2t/hs1.ncbiRefSeq.gtf', 
                bowtie2_index = '/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/WenhuoHu/genomes/t2t/bowtie2_index/', 
                star_index    = '/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/WenhuoHu/genomes/t2t/star_indx/', 
                cellranger_index = '/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/WenhuoHu/genomes/t2t/t2t_cellranger/', 
                bismark_index = '/research_jude/rgs01_jude/groups/jxugrp/home/common/Lab_Members/WenhuoHu/genomes/t2t/Bisulfite_Genome/')
    genome_list = list(hg38 = hg38, GRCh38 = GRCh38, hg19 = hg19, mm10 = mm10, t2t = t2t)
    if(genome %in% names(genome_list)){
        genome_list[[genome]]
    }else{
        genome_list
    }
}
                                     
sif = list( 
    bedtool = paste0('/home/whu78/images/depot.galaxyproject.org-singularity-bedtools-2.30.0--hc088bd4_0.img'), 
    bowtie2 = paste0('/home/whu78/images/bowtie2_a2.5.4.sif'),
    velo = paste0('/home/whu78/images/velocyto_a1.17.17.sif'),
    bismark = paste0('/home/whu78/images/bismark_a0.24.2.sif'),
    samtools= paste0('/home/whu78/images/depot.galaxyproject.org-singularity-samtools-1.20--h50ea8bc_0.img'),
    sratool = paste0('/home/whu78/images/sratoolkit_latest.sif'), 
    deeptool = paste0('/home/whu78/images/depot.galaxyproject.org-singularity-deeptools-3.5.5--pyhdfd78af_0.img')
)

surv_fun = function(data, vars, n_cut = 2, n_names = NULL){
    surv_list = lapply(vars, function(var){
        plotdat = data[, c(var, 'OS_MONTHS', 'OS_STATUS'), with = F]
        if(is.null(n_names)){n_names = paste0('Q', 1:n_cut)}
        plotdat$lvl = cut_quantile(unlist(data[, var, with = F]), n = n_cut, n_names = n_names)

        # Extract p-value
        test <- survdiff(Surv(OS_MONTHS, OS_STATUS) ~ lvl, data = plotdat)
        p_value <- 1 - pchisq(test$chisq, length(test$n) - 1)

        # Summary of the survfit object
        fit <- survfit(Surv(OS_MONTHS, OS_STATUS) ~ lvl, data = plotdat)
        summary_fit <- summary(fit)
        tmp = as.data.table(summary_fit$table, keep.rownames = T)
        tmp$pval = p_value
        tmp$com = var
        tmp
    }) 
    surv_dt = rbindlist(surv_list)
    surv_dt = surv_dt[order(pval), ]
    surv_dt
}
                                     
surv_sum_fun = function(data, group = 'lvl', os_time = 'OS_MONTHS', os_status = 'OS_STATUS'){
    # Extract p-value
    fm = as.formula(paste0("Surv(", os_time, ', ', os_status, ') ~', group))
    test <- survdiff(fm, data = data)
    p_value <- 1 - pchisq(test$chisq, length(test$n) - 1)
    # Summary of the survfit object
    fit <- survfit(fm, data = data)
    summary_fit <- summary(fit)
    tmp = as.data.table(summary_fit$table, keep.rownames = T)
    tmp$pval = p_value
    tmp
}

size_fun = function(hh = 4, ww = 4, res = 300){
    options(repr.plot.height = hh, repr.plot.width = ww, repr.plot.res = res)
}

size = size_fun

combn_fun = function(groups){
    comb = combn(groups, 2, simplify = F)
    comb = as.data.table(t(as.data.table(comb)), keep.rownames = T)
    setnames(comb, 'V1', 'grp1')
    setnames(comb, 'V2', 'grp2')
    comb[, tag := paste0(grp1, '__', grp2) ]
    comb$rn = 1:nrow(comb)
    (comb)
}

# Function to find the latest version file, read, modify, and save with the new version
read_file <- function(file_prefix, path = '~/WenhuoHu/collab_Farhan/') {
    # Step 1: List all .rds files with the prefix
    files <- list.files(pattern = paste0("^", file_prefix, "_v\\d+\\.rds$"), path = path)

    # If no files exist, return an error
    if (length(files) == 0) {
        cat("No RDS files with the given prefix found.\n")
        return(NULL)
    }

    # Step 2: Extract version numbers from the filenames
    version_numbers <- gsub(paste0("^", file_prefix, "_v(\\d+)\\.rds$"), "\\1", files)
    version_numbers <- as.numeric(version_numbers)

    # Find the maximum version number
    max_version <- max(version_numbers)

    # Step 3: Read the RDS file with the maximum version
    latest_file <- paste0(path, '/', file_prefix, "_v", max_version, ".rds")
    cat("Current maximum version:", max_version, "\t", latest_file, "\n")

    loaded_data <- readRDS(latest_file)
}

save_file = function(val, file_prefix = NULL, path = '~/WenhuoHu/collab_Farhan/'){
    if(is.null(file_prefix)){ file_prefix = substitute(val) }
    # Step 1: List all .rds files with the prefix
    files <- list.files(pattern = paste0("^", file_prefix, "_v\\d+\\.rds$"), path = path)

    # If no files exist, return an error
    if (length(files) == 0) {
        cat("No RDS files with the given prefix found.\n")
        return(NULL)
    }

    # Step 2: Extract version numbers from the filenames
    version_numbers <- gsub(paste0("^", file_prefix, "_v(\\d+)\\.rds$"), "\\1", files)
    version_numbers <- as.numeric(version_numbers)

    # Find the maximum version number
    max_version <- max(version_numbers)

    # Step 3: Read the RDS file with the maximum version
    latest_file <- paste0(path, '/', file_prefix, "_v", max_version, ".rds")

    # Step 4: Increment the version
    new_version <- max_version + 1
    cat("New version:", new_version, "\n")

    # Step 6: Save the updated data with the new version
    new_file <- paste0(path, '/', file_prefix, "_v", new_version, ".rds")
    saveRDS(val, new_file)
    cat("File saved as:", new_file, "\n")
}

add_gene_list = function(val, comment){
    param_name <- substitute(val)
    if (!exists("gene_list", envir = .GlobalEnv)) {
        gene_list = read_file('gene_list')
    } 
    gene_list[[param_name]] = val
    comment(gene_list[[param_name]]) = comment
    save_file(gene_list)
    gene_list
}

srr_id_fun = function(dsn){
	exe_dir = '/cm/shared/apps/python/3.12.x-mamba/bin/'
		dsn[, srr_id := {
			out <- system(paste0('esearch -db sra -query "', gsm_id, '" | efetch -format docsum | grep -o "SRR[0-9]*" | tr "\n" "," '), intern = TRUE)
			if (length(out) == 0) NA_character_ else sub(",$", "", out)
		}, by = 1:nrow(dsn)]

	dsn
}

seurat_process <- function(seu, dims = 1:20, resolution = 0.5, verbose = F) {
  
  if(verbose) {
    msg <- function(x) x  # just run normally
  } else {
    msg <- function(x) suppressMessages(suppressWarnings(x))
  }
  
  # Normalize
  seu <- msg(NormalizeData(seu))
  
  # Identify highly variable features
  seu <- msg(FindVariableFeatures(seu))
  
  # Scale the data
  seu <- msg(ScaleData(seu))
  
  # PCA
  seu <- msg(RunPCA(seu, features = VariableFeatures(seu)))
  
  # UMAP
  seu <- msg(RunUMAP(seu, dims = dims))
  
  # Find neighbors and clusters
  seu <- msg(FindNeighbors(seu, dims = dims))
  seu <- msg(FindClusters(seu, resolution = resolution))
  
  # Return processed Seurat object
  return(seu)
}


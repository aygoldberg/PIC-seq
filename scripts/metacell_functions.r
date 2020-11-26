require("KernSmooth")
require("reshape2")
require("RANN")
require("plyr")
require("plotrix")
require("gplots")
require("parallel")
require("compositions")
require(RCurl)

vecsplit = function(strvec, del, i) {
	 unlist(lapply(sapply(strvec, strsplit, del), "[[", i))
}

image.2 = function(X, col = colorRampPalette(c("blue", "white", "red"))(1000), balance=F, annotate="both", zlim = NULL,
        hct = NULL, vct = NULL, lwd=1, lty=1, cex = 1, text=F, text_mat = X) {
        if (is.null(zlim)) {
                if (balance) {
                                zlim = max(abs(X), na.rm=T) * c(-1,1)
                } else {
                        zlim = quantile(X, c(0,1), na.rm=T)
                }
	}
        hcls = NULL; vcls = NULL; hticks = NULL; vticks = NULL
        nrow = nrow(X); ncol = ncol(X)
        if (!is.null(hct)) {
                X = X[order(hct),]; text_mat = text_mat[order(hct),]
                hticks = seq(-1 / (2 * (nrow-1)),1 + 1 / (2 * (nrow-1)), length.out = nrow + 1)
		hcls = cumsum(table(hct)) + 1
        }
        if (!is.null(vct)) {
                X = X[,order(vct)]; text_mat = text_mat[,order(vct)]
                vticks = seq(-1 / (2 * (ncol-1)),1 + 1 / (2 * (ncol-1)), length.out = ncol + 1)
		vcls = cumsum(table(vct)) + 1
        }
        message("zlim: ", zlim[1], "<>", zlim[2])
        image(t(X), axes = F, col = col, zlim = zlim)
        abline(h = hticks[hcls], v = vticks[vcls], lwd=lwd, lty=lty)
        if (annotate %in% c("rows", "both")) {
	        mtext(rownames(X), las = 2, side=2, at = (1 - seq_len(nrow(X))) / (1 - nrow(X)), cex = cex)
                if (!is.null(hct)) {
                        mtext(names(hcls), side = 4, las = 2, at = rowMeans(cbind(hticks[c(1,hcls[-length(hcls)])], hticks[hcls])), cex = cex)}
	}
	if (annotate %in% c("columns", "both")) {
	        mtext(colnames(X), las = 2, side=1, at = (1 - seq_len(ncol(X))) / (1 - ncol(X)), cex = cex)
                if (!is.null(vct)) {
                                   mtext(names(vcls), side = 3, las = 2, at = rowMeans(cbind(vticks[c(1,vcls[-length(vcls)])], vticks[vcls])), cex = cex)}
        }
	if (text) {
           hmed = seq(0,1,length.out=nrow); vmed = seq(0,1,length.out=ncol)
           text(rep(vmed, each = nrow), rep(hmed, ncol), text_mat)
	}
}

xy_scatter_genes = function(x,y, bad_genes = c(), disp_genes = c(), cols = c("navyblue", "chocolate2"), text = F, fc = 1, reg = 10, lwd = 1) {
	if (length(disp_genes) == 0) {
		z = log2((y + reg) / (x + reg))
		disp_genes = names(which(abs(z) > fc))
	}
        good_genes = setdiff(names(x), bad_genes)
	lim = log2(c(reg, reg + max(c(x,y))))
        plot(log2(x[good_genes] + reg), log2(y[good_genes] + reg), pch = 20, cex = 2, col = cols[1 + good_genes %in% disp_genes],
        xlim = lim, ylim = lim, axes = F, xlab = "", ylab = "")
        if (text & length(disp_genes) > 0) { text(log2(x[disp_genes] + reg), log2(y[disp_genes] + reg), disp_genes)}
        abline(coef = c(fc,1), lty = 2); abline(coef = c(-fc,1), lty = 2)
        axis(1); axis(2)
}

choose_genes_from_clust = function(mc_id, mat_id, good_clusts = colnames(sc_cl@mc_fp), 
	nms_per_clust = 5, nms_thresh = 5, max_num = Inf, bad_genes = c(), must_haves = c(),
	ord = "none") {

	sc_cl = scdb_mc(mc_id); sc_mat = scdb_mat(mat_id)
	lfp = log2(sc_cl@mc_fp[,good_clusts])
	nms = unique(as.vector(apply(lfp,2,function(x){ head(names(sort(x, decreasing = T)),nms_per_clust)})))
	nms = setdiff(nms, c(bad_genes, names(which(apply(lfp[nms,],1,max) < nms_thresh))))
	nms = union(must_haves, head(names(sort(apply(lfp[nms, ], 1, max),T)), max_num - length(must_haves)))
	if (ord == "hc") {
		nms = nms[ hclust(dist(cor(t(lfp[nms,]))), "ward.D2")$order]
	} else if (ord == "max.col") {
		nms = nms[ order(max.col(lfp[nms,]), rowSums(as.matrix(sc_mat@mat[nms,])))]
	}
	nms
}

sc_to_bulk = function(mc_id, mat_id, comb=NULL, bad_genes = c(), cells = names(comb), min_comb = 0, choose_genes = T, normalize = T,
	g1 = NULL, g2 = NULL) {

	if (is.null(comb)) { cells = union(g1,g2)}	
        sc_cl = scdb_mc(mc_id); sc_mat = scdb_mat(mat_id)
	umis = as.matrix(sc_mat@mat[,cells])
	umis_n = sweep(umis,2,colSums(umis),"/") * 1000
	
	if (!normalize) { umis_n = umis}
	if (is.null(comb)) {
		comb = 1 + (cells %in% g1)
		names(comb) = cells
	}
	MAP = as.numeric(factor(comb[cells])); names(MAP) = cells
	if (choose_genes) {
		genes = setdiff(scr_chi_square_diff_genes(umis[,cells], MAP = MAP[cells], fdr = T, pval = 1e-3), bad_genes)
	} else {
		genes = rownames(umis)
	}
	m = t(apply(umis_n[genes, cells],1,tapply,comb[cells],sum))
	sizes = table(comb[cells])
	good_clusts = names(which(sizes >= min_comb))
	m = m[,good_clusts]; sizes = sizes[good_clusts]
	m = sweep(m,2,as.vector(sizes),"/") * min(sizes)
	m
}

create_batch_matrix = function(mc_id, mat_id, comb, fname, batch_iden = "amp_batch_id", cord = colnames(sc_cl@mc_fp), clusts = NULL,
        batch_shades = colorRampPalette(c("white", "navyblue"))(1000), color_batches = F, txt_file = F, norm_by = "col",
        mar = c(3,10,0,20)) {

        sc_cl = scdb_mc(mc_id); sc_mat = scdb_mat(mat_id)
        cells = names(comb)
	cell_stats = sc_mat@cell_metadata[cells,]
        if (is.null(clusts)) { clusts = sc_cl@mc[cells]}
        comb_levels = names(table(comb))
        if (is.null(levels(comb))) { comb = factor(comb); names(comb) = cells}
        batch_meta = unique(cbind(as.vector(cell_stats[cells,batch_iden]), as.vector(comb)))
        B = batch_meta[,1]
        batch_dist = table(cell_stats[cells,batch_iden], clusts[cells])
        batch_dist = batch_dist[B,]
        bct = factor(batch_meta[,2], levels = comb_levels); names(bct) = batch_meta[,1]
	if (color_batches) {
                cols = sc_cl@colors[ as.numeric(cord)]
        } else {cols = rep("black", length(cord))}
        batch_dist = batch_dist[, as.character(cord)]
        if (!is.null(fname)) {
                png(fname, height = max(nrow(batch_dist) * 12, 1000), width = max(ncol(batch_dist) * 17, 1500))
		par(mar = mar)
        }
        if (norm_by == "col") {
                dist_n = sweep(batch_dist,2,colSums(batch_dist),"/")
        } else if (norm_by == "row") {
		dist_n = sweep(batch_dist,1,rowSums(batch_dist),"/")
        } else {
                dist_n = sweep(batch_dist,1,rowSums(batch_dist),"/")
                dist_n = sweep(dist_n,2,colSums(dist_n),"/")
        }
        image.2(dist_n, col = batch_shades, , zlim = c(0, 1), hct = bct, text = F)
	if (!is.null(fname)) {
                dev.off()
        }
	rownames(batch_meta) = batch_meta[,1]
	df = cbind(batch = rownames(batch_dist), source = batch_meta[rownames(batch_dist),2], as.matrix(batch_dist))
	if (txt_file) {
                write.table(df, sep = "\t", quote = F, row.names=F, file = gsub(".png", ".txt", fname))
        }
}

plot_sc_heatmap = function(mc_id, mat_id, nms, clusts = NULL, good_clusts = NULL, cells = NULL, fname = NULL, annotate = F, mar = rep(0,4),
	genes_shades = colorRampPalette(c("white", "orange", "tomato","mediumorchid4", "midnightblue"))(1000), draw_cls = T, normalize = T,
	lwd=1, lty=2) {

	sc_cl = scdb_mc(mc_id); sc_mat = scdb_mat(mat_id)
	if (is.null(good_clusts)) {
		good_clusts = colnames(sc_cl@mc_fp)
	}
	if (is.null(cells)) {
		cells = names(which(sc_cl@mc > 0 & sc_cl@mc %in% good_clusts))
	}
	if (is.null(clusts)) {
		clusts = sc_cl@mc
	}
	umis = as.matrix(sc_mat@mat[,cells])
	umis_n = sweep(umis,2,colSums(umis),"/") * 1000
	if (normalize) {
		foc = log(1 + 7 * umis_n)
	} else {
	        foc = log(1 + 7 * umis)
	}
	cls = cumsum(table(factor(clusts[cells], levels = good_clusts))) / length(cells)	
	if (!is.null(fname)) { 
		png(fname, height = 2000, width = 3300)
		par(mar = mar)
	}
	cell_ord = cells[order(factor(clusts[cells], levels = good_clusts), sample(length(cells)))]
	IM = foc[nms, cell_ord]
	image(t(IM), col = genes_shades, axes = F)
	zlim = quantile(IM, c(0,1))
	message("zlim: ", zlim[1], " - ", zlim[2])
	if (annotate) {
		mtext(nms, side = 2, las = 2, at = (1 - seq_along(nms)) / (1 - length(nms)))
		mtext(names(cls), side = 1, las = 2, at = rowMeans(cbind(c(0,cls[-length(cls)]), cls)))
	}
	if (draw_cls) {
		abline(v = cls, lty = lty, lwd=lwd)
	}
        if (!is.null(fname)) {
		dev.off()
	}
	cell_ord
}

scr_write_models_file = function(mc_id, mat_id, filename = "models.txt") {
	sc_cl = scdb_mc(mc_id); sc_mat = scdb_mat(mat_id)	
	umis = as.matrix(sc_mat@mat[, names(sc_cl@mc)])
        genes = rownames(sc_cl@mc_fp)
        write.table(data.frame(log2(sc_cl@mc_fp), umicount = rowSums(umis[genes,])), col.names = NA, quote = F, sep = "\t", file = filename)
}

scr_chi_square_diff_genes = function(umis, MAP = NULL, g1 = NULL, g2 = NULL, pval, fdr = F) {

  if (is.null(MAP)) {
    MAP = c(rep(1,length(g1)), rep(2, length(g2)))
    names(MAP) = c(g1, g2)
  }
  cells = names(MAP)
  umis = umis[,cells]
  uniform_a = rowSums(umis)/sum(umis)
  exp_count = matrix(uniform_a, ncol = 1) %*% matrix(colSums(umis),1) # exp_counts per cell
  dimnames(exp_count)  = dimnames(umis)
  ex = t(daply(.data= data.frame(cbind(V1 = MAP, t(exp_count)), check.names = F), .(V1), colSums))[-1,]
  obs = t(daply(.data= data.frame(cbind(V1 = MAP, t(umis)), check.names = F), .(V1), colSums))[-1,]

  x2 = rowSums(((obs-ex)^2 )/ex ) # chi^2 with df = ncluster-1

  if (!fdr) {
    sig_genes = x2 > qchisq(1-pval,df= length(unique(MAP)) - 1)
  } else {
    pvals = p.adjust(1 - pchisq(x2, df = length(unique(MAP)) - 1), "fdr")
    sig_genes = pvals < pval
  }
  sig_genes[ is.na(sig_genes)] = F
  return (names(sig_genes)[sig_genes])

}

.downsamp = function (umis, n, replace = F) {
        m = nrow(umis)
        .downsamp_one=function(v,n, replace = F){
                a = tabulate(sample(rep(1:length(v),times=v),replace=replace,size=n),nbins=m)
                return (a)
        }
	ret = apply(umis[, colSums(umis) >= n], 2, .downsamp_one, n)
	rownames(ret) = rownames(umis)
	return(ret)
}

mc_compute_unnorm_fp = function(mc, us, mc_cores = 16)
{
        f_g_cov = rowSums(us) > 10
	cells = intersect(names(mc), colnames(us))
        all_gs = rownames(us[f_g_cov, cells])
        n_g = length(all_gs)
        g_splts = split(all_gs, 1+floor(mc_cores*(1:n_g)/(n_g+1)))
        fnc = function(gs) {
                                        .row_stats_by_factor(us[gs, cells],
                                                                        mc[cells],
                                                                        function(y) {exp(rowMeans(log(1+y)))-1}) }

        clust_geomean = do.call(rbind, mclapply(g_splts, fnc, mc.cores = mc_cores))

        mc_meansize = tapply(colSums(us[,cells]), mc[cells], mean)
        ideal_cell_size = pmin(1000, median(mc_meansize))
        g_fp = t(ideal_cell_size*t(clust_geomean)/as.vector(mc_meansize))
        return(g_fp)
}

.row_stats_by_factor = function (data, fact, rowFunction = rowMeans) {
        u = as.character(sort(unique(fact)))
        fact[is.na(fact)] = F
        n=length(u)
        centers = matrix(NA,dim(data)[1], n, dimnames = list(rownames(data), u))
        for (i in u) {
                if(sum(fact==i, na.rm=T)>1) {
                        centers[,i] = rowFunction(data[,fact==i,drop=F])
                } else {
                        centers[,i] = data[,fact==i]
                }
        }
        return(centers)
}

mc_compute_fp_abs = function(mc, us, norm_to_med = T)
{
        f_g_cov = rowSums(us) > 10

        mc_cores = 16
        doMC::registerDoMC(mc_cores)
        all_gs = rownames(us[f_g_cov,])
        n_g = length(all_gs)
        g_splts = split(all_gs, 1+floor(mc_cores*(1:n_g)/(n_g+1)))
        fnc = function(gs) {
                                        .row_stats_by_factor(us[gs,],
                                                                        mc,
                                                                        function(y) {exp(rowMeans(log(1+y)))-1}) }
        clust_geomean = do.call(rbind, lapply(g_splts, fnc))
        mc_meansize = tapply(colSums(us), mc, mean)
        ideal_cell_size = pmin(1000, median(mc_meansize))
        g_fp = clust_geomean #t(ideal_cell_size*t(clust_geomean)/as.vector(mc_meansize))
        #normalize each gene
        fp_reg = 0.1
	if (norm_to_med) {
	        g_fp = (fp_reg+g_fp)/apply(fp_reg+g_fp, 1, median)
	}
        return(g_fp)
}

plot_two_genes_fp = function(mc_id, ga, gb, log = T) {
	sc_cl = scdb_mc(mc_id)
	fp = sc_cl@mc_fp
	if (log) { fp = log2(fp)}
	a = fp[ga,]; b = fp[gb,]
	plot(a,b, xlab = ga, ylab = gb, pch = 21, cex = 3, bg = sc_cl@colors)
	text(a,b, names(a))
	return(data.frame(a = a, b = b))
}

stirlings_lmulti = function(vec) { 
	vec = vec[vec > 0]
	n = sum(vec)
	log_stirl = 0.5 * (log(n) - (length(vec) - 1) * log(2 * pi) - sum(log(vec))) + 
		n * log(n) - sum(vec * log(vec))
	log_stirl
}

explore_gene_set = function(mc_id, mat_id, cells, ct, outdir, prefix, comb,
	downsamp_n = 500, pos_gene = NULL, neg_gene = NULL, batch_attr = "amp_batch_id",
	gradient = colorRampPalette(c("blue", "white", "red"))(1000)) {

        dir.create(outdir)
	sc_cl = scdb_mc(mc_id); sc_mat = scdb_mat(mat_id)
	umis = as.matrix(sc_mat@mat[,cells])
	nms = names(ct)
	mat_ds = as.matrix(scm_downsamp(umis, downsamp_n))

	C = cor(t(log2(1 + mat_ds[nms, ])), m = "spearman")
	diag(C) = NA
	ord = order(ct, -rowSums(umis[nms, ]))
	cls = cumsum(table(ct)) / length(ct)

	modules_n = apply(mat_ds[nms,], 2, tapply, as.vector(ct), sum)
	modules_c = t(apply(modules_n,1,tapply, sc_cl@mc[colnames(modules_n)], mean))
	modules_c = log2((1 + modules_c) / (1 + rowMedians(modules_c)))
	modules_c[ is.infinite(modules_c)] = NA

	chc = hclust(dist(cor(modules_c))); 
	if (!is.null(pos_gene)) {
		chc = as.hclust(reorder(as.dendrogram(chc),
                      modules_c[ct[pos_gene],] - modules_c[ct[neg_gene],],
                      agglo.FUN=mean))
	}
	cord = chc$order; ccls = NA
	modules_c = modules_c[,cord]
	modules_c = modules_c[ order(max.col(modules_c)),]

	message("Plotting modules on clusters")
	png(paste0(outdir, "/", prefix, "_modules_on_clusts.png"), height=2000, width=2000)
	par (mar = c(1,20,0,10), fig = c(0,1,0.8,1))	
	image.2(modules_c, col = gradient, balance = T)
	abline(v = ccls)
	par(mar = c(3,20,0,10), fig = c(0,1,0,0.8), new=T)
	create_batch_matrix(mc_id, mat_id, comb, NULL, batch_attr,
		cord = colnames(modules_c), color_batches=T)
	abline(v = ccls)
	dev.off()

	df = data.frame(ct = ct, umicount = rowSums(umis[nms, cells]))
	df = df[ order(factor(df$ct, levels = rownames(modules_c)), -df$umicount),]
	write.table(df, sep = "\t", quote = F, col.names = NA, file = paste0(outdir, "/", prefix, "_modules.txt"))

	message("Plotting Truth matrix")
	cls = cumsum(table(factor(ct, levels = rownames(modules_c)))) / length(ct)
	IM = matrix(NA, nrow = length(ct), ncol = ncol(modules_c), dimnames = list(names(ct), colnames(modules_c)))
	shared_genes = intersect(rownames(sc_cl@mc_fp), names(ct))
	IM[shared_genes,] = log2(sc_cl@mc_fp[shared_genes, colnames(modules_c)])
	IM = IM[rownames(df), colnames(modules_c)]
	height = 10 * nrow(IM)
	png(paste0(outdir, "/", prefix, "_trurh.png"), height=height, width=2000)
	image.2(IM, col = gradient, balance = T)
	abline(v=ccls, h = cls)
	mtext(names(cls), las = 2, side = 4, at = rowMeans(cbind(c(0,cls[-length(cls)]),cls)))
	dev.off()

	message("Plotting CT correlation")
	png(paste0(outdir, "/", prefix, "_cor.png"), height=height, width=height)
	par(mar=c(10,10,3,3))
	image.2(C[rownames(df),rownames(df)], col = gradient, balance = T)
	abline(h=cls,v=cls)
	mtext(names(cls), las = 2, side = 3, at = rowMeans(cbind(c(0,cls[-length(cls)]),cls)))
	mtext(names(cls), las = 2, side = 4, at = rowMeans(cbind(c(0,cls[-length(cls)]),cls)))
	dev.off()
	
	colnames(modules_c)
}


import_metacell_structure = function(id, folder=id, all_id="all", bad_genes = c(), url=NULL) {

	if (is.null(url)) {
		sin_tab = read.delim(paste0(folder, "/mc.txt"), stringsAsFactor=F)
		db_tab = read.delim(paste0(folder, "/db_mc.txt"), stringsAsFactor=F)
		col_tab = read.delim(paste0(folder, "/colors.txt"), stringsAsFactor=F)
		color_key = read.delim(paste0(folder, "/color_key.txt"), stringsAsFactor=F)
		mc2d_mc = read.delim(paste0(folder, "/mc2d_mc.txt"), stringsAsFactor=F, row.names=1)
		mc2d_sc	= read.delim(paste0(folder, "/mc2d_sc.txt"), stringsAsFactor=F, row.names=1)
		mc2d_graph = read.delim(paste0(folder, "/mc2d_graph.txt"), stringsAsFactor=F)
		cgraph = read.delim(paste0(folder, "/cgraph.txt"), stringsAsFactor=F)
		gset_tab = read.delim(paste0(folder, "/gset.txt"), stringsAsFactor=F)
	} else {
		sin_tab = read.delim(text = getURL(paste0(url, "/mc.txt")), stringsAsFactor=F)
		db_tab = read.delim(text = getURL(paste0(url, "/db_mc.txt")), stringsAsFactor=F)
		col_tab = read.delim(text = getURL(paste0(url, "/colors.txt")), stringsAsFactor=F)
		color_key = read.delim(text = getURL(paste0(url, "/color_key.txt")), stringsAsFactor=F)
		mc2d_mc = read.delim(text = getURL(paste0(url, "/mc2d_mc.txt")), stringsAsFactor=F)
		mc2d_sc = read.delim(text = getURL(paste0(url, "/mc2d_sc.txt")), stringsAsFactor=F)
		mc2d_graph = read.delim(text = getURL(paste0(url, "/mc2d_graph.txt")), stringsAsFactor=F)
		cgraph = read.delim(text = getURL(paste0(url, "/cgraph.txt")), stringsAsFactor=F)
		gset_tab = read.delim(text = getURL(paste0(url, "/gset.txt")), stringsAsFactor=F)
	}

	sin_mc = sin_tab[,2]; names(sin_mc) = sin_tab[,1]
	db_mc = db_tab[,2]; names(db_mc) = db_tab[,1]

	all_mat = scdb_mat(all_id)
	id_s = paste0(id, "_singlets")
	id_d = paste0(id, "_PIC")

	all_cells = union(names(sin_mc), names(db_mc))
	mcell_mat_ignore_genes(new_mat_id=id, mat_id="all", bad_genes, reverse=F)
	
	bad_cells = setdiff(all_mat@cells, all_cells)
	mcell_mat_ignore_cells(id, id, union(bad_cells, mat@ignore_cells))
	mcell_mat_ignore_cells(id_s, id, unique(c(bad_cells, names(db_mc), mat@ignore_cells)))
	mcell_mat_ignore_cells(id_d, id, unique(c(bad_cells, names(sin_mc), mat@ignore_cells)))

	sin_mat = scdb_mat(id_s)
	sin_cl = tgMCCov(sin_mc, setdiff(sin_mat@cells, names(sin_mc)), sin_mat)
	colors = rep(NA, ncol(sin_cl@mc_fp))
	col_vec = col_tab[,2]; names(col_vec) = col_tab[,1]
	colors[as.numeric(names(col_vec))] = col_vec
	sin_cl@colors = colors
	sin_cl@color_key = rbind(sin_cl@color_key, color_key)
	scdb_add_mc(id_s, sin_cl)

	db_mat = scdb_mat(id_d)
        db_cl = tgMCCov(db_mc, setdiff(db_mat@cells, names(db_mc)), db_mat)
        scdb_add_mc(id_d, db_cl)

	mc_x = mc2d_mc$x; names(mc_x) = rownames(mc2d_mc)
	mc_y = mc2d_mc$y; names(mc_y) = rownames(mc2d_mc)
	sc_x = mc2d_sc$x; names(sc_x) = rownames(mc2d_sc)
	sc_y = mc2d_sc$y; names(sc_y) = rownames(mc2d_sc)
	sin_2d = tgMC2D(mc_id=id_s, mc_x=mc_x, mc_y=mc_y, sc_x=sc_x, sc_y=sc_y, graph=mc2d_graph)
	scdb_add_mc2d(id_s, sin_2d)
	
	sin_cgraph = tgCellGraph(cgraph, sin_mat@cells)
	scdb_add_cgraph(id_s, sin_cgraph)

	gset = gset_tab[,2]; names(gset) = gset_tab[,1]
	sc_gset = tgGeneSets(gset)
	scdb_add_gset(id, sc_gset)
}
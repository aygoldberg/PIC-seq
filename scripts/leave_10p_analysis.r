#########################
#                       #
# Reproduce Figure S2g  #
#                       #
#########################

message("Generating Figure S2g")

##########
# leave 10% out analysis

cvdir = paste0(supdir, "/mle_cv")
dir.create(cvdir)

genes = rownames(umis); 
k = 10
gene_ct = seq_along(genes);
gene_ct = as.numeric(cut(gene_ct, quantile(gene_ct, (0:k)/k), include.lowest = T))
gene_ct = sample(gene_ct)
names(gene_ct) = genes

#write.table(gene_ct, sep = "\t", quote = F, col.names = F, file = paste0(cvdir, "/gene_ct.txt"))
#gene_table = read.table(paste0(cvdir, "/gene_ct.txt"), sep = "\t", row.names = 1, stringsAsFactor=F)

gene_ct = gene_table[,1]; names(gene_ct) = rownames(gene_table)

comb = paste0(cell_stats$timepoint, ".", cell_stats$date);
names(comb) = rownames(cell_stats)

mles = list()
for (i in seq_len(k)) {
	naughty_genes = names(which(gene_ct == i))
	message("Iteration: ", i)

	mles[[i]] = run_pic_seq(id_s, id_s, ds[,good_pics], t_cells, dc_cells,
		lr_features, mle_features, paste0(cvdir, "/temp.png"), numis = 1000, bad_genes = union(bad_genes, naughty_genes),
		comb = comb, reg = 1e-4,  downsample=F)
}

marker_cor = matrix(NA, nrow = length(mle_features), ncol = 4, dimnames = list(mle_features, c("full", "alpha", "mle", "empty")))
for (i in seq_len(k)) {
	mle = mles[[i]]
	exp_us = generate_expected_pics_from_mle(id_s, mle[,c("a_mc", "b_mc")], mle$alpha, colSums(ds[rownames(sin_cl@e_gc), rownames(mle)]), bad_genes = bad_genes)
	naughty_genes = names(which(gene_ct == i))
	dus = ds[ rownames(exp_us), colnames(exp_us)]
	markers2 = intersect(naughty_genes, intersect(mle_features, rownames(dus)))
	C = cor(t(log(1 + exp_us[markers2,])), t(log(1 + dus[markers2,])))
	marker_cor[markers2, "full"] = diag(C)

	mle_s1 = mle
	mle_s1$a_mc = sample(mle_s1$a_mc); mle_s1$b_mc = sample(mle_s1$b_mc)
	exp_us = generate_expected_pics_from_mle(id_s, mle_s1[,c("a_mc", "b_mc")], mle_s1$alpha, colSums(ds[rownames(sin_cl@e_gc), rownames(mle)]), bad_genes = bad_genes)
	C = cor(t(log(1 + exp_us[markers2,])), t(log(1 + dus[markers2,])))
	marker_cor[markers2, "alpha"] = diag(C)	

	mle_s2 = mle
	mle_s2$alpha = 0.5
	exp_us = generate_expected_pics_from_mle(id_s, mle_s2[,c("a_mc", "b_mc")], mle_s2$alpha, colSums(ds[rownames(sin_cl@e_gc), rownames(mle)]), bad_genes = bad_genes)
	C = cor(t(log(1 + exp_us[markers2,])), t(log(1 + dus[markers2,])))
	marker_cor[markers2, "mle"] = diag(C)

	mle_s3 = mle_s1
	mle_s3$alpha = 0.5
	exp_us = generate_expected_pics_from_mle(id_s, mle_s3[,c("a_mc", "b_mc")], mle_s3$alpha, colSums(ds[rownames(sin_cl@e_gc), rownames(mle)]), bad_genes = bad_genes)
	C = cor(t(log(1 + exp_us[markers2,])), t(log(1 + dus[markers2,])))
	marker_cor[markers2, "empty"] = diag(C)
}

high_markers = names(which(rowSums(umis[mle_features, good_pics]) > 1000))
IM = marker_cor[ rowSums(is.na(marker_cor)) == 0,]
IM = IM[ intersect(high_markers, rownames(IM)),]
IM = IM[ order(IM[,"full"]), c(1,3,2,4)]
IM2 = tail(IM,100)
grad = rev(colorRampPalette(c("#293377", "#008783", "#54B582", "#96C882", "#ECE986", "white"))(1000))
png(paste0(supdir, "/FigS2g.png"), height=1500, width = 1000)
par(mar = c(0,0,0,0))
image.2(IM2, col = grad, annotate = "none")
dev.off()
write.table(rev(rownames(IM2)), row.names = F, quote = F, col.names = F, file = paste0(supdir, "/leave_p10_models_cor.txt"))
png(paste0(supdir, "/FigS2g_colorbar.png"), height=100, width=1000)
par(mar = rep(0,4))
image(matrix(1:1000), col = grad, axes = F)
dev.off()

#########

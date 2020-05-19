#########################
#                       #
# Reproduce Figure 2    #
#                       #
#########################

message("Generating Figure 2")

good_pics = rownames(mle_res)
outdir = "figures/figure2"
dir.create(outdir)

supdir = "figures/supp_figures_1-4/"

###########
# Compositional changes

alpha = mle_res[good_pics, "alpha"]; names(alpha) = good_pics
t_mc = mle_res[good_pics, "a_mc"]; names(t_mc) = good_pics
dc_mc = mle_res[good_pics, "b_mc"]; names(dc_mc) = good_pics
parser_t = color2name[ sin_cl@colors[ t_mc]]; names(parser_t) = good_pics
parser_dc = color2name[ sin_cl@colors[ dc_mc]]; names(parser_dc) = good_pics

nice_cells = good_pics
bad_treat = "OVA + LPS transwell"
bad_cells = rownames(cell_stats)[ cell_stats$treatment == bad_treat]
t_cells = setdiff(t_cells, bad_cells)
dc_cells = setdiff(dc_cells, bad_cells)
nice_cells = setdiff(nice_cells, bad_cells)

tp_ord = c("3h", "20h", "48h")
comb = with(cell_stats, paste0(timepoint, "@", sorting.scheme, ".", treatment)); names(comb) = rownames(cell_stats)

t_dist = rbind(table(comb[nice_cells], parser_t[nice_cells]), 
	table(comb[t_cells], sin_names[t_cells]))

t_dist = t_dist[ rowSums(t_dist) > 20,]
t_dist = t_dist[ rev(order(factor(vecsplit(rownames(t_dist), "@", 1), levels = tp_ord),
	factor(vecsplit(rownames(t_dist), "@", 2), 
		levels = c("Trbc+.OVA + LPS T only", "Trbc+.OVA + LPS", "doublets.OVA + LPS")))),]

dc_dist = rbind(table(comb[nice_cells], parser_dc[nice_cells]), 
        table(comb[dc_cells], sin_names[dc_cells]))

dc_dist = dc_dist[ rowSums(dc_dist) > 20,]
dc_dist = dc_dist[ rev(order(factor(vecsplit(rownames(dc_dist),	"@", 1), levels = tp_ord),
	factor(vecsplit(rownames(dc_dist),"@", 2), 
		levels = c("Cd11c+.OVA + LPS DC only", "Cd11c+.OVA + LPS", "doublets.OVA + LPS")))),]

t_dist = t_dist[, intersect(lin_ord, colnames(t_dist))]
dc_dist = dc_dist[, intersect(lin_ord, colnames(dc_dist))]

t_n = t_dist / rowSums(t_dist); dc_n = dc_dist / rowSums(dc_dist)

png(paste0(outdir, "/Fig2a.png"), height=1000, width=1000)
par(lwd=6)
barplot(t(t_n), horiz = T, col = name2color[ colnames(t_n)], space = c(1.5,0.3,0.3), names.arg = rep("", nrow(t_n)))
dev.off()

png(paste0(outdir, "/Fig2b.png"), height=1000, width=1000)
par(lwd=6)
barplot(t(dc_n), horiz = T, col = name2color[ colnames(dc_n)], space = c(1.5,0.3,0.3), names.arg = rep("", nrow(dc_n)))
dev.off()

##############
# Compute the expected transcription of PICs

numis = 1000
ds = .downsamp(umis[,good_pics], numis)
ds_f = ds[setdiff(rownames(sin_cl@e_gc), bad_genes),]
us = umis[ rownames(ds_f), colnames(ds_f)]
cells = intersect(colnames(ds_f), rownames(mle_res))
exp_us = generate_expected_pics_from_mle(id_s, mle_res[good_pics, c("a_mc", "b_mc")], mle_res[ good_pics, "alpha"], 
	colSums(us), bad_genes = bad_genes)
exp_n = generate_expected_pics_from_mle(id_s, mle_res[good_pics, c("a_mc", "b_mc")], mle_res[ good_pics, "alpha"], 
	colSums(ds_f), bad_genes = bad_genes)
t_n = generate_expected_pics_from_mle(id_s, mle_res[good_pics, c("a_mc", "b_mc")], rep(1, length(good_pics)),
	colSums(ds_f) * alpha, bad_genes = bad_genes)
dc_n = generate_expected_pics_from_mle(id_s, mle_res[good_pics, c("a_mc", "b_mc")], rep(0, length(good_pics)),
        colSums(ds_f) * (1 - alpha), bad_genes = bad_genes)
genes = rownames(exp_us)
y = rowSums(umis[genes,good_pics]); x = rowSums(exp_us[genes,good_pics]);

# generate PICs expected UMIs when modeled as 2T+DC triplets
tr_res = triplet_res$mle_res
triplet_mc = tr_res[cells, c("forward_a1_mc", "forward_a2_mc", "b_mc")]
triplet_alpha = with(tr_res[cells,], cbind(forward_alpha, 1 - forward_alpha) *  alpha); 
exp_tr_us = generate_expected_pics_from_mle(id_s, triplet_mc, triplet_alpha,
        colSums(us), bad_genes = bad_genes)
x2 = rowSums(exp_tr_us[genes, cells])

################

reg = 40
sum_t = rowSums(t_n[genes,]); sum_dc = rowSums(dc_n[genes,])
z = log2((10 + sum_dc) / (10 + sum_t))
grad = colorRampPalette(c("limegreen", "gray40", "firebrick3"))(101)
disp_genes = names(which(abs(log2((y+reg)/(x+reg))) > 1))
val = z[disp_genes]; zlim = max(abs(val))
val_n = round((val + zlim) / (2 * zlim) * 100) + 1
lx = log2(reg+x); ly = log2(reg+y)
lim = quantile(c(lx,ly), c(0,1))
png(paste0(outdir, "/Fig2c.png"), height=1000, width=1000)
plot(lx, ly, pch = 20, col = "gray", cex = 2,
        xlim = lim, ylim = lim, axes = F, xlab = "", ylab = "")
axis(1); axis(2);
abline(coef = c(1,1), lty = 2); abline(coef = c(-1,1), lty = 2)
points(lx[ disp_genes], ly[disp_genes], cex = 3, pch = 21, bg = grad[val_n[disp_genes]])
dev.off()

png(paste0(outdir, "/Fig2c_text.png"), height=1000, width=1000)
xy_scatter_genes(x,y, col = c(NA,NA), text = T, reg = reg)
dev.off()

z1 = log2((y+reg)/(x+reg)); z2 = log2((y+reg)/(x2+reg))
png(paste0(supdir, "/FigS4b.png"), height=1000, width=1000)
plot(z1, z2, type = "n"); grid(col = "black"); points(z1, z2, cex = 2, pch = 20, col = "navyblue")
dev.off()
png(paste0(supdir, "/FigS4b_text.png"), height=1000, width=1000)
plot(z1, z2, type = "n"); grid(col = "black"); text(z1, z2, names(z1), col = "navyblue")
dev.off()

#########
# Compute significant changes between PICs observed and expected gene expression

comb = with(cell_stats, 
	paste0(sorting.scheme, ".", treatment, ".", timepoint)); names(comb) = rownames(cell_stats)
pic_joint = interaction(factor(parser_dc, levels = lin_ord),
        factor(parser_t, levels = lin_ord)); names(pic_joint) = names(parser_t)

pic_comb = interaction(pic_joint[good_pics], comb[good_pics])
names(pic_comb) = good_pics

good_combs = names(which(table(pic_comb) > 30))
analyzed_pics = good_pics[ pic_comb %in% names(which(table(pic_comb) > 30))]
pic_comb = factor(pic_comb[analyzed_pics], levels = good_combs)

real_m = t(apply(ds_f[,analyzed_pics], 1, tapply, pic_comb[analyzed_pics], sum))
exp_m =  t(apply(exp_n[,analyzed_pics], 1, tapply, pic_comb[analyzed_pics], sum))
x2 = rowSums(((real_m - exp_m)^2 )/exp_m)
qvals = p.adjust(1 - pchisq(x2, df = ncol(real_m) - 1), "fdr")
z = log2(real_m / exp_m);
z[real_m==0] = NA

y = log2(rowSums(us[names(qvals),analyzed_pics]));

genes = union(c(), setdiff(names(which(qvals < 1e-6 & apply(us[ names(qvals), analyzed_pics], 1, function(x) sort(x,T)[3]) > 3 & y > 6)), bad_genes))

z_reg = log2((real_m + 5) / (exp_m + 5))
IM = z_reg[genes, ]#grep(paste0("^|\\.", pop, "\\.|$"), colnames(z))]
IM = IM[rowSums(!is.na(IM)) > 0,]

library(tglkmeans)

tp_comp = factor(vecsplit(colnames(IM), "\\.", 5), levels = c("3h", "20h", "48h"))
t_comp = vecsplit(colnames(IM), "\\.", 2)
dc_comp = vecsplit(colnames(IM), "\\.", 1)
samp_ord = colnames(IM)[order(tp_comp, factor(t_comp, levels = lin_ord), factor(dc_comp, levels = lin_ord))]

tp_cols = c("lightskyblue1", "lightslateblue", "blue")
IM2 = IM[,samp_ord]
k = 30
data = as.data.frame(IM2)
data$id = rownames(data)
data = data[,c(ncol(data), 1:(ncol(data) - 1))]
km <- TGL_kmeans_tidy(data, k=k, metric='euclid', verbose=TRUE, seed = 18)
centers = as.matrix(km$centers[,-1]); rownames(centers) = seq_len(k)
x_min = apply(centers,1,min); x_max = apply(centers,1,max)
centers = centers[ order(abs(x_min) < x_max, max.col(centers)),]
km_clusts = as.numeric(factor(km$cluster$clust, levels = rownames(centers))); names(km_clusts) = rownames(data)

png(paste0(supdir, "/FigS4d.png"), height=max(2000, nrow(IM2) * 12), width=1000)
par(mar = c(15,15,3,3))
image.2(IM2, balance = T, annotate = "both", hct = km_clusts)
dev.off()

############



############

marker_genes = c("Foxp3", "Ikzf2", "Il22", "Tnfrsf9", "Hopx", "Cxcr6")

real_m = t(apply(ds_f[marker_genes,analyzed_pics], 1, tapply, pic_comb[analyzed_pics], mean))
exp_m =  t(apply(exp_n[marker_genes,analyzed_pics], 1, tapply, pic_comb[analyzed_pics], mean))
t_m =  t(apply(t_n[marker_genes,analyzed_pics], 1, tapply, pic_comb[analyzed_pics], mean))
dc_m =  t(apply(dc_n[marker_genes,analyzed_pics], 1, tapply, pic_comb[analyzed_pics], mean))
reg = 0.1
z_m = log2((dc_m + reg) / (t_m + reg));
zlim = max(abs(z_m))

library(Hmisc)
m = t(apply(ds_f[marker_genes, analyzed_pics],1,tapply, pic_comb[analyzed_pics], sum))
n = tapply(colSums(ds_f[,analyzed_pics]),pic_comb[analyzed_pics],sum)
Y = sweep(t(apply(m,1,binconf,n)), 2, rep(tapply(colSums(ds_f[,analyzed_pics]), pic_comb[analyzed_pics], mean), 3), "*")

circ_frac = 0.1
cls = cumsum(table(tp_comp))
cls_s = cls[-length(cls)]

tp_comp = factor(vecsplit(colnames(real_m), "\\.", 5), levels = c("3h", "20h", "48h"))
t_comp = factor(vecsplit(colnames(real_m), "\\.", 2), levels = lin_ord)
dc_comp = factor(vecsplit(colnames(real_m), "\\.", 1), levels = lin_ord)
ord = order(tp_comp, t_comp, dc_comp);

gene = marker_genes[1]
exp_tab = rbind(t_m[gene,ord], dc_m[gene,ord], 0, 0)
real_tab = rbind(0, 0, real_m[gene,ord], 0)
tab = cbind(exp_tab, real_tab)[, rep(seq_along(ord), each = 2) + rep(c(0, length(ord)), length(ord))]
X = barplot(tab)

png(paste0(outdir, "/Fig2d.png"), height=3000, width=2400)
x0_fig = rep(c(0,0.5), each = 3)
x1_fig = rep(c(0.5,1), each = 3)
y0_fig = rep(c(0.1,0.4,0.7), 2)
y1_fig = rep(c(0.4,0.7,1), 2)
for (i in 1:6) {
        gene = marker_genes[i]
        par(fig = c(x0_fig[i], x1_fig[i], y0_fig[i], y1_fig[i]), new=(i>1), lwd=5)
        exp_tab = rbind(t_m[gene,ord], dc_m[gene,ord], 0, 0)
        real_tab = rbind(0, 0, real_m[gene,ord], 0)
        tab = cbind(exp_tab, real_tab)[, rep(seq_along(ord), each = 2) + rep(c(0, length(ord)), length(ord))]
        mtab = round(max(Y[gene,], na.rm=T),2)
        y_dc = -circ_frac * mtab; y_t = y_dc * 2
        ylim = c(0, mtab)
        plot(X, rep(0, length(X)), ylim = ylim, axes=F, xlab = "", ylab = "", type = "n")
        X = barplot(tab, col = c("limegreen", "firebrick3", "gray70", "white"), ylim = ylim, las = 2,
                names.arg = rep("", ncol(tab)), space = c(0.5,0), axes = F, add=T, main = gene, cex.main=2)
        obs_coords = X[seq(2,length(X),2)]; ci.l = Y[gene, ord + length(n)]; ci.u = Y[gene, ord + 2 * length(n)]
        segments(obs_coords, ci.l, y1 = ci.u);
        segments(obs_coords-0.2, ci.l, x1 = obs_coords + 0.2); segments(obs_coords-0.2, ci.u, x1 = obs_coords + 0.2);
        axis(2, at = c(0, mtab), las = 2)
        #X = rowMeans(cbind(X[seq(1,length(X),3)], X[seq(2,length(X),3)]))
}

X2 = rowMeans(cbind(X[seq(1,length(X),2)], X[seq(2,length(X),2)]))
par(fig = c(0,0.5,0,0.1), new=T)
plot(X2, rep(0, length(X2)), type="n", axes=F, xlab="", ylab="", ylim = c(-0.3,1.3), xlim = quantile(X,c(0,1)))
segments(X2, 0, y1 = 1)
points(X2, rep(0, length(X2)), pch = 21, bg = name2color[as.vector(dc_comp[ord])], cex = 6)
points(X2, rep(1, length(X2)), pch = 21, bg = name2color[as.vector(t_comp[ord])], cex = 6)
box()
par(fig = c(0.5,1,0,0.1), new=T)
plot(X2, rep(0, length(X2)), type="n", axes=F, xlab="", ylab="", ylim = c(-0.3,1.3), xlim = quantile(X, c(0,1)))
segments(X2, 0, y1 = 1)
points(X2, rep(0, length(X2)), pch = 21, bg = name2color[as.vector(dc_comp[ord])], cex = 6)
points(X2, rep(1, length(X2)), pch = 21, bg = name2color[as.vector(t_comp[ord])], cex = 6)
box()
dev.off()

#############################

marker_genes = c("Pdcd1", "Tnfsf4", "Il2", "Icos", "Cd40lg", "Top2a", "Gzmb", "Tigit", "Ifngr1", "Nr4a2")

real_m = t(apply(ds_f[marker_genes,analyzed_pics], 1, tapply, pic_comb[analyzed_pics], mean))
exp_m =  t(apply(exp_n[marker_genes,analyzed_pics], 1, tapply, pic_comb[analyzed_pics], mean))
t_m =  t(apply(t_n[marker_genes,analyzed_pics], 1, tapply, pic_comb[analyzed_pics], mean))
dc_m =  t(apply(dc_n[marker_genes,analyzed_pics], 1, tapply, pic_comb[analyzed_pics], mean))
reg = 0.1
z_m = log2((dc_m + reg) / (t_m + reg));
zlim = max(abs(z_m))

library(Hmisc)
m = t(apply(ds_f[marker_genes, analyzed_pics],1,tapply, pic_comb[analyzed_pics], sum))
n = tapply(colSums(ds_f[,analyzed_pics]),pic_comb[analyzed_pics],sum)
Y = sweep(t(apply(m,1,binconf,n)), 2, rep(tapply(colSums(ds_f[,analyzed_pics]), pic_comb[analyzed_pics], mean), 3), "*")

dir.create(paste0(supdir, "/FigS4c/"))
for (gene in marker_genes) {
	png(paste0(supdir, "/FigS4c/", gene, ".png"), height=1000, width=1800)
	par(lwd=5)
	exp_tab = rbind(t_m[gene,ord], dc_m[gene,ord], 0, 0)
	real_tab = rbind(0, 0, real_m[gene,ord], 0)
	tab = cbind(exp_tab, real_tab)[, rep(seq_along(ord), each = 2) + rep(c(0, length(ord)), length(ord))]
	mtab = round(max(Y[gene,], na.rm=T),3)
	y_dc = -circ_frac * mtab; y_t = y_dc * 2
	ylim = c(y_t + y_dc, mtab)
	X = barplot(tab, col = c("limegreen", "firebrick3", "gray70", "white"), ylim = ylim, xaxs = "i", las = 2,
		names.arg = rep("", ncol(tab)), space = c(0.5,0), axes = F)
	obs_coords = X[seq(2,length(X),2)]; ci.l = Y[gene, ord + length(n)]; ci.u = Y[gene, ord + 2 * length(n)]
	segments(obs_coords, ci.l, y1 = ci.u); 
	segments(obs_coords-0.2, ci.l, x1 = obs_coords + 0.2); segments(obs_coords-0.2, ci.u, x1 = obs_coords + 0.2);
	axis(2, at = c(0, mtab), las = 2)
	X2 = rowMeans(cbind(X[seq(1,length(X),2)], X[seq(2,length(X),2)]))
	segments(X2, y_dc, y1 = y_t)
	points(X2, rep(y_dc, length(X2)), pch = 21, bg = name2color[ as.vector(dc_comp[ord])], cex = 10)
	points(X2, rep(y_t, length(X2)), pch = 21, bg = name2color[ as.vector(t_comp[ord])], cex = 10)
	dev.off()
}

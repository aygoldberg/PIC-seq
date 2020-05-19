#################################
#                       	#
# Reproduce Supp Figures 1-4    #
#                       	#
#################################

message("Generating supplementary figures 1-4")

b = 15
dir.create(paste0(supdir, "/FigS1b"))
palette_t = colorRampPalette(c("white", "orange", "tomato", "mediumorchid4", "midnightblue"))
palette = palette_t(b)
nms = c("Trbc2", "Fscn1")
for (val in nms) {
        vals = sin_n[val,cells]
        norm_val = rep(1, length(vals)); names(norm_val) = names(vals)
        norm_val[ vals != 0] = as.numeric(cut(vals[ vals != 0], unique(quantile(vals[ vals != 0], (0:b)/b)), include.lowest = T)) + 1;
        cols = colorRampPalette(palette)(max(norm_val))
        png(paste0(supdir, "/FigS1b/", val, ".png"), height = 1000, width = 1000)
        plot(sin_2d@sc_x[cells], sin_2d@sc_y[cells], pch = 20, col = "gray80", cex=1, axes = F, xlab = "", ylab = "")
        exp_cells = intersect(cells, names(which(vals > 0)))
        points(sin_2d@sc_x[exp_cells], sin_2d@sc_y[exp_cells], cex = 1 + 0.4 * round((norm_val[exp_cells] - 1) / max(norm_val) * 5),
                pch = 21, bg = cols[norm_val[exp_cells]])
        axis(1); axis(2)
        dev.off()
}

png(paste0(outdir, "/FigS1b_colorbar.png"), height=100, width = 1000)
par(mar = rep(0,4))
image(matrix(1:b), axes = F, col = palette)
dev.off()


############
# batch distribution

supdir = "figures/supp_figures_1-4/"
sin_vec = sin_names[cells]
sin_ord = cells[ order(factor(sin_vec, levels = lin_ord))]

cord = unique(sin_cl@mc[ sin_ord])
batch_comb = unique(with(sin_stats[sin_ord,], paste0(amp_batch_id, "!", sorting.scheme, "@", timepoint, "@", treatment, "@", date)))
B = vecsplit(batch_comb, "!", 1); batch_comb = vecsplit(batch_comb, "!", 2); names(batch_comb) = B
comb_sort    = factor(vecsplit(batch_comb, "@", 1), levels = c("Trbc+", "Cd11c+"))
comb_tp      = factor(vecsplit(batch_comb, "@", 2), levels = c("3h", "20h", "48h"))
comb_treat   = factor(vecsplit(batch_comb, "@", 3), levels = c("OVA + LPS T only", "OVA + LPS DC only", "OVA + LPS"))
bord = B[ order(comb_sort, comb_tp, comb_treat)]

cell_comb = with(sin_stats[sin_ord,], paste0(amp_batch_id, "!", sorting.scheme, "@", timepoint, "@", treatment, "@", date))
cell_comb = factor(cell_comb, levels = sapply(bord, grep, names(table(cell_comb)), v=T))
names(cell_comb) = sin_ord
create_batch_matrix(id_s, id_s, cell_comb, paste0(supdir, "/FigS1c.png"), cord = cord)

png(paste0(supdir, "/FigS1c_colorbar.png"), height=100, width=1000)
par(mar = rep(0,4))
image(matrix(1:1000), col = colorRampPalette(c("white", "navyblue"))(1000), axes = F)
dev.off()

png(paste0(supdir, "/FigS1c_metacell_bar.png"), height=100, width=1000)
par(mar = rep(0,4))
image(matrix(seq_along(cord)), col = sin_cl@colors[ cord])
dev.off()

###########
# joint expression

dist.grad = colorRampPalette(rev(c("#293377", "#008783", "#54B582", "#96C882", "#ECE986", "white")))
y = rowMeans(sin_n[,t_cells]) * min(length(t_cells), length(dc_cells))
x = rowMeans(sin_n[,dc_cells]) * min(length(t_cells), length(dc_cells))

reg = 100; z_big = log2((y + reg) / (x + reg))
num_genes = 10
t_genes = names(head(sort(z_big,T),num_genes)); dc_genes = names(head(sort(z_big,F),num_genes))
all_ds = cbind(ds, sin_ds)
x = colSums(all_ds[t_genes,]); y = colSums(all_ds[dc_genes,])

xlim = c(0,max(x)); ylim = c(0,max(y))
cells = intersect(colnames(all_ds), union(t_cells, dc_cells)); 
png(paste0(supdir, "/FigS2b.png"), height=1000, width=2000)
par(mfrow = c(1,2))
smoothScatter(x[cells], y[cells], colramp=dist.grad, xlim = xlim, ylim = ylim, axes = F, xlab = "", ylab = "")
grid(col = "black"); axis(1); axis(2)
smoothScatter(x[good_pics], y[good_pics], colramp=dist.grad, xlim = xlim, ylim = ylim, axes = F, xlab = "", ylab = "")
grid(col = "black"); axis(1); axis(2)
dev.off()

############

umis_n = sweep(umis, 2, colSums(umis), "/") * 1000
t_pop = "act"
sub_cells = names(which(parser_t == t_pop))
sub_cells = sub_cells[ parser_dc[sub_cells] %in% names(which(table( parser_dc[sub_cells]) > 50))]
pops = names(table(parser_dc[sub_cells]))
sin_cells = dc_cells[ sin_names[ dc_cells] %in% pops]
dir.create(paste0(supdir, "/FigS2h"))
reg = 20
for (dc_pop in pops) {
	sin_g1 = sin_cells[ sin_names[ sin_cells] == dc_pop]
	sin_g2 = setdiff(sin_cells, sin_g1)
	db_g1 = sub_cells[ parser_dc[ sub_cells] == dc_pop]
	db_g2 = setdiff(sub_cells, db_g1)
	sin_x = rowMeans(umis_n[,sin_g2]) * min(length(sin_g1), length(sin_g2))
	sin_y = rowMeans(umis_n[,sin_g1]) * min(length(sin_g1), length(sin_g2))
	db_x  = rowMeans(umis_n[,db_g2])  * min(length(db_g1), length(db_g2))
	db_y  = rowMeans(umis_n[,db_g1])  * min(length(db_g1), length(db_g2))
	sin_z = log2((reg + sin_y) / (reg + sin_x))	
	db_z = log2((reg + db_y) / (reg + db_x))
	xlim = c(-4,4); ylim = c(-2.5,2.5)
	png(paste0(supdir, "/FigS2h/", dc_pop, ".png"), height=1000, width=1000)
	plot(sin_z, db_z, pch = 20, col = rgb(0,0,0.5,0.5), xlim = xlim, ylim = ylim, cex = 2.5, axes = F, xlab = "", ylab = "")
	axis(1); axis(2)
	grid(col = "black"); abline(h=0,v=0,lwd=2)
	dev.off()
	png(paste0(supdir, "/", t_pop, "_pic_vs_sin/", dc_pop, "_text.png"), height=1000, width=1000)
	plot(sin_z, db_z, type = "n", xlim = xlim, ylim = ylim, axes = F)
	axis(1); axis(2)
	grid(col = "black"); abline(h=0,v=0,lwd=2)
	text(sin_z, db_z, names(sin_z), col = rgb(0,0,0.5,0.5))
	dev.off()
}

############
# Figure S3

t_cells = setdiff(names(sin_cl@mc)[ color2name[ sin_cl@colors[ sin_cl@mc]] %in% lin_ord[1:5]], c())
dc_cells = names(sin_cl@mc)[ color2name[ sin_cl@colors[ sin_cl@mc]] %in% lin_ord[6:9]]

t_modules = read.delim("annotations/t_modules.txt", stringsAsFactor=F, row.names=1)
annotations = c("Th precursor", "Type I IFN", "Early activation", "Activation", "Proliferation")
mat_ds = as.matrix(scm_downsamp(umis[,t_cells], 500))

t_nms = rownames(t_modules)
t_ct = factor(t_modules$annotation, levels = annotations); names(t_ct) = t_nms
t_nms = names(sort(t_ct))
C = cor(t(log2(1 + mat_ds[t_nms, ])), m = "spearman"); diag(C) = NA
comb = with(sin_stats, paste0(sorting.scheme, ".", treatment, ".", timepoint)); names(comb) = rownames(sin_stats)
t_m = sc_to_bulk(paste0(id_s, "_f"), id_s, comb, cells = t_cells, choose_genes = F, min_comb = 20)
IM = log(10 + t_m[t_nms,])
IM2 = t(apply(IM,1,scale))
dimnames(IM2) = dimnames(IM)
mono = c(2,1,3); coc = c(6,5,7); tw = 4;
bars_width = 0.333; bar_width = bars_width / ncol(IM2)
t_cls = cumsum(table(t_ct)) / length(t_ct)
cor_grad = colorRampPalette(c("lightseagreen", "white", "lightsalmon4"))(1000)
zlim = max(abs(IM2)) * c(-1,1)
png(paste0(supdir, "/FigS3a.png"), height = 1000, width = 1500)
par(mar = c(0.5,0.5,0.5,0.5), fig = c(0,1 - bars_width,0,1), lwd = 2)
image.2(C, balance = T, annotate = "none", col = cor_grad); box()
par(fig = c(1 - bars_width, 1 - bars_width + bar_width * 3,0,1), new = T)
image.2(IM2[,mono], zlim = zlim, annotate = "none"); box()
par(fig = c(1 - bars_width + bar_width * 3, 1 - bars_width + bar_width * 4,0,1), new = T)
image.2(IM2[,tw], zlim = zlim, annotate = "none"); box()
par(fig = c(1 - bars_width + bar_width * 4, 1 - bars_width + bar_width * 7,0,1), new = T)
image.2(IM2[,coc], zlim = zlim, annotate = "none"); box()
dev.off()

############
cells = union(t_cells, dc_cells)
sin_foc = log(1 + 7 * sin_n)
sin_mn = apply(sin_foc[names(t_ct),cells], 2, tapply, t_ct, sum)

t_stats = sin_stats[ t_cells,]
t_stats = t_stats[ t_stats$sorting.scheme != "Cd11c+" & t_stats$treatment %in% c("OVA + LPS", "OVA + LPS T only"),]
t_ass = factor(ifelse(t_stats$treatment == "OVA + LPS T only", "Mono", as.vector(t_stats$timepoint)), c("Mono", "3h", "20h", "48h"))
disp_cells = rownames(t_stats); names(t_ass) = disp_cells
IM = sin_mn[, disp_cells]
cell_ord = colnames(IM)[ order(t_ass, 
	4 * IM["Proliferation",] + IM["Activation",] - 2 * IM["Th precursor",])]
tp_cols = c("gray60", "#96BBCC", "lightslateblue", "blue")
x = rep(NA, length(disp_cells)); names(x) = disp_cells
for (i in seq_along(levels(t_ass))) {
        ass = levels(t_ass)[i]
	x[ intersect(cell_ord, names(t_ass)[t_ass == ass])] = seq(i, i+1, length.out = sum(t_ass == ass))
}
cls = seq_along(levels(t_ass))
dir.create(paste0(supdir, "/FigS3b"))
for (val in rownames(IM)) {
        png(paste0(supdir, "/FigS3b/", val, ".png"), height=500, width=1500)
        plot(x, IM[val,], pch = 20,
                col = sin_cl@colors[sin_cl@mc[ names(x)]],
                axes = F, xlab = "", ylab = "", cex = 4.5, xaxs = "i");
        axis(2)
        abline(v = cls); grid(col = "black")
        dev.off()
}

############
dc_modules = read.delim("annotations/dc_modules.txt", stringsAsFactor=F, row.names=1)
annotations = c("Type I IFN", "LPS response", "Irf8 program", "Costimulatory program")

mat_ds = as.matrix(scm_downsamp(umis[,dc_cells], 500))
dc_nms = rownames(dc_modules)
dc_ct = factor(dc_modules$annotation, levels = annotations); names(dc_ct) = dc_nms
dc_nms = names(sort(dc_ct))
C = cor(t(log2(1 + mat_ds[dc_nms, ])), m = "spearman"); diag(C) = NA
comb = with(sin_stats, paste0(sorting.scheme, ".", treatment, ".", timepoint)); names(comb) = rownames(sin_stats)
dc_m = sc_to_bulk(paste0(id_s, "_f"), id_s, comb, cells = dc_cells, choose_genes = F, min_comb = 20)
IM = log(10 + dc_m[dc_nms,])
IM2 = t(apply(IM,1,scale))
dimnames(IM2) = dimnames(IM)
mono = c(2,1,3); coc = c(6,5,7); tw = 4;
bars_width = 0.333; bar_width = bars_width / ncol(IM2)
dc_cls = cumsum(table(dc_ct)) / length(dc_ct)
cor_grad = colorRampPalette(c("lightseagreen", "white", "lightsalmon4"))(1000)
zlim = max(abs(IM2)) * c(-1,1)
png(paste0(supdir, "/FigS3c.png"), height = 1000, width = 1500)
par(mar = c(0.5,0.5,0.5,0.5), fig = c(0,1 - bars_width,0,1), lwd = 2)
image.2(C, balance = T, annotate = "none", col = cor_grad); box()
par(fig = c(1 - bars_width, 1 - bars_width + bar_width * 3,0,1), new = T)
image.2(IM2[,mono], zlim = zlim, annotate = "none"); box()
par(fig = c(1 - bars_width + bar_width * 3, 1 - bars_width + bar_width * 4,0,1), new = T)
image.2(IM2[,tw], zlim = zlim, annotate = "none"); box()
par(fig = c(1 - bars_width + bar_width * 4, 1 - bars_width + bar_width * 7,0,1), new = T)
image.2(IM2[,coc], zlim = zlim, annotate = "none"); box()
dev.off()

##########

dc_stats = sin_stats[dc_cells,]
dc_stats = dc_stats[ dc_stats$sorting.scheme == "Cd11c+" & dc_stats$treatment %in% c("OVA + LPS", "OVA + LPS DC only"),]
dc_ass = factor(ifelse(dc_stats$treatment == "OVA + LPS DC only", "Mono", as.vector(dc_stats$timepoint)), c("Mono", "3h", "20h", "48h"))
disp_cells = rownames(dc_stats); names(dc_ass) = disp_cells

sin_mn = apply(sin_foc[names(dc_ct),cells], 2, tapply, dc_ct, sum)
IM = sin_mn[, disp_cells]
tp_cols = c("gray60", "#96BBCC", "lightslateblue", "blue")

dir.create(paste0(supdir, "/FigS3d"))

for (val in rownames(IM)) {
	png(paste0(supdir, "/FigS3d/", val, ".png"), height=700, width=1000)
	boxplot(IM[val, disp_cells] ~ dc_ass[ disp_cells], outline = F, col = tp_cols, boxwex = 0.3, axes = F)
	grid(col = "black", lwd=1.5)
	boxplot(IM[val, disp_cells] ~ dc_ass[ disp_cells], outline = F, col = tp_cols, boxwex = 0.3, axes = F, add = T)
	axis(2)
	dev.off()
}

###########


tps = factor(cell_stats$timepoint, levels = c("3h", "20h", "48h"))
names(tps) = rownames(cell_stats)
tp_cols = c("lightskyblue1", "lightslateblue", "blue")
dir.create(paste0(supdir, "/FigS3e"))
for (treatment in names(table(factor(cell_stats[t_cells, "treatment"])))) {
        sub_cells = sample(t_cells[ cell_stats[t_cells, "treatment"] == treatment & cell_stats[t_cells, "timepoint"] == "20h"])
	png(paste0(supdir, "/FigS3e/T_", treatment, ".png"), height = 1000, width = 1000)
	plot(sin_2d@sc_x[cells], sin_2d@sc_y[cells], pch = 20, col = "gray80", cex=1, axes = F, xlab = "", ylab = "")
        points(sin_2d@sc_x[sub_cells], sin_2d@sc_y[sub_cells], cex = 2,
                pch = 21, bg = sin_cl@colors[sin_cl@mc[sub_cells]])
	dev.off()
}

for (treatment in names(table(factor(cell_stats[dc_cells, "treatment"])))) {
        sub_cells = sample(dc_cells[ cell_stats[dc_cells, "treatment"] == treatment & cell_stats[dc_cells, "timepoint"] == "20h"])
        png(paste0(supdir, "/FigS3e/DC_", treatment, ".png"), height = 1000, width = 1000)
        plot(sin_2d@sc_x[cells], sin_2d@sc_y[cells], pch = 20, col = "gray80", cex=1, axes = F, xlab = "", ylab = "")
        points(sin_2d@sc_x[sub_cells], sin_2d@sc_y[sub_cells], cex = 2,
		pch = 21, bg = sin_cl@colors[sin_cl@mc[sub_cells]])
        dev.off()
}


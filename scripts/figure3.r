#########################
#                       #
# Reproduce Figure 3    #
#                       #
#########################

message("Generating Figure 3")

###########
# Compositional changes

good_pics = rownames(mle_res)
outdir = "figures/figure3"
dir.create(outdir)

alpha = mle_res[good_pics, "alpha"]; names(alpha) = good_pics
t_mc = mle_res[good_pics, "a_mc"]; names(t_mc) = good_pics
dc_mc = mle_res[good_pics, "b_mc"]; names(dc_mc) = good_pics
parser_t = color2name[ sin_cl@colors[ t_mc]]; names(parser_t) = good_pics
parser_dc = color2name[ sin_cl@colors[ dc_mc]]; names(parser_dc) = good_pics

###########


cells = union(t_cells, dc_cells)

png(paste0(outdir, "/Fig3g.png"), height = 1000, width = 1000)
plot(sin_2d@sc_x[cells], sin_2d@sc_y[cells], pch = 21, bg = sin_cl@colors[ sin_cl@mc[ cells]],
        axes = F, xlab = "", ylab = "", cex = 1.5)
dev.off()

############

cells = union(t_cells, dc_cells)
sin_ds = .downsamp(umis[,cells], 500)
good_pics = rownames(mle_res) #[ !(mle_res$forward_triplet | mle_res$reverse_triplet)]

alpha = mle_res[good_pics, "alpha"]; names(alpha) = good_pics
t_mc = mle_res[good_pics, "a_mc"]; names(t_mc) = good_pics
dc_mc = mle_res[good_pics, "b_mc"]; names(dc_mc) = good_pics
parser_t = color2name[ sin_cl@colors[ t_mc]]; names(parser_t) = good_pics
parser_dc = color2name[ sin_cl@colors[ dc_mc]]; names(parser_dc) = good_pics

t_col = "limegreen"; dc_col =  "firebrick3"

t_clusts = names(table(sin_cl@mc[ t_cells]))
t_ord = t_clusts[order(factor(color2name[ sin_cl@colors[ as.numeric(t_clusts)]], levels = lin_ord))]

dc_clusts = names(table(sin_cl@mc[ dc_cells]))
dc_ord = dc_clusts[order(factor(color2name[ sin_cl@colors[ as.numeric(dc_clusts)]], levels = lin_ord))]

clust_ord = c(t_ord, dc_ord); 
cells = setdiff(union(t_cells, dc_cells), rownames(cell_stats)[ cell_stats$treatment == "OVA + LPS transwell"])

umicount = colSums(umis)
ylim = quantile(log2(umicount[union(cells, good_pics)]), c(0,1))
t_nms = rev(read.table(paste0(outdir, "/t_nms.txt"), stringsAsFactor=F)[[1]])
dc_nms = rev(read.table(paste0(outdir, "/dc_nms.txt"), stringsAsFactor=F)[[1]])

png(paste0(outdir, "/Fig3h.png"), height=2000, width=2000)
sin_vec = sin_names[cells]
par(mar = c(0.5,3,0.5,0), fig = c(0,1,0,0.4), lwd = 2)
sin_ord = plot_sc_heatmap(paste0(id_s, "_f"), id, t_nms, clusts = sin_vec, good_clusts = lin_ord, cells = cells, annotate=F, normalize=T, lty=1, lwd=2); box()
par(fig = c(0,1,0.4,0.8), new = T)
sin_ord = plot_sc_heatmap(paste0(id_s, "_f"), id, dc_nms, clusts = sin_vec, good_clusts = lin_ord, cells = cells, annotate=F, normalize=T, lty=1, lwd=2); box()
par(fig=c(0,1,0.85,0.9), new=T)
image(matrix(1 * (cell_stats[ sin_ord, "treatment"] == "helminths")), axes = F, col = c("gray80", "gray20")); box()
par(fig=c(0,1,0.8,0.85), new=T)
image(matrix(1 * (cell_stats[ sin_ord, "sorting.scheme"] == "Ag+ Cd11c+")), axes = F, col = c("white", "green2")); box()
par(fig=c(0,1,0.9,1), new=T)
plot(log2(umicount[sin_ord]), col = rgb(0.2,0.2,0.2,0.4), pch = 20, xaxs = "i", axes = F, xlab = "", ylab = "", cex = 1.5, ylim = ylim); 
axis(2, las = 2); box()
dev.off()

png(paste0(outdir, "/Fig3i.png"), height=200, width=2000)
par(mar = c(0.5,3,0.5,0), fig=c(0,1,0,0.5))
image(matrix(seq_along(sin_ord)), col = ifelse(sin_ord %in% t_cells, sin_cl@colors[ sin_cl@mc[ sin_ord]], NA), axes=F)
par(fig=c(0,1,0.5,1), new=T)
image(matrix(seq_along(sin_ord)), col = ifelse(sin_ord %in% dc_cells, sin_cl@colors[ sin_cl@mc[ sin_ord]], NA), axes=F)
dev.off()

png(paste0(outdir, "/Fig3j.png"), height=2000, width=2000)
db_clusts = interaction(factor(parser_dc, levels = lin_ord),
        factor(parser_t, levels = lin_ord))
db_vec = as.vector(db_clusts); names(db_vec) = good_pics
par(mar = c(0.5,3,0.5,0), fig = c(0,1,0,0.4), lwd = 2)
db_ord = plot_sc_heatmap(id_d, id, t_nms, clusts = db_vec, cells = good_pics, good_clusts = names(table(db_clusts)),
        annotate=F, normalize=T, lty=1, lwd=2); box()
par(fig = c(0,1,0.4,0.8), new = T)
db_ord = plot_sc_heatmap(id_d, id, dc_nms, clusts = db_vec, cells = good_pics, good_clusts = names(table(db_clusts)),
        annotate=F, normalize=T, lty=1, lwd=2); box()
par(fig=c(0,1,0.85,0.9), new=T)
image(matrix(1 * (cell_stats[ db_ord, "treatment"] == "helminths")), axes = F, col = c("gray80", "gray20")); box()
par(fig=c(0,1,0.8,0.85), new=T)
image(matrix(1 * (cell_stats[ db_ord, "sorting.scheme"] == "Ag+ doublets")), axes = F, col = c("white", "green2")); box()
par(fig=c(0,1,0.9,1), new=T)
plot(log2(umicount[db_ord]), col = rgb(0.2,0.2,0.2,0.4), pch = 20, xaxs = "i", axes = F, xlab = "", ylab = "", cex = 1.5, ylim = ylim); 
axis(2, las = 2); box()
dev.off()

png(paste0(outdir, "/Fig3k.png"), height=200, width=2000)
par(mar = c(0.5,3,0.5,0), fig=c(0,1,0,0.5))
image(matrix(seq_along(db_ord)), col = name2color[ parser_t[ db_ord]], axes=F)
par(fig=c(0,1,0.5,1), new=T)
image(matrix(seq_along(db_ord)), col = name2color[ parser_dc[ db_ord]], axes=F)
dev.off()

png(paste0(outdir, "/Fig3l.png"), height=200, width=2000)
par(mar = c(0.5,3,0.5,0))
split_count = cbind(alpha[good_pics], 1 - alpha[good_pics])# * umicount[good_pics]
#split_count = split_count / rowSums(split_count)
barplot(t(split_count[db_ord,]), col = c(t_col, dc_col), border = NA, xaxs = "i", space = 0, names.arg = rep("", length(db_ord)), las = 2)
box()
dev.off()

############

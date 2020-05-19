#################################
#                               #
# Reproduce Supp Figures 6-8    #
#                               #
#################################

message("Generating supplementary figures 6-8")

umis_n = sweep(umis, 2, colSums(umis), "/") * 1000
dc_pop = "MigDC"
sub_cells = names(which(parser_dc == dc_pop))
sub_cells = sub_cells[ parser_t[ sub_cells] %in% names(which(table( parser_t[sub_cells]) > 50))]
pops = names(table(parser_t[sub_cells]))
sin_cells = t_cells[ color2name[ sin_cl@colors[ sin_cl@mc[t_cells]]] %in% pops]
dir.create(paste0(supdir, "/FigS7f")) 
reg = 20
for (t_pop in pops) {
        sin_g1 = sin_cells[ sin_cl@mc[ sin_cells] %in% which(sin_cl@colors == name2color[t_pop])]
        sin_g2 = setdiff(sin_cells, sin_g1)
        db_g1 = names(which(parser_t[ sub_cells] == t_pop))
        db_g2 = setdiff(sub_cells, db_g1)
        sin_x = rowSums(umis_n[,sin_g2]) / length(sin_g2) * min(length(sin_g1), length(sin_g2))
        sin_y = rowSums(umis_n[,sin_g1]) / length(sin_g1) * min(length(sin_g1), length(sin_g2))
        db_x  = rowSums(umis_n[,db_g2]) / length(db_g2) * min(length(db_g1), length(db_g2))
        db_y  = rowSums(umis_n[,db_g1]) / length(db_g1) * min(length(db_g1), length(db_g2))
        sin_z = log2((reg + sin_y) / (reg + sin_x))
        db_z = log2((reg + db_y) / (reg + db_x))
        xlim = max(abs(sin_z)) * c(-1,1); ylim = max(abs(db_z)) * c(-1,1)
        png(paste0(supdir, "/FigS7f/", t_pop, ".png"), height=1000, width=1000)
        plot(sin_z, db_z, pch = 20, col = rgb(0,0,0.5,0.5), xlim = xlim, ylim = ylim, cex = 2.5, axes = F, xlab = "", ylab = "")
        axis(1); axis(2)
        grid(col = "black"); abline(h=0,v=0,lwd=2)
        dev.off()
        png(paste0(supdir, "/FigS7f/", t_pop, "_text.png"), height=1000, width=1000)
        plot(sin_z, db_z, type = "n", xlim = xlim, ylim = ylim, axes = F)
        axis(1); axis(2)
        grid(col = "black"); abline(h=0,v=0,lwd=2)
        text(sin_z, db_z, names(sin_z), col = rgb(0,0,0.5,0.5))
        dev.off()
}

################

sin_vec = sin_names[cells]
sin_ord = cells[ order(factor(sin_vec, levels = lin_ord))]

cord = unique(sin_cl@mc[ sin_ord])
batch_comb = unique(with(sin_stats[sin_ord,], paste0(amp_batch_id, "!", sorting.scheme, "@", treatment, "@", replicate)))
B = vecsplit(batch_comb, "!", 1); batch_comb = vecsplit(batch_comb, "!", 2); names(batch_comb) = B
comb_sort = factor(vecsplit(batch_comb, "@", 1),  c("Trbc+", "CD160+ Trbc+", "Icos+ Trbc+", "Tigit+ Trbc+",
                "Cd11c+", "Ag+ Cd11c+", "doublets", "Ag+ doublets"))
comb_treat   = factor(vecsplit(batch_comb, "@", 2), levels = c("PBS", "helminths"))
comb_rep   = factor(vecsplit(batch_comb, "@", 3))
bord = B[ order(comb_sort, comb_treat, comb_rep)]

cell_comb = with(sin_stats[sin_ord,], paste0(amp_batch_id, "!", sorting.scheme, "@", treatment, "@", replicate))
cell_comb = factor(cell_comb, levels = sapply(bord, grep, names(table(cell_comb)), v=T))
names(cell_comb) = sin_ord
create_batch_matrix(id_s, id_s, cell_comb, paste0(supdir, "/FigS6e.png"), cord = cord)

png(paste0(supdir, "/FigS6e_colorbar.png"), height=100, width=1000)
par(mar = rep(0,4))
image(matrix(seq_along(cord)), col = sin_cl@colors[cord])
dev.off()

##############

comb = with(cell_stats, paste0(sorting.scheme, ".", treatment)); names(comb) = rownames(cell_stats)
sample_dist = table(comb[t_cells], color2name[ sin_cl@colors[ sin_cl@mc[t_cells]]])
sample_dist = sample_dist[rev(c("Trbc+.PBS", "Trbc+.helminths", "Icos+ Trbc+.helminths", "CD160+ Trbc+.helminths")),
        intersect(lin_ord, colnames(sample_dist))]
t_n = sample_dist / rowSums(sample_dist)

sample_dist = table(comb[dc_cells], color2name[ sin_cl@colors[ sin_cl@mc[dc_cells]]])
sample_dist = sample_dist[rev(c("Cd11c+.PBS", "Cd11c+.helminths", "Ag+ Cd11c+.helminths")),
        intersect(lin_ord, colnames(sample_dist))]
dc_n = sample_dist / rowSums(sample_dist)

png(paste0(supdir, "/FigS6f.png"), height=1400, width=1000)
par(lwd=4, mfrow=c(2,1))
barplot(t(t_n), horiz = T, col = name2color[ colnames(t_n)])

barplot(t(dc_n), horiz = T, col = name2color[ colnames(dc_n)])
dev.off()

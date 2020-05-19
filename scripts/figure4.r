#########################
#                       #
# Reproduce Figure 4    #
#                       #
#########################

message("Generating Figure 4")

#######
# population analysis

good_pics = rownames(mle_res)
alpha = mle_res[good_pics, "alpha"]; names(alpha) = good_pics
t_mc = mle_res[good_pics, "a_mc"]; names(t_mc) = good_pics
dc_mc = mle_res[good_pics, "b_mc"]; names(dc_mc) = good_pics
parser_t = color2name[ sin_cl@colors[ t_mc]]; names(parser_t) = good_pics
parser_dc = color2name[ sin_cl@colors[ dc_mc]]; names(parser_dc) = good_pics

outdir = "figures/figure4"
dir.create(outdir)

bad_cells = c()
t_cells = setdiff(t_cells, bad_cells)
dc_cells = setdiff(dc_cells, bad_cells)

comb = with(cell_stats, paste0(sorting.scheme, ".", treatment, "@", timepoint, ".", replicate)); names(comb) = rownames(cell_stats)
t_dist = rbind(table(comb[good_pics], color2name[ sin_cl@colors[t_mc[good_pics]]]),
	table(comb[t_cells], color2name[ sin_cl@colors[ sin_cl@mc[t_cells]]]))
t_dist = t_dist[,intersect(lin_ord, colnames(t_dist))]
t_dist = t_dist[ rowSums(t_dist) > 20,]
t_dist = t_dist[ order(factor(vecsplit(rownames(t_dist), "\\.", 2), levels = c("helminths", "PBS")),
        factor(vecsplit(rownames(t_dist),"\\.", 1), levels = c("Trbc+", "Cd11c+", "doublets",  "Ag+ Cd11c+", "Ag+ doublets"))),]

dc_dist = rbind(table(comb[good_pics], color2name[ sin_cl@colors[dc_mc[good_pics]]]),
	table(comb[dc_cells], color2name[ sin_cl@colors[ sin_cl@mc[dc_cells]]]))
dc_dist = dc_dist[,intersect(lin_ord, colnames(dc_dist))]
dc_dist = dc_dist[ rowSums(dc_dist) > 20,]
dc_dist = dc_dist[ order(factor(vecsplit(rownames(dc_dist), "\\.", 2), levels = c("helminths", "PBS")),
        factor(vecsplit(rownames(dc_dist),"\\.", 1), levels = c("Trbc+", "Cd11c+", "doublets", "Ag+ Cd11c+", "Ag+ doublets"))),]

t_n = t_dist / rowSums(t_dist); dc_n = dc_dist / rowSums(dc_dist)

t_melt = melt(t_n)
t_melt$tp = vecsplit(as.vector(t_melt$Var1), "\\.", 2)
t_melt$group = vecsplit(as.vector(t_melt$Var1), "\\.", 1)
t_melt$rep = vecsplit(as.vector(t_melt$Var1), "\\.", 3)

dc_melt = melt(dc_n)
dc_melt$tp = vecsplit(as.vector(dc_melt$Var1), "\\.", 2)
dc_melt$group = vecsplit(as.vector(dc_melt$Var1), "\\.", 1)
dc_melt$rep = vecsplit(as.vector(dc_melt$Var1), "\\.", 3)

##########
# enrichment graph

cells = union(t_cells, dc_cells)
db_cells = intersect(good_pics, rownames(cell_stats)[ cell_stats$sorting.scheme == "doublets" & cell_stats$timepoint == "48h"])
sin_cells = intersect(cells, rownames(cell_stats)[ cell_stats$timepoint == "48h" & cell_stats$sorting.scheme %in% c("Trbc+", "Cd11c+")])

nb_cells =  intersect(sin_cells, rownames(sin_stats)[ sin_stats$treatment == "helminths"])
pbs_cells = intersect(sin_cells, rownames(sin_stats)[ sin_stats$treatment == "PBS"])

nb_doublets = intersect(db_cells, rownames(cell_stats)[ cell_stats$treatment == "helminths"])
pbs_doublets = intersect(db_cells, rownames(cell_stats)[ cell_stats$treatment == "PBS"])

pbs_obs = c(table(factor(t_mc[pbs_doublets], levels = t_clusts)),
         table(factor(dc_mc[pbs_doublets], levels = dc_clusts)))
pbs_exp = c(table(factor(sin_cl@mc[intersect(t_cells, pbs_cells)], levels = t_clusts)),
        table(factor(sin_cl@mc[intersect(dc_cells, pbs_cells)], levels = dc_clusts)))
pbs_obs_n = pbs_obs / sum(pbs_obs)
pbs_exp_n = pbs_exp / sum(pbs_exp)
reg = 0.01
pbs_enr = log2((reg + pbs_obs_n) / (reg + pbs_exp_n))

nb_obs = c(table(factor(t_mc[nb_doublets], levels = t_clusts)),
         table(factor(dc_mc[nb_doublets], levels = dc_clusts)))
nb_exp = c(table(factor(sin_cl@mc[intersect(t_cells, nb_cells)], levels = t_clusts)),
        table(factor(sin_cl@mc[intersect(dc_cells, nb_cells)], levels = dc_clusts)))
nb_obs_n = nb_obs / sum(nb_obs)
nb_exp_n = nb_exp / sum(nb_exp)
reg = 0.01
nb_enr = log2((reg + nb_obs_n) / (reg + nb_exp_n))

clust_ord = names(pbs_enr)[ order(factor(color2name[ sin_cl@colors[ as.numeric(names(pbs_enr))]], levels = lin_ord), rowMeans(cbind(pbs_enr, nb_enr)))]
t_ord = intersect(clust_ord, names(table(sin_cl@mc[t_cells])))
t_ord = intersect(t_ord, names(which(pbs_exp[t_ord] > 5 & nb_exp[t_ord] > 5 & pbs_obs[t_ord] > 5 & nb_obs[t_ord] > 5)))
dc_ord = intersect(clust_ord, names(table(sin_cl@mc[dc_cells])))
dc_ord = intersect(dc_ord, names(which(pbs_exp[dc_ord] > 5 & nb_exp[dc_ord] > 5 & pbs_obs[dc_ord] > 5 & nb_obs[dc_ord] > 5)))
ord = c(t_ord, dc_ord)
t_share = length(t_ord) / length(ord)
IM = cbind(pbs_enr[ord], nb_enr[ord])
ylim = max(abs(IM)) * c(-1,1)

png(paste0(outdir, "/Fig4a.png"), height = 700, width = 1000)
barplot(pbs_enr[t_ord], border = NA, col = sin_cl@colors[ as.numeric(t_ord)], ylim = ylim, names.arg = rep("", length(t_ord)))
dev.off()

png(paste0(outdir, "/Fig4b.png"), height = 700, width = 1000)
barplot(pbs_enr[dc_ord], border = NA, col = sin_cl@colors[ as.numeric(dc_ord)], ylim = ylim, names.arg = rep("", length(dc_ord)))
dev.off()

png(paste0(outdir, "/Fig4c.png"), height = 700, width = 1000)
barplot(nb_enr[t_ord], border = NA, col = sin_cl@colors[ as.numeric(t_ord)], ylim = ylim, names.arg = rep("", length(t_ord)))
dev.off()

png(paste0(outdir, "/Fig4d.png"), height = 700, width = 1000)
barplot(nb_enr[dc_ord], border = NA, col = sin_cl@colors[ as.numeric(dc_ord)], ylim = ylim, names.arg = rep("", length(dc_ord)))
dev.off()

#############
# T cells

library(Hmisc)
t_melt$abs = with(t_melt, t_dist[cbind(as.vector(Var1), as.vector(Var2))])
t_melt$bar_id = with(t_melt, paste0(Var2, ":", group, "#", tp))
m = with(t_melt, tapply(abs, bar_id, sum))
tot_df = data.frame(m = m, bar_id = vecsplit(names(m), ":", 2), pop = vecsplit(names(m), ":", 1))
tot_df$bar_id = as.vector(tot_df$bar_id)
tot_df$group = vecsplit(tot_df$bar_id, "#", 1); tot_df$tp = vecsplit(tot_df$bar_id, "#", 2)
totals = with(tot_df, tapply( m, bar_id, sum))
tot_df$n = totals[ tot_df$bar_id]

treats = c("PBS@48h", "helminths@48h"); groups = c("Trbc+", "doublets")
tot_df = tot_df[ vecsplit(as.vector(tot_df$bar_id), "#", 1) %in% groups & vecsplit(as.vector(tot_df$bar_id), "#", 2) %in% treats,]

conf = with(tot_df, binconf(m,n))
rownames(conf) = rownames(tot_df)
#Y = sweep(t(apply(m,1,binconf,n)), 2, rep(tapply(colSums(ds_f[,cells2]), db_int[cells2], mean), 3), "*")
tot_df = cbind(tot_df, conf)

ylim = c(0,0.63)
X = t_melt[t_melt$tp %in% treats & t_melt$group %in% groups,]
Y = dcast(tot_df[,c("tp", "group", "pop", "PointEst")], factor(group, levels = groups) ~ factor(pop, levels = lin_ord) * factor(tp, levels = treats), mean)
Y.m = as.matrix(Y[,-1]); rownames(Y.m) = Y[,1]
Y = dcast(tot_df[,c("tp", "group", "pop", "Lower")], factor(group, levels = groups) ~ factor(pop, levels = lin_ord) * factor(tp, levels = treats), mean)
Y.l = as.matrix(Y[,-1]); rownames(Y.l) = Y[,1]
Y = dcast(tot_df[,c("tp", "group", "pop", "Upper")], factor(group, levels = groups) ~ factor(pop, levels = lin_ord) * factor(tp, levels = treats), mean)
Y.u = as.matrix(Y[,-1]); rownames(Y.u) = Y[,1]

fill_cols = rep(c("white", "slateblue2"),2)
X2 = X[ X$group %in% rownames(Y.m) & X$tp %in% treats,]

Z = barplot(Y.m, beside=T, las = 2, col = fill_cols, space = rep(c(1,0.2,1,0.2), ncol(Y.m) / 2), ylim = ylim) #c(0, max(X2$value) * 1.05))
dimnames(Z) = dimnames(Y.m)
X2$title = paste0(X2$Var2, "_", X2$tp)
X2$coord = Z[cbind(X2$group, X2$title)]
X2$offset = X2$coord + runif(nrow(X2), -0.4, 0.4)
X2$side = with(X2, paste0(tp, "#", rep))
two_sides = names(which(table(X2$side) == ncol(Y.m)))
X3 = X2[ X2$side %in% two_sides,]
X3 = X3[ order(X3$group, X3$Var2, X3$side),]
X3_singlets = X3[ X3$group == "Trbc+",]; X3_doublets = X3[X3$group == "doublets",]

paired_pvals = rep(NA, nrow(X3_singlets));
for (i in seq_along(paired_pvals)) {
        pop = as.vector(X3_singlets[i, "Var2"])
        sin_cells = names(which(comb == as.vector(X3_singlets[i, "Var1"])))
        db_cells = names(which(comb == as.vector(X3_doublets[i, "Var1"])))
        dist = rbind(table(color2name[ sin_cl@colors[ sin_cl@mc[ sin_cells]]] == pop), table(parser_t[db_cells] == pop))
        paired_pvals[i] = fisher.test(dist)$p.value
}

paired_qvals = p.adjust(paired_pvals, "fdr")
line_col = ifelse(paired_qvals < 0.05, ifelse(X3_doublets$value > X3_singlets$value, "red2", "blue2"), "gray20")
png(paste0(supdir, "/FigS8a.png"), height=700, width = 1800)
Z = barplot2(Y.m, beside=T, las = 2, col = "gray80", space = rep(c(1,0.2,1,0.2), ncol(Y.m) / 2), ylim = ylim, axes = F, border = NA,
	plot.ci = T, ci.l = Y.l, ci.u = Y.u, ci.lwd=4)
axis(2); axis(1, at = colMeans(Z), labels = colnames(Z), las = 2)
segments(X3_singlets$offset, X3_singlets$value, X3_doublets$offset, X3_doublets$value, col = line_col, lwd = ifelse(line_col != "gray20", 3, 1.5))
with(X2, points(offset, value, pch = ifelse(group == "doublets", 21, 23), cex = 3, bg = name2color[ as.vector(X2$Var2)]))
dev.off()

#######
# DC

dc_melt$abs = with(dc_melt, dc_dist[cbind(as.vector(Var1), as.vector(Var2))])
dc_melt$bar_id = with(dc_melt, paste0(Var2, ":", group, "#", tp))
m = with(dc_melt, tapply(abs, bar_id, sum))
tot_df = data.frame(m = m, bar_id = vecsplit(names(m), ":", 2), pop = vecsplit(names(m), ":", 1))
tot_df$bar_id = as.vector(tot_df$bar_id)
tot_df$group = vecsplit(tot_df$bar_id, "#", 1); tot_df$tp = vecsplit(tot_df$bar_id, "#", 2)
totals = with(tot_df, tapply( m, bar_id, sum))
tot_df$n = totals[ tot_df$bar_id]

treats = c("PBS@48h", "helminths@48h"); groups = c("Cd11c+", "doublets") #, "Ag+ Cd11c+", "Ag+ doublets")
tot_df = tot_df[ vecsplit(as.vector(tot_df$bar_id), "#", 1) %in% groups & vecsplit(as.vector(tot_df$bar_id), "#", 2) %in% treats,]

conf = with(tot_df, binconf(m,n))
rownames(conf) = rownames(tot_df)
tot_df = cbind(tot_df, conf)

ylim = c(0,0.63)
X = dc_melt[dc_melt$tp %in% treats & dc_melt$group %in% groups,]
Y = dcast(tot_df[,c("tp", "group", "pop", "PointEst")], factor(group, levels = groups) ~ factor(pop, levels = lin_ord) * factor(tp, levels = treats), mean)
Y.m = as.matrix(Y[,-1]); rownames(Y.m) = Y[,1]
Y = dcast(tot_df[,c("tp", "group", "pop", "Lower")], factor(group, levels = groups) ~ factor(pop, levels = lin_ord) * factor(tp, levels = treats), mean)
Y.l = as.matrix(Y[,-1]); rownames(Y.l) = Y[,1]
Y = dcast(tot_df[,c("tp", "group", "pop", "Upper")], factor(group, levels = groups) ~ factor(pop, levels = lin_ord) * factor(tp, levels = treats), mean)
Y.u = as.matrix(Y[,-1]); rownames(Y.u) = Y[,1]

X2 = X[ X$group %in% rownames(Y.m) & X$tp %in% treats,]
Z = barplot(Y.m, beside=T, las = 2, col = fill_cols, space = rep(c(1,0.2,1,0.2), ncol(Y.m) / 2), ylim = ylim) #c(0, max(X2$value) * 1.05))
dimnames(Z) = dimnames(Y.m)
X2$coord = Z[cbind(X2$group, paste0(X2$Var2, "_", X2$tp))]
with(X2, points(coord + runif(nrow(X2), -0.2, 0.2), value, pch = 21, cex = 4, bg = name2color[ as.vector(X2$Var2)]))

X2$offset = X2$coord + runif(nrow(X2), -0.4, 0.5)
X2$side = with(X2, paste0(tp, "#", rep))
two_sides = names(which(table(X2$side) == ncol(Y.m)))
X3 = X2[ X2$side %in% two_sides,]
X3 = X3[ order(X3$group, X3$Var2, X3$side),]
X3_singlets = X3[ X3$group == "Cd11c+",]; X3_doublets = X3[X3$group == "doublets",]

paired_pvals = rep(NA, nrow(X3_singlets)); 
for (i in seq_along(paired_pvals)) {
	pop = as.vector(X3_singlets[i, "Var2"])
	sin_cells = names(which(comb == as.vector(X3_singlets[i, "Var1"])))
	db_cells = names(which(comb == as.vector(X3_doublets[i, "Var1"])))
	dist = rbind(table(color2name[ sin_cl@colors[ sin_cl@mc[ sin_cells]]] == pop), table(parser_dc[db_cells] == pop))
	paired_pvals[i] = fisher.test(dist)$p.value
}

paired_qvals = p.adjust(paired_pvals, "fdr")
line_col = ifelse(paired_qvals < 0.05, ifelse(X3_doublets$value > X3_singlets$value, "red2", "blue2"), "gray20")
png(paste0(supdir, "/FigS8b.png"), height=700, width = 1800)
Z = barplot2(Y.m, beside=T, las = 2, col = "gray80", space = rep(c(1,0.2,1,0.2), ncol(Y.m) / 2), ylim = ylim, axes = F, border = NA,
	plot.ci = T, ci.l = Y.l, ci.u = Y.u, ci.lwd=4)
axis(2); axis(1, at = colMeans(Z), labels = colnames(Z), las = 2)
segments(X3_singlets$offset, X3_singlets$value, X3_doublets$offset, X3_doublets$value, col = line_col, lwd = ifelse(line_col != "gray20", 3, 1.5))
with(X2, points(offset, value, pch = ifelse(group == "doublets", 21, 23), cex = 3, bg = name2color[ as.vector(X2$Var2)]))
dev.off()

############
# Joint and Antigen

treats = c("helminths@48h"); groups = c("doublets", "Ag+ doublets")
X = t_melt[t_melt$tp %in% treats & t_melt$group %in% groups,]
Y = dcast(X[,c(5,2,3)], factor(group, levels = groups) ~ factor(Var2, levels = lin_ord), mean)
Y2 = as.matrix(Y[,-1]); rownames(Y2) = Y[,1]
fill_cols = rep(c("gray70", "gray30"),2)
X2 = X[ X$group %in% rownames(Y2) & X$tp %in% treats,]
png(paste0(supdir, "/FigS8d.png"), height=700, width=1500)
par(lwd=6)
Z = barplot(Y2, beside=T, las = 1, col = fill_cols, ylim = c(0, max(X2$value) * 1.05), space = rep(c(2,0), ncol(Y2)))
dimnames(Z) = dimnames(Y2)
X2$coord = Z[cbind(X2$group, as.vector(X2$Var2))]
with(X2, points(coord + runif(nrow(X2), -0.2, 0.2), value, pch = 21, cex = 4, bg = name2color[ as.vector(X2$Var2)]))
dev.off()

#################
# gene analysis

numis = 1000
ds = .downsamp(umis[,good_pics], numis)
ds_f = ds[setdiff(rownames(sin_cl@e_gc), bad_genes), good_pics]
us = umis[ rownames(ds_f), colnames(ds_f)]

exp_us = generate_expected_pics_from_mle(id_s, mle_res[good_pics, c("a_mc", "b_mc")], mle_res[ good_pics, "alpha"],
        colSums(us), bad_genes = bad_genes)
exp_n = generate_expected_pics_from_mle(id_s, mle_res[good_pics, c("a_mc", "b_mc")], mle_res[ good_pics, "alpha"],
	colSums(ds_f), bad_genes = bad_genes)
t_n = generate_expected_pics_from_mle(id_s, mle_res[good_pics, c("a_mc", "b_mc")], rep(1, length(good_pics)),
        colSums(ds_f) * alpha, bad_genes = bad_genes)
dc_n = generate_expected_pics_from_mle(id_s, mle_res[good_pics, c("a_mc", "b_mc")], rep(0, length(good_pics)),
	colSums(ds_f) * (1 - alpha), bad_genes = bad_genes)
genes = rownames(exp_us)

###############

sub_pics = intersect(good_pics, rownames(cell_stats)[ cell_stats$sorting.scheme == "doublets" & cell_stats$timepoint == "48h"])

comb = with(cell_stats, paste0(sorting.scheme, ".", treatment, ".", timepoint, ".", replicate)); names(comb) = rownames(cell_stats)
pic_joint = interaction(factor(color2name[sin_cl@colors[t_mc[sub_pics]]], levels = lin_ord),
        factor(color2name[sin_cl@colors[dc_mc[sub_pics]]], levels = lin_ord))
names(pic_joint) = sub_pics

t_ord = intersect(lin_ord, names(table(parser_t))); dc_ord = intersect(lin_ord, names(table(parser_dc)))
pop = "Treg"; prefix = "treg"; side = 2
cell_thresh = 10
analyzed_pics = intersect(sub_pics, names(which(parser_t == pop)))
pic_comb = factor(pic_joint[analyzed_pics]); 
good_joint = names(which(rowSums(table(pic_joint[analyzed_pics], as.vector(cell_stats[analyzed_pics, "treatment"])) < 15) == 0))
names(pic_comb) = analyzed_pics
real_m = t(apply(us[,analyzed_pics], 1, tapply, pic_comb[analyzed_pics], sum))
exp_m =  t(apply(exp_us[,analyzed_pics], 1, tapply, pic_comb[analyzed_pics], sum))
x2 = rowSums(((real_m - exp_m)^2 )/exp_m)
qvals = p.adjust(1 - pchisq(x2, df = ncol(real_m) - 1), "fdr")
z = log2(real_m / exp_m); 
z[real_m==0] = NA

y = log2(rowSums(us[names(qvals),analyzed_pics]));
genes = names(which(qvals < 1e-6 & apply(us[ names(qvals), analyzed_pics], 1, function(x) sort(x,T)[3]) > 3 & y > 6))

z_reg = log2((real_m + 5) / (exp_m + 5))
IM = z_reg[genes, ]#grep(paste0("^|\\.", pop, "\\.|$"), colnames(z))]
IM = IM[rowSums(!is.na(IM)) > 0,]

x_min = apply(IM,1,min); x_max = apply(IM,1,max)
x = ifelse(abs(x_min) > x_max, x_min, x_max)

library(tglkmeans)

k = 15
data = as.data.frame(IM)
data$id = rownames(data)
data = data[,c(ncol(data), 1:(ncol(data) - 1))]
km <- TGL_kmeans_tidy(data, k=k, metric='euclid', verbose=TRUE, seed = 18)
centers = as.matrix(km$centers[,-1]); rownames(centers) = seq_len(k)
x_min = apply(centers,1,min); x_max = apply(centers,1,max)
centers = centers[ order(abs(x_min) < x_max, max.col(centers)),]
km_clusts = as.numeric(factor(km$cluster$clust, levels = rownames(centers))); names(km_clusts) = rownames(IM)
ord = order(km_clusts)
cls = cumsum(table(km_clusts)) / length(km_clusts)

png(paste0(supdir, "/FigS8c.png"), height = nrow(IM) * 12, width=1000)
par(mar = c(10,15,5,5))
image.2(IM, balance = T, annotate = "both", hct = km_clusts)
dev.off()

##############

marker_genes = c("Il12b", "Ccl6")
real_m = t(apply(ds_f[marker_genes,analyzed_pics], 1, tapply, pic_comb[analyzed_pics], mean)); real_m = real_m[,colnames(IM)]
exp_m =  t(apply(exp_n[marker_genes,analyzed_pics], 1, tapply, pic_comb[analyzed_pics], mean)); exp_m = exp_m[,colnames(IM)]
t_m =  t(apply(t_n[marker_genes,analyzed_pics], 1, tapply, pic_comb[analyzed_pics], mean)); t_m = t_m[,colnames(IM)]
dc_m =  t(apply(dc_n[marker_genes,analyzed_pics], 1, tapply, pic_comb[analyzed_pics], mean)); dc_m = dc_m[,colnames(IM)]

library(Hmisc)
m = t(apply(ds_f[marker_genes, analyzed_pics],1,tapply, pic_comb[analyzed_pics], sum))
n = tapply(colSums(ds_f[,analyzed_pics]),pic_comb[analyzed_pics],sum)
Y = sweep(t(apply(m,1,binconf,n)), 2, rep(tapply(colSums(ds_f[,analyzed_pics]), pic_comb[analyzed_pics], mean), 3), "*")
ord = seq_along(n)

png(paste0(outdir, "/FigS4g.png"), height=700, width=2000)
par(mfrow = c(1,2))
for (gene in marker_genes) {
        exp_tab = rbind(t_m[gene,ord], dc_m[gene,ord], 0, 0)
        real_tab = rbind(0, 0, real_m[gene,ord], 0)
        tab = cbind(exp_tab, real_tab)[, rep(seq_along(ord), each = 2) + rep(c(0,length(ord)), length(ord))]
        mtab = max(Y[gene,], na.rm=T)
        X = barplot(tab, col = c("limegreen", "firebrick3", "gray70", "white"), xaxs = "i", las = 2, space = c(1,0), ylim = c(0,mtab), main = gene, cex.main=2)
        obs_coords = X[seq(2,length(X),2)]; ci.l = Y[gene, ord + length(n)]; ci.u = Y[gene, ord + 2 * length(n)]
        segments(obs_coords, ci.l, y1 = ci.u);
        segments(obs_coords-0.2, ci.l, x1 = obs_coords + 0.2); segments(obs_coords-0.2, ci.u, x1 = obs_coords + 0.2);
}
dev.off()

########################################
# analyze pathogen-specific interactions

dc_pop = "MigDC";
sub_pics = intersect(good_pics, rownames(cell_stats)[ cell_stats$timepoint == "48h" & cell_stats$sorting.scheme %in% c("doublets", "Ag+ doublets") & 
	cell_stats$treatment %in% c("PBS","helminths")])
analyzed_pics = intersect(sub_pics, names(parser_dc)[parser_dc == dc_pop])
pic_joint = factor(paste0(cell_stats[analyzed_pics, "treatment"], "@", cell_stats[analyzed_pics, "sorting.scheme"]), 
	levels = c("PBS@doublets", "helminths@doublets", "helminths@Ag+ doublets"))
names(pic_joint) = analyzed_pics

good_joint = names(which(rowSums(table(pic_joint[analyzed_pics], as.vector(cell_stats[analyzed_pics, "treatment"])) < 15) == 0))
pic_comb = pic_joint
names(pic_comb) = analyzed_pics
real_m = t(apply(us[,analyzed_pics], 1, tapply, pic_comb[analyzed_pics], sum))
exp_m =  t(apply(exp_us[,analyzed_pics], 1, tapply, pic_comb[analyzed_pics], sum))
x2 = rowSums(((real_m - exp_m)^2 )/exp_m)
qvals = p.adjust(1 - pchisq(x2, df = ncol(real_m) - 1), "fdr")
z = log2(real_m / exp_m); 
z[real_m==0] = NA

y = log2(rowSums(us[names(qvals),analyzed_pics]));
genes = union(names(which(qvals < 1e-5 & apply(us[ names(qvals), analyzed_pics], 1, function(x) sort(x,T)[3]) > 3 & y > 5)), c())

z_reg = log2((real_m + 5) / (exp_m + 5))
IM = z_reg[genes, ]#grep(paste0("^|\\.", pop, "\\.|$"), colnames(z))]
IM = IM[rowSums(!is.na(IM)) > 0,]

x_min = apply(IM,1,min); x_max = apply(IM,1,max)
x = ifelse(abs(x_min) > x_max, x_min, x_max)

library(tglkmeans)

k = 12
data = as.data.frame(IM)
data$id = rownames(data)
data = data[,c(ncol(data), 1:(ncol(data) - 1))]
km <- TGL_kmeans_tidy(data, k=k, metric='euclid', verbose=TRUE, seed = 18)
#gc = km$cluster$clust; names(gc) = km$cluster$id
#cls = cumsum(table(gc)) / length(gc)
centers = as.matrix(km$centers[,-1]); rownames(centers) = seq_len(k)
x_min = apply(centers,1,min); x_max = apply(centers,1,max)
centers = centers[ order(abs(x_min) < x_max, max.col(centers)),]
km_clusts = as.numeric(factor(km$cluster$clust, levels = rownames(centers))); names(km_clusts) = rownames(IM)
ord = order(km_clusts)
cls = cumsum(table(km_clusts)) / length(km_clusts)

png(paste0(supdir, "/FigS8e.png"), height = nrow(IM) * 12, width=1000)
par(mar = c(5,15,3,3))
image.2(IM, balance = T, annotate = "both", hct = km_clusts)
dev.off()

####################

genes = c("Ccl22", "Ccl17", "Cd40", "Ebi3", "Dll4")
real_m = t(apply(ds_f[genes,analyzed_pics], 1, tapply, pic_joint[analyzed_pics], mean))
exp_m =  t(apply(exp_n[genes,analyzed_pics], 1, tapply, pic_joint[analyzed_pics], mean))
t_m =  t(apply(t_n[genes,analyzed_pics], 1, tapply, pic_joint[analyzed_pics], mean))
dc_m =  t(apply(dc_n[genes,analyzed_pics], 1, tapply, pic_joint[analyzed_pics], mean))

library(Hmisc)
m = t(apply(ds_f[genes, analyzed_pics],1,tapply, pic_comb[analyzed_pics], sum))
n = tapply(colSums(ds_f[,analyzed_pics]),pic_comb[analyzed_pics],sum)
Y = sweep(t(apply(m,1,binconf,n)), 2, rep(tapply(colSums(ds_f[,analyzed_pics]), pic_comb[analyzed_pics], mean), 3), "*")
ord = seq_along(n)

png(paste0(outdir, "/Fig4h.png"), height=700, width=2500)
par(lwd=2, mfrow = c(1,length(genes)))
for (gene in genes) {
        exp_tab = rbind(t_m[gene,], dc_m[gene,], 0, 0)
        real_tab = rbind(0, 0, real_m[gene,], 0)
        tab = cbind(exp_tab, real_tab)[, rep(seq_len(ncol(real_m)), each = 2) + rep(c(0,ncol(real_m)), ncol(real_m))]
        mtab = max(Y[gene,], na.rm=T)
        X = barplot(tab, col = c("limegreen", "firebrick3", "gray70", "white"), xaxs = "i", las = 2, space = c(1,0), ylim = c(0,mtab), main = gene, cex.main=2)
        obs_coords = X[seq(2,length(X),2)]; ci.l = Y[gene, ord + length(n)]; ci.u = Y[gene, ord + 2 * length(n)]
        segments(obs_coords, ci.l, y1 = ci.u);
	segments(obs_coords-0.2, ci.l, x1 = obs_coords + 0.2); segments(obs_coords-0.2, ci.u, x1 = obs_coords + 0.2);
}
dev.off()

#########

genes = c("Serpinb9", "Serpinb9b", "Fabp4", "Fabp5", "Ccl9", "Nrp2",
	"Bcl2a1d", "Bcl2l1", "Tnfsf9", "Cd86")

real_m = t(apply(ds_f[genes,analyzed_pics], 1, tapply, pic_joint[analyzed_pics], mean))
exp_m =  t(apply(exp_n[genes,analyzed_pics], 1, tapply, pic_joint[analyzed_pics], mean))
t_m =  t(apply(t_n[genes,analyzed_pics], 1, tapply, pic_joint[analyzed_pics], mean))
dc_m =  t(apply(dc_n[genes,analyzed_pics], 1, tapply, pic_joint[analyzed_pics], mean))

library(Hmisc)
m = t(apply(ds_f[genes, analyzed_pics],1,tapply, pic_comb[analyzed_pics], sum))
n = tapply(colSums(ds_f[,analyzed_pics]),pic_comb[analyzed_pics],sum)
Y = sweep(t(apply(m,1,binconf,n)), 2, rep(tapply(colSums(ds_f[,analyzed_pics]), pic_comb[analyzed_pics], mean), 3), "*")
ord = seq_along(n)

dir.create(paste0(supdir, "/FigS8f/"))
for (gene in genes) {
	png(paste0(supdir, "/FigS8f/", gene, ".png"), height=700, width=1200)
	par(lwd=2)
        exp_tab = rbind(t_m[gene,], dc_m[gene,], 0, 0)
        real_tab = rbind(0, 0, real_m[gene,], 0)
	tab = cbind(exp_tab, real_tab)[, rep(seq_len(ncol(real_m)), each = 2) + rep(c(0,ncol(real_m)), ncol(real_m))]
        mtab = max(Y[gene,], na.rm=T)
        X = barplot(tab, col = c("limegreen", "firebrick3", "gray70", "white"), xaxs = "i", las = 2, space = c(1,0), ylim = c(0,mtab))
        obs_coords = X[seq(2,length(X),2)]; ci.l = Y[gene, ord + length(n)]; ci.u = Y[gene, ord + 2 * length(n)]
	segments(obs_coords, ci.l, y1 = ci.u);
        segments(obs_coords-0.2, ci.l, x1 = obs_coords + 0.2); segments(obs_coords-0.2, ci.u, x1 = obs_coords + 0.2);
	dev.off()
}

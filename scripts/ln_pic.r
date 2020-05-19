#################################################################
#                                                               #
# Run PIC-seq on in vivo infection model (Figures 3 & 4)	#
#                                                               #
#################################################################

message("Running PIC-seq on in vivo infection model (Figures 3 & 4)")
id = "ln"

############

id_s = paste0(id, "_singlets")
id_d = paste0(id, "_PIC")
sc_mat = scdb_mat(id)
sin_2d = scdb_mc2d(id_s); sin_cl = scdb_mc(id_s); sin_mat = scdb_mat(id_s)
db_2d = db_cl = scdb_mc(id_d); db_mat = scdb_mat(id_d)

all_cells = union(names(sin_cl@mc), names(db_cl@mc))
umis = as.matrix(sc_mat@mat[,all_cells])
cell_stats = sc_mat@cell_metadata[all_cells,]
fp = sin_cl@mc_fp
lfp = log2(sin_cl@mc_fp)
supdir = "figures/supp_figures_6-8/"
dir.create(supdir)
scfigs_init(supdir)

sin_stats = sin_mat@cell_metadata[names(sin_cl@mc),]
sin_comb = with(sin_stats, paste0(sorting.scheme, ".", treatment, ".", timepoint, ".", date))
names(sin_comb) = rownames(sin_stats)
db_stats = db_mat@cell_metadata[names(db_cl@mc),]
db_comb = with(db_stats, paste0(sorting.scheme, ".", treatment, ".", timepoint, ".", date))
names(db_comb) = rownames(db_stats)
sin_umis = as.matrix(sin_mat@mat[, names(sin_cl@mc)])
sin_n = sweep(sin_umis,2,colSums(sin_umis),"/") * 1000

color_scheme = sin_cl@color_key
color2name = color_scheme$group; names(color2name) = color_scheme$color
name2color = color_scheme$color; names(name2color) = color_scheme$group
sin_names = color2name[ sin_cl@colors[ sin_cl@mc]]; names(sin_names) = names(sin_cl@mc)

lin_ord = c("T", "CD8T", "CD8mem", "ActT", "Treg",
	"MigDC", "cDC1", "cDC2", "Monocytes", "pDC", "NK", "B")

t_cells = setdiff(names(sin_cl@mc)[ color2name[ sin_cl@colors[ sin_cl@mc]] %in% lin_ord[1:5]], c())
dc_cells = names(sin_cl@mc)[ color2name[ sin_cl@colors[ sin_cl@mc]] %in% lin_ord[6:10]]
cells = union(t_cells, dc_cells)

############
# Discard irrelevant PIC populations

db_col = "blue"
clusts = c(paste0("S.", sin_cl@mc), paste0("P.", db_cl@mc))
names(clusts) = c(names(sin_cl@mc), names(db_cl@mc))
num_clusts = as.numeric(factor(clusts)); names(num_clusts) = names(clusts)
cells = names(clusts)

afp = log2(mc_compute_fp_abs(num_clusts[cells], umis[,cells]))
clust_title = names(table(clusts))
colnames(afp) = clust_title
clust_num = as.numeric(vecsplit(clust_title, "\\.", 2)) #unlist(lapply(sapply(clust_title, strsplit, "\\."), "[[", 2)))
clust_orig = vecsplit(clust_title, "\\.", 1)

ga = "Cst3"; gb = "Trbc2"; a = afp[ga,]; b = afp[gb,];
gc = "Xcl1"; gd = "Igkc"; c = afp[gc,]; d = afp[gd,];
#ga = "Fcer1g"; gb = "Gzma"; a = afp[ga,]; b = afp[gb,]; 
db_clusts = grep("P", clust_title, v = T)
sin_clusts = grep("S", clust_title, v = T)
good_clusts = names(which(c < 2 & d < 3))
clust_cols = ifelse(clust_orig == "P", ifelse(clust_title %in% good_clusts, db_col, "gray60"), sin_cl@colors[ clust_num])

color2name[db_col] = "PIC"
name2color["PIC"] = db_col
good_clusts = names(which(clust_cols == db_col))
good_pics = names(clusts)[ clusts %in% good_clusts]

############

anchor_genes = c("Top2a", "Mki67", "Pcna", "Mcm4", "Cdk1", "Ube2c")
mcell_mat_rpt_cor_anchors(mat_id=id, gene_anchors = anchor_genes, cor_thresh = 0.1, gene_anti = c(),
        tab_fn = paste0(supdir, "/lateral_gmods.txt"), sz_cor_thresh = 0.2)
gcor_mat = read.delim(paste0(supdir, "/lateral_gmods.txt"), stringsAsFactor=F, row.names=1)
foc_genes = apply(gcor_mat[, gsub("-", ".", anchor_genes)], 1, which.max)
gset = gset_new_gset(sets = foc_genes, desc = "Cell cycle and MHC-II genes")
scdb_add_gset(paste0(id, "_lateral"), gset)

mcell_mat_ignore_genes(new_mat_id = paste0(id, "_lateral"), mat_id = id, ig_genes = names(foc_genes), reverse = T)
mcell_gset_split_by_dsmat(gset_id = paste0(id, "_lateral"), mat_id = paste0(id, "_lateral"), K = 7)
gset = scdb_gset(paste0(id, "_lateral"))
mcell_plot_gset_cor_mats(gset_id = paste0(id, "_lateral"), scmat_id = paste0(id, "_lateral"))
lateral_clusts = unique(gset@gene_set[ anchor_genes])
mcell_gset_remove_clusts(gset_id = paste0(id, "_lateral"), filt_clusts = lateral_clusts, new_id = paste0(id, "_lateral_f"), reverse=T)

############
# feature selection

bad_genes = grep("Gm[0-9].|Mir|-ps|Rpl|Rps|Ig|Jchain", rownames(umis), v=T)
cells = union(t_cells, dc_cells)
rel_cells = cells

lateral_set = scdb_gset(paste0(id, "_lateral_f"))@gene_set
cc_genes = setdiff(names(lateral_set)[ lateral_set %in% c(lateral_set["Top2a"], lateral_set["Ube2c"])], "Ldha")
cc_umis = colSums(umis[cc_genes,rel_cells])
cc_cells = names(which(cc_umis >= 16))

sub_t_cells = setdiff(intersect(rel_cells, t_cells), cc_cells)
sub_dc_cells = setdiff(intersect(rel_cells, dc_cells), cc_cells)

lr_features = choose_lr_features(id_s, sub_t_cells, sub_dc_cells, bad_genes=bad_genes, cor_n=100, must_haves=names(scdb_gset(id)@gene_set))
mle_features = choose_mle_features(id_s, id_s,  t_cells, dc_cells, bad_genes=c("Ftl1", bad_genes), shared_thresh=1,
	existing_list = names(scdb_gset(id)@gene_set))

##############
# Run PIC-seq

comb = paste0(cell_stats$treatment, ".", cell_stats$date); names(comb) = rownames(cell_stats)
names(comb) = rownames(cell_stats)

bad_cells = names(which(comb == "helminths.20190326"))

numis=1000
ds = .downsamp(umis[,good_pics], numis)
good_pics = intersect(good_pics, colnames(ds))

mle_res = run_pic_seq(id_s, id_s, ds[,good_pics], setdiff(t_cells, bad_cells), setdiff(dc_cells, bad_cells),
	lr_features, mle_features, fname=paste0(supdir, "/FigS7a.png"), bad_genes = bad_genes,
        comb = comb, reg = 1e-4, numis = 1000, downsample=F)

###############
# Analyze triplets

triplet_dir = paste0(supdir, "/FigS7d-e/")
dir.create(triplet_dir)

triplet_res = analyze_triplets(id_s, id_s, id_s, ds[,good_pics], setdiff(t_cells, bad_cells), setdiff(dc_cells, bad_cells),
	mle_res, mle_features, bad_genes=bad_genes, reg=1e-4, outdir=paste0(supdir, "/FigS7d-e/"), downsample=F)

forward_pr = triplet_res$tr_p
reverse_pr = triplet_res$rev_p
tr_res = triplet_res$mle_res

#############
# run mle on simulated PICs

numis = 1000; k = 5000
res = simulate_doublets(id, setdiff(t_cells, bad_cells), setdiff(dc_cells, bad_cells), k, comb, numis = rep(numis, k))

sim_umis = res$sim_umis; sim_info = res$info
sim_cells = names(which(colSums(sim_umis) == numis))
sim_umis = sim_umis[,sim_cells]; sim_info = sim_info[sim_cells,]

sim_alpha = sim_info$alpha.1; names(sim_alpha) = rownames(sim_info)
sim_mle_res = assign_pics_to_singlets(id_s, id_s, sim_umis, setdiff(t_cells, bad_cells), setdiff(dc_cells, bad_cells), sim_alpha,
	verbose=T, bad_genes = bad_genes, markers = mle_features, reg = 1e-4)

t_confu = table(sin_cl@mc[ as.vector( sim_info$sim.1)], sim_mle_res$a_mc)
t_n = t_confu / rowSums(t_confu)
dc_confu = table(sin_cl@mc[ as.vector( sim_info$sim.2)], sim_mle_res$b_mc)
dc_n = dc_confu / rowSums(dc_confu)

grad = colorRampPalette(c("white", "#FDC51D", "#CA531C", "#951851", "#36277A", "black"))(1000)
t_cls = factor(color2name[ sin_cl@colors[ as.numeric(rownames(t_n))]], levels = lin_ord)
png(paste0(supdir, "/FigS7b.png"), height = 1500, width = 1500)
par(mar = rep(1,4), lwd = 3, fig = c(0.05,1,0.05,1))
image.2(t_n, zlim = c(0,1), col = grad, annotate = "none", hct = t_cls, vct = t_cls); box()
par(fig = c(0.05,1,0,0.05), new = T)
image(matrix(seq_along(t_cls)), axes = F, col = name2color[ as.vector(sort(t_cls))]); box()
par(fig = c(0,0.05,0.05,1), new	= T)
image(t(seq_along(t_cls)), axes = F, col = name2color[ as.vector(sort(t_cls))]); box()
dev.off()

dc_cls = factor(color2name[ sin_cl@colors[ as.numeric(rownames(dc_n))]], levels = lin_ord)
png(paste0(supdir, "/FigS7c.png"), height = 1500, width = 1500)
par(mar = rep(1,4), lwd = 3, fig = c(0.05,1,0.05,1))
image.2(dc_n, zlim = c(0,1), col = grad, annotate = "none", hct = dc_cls,vct = dc_cls); box()
par(fig = c(0.05,1,0,0.05), new	= T)
image(matrix(seq_along(dc_cls)), axes = F, col = name2color[ as.vector(sort(dc_cls))]);box()
par(fig = c(0,0.05,0.05,1), new = T)
image(t(seq_along(dc_cls)), axes = F, col = name2color[ as.vector(sort(dc_cls))]); box()
dev.off()

png(paste0(supdir, "/FigS7b-c_colorbar.png"), height = 100, width = 1000)
par(mar = c(3,0,0,0))
image(matrix(1:1000), col = grad)
dev.off() 

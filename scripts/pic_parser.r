require(glmnet)

###############################################################################################################
#
# choose_lr_features
# ------------------
#
# Returns a set of genes to be used as features for inferring the mixing coefficient (alpha)
# Arguments:
#	mat_id:		a MetaCell ID for a count matrix
#	a_cells:	A set of cell identifiers for population A
#	b_cells:	A set of cell identifiers for population B
#	bad_genes:	a black list of genes not to be used as features
#	cor_n:		Number of genes correlated with cell size per population to be added to features list
#	must_haves:	A list of predefined genes to be included in the features list
#
###############################################################################################################

choose_lr_features = function(mat_id, a_cells, b_cells, bad_genes, cor_n = 100, must_haves = c()) {

	sc_mat = scdb_mat(mat_id)
	mcell_mat_ignore_cells(paste0(mat_id, "_a"), mat_id, unique(c(setdiff(colnames(sc_mat@mat), a_cells), sc_mat@ignore_cells)))
	mcell_mat_ignore_cells(paste0(mat_id, "_b"), mat_id, unique(c(setdiff(colnames(sc_mat@mat), b_cells), sc_mat@ignore_cells)))
	mcell_add_gene_stat(gstat_id=paste0(mat_id, "_a"), mat_id=paste0(mat_id, "_a"));
	mcell_add_gene_stat(gstat_id=paste0(mat_id, "_b"), mat_id=paste0(mat_id, "_b"));

	a_gstat = scdb_gstat(paste0(mat_id, "_a")); b_gstat = scdb_gstat(paste0(mat_id, "_b"))
	a_genes = with(a_gstat[setdiff(rownames(a_gstat)[ a_gstat$sz_cor_norm > 0], bad_genes),], tail(name[ order(sz_cor)],cor_n))
	b_genes = with(b_gstat[setdiff(rownames(b_gstat)[ b_gstat$sz_cor_norm > 0], bad_genes),], tail(name[ order(sz_cor)],cor_n))
	genes = unique(c(a_genes, b_genes, must_haves))
}

###############################################################################################################
#
# choose_mle_features
# ------------------
#
# Returns a set of genes to be used as features for inferring the metacell identities of the two contributing
# 	cell types (a_mc, b_mc)
# Arguments:
#       mc_id:		a MetaCell ID for a clustering model
#       mat_id:         a MetaCell ID for a count matrix
#       a_cells:        A set of cell identifiers for population A
#       b_cells:        A set of cell identifiers for population B
#       bad_genes:      a black list of genes not to be used as features
#       nms_per_clust:	Maximum number of differential genes to extract per metacell
#	nms_thresh:	A specificity (log2 fold change) threshold above which a gene is considered differential
#	shared_thresh:	A specificity (log2 fold change) threshold, such that genes differentially expressed
#			in BOTH populations will be excluded from the list
#       existing_list:  A list of predefined genes to be included in the features list (unless differentially
#			expressed in both populations
#
###############################################################################################################

choose_mle_features = function(mc_id, mat_id, a_cells, b_cells, bad_genes, nms_per_clust=20, nms_thresh=1, shared_thresh=1,
	existing_list = c()) {
	
	sc_cl = scdb_mc(mc_id); sc_mat = scdb_mat(mat_id)
	a_num = as.numeric(factor(sin_cl@mc[a_cells])); names(a_num) = a_cells
	scdb_add_mc(paste0(mc_id, "_a"), tgMCCov(a_num, setdiff(colnames(sc_mat@mat), names(a_num)), sc_mat))

	b_num = as.numeric(factor(sin_cl@mc[b_cells])); names(b_num) = b_cells
	scdb_add_mc(paste0(mc_id, "_b"), tgMCCov(b_num, setdiff(colnames(sc_mat@mat), names(b_num)), sc_mat))

	a_cl = scdb_mc(paste0(mc_id, "_a"))
	a_fp = log2(a_cl@mc_fp)
	b_cl = scdb_mc(paste0(mc_id, "_b"))
	b_fp = log2(b_cl@mc_fp)

	a_nms = choose_genes_from_clust(paste0(mc_id, "_a"), mat_id, nms_per_clust=nms_per_clust, nms_thresh=nms_thresh, ord = "max.col", bad_genes = bad_genes)
	b_nms = choose_genes_from_clust(paste0(mc_id, "_b"), mat_id, nms_per_clust=nms_per_clust, nms_thresh=nms_thresh, ord = "max.col", bad_genes = bad_genes)
	nms = union(existing_list, union(a_nms, b_nms))

	shared_markers = intersect(nms, intersect(rownames(a_fp), rownames(b_fp)))
	x = apply(abs(a_fp[shared_markers,]),1,max); y = apply(abs(b_fp[shared_markers, ]),1,max)
	bad_shared = names(which(x > shared_thresh & y > shared_thresh))
	setdiff(nms, bad_shared)
}

###############################################################################################################
#
# run_pic_seq
# ------------------
#
# General method to run the basic PIC-seq algorithm
# Returns a data frame with information per PIC.
# Arguments:
#       mc_id:          a MetaCell ID for a clustering model
#       mat_id:         a MetaCell ID for a count matrix
#	pic_umis:	a count matrix of the PIC population.
#       a_cells:        a set of cell identifiers for population A
#       b_cells:        a set of cell identifiers for population B
#	lr_features:	a set of genes to be used as features for inferring the mixing coefficient
#	mle_features:	a set of genes to be used as features for inferring the metacell identities 
#			of the two contributing cell types
#	fname:		A file name for output figure
#	numis:		Total number of UMI for simulated PIC
#	lr_k:		Number of synthetic PIC to simulate
#       bad_genes:      a black list of genes not to be used as features
#	comb:		a vector grouping cells from populations A and B. Simulated PIC will be
#			composed of cells derived from the same comb group.
#	reg:		A regularization coefficient used to regularize the multinomial distributions
#	downsample:	Whether to downsample the PIC count matrix to library sizes of "numis"
#
# Return values (for each PIC):
#	a_mc:		metacell identity of the contributing cell of population A
#	b_mc:		metacell identity of the contributing cell of population B
#	ll:		log likelihood score of the PIC-seq assignment
#	alpha:		Inferred mixing coefficient (alpha)
#
###############################################################################################################

run_pic_seq = function(mc_id, mat_id, pic_umis, a_cells, b_cells, 
	lr_features, mle_features, fname, numis = 1000, 
	lr_k = 2e4, bad_genes = c(), comb = NULL, reg = 1e-4, downsample=T) {
	
	sc_cl = scdb_mc(mc_id); sc_mat = scdb_mc(mat_id); 
	cells = union(a_cells, b_cells)
	lr_features = setdiff(lr_features, bad_genes)
	mle_features = setdiff(mle_features, bad_genes)
	if (is.null(comb)) {comb = rep("X", length(cells)); names(comb) = cells}

	message("Simulating ", lr_k, " PIC (numis = ", numis, ") for alpha estimation")
	lr_res = simulate_doublets(mat_id, a_cells, b_cells, lr_k, comb, numis = rep(numis, lr_k))
	sim_umis = lr_res$sim_umis; sim_info = lr_res$info
	sim_cells = names(which(colSums(sim_umis) == numis))
	sim_umis = sim_umis[,sim_cells]; sim_info = sim_info[sim_cells,]

        message("Building LR for alpha")
	alpha_fit = estimate_mixing(sim_umis, sim_info$alpha.1, lr_features, fname)

	message("Running MLE on PIC")
	
	if (downsample) {
		ds = .downsamp(pic_umis, numis)
	} else {
		ds = pic_umis
	}
	real_t_frac = predict(alpha_fit, newx = t(ds[lr_features,]), s = "lambda.min")
	alpha = real_t_frac[colnames(ds),]
	good_cells = intersect(colnames(ds), names(which(alpha >= 0 & alpha <= 1)))
	mle_res = assign_pics_to_singlets(mc_id, mat_id, ds[,good_cells], a_cells, b_cells, alpha[good_cells],
		verbose=T, bad_genes = bad_genes, markers = mle_features, reg = reg)

	t_mc = mle_res$a_mc; dc_mc = mle_res$b_mc
	names(t_mc) = rownames(mle_res); names(dc_mc) = rownames(mle_res)
	mle_res$alpha = alpha[good_cells]

        mle_res
}

###############################################################################################################
#
# confusion_matrix 
# ------------------
#
# Returns a matrix depicting for each pair of metacells the log ratio between the observed and expected
# numbers of shared neighbors between these metacells
# 
# Arguments:
#       mc_id:          a MetaCell ID for a clustering model
#       graph_id:       a MetaCell ID for a kNN graph structure
#       cells:          A set of cell identifiers over which to compute confusion
#
###############################################################################################################

confusion_matrix = function(mc_id, graph_id, cells) {
	sc_cl = scdb_mc(mc_id); sc_graph = scdb_cgraph(graph_id)
	edges = sc_graph@edges[ as.vector(sc_graph@edges$mc1) %in% cells & as.vector(sc_graph@edges$mc2) %in% cells,]
	edges$mc1 = as.vector(edges$mc1)
	edges$mc2 = as.vector(edges$mc2)
	edges$mc_in = sc_cl@mc[ edges$mc1]
	edges$mc_out = sc_cl@mc[ edges$mc2]
	confusion = table(edges$mc_in, edges$mc_out)
	exp_confu = outer(rowSums(confusion), colSums(confusion)) / sum(confusion)
	log2((confusion + 1) / (exp_confu + 1))
}

###############################################################################################################
#
# analyze_triplets
# ------------------
#
# General method to estimate the triplet frequency in a PIC population.
# Returns a list containing a dataframe with triplets information per PIC + two numbers depicting the estimated
# triplet frequencies (one for 2A+B triplets and one for A+2B triplets)
#
# Arguments:
#       mc_id:          a MetaCell ID for a clustering model
#       mat_id:         a MetaCell ID for a count matrix
#	graph_id:       a MetaCell ID for a kNN graph structure
#       pic_umis:       a count matrix of the PIC population.
#       a_cells:        a set of cell identifiers for population A
#       b_cells:        a set of cell identifiers for population B
#	mle_res:	A PIC-seq data frame (return value of run_pic_seq)
#       mle_features:   a set of genes to be used as features for inferring the metacell identities
#                       of the two contributing cell types
#       numis:          Total number of UMI for simulated PIC
#       bad_genes:      a black list of genes not to be used as features
#       outdir:         A directory name for output figures
#       prefix:         A prefix for output file names
#       reg:            A regularization coefficient used to regularize the multinomial distributions
#	confu_thresh:
#       downsample:     Whether to downsample the PIC count matrix to library sizes of "numis"
#	tr_k:		Number of synthetic doublets, 2A+B & A+2B triplets to simulate
#	tr_k2:		Maximum number of doublets, 2A+B & A+2B triplets to run PIC-seq on
#	tr_k3:		Total number of mixed synthetic doublets and triplets to use for testing triplet filtering
#
# Return values (as a list):
#       tr_p:           Estimated frequency of 2A+B triplets
#       rev_p:          Estimated frequency of A+2B triplets
#       mle_res:        An extended PIC-seq data frame with additional information per PIC:
#		forward_diff: 	Triplet score (likelihood ratio) when modeled as a 2A+B triplet
#		reverse_diff: 	Triplet score (likelihood ratio) when modeled as a A+2B triplet
#		forward_a1_mc:	Metacell identity of the first A cell contributing to a 2A+B triplet
#		forward_a2_mc:	Metacell identity of the second A cell contributing to a 2A+B triplet
#		forward_ll:	Log likelihood score of the best 2A+B assignment
#		forward_alpha:	The intra-mixing coefficient of the two A cells contributing to the 2A+B triplet
#		reverse_a1_mc:	Metacell identity of the first B cell contributing to a A+2B triplet
#		reverse_a2_mc:	Metacell identity of the second B cell contributing to a A+2B triplet
#		reverse_ll:	Log likelihood score of the best A+2B assignment
#		reverse_alpha:	The intra-mixing coefficient of	the two	B cells	contributing to	the A+2B triplet
#
###############################################################################################################

analyze_triplets = function(mc_id, mat_id, graph_id, pic_umis, a_cells, b_cells, mle_res,
	mle_features, numis=1000, bad_genes = c(), outdir = "./", prefix = "pic-seq",
	reg=1e-6, confu_thresh = 0, downsample=T,
	tr_k = 5e3, tr_k2 = 1e3, tr_k3 = 1e3) {

	
	t_confu = confusion_matrix(mc_id, graph_id, a_cells)
	dc_confu = confusion_matrix(mc_id, graph_id, b_cells)

        message("Estimate triplet frequency")

	sc_cl = scdb_mc(mc_id); sc_mat = scdb_mat(mat_id)
	if (downsample) {
	        ds = .downsamp(pic_umis, numis)
	} else {
		ds = pic_umis
	}
	good_cells = intersect(colnames(ds), rownames(mle_res))
        message("1. Calculate PIC triplet scores")
	
	mle_res = mle_res[good_cells,]
	alpha = mle_res$alpha; names(alpha) = rownames(mle_res)
	a_mc = mle_res$a_mc; names(a_mc) = rownames(mle_res)
	b_mc = mle_res$b_mc; names(b_mc) = rownames(mle_res)

	triplet_res = .test_for_triplets(mc_id, mat_id, graph_id, ds[,good_cells], a_cells, b_cells, alpha[good_cells], b_mc[good_cells],
		bad_genes = bad_genes, confu_thresh = confu_thresh, markers = mle_features, reg = reg)

	reverse_res = .test_for_triplets(mc_id, mat_id, graph_id, ds[,good_cells], b_cells, a_cells, 1 - alpha[good_cells], a_mc[good_cells],
        	bad_genes = bad_genes, confu_thresh = confu_thresh, markers = mle_features, reg = reg)
	
	message("2. Simulate triplets")

	db_res =  simulate_doublets(mat_id, a_cells, b_cells, tr_k, comb, numis = rep(numis, tr_k))
	tr_res =  simulate_doublets(mat_id, a_cells, b_cells, tr_k, comb, numis = rep(numis, tr_k), composition = c(1,1,2))
	rev_res = simulate_doublets(mat_id, a_cells, b_cells, tr_k, comb, numis = rep(numis, tr_k), composition = c(1,2,2))

	sim_db_umis = db_res$sim_umis; sim_db_info = db_res$info; sim_db_cells = names(which(colSums(sim_db_umis) == numis))
	sim_tr_umis = tr_res$sim_umis; sim_tr_info = tr_res$info;
	sim_rev_umis = rev_res$sim_umis; sim_rev_info = rev_res$info;
	
	colnames(sim_tr_umis) = paste0("Tr", seq_len(ncol(sim_tr_umis))); 
	rownames(sim_tr_info) = colnames(sim_tr_umis); 
	sim_tr_cells = names(which(colSums(sim_tr_umis) == numis))
	
	colnames(sim_rev_umis) = paste0("Rev", seq_len(ncol(sim_rev_umis))); 
	rownames(sim_rev_info) = colnames(sim_rev_umis); 
	sim_rev_cells = names(which(colSums(sim_rev_umis) == numis))
	
	sim_tr_info$mc.1 = sc_cl@mc[ as.vector(sim_tr_info$sim.1)]
	sim_tr_info$mc.2 = sc_cl@mc[ as.vector(sim_tr_info$sim.2)]
	sim_tr_info$enr1 = with(sim_tr_info, t_confu[cbind(as.character(mc.1), as.character(mc.2))])
	sim_tr_info$enr2 = with(sim_tr_info, t_confu[cbind(as.character(mc.2), as.character(mc.1))])
	sim_tr_info$good = with(sim_tr_info, enr1 < confu_thresh & enr2 < confu_thresh)
	sim_tr_cells = intersect(sim_tr_cells, rownames(sim_tr_info)[ sim_tr_info$good])
	
	sim_rev_info$mc.2 = sc_cl@mc[ as.vector(sim_rev_info$sim.2)]
	sim_rev_info$mc.3 = sc_cl@mc[ as.vector(sim_rev_info$sim.3)]
	sim_rev_info$enr1 = with(sim_rev_info, dc_confu[cbind(as.character(mc.2), as.character(mc.3))])
	sim_rev_info$enr2 = with(sim_rev_info, dc_confu[cbind(as.character(mc.3), as.character(mc.2))])
	sim_rev_info$good = with(sim_rev_info, enr1 < confu_thresh & enr2 < confu_thresh)
	sim_rev_cells = intersect(sim_rev_cells, rownames(sim_rev_info)[ sim_rev_info$good])

	tr_k2 = min(c(tr_k2, length(sim_db_cells), length(sim_tr_cells), length(sim_rev_cells)))
        message("3. Calculate simulation triplet score (k = ", tr_k2, ")")
	
	sim_db_alpha = sim_db_info$alpha.1; names(sim_db_alpha) = rownames(sim_db_info)
	sim_tr_alpha = sim_tr_info$alpha.1 + sim_tr_info$alpha.2; names(sim_tr_alpha) = rownames(sim_tr_info)
	sim_rev_alpha = sim_rev_info$alpha.1; names(sim_rev_alpha) = rownames(sim_rev_info)
	sim_cells = c(sample(sim_db_cells, tr_k2), sample(sim_tr_cells, tr_k2), sample(sim_rev_cells, tr_k2))
	sim_umis = cbind(sim_db_umis, sim_tr_umis, sim_rev_umis)[,sim_cells]
	sim_alpha = c(sim_db_alpha, sim_tr_alpha, sim_rev_alpha)[sim_cells]
	sim_mle_res = assign_pics_to_singlets(mc_id, mat_id, sim_umis, a_cells, b_cells, sim_alpha, verbose=T, bad_genes = bad_genes, markers = mle_features, reg = reg)
	
	sim_b = as.vector(sim_mle_res$b_mc); names(sim_b) = names(sim_alpha)
	sim_a = as.vector(sim_mle_res$a_mc); names(sim_a) = names(sim_alpha)

	sim_triplet_res = .test_for_triplets(mc_id, mat_id, graph_id, sim_umis, a_cells, b_cells, sim_alpha, sim_b,
        	bad_genes = bad_genes, confu_thresh = confu_thresh, markers = mle_features, reg = reg)
	sim_reverse_res = .test_for_triplets(mc_id, mat_id, graph_id, sim_umis, b_cells, a_cells, 1 - sim_alpha, sim_a,
	        bad_genes = bad_genes, confu_thresh = confu_thresh, markers = mle_features, reg = reg)

        message("4. Estimae triplet rate")
	#write.table(sim_triplet_res, sep = "\t", quote=F, col.names=NA, file = paste0(outdir, "/", prefix, "_sim_triplets.txt"))
	#write.table(sim_mle_res, sep = "\t", quote=F, col.names=NA, file = paste0(outdir, "/", prefix, "_sim_doublets.txt"))
	#write.table(sim_reverse_res, sep = "\t", quote=F, col.names=NA, file = paste0(outdir, "/", prefix, "_sim_reverse_triplets.txt"))

	sim_diff  = sim_triplet_res$ll - sim_mle_res$ll; names(sim_diff) = rownames(sim_triplet_res)
	real_diff = triplet_res$ll - mle_res$ll; names(real_diff) = rownames(triplet_res)

	real_reverse = reverse_res$ll - mle_res$ll;
	sim_reverse  = sim_reverse_res$ll - sim_mle_res$ll
	names(real_reverse) = rownames(mle_res); names(sim_reverse) = rownames(sim_mle_res)

	tr_p  = .estimate_triplets_rate(real_diff, sim_diff[intersect(names(sim_diff), sim_tr_cells)], sim_diff[intersect(names(sim_diff), sim_db_cells)], 
		prefix = paste0(prefix, "_forward"), outdir=outdir)

	rev_p = .estimate_triplets_rate(real_reverse, sim_reverse[intersect(names(sim_reverse), sim_rev_cells)], 
		sim_reverse[intersect(names(sim_reverse), sim_db_cells)], prefix = paste0(prefix, "_reverse"), outdir=outdir)
	mle_res$forward_diff = real_diff
	mle_res$reverse_diff = real_reverse
	colnames(triplet_res) = paste0("forward_", colnames(triplet_res))
	mle_res = cbind(mle_res, triplet_res)
        colnames(reverse_res) = paste0("reverse_", colnames(reverse_res))
	mle_res = cbind(mle_res, reverse_res)

	message("5. Simulate triplet filtration")
	.simulate_triplets_filtration(db_diff=sim_diff[sample(grep("Sim", names(sim_diff), v=T), round((1-tr_p) * tr_k3))], 
		tr_diff = sim_diff[sample(grep("Tr", names(sim_diff), v=T),round(tr_p * tr_k3))], outdir=outdir, prefix=paste0(prefix, "_forward"))

	.simulate_triplets_filtration(db_diff=sim_diff[sample(grep("Sim", names(sim_diff), v=T), round((1-rev_p) * tr_k3))], 
		tr_diff = sim_diff[sample(grep("Rev", names(sim_diff), v=T),round(rev_p * tr_k3))], outdir=outdir, prefix=paste0(prefix, "_reverse"))

	list(mle_res = mle_res, tr_p = tr_p, rev_p = rev_p)
}

###############################################################################################################
#
# generate_expected_pics_from_mle
# ------------------
#
# Returns an expected count matrix for PICs given their inferred contributing cells and mixing coefficients.
#
# Arguments:
#       mc_id:          a MetaCell ID for a clustering model
#	mc_tab:		A n*2 matrix (or n*3 when modeling triplets) containing the metacell assignments of the
#			contributing cells for each PIC.
#	alpha:		A vector (or n*2 matrix when modeling triplets) containing the mixing coefficients for
#			each PIC.
#	numis:		A vector depicting each PIC's cell size (or a fixed number for a downsampled matrix).
#       reg:            A regularization coefficient used to regularize the multinomial distributions.
#       bad_genes:      a black list of genes to be excluded from the expected count matrix.
###############################################################################################################

generate_expected_pics_from_mle = function(mc_id, mc_tab, alpha, numis, reg = 1e-6, bad_genes = c())  {

	sc_cl = scdb_mc(mc_id)
	mc_prof = sc_cl@e_gc[setdiff(rownames(sc_cl@e_gc), bad_genes),] + reg
	mc_prof = sweep(mc_prof, 2, colSums(mc_prof), "/")

	cells = rownames(mc_tab)
	if (is.null(ncol(alpha))) {
		alpha = cbind(alpha, 1-alpha)
	}
	if (ncol(alpha) < ncol(mc_tab)) {
		alpha = cbind(alpha, 1 - rowSums(alpha))
	}

	p_sum = matrix(0, nrow = nrow(mc_prof), ncol = length(cells), dimnames = list(rownames(mc_prof), cells))
	for (i in seq_len(ncol(mc_tab))) {
		p_sum = p_sum + sweep(mc_prof[, as.vector(mc_tab[,i])], 2, alpha[,i], "*")
	}
	exp_pic = sweep(p_sum, 2, numis, "*")
	exp_pic
}

###############################################################################################################
#
# simulate_doublets
# ------------------
#
# Simulates synthetic doublets or triplets composed of singlets from populations A and B.
# Returns a list containing a count matrix and a data frame.
#
# Arguments:
#       mat_id:         A MetaCell ID for a count matrix       
#       a_cells:        A set of cell identifiers for population A
#       b_cells:        A set of cell identifiers for population B
#       k:	       	Total number of synthetic PIC to simulate
#       comb:           a vector grouping cells from populations A and B. Simulated PIC will be
#                       composed of cells derived from the same comb group.
#       numis:          A vector depicting each synthetic PIC's cell size (use a fixed number for a downsampled matrix).
#			If value is Inf - each PIC's cell size will be determined by its contributing cells.
#       composition:    A vector of 1's and 2's denoting the PIC composition. (1;2) is a doublet,
#			(1;1;2) is a 2A+B triplet and (1;2;2) is a A+2B triplet.
#	sample:		If ="sc", each synthetic PIC is created by summing the gene expression vectors of its contributors.
#			If ="mc", each synthetic PIC is created by sampling transcripts from the multinomial distributions
#			of its contributors' metacell identities.
#       mc_id:          a MetaCell ID for a clustering model (used only ifs sample="mc")
#       reg:            A regularization coefficient used to regularize the multinomial distributions (used only if sample="mc")
#
#	Return values (as a list):
#	sim_umis:	A count matrix of the synthetic PICs
#	sim_info:	A data frame containing for each synthetic PIC:
#		sim.[1-n]:	the i'th contributing cell
#		alpha.[1-n]:	the i'th mixing coefficient (fraction of transcripts derived from the i'th contributor)
###############################################################################################################

simulate_doublets = function(mat_id, a_cells, b_cells, k, comb = NULL, numis = Inf, composition = c(1,2), sample = "sc",
	mc_id = NULL, reg = 1e-6) {

	cells = union(a_cells, b_cells)
	sc_mat = scdb_mat(mat_id)
	umis = as.matrix(sc_mat@mat[,cells])
	if (is.null(comb)) { 
		comb = rep(1, length(cells)); names(comb) = cells
	}
	comb_dist = table(comb[cells])
	k_dist = round(comb_dist / sum(comb_dist) * k)
	k = sum(k_dist); numis = numis[seq_len(k)]

	sim_a = NULL; sim_b = NULL
	sim = matrix(NA, nrow = length(composition), ncol = k, dimnames = list(composition, rep(names(k_dist), k_dist)))
	for (c in names(k_dist)) {
		pos = which(colnames(sim) == c)
		sim[ composition == 1, pos] = sample(a_cells[ comb[a_cells] == c], length(pos) * sum(composition == 1), replace = T)
		sim[ composition == 2, pos] = sample(b_cells[ comb[b_cells] == c], length(pos) * sum(composition == 2), replace = T)
	}
	umicount = colSums(umis)
	uc = matrix(umicount[sim], nrow = nrow(sim))
	numis = pmin(colSums(uc), numis)
	nuc = round(sweep(uc, 2, numis / colSums(uc), "*"))
	nuc = pmin(nuc, uc); numis = colSums(nuc)
	if (sample == "sc") {
		sim_umis = matrix(0, nrow = nrow(umis), ncol = k, dimnames = list(rownames(umis), paste0("Sim", seq_len(k))))
		for (i in seq_len(nrow(sim))) {
			sim_umis = sim_umis + .downsamp_var(umis[,sim[i,]], nuc[i,])
		}
	} else if (sample == "mc") {
        	sc_cl = scdb_mc(mc_id); #sc_mat = scdb_mat(mat_id)
	        mc_prof = sc_cl@e_gc[setdiff(rownames(sc_cl@e_gc), bad_genes),]
	        mc_prof = sweep(mc_prof, 2, colSums(mc_prof), "/")
		sim_umis = matrix(0, nrow = nrow(mc_prof), ncol = k, dimnames = list(rownames(mc_prof), paste0("Sim", seq_len(k))))
		for (i in seq_len(nrow(sim))) {
			p = mc_prof[, sc_cl@mc[sim[i,]]] + reg
                        sim_umis = sim_umis + apply(rbind(nuc[i,], p), 2, function(x) .sample_doublets_from_prof(x[-1], x[1]))
		}
	}
	colnames(sim) = colnames(sim_umis); rownames(sim) = seq_len(nrow(sim))
	list(sim_umis = sim_umis, info = data.frame(sim = t(sim), alpha = t(nuc) / numis))
}

###############################################################################################################
#
# estimate_mixing
# ------------------
#
# Constructs a linear regression (LR) model to infer the mixing coefficient (alpha) from a set of synthetic PICs
#
# Arguments:
#       sim_umis:       A count matrix of synthetic PICs.
#       alpha:        	A vector denoting the mixing coefficient of each synthetic PIC.
#       genes:		A set of genes used as features for the LR model.
#       fname:		A file name for output figure.
#       normalize:	Whether to size-normalize the synthetic UMIs prior to building the model.
#
# Return value:
#	a glmnet object containing the linear regression model.
#
###############################################################################################################

estimate_mixing = function(sim_umis, alpha, genes, fname, normalize = T) {
	
	sim_n = sweep(sim_umis, 2, colSums(sim_umis), "/") * 1000
	if (normalize) { us = sim_n} else {us = sim_umis}
	alpha_fit = cv.glmnet(t(us[genes,]), alpha)
	alpha_tag = predict(alpha_fit, newx = t(us[genes,]), s = "lambda.min")
	i = which(alpha_fit$lambda == alpha_fit$lambda.min)
	if (!is.null(fname)) {
		png(fname, height=1000, width=1000)
		#par(bg = bg, fg = fg, col.axis = col.axis, col.lab = col.lab)
	}
	lim = c(0,1)
	plot(alpha_tag, alpha, pch = 20, , col = "gray", axes = F, xlab = "", ylab = "", cex = 2, xlim = lim, ylim = lim,
		main = round(1 - alpha_fit$cvm[i] / var(alpha),4))
	axis(1); axis(2); abline(coef = c(0,1), lwd = 2)
        if (!is.null(fname)) {dev.off()}
	alpha_fit
}

###############################################################################################################
#
# assign_pics_to_singlets
# ------------------
#
# Find the maximum likelihood assignment of a PIC to two contributing metacells, given alpha
#
# Arguments:
#       mc_id:          a MetaCell ID for a clustering model
#       mat_id:         a MetaCell ID for a count matrix
#       pic_umis:       a count matrix of the PIC population.
#       a_cells:        a set of cell identifiers for population A
#       b_cells:        a set of cell identifiers for population B
#       alpha:          A vector denoting the inferred mixing coefficient of each PIC.
#	bins:		Used to stratitfy PICs with similar alpha values together by binning alpha, to reduce computing time.
#			If =0 each PIC likelihoods are computed separately.
#       reg:            A regularization coefficient used to regularize the multinomial distributions
#	verbose:	If TRUE, a counter of progress is displayed.
#	likelihoods_fn:	A text filename to store likelihood values.
#       bad_genes:      Genes to exclude and not use a features
#       markers:      	A set of genes to use a features. If =NULL, likelihoods are computed over all genes not in bad_genes.
#
# Return value (for each PIC):
#	a_mc:		Metacell assignment of the contributing cells from population A.
#	b_mc:		Metacell assignment of the contributing cells from population B.
#	ll:		Maximum likelihood value.
#
###############################################################################################################

assign_pics_to_singlets = function(mc_id, mat_id, pic_umis, a_cells, b_cells, alpha, bins = 0,
	reg = 1e-6, verbose = T, likelihoods_fn = NULL, bad_genes = bad_genes, markers = NULL) {

	sc_cl = scdb_mc(mc_id); sc_mat = scdb_mat(mat_id)
	umis = as.matrix(sc_mat@mat[,names(sc_cl@mc)])
        if (is.null(markers)) {
                mc_prof = sc_cl@e_gc[setdiff(rownames(sc_cl@e_gc), bad_genes),]
	} else {
                mc_prof = mc_compute_unnorm_fp(sc_cl@mc, umis[setdiff(markers, bad_genes),])
        }
	mc_prof = sweep(mc_prof, 2, colSums(mc_prof), "/")

	sinsqr = melt(outer(table(sc_cl@mc[a_cells]), table(sc_cl@mc[b_cells])))
	sinsqr$Var1 = as.vector(sinsqr$Var1); sinsqr$Var2 = as.vector(sinsqr$Var2)
	sinsqr$id = paste0(sinsqr$Var1, "-", sinsqr$Var2)
	rownames(sinsqr) = sinsqr$id

	p_a = mc_prof[,as.character(unique(sinsqr$Var1))] + reg
	p_b = mc_prof[,as.character(unique(sinsqr$Var2))] + reg

	db_likelihoods = matrix(NA, nrow = ncol(pic_umis), ncol = nrow(sinsqr), dimnames = list(colnames(pic_umis), rownames(sinsqr)))
	n = nrow(db_likelihoods); processed = 0; range = seq(0,1,0.05)
	if (bins > 0) {
		h = hist(alpha, bins, plot = F)
		alpha_cut = h$mids[as.numeric(cut(alpha, h$breaks, include.lowest = T))]
		names(alpha_cut) = names(alpha)
		for (a in unique(alpha_cut)) {
			cells = names(which(alpha_cut == a))
	                u = pic_umis[rownames(mc_prof), cells];
	                p = with(sinsqr, p_a[, as.character(Var1)] * a +
        	                p_b[, as.character(Var2)] * (1 - a))
	                db_likelihoods[cells,] = t(u) %*% log(p)
			new_processed = processed + (length(cells) / n)
			if (sum(range > processed & range <= new_processed) > 0) {message(round(new_processed * 100), "%")}
			processed = new_processed			
		}
	} else {
		for (cell in colnames(pic_umis)) {
			u = pic_umis[rownames(mc_prof), cell]; a = alpha[cell]
			p = with(sinsqr, p_a[, as.character(Var1)] * a + 
				p_b[, as.character(Var2)] * (1 - a))
		        db_likelihoods[cell,] = t(u) %*% log(p)
			new_processed = processed + 1
	                if (sum(range > processed & range <= new_processed) > 0) {message(round(new_processed * 100), "%")}
			processed = new_processed
		}		
	}
	lmulti = apply(pic_umis[rownames(mc_prof), ], 2, stirlings_lmulti)
	cor_likelihoods = db_likelihoods + lmulti
	if (!is.null(likelihoods_fn)) {
		write.table(cor_likelihoods, sep = "\t", quote = F, col.names = NA, file = likelihoods_fn)
	}
	cell_likelihood = apply(cor_likelihoods[ rowSums(!is.na(db_likelihoods)) > 0,],1,max, na.rm = T)
	mle = apply(db_likelihoods[ rowSums(!is.na(db_likelihoods)) > 0,],1,which.max)

	a_mc = sinsqr[mle, "Var1"]; b_mc = sinsqr[mle, "Var2"]
	names(a_mc) = names(mle); names(b_mc) = names(mle)
	data.frame(a_mc = a_mc, b_mc = b_mc, ll = cell_likelihood)	
}

.test_for_triplets = function(mc_id, mat_id, graph_id, pic_umis, a_cells, b_cells, a_alpha, b_mc, #homo_alpha = 0.5, #bins = 200,
        reg = 1e-6, verbose = T, likelihoods_fn = NULL, bad_genes = c(), confu_thresh = Inf, exclusion_matrix = NULL, markers = NULL) {

	sc_cl = scdb_mc(mc_id); sc_mat = scdb_mat(mat_id)
	umis = as.matrix(sc_mat@mat[,names(sc_cl@mc)])
	if (is.null(markers)) {
	        mc_prof = sc_cl@e_gc[setdiff(rownames(sc_cl@e_gc), bad_genes),]	
	} else {
		mc_prof = mc_compute_unnorm_fp(sc_cl@mc, umis[setdiff(markers, bad_genes),])
	}
        mc_prof = sweep(mc_prof, 2, colSums(mc_prof), "/")

	confu_enr = confusion_matrix(mc_id, graph_id, a_cells)

	# create all homo_alpha combinations of mc from a
	umicount = colSums(umis)
	uc_med = tapply(umicount[a_cells], sc_cl@mc[a_cells], median)
        sinsqr = melt(outer(table(sc_cl@mc[a_cells]), table(sc_cl@mc[a_cells])))
	sinsqr = sinsqr[ sinsqr$Var2 >= sinsqr$Var1,]
        sinsqr$Var1 = as.vector(sinsqr$Var1); sinsqr$Var2 = as.vector(sinsqr$Var2)
	sinsqr$enr1 = with(sinsqr, confu_enr[cbind(as.character(Var1), as.character(Var2))])
	sinsqr$enr2 = with(sinsqr, confu_enr[cbind(as.character(Var2), as.character(Var1))])
	sinsqr = sinsqr[ (sinsqr$enr1 < confu_thresh & sinsqr$enr2 < confu_thresh) | sinsqr$Var2 == sinsqr$Var1,]
	sinsqr$uc1 = uc_med[as.character(sinsqr$Var1)]
	sinsqr$uc2 = uc_med[as.character(sinsqr$Var2)]
	sinsqr$homo_alpha = with(sinsqr, uc1 / (uc1+uc2))
	if (!is.null(exclusion_matrix)) {
		sinsqr$good = with(sinsqr, exclusion_matrix[ cbind(as.character(Var1), as.character(Var2))])
		sinsqr = sinsqr[ sinsqr$good,]
	}
        sinsqr$id = paste0(sinsqr$Var1, "-", sinsqr$Var2)
        rownames(sinsqr) = sinsqr$id

        p_a = mc_prof[,names(table(sc_cl@mc[a_cells]))] + reg
	p_a = with(sinsqr, sweep(p_a[, as.character(Var1)],2,homo_alpha,"*") + sweep(p_a[, as.character(Var2)],2,1 - homo_alpha,"*"))
	colnames(p_a) = sinsqr$id
	
	p_b = mc_prof + reg
        db_likelihoods = matrix(NA, nrow = ncol(pic_umis), ncol = nrow(sinsqr), dimnames = list(colnames(pic_umis), rownames(sinsqr)))

	for (cell in colnames(pic_umis)) {
                u = pic_umis[rownames(mc_prof), cell]
		a = a_alpha[cell]; b = as.character(b_mc[cell])
		#if (!(b %in% colnames(p_b))) {
		#message(cell, ": ", b)}
		p = p_a * a + p_b[,b] * (1 - a)
                db_likelihoods[cell, colnames(p)] = t(u) %*% log(p)
	}
        lmulti = apply(pic_umis[rownames(mc_prof), ], 2, stirlings_lmulti)
        cor_likelihoods = db_likelihoods + lmulti
        cell_likelihood = apply(cor_likelihoods[ rowSums(!is.na(db_likelihoods)) > 0,],1,max, na.rm = T)
	mle = apply(db_likelihoods[ rowSums(!is.na(db_likelihoods)) > 0,],1,which.max)
        a1_mc = sinsqr[mle, "Var1"]; a2_mc = sinsqr[mle, "Var2"]
        names(a1_mc) = names(mle); names(a2_mc) = names(mle)
	res = data.frame(a1_mc = a1_mc, a2_mc = a2_mc, ll = cell_likelihood)
	res$alpha = sinsqr[ paste0(a1_mc, "-", a2_mc), "homo_alpha"]
	res
}

.sample_doublets_from_prof = function(prof, nu) {
	cdf = c(0, cumsum(prof)); names(cdf) = c(names(prof), "XXX")
	samps = runif(nu)
	X = as.numeric(table(factor(names(cdf[ as.numeric(cut(samps, cdf, include.lowest=T))]), levels = names(prof))))
	names(X) = names(prof)
	X
}

.simulate_triplets_filtration = function(db_diff, tr_diff, outdir="./", prefix="pic-seq") {
	df 		= data.frame(source = c(rep("Doublet", length(db_diff)), rep("Triplet", length(tr_diff))), diff = c(db_diff, tr_diff))
	ord_df 		= df[ order(-df$diff, sample(nrow(df))),]
	ord_df$tp 	= cumsum(ord_df$source == "Triplet"); 
	ord_df$fp 	= cumsum(ord_df$source == "Doublet")
	ord_df$fn 	= length(tr_diff) - ord_df$tp
	ord_df$retained = with(ord_df, nrow(df) - tp-fp)
	ord_df$fnr 	= with(ord_df, fn/retained)
	ord_df$discarded_r = with(ord_df, 1 - retained / nrow(df))
	png(paste0(outdir, "/", prefix, "_filtration_rate.png"), height=1000, width=1000)
	with(ord_df, plot(discarded_r, fnr, type = "l", lwd = 4, col = "red", axes=F, xlab="", ylab=""))
	grid(col = "black"); axis(1); axis(2)
	abline(v=length(tr_diff)/nrow(df), lwd=2)
	dev.off()
}

.tr_posterior = function(x, tr_prior = tr_prior, p_tr, p_db) {
	x = as.character(x); p_tr[x] * tr_prior / (p_tr[x] * tr_prior + p_db[x] * (1 - tr_prior))
}

.estimate_triplets_rate = function(real_diff, tr_diff, db_diff, outdir = ".", prefix = "triplet",
	cols = c("green3", "navyblue", "gray"), reg=0.1) {

	png(paste0(outdir, "/", prefix, "_ecdf.png"), height=1000,   width = 1000)
	all_diff = c(real_diff, tr_diff, db_diff)
	plot(ecdf(-all_diff), col = "white"); grid(col = "black")
	lines(ecdf(-real_diff), cex = 0, col = cols[1], lwd = 6)
	lines(ecdf(-tr_diff), cex = 0, col = cols[2], lwd = 6)
	lines(ecdf(-db_diff), cex = 0, col = cols[3], lwd = 6)	
	dev.off()

	all_diff = pmax(all_diff, 0)
	all_round = round(all_diff)
	all_ass = c(rep("Real", length(real_diff)), rep("Triplet", length(tr_diff)), rep("Doublet", length(db_diff)))
	#write.table(data.frame(names(all_round), all_ass, all_round), sep = "\t", quote = F)
	bins = seq(0, max(all_round))

	db_dist = table(factor(all_round[all_ass == "Doublet"], levels = bins))
	tr_dist = table(factor(all_round[all_ass == "Triplet"], levels = bins))
	p_db = (reg + db_dist) / sum(reg + db_dist)
	p_tr = (reg + tr_dist) / sum(reg + tr_dist)
	vals = all_round[ all_ass == "Real"]

	priors = seq(0,1,0.01)
	posteriors = sapply(priors, function(x) sum(.tr_posterior(vals, x, p_tr,  p_db))) / length(vals)
	diff = abs(posteriors - priors); names(diff) = priors
	prior = as.numeric(names(sort(diff)[3]))
	message(prefix, " triplet rate: ", prior)
	prior
}

.downsamp_var = function (umis, n, replace = F) {
        m = nrow(umis)
        .downsamp_one=function(v,n, replace = F){
                a = tabulate(sample(rep(1:length(v),times=v),replace=replace,size=n),nbins=m)
                return (a)
        }
	ret = umis[,colSums(umis) >= n]; n = n[colSums(umis) >= n]
	idx = which(colSums(ret) > n)
	if (length(idx) > 1) {
		ret[,idx] = apply(rbind(n, ret)[, idx], 2, function(x) .downsamp_one(x[-1], x[1]))
	} else if (length(idx) == 1) {
		ret[,idx] = .downsamp_one(ret[,idx], n[idx])
	}
	rownames(ret) = rownames(umis)
	return(ret)
}

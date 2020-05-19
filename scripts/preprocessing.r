set.seed(1111)
library(metacell)
source("scripts/metacell_functions.r")
source("scripts/pic_parser.r")

dir.create("saved_work")
dir.create("saved_work/metacell")
scdb_init("saved_work/metacell", force_reinit=F)

##################################
# Import the entire count matrix

mcell_import_multi_mars("all", "annotations/metadata.txt", "output/umi.tab/", force = F)

mat = scdb_mat("all")
nms = c(rownames(mat@mat), rownames(mat@ignore_gmat))
bad_genes = unique(c(grep("ERCC", nms, v=T), grep("^mt-", nms, v=T), "Gm42418", "Gm26917", "Gm29216"))

#########################################################
#							#
# Import the MetaCell structures from the paper		#
#							#
#########################################################

import_metacell_structure("in_vitro", "import/in_vitro/", "all", bad_genes)
import_metacell_structure("ln", "import/ln/", "all", bad_genes)

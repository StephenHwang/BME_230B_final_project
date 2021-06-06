# Running cibersort

home.dir <- '/home/stephen/Documents/classes/bme/230B/decon/'
setwd(home.dir)

# .libPaths(c("/soe/vfriedl/R/x86_64-redhat-linux-gnu-library/3.5"
            # ,"/projects/sysbio/apps/x86_64/Rlib"
            # ,"/usr/lib64/R/library"
            # ,"/usr/share/R/library"))
source('bin/CIBERSORT.R')

# run Cibersort w/
# (1) Signature file (for example LM22.txt)
# (2) Mixture file (the expression matrix for example ExampleMixtures-GEPs.txt)
# (3) 100 permutations (to get a p-value)
# (4) quantile normalization = FALSE (recommended for RNA-Seq data)
# rest default, which basically means to not run 'absolute' mode
# (5) absolute (default=FALSE): run Cibersort in absolute mode (not fractions)
# (6) abs_method: if set absolute=TRUE, choose 'no.sumto1' or 'sig.score'. 'no.sumto1' removes constraints

#run Cibersort, remember the first file should be the signature matrix, the second file is mixture matrix
# example
# sig.matrix <-"PBMC_signature_matrix/PBMC_example_signature_matrix_v1.tsv"
# mixture.matrix <- "train_bulk_mixture_data/CellLineA_Tumor_bulk_composition_noise_0.tsv"

sig.matrix <-"signature_matrix/LM5.txt"
mixture.matrix <- "train_bulk_mixture_data/CellLineA_Tumor_bulk_composition_noise_0.tsv"

results <- CIBERSORT(sig.matrix, mixture.matrix, perm=100, QN=FALSE, absolute=TRUE, abs_method='no.sumto1')
write.table(results,file="./Cibersort_ExampleOutput.tsv",sep = "\t", quote = F,row.names = T, col.names = T) # write result table to file
# 1:21 -

# View results
results


test.mixture.matrix.7 <- 'test_bulk_mixture_data/bulk_composition_from_scRNA-seq/bulk_composition_7.tsv'
test.mixture.matrix.8 <- 'test_bulk_mixture_data/bulk_composition_from_scRNA-seq/bulk_composition_8.tsv'

tmp <- results[, c("B cells naive", 'B cells memory')]
tmp <- results[, c("T cells CD4 naive",  "T cells CD4 memory resting")]

colSums(tmp)

################################################################################
# Deconvolution results evaluation, there are two metrics we use for evaluation:
#  pearson correlation
#  RMSE

#function to annotate RMSE and correlation score in plots
Corner_text <- function(text, location="topright", text.col="black"){
  legend(location,legend=text, bty ="n", pch=NA, text.col=text.col)
}

#function of RMSE
RMSE <- function(m, o){
  sqrt(mean((m - o)^2))}
file <- "Cibersort_ExampleOutput.tsv" #load cibersort results
cells <- c("T","B","M", "NK", "Dendritic")

#pdf("Cibersort_ExampleOutput_plot.pdf", width=6, height=12)
#par(mfrow=c(6,3), mar=c(3,3,3,3))  #arrange the plot in pdf and set the margin
deconv.results <- read.delim(file)

#load cell proportions from mixture file
true.T.prop <- unlist(lapply(rownames(deconv.results), function(x) as.numeric(unlist(strsplit(x, split="::", fixed=T))[12])/100))
true.B.prop <- unlist(lapply(rownames(deconv.results), function(x) as.numeric(unlist(strsplit(x, split="::", fixed=T))[2])/100))
true.M.prop <- unlist(lapply(rownames(deconv.results), function(x) as.numeric(unlist(strsplit(x, split="::", fixed=T))[6])/100))
true.NK.prop <- unlist(lapply(rownames(deconv.results), function(x) as.numeric(unlist(strsplit(x, split="::", fixed=T))[8])/100))
true.Dendr.prop <- unlist(lapply(rownames(deconv.results), function(x) as.numeric(unlist(strsplit(x, split="::", fixed=T))[4])/100))

## Plot results
#load results into a matrix, change colnames according
pbmc.results <- do.call(cbind, list(true.T.prop, true.B.prop, true.M.prop, true.NK.prop, true.Dendr.prop,
                                    as.numeric(deconv.results[,"avg_T.cell"]),
                                    as.numeric(deconv.results[,"avg_B.cell"]),
                                    as.numeric(deconv.results[,"avg_Monocyte"]),
                                    as.numeric(deconv.results[,"avg_NK.cell"]),
                                    as.numeric(deconv.results[,"avg_Dendritic.cell"])))
colnames(pbmc.results) <- c("true.T.prop", "true.B.prop", "true.M.prop", "true.NK.prop", "true.Dendritic.prop", "T.result", "B.result", "M.result", "NK.result", "Dendritic.result")
for (cell in cells){
  plot(x=pbmc.results[,paste(c("true.",cell,".prop"), collapse="")], y=pbmc.results[,paste(c(cell,".result"), collapse="")],
       xlab=paste(c("true.",cell,".prop"), collapse=""), ylab=paste(c("cibersort", cell, "result"), collapse="."),
       main=paste(c(cell, "cell", "in", "Example"), collapse=" "), cex=0.5)
  mtext(paste(c("true.",cell,".prop"), collapse=""), side=1, line=2, font=1, cex=0.8)
  mtext(paste(c("cibersort", cell, "result"), collapse="."), side=2, line=2, font=1, cex=0.8)
  Corner_text(paste("corr=", round(cor.test(pbmc.results[,paste(c("true.",cell,".prop"), collapse="")],  #calculate pearson correlation score
                                            pbmc.results[,paste(c(cell,".result"), collapse="")], method=c("pearson"))$estimate, 5), sep=" "),location="topleft")
  Corner_text(paste("RMSE=", round(mean(na.omit(unlist(lapply(1:200, function(i) RMSE(pbmc.results[,paste(c(cell,".result"), collapse="")][i], pbmc.results[,paste(c("true.",cell,".prop"), collapse="")][i]))))),3), sep=" "),location="bottomright")
}
#dev.off()



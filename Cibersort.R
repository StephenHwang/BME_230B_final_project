# Running cibersort

home.dir <- '/home/stephen/Documents/classes/bme/230B/decon/'
setwd(home.dir)

source('bin/CIBERSORT.R')

# run Cibersort w/
# (1) Signature file (for example LM22.txt)
# (2) Mixture file (the expression matrix for example ExampleMixtures-GEPs.txt)
# (3) 100 permutations (to get a p-value)
# (4) quantile normalization = FALSE (recommended for RNA-Seq data)
# rest default, which basically means to not run 'absolute' mode
# (5) absolute (default=FALSE): run Cibersort in absolute mode (not fractions)
# (6) abs_method: if set absolute=TRUE, choose 'no.sumto1' or 'sig.score'. 'no.sumto1' removes constraints

# example datasets:
# sig.matrix <-"signature_matrix/PBMC_example_signature_matrix_v1.tsv"
# mixture.matrix <- "train_bulk_mixture_data/CellLineA_Tumor_bulk_composition_noise_0.tsv"

# Cell types:
#   B cells
#   T cells
#   Monocytes
#   Nature killer cells
#   Dendritic cells


#sig.matrix <-"signature_matrix/LM5.txt"        # 550 genes
#sig.matrix <-"signature_matrix/holiday.txt"    # too many genes, takes too long to run (12366)
#sig.matrix <-"signature_matrix/LM_joint.txt"    # not as good at LM5

sig.matrix <-"signature_matrix/signature_1.txt"
View(read.table(sig.matrix, header=T, sep="\t", row.names=1, check.names=F))

mixture.matrix <- "train_bulk_mixture_data/CellLineA_Tumor_bulk_composition_noise_0.tsv"

results <- CIBERSORT(sig.matrix, mixture.matrix, perm=100, QN=FALSE, absolute=TRUE, abs_method='no.sumto1')
write.table(results,file="./Cibersort_ExampleOutput.tsv",sep = "\t", quote = F,row.names = T, col.names = T) # write result table to file
# 1:21 -

## LM22_joint: join columns and re-name
# merged_results <- cbind(
# rowSums(results[,  c("B cells naive", "B cells memory" )]),
# rowSums(results[, c("T cells CD8",
#             "T cells CD4 naive",
#             "T cells CD4 memory resting",
#             "T cells CD4 memory activated",
#             "T cells follicular helper",
#             "T cells regulatory (Tregs)",
#             "T cells gamma delta")]),
# rowSums(results[, c("NK cells resting", "NK cells activated")]),
# results[, c("Monocytes")],
# rowSums(results[, c("Dendritic cells resting", "Dendritic cells activated")])
# )
# colnames(merged_results) <- c("avg_B.cell", "avg_T.cell", "avg_NK.cell", "avg_Monocyte", "avg_Dendritic.cell")
# deconv.results <- merged_results

# View results
results






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



# Overall results
matrix1 <- pbmc.results[, 1:5]
matrix2 <- pbmc.results[, 6:10]
results.cor <- cor(c(matrix1), c(matrix2), method='pearson')
results.RMSE <- RMSE(c(matrix1), c(matrix2))

cat(paste(paste('Sig:', gsub('signature_matrix/', '', sig.matrix)),
      paste('cor:', results.cor),
      paste('RMSE:', results.RMSE),
      sep='\n'))








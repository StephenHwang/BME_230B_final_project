# Running cibersort on test data
home.dir <- '/home/stephen/Documents/classes/bme/230B/decon/'
setwd(home.dir)

source('bin/CIBERSORT.R')

# sig.matrix <-"signature_matrix/signature_5.txt"
sig.matrix <-"signature_matrix/signature_2.txt"

# make a list of files to run thru
# gsub to correct output name
results.path <- '/home/stephen/Documents/classes/bme/230B/decon/test_bulk_mixture_data/results/'
data.dir <- '/home/stephen/Documents/classes/bme/230B/decon/test_bulk_mixture_data/bulk_composition_from_cell_lines/'
datasets.paths <- list.files(data.dir)

dataset <- datasets.paths[1]
# "Running: CellLineE_Tumor_bulk_composition_noise_0.3.tsv"

# for (dataset in datasets.paths[12:33]) {
for (dataset in datasets.paths) {
  print(paste0('Running: ', dataset))

  dataset.input.path <- paste0(data.dir, dataset)
  dataset.output.path <- paste0(results.path, gsub('.tsv', '', dataset), '_Output.tsv')
  print(dataset.input.path)
  print(dataset.output.path)

  results <- CIBERSORT(sig.matrix, dataset.input.path, perm=100, QN=FALSE, absolute=TRUE, abs_method='no.sumto1')

  write.table(results, file=dataset.output.path, sep = "\t", quote = F,row.names = T, col.names = T) # write result table to file

  rm(dataset.path)
  rm(dataset.name)
  rm(results)
}

print('done')

library(scHCL)
data(hcl_lung)
data(ref.expr)

write_csv(ref.expr, gzfile('/tmp/scHCL_ref.expr.csv.gz'))
write_csv(hcl_lung, gzfile('/tmp/scHCL_lung.expr.csv.gz'))


rownames(hcl_lung) = toupper(rownames(hcl_lung))

hcl_result <- scHCL(scdata = hcl_lung, numbers_plot = 3)

cor = hcl_result$cors_matrix

# what the correlation of one specific cell with all the reference
cellid = "E18_2_C06"
cor[,cellid]

# lets see if we can recreate that
expr_vector = hcl_lung[cellid]

hcl_result$cors_matrix[,cellid]

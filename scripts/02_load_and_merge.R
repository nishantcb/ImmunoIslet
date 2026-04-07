library(Seurat)
library(dplyr)

data_dir <- ""

files <- c(
  paste0(data_dir, "GSM2230757_human1_umifm_counts.csv.gz"),
  paste0(data_dir, "GSM2230758_human2_umifm_counts.csv.gz"),
  paste0(data_dir, "GSM2230759_human3_umifm_counts.csv.gz"),
  paste0(data_dir, "GSM2230760_human4_umifm_counts.csv.gz")
)

library(Seurat)
library(dplyr)
#############################################################################
# For all 4 Files
##############################################################################
data_list <- list()

for (i in seq_along(files)) {
  
  df <- read.csv(gzfile(files[i]), header = TRUE, check.names = FALSE)
  df <- df[!(df[,1] %in% c("barcode", "assigned-cluster")), ]
  # Set gene names
  rownames(df) <- df[,1]
  df <- df[,-1]
  
  # Remove duplicates
  rownames(df) <- make.unique(rownames(df))
  
  df[df == ""] <- 0
  df[df == "NA"] <- 0
  
  mat <- as.matrix(df)
  mat <- apply(mat, 2, function(x) as.numeric(as.character(x)))
  mat <- as.matrix(mat)
  
  mat[is.na(mat)] <- 0
  
  mat <- t(mat)
  obj <- CreateSeuratObject(counts = mat, project = paste0("Sample_", i))
  
  data_list[[i]] <- obj
}

# Merge
islet <- merge(
  x = data_list[[1]],
  y = data_list[2:length(data_list)],
  add.cell.ids = paste0("D", seq_along(data_list))
)

print(islet)

# Save
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
saveRDS(islet, "data/processed/islet_merged.rds")




#####################################################################

library(Seurat)
library(dplyr)

library(Seurat)

file <- "GSM2230757_human1_umifm_counts.csv.gz"

df <- read.csv(gzfile(file), header = TRUE, check.names = FALSE)

# 🔍 Inspect
print(head(df[,1:5]))

# Remove extra column
df <- df[, !(colnames(df) %in% c("Unnamed: 0"))]

df_clean <- df[, !(colnames(df) %in% c(
  "Unnamed: 0",
  "barcode",
  "assigned-cluster",
  "assigned_cluster",
  ""
))]

# ✅ Ensure unique gene names
colnames(df_clean) <- make.unique(colnames(df_clean))


# ✅ Extract cell names (make unique)
cell_names <- make.unique(as.character(df$barcode))

# 🔍 Check gene columns
cat("Gene columns preview:\n")
print(head(colnames(df_clean)))

# 🧹 Clean values BEFORE conversion
#df_clean[df_clean == ""] <- NA
#df_clean[df_clean == "NA"] <- NA

# 🔢 Convert to numeric safely (column-wise)
df_clean[] <- lapply(df_clean, function(x) {
  suppressWarnings(as.numeric(as.character(x)))
})

# Replace NA with 0
df_clean[is.na(df_clean)] <- 0

# ✅ Convert to matrix (cells x genes at this stage)
mat <- as.matrix(df_clean)

# 🔍 Check before transpose
cat("Before transpose (cells x genes):\n")
print(dim(mat))

# 🔁 TRANSPOSE → genes x cells
mat <- t(mat)

# ✅ Assign proper names
rownames(mat) <- colnames(df_clean)   # genes
colnames(mat) <- cell_names           # cells

# 🔍 Final checks
cat("After transpose (genes x cells):\n")
print(dim(mat))

cat("Gene names:\n")
print(head(rownames(mat)))

cat("Cell names:\n")
print(head(colnames(mat)))

# ✅ Final safety check
stopifnot(is.numeric(mat))
stopifnot(!any(is.na(mat)))

# ✅ Create Seurat object
obj <- CreateSeuratObject(counts = mat, project = "Sample_1")

print(obj)
head(obj)


###################################################
#Final
#####################3

library(Seurat)


library(dplyr)

# 📂 File paths
data_dir <- ""

files <- c(
  paste0(data_dir, "GSM2230757_human1_umifm_counts.csv.gz"),
  paste0(data_dir, "GSM2230758_human2_umifm_counts.csv.gz"),
  paste0(data_dir, "GSM2230759_human3_umifm_counts.csv.gz"),
  paste0(data_dir, "GSM2230760_human4_umifm_counts.csv.gz")
)

# 📦 Store Seurat objects
data_list <- list()

# 🔁 Loop through each file
for (i in seq_along(files)) {
  
  cat("\n============================\n")
  cat("Processing file:", files[i], "\n")
  
  # 📥 Read CSV
  df <- read.csv(gzfile(files[i]), header = TRUE, check.names = FALSE)
  
  # 👀 Inspect first file only
  if (i == 1) {
    print(head(df[,1:5]))
  }
  
  # Remove extra column
  df <- df[, !(colnames(df) %in% c("Unnamed: 0"))]
  df_clean <- df[, !(colnames(df) %in% c(
    "Unnamed: 0",
    "barcode",
    "assigned-cluster",
    "assigned_cluster",
    ""
  ))]
  
  
  # 🔤 Ensure unique gene names
  colnames(df_clean) <- make.unique(colnames(df_clean))
  
  # 🧬 Extract cell names (make unique)
  cell_names <- make.unique(as.character(df$barcode))
  
  # 🔢 Convert all columns to numeric safely
  df_clean[] <- lapply(df_clean, function(x) {
    suppressWarnings(as.numeric(as.character(x)))
  })
  
  # Replace NA with 0
  df_clean[is.na(df_clean)] <- 0
  
  # 🔁 Convert to matrix (cells x genes)
  mat <- as.matrix(df_clean)
  
  # 🔍 Check dimensions before transpose
  cat("Before transpose (cells x genes):", dim(mat), "\n")
  
  # 🔁 TRANSPOSE → genes x cells
  mat <- t(mat)
  
  # 🧾 Assign names
  rownames(mat) <- colnames(df_clean)   # genes
  colnames(mat) <- cell_names           # cells
  
  # 🔍 Final checks
  cat("After transpose (genes x cells):", dim(mat), "\n")
  cat("Gene preview:", head(rownames(mat)), "\n")
  
  # ✅ Safety checks
  stopifnot(is.numeric(mat))
  stopifnot(!any(is.na(mat)))
  stopifnot(!any(rownames(mat) == ""))
  
  # 🧪 Create Seurat object
  obj <- CreateSeuratObject(counts = mat, project = paste0("Sample_", i))
  
  # 📦 Store
  data_list[[i]] <- obj
}

# 🔗 Merge all samples
islet <- merge(
  x = data_list[[1]],
  y = data_list[2:length(data_list)],
  add.cell.ids = paste0("D", seq_along(data_list)))

# 👀 Check merged object
print(islet)

# 💾 Save
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
saveRDS(islet, "data/processed/islet_merged.rds")

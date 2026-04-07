dir()
setwd(dir = '/home/nishant/Desktop/RNA_linux/ImmunoIslet/')
# Create directory (if not exists)
dir.create("data/raw/islet_dataset", recursive = TRUE, showWarnings = FALSE)
setwd("data/raw/islet_dataset")
url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE84nnn/GSE84133/suppl/GSE84133_RAW.tar"
download.file(url, destfile = "GSE84133_RAW.tar", mode = "wb")
# Extract files
untar("GSE84133_RAW.tar")
# List extracted files
files <- list.files(pattern = "*.csv")
print(files)


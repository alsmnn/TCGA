# download TCGA data
#

download_tcga <- function(project_id) {
  # save a query in order to download all TCGA-data
  #
  query <- GDCquery(
    project = project_id,  # project-code(s) you want to download
    data.category = "Transcriptome Profiling",  # catogry for RNA-Seq Data
    data.type = "Gene Expression Quantification",  # raw reads or aligned data
    workflow.type = "HTSeq - Counts"  # FPKM / HTSeq-Count / etc.
  )
  ## download selected data
  #
  GDCdownload(
    query,  # name of the filtered dataset assigned above
    method = "api",  # needs to be set in order to download from the api
    files.per.chunk = 10,  # this should minimise prob. of corruption
    directory = "input/GDCdata"  # save files to a seperate directory
  )
}


# create a SummarizedExperiment from downloaded files
#

prepare_tcga <- function(project_id) {
  # save a query in order to download all TCGA-data
  #
  query <- GDCquery(
    project = project_id,  # project-code(s) you want to download
    data.category = "Transcriptome Profiling",  # catogry for RNA-Seq Data
    data.type = "Gene Expression Quantification",  # raw reads or aligned data
    workflow.type = "HTSeq - Counts"  # FPKM / HTSeq-Count / etc.
  )
  # create a SummarizedExperiment with counts and clinical data
  #
  tcga_data <- GDCprepare(
    query,
    directory = "input/GDCdata",
    summarizedExperiment = TRUE  # Downloads also clinical data
  )
  tcga_data
}


# create another round function (round up at 0.5)
#

round2 <- function(x, n = 0) {
  posneg <- sign(x)
  z <- abs(x) * 10^n
  z <- z + 0.5
  z <- trunc(z)
  z <- z / 10^n
  z * posneg
}


# create DSeqDataSet of SummarizedExperiments
#

dds_conversion <- function(tcga) {
  dds <- DESeqDataSet(
    tcga,
    design = ~race_wo
  )
  dds
}


# remove duplicate entries of genes, created by the conversion to a DDS
#

remove_duplicates <- function(dds1) {
  duplicates <- grep("\\.", rownames(dds1))
  dds <- dds1[-duplicates]
  dds
}


# we need a factor without spaces for DESeq analysis to work
#

race_wo <- function(raw_data) {
  race_wo_na <- raw_data[, !is.na(raw_data$race)]

  race_wo_na$race_wo <- NA

  if ("not reported" %in% race_wo_na$race) {
    nr <- race_wo_na$race == "not reported"
    race_wo_na$race_wo[nr] <- "not_reported"
  }

  if ("black or african american" %in% race_wo_na$race) {
    aa <- race_wo_na$race == "black or african american"
    race_wo_na$race_wo[aa] <- "black_or_african_american"
  }

  if ("white" %in% race_wo_na$race) {
    wh <- race_wo_na$race == "white"
    race_wo_na$race_wo[wh] <- "white"
  }

  if ("american indian or alaska native" %in% race_wo_na$race) {
    ai <- race_wo_na$race == "american indian or alaska native"
    race_wo_na$race_wo[ai] <- c("american_indian_or_alaska_native")
  }

  if ("asian" %in% race_wo_na$race) {
    as <- race_wo_na$race == "asian"
    race_wo_na$race_wo[as] <- "asian"
  }

  if ("native hawaiian or other pacific islander" %in% race_wo_na$race) {
    nh <- race_wo_na$race == "native hawaiian or other pacific islander"
    race_wo_na$race_wo[nh] <- c("native_hawaiian_or_other_pacific_islander")
  }
  race_wo_na
}


# create a DGEList element out of DESeqDataSet
#

normalize_tcga <- function(dds) {
  dgelist <- as.DGEList(dds)
  # Keep only genes expressed in >= 50% of the samples and recalculate libsize
  keep <- rowSums(cpm(dgelist) > 0.1) >= round2(nrow(colData(dds)) / 2)

  dgelist <- dgelist[keep, ,
    keep.lib.sizes = FALSE
  ]

  dgelist <- calcNormFactors(
    dgelist,
    method = "upperquartile",
    p = 0.75
  )

  dgelist <- estimateCommonDisp(
    dgelist,
    verbose = TRUE
  )
  dgelist
}

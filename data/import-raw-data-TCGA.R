# Setup -------------------------------------------------------------------

# setup your environment with all needed packages and configs
#

library(ProjectTemplate)
load.project(list(munging = FALSE, data_loading = FALSE, cache_loading = FALSE))


# Download and prepare tcga data ------------------------------------------

# get project ids for TCGA data
#

project_ids <- TCGAbiolinks:::getGDCprojects()$project_id

project_ids <- project_ids[grep("TCGA", project_ids)]

project_ids <- project_ids[order(project_ids)]

project_ids <- as.list(project_ids)


# download, save, and prepare summarized experiment for each entity
#

lapply(
 project_ids,
 download_tcga
 )


ProjectTemplate::cache("raw_tcga", depends = c("project_ids"), {
    raw_tcga <- lapply(
      project_ids,
      prepare_tcga
    )
    # assign names
    #
    names(raw_tcga) <- gsub("-", "_", project_ids)
  }
)

# Filtering and cleaning --------------------------------------------------

# we need a factor without spaces for DESeq analysis to work
#
# loop over every element in the tcga data and make the column for the DESeq
# analysis
#
tcga_raw <- lapply(raw_tcga, race_wo)


# DESeqDataSet (DDS) ------------------------------------------------------

# create names for the normalized data sets
#

tcga_names <- vector("character", length(tcga_raw))

for (i in seq_len(names(tcga_raw))) {
  tcga_names[i] <- paste(strsplit(names(tcga_raw[i]),
                                  "_")[[1]][2],
                         "se", sep = "_")
}


# create DSeqDataSet of SummarizedExperiments
#

tcga_raw_dds1 <- lapply(
  tcga_raw,
  dds_conversion
)


# remove duplicate entries of genes, created by the conversion to a DDS
#

tcga_raw_dds <- lapply(tcga_raw_dds1, remove_duplicates)


# assign names
#

names(tcga_raw_dds) <- paste(names(tcga_raw), "dds", sep = "_")



# DGEList -----------------------------------------------------------------

# create a DGEList element out of DESeqDataSet
#

tcga_raw_dgelist <- lapply(tcga_raw_dds, normalize_tcga)


names(tcga_raw_dgelist) <- paste(names(tcga_raw), "dgelist", sep = "_")



# SummarizedExperiment (se) -----------------------------------------------

# transforming DGEList object to RangedSummarizedExperiment (RSE) via
# DESeqDataSet
#

tcga <- vector("list", length(tcga_raw_dgelist))

for (i in seq_len(tcga_raw_dgelist)) {
  tcga[[i]] <- as.DESeqDataSet(tcga_raw_dgelist[[i]])
  # log2 transformation of cpm
  #
  assays(tcga[[i]])$log2ps_counts <- cpm(
    tcga_raw_dgelist[[i]]$pseudo.counts,
    log = TRUE
  )
  # if patients are not dead assign "days_to_last_follow_up" to the $surv_time
  #
  not_dead <- is.na(tcga[[i]]$days_to_death)
  dead <- !is.na(tcga[[i]]$days_to_death)

  if (any(notDead == TRUE)) {
    tcga[[i]]$surv_times[not_dead] <- tcga[[i]]$days_to_last_follow_up[not_dead]
    tcga[[i]]$surv_times[dead] <- tcga[[i]]$days_to_death[dead]
  }
  # assign TRUE and FALSE values to surv_events in order to make survival plots
  #
  tcga[[i]]$surv_events <- grepl(
    "dead|deceased",
    tcga[[i]]$vital_status,
    ignore.case = TRUE
  )
}


# assign names earlier created
#

names(tcga) <- tcga_names


# cache the normalized data sets.
ProjectTemplate::cache("tcga")
ProjectTemplate::cache("tcga_names")

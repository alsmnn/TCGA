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
    raw_tcga
  }
)


# create names for the normalized data sets
#

tcga_names <- vector("character", length(raw_tcga))

for (i in seq_along(names(raw_tcga))) {
  tcga_names[i] <- paste(strsplit(names(raw_tcga[i]),
                                  "_")[[1]][2],
                         "se", sep = "_")
}


# Preprocess TCGA data ----------------------------------------------------

# cache the normalized data sets.
ProjectTemplate::cache("tcga", depends = c("raw_tcga"), {
  tcga <- lapply(
    raw_tcga,
    munge_tcga
  )
  names(tcga) <- tcga_names
  tcga
}
)

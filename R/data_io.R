#' Load datasets in the demohaz R package
#'
#' Load datasets in the demohaz R package. This is done with a function to
#' support reformating of the raw data text files and to make the associated
#' tests clearer. The following datasets can be loaded:
#'
#' Ache       The Ache dataframe or lifetable
#' Hadza      The Hadza dataframe or lifetable
#' Hiwi       The Hiwi dataframe or lifetable
#' Tsimane    The Tsimane dataframe or lifetable
#'
#' @param dataset_name Name of the dataset to load
#' @param lifetable Whether to return the raw data (a dataframe) or the
#'   lifetable
#' @return The requested dataset
#' @export
#'
load_demohaz_data <- function(dataset_name,raw=FALSE) {
  if(tolower(dataset_name) == "ache") {
    file_path <- system.file("extdata",
                           "HillHurtadoLifeTable.txt",
                           package="demohaz")
    dataframe <- read.table(file_path,header=TRUE)
    if(raw) {
      return(dataframe)
    }
    lifetable <- demogR::life.table(x=dataframe$Age,
                                    nDx=dataframe$nDx,
                                    nKx=dataframe$nKx,
                                    type="cohort",
                                    iwidth=1,
                                    width12=c(1,1))
    return(lifetable)
  } else if (tolower(dataset_name) == "hadza") {
    file_path <- system.file("extdata",
                             "hadza.txt",
                             package="demohaz")
    dataframe <- read.table(file_path,header=TRUE)
    if(raw) {
      return(dataframe)
    }
    lifetable <- demogR::life.table(x=dataframe$Age,
                                    nDx=dataframe$nDx,
                                    nKx=dataframe$nKx,
                                    type="cd")
    return(lifetable)
  } else if (tolower(dataset_name) == "hiwi") {
    file_path <- system.file("extdata",
                             "hiwi.txt",
                             package="demohaz")
    dataframe <- read.table(file_path,header=TRUE)
    if(raw) {
      return(dataframe)
    }
    lifetable <- demogR::life.table(x=dataframe$Age,
                                    nDx=dataframe$nDx,
                                    nKx=dataframe$nKx,
                                    type="cd")
    return(lifetable)
  } else if (tolower(dataset_name) == "tsimane") {
    file_path <- system.file("extdata",
                             "newtsimane.txt",
                             package="demohaz")
    dataframe <- read.table(file_path,header=TRUE,skip=1)
    if(raw) {
      return(dataframe)
    }
    lifetable <- demogR::life.table(x=dataframe$Age,
                                    nDx=dataframe$nDx,
                                    nKx=dataframe$nKx,
                                    type="cd")
    return(lifetable)
  } else {
    stop(paste0("Unrecognized dataset = ",dataset_name))
  }

  stop("This point should never be reached.")
}
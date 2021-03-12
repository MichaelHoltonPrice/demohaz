#' Load datasets in the demohaz R package
#'
#' Load datasets in the demohaz R package. This is done with a function to
#' support reformating of the raw data text files and to make the associated
#' tests clearer. The following datasets can be loaded:
#'
#' Ache       The Ache lifetable (see demogR::life.table)
#' Hadza      The Hadza lifetable (see demogR::life.table)
#' Hiwi       The Hiwi lifetable (see demogR::life.table)
#' Tsimane    The Tsimane lifetable (see demogR::life.table)
#'
#' @param dataset_name Name of the dataset to load
#' @return The requested dataset
#' @export
#'
load_demohaz_data <- function(dataset_name) {
  if(tolower(dataset_name) == "ache") {
    file_path <- system.file("extdata",
                           "HillHurtadoLifeTable.txt",
                           package="demohaz")
    dataframe <- read.table(file_path,header=TRUE)
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
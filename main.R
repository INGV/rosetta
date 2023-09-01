#!/usr/bin/env Rscript

# This is an R implementation of the "Rosetta" algorithm described
# by Sbarra et al. (https://doi.org/10.1038/s41598-019-51966-4) to infer
# depth and magnitude of pre-instrumental earthquakes for Italy

# It takes in input a csv file with a list of DBMI earthquake IDs
# (https://doi.org/10.13127/DBMI/DBMI15.4), downloads related data
# and calculates hypocenter depth plotting mobile averaged macroseismic
# intensities in an epicentral distance of 0 to 55 kms, and saves both plot
# and average intensity tables.

### License Notice ###
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.#
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
### END License Notice ###

#
# Developed by Roberto Vallone <roberto.vallone@ingv.it>
# orcid:0000-0003-1208-9412
#

# This script is an attempt to rewrite in a cleaner and more rational way
# the original code used for the aforementioned publication

# version 2023-08-29


# Establishing filters -----------------------------------------------
eqIntNumb55Min <- 30 # Minimum number of intensities in 55 km from epicenter
maxStdError <- 0.01 # Maximum standard error admitted
mobAvgNumbMin <- 6 # Minimum total mobile averages number
meanv0_10IntMin <- 4 # Minimum average value of MCS
# between 0 and 10 kms (or fault R) from epicenter
resultAzMin <- 18 # Minimum number of azimuth slices of 10 degs to consider
# the macroseismic field spatially "homogeneous"

# Establishing whether to PLOT and/OR SAVE graphics ------------------
liPlotto <- T
liSalvo <- T
histYES <- T

# USER input/output preferences --------------------------------------
## Input file with list of earthquakes to download -------------------
# File must be a single column txt file with DBMI IDs one per row
# see file eventListTemplate.csv for an example
eqEventsFile <- "~/Dropbox/Projects/Rosetta/dbmiEventsList2.csv"

## output directory declaration --------------------------------------
#outdir <- "/home/roberto/Dropbox/Rosetta/outputDirectory"
outdir <-
  paste("../../",
        "OUTDIR_",
        as.character(format(Sys.time(),
                            "%Y_%m_%d_%H_%M_%S")),
        sep = "")
dir.create(outdir)

# Loading required libraries -----------------------------------------
#library(yesno)
library(tools)
library(ggplot2)
library(readr)
library(dplyr)
library(geosphere)
#library(plotrix)
library(REdaS)

# ASMI services URL
asmiServUrl <- "https://emidius.mi.ingv.it/services/macroseismic/query"

# Read event list
eqEventsList <- read.csv(eqEventsFile, header = F)

# set and create download directory
# downloadPath <- paste(dirname(eqEventsFile), "/tmpFILEZdbmi2", sep = "" )
downloadPath <- paste("/tmp/rrosettaTmpFilezDBMI", as.character(format(
  Sys.time(),
  "%Y_%m_%d_%H_%M_%S"
)), sep = "")

dir.create(downloadPath)

# download requested files ---------------------------------------------
scaricali <- function(inputID) {download.file(
  paste(asmiServUrl,
        "?eventid=",
        inputID,
        "&includemdps=true&format=textmacro", sep = ""),
  destfile = paste(downloadPath, "/", inputID, "_DBMI.csv", sep = ""))
  Sys.sleep(0.5)
}

sapply(eqEventsList, scaricali)

eqFiles <-
  list.files(path = downloadPath,
             pattern = "\\_DBMI.csv",
             full.names = TRUE)

# Create output CSV file in the output directory
outFileName <-
  paste(outdir,
        "/rosetta-DBMI_",
        as.character(format(Sys.time(), "%Y_%m_%d_%H_%M_%S")),
        ".csv",
        sep = "")

file.create(outFileName)

# Inserting filter variable in a pre-header
preHeader <- c(
  "## Filters set up ##",
  paste(
    "# Minimum number of intensities in 55 kms from epicenter",
    eqIntNumb55Min,
    sep = ": "
  ),
  paste("# Maximum standard error admitted",
        maxStdError,
        sep = ": "),
  paste("# Minimum total mobile averages number",
        mobAvgNumbMin,
        sep = ": "),
  paste(
    "# Minimum average value of intensities in the first mobile window from epicenter",
    meanv0_10IntMin,
    sep = ": "
  )
)

# Creating header
header <- list(
  "DBMI Id",
  "Epicentral Area",
  "Date",
  "Mw",
  "Instrumental Depth",
  "Total MDPs",
  "MDPs within 55 km",
  "Mobile avgs within 55 km",
  "Epicenter Lon",
  "Epicenter Lat",
  "Steepness",
  "Calculated Depth",
  "Std. error",
  "Y-axis intercept",
  "Y-axis corrected intercept",
  "Mw_Y intercept",
  "Mw_I corrected intercept",
  "N. azimuth slices",
  "Fault radius (W&C)",
  "Plot file"
)
writeLines(preHeader, outFileName)
# Insert blank line between pre-header and header
cat("\n", file = outFileName, append = TRUE)
write.table(
  header,
  file = outFileName,
  sep = ",",
  quote = FALSE,
  col.names = FALSE,
  row.names = FALSE,
  append = TRUE,
)

# Load rosetta cycle source -------------------------------------------
source("rosettaCycle.R")

# (l)Apply ROSETTA FUNCTION to the dataset
lapply(eqFiles, exe.rosetta.cycle)
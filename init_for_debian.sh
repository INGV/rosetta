#!/bin/sh

CRAN_MIRROR='http://cran.mirror.garr.it/CRAN/'

sudo apt-get update
sudo apt-get install --yes r-base r-base-dev
for pkg in tools ggplot2 readr dplyr geosphere REdaS ;
do
	Rscript -e "install.packages('${pkg}', repo='${CRAN_MIRROR})"
done

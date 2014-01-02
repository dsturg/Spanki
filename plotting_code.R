
################################################
################################################
#  Splicing Analysis Kit (Spanki)
#
#  Code for producing plots for visualizing splicing differences
#  Requires the 'vcd' library (Author/Maintainer - David Myer):
#	http://cran.r-project.org/web/packages/vcd/
#
################################################
################################################


#  Load the library
#  Make sure you have "vcd" installed
#  See:
#  http://cran.r-project.org/web/packages/vcd/index.html
#  David Meyer, Achim Zeileis, and Kurt Hornik (2012). vcd: Visualizing Categorical Data. R package version 1.2-13.

library(vcd)


################################################
################################################
#
#  Fourfold. mosaic plots for 2x2's
#
#
################################################
################################################


#///////////////////
# Basic fourfold and mosaic plots
#//////////////////

# Consider a splicing event with the following data:
# Sample 1:
#	Inclusion = 9
#   Exclusion= 30
# Sample 2:
#	Inclusion = 132
#   Exclusion = 157

# Put the data in a 2x2 matrix:
mymat <- matrix(c(20,50,45,40),ncol=2,byrow=T)

# Perform Fisher's Exact Test:
fisher.test(mymat)

# Make a fourfold plot:
fourfoldplot(mymat)

# Make a mosaic plot
mosaic(mymat)

# Alternative fill colors for mosaics:

fill_colors <- matrix(c("dark cyan", "gray", "gray", "dark magenta"), ncol = 2)
mosaic(mymat, gp = gpar(fill = fill_colors, col = 0))
fill_colors <- matrix(c("blue", "blue", "red", "red"), ncol = 2)
mosaic(mymat, gp = gpar(fill = fill_colors, col = 0))


#///////////////////
# Example mosaic plot 
# using sample data output
#//////////////////

# Set the directory to 
setwd("~/Desktop/spankitest/")
library(vcd)

splicecomp <- read.delim(file="F_vs_M_splicecomp/event_compare.out",stringsAsFactors=F)

i <- 1

mymat <- matrix(c(splicecomp$inc1[i],splicecomp$exc1[i],splicecomp$inc2[i],splicecomp$exc2[i]),ncol=2,byrow=T)

# Make a mosaic plot
fill_colors <- matrix(c("blue", "blue", "red", "red"), ncol = 2)
mosaic(mymat, gp = gpar(fill = fill_colors, col = 0))






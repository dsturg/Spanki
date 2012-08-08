
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
Example data
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

#\\\\\\\\\\\\\\\\\\
Variations
#\\\\\\\\\\\\\\\\\\\


# Alternative fill colors for mosaics:

fill_colors <- matrix(c("dark cyan", "gray", "gray", "dark magenta"), ncol = 2)
fill_colors <- matrix(c("blue", "blue", "red", "red"), ncol = 2)

mosaic(mymat, gp = gpar(fill = fill_colors, col = 0))





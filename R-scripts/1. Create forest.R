# R-script to create a spatially explicit representation of a forest, 
# based on rank abundance data provided at specific locations. 
# The modelled area size can be adjusted with the parameters area_length
# and area_width (both in m). Each 1 x 1 m patch can hold up to one tree. 

# We start with an empty space containing area_length x area_width patches.
# At specified locations (specified in the input file 'plot_locations.txt'), 
# sources of tree seeds can be found. These sources contain a number of 
# individuals of each species, that is specified in the input file 
# 'RAI_per_plot.txt'. 
# Each patch can become occupied by a tree with probability (Pocc). 
# The probability that a seed comes from a specific source population
# depends on the relative distance of the source to the patch (i.e. the 
# probability is highest that a seed comes from the nearest source population).
# Every tree in a source population has an equal change of colonizing the
# empty patch. 




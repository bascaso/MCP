#  Begoña Ascaso. January 2017

#  DESCRIPTION:

#  MCP (Matching Catalogues and Computing Completeness and Purity rates)
#  The MCP is a python package that matches one list of structures (e.g dark matter haloes) to other (e.g. cluster detections) using Cylindrical Matching and FoF+Cylindrical Matching. 
#  Additionally, if requested, the program will compute completeness and purity rates.
#  Note that before producing the completeness and purity plots, we need to match two ways.



#  CONTENT:

#  The package includes the following programs:

#  mainMCP.py (main program to be called with the MCP.pars file)
#  MatchB.py (program containing subroutines to match using the CM and the FCM methods)
#  CompPurB.py (program that computes the Completeness and Purity rate as a function of redshift and mass/richness and creates plots)
#  funtionsB.py (auxiliar program containing help functions)

#  Additionally:

#  PurityvsCompPlot.py (program that will read different matched files and will compute Purity vs Completeness for all of them, to be modified manually).


# USE

1. python mainMCP.py MCP.pars

Need to modify MCP.pars according to our needs.


if MATCHING=yes -> Will match Haloes to Detections (H2D) or Detections to Haloes (D2H) with a Cylindrical Matching (CM) method or FoF+CM (FCM) method. It will save files with the matches in both directions. A photo-z scatter needs to be specified in the MCP.pars parameter

The input of the files need to be 

Id, ra, dec, z, number of galaxies (or richness), halo mass (or SNR), radius


if COMPPUR=yes -> It will compute completeness and purity curves as a function of redshift and mass/richness and it will make plots of them. Note that to do this you need to have matched the files both ways first.




2. python PurityvsCompPlot.py 

It will create PurityvsComp plots using the matchings previsouly made with a particular method (FCM) or CM and making a rough mass-observable fit.

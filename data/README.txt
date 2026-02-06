                            CONTENTS OF THIS ARCHIVE
                            ========================
                         	
  This archive contains observed values and projections of precipitation
indices for 16 regions of the UK. Filenames are as follows:
	
  - HadUK.csv: observed values, derived from the HadUK gridded dataset. 
    Columns are:
    . Year, from 1892 to 2021. Years run from December to November, so 
      that (for example) "1970" refers to the period from December 1969
      to November 1970.
    . Season: this is either "Annual", "DJF", "MAM", "JJA" or "SON".
    . Region: this is the name of one of the 16 UK regions. Note, 
      however, that the values for the Channel Islands are all NA  
      (this is a problem with the HadUK dataset itself, presumably
      because it's a land surface dataset and someone at the Met 
      Office forgot about the Channel Islands when preparing it). 
    . fwd: the mean fraction of wet days for the period, over all grid 
      cells within the region. For an individual grid cell, a wet day 
      is defined as a day with at least 1mm of rainfall. 
    . wsmax: the mean, over all grid cells within the region, of the 
      maximum duration of wet spells starting within the given period.
      For an individual grid cell, a wet spell is a sequence of 
      consecutive wet days preceded and followed by a dry day.

  - UKCP18_##.csv: projections from the UKCP18 ensemble for the period
    1981-2080. These are generated using historical greenhouse gas emissions 
    until 2015 and the RCP8.5 scenario thereafter. "##" in the filename
    is the number of the ensemble member. The file format is the same as  
    for HadUK.csv, except that data for the Channel Islands are present
    here. 

  - CORDEX_$$$_%%%_@@@.csv: projections from the EuroCORDEX ensemble for the 
    period 1981-2080. Again, these are generated using historical greenhouse
    gas emissions until 2015 and the RCP8.5 scenario thereafter. In the 
    filenames, "$$$" and "%%%" denote the GCM and RCM used respectively
    while "@@@" denotes the run identifier for that particular GCM:RCM
    combination. The file format is the same as for HadUK.csv. Note that
    the different RCMs produce outputs on different grids: in each case, 
    the fractions of wet days and maximum wet spell durations are 
    initially calculated on the RCMs' "native" grids before averaging over
    regions, with the corresponding "region masks" being defined separately
    for each RCM. 

  - regionmask-region_osgb.nc: this contains the "mask" information for the
    16 regions, in the grid used by the HadUK gridded dataset. To read the
    information, do something like the following:
    
    	library(RNetCDF)
    	ncFile <- open.nc("regionmask-region_osgb.nc")
    	RegNames <- var.get.nc(ncFile, "geo_region") # Region names 
    	RegMasks <- var.get.nc(ncFile, "region_mask")
    	close.nc(ncFile)
    	
    RegNames is now a character vector of length 16, while RegMasks is an array
    of dimension 82x112x16 containing zeroes and 1s: each of the 16 "slices"
    is thus a matrix of dimension 82x112, containing 1s in the grid cells 
    in the corresponding region. A quick way to plot the regions is
    
    	MapData <- apply(RegMasks, MARGIN=1:2, FUN=function(x) sum(x*(1:16)))
    	MapData[MapData==0] <- NA
    	image(MapData, col=hcl.colors(16, palette="Set 2"))
    	
    							Richard Chandler
    							19th November 2025



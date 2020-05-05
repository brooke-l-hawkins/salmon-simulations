#*******************************************************************************
#
# An individual-based model for evaluating 
#   how phenology, growth and survival of juvenile Chinook salmon may respond 
#   to altered thermal regimes and the predation by largemouth bass in the
#   Snoqualmie River watershed, WA, USA.
#
# Coded by:
#    B.L. Hawkins (general, second species, predation), A.H. Fullerton (general), 
#    B.J. Burke (movement), N. Som (network shapes) & M. Nahorniak (bioenergetics) 
#    Last Updated 30 Apr 2020
#*******************************************************************************

# clear workspace
rm(list=ls())
gc()

#=== SCENARIOS =================================================================
climate.year = 2014 # 2014 for Baseline, 2015 for Warm
predation.flag = FALSE # FALSE for Baseline, TRUE for Predator

#=== SENSITIVITY ANALYSIS ======================================================
SA = FALSE # TRUE to run sensitivity analysis, FALSE to run scenarios only
iter.list = 1 # how many replicates (iterations) to run, 1:300 for sensitivity analysis

#=== SETUP =====================================================================

start.time = proc.time()

# Libraries
library(SSN)
library(plyr)
library(rgdal)
library(raster)
library(RColorBrewer)
library(GLDEX)

# Load Functions
source("code/functions.R")

# Global Variables
salmon.nm = "chinook"
fish_other.nm = "lmb"
netnm = "sno" # used for plotting and some functions
first.date = as.Date(paste0(climate.year-1, "-", "09-01")) #starting date for simulation and for spawning
last.date = as.Date(paste0(climate.year, "-", "08-31")) #last date of the simulation
dat.idx = seq(from = first.date, to = last.date, by = 1) # list of all dates to model
day1 = strptime(dat.idx[1], format="%Y-%m-%d")
spawn.date.variable = TRUE # variable spawn date? (set to F to have all fish spawn on first day)
ifelse(predation.flag == TRUE, tag <- "P", tag <- "A") #Predation or Alone
ifelse(tag == "P", SecondSpecies<- TRUE, SecondSpecies<- FALSE)
iter = 1
plot.iter = 1 # which iteration will have maps plotted for each time step
run = fncGetRun()
show.progress = TRUE # send statements to the console showing progress


# Network-specific Fields
  length.field = "LENGTHKM"
  so.field = "SO"
  width.field = "WIDTH_M"
  wt.field = "WT"
  xlat = "NEAR_X"
  ylon = "NEAR_Y"
  lng2b.field = "lngth2B" #this is in fish shapefiles

# Directories
mydir = getwd()
loadDir = "data.in"
ssn.folder = "sno.ssn"
outputDir = paste0("data.out/", tag, ".", climate.year) 
# make directories, if they don't already exist
if (! dir.exists(file.path(paste0(mydir, "/", outputDir)))) {
  dir.create(file.path(paste0(mydir, "/", outputDir)))
}

# Plotting
# If plot.flag set to TRUE, first iteration will plot fish locations throughout time
plot.flag = TRUE
track.individuals = FALSE
if (plot.flag == TRUE) {
  imageDir= paste0(outputDir, "/images")
  
  if (!dir.exists(file.path(imageDir))) {
    dir.create(file.path(imageDir))
  }
  
  # Color palette for plotting
  grays = brewer.pal(9, "Greys")
  grays = grays[c(9,8,7,6,5,4,3,2,1)]
  
  # water temperature
  cb=  fncColBrewPlus(n = 14, paint = F)
  left = c(0,2,4,6,8,10,12,14,16,18,20,22,24)
  rght = c(2,4,6,8,10,12,14,16,18,20,22,24,27)
  
  # Read in basin outline or hillshade for plotting
  outline = TRUE
  if (outline) {
    basin = readOGR(paste0(mydir,"/",loadDir,"/shapefiles"),"Basin_snq")
  } else {
    hs = raster(paste0(mydir, "/", loadDir, "/shapefiles/snoq_hs"))
    units = km
    projstr = CRS("+proj=aea +lat_1=47.66564 +lat_2=47.50781 +lat_0=23.0 +lon_0=-122.4525 +x_0=0 +y_0=0")
    basin = projectRaster(hs, crs = projstr)
  }
  
  streamlines = TRUE
  if(streamlines) streams<- readOGR(paste0(getwd(), "/", loadDir, "/", ssn.folder), "edges")
  
  # Set extent for plotting (will be updated later once ssn is loaded)
  ex = extent(basin)
}

# Load model parameters
parameters.all = read.csv(paste0(loadDir, "/parameters/parameters.csv"), as.is = TRUE)
parameters = cbind.data.frame(rbind(salmon.nm, fish_other.nm), t(parameters.all[, c(salmon.nm, fish_other.nm)]), stringsAsFactors = FALSE)
colnames(parameters) = c("species",parameters.all$parameter); row.names(parameters) = c("salmon", "fish_other")
for(c in 7:ncol(parameters)){ parameters[, c] = as.numeric(parameters[,c])}
for(i in c("spawn.date.begin", "om.date.taper", "om.date.end")){ parameters[, i] = gsub("x", "", parameters[, i])}
parameters.static = parameters

# Species-specific constants used in Wisconsin bioenergetics model
salmon.constants = fncReadConstants("salmon")
if(SecondSpecies == TRUE) fish_other.constants = fncReadConstants("fish_other")

# Load Attribute, Network, and Fish growth data
  # Load pre-calculated SSN water temperature data for the climate scenario
  WT.df = read.csv(paste0(loadDir, "/", ssn.folder, "/WT.df", climate.year, ".csv"), header = TRUE, stringsAsFactors = FALSE)

# Load copy of SSN (will be kept as loaded backup)
sno.ssn = importSSN(paste0(loadDir, "/", ssn.folder) , predpts='preds')
  # There are 15 points where ratio is >1; checked in GIS, this changes them to 1:
  sno.ssn@predpoints@SSNPoints[[1]]@point.data$ratio[sno.ssn@predpoints@SSNPoints[[1]]@point.data$ratio > 1] = 1
  # Reaches with no prediction points
  missed = c(0, which(!1:nrow(sno.ssn@data)%in%sort(unique(sno.ssn@predpoints@SSNPoints[[1]]@point.data$rid))))
  network.base.segs = 1:5
# Create vector of all reaches
allsegs = sno.ssn@data$rid

# Load topology lists for the network (previously created)
# all segments downstream from a given seg
load(paste0(mydir, "/", loadDir, "/", ssn.folder, "/dnsegs.RData"))
# all segments upstream from a given seg
load(paste0(mydir, "/", loadDir, "/", ssn.folder, "/upsegs.RData"))
# all segments at the upper end/confluence of a given seg
load(paste0(mydir, "/", loadDir, "/", ssn.folder, "/jct.list.RData"))

# Load pre-calculated growth for each temperature/ration/size combination 
# that can be looked up instead of having to run the bioenergetics model on the fly each time
# dims = no. of WT (500; 0 to 25 C), ration (151; 0.01 to 0.3), weight (500, 1 to 500 g)
load(paste0(mydir, "/", loadDir, "/fish.growth.lookup/wt.growth.array.", parameters["salmon","species"], ".RData"))
assign("wt.growth.salmon", wt.growth)
if(SecondSpecies == TRUE){
  load(paste0(mydir, "/", loadDir, "/fish.growth.lookup/wt.growth.array.", parameters["fish_other","species"], ".RData"))
  assign("wt.growth.fish_other", wt.growth)
}
rm(wt.growth)

  
# Set up arrays that will store data
# List of tracked fields to store permanently
array.cols2keep = c("pid", "seg", "xloc", "yloc", "upDist", "SO", "WT", "TU","emrg", "survive", 
                    "conspecificdensity", "otherdensity", "totaldensity", "direction", "movedist", 
                    "pvals", "consInst", "consCum", "growth", "weight", "ration", "om.prob",
                    "pred.prob", "num.prey", "num.eaten")
output.cols2keep = c("pid", "TU","emrg", "survive", "consCum", "weight", "dateSp", "dateEm", "dateOm", "datePr", "dateDi")

# Salmon arrays
nFish = parameters["salmon", "nFish"]
# no. of fish, no. variables, no. time steps, no. iterations
salmon.array = array(NA, dim = c(nFish, length(array.cols2keep), length(dat.idx) * 2, length(iter.list))) 
# no. of fish, no. variables, no. iterations; final values so no time component
salmon.finalstep = array(NA, dim = c(nFish, length(output.cols2keep), length(iter.list))) 

# fish_other arrays
if(SecondSpecies == TRUE){
  nFish = parameters["fish_other", "nFish"]
  # no. of fish, no. variables, no. time steps, no. iterations
  fish_other.array = array(NA, dim = c(nFish, length(array.cols2keep), length(dat.idx) * 2, length(iter.list))) 
  # no. of fish, no. variables, no. iterations; final values so no time component
  fish_other.finalstep = array(NA, dim = c(nFish, length(output.cols2keep), length(iter.list))) 
}
  
#=== ITERATE (if running more than one replicate) ==============================

# Start looping through iterations
for(iter in iter.list){
  start.time = proc.time()
  (scenario = paste0(climate.year, ".", iter))
  
  # Seed for reproducible results, change for each iteration
  set.seed(iter)
  
  # if running sensitivity analyses:
  if(SA == TRUE){
    parameters = parameters.static
    parameters = fncGetParametersSA(vary_nFish = TRUE)
  }
  
  # Salmon arrays
  nFish = parameters["salmon","nFish"]
  # no. of fish, no. variables, no. time steps, no. iterations
  salmon.array = array(NA, dim = c(nFish, length(array.cols2keep), length(dat.idx) * 2)) 
  # no. of fish, no. variables, no. iterations; final values so no time component
  salmon.finalstep = array(NA, dim = c(nFish, length(output.cols2keep))) 
  
  # fish_other arrays
  if(SecondSpecies == TRUE){
    nFish = parameters["fish_other", "nFish"]
    # no. of fish, no. variables, no. time steps, no. iterations
    fish_other.array= array(NA, dim = c(nFish, length(array.cols2keep), length(dat.idx) * 2)) 
    # no. of fish, no. variables, no. iterations; final values so no time component
    fish_other.finalstep = array(NA, dim = c(nFish, length(output.cols2keep))) 
  }
  
#=== INITIALIZE HABITAT ========================================================

  # Load the Snoqualmie Spatial Stream Network (SSN)
  ssn = sno.ssn

  if(plot.flag == TRUE){ #set plotting extent based on SSN bounding box
    ex@xmin = ssn@bbox[1]
    ex@xmax = ssn@bbox[3]
    ex@ymin = ssn@bbox[2]
    ex@ymax = ssn@bbox[4]
  }
  
  # Get reach widths, adjusted to represent the proportion useable by fish, which decreases as width increases
  useable.widths = fncUseableWidths(dat = ssn@data[,c("rid", width.field)])
  ssn@data$UseableWidth = useable.widths[useable.widths[,1] == ssn@data$rid,3]
  
  # Extract dataset (NHD networks have variable numbers)
  wq.df = getSSNdata.frame(ssn, "preds")
  # change factor fields to numeric, if necessary
  if(is.factor(wq.df$rid)) wq.df$rid = as.numeric(levels(wq.df$rid)[as.numeric(wq.df$rid)])

  # Get lists of network-specific segments above and below Snoqualmie Falls (a natural anadromous barrier) and Tolt Reservoir dam
  segsAbvFalls= upsegs[[53]]
  segsAbvRes = upsegs[[308]]
  segsBlwFalls = unique(ssn@data$rid)[! unique(ssn@data$rid) %in% segsAbvFalls]
  segsBlwRes = unique(segsBlwFalls)[! unique(segsBlwFalls) %in% segsAbvRes]
  segsAccessible = intersect(segsBlwFalls, segsBlwRes)
  
  # Add column to denote which reaches are accessible to fish (for blocking movement)
  ssn@data$accessible = 0
  ssn@data$accessible[ssn@data$rid %in% segsAccessible] = 1 
  ssn@data$Chinook_rg[is.na(ssn@data$Chinook_rg)] = 0
  ssn@data$Chinook_sp[is.na(ssn@data$Chinook_sp)] = 0
  if(SecondSpecies == TRUE){
    ssn@data$LMB_range[is.na(ssn@data$LMB_range)] = 0
    ssn@data$LMB_sp[is.na(ssn@data$LMB_sp)] = 0
  }
  
  # Calculate ration (g/g/d) available to fish; More productive food webs in lower mainstem habitats
  #linearly relate to log(stream order) and position within reach
  ration_1 = fncRescale(log(wq.df[,so.field] + (1 - wq.df$ratio)), c(parameters["salmon","ration.lo"], parameters["salmon","ration.hi"]))
  #inverse-linearly relate to distance from the mouth (upDist)
  ration_2 = fncRescale(wq.df$upDist, c(parameters["salmon","ration.hi"], parameters["salmon","ration.lo"])) 
  wq.df$ration = ssn@predpoints@SSNPoints[[1]]@point.data$ration = (ration_1 + ration_2) / 2

  # Save this 'base case' ration because we will be updating ration later
  wq.df$ration_base = ssn@predpoints@SSNPoints[[1]]@point.data$ration_base = wq.df$ration

  if(SecondSpecies == TRUE){
    # Add ration field for second species
    ration_1 = fncRescale(log(wq.df[,so.field] + (1 - wq.df$ratio)), c(parameters["fish_other","ration.lo"], parameters["fish_other","ration.hi"]))
    ration_2 = fncRescale(wq.df$upDist, c(parameters["fish_other","ration.hi"], parameters["fish_other","ration.lo"])) 
    wq.df$ration_ss = ssn@predpoints@SSNPoints[[1]]@point.data$ration_ss = (ration_1 + ration_2) / 2
    wq.df$ration_base_ss = ssn@predpoints@SSNPoints[[1]]@point.data$ration_base_ss = wq.df$ration_ss
  }
  
  # Initialize fish density
  wq.df$density = ssn@predpoints@SSNPoints[[1]]@point.data$density = 0
  
  # Reduce field set
  if(SecondSpecies == TRUE) ration.cols = c("ration","ration_ss","ration_base","ration_base_ss") else ration.cols = c("ration","ration_base")
  wq.df = wq.df[, c("pid", "rid", so.field, xlat, ylon, "ratio", "upDist", "afvFlow", "locID", "netID", width.field, ration.cols, "density")]
  nwq = nrow(wq.df)
  
  # Put wq.df back into ssn
  ssn = putSSNdata.frame(wq.df, ssn, "preds")
  
#=== INITIALIZE FISH ===========================================================

  # 1. Import the fish locations into the SSN
  
  # Salmon
  fish.shp = parameters["salmon", "fish.shp"] # get name of shapefile to load
  ssn = importPredpts(ssn, fish.shp, "ssn") # load into SSN
  salmon.id = 2 # this is the 2nd 'preds' file loaded in the SSN
  ssn@predpoints@ID[salmon.id] = parameters["salmon","species"] # name it

  # Update fish data files with some basic tracking info and fields that will hold scenario data
  # rename fish ID sequentially (critical for tracking movement later)
  ssn@predpoints@SSNPoints[[salmon.id]]@point.data$pid = 
    rownames(ssn@predpoints@SSNPoints[[salmon.id]]@point.data) = 1:nrow(ssn@predpoints@SSNPoints[[salmon.id]]@point.data) 
  
  # create 'seg' field and set equal to 'rid' field
  ssn@predpoints@SSNPoints[[salmon.id]]@point.data$seg = ssn@predpoints@SSNPoints[[salmon.id]]@point.data$rid
  
  # create 'wt.field'
  ssn@predpoints@SSNPoints[[salmon.id]]@point.data[,wt.field] = 0 # water temperature field


  # Other fish
  if(SecondSpecies == TRUE){
    fish.shp = parameters["fish_other", "fish.shp"] # get name of shapefile to load
    ssn = importPredpts(ssn, fish.shp, "ssn") # load into SSN
    fish_other.id = 3 # this is the 3rd 'preds' file loaded in the SSN
    ssn@predpoints@ID[fish_other.id] = parameters["fish_other","species"] # name it

    rm(fish.shp)
    
    # Update fish data files with some basic tracking info and fields that will hold scenario data
    # rename fish ID sequentially (critical for tracking movement later)
    ssn@predpoints@SSNPoints[[fish_other.id]]@point.data$pid = 
      rownames(ssn@predpoints@SSNPoints[[fish_other.id]]@point.data) = 1:nrow(ssn@predpoints@SSNPoints[[fish_other.id]]@point.data) 
    
    # create 'seg' field and set equal to 'rid' field
    ssn@predpoints@SSNPoints[[fish_other.id]]@point.data$seg = ssn@predpoints@SSNPoints[[fish_other.id]]@point.data$rid
    
    # create 'wt.field'
    ssn@predpoints@SSNPoints[[fish_other.id]]@point.data[,wt.field] = 0 # water temperature field

  }

  # 2. Create fish data frame for manipulating during simulation
  # pid: unique identification for each fish
  # seg: which network segment the fish is in (same as rid)
  # ratio: relative distance along the segment (from bottom to top [0-1])
  # xloc: x coordinate of fish
  # yloc: y coordinate of fish
  # length2segBase is the length from a fish's position in the current segment to its base
  # upDist is the distance from a fish's position to the base of the whole network
  # emrg values: 0 = egg placed (initial value), -1 = egg spawned, 1 = egg emerged
  # survive values: 1 = alive (initial value), 0 = dead, 2 = outmigrated/smolted, -1 = eaten by predator
  
  # Get data frame from SSN
  salmon.df = getSSNdata.frame(ssn, "chinook")[c(xlat, ylon, "pid", "rid", "ratio", "upDist", lng2b.field, wt.field)]

  # Create new dataframe
  salmon = data.frame(pid = 1:nrow(salmon.df), seg = 0, xloc = 0, yloc = 0, ratio = 0, length2segBase = 0, upDist = 0, 
          SO = 0, WT = 0, TU = 0, emrg = 0, survive = 1, conspecificdensity = parameters["salmon", "initial.mass"], 
          otherdensity = 0, totaldensity = 0, direction = 0, movedist = 0, pvals = 0, ration = 0, consInst = 0, consCum = 0,  
          growth = 0, weight = parameters["salmon", "initial.mass"], dateSp = day1, dateEm = day1, dateOm = day1,
          datePr = day1, dateDi = day1, om.prob = 0, pred.prob = 0, num.prey = 0, num.eaten = 0, stringsAsFactors = FALSE)

  # Fill in what we can from the ssn
  salmon$seg = salmon.df$rid
  salmon$ratio = salmon.df$ratio
  salmon$xloc = salmon.df[,xlat]
  salmon$yloc= salmon.df[,ylon]
  salmon$length2segBase = salmon.df[,lng2b.field]
  salmon$upDist = salmon.df$upDist
  rm(salmon.df)
  
  if(SecondSpecies == TRUE){
    fish_other.df = getSSNdata.frame(ssn, "lmb")[c(xlat, ylon, "pid", "rid", "ratio", "upDist", lng2b.field, wt.field)]
    fish_other = data.frame(pid = 1:nrow(fish_other.df), seg = 0, xloc = 0, yloc = 0, ratio = 0, length2segBase = 0, upDist = 0, 
            SO = 0, WT = 0, TU = 0, emrg = 0, survive = 1, conspecificdensity = parameters["fish_other", "initial.mass"], 
            otherdensity = 0, totaldensity = 0, direction = 0, movedist = 0, pvals = 0, ration = 0, consInst = 0, consCum = 0,  
            growth = 0, weight = 0, dateSp = day1, dateEm = day1, dateOm = day1, datePr = day1, dateDi = day1, om.prob = 0, 
            pred.prob = 0, num.prey = 0, num.eaten = 0, stringsAsFactors = FALSE)
    fish_other$seg = fish_other.df$rid
    fish_other$ratio = fish_other.df$ratio
    fish_other$xloc = fish_other.df[,xlat]
    fish_other$yloc= fish_other.df[,ylon]
    fish_other$length2segBase = fish_other.df[,lng2b.field]
    fish_other$upDist = fish_other.df$upDist
    set.seed(iter)
    fish_other$weight = sample(x = seq(1:parameters["fish_other","max.initial.mass"]), size = parameters["fish_other","nFish"], 
              prob = sort(fncRescale(rlnorm(parameters["fish_other","max.initial.mass"],1,parameters["fish_other", "mass.shape"])), decreasing = TRUE))
    rm(fish_other.df)
  }
  
  # 3. Initialize variables and determine timing for salmon spawning
  # initialize spawning start date and counter
  spawn.date.begin = parameters["salmon", "spawn.date.begin"]
  spawn.init = as.Date(paste0(climate.year - 1, "-", spawn.date.begin))
  spawn.counter = 1
  # determine spawning window and distribution
  # this only applies if salmon spawning is variable (over multiple days) instead of fixed (on one day).
  # spawning window lasts a mean of 51 days (31-71 days range) based on empirical data from WDFW spawner surveys 2005-2015.
  # spawning distribution is modeled as a normal distribution:
  # - number of observations is the number of simulated salmon
  # - mean is zero
  # - standard deviation is parameter "spawn.date.shape"
  # spawning times are the standard distribution broken up into parameter "nSpawnDays" * 2,
  # to represent two modeled time steps per calendar day.
  # salmon.spawn.times is a vector of numbers, where length is the number of timesteps in the spawning window.
  # each element represents the number of salmon that will spawn for a given timestep,
  # the first element is the number of salmon that will sapwn on the first timestep of spawn.init
  set.seed(iter)
  salmon.spawn.times = histsu(rnorm(nrow(salmon), mean = 0, sd = parameters["salmon", "spawn.date.shape"]), breaks = (parameters["salmon", "nSpawnDays"] * 2), plot = FALSE)$counts
  
#=== RUN SIMULATION ============================================================

  # Loop over each day (dd) and time step (tt) in an iteration (iter)
  dt = 1 # day/time counter
  for(dd in 1:length(dat.idx)){
    for(tt in c(6, 18)){ #6am, 6pm
      set.seed(iter)
      thetitle = fncGetTitle(dat.idx[dd], tt)
      cat(iter, ": ", thetitle, "runtime (h): ", proc.time()[3]/60/60, "\n")
      if(show.progress == TRUE){
        cat(length(salmon$emrg[salmon$emrg == -1]), " eggs incubating", "\n")
        cat(length(salmon$weight[salmon$survive == 2]), " subyearling smolts: ", quantile(salmon$weight[salmon$survive == 2]), "\n")
        cat(length(salmon$weight[salmon$survive == 1 & salmon$emrg == 1]), " rearing in stream: ", quantile(salmon$weight[salmon$survive == 1]), "\n")
        cat("proportion yearlings: ", length(salmon$weight[salmon$survive == 1]) / (length(salmon$weight[salmon$survive == 2]) + length(salmon$weight[salmon$survive == 1])), "\n")
        cat("fry-to-smolt survival: ", length(salmon$weight[salmon$survive > 0]) / parameters["salmon","nFish"], "\n")
        if(SecondSpecies == TRUE) cat(length(fish_other$weight[fish_other$survive == 1]), " fish_other: ", 
                                      quantile(fish_other$weight[fish_other$survive == 1]), "\n")
      }
      
  # 1. UPDATE HABITAT: Get predicted stream temperatures
      
      # Load temperature for this date and time
      if("WT" %in% colnames(ssn@predpoints@SSNPoints[[1]]@point.data)){
        ssn = fncUnloadWQ("WT",ssn) # unload if it's already loaded
      }
      ssn = fncLoadWQdata(type = "WT", dat.df = WT.df, thedate = dat.idx[dd], thetime = tt, ssn = ssn, plotit = FALSE)
      
      # Extract data frame that has had temperautre loaded
      wq.df = getSSNdata.frame(ssn,"preds")
      
      # changing factor/character fields to numeric, if necessary
      for(field in c("rid", so.field, wt.field)){
        if (is.factor(wq.df[,field])) {
          wq.df[,field] = as.numeric(levels(wq.df[,field])[as.numeric(wq.df[,field])])
        }
        if (is.character(wq.df[,field])) {
          wq.df[,field] = as.numeric(wq.df[,field])
        }
      }
      
      # Plot updated water temperature and fish locations
      # (only for this iteration due to storage constraints)
      if(plot.flag == TRUE & iter == plot.iter){
        
        png(paste0(mydir, "/", imageDir, "/", dat.idx[dd], "-", sprintf("%03d", tt), ".png"), width = 9, height = 6, units = "in", res = 150)
        par(mar = c(1, 0, 4, 0))
  
        # plot background
        if (outline) {
          plot(basin,col="gray80",border=NA)
        } else {
          plot(basin,col=grays,ext=ex,axes=FALSE, xlab = "", ylab = "", box = FALSE, legend = FALSE)
        }
        
        plot(streams, col="darkgray", add = TRUE)

        #plot water temperature for dd & tt as colored stream lines:
          datap = wq.df[,c("rid", wt.field)]
          datap2 = tapply(datap[,wt.field], datap$rid, mean, na.rm = TRUE)
          datap3 = cbind.data.frame("rid" = as.numeric(names(datap2)), "WT" = datap2)
          
          linedata = ssn@data
          rownames(linedata) = NULL
          linedata$sort.order = as.integer(rownames(linedata))
          if("WT" %in% colnames(linedata)) linedata = linedata[,-which(colnames(linedata) == "WT")] #remove WT field
          linedata2 = merge(linedata, datap3, by = "rid", all.x = TRUE) 
          linedata2 = linedata2[order(linedata2$sort.order),]
          for(n in 1:length(cb)) {
            linedata2$col.class[linedata2[,wt.field] >= left[n] & linedata2[,wt.field] <= rght[n]] = n
          }
          ssn@data = linedata2
          for (i in 1:length(ssn@lines)) {
            for (j in 1:length(ssn@lines[[i]])) {
              lines(ssn@lines[[i]]@Lines[[j]]@coords, col = cb[ssn@data[i,"col.class"]], lwd = 5*(ssn@data[i,"afvFlow"] + 0.4))
            }
          }
        
        # Plot fish locations for dd & tt
        open.salmon = salmon$pid[salmon$emrg == 0 & salmon$survive == 1]
        spawned.salmon = salmon$pid[salmon$emrg == -1 & salmon$survive == 1]
        emerged.salmon = salmon$pid[salmon$emrg == 1 & salmon$survive == 1]
        
        if(SecondSpecies == TRUE){
          alive.fish_other = fish_other$pid[fish_other$survive == 1]
          points(ssn@predpoints@SSNPoints[[fish_other.id]]@point.coords[,1][ssn@predpoints@SSNPoints[[fish_other.id]]@point.data$pid %in% alive.fish_other],
                 ssn@predpoints@SSNPoints[[fish_other.id]]@point.coords[,2][ssn@predpoints@SSNPoints[[fish_other.id]]@point.data$pid %in% alive.fish_other],
                 col = "gray50", pch = 15, cex = 0.5)
        }
        
        # salmon spawning sites
        points(ssn@predpoints@SSNPoints[[salmon.id]]@point.coords[,1][ssn@predpoints@SSNPoints[[salmon.id]]@point.data$pid %in% open.salmon],
              ssn@predpoints@SSNPoints[[salmon.id]]@point.coords[,2][ssn@predpoints@SSNPoints[[salmon.id]]@point.data$pid %in% open.salmon],
              col = 1, pch = 1, cex =0.5)
        
        # incubating salmon eggs
        points(ssn@predpoints@SSNPoints[[salmon.id]]@point.coords[,1][ssn@predpoints@SSNPoints[[salmon.id]]@point.data$pid %in% spawned.salmon],
              ssn@predpoints@SSNPoints[[salmon.id]]@point.coords[,2][ssn@predpoints@SSNPoints[[salmon.id]]@point.data$pid %in% spawned.salmon],
              col = 2, pch = 20, cex = 0.8)
        
        # emerged salmon
        points(ssn@predpoints@SSNPoints[[salmon.id]]@point.coords[,1][ssn@predpoints@SSNPoints[[salmon.id]]@point.data$pid %in% emerged.salmon],
              ssn@predpoints@SSNPoints[[salmon.id]]@point.coords[,2][ssn@predpoints@SSNPoints[[salmon.id]]@point.data$pid %in% emerged.salmon],
              col = 4, pch = 20, cex = 0.8)
        
        # all points: points(ssn@predpoints@SSNPoints[[id]]@point.coords[,1],ssn@predpoints@SSNPoints[[id]]@point.coords[,2],col=1,pch=20,cex=0.8)
    
        # track a few individuals:
        if(track.individuals == TRUE){
          fish.to.track = 1:5 #sample(salmon$pid, 5)
          i = 10
          for(f in fish.to.track){
            points(ssn@predpoints@SSNPoints[[salmon.id]]@point.coords[,1][ssn@predpoints@SSNPoints[[salmon.id]]@point.data$pid == f],
                   ssn@predpoints@SSNPoints[[salmon.id]]@point.coords[,2][ssn@predpoints@SSNPoints[[salmon.id]]@point.data$pid == f],
                   col = i, pch = 20, cex = 0.8)
            i= i + 1
          }
          rm(i)
        }
        
        # Add legend
        leglabs = paste(left, "to", rght)
        legend("right", legend = leglabs, title=expression("Temperature "~degree(C)),bty = "n", pch = 19, col = cb, cex = 0.8)
        
        # Add scale bar
        rect(xleft = ex[1] + 5000, ybottom = ex[3] + 3000, xright = ex[1] + 7500, ytop = ex[3] + 3500)
        rect(xleft = ex[1] + 7500, ybottom = ex[3] + 3000, xright = ex[1] + 10000, ytop= ex[3] + 3500, col = 1)
        rect(xleft = ex[1] + 10000, ybottom = ex[3] + 3000, xright = ex[1] + 12500, ytop = ex[3] + 3500)
        rect(xleft = ex[1] + 12500, ybottom = ex[3] + 3000, xright = ex[1] + 15000, ytop = ex[3] + 3500, col = 1)
        segments(x0 = ex[1] + 5000, y0 = ex[3] + 3000, x1 = ex[1] + 5000, y1 = ex[3] + 2500)
        segments(x0 = ex[1] + 10000, y0=ex[3] + 3000, x1 = ex[1] + 10000, y1 = ex[3] + 2500)
        segments(x0 = ex[1] + 15000, y0 = ex[3] + 3000, x1 = ex[1] + 15000, y1 = ex[3] + 2500)
        text(x = ex[1] + 5000, y = ex[3] + 1500, "0", cex = 0.8)
        text(x = ex[1] + 10000, y = ex[3] + 1500, "5", cex = 0.8)
        text(x = ex[1] + 15000, y = ex[3] + 1500, "10", cex = 0.8)
        text(x = ex[1] + 18000, y = ex[3] + 1500, "km", cex = 0.8) 
        
        # Add north arrow
        arrows(ex[1] + 2000, ex[3] + 2300, ex[1] + 2000, ex[3] + 4000, length = 0.1, lwd = 5)
        text(ex[1] + 2000, ex[3] + 1000, "N")
        
        # Add title
        mtext(thetitle, side = 3, line=2, cex = 1.3)
    
        dev.off()
      } # end plotting for the selected iter

      
  # 2. SPAWN SALMON
      # if there are salmon left to spawn,...
      if (any(salmon$emrg == 0)){
        # if today is on or past the day spawning begins,...
        if(dat.idx[dd] >= spawn.init){
          # if spawning occurs over multiple days,...
          if(spawn.date.variable){
            # if there are many salmon left to spawn,...
            if (sum(salmon$emrg == 0) > 1) {
              # sample x = indices of salmon that haven't spawned yet
              # sample n = number of salmon at species.spawn.times[spawn.counter] or number of salmon left unemerged, whichever is smaller
              set.seed(iter)
              spawn.index = sample(x = salmon$pid[salmon$emrg == 0], size = min(sum(salmon$emrg == 0), salmon.spawn.times[spawn.counter]))
            } else { # else there is one salmon left to spawn,...
              # if one salmon should spawn today, spawn the last salmon.
              if (salmon.spawn.times[spawn.counter] == 1) spawn.index = which(salmon$emrg == 0)
              # (can't use sample because function behavior is different when x has length of 1)
            }
            # spawn salmon and record date
            salmon$emrg[salmon$pid %in% spawn.index] = -1
            salmon$dateSp[salmon$dateSp == day1 & salmon$emrg == -1] = strptime(dat.idx[dd], format = "%Y-%m-%d")
            # increment spawning counter
            spawn.counter = spawn.counter + 1 
          } else { # else spawning is fixed,...
            # all salmon spawn now
            salmon$emrg = -1
            salmon$dateSp = strptime(dat.idx[dd], format = "%Y-%m-%d")
          } # end checking variable spawning
        } # end checking spawn start date
      } # end checking for unspawned salmon

  # 3. INCUBATE salmon & accumulate thermal exposure
      
      # Assign closest water temperature to salmon
      salmon$WT[salmon$survive == 1] = fncGetNearestAttribute(wq.df[,c("pid", "rid", "ratio", wt.field)], salmon[salmon$survive == 1, c("pid", "seg", "ratio")], ssn)[,wt.field]
      
      # Accumulate thermal units (TUs) While salmon are spawned or emerged
      salmon$TU[salmon$survive == 1 & salmon$emrg != 0] = salmon$TU[salmon$survive == 1  & salmon$emrg != 0] + salmon$WT[salmon$survive == 1  & salmon$emrg != 0] * 0.5 #halved to account for half-day time step
      
      if(SecondSpecies == TRUE){
        fish_other$WT[fish_other$survive == 1] = fncGetNearestAttribute(wq.df[,c("pid", "rid", "ratio", wt.field)], fish_other[fish_other$survive == 1, c("pid", "seg", "ratio")], ssn)[,wt.field]
        fish_other$TU[fish_other$survive == 1] = fish_other$TU[fish_other$survive == 1] + fish_other$WT[fish_other$survive == 1] * 0.5 #halved to account for half-day time step
      }
        
    
  # 4. EMERGE salmon
      # Alevins emerge after reaching an accumulated thermal unit (TU) threshold
      if (any(salmon$emrg == -1)) {
        salmon$emrg[salmon$emrg == -1 & salmon$TU >= parameters["salmon", "ATU.crit"]] = 1
        salmon$dateEm[salmon$dateEm == day1 & salmon$emrg == 1] = strptime(dat.idx[dd], format="%Y-%m-%d")
      }
      
      
  # 5. MOVE FISH
      # After emergence, assess probability of movement and move fish if appropriate

      # 5a. Move salmon individually

      # index for which salmon have emerged and can grow, and are still alive
      salmon.alive.emerge.index = salmon$pid[salmon$emrg == 1 & salmon$survive == 1]
      
      if(length(salmon.alive.emerge.index) > 0){
        
        # Set which growth lookup table to use
        sp.idx = 1
        wt.growth = wt.growth.salmon
        fish = salmon[salmon.alive.emerge.index,]
        
        # Determine a fish's drive to move downstream
        set.seed(iter)
        om.prob = fncDownstreamDrive(w = salmon[salmon.alive.emerge.index, "weight"], om.mass = parameters["salmon", "om.mass"], om.date.taper = parameters["salmon", "om.date.taper"])
        fish$direction = rbinom(n = nrow(fish), size = 1, prob = om.prob)
        fish$om.prob = om.prob

        # Run movement decision functions
          #1. Get movement distance for fish that can move during this time step
          moveDist = fncMoveDistance(fish, sp.idx, ssn)
          
          #2. Move fish individually (one fish, one time step)
          moved = lapply(fish[,"pid"], function(x) fncMoveIndividual(fpid = x, fish, sp.idx, ssn))
          fish = do.call(rbind, moved)
          
          #3. Update x and y coordinates
          # if this was a generated network, use straight-line distance calcs
          # if this network was created in GIS, use different approach for moving along a potentially curvy reach
          ifelse(length(grep("network-", netnm)) > 0, 
                 locs<- ddply(fish, .(pid), function(x) data.frame(fncGetXY(x$seg, x$ratio, ssn))), 
                 locs<- ddply(fish, .(pid), function(x) data.frame(fncGetXY.arc(x$seg, x$ratio, ssn))))
          fish$xloc<- locs$xloc
          fish$yloc<- locs$yloc
          
          #4. Update data in SSN
          # data table
          ssn@predpoints@SSNPoints[[sp.idx + 1]]@point.data$upDist[salmon.alive.emerge.index] = fish$upDist
          ssn@predpoints@SSNPoints[[sp.idx + 1]]@point.data$ratio[salmon.alive.emerge.index] = fish$ratio
          ssn@predpoints@SSNPoints[[sp.idx + 1]]@point.data[,xlat][salmon.alive.emerge.index] = fish$xloc
          ssn@predpoints@SSNPoints[[sp.idx + 1]]@point.data[,ylon][salmon.alive.emerge.index] = fish$yloc
          ssn@predpoints@SSNPoints[[sp.idx + 1]]@point.data$rid[salmon.alive.emerge.index] = fish$seg
          # note: sp.idx + 1 because 'preds' is in the 1st position, so the first fish is in the second position, etc.
          
          # coordinates
          ssn@predpoints@SSNPoints[[sp.idx + 1]]@point.coords[,1][salmon.alive.emerge.index] = fish$xloc 
          ssn@predpoints@SSNPoints[[sp.idx + 1]]@point.coords[,2][salmon.alive.emerge.index] = fish$yloc
          
          salmon[salmon.alive.emerge.index,] = fish
          rm(fish)
          
        # Filter out fish that have outmigrated as subyearlings
          #unless the date is past when salmon have been observed in smolt trap, then turn them around
          if(dat.idx[dd] < as.Date(paste0(climate.year,"-", parameters["salmon","om.date.end"]))){
            salmon$survive[salmon$seg %in% network.base.segs] = 2
          } else{
            salmon$direction[salmon$seg == network.base.segs[1] & salmon$direction == 1] = 0
          }
          
        # Reset salmon.alive.emerge.index
          salmon.alive.emerge.index = salmon$pid[salmon$survive == 1 & salmon$emrg == 1]
      }
      
      # 5b. Move fish_other individually
      
      if(SecondSpecies == TRUE){
        
        # index for which fish_other are alive
        fish_other.alive.index = fish_other$pid[fish_other$survive == 1]
        
        # Jitter initial position
        if(dd == 1 & tt == 6){
          fish_other[fish_other.alive.index,] = fncJitterPosition(wt.growth = wt.growth.fish_other, 
                    fae = fish_other.alive.index, fish = fish_other[fish_other.alive.index,], sp.idx = 2)
        }
        
        if(length(fish_other.alive.index) > 0){
        
        # set which growth lookup table to use based on species 'spp'
        sp.idx = 2
        wt.growth = wt.growth.fish_other
        fish = fish_other[fish_other.alive.index,]
        
        # Run movement decision functions
        # Get movement distance for fish that can move during this time step
        moveDist = fncMoveDistance(fish, sp.idx, ssn)
        
        # Move fish individually (one fish, one time step)
        moved = lapply(fish[,"pid"], function(x) fncMoveIndividual(fpid = x, fish, sp.idx, ssn))
        fish = do.call(rbind, moved)
        
        # Update x and y coordinates
        # if this was a generated network, use straight-line distance calcs
        # if this network was created in GIS, use different approach for moving along a potentially curvy reach
        ifelse(length(grep("network-", netnm)) > 0, 
               locs<- ddply(fish, .(pid), function(x) data.frame(fncGetXY(x$seg, x$ratio, ssn))), 
               locs<- ddply(fish, .(pid), function(x) data.frame(fncGetXY.arc(x$seg, x$ratio, ssn))))
        fish$xloc<- locs$xloc
        fish$yloc<- locs$yloc
        
        # Update data in SSN
        # data table
        ssn@predpoints@SSNPoints[[sp.idx + 1]]@point.data$upDist[fish_other.alive.index] = fish$upDist
        ssn@predpoints@SSNPoints[[sp.idx + 1]]@point.data$ratio[fish_other.alive.index] = fish$ratio
        ssn@predpoints@SSNPoints[[sp.idx + 1]]@point.data[,xlat][fish_other.alive.index] = fish$xloc
        ssn@predpoints@SSNPoints[[sp.idx + 1]]@point.data[,ylon][fish_other.alive.index] = fish$yloc
        ssn@predpoints@SSNPoints[[sp.idx + 1]]@point.data$rid[fish_other.alive.index] = fish$seg
        # note: sp.idx + 1 because 'preds' is in the 1st position, so the first fish is in the second position, etc.
        
        # coordinates
        ssn@predpoints@SSNPoints[[sp.idx + 1]]@point.coords[,1][fish_other.alive.index] = fish$xloc 
        ssn@predpoints@SSNPoints[[sp.idx + 1]]@point.coords[,2][fish_other.alive.index] = fish$yloc
        
        fish_other[fish_other.alive.index,] = fish
        rm(fish)
        
        }
       
       # Turn fish around if they've reached the base of the network 
        fish_other$direction[fish_other$seg == network.base.segs[1] & fish_other$direction == 1] = 0
        
       # Reset fish_other.alive.index
       fish_other.alive.index = fish_other$pid[fish_other$survive == 1]
      }
      

  # 6. RE-ASSESS HABITAT QUALITY (to incorporate changes after movement into growth calcs)
        
      # 6a. Get nearest water temperature
      # Update these fields in fish dataframes
      if(length(salmon.alive.emerge.index) > 0){
      salmon$WT[salmon.alive.emerge.index] = fncGetNearestAttribute(wq.df[,c("pid", "rid", "ratio", wt.field)], salmon[salmon.alive.emerge.index,c("pid", "seg", "ratio")], ssn)[,wt.field]
      }
      
      if(SecondSpecies == TRUE){
        if(length(fish_other.alive.index) > 0){
          fish_other$WT[fish_other.alive.index] = fncGetNearestAttribute(wq.df[,c("pid", "rid", "ratio", wt.field)], fish_other[fish_other.alive.index,c("pid","seg","ratio")], ssn)[,wt.field] 
        }
      }
        
      # 6b. Update fish densities to incorporate fish that just moved into density estimates
      if(length(salmon.alive.emerge.index) > 0){
        if(SecondSpecies == FALSE) {fish_other = NULL; fish_other.alive.index = NULL}
        density = as.data.frame(fncFishDensity(salmon, fish_other, salmon.alive.emerge.index, fish_other.alive.index, ssn = ssn))
        namedvec = density$conspecific.bio.density; names(namedvec) = paste0("x",density$seg)
        salmon$conspecificdensity = fncGetValue(mykey = paste0("x",salmon$seg), mylookupvector = namedvec)
        #reset eggs to zero density
        salmon$conspecificdensity[salmon$emrg!=1] = 0
        namedvec = density$other.bio.density; names(namedvec) = paste0("x",density$seg)
        salmon$otherdensity = fncGetValue(mykey = paste0("x",salmon$seg), mylookupvector = namedvec)
        salmon$totaldensity = salmon$conspecificdensity + salmon$otherdensity
        rm(density, namedvec)
      }
      
      if(SecondSpecies == TRUE){
        if(length(fish_other.alive.index) > 0){
          density = as.data.frame(fncFishDensity(fish_other, salmon, fish_other.alive.index, salmon.alive.emerge.index))
          namedvec = density$conspecific.bio.density; names(namedvec) = paste0("x",density$seg)
          fish_other$conspecificdensity = fncGetValue(mykey = paste0("x",fish_other$seg), mylookupvector = namedvec)
          namedvec = density$other.bio.density; names(namedvec) = paste0("x",density$seg)
          fish_other$otherdensity = fncGetValue(mykey = paste0("x",fish_other$seg), mylookupvector = namedvec)
          fish_other$totaldensity = fish_other$conspecificdensity + fish_other$otherdensity
          rm(density, namedvec)
        }
      }
        
      # 6c. Update available ration after fish have moved and after accounting for density
        # (ration linearly decreases with fish density at the new location)
      
      # For reaches with salmon in them:
      if(length(salmon.alive.emerge.index) > 0){
        
        # get salmon density
        fdens <- salmon$conspecificdensity
        
        # force high densities to cap out at this parameter value
        fdens[fdens > parameters["salmon", "max.density"]] = parameters["salmon", "max.density"]
        
        # calculate density effect
        density.effect = fncRescale((1 - c(fdens, 0.01, parameters["salmon", "max.density"])), c((0.5), 1))
        density.effect = density.effect[-c(length(density.effect), (length(density.effect) -1 ))] #remove the last 2 values that were used to standardize the range
        density.effect = cbind("rid" = salmon$seg, "pid" = salmon$pid, density.effect)
  
        # Update ration in SSN & wq.df based on fish density from above
        namedvec = wq.df$ration_base; names(namedvec) = paste0("x", wq.df$rid)
        ration1 = fncGetValue(paste0("x", wq.df$rid), namedvec)
        namedvec = tapply(density.effect[,"density.effect"], density.effect[,"rid"], mean, na.rm = TRUE); names(namedvec) = paste0("x",names(namedvec))
        dns.eff1 = fncGetValue(paste0("x", wq.df$rid), namedvec)
        updated.ration = ration1 * dns.eff1
        #updated.ration = updated.ration[1:nwq]
        ssn@predpoints@SSNPoints[[1]]@point.data$ration[!is.na(updated.ration)] = 
          wq.df$ration[!is.na(updated.ration)] = updated.ration[!is.na(updated.ration)]
      }
      
      # For reaches with fish_other in them:
      if(SecondSpecies == TRUE){
        if(length(fish_other.alive.index) > 0){
          
          # get fish_other density
          fdens <- fish_other$conspecificdensity
          
          # force high densities to cap out at this parameter value
          fdens[fdens > parameters["fish_other", "max.density"]] = parameters["fish_other", "max.density"]
          
          # calculate density effect
          density.effect = fncRescale((1 - c(fdens, 0.01, parameters["fish_other", "max.density"])), c((0.5), 1))
          density.effect = density.effect[-c(length(density.effect), (length(density.effect) -1 ))] #remove the last 2 values that were used to standardize the range
          density.effect = cbind("rid" = fish_other$seg, "pid" = fish_other$pid, density.effect)
          
          # Update ration in SSN & wq.df based on fish density from above
          namedvec = wq.df$ration_base_ss; names(namedvec) = paste0("x", wq.df$rid)
          ration1 = fncGetValue(paste0("x", wq.df$rid), namedvec)
          namedvec = tapply(density.effect[,"density.effect"], density.effect[,"rid"], mean, na.rm = TRUE); names(namedvec) = paste0("x",names(namedvec))
          dns.eff1 = fncGetValue(paste0("x", wq.df$rid), namedvec)
          updated.ration = ration1 * dns.eff1
          #updated.ration = updated.ration[1:nwq]
          ssn@predpoints@SSNPoints[[1]]@point.data$ration_ss[!is.na(updated.ration)] = 
            wq.df$ration_ss[!is.na(updated.ration)] = updated.ration[!is.na(updated.ration)]
        }
      }
      
      # Update ration in fish dataframes
      if(length(salmon.alive.emerge.index) > 0){
        salmon$ration[salmon.alive.emerge.index] = fncGetNearestAttribute(wq.df[,c("pid", "rid", "ratio", "ration")], salmon[salmon.alive.emerge.index,c("pid","seg","ratio")], ssn)[,"ration"] 
      }
      
      if(SecondSpecies == TRUE){
        if(length(fish_other.alive.index) > 0){
          fish_other$ration[fish_other.alive.index] = fncGetNearestAttribute(wq.df[,c("pid", "rid", "ratio", "ration_ss")], fish_other[fish_other.alive.index,c("pid","seg","ratio")], ssn)[,"ration_ss"] 
        }
      }
         
  # 7. GROW FISH. (calculated using the Wisconsin Bioenergetics model equations)

      # salmon:
      if(length(salmon.alive.emerge.index) > 0){
        # Get parameters and inputs for Bioenergetics model
        # provide nearest water temperature, fish weight (grams), and pvals or ration
        salmon.input= fncGetBioEParms(parameters["salmon","spp"], 
              parameters["salmon","pred.en.dens"], parameters["salmon","prey.en.dens"],
              wt.nearest = salmon[salmon.alive.emerge.index,c("pid",wt.field)], 
              startweights = salmon$weight[salmon.alive.emerge.index],
              pvals = rep(-9, nrow(salmon[salmon.alive.emerge.index,])), 
              ration = salmon$ration[salmon.alive.emerge.index])
        
        # Run bioenergetics for all fish at this 12-hour time step
        salmon.results = BioE(salmon.input, salmon.constants)
        
        # Put results back into fish dataframe
        # note: splitting results in half because bioenergetics model operates on 
        # a daily time step and we are using a half-day time step
        
        # instantaneous growth (units are g/d; halved and divided by fish weight to get g/g/12h)
        salmon$growth[salmon.alive.emerge.index]= t(salmon.results$Growth / 2) / salmon$weight[salmon.alive.emerge.index]
        # cumulative growth
        salmon$weight[salmon.alive.emerge.index] = salmon$weight[salmon.alive.emerge.index] + t(salmon.results$Growth / 2)
        # instantaneous consumption (units are g/g/d; halved to get g/g/12h)
        salmon$consInst[salmon.alive.emerge.index] = t(salmon.results$Consumption / 2)
        # cumulative amount consumed
        salmon$consCum[salmon.alive.emerge.index] = salmon$consCum[salmon.alive.emerge.index] + t(salmon.results$Consumption / 2)
        # bioenergetic p-values
        pp = salmon.input$ration/salmon.results$CMAX; pp[pp > 1] = 1
        salmon$pvals[salmon.alive.emerge.index] = pp
        
        # mark salmon with zero or negative weights as dead and reset index
        salmon$survive[salmon$weight <= 0] = 0
        salmon.alive.emerge.index = salmon$pid[salmon$emrg == 1 & salmon$survive == 1]
      }
      
      # fish_other:
      if(SecondSpecies == TRUE){
        if(length(fish_other.alive.index) > 0){
          # Get parameters and inputs for Bioenergetics model
          # provide nearest water temperature, fish weight (grams), and pvals or ration
          fish_other.input= fncGetBioEParms(parameters["fish_other","spp"], 
                parameters["fish_other","pred.en.dens"], parameters["fish_other","prey.en.dens"],
                wt.nearest = fish_other[fish_other.alive.index,c("pid",wt.field)], 
                startweights = fish_other$weight[fish_other.alive.index],
                pvals = rep(-9, nrow(fish_other[fish_other.alive.index,])), 
                ration = fish_other$ration[fish_other.alive.index])
        
          # Run bioenergetics for all fish at this 12-hour time step
          fish_other.results = BioE(fish_other.input, fish_other.constants)
  
          # Put results back into fish dataframe
          # note: splitting results in half because bioenergetics model operates on 
          # a daily time step and we are using a half-day time step
          
          # instantaneous growth (units are g/d; halved and divided by fish weight to get g/g/12h)
          fish_other$growth[fish_other.alive.index] = t(fish_other.results$Growth / 2) / fish_other$weight[fish_other.alive.index]
          # cumulative growth
          fish_other$weight[fish_other.alive.index] = fish_other$weight[fish_other.alive.index] + t(fish_other.results$Growth / 2)
          # instantaneous consumption (units are g/g/d; halved to get g/g/12h)
          fish_other$consInst[fish_other.alive.index] = t(fish_other.results$Consumption / 2)
          # cumulative amount consumed
          fish_other$consCum[fish_other.alive.index] = fish_other$consCum[fish_other.alive.index] + t(fish_other.results$Consumption / 2)
          # bioenergetic p-values
          pp = fish_other.input$ration/fish_other.results$CMAX; pp[pp > 1] = 1
          fish_other$pvals[fish_other.alive.index] = pp
          
          # mark fish with zero or negative weights as dead and reset index
          fish_other$survive[fish_other$weight <= 0] = 0
          fish_other.alive.index = fish_other$pid[fish_other$survive == 1]
        }
      }
      
      # Record the date fish smolt (oumigrate as subyearling)
      salmon$dateOm[salmon$dateOm == day1 & salmon$survive == 2] = strptime(dat.idx[dd], format = "%Y-%m-%d")

      
  # 8. PREDATION
      # Bass are predators if they are alive and above a critical temperature.
      # Salmon are potentialy prey if they are alive, emerged, in the same segment 
      # as a predator, within a predator's movement distance, and they're smaller
      # than a predator by a critical mass ratio.
      # Probability of predation increases linearly with temperature.
      if (predation.flag == TRUE) {
        for (pred.index in 1:nrow(fish_other)) {
          pred.output <- fncOnePredator(pred.survive = fish_other$survive[pred.index],
                                        pred.temp   = fish_other$WT[pred.index],
                                        pred.seg    = fish_other$seg[pred.index],
                                        pred.weight = fish_other$weight[pred.index],
                                        pred.dist   = fish_other$length2segBase[pred.index],
                                        pred.move   = parameters["fish_other", "pred.move"],
                                        pred.temp.crit = parameters["fish_other", "pred.temp.crit"],
                                        pred.mass.crit = parameters["fish_other", "pred.mass.crit"],
                                        max.pred.prob = parameters["fish_other", "max.pred.prob"],
                                        prey.id = salmon$pid,
                                        prey.survive = salmon$survive,
                                        prey.emerge = salmon$emrg,
                                        prey.seg    = salmon$seg,
                                        prey.weight = salmon$weight,
                                        prey.dist   = salmon$length2segBase)
          fish_other$pred.prob[pred.index] <- pred.output$predation.probability
          fish_other$num.prey[pred.index] <- pred.output$number.potential.prey
          fish_other$num.eaten[pred.index] <- pred.output$number.prey.eaten
          salmon$survive[which(salmon$pid %in% pred.output$prey.eaten.id)] <- -1
        }
        rm(pred.output)
        
        # Update living emerged fish indices
        salmon.alive.emerge.index = salmon$pid[salmon$emrg == 1 & salmon$survive == 1]

        # Record the date fish were eaten by predator
        salmon$datePr[salmon$datePr == day1 & salmon$survive == -1] = strptime(dat.idx[dd], format = "%Y-%m-%d")
      }
        
      # 9. OTHER MORTALITY
      set.seed(iter)
      salmon$survive[salmon.alive.emerge.index] = fncSurvive(df = salmon[salmon.alive.emerge.index, c("weight","growth")], minprob = parameters["salmon","survival.min"], maxprob = parameters["salmon","survival.max"], b = parameters["salmon","survival.shape"])
      salmon.alive.emerge.index = salmon$pid[salmon$emrg == 1 & salmon$survive == 1]
      # Record the date fish died
      salmon$dateDi[salmon$dateDi == day1 & salmon$survive == 0] = strptime(dat.idx[dd], format = "%Y-%m-%d")
        
      if(SecondSpecies == TRUE){
        if(length(fish_other.alive.index) > 0){
          fish_other$survive[fish_other.alive.index] = fncSurvive(df = fish_other[fish_other.alive.index, c("weight","growth")], minprob = parameters["fish_other","survival.min"], maxprob = parameters["fish_other","survival.max"], b = parameters["fish_other","survival.shape"])
          fish_other.alive.index = fish_other$pid[fish_other$survive ==1]
        }
        # Record the date fish died
        fish_other$dateDi[fish_other$dateDi == day1 & fish_other$survive == 0] = strptime(dat.idx[dd], format = "%Y-%m-%d")
      }

      # Remove WT 
      ssn = fncUnloadWQ("WT",ssn) 
    
      # Force dates to numeric for storage in arrays
      # Can be reformatted as date again later as follows:
      # as.POSIXct(salmon$dateSp, origin = "1970-01-01")
      salmon.numeric = salmon
      salmon.numeric[,"dateSp"][salmon.numeric[,"dateSp"] == day1] = NA #unemerged fish
      salmon.numeric[,"dateEm"][salmon.numeric[,"dateEm"] == day1] = NA #unemerged fish
      salmon.numeric[,"dateOm"][salmon.numeric[,"dateOm"] == day1] = NA #unemerged fish
      salmon.numeric[,"datePr"][salmon.numeric[,"datePr"] == day1] = NA #unemerged fish
      salmon.numeric[,"dateDi"][salmon.numeric[,"dateDi"] == day1] = NA #unemerged fish
      salmon.numeric[,"dateEm"] = as.numeric(salmon.numeric[,"dateEm"])
      salmon.numeric[,"dateSp"] = as.numeric(salmon.numeric[,"dateSp"])
      salmon.numeric[,"dateOm"] = as.numeric(salmon.numeric[,"dateOm"])
      salmon.numeric[,"datePr"] = as.numeric(salmon.numeric[,"datePr"])
      salmon.numeric[,"dateDi"] = as.numeric(salmon.numeric[,"dateDi"])
      
      if(SecondSpecies == TRUE){
        fish_other.numeric = fish_other
        fish_other.numeric[,"dateSp"] = NA
        fish_other.numeric[,"dateEm"] = NA
        fish_other.numeric[,"dateOm"] = NA
        fish_other.numeric[,"datePr"] = NA
        fish_other.numeric[,"dateDi"] = as.numeric(fish_other.numeric[,"dateDi"])
      }
      
      # Store multidimensional array of results for this date/time
      salmon.array[,,dt] = as.matrix(salmon.numeric)[,array.cols2keep]; colnames(salmon.array) = array.cols2keep
      if(SecondSpecies == TRUE) {fish_other.array[,,dt] = as.matrix(fish_other.numeric)[,array.cols2keep]; colnames(fish_other.array) = array.cols2keep}
  
      dt = dt + 1
    } #end time steps (tt)
  } #end dates (dd)

  # Mark any unemerged fish as dead
  salmon.array[,"survive",dim(salmon.array)[3]][salmon.array[,"emrg",dim(salmon.array)[3]] != 1] = 0
  
  # Store final step output
  salmon.finalstep = as.matrix(salmon.numeric)[,output.cols2keep]; colnames(salmon.finalstep) = output.cols2keep
  if(SecondSpecies == TRUE){fish_other.finalstep = as.matrix(fish_other.numeric)[,output.cols2keep]; colnames(fish_other.finalstep) = output.cols2keep}


# Save data for the iteration
  save("salmon.array", file = paste0(outputDir,"/salmon.array.",climate.year,".",iter,".RData"))
  save("salmon.finalstep",file = paste0(outputDir,"/salmon.finalstep.",climate.year,".",iter,".RData"))
  if(SecondSpecies == TRUE){
    save("fish_other.array", file = paste0(outputDir,"/fish_other.array.",climate.year,".",iter,".RData"))
    save("fish_other.finalstep",file = paste0(outputDir,"/fish_other.finalstep.",climate.year,".",iter,".RData"))
  }

  end.time = proc.time()
  runtime = (end.time[3] - start.time[3]) / 60 / 60
  
  # Store info on no. replicates, climate scenarios, and runtime
  textDir = paste0(outputDir, "/run.info.", climate.year, ".", iter, ".txt")
  file.create(textDir)
  
  run.info<- c(paste0("netnwork: ",netnm),
          paste0("year: ", climate.year), 
          paste0("runtime (h): ", runtime),
          paste0("seed: ", iter),
          paste0("run: ", run),
          paste0("predation.flag: ", predation.flag),
          paste0("spawn.date.variable: ",spawn.date.variable),
          paste0("plot.flag: ",plot.flag),
          paste0("first date:",first.date),
          paste0("last.date:",last.date))
  salmon.parms = NULL; for(i in 1:ncol(parameters)) {var = colnames(parameters)[i]; salmon.parms[i] = paste0(var,": ",parameters["salmon",var]); fish_other.parms = NULL}
  if(SecondSpecies == TRUE) for(i in 1:ncol(parameters)) {var = colnames(parameters)[i]; fish_other.parms[i] = paste0(var,": ",parameters["fish_other",var])}
  
  write(c("RUN INFO:", run.info," ","SALMON PARAMETERS:",salmon.parms," ","FISH_OTHER PARAMETERS:",fish_other.parms),file = textDir)

  # Produce quick summary
  source("code/summary.R")
  
} #end iteration

#=== END OF FILE ===============================================================

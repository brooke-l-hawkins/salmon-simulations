#*******************************************************************************
#
# Functions for SnoIBM
#   An individual-based model for evaluating 
#   how phenology, growth and survival of juvenile Chinook salmon may respond 
#   to altered thermal regimes and predation by largemouth bass in the
#   Snoqualmie River watershed, WA, USA.
#
# Coded by:
#    B.L. Hawkins (general, second species, predation), A.H. Fullerton (general), 
#    B.J. Burke (movement), N. Som (network shapes) & M. Nahorniak (bioenergetics) 
#    Last Updated 17 Apr 2020
#
#*******************************************************************************


#=== FUNCTIONS TO LOAD DATA INTO SSN ========================================================

#Load temperature data for a particular date & time
fncLoadWQdata<- function(type,  dat.df, thedate, thetime, ssn = parent.frame()$ssn, plotit = F){
  
  idnm = "pid"
  obsdat = getSSNdata.frame(ssn, "preds")

  arcID = as.numeric(gsub("X", "", colnames(dat.df)[-c(1:2)]))
  dat = t(dat.df[dat.df[,"Date"] == as.Date(thedate) & dat.df[,"Time"] == thetime, 3:ncol(dat.df)])
  dat = cbind(arcID, dat)
  row.names(dat) = NULL
  if(type == "WT") colnames(dat) = c(idnm, wt.field)
  
  #dat can now be merged to wq.df and has the temperature values for each reach for that date & time
  obsdat = merge(obsdat, dat, by.x = idnm, by.y= idnm, all.x = TRUE) 
  obsdat = obsdat[order(obsdat$pid),]; row.names(obsdat) = obsdat$pid #super important step: needs to be ordered by pid and row.names = pid
  
  #put data back in SSN
  ssn = putSSNdata.frame(obsdat, ssn, "preds")
  
  if(plotit == TRUE) plot(ssn, PredPointsID = "preds", VariableName = type)
  
  return(ssn)
}

#Unload water temperature data from SSN
fncUnloadWQ<- function(type, ssn = parent.frame()$ssn){
  
  dat = getSSNdata.frame(ssn, "preds")
  
  if(type == "WT"){
    if(length(which(colnames(dat) == "WT")) > 0) dat = dat[,-which(colnames(dat) == "WT")]
    if(length(which(colnames(ssn@data) == "WT")) > 0) ssn@data = ssn@data[,-which(colnames(ssn@data) == "WT")]
    if(length(which(colnames(ssn@data) == "WT.color")) > 0) ssn@data = ssn@data[,-which(colnames(ssn@data) == "WT.color")]
  }
  ssn = putSSNdata.frame(dat, ssn, "preds")
  
  return(ssn)
}

#Subsample fish shapefile
fncFishShp<- function(species, nSpawners, nStart, ssn = parent.frame()$ssn, plotit = FALSE){
  
  # Load shapefile with largest number of possible spawning locations
  fish.id = which(parameters[,"species"] == species)[1] 
  shpnm = substr(parameters[fish.id, "fish.shp"], 1, (nchar(parameters[fish.id, "fish.shp"])-4))
  shpnm2 = paste0(shpnm, nStart)
  shp = readOGR(paste0(mydir,"/",loadDir,"/", ssn.folder), shpnm2)

  # Subsample how many fish to keep ('nSpawners')
  set.seed(1)
  fish2keep = sample(shp@data$pid, nSpawners)
  fish_spawn = subset(shp, shp@data$pid %in% fish2keep)
  row.names(fish_spawn@data) = NULL

  # Save shapefile for later use
  writeOGR(fish_spawn, paste0(mydir,"/",loadDir,"/", ssn.folder), substr(shpnm, 1, nchar(shpnm) - 1), driver="ESRI Shapefile", overwrite_layer = T)

  if(plotit == TRUE){
    plot(streams, lwd = streams@data$meanmsq / 100000000, col = "darkgray")
    plot(shp, add = 2, col = 2, cex = 0.5)
    plot(shp[shp@data$pid %in% fish2keep,], add = 2, col = 4, cex = 0.5)
  }
  #\\ Read back in subsetted shapefile [this how a file will be loaded for modeling]
  # ssn = importSSN("data.in/sno.ssn", predpts = 'preds')
  # ssn = importPredpts(ssn, paste0(shpnm, nSpawners), "ssn")
  # salmon.id = 2 #this is the 2nd preds file loaded in the SSN (see ssn@predpoints@ID)
  # ssn@predpoints@ID[salmon.id]="chinook"
  # td = getSSNdata.frame(ssn,"chinook")
  
}


#=== BASIC FUNCTIONS ===========================================================

# Look up value based on named vector
fncGetValue <- function(mykey, mylookupvector){
  myvalue <- mylookupvector[mykey]
  myvalue <- unname(myvalue)
  return(myvalue)
}

# Get the four coordinates for a single segment
fncGetSegCoords<- function(seg = seg, ssn = parent.frame()$ssn) {
  idx = which(ssn@data$rid == seg)
  return(ssn@lines[[idx]]@Lines[[1]]@coords)
}

# Distance to network base for a segment (including the length of this segment)
fncDist2base4seg<-  function(seg = seg, ssn = parent.frame()$ssn) {
  return(ssn@data$upDist[ssn@data$rid == seg])
}

# Distance to network base for a fish
fncDist2base4fish<- function(fishID = pid, fsh = NULL, id = salmon.id, ssn = parent.frame()$ssn){
  if(is.null(fsh)) fsh = ssn@predpoints@SSNPoints[[id]]@point.data #get fish data
  seg = fsh$seg[fsh$pid == fishID] #which seg is point sitting on
  seg.updist = ssn@data$upDist[ssn@data$rid == seg] #distance to top of segment from base of network
  seg.length = ssn@data[,length.field][ssn@data$rid == seg] #length of seg
  fish.lng2base = ssn@data[ssn@data$rid == seg, length.field] * fsh$ratio[fsh$pid == fishID] #length from point to base of its seg
  ans = seg.updist -seg.length + fish.lng2base #distance to point from base of network
  return(ans)
  # distance from base to top of seg - length of seg + distance from base of seg up to point
}

# Convert segment ratio (proprtional length) to length downstream
fncSegRatio2length<- function(seg = seg, segRatio = 0.5, ssn = parent.frame()$ssn) {
  return (ssn@data[ssn@data[,"rid"] == seg,length.field] * segRatio)
}

# Convert length downstream to segment ratio
fncLength2segRatio<- function(seg = seg, lengthDownstream = 0.5, ssn = parent.frame()$ssn) {
  return (lengthDownstream/ssn@data[ssn@data[,"rid"] == seg,length.field])
}

# Returns the downstream segment
fncGetSegDownstream<- function(seg = seg, ssn = parent.frame()$ssn) {
  
  # bottom seg already
  if(seg == min(ssn@data$rid)) return(seg)
  
  # loop over all segs to identify the seg whose first coordinates match the original seg's last coordinates
  for (x in min(ssn@data$rid):max(ssn@data$rid)) {
    coords.x = fncGetSegCoords(seg = x, ssn)
    coords.seg = fncGetSegCoords(seg = seg, ssn)
    nrows.seg = nrow(coords.seg)
    if(coords.x[1,1] == coords.seg[nrows.seg,1] & coords.x[1,2]==coords.seg[nrows.seg,2]) seg = x
  }
  if(!seg %in% missed){
    return(seg)
  } else{
    # call function recursively until find a seg not in 'missed'
    return(fncGetSegDownstream(seg))
  }
}

# Called by fncGetSegUpstream
fncGetNextUpSeg<- function(seg = seg, ssn = parent.frame()$ssn){
  # loop over all segs to identify the seg whose last coordinates match the original seg's first coordinates
  mySegs = c()
  for (x in min(ssn@data$rid):max(ssn@data$rid)) {
    coords.x = fncGetSegCoords(seg = x, ssn)
    coords.seg = fncGetSegCoords(seg = seg, ssn)
    nrows.x = nrow(coords.x)
    if(coords.x[nrows.x,1] == coords.seg[1,1] & coords.x[nrows.x,2] == coords.seg[1,2]) mySegs = c(mySegs,x)
  }
  if(is.null(mySegs)) mySegs = NA
  return(mySegs)
}

# Returns the upstream segment(s)
fncGetSegUpstream<- function(seg = seg, ssn = parent.frame()$ssn) {
  
  # headwater seg already
  if(is.na(seg)) return(seg)
  
  upSegs = fncGetNextUpSeg(seg)

  # check for 'missed' segs
  while(length(upSegs) == 1 & upSegs[1] %in% missed){
    upSegs = fncGetNextUpSeg(upSegs)
  }
  while(any(upSegs %in% missed)){
    for(s in upSegs){
      if(s %in% missed){
        new = fncGetNextUpSeg(s)
        upSegs = c(new, upSegs)
        if(s %in% missed) upSegs = upSegs[!upSegs == s]

      }
    }
  }

  # remove any NAs caused by headwater reaches
  upSegs = upSegs[!is.na(upSegs)]
  
  return(upSegs)
}

# Get the x,y coordinate for a particular location along a straight-line segment in generated stream networks
#    For line coordinates, the higher values (upstream) are the first row in the coordinate set, as if going
#    from x1 to x2 is moving downstream, but for other things, like ratio, it's more like moving upstream
#    Called automatically when 'netnm' begins with "network-".
fncGetXY<- function(seg = seg, ratio = 0.5, ssn = parent.frame()$ssn) {
  # Function is (x1-x2) * proportion + x2
  rowtoget = nrow(fncGetSegCoords(seg, ssn))
  idx = which(ssn@data$rid == seg)
  xloc = (ssn@lines[[idx]]@Lines[[1]]@coords[1,1] - ssn@lines[[idx]]@Lines[[1]]@coords[rowtoget,1]) * ratio + ssn@lines[[idx]]@Lines[[1]]@coords[rowtoget,1]
  yloc = (ssn@lines[[idx]]@Lines[[1]]@coords[1,2] - ssn@lines[[idx]]@Lines[[1]]@coords[rowtoget,2]) * ratio + ssn@lines[[idx]]@Lines[[1]]@coords[rowtoget,2]
  return(cbind(xloc,yloc))
}

# This is the curvy line version for segments in stream networks generated in GIS
#   It operates by moving fish to closest existing node within a reach (i.e., coordinates of the points along
#   a line that was created in GIS); This version is called if 'netnm' does not begin with "network-".
fncGetXY.arc<- function(seg = seg, ratio = 0.5, ssn = parent.frame()$ssn) {
  idx = which(ssn@data$rid == seg)
  ratio = 1 - ratio #may need to check for your network if it was generated in GIS
  if(ratio == 0) ratio = 0.0001 #can't be zero for while statement below
  ll = LineLength(ssn@lines[[idx]]@Lines[[1]]@coords,sum = FALSE)
  lL = LineLength(ssn@lines[[idx]]@Lines[[1]]@coords,sum = TRUE)
  l.idx = 0; counter = 1
  while(l.idx < ratio * lL){l.idx = l.idx + ll[counter]; counter = counter + 1}
  foo = c(lL - sum(ll[1:(counter - 1)]), lL - sum(ll[1:(counter - 2)])); l.new = which.min(foo)
  ifelse(l.new == 1, new.idx<- counter, new.idx<- counter - 2)
  xloc = ssn@lines[[idx]]@Lines[[1]]@coords[new.idx,1]
  yloc = ssn@lines[[idx]]@Lines[[1]]@coords[new.idx,2]
  return(cbind(xloc,yloc))
}

# Get stream order for a seg
fncGetSO<- function(s, ssn = parent.frame()$ssn) {
  return(wq.df[,so.field][wq.df$rid %in% s][1])
}

# Trace upstream
fncTraceUp<- function(seg, ssn = parent.frame()$ssn, plotit = TRUE, col = "turquoise", no2get = 50000){
  
  ups = fncGetSegUpstream(seg)
  set = c(seg,ups)
  new = ups
  
    while(length(ups) > 0 & length(set) < no2get){
    for(u in 1:length(ups)){
        new = c(new, fncGetSegUpstream(seg = ups[u]))
      }
      new = new[! new %in% set]
      ups = new
      set = c(set, new)
    }

  set = set[set != seg]
  
  if(plotit == TRUE){ for(x in set){ fncHighlightSeg(x, col = col) } }
  
  return(set)
}    

# Trace downstream
fncTraceDn<- function(seg, ssn = parent.frame()$ssn, plotit = TRUE, col = "yellow", no2get = 50000){
  
  dn = fncGetSegDownstream(seg)
  set = c(seg, dn)
  new = dn
  
  while(length(dn) > 0 & dn > 1 & length(set) < no2get){
    for(u in 1:length(dn)){
      new = c(new, fncGetSegDownstream(seg = dn[u]))
    }
    new = new[! new %in% set]
    dn = new
    set = c(set, new)
  }
  set = set[set != seg]
  
  if(plotit == TRUE){ for(x in set){fncHighlightSeg(x, col = col)} }
  
  return(set)
}

# Highlight a particular fish
fncHighlightFish<- function(fish, fishID, col = "turquoise", cex = 2, pch = 20) {
  points(fish$xloc[fish$pid == fishID], fish$yloc[fish$pid == fishID], col = col, cex = cex, pch = pch)
}

# Highlight certain segments
fncHighlightSeg<- function(seg, ssn = parent.frame()$ssn, col = "turquoise", lwd = 2) {
  lines(fncGetSegCoords(seg)[,1], fncGetSegCoords(seg)[,2], col = col, lwd = lwd)
}

# Label chosen segment
fncLabelSeg<- function(seg, ssn = parent.frame()$ssn, col = 1, lwd = 1, label = "seg") {
  lines(fncGetSegCoords(seg)[,1], fncGetSegCoords(seg)[,2], col = col, lwd = lwd)
  #label seg:
  if(label == "seg") text(mean(fncGetSegCoords(seg)[,1]), mean(fncGetSegCoords(seg)[,2]), seg, cex = 0.8)
  #label stream order:
  if(label == "so") text(mean(fncGetSegCoords(seg)[,1]), mean(fncGetSegCoords(seg)[,2]), ssn@data[ssn@data$rid == seg, so.field], cex = 0.8)
}

# Get closest pred points for all branches of a junction at the upper end of a given segment where water temperature values exist
  # only needed for stream networks created in GIS, not generated networks
  # for junctions at confluences, function will return the top point on the segment itself and lowest points on each of the two upstream segments for which predictions exist
  # for a junction of two arcs (no branches), function will return the the top point on the segment and the lowest point on the upstream segment for which predictions exist
  # if a segment is a headwater segment, function will return the top point on the segment itself for which predictions exist
  # if provided a 'missed' seg, the function will remove that missed seg from the result and begin with the next upstream seg
fncGetJunctionPreds<- function(seg, ssn = parent.frame()$ssn){
  mySegs = c(seg, fncGetSegUpstream(seg))
  out = matrix(NA,length(mySegs),2)
  nn = 1
  for(r in mySegs){
      # for the trunk, we want the upstream-most point, for branches we want the downstream-most point
      tmp.df = ssn@predpoints@SSNPoints[[1]]@point.data[ssn@predpoints@SSNPoints[[1]]@point.data$rid == r, c("rid","pid","ratio")]
      ifelse(r == mySegs[1], pt<- tmp.df[which.max(tmp.df$ratio), "pid"], pt<- tmp.df[which.min(tmp.df$ratio), "pid"])
      out[nn,] = cbind(r,pt)
      nn = nn + 1
  }
  colnames(out) = c("r","pt")
  idx = which(out[,"r"] %in% missed)
  if(length(idx) > 0) out = out[-idx,]
  if(!is.matrix(out)) out = as.matrix(t(out))
  
  return(out)
}

# Get list of downstream segs for each seg
fncDnSegs<- function(ssn = parent.frame()$ssn, path){
  seglist = sort(unique(ssn@data$rid))
  if(min(seglist) == 0) seglist = seglist[-1] #remove the base one because can't index a list on zero
  
  dnsegs = list()
  for(seg in seglist){
    cat(seg, "\n")
    sg = fncTraceDn(seg, plotit = FALSE)
    dnsegs[[seg]] = sg
  }
  
  save("dnsegs",file = paste0(path, "/dnsegs.RData"))
  return(dnsegs)
}

# Get list of Upstream segs for each seg
fncUpSegs<- function(ssn = parent.frame()$ssn, path){
  seglist = sort(unique(ssn@data$rid))
  if(min(seglist) == 0) seglist = seglist[-1] #remove the base one because can't index a list on zero
  
  upsegs = list()
  for(seg in seglist){
    cat(seg,"\n")
    sg = fncTraceUp(seg, plotit = FALSE)
    upsegs[[seg]] = sg
  }
  
  save("upsegs",file=paste0(path, "/upsegs.RData"))
  return(upsegs)
}

# Get list of segs at junctions for each seg
fncJctLst<- function(ssn = parent.frame()$ssn, path){
  seglist = sort(unique(ssn@data$rid))
  if(min(seglist) == 0) seglist = seglist[-1] #remove the base one because can't index a list on zero
  
  jct.list = list()
  for(seg in seglist){
    jct.list[[seg]] = fncGetJunctionPreds(seg)
  }
  
  save("jct.list",file=paste0(path, "/jct.list.RData"))
  return(jct.list)
}

# Rescale a vector to a specific range
fncRescale<- function(x, to = c(0, 1), from = range(x, na.rm = TRUE, finite = TRUE)) {
    (x - from[1]) / diff(from) * diff(to) + to[1]
}

# Get Julian month from Julian day
fncJulianMonth<- function(x){
  if(x >= 0){
    if(x <= 31) return(1)
    else if(x >= 32 & x < 60) return (2)
    else if(x >= 60 & x< 91) return (3)
    else if(x >= 91 & x< 121) return (4)
    else if(x >= 121 & x< 152) return (5)
    else if(x >= 152 & x< 182) return (6)
    else if(x >= 182 & x< 213) return (7)
    else if(x >= 213 & x< 244) return (8)
    else if(x >= 244 & x< 274) return (9)
    else if(x >= 274 & x< 305) return (10)
    else if(x >= 305 & x< 335) return (11)
    else if(x >= 335 & x< 366) return (12)
  } else if(x < 0){
    if(-x <= 31) return(12)
    else if(-x >= 32 & - x < 60) return (11)
    else if(-x >= 60 & - x < 91) return (10)
    else if(-x >= 91 & - x < 121) return (9)
    else if(-x >= 121 & - x  <152) return (8)
    else if(-x >= 152 & - x  <182) return (7)
    else if(-x >= 182 & - x  <213) return (6)
    else if(-x >= 213 & - x  <244) return (5)
    else if(-x >= 244 & - x  <274) return (4)
    else if(-x >= 274 & - x  <305) return (3)
    else if(-x >= 305 & - x  <335) return (2)
    else if(-x >= 335 & - x  <366) return (1)
    
  }
}

# Get iteration
fncGetRun<- function(){
  run= as.character(Sys.time())
  run= gsub("-", ".", run)
  run= gsub(":", ".", run)
  run= gsub(" ", ".", run)
  return(run)
}

#*******************************************************************************
#=== HABITAT FUNCTIONS ==========================================================
# Calculate useable width of a stream reach based on width
  # scales down usability of wide mainstems but leaves tribs as 100% useable
fncUseableWidths<- function(dat = parent.frame()$ssn@data[,c("rid", "WIDTH_M")]){
  widths = dat[,2]
  rid = dat[,1]
  prop.useable = 1 - fncRescale(log(widths), to = c(0, 0.9)) #assuming 10% of really wide rivers is still useable, and 100% of tribs is useable habitat
  useable.widths = widths * prop.useable
  return(cbind(rid, widths, useable.widths))
}

# Get attribute nearest to a fish
  # 'field' can be "WT", "ration", "SO" etc. 
fncGetNearestAttribute<- function(dat.df = dat.df[,c("pid", "rid", "ratio", field)],
                  fish = fish[,c("pid", "seg", "ratio")], ssn = parent.frame()$ssn){
  field = colnames(dat.df)[4]
  fish = as.matrix(fish)
  dat = as.matrix(dat.df)
  attribs = matrix(NA, length(unique(fish[,"pid"])), 4)
  thesegs = unique(fish[,"seg"]) #each segment that has fish in it
  nn = 1
  
  for(r in thesegs){
    n.fish = length(fish[fish[,"seg"] == r, "pid"])
    tmp.fish = matrix(fish[fish[,"seg"] == r], n.fish, 3) #get fish in that seg
    colnames(tmp.fish) = c("pid", "seg", "ratio")
    tmp.dat = dat[dat[,"rid"] == r,] #get attribute data for the seg
    if(is.vector(tmp.dat)) tmp.dat = matrix(tmp.dat, 1, 4)
    colnames(tmp.dat) = c("pid", "rid", "ratio", field)
    for(i in 1:nrow(tmp.fish)){ #for each of these fish
      idx = which.min(abs(tmp.fish[i,"ratio"] - tmp.dat[,"ratio"]))
      wt = tmp.dat[idx, field] #record the attribute value in the seg that is closest to the fish's position
      attribs[nn,] = cbind(tmp.fish[i, "pid"], r, tmp.fish[i, "ratio"], wt)
      nn = nn + 1
    }
  }
  
  colnames(attribs) = c("pid", "seg", "ratio", field)
  attribs = attribs[order(attribs[,"pid"]),,drop = FALSE]
  return(attribs)
}

# Calculate the density of species1, species2, and total fish in all accessible reaches based on fish biomass or fish counts
fncFishDensity<- function(fish1, fish2, fish1.idx, fish2.idx, ssn = parent.frame()$ssn){
  
  rid = ssn@data$rid[ssn@data$accessible ==1]
  length = ssn@data[,length.field][ssn@data$accessible ==1]
  width = ssn@data$UseableWidth[ssn@data$accessible ==1]
  area = length * 1000 * width
  # multiply by 1000 because length is in km, width is in m
  # area is in m^2
  
  density.species1.num = density.species2.num =  density.total.num = vector("numeric", length(rid))
  density.species1.bio = density.species2.bio =  density.total.bio = vector("numeric", length(rid))
  for(i in 1:length(rid)){ 
    r = rid[i]
    seg.area = area[i]
    #only include fish that are alive and emerged and in this reach
    fsh1.idx = fish1.idx[fish1[fish1.idx,"seg"] == r]
    fsh2.idx = fish2.idx[fish2[fish2.idx,"seg"] == r]
    num.fish1.in.seg = length(fish1[fsh1.idx,1])
    num.fish2.in.seg = length(fish2[fsh2.idx,1])
    density.species1.num[i] = num.fish1.in.seg / seg.area
    density.species2.num[i] = num.fish2.in.seg / seg.area
    biomass.fish1.in.seg = sum(fish1[fsh1.idx,"weight"]) 
    biomass.fish2.in.seg = sum(fish2[fsh2.idx,"weight"]) 
    density.species1.bio[i] = biomass.fish1.in.seg / seg.area
    density.species2.bio[i] = biomass.fish2.in.seg / seg.area
  }
  
  density.total.num = density.species1.num + density.species2.num
  density.total.bio = density.species1.bio + density.species2.bio
  
  density.num = matrix(c(rid, density.species1.num, density.species2.num, density.total.num), nrow = length(rid), ncol = 4)
  density.bio = matrix(c(rid, density.species1.bio, density.species2.bio, density.total.bio), nrow = length(rid), ncol = 4)
  density = cbind(density.bio, density.num[,c(2:4)])
  
  #apply the "superindividual" correction factor that makes densities more representative of observed fish densities
  density[,c(2:7)] = density[,c(2:7)] * parameters[1,"corr.factor"] 

  colnames(density) = c("seg", "conspecific.bio.density", "other.bio.density", "total.bio.density", "conspecific.num.density", "other.num.density", "total.num.density")
  rownames(density) = rid
  
  return(density)
}

#=== MOVEMENT FUNCTIONS ========================================================

# Calculates movement distance potentials based on conditions at current location vs. elsewhere
# returns data frame of max movement distances for each fish
fncMoveDistance<- function(fish = parent.frame()$fish, sp.idx, ssn = parent.frame()$ssn){
  
  # get global variables
  mvmt.scalar = parameters[sp.idx, "mvmt.scalar"]
  mvdist.shape = parameters[sp.idx, "mvdist.shape"]
  wt.field = get("wt.field")
  fpids = fish[,"pid"]
  nFish = length(fpids)
  
  # look up growth at fish’s location based on its WT, ration, and weight
  growths = matrix(NA, nFish, 3, dimnames = list(fpids, c("growth", "growth.pot", "growth.min")))
  growths[,"growth"] = fncGrowthFish(fpids, fish, sp.idx)
  
  # look up max growth possible for each fish across all accessible reaches during this time step
  tmp.lst = list()
  tmp.lst$result = lapply(fpids, function(x) fncGrowthPossible(fweight = fish[,"weight"][fish[,"pid"] == x], sp.idx, ssn))
  tmp.range = do.call(rbind,lapply(tmp.lst$result, range))
  growths[,"growth.pot"] = tmp.range[,2]
  growths[,"growth.min"] = tmp.range[,1]

  #movement probability and distance is influenced by how close fish's current growth is to max growth potential
  #bigger values mean higher probability of longer movements
  gr.diff = apply(growths[,c("growth","growth.pot"), drop = FALSE], 1, diff) / apply(growths[,c("growth.min","growth.pot"), drop = FALSE], 1, diff)

  # draw from lognormal distribution where mean is the growth difference and sd is a parameter
  # movement multiplier parameter is used to increase movement in certain scenarios
  moveDist<- rlnorm(nFish, meanlog = (gr.diff * mvmt.scalar), sdlog = mvdist.shape)
  
  if(SecondSpecies == TRUE) moveDist = moveDist * mvmt.scalar #high site fidelity
  
  moveDist <- cbind(fpids, moveDist)
  colnames(moveDist) = c("pid", "moveDist")
  
  return(moveDist)
}

# Determines if a fish can stop before its assigned movement distance if it encounters good habitat; returns list(T/F, pStop, dist2move, movedirection)
fncStopEarly<- function(fpid, seg = parent.frame()$seg, remainingDist, moveLength, length2segBase, om.prob2 = 1, fish = parent.frame()$fish, ssn = parent.frame()$ssn){
  # the probability of stopping increases as growth potential increases
  # but probability of stopping decreases if fish are motivated to outmigrate
  
    length.field = get("length.field")
  
    # look up growth possible for the fish at all locations in this reach
    result = fncGrowthPossInReach(fish[,"weight"][fish[,"pid"] == fpid], sp.idx, seg = seg)
    growth.here = max(result)
    # for an SSN with multiple points per reach, this syntax matters
    best.obs = as.numeric(names(result)[which.max(result)])
    
    # look up growth possible for each fish across accessible reaches during this time step
    result2 = fncGrowthPossible(fish[,"weight"][fish[,"pid"] == fpid], sp.idx)
     
    
    # probability of stopping due to good habitat (growth in reach relative to max growth potential in stream network)
      # bigger values mean higher probability of stopping
    probStop1 = 1 - diff(c(growth.here, max(result2))) / diff(range(result2))
    if(probStop1 > 1) probStop1 = 1 # can't be higher than 1
    
    # probability of stopping due to fish's drive to outmigrate
      # bigger values mean higher probability of stopping
    probStop2 = 1 - om.prob2
    
    # cumulative probability of stopping
    probStop = probStop1 * probStop2
    
    #prob. of stopping, prob. of continuing
    pStop<- c(probStop, (1 - probStop)) 
  
    # determine position of fish relative to obs to stop at & update distance to move
    obs.seg.ratio = wq.df$ratio[wq.df$rid == seg & wq.df$pid == best.obs] #seg ratio at best.obs
    obs.dist2Base = ssn@data[ssn@data[,"rid"] == seg, length.field] * obs.seg.ratio
    ifelse(length2segBase >= obs.dist2Base, moveDirection<- 1, moveDirection<- 0) 
    dist2move = abs(length2segBase - obs.dist2Base)
    
    StopEarly<- sample(c(T, F), size = 1, prob = pStop) #T=stop, F=continue
    
    return(list(StopEarly, pStop, dist2move, moveDirection))
}

# Determines if there's room within the segment to move; returns T/F
fncRoom2Move<- function(seg = parent.frame()$seg, moveDirection, remainingDist, length2segBase, ssn = parent.frame()$ssn) {
  
  length.field = get("length.field")
  if (moveDirection == 1 & remainingDist <= length2segBase) {
    return(TRUE)
  } else if (moveDirection == 0 & remainingDist <= (ssn@data[ssn@data$rid == seg, length.field] - length2segBase))  {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

# Determines which segment a fish moves into as it encounters junctions; returns 'seg'
fncJunctionDecisions<- function(fpid, sp.idx, seg = parent.frame()$seg, currMoveDirection, stochastic = TRUE, fish = parent.frame()$fish, ssn = parent.frame()$ssn) {
  # As a fish gets to a junction, whether heading up or downstream,
  # it looks at all possible options for movement and chooses the best one
  
  # Always look at the junction upstream of a particular segment,
  # So, if a fish is heading downstream
  # it's easiest to go to the downstream segment and look at all options from there
  if(currMoveDirection == 1){
    if(seg > min(ssn@data$rid)){ seg = fncGetSegDownstream(seg)}
  }
  
  if(seg == 0){ return(bestSeg = 1)} #can't go farther downstream
  
  # Get conditions at downstream segment and at each of the upstream segments
  
  # First get which reaches participate in the junction
  rids = unique(jct.list[[seg]][,"r"])

  # Check that these reaches are accessible to fish & remove any that aren't
  rg.col = parameters[sp.idx, "spp.range"]
  idx = ssn@data[,rg.col][ssn@data$rid %in% rids]
  rids = rids[idx == 1]
  
  # If any accessible segs in 'rids'
  if(length(rids) > 1){
    # Look up previously calculated growth values
    growth.confluence = fncGrowthPossible(fish[,"weight"][fish[,"pid"] == fpid], sp.idx, ssn, segs = rids)
  
  # Choose based on growth in each seg option
    gr.segs = tapply(growth.confluence, names(growth.confluence), mean)

    if(any(gr.segs < 0)) gr.segs = gr.segs + 2 * abs(min(gr.segs)) #ensure no negative probabilities
    
    # choose best seg
    if(stochastic == TRUE){
      # base probabilities on relative growth potential in each reach
      ifelse(sum(gr.segs) > 0, probs<- gr.segs / sum(gr.segs), probs<- rep(0.5, length(gr.segs)))
      probs<- gr.segs / sum(gr.segs)
      ifelse(length(gr.segs) > 1, bestSeg<- sample(as.integer(names(gr.segs)), size = 1, prob = probs), bestSeg<- seg)
    } else {
      # choose seg deterministically
      bestSeg = as.integer(names(which.max(gr.segs)))
    }
  
  } else {
    bestSeg = seg
  }
  return(bestSeg)
}

# Determines a fish's position based on it's current circumstances; returns list(dist2move,remainingDist,length2segBase,segRatio, moveDirection)
  #where seg was original seg and newSeg is the result of fncJunctionDecisions
fncUpdatePosition<- function(newSeg, seg = parent.frame()$seg, sp.idx, moveDirection, length2segBase, dist2move, remainingDist, 
                             SE = StopEarly, RM = Rm2Move, ssn = parent.frame()$ssn){
  
  if(newSeg == 0) return(list(dist2move = 0, remainingDist = 0, length2segBase = 0, segRatio = 0, moveDirection = 1)) #fish is at base of network
  
  if(SE == TRUE | RM == TRUE){ #if fish stopped in the reach, update its position
    
    if (moveDirection == 1) {
      length2segBase = length2segBase - dist2move 
    } else if (moveDirection == 0) {
      length2segBase = length2segBase + dist2move 
    }
    # Then update it's segRatio
    segRatio = fncLength2segRatio(seg, length2segBase)
    remainingDist = 0
    #seg and moveDirection stay the same
    
  } else { #move fish to next upstream or downstream reach and update its position
    
    if(moveDirection == 0){ # if heading upstream
      
      dist2move = ssn@data[ssn@data$rid == seg,length.field] - length2segBase #travel the rest of the way up the original reach
      
      # if no accessible upstream reaches or if best choice is to stay in same reach,
      # then locate fish at top of existing reach facing downstream
      if(sum(ssn@data$accessible[ssn@data$rid %in% fncGetSegUpstream(seg)]) == 0 | newSeg == seg){
        segRatio = 1
        length2segBase = ssn@data[ssn@data$rid == seg, length.field] #length of current reach
        moveDirection = 1
        
      } else{ #otherwise, 
        
        #if newSeg is upstream of seg or if both seg and newSeg are upper branches
        ifelse (fncGetSegDownstream(seg) > 0, test.seg<- fncGetSegUpstream(fncGetSegDownstream(seg)), test.seg<- fncGetSegUpstream(seg)) #(no seg zero)
        if(newSeg %in% fncGetSegUpstream(seg) | newSeg %in% test.seg){ 
          # then locate fish at the base of the new reach and maintain original direction
          segRatio = 0 
          length2segBase = 0
        } 
      }
      
    } else if(moveDirection == 1){ # if heading downstream
      
      dist2move = length2segBase #travel to base of the reach
      
      # if no accessible downstream reaches or if best choice is to stay in same reach,
      # then locate fish at bottom of existing reach facing upstream
      if(newSeg <= min(ssn@data$rid) | seg == newSeg){
        segRatio = 0 
        length2segBase = 0
        moveDirection = 0
        if(newSeg < min(ssn@data$rid)) newSeg = min(ssn@data$rid)
        
      } else{ #otherwise, 
        
        if(newSeg %in% fncGetSegDownstream(seg)){ #if newSeg is downstream of seg
          # then locate fish at the top of the new reach and maintain original direction
          segRatio = 1
          length2segBase = ssn@data[ssn@data$rid == newSeg, length.field]
          
        } else { #i.e., newSeg and seg are both upper branches
          # then locate fish at the bottom of the new reach facing upstream
          segRatio = 0
          length2segBase = 0
          moveDirection = 0
        }
        #seg = newSeg
      }
    } #end move up or down
    
    #update distance remaining (i.e., distance a fish can still move this time step)
    remainingDist = remainingDist - dist2move

  } #end stop early or continue
  
  return(list(dist2move, remainingDist, length2segBase, segRatio, moveDirection)) 
}

# Move fish individually (one fish, one time step)
fncMoveIndividual<- function(fpid, fish = parent.frame()$fish, sp.idx, ssn = parent.frame()$ssn) {
    set.seed(fpid) 
    
    #Get starting reach and location info
    seg = fish$seg[fish$pid == fpid]
    segLength = ssn@data[ssn@data$rid == seg, length.field]
    length2segBase = fish$length2segBase[fish$pid == fpid]
    segRatio = fncLength2segRatio(seg, length2segBase)
    initDirection = moveDirection = fish$direction[fish$pid == fpid]
    ifelse(sp.idx == 1, om.prob2<- get("om.prob")[fish$pid == fpid], om.prob2<- 0)
    remainingDist = moveLength = moveDist[,"moveDist"][moveDist[,"pid"] == fpid]
    dist2move = dist.moved.total = 0
    
    while(remainingDist > 0){ #Continue as long as the fish has not yet moved its allocated distance
      
      if(seg == 0){
        seg = network.base.segs[1]
      }
      
      #Determine if there is room within the reach to move
      Rm2Move = fncRoom2Move(seg, moveDirection, remainingDist, length2segBase, ssn)
      if(Rm2Move == TRUE){
        dist2move = remainingDist
      }
      
      # Determine if fish will stop early due to good thermal habitat somewhere in the reach
      # can happen even if room to move is false
        StopResult = fncStopEarly(fpid, seg, remainingDist, moveLength, length2segBase, om.prob2, fish, ssn)
        StopEarly = StopResult[[1]] #Did the fish stop early in good habitat? T=stop, F=continue
        if(StopEarly == TRUE) {
          dist2move = StopResult[[3]]
          moveDirection = StopResult[[4]]
        }

      # If fish is not stopping in the reach, determine which will be the next reach it enters
      # (only allowed to enter accessible reaches; allows fish to change direction at any segment edge, not just confluences)
      if (StopEarly == FALSE & Rm2Move == FALSE) {
          newSeg = fncJunctionDecisions(fpid, sp.idx, seg, moveDirection, stochastic = FALSE, fish, ssn)
        } else {newSeg = seg}

      # Update fish's position
      # to either somewhere within the reach, or the rest of the way up/down the reach to the junction   
      Position = fncUpdatePosition(newSeg, seg, sp.idx, moveDirection, length2segBase, dist2move, remainingDist, StopEarly, Rm2Move, ssn)
      dist2move = Position[[1]]
      remainingDist = Position[[2]]
      length2segBase = Position[[3]]
      segRatio = Position[[4]]
      seg = newSeg
      moveDir.orig = moveDirection
      moveDirection = Position[[5]]
      dist.moved.total = dist.moved.total + dist2move
      
      # Stop if errors
      if(remainingDist < 0) browser()
      if(segRatio > 1 | segRatio < 0) browser()
      if(seg == 0) browser()
      
    } # end while remaining dist loop
    
    # Put info back in fish table
    fish$seg[fish$pid == fpid] = seg
    fish$length2segBase[fish$pid == fpid] = length2segBase
    fish$ratio[fish$pid == fpid] = segRatio
    fish$direction[fish$pid == fpid] = moveDirection
    if(seg > 0){
      fish$upDist[fish$pid == fpid] = fncDist2base4fish(fpid, fish)
      fish$SO[fish$pid == fpid] = fncGetSO(seg)
    } else { 
      fish$upDist[fish$pid == fpid] = 0
      fish$SO[fish$pid == fpid] = fncGetSO(1)
    }
      fish$movedist[fish$pid == fpid] = dist.moved.total
      
    return(fish[fish$pid == fpid,,drop = FALSE])
}

# Bias movement direction downstream as fish grows larger 
  # but after a date threshold, fish loses the outmigration drive and becomes a yearling
fncDownstreamDrive<- function(w, om.mass = 1.5, om.date.taper = "04-01"){
  om.date.taper = as.Date(paste0(climate.year, "-", om.date.taper))
  if(dat.idx[dd] > om.date.taper){
    om.prob = (w / om.mass) * (1 - as.numeric(dat.idx[dd] - om.date.taper) / as.numeric(dat.idx[length(dat.idx)] - om.date.taper))
  } else{
    om.prob = w / om.mass
  }
  # constrain between 0.5 (equal chance of moving up or down) and nearly 1 (maximum downstream drive)
  om.prob[om.prob >= 1] = 0.999; om.prob[om.prob < 0.5] = 0.5
  
  return(om.prob)
}

# Jitter a fish's position once
fncJitterPosition<- function(wt.growth, fae, fish = parent.frame()$fish, sp.idx, ssn = parent.frame()$ssn, plotit = "none"){
  
  # Create optional plot showing original fish locations [part 1]
  if(plotit != "none"){
    if(plotit != "screen"){png(paste0(plotDir,"/",plotit,".png"), width = 6, height = 6, units = "in", res = 300)}
    plot(basin,col="gray80",border=NA, main = "")
    plot(streams, col="darkgray", add = TRUE)
    points(ssn@predpoints@SSNPoints[[sp.idx + 1]]@point.coords, col = 1)
  }
  
  # set global variables that will be called inside fncMoveIndividual
  sp.idx <<- sp.idx
  wt.growth <<- wt.growth

  # get uniform random movement distances
  moveDist <<- cbind("pid" = fish[,"pid"], "moveDist" = runif(nrow(fish), min = 0.01, max = 0.5))
  
  # call movement function once
  moved = lapply(fish[,"pid"], function(x) fncMoveIndividual(fpid = x, fish, sp.idx, ssn))
  fish = do.call(rbind, moved)
  
  # put results back into data frame and SSN
  ifelse(length(grep("network-", netnm)) > 0, 
         locs<- ddply(fish, .(pid), function(x) data.frame(fncGetXY(x$seg, x$ratio, ssn))), 
         locs<- ddply(fish, .(pid), function(x) data.frame(fncGetXY.arc(x$seg, x$ratio, ssn))))
  fish$xloc<- locs$xloc; fish$yloc<- locs$yloc
  ssn@predpoints@SSNPoints[[sp.idx + 1]]@point.data$upDist[fae] = fish$upDist
  ssn@predpoints@SSNPoints[[sp.idx + 1]]@point.data$ratio[fae] = fish$ratio
  ssn@predpoints@SSNPoints[[sp.idx + 1]]@point.data[,xlat][fae] = fish$xloc
  ssn@predpoints@SSNPoints[[sp.idx + 1]]@point.data[,ylon][fae] = fish$yloc
  ssn@predpoints@SSNPoints[[sp.idx + 1]]@point.data$rid[fae] = fish$seg
  ssn@predpoints@SSNPoints[[sp.idx + 1]]@point.coords[,1][fae] = fish$xloc 
  ssn@predpoints@SSNPoints[[sp.idx + 1]]@point.coords[,2][fae] = fish$yloc
  
  # turn fish back upstream if it's reached the base of the network
  fish$direction[fish$seg == network.base.segs[1] & fish$direction == 1] = 0
  
  # Create optional plot showing new fish locations [part 2]
  if(plotit != "none"){
    points(ssn@predpoints@SSNPoints[[sp.idx + 1]]@point.coords[fish_other.alive.index,], col = 2)
    if(plotit != "screen"){dev.off()}
  }
  
  return(fish)
}

#*******************************************************************************

#=== BIOENERGETICS FUNCTIONS ===================================================

# Look up growth based on actual WT, ration, and weight of fish
# PIDs is vector of fish or a single fish
fncGrowthFish<- function(PIDs = NA, fish = parent.frame()$fish, sp.idx = 1){
  
  #these need to match what was used when pre-calculating growth 
  #(they define the dimensions of the array that the precalculated data are stored in)
  wt.seq = seq(0.05, 25, 0.05) #water temperature
  ra.seq = seq(0.01, 0.3, 0.001) #ration
  ma.seq = seq(1, 500, 1) #fish mass
  if(sp.idx == 2) ma.seq = c(1:100,seq(105,995,5), seq(1000,6000,10)) else ma.seq = seq(1, 500, 1) #fish mass

  #get data
  if(any(is.na(PIDs))){
    td = fish[, c("pid", "weight", ration.field, wt.field)] #all fish
  } else {
    td = fish[fish$pid %in% PIDs, c("pid", "weight", "ration", wt.field)] #specific fish
  }
  td[,"weight"][td[,"weight"] < 1] = 1 #minimum value that can be looked up is 1 g
  
  if(nrow(td) > 0){ #ensure there are data before proceeding
    growth = watemp = vector(length = length(PIDs)) #create empty vectors
    #for each fish, lookup pre-calculated growth
    for(x in 1:length(PIDs)){
      wt.idx = which.min(abs(wt.seq - td[,wt.field][x])) #which water temperature is closest to fish's?
      ra.idx = which.min(abs(ra.seq - td[,"ration"][x])) #which ration is closest to fish's?
      ma.idx = which.min(abs(ma.seq - td[,"weight"][x])) #which weight is closest to fish's?
      growth[x] = wt.growth[wt.idx, ra.idx, ma.idx] #use these indices to look up pre-calculated growth
      watemp[x] = wt.seq[wt.idx]
    }
    if(any(is.na(PIDs))) PIDs = fish[,"pid"]

    return(growth)
    
  } else{
    cat(fpid, "\n")
    browser()
    return(NA)
  }
  
}

# Look up best possible location for growing based on all current temperature/ration across stream network
fncGrowthPossible<- function(fweight = 1, sp.idx = 1, ssn = parent.frame()$ssn, segs = NA){

  #these need to match what was used when pre-calculating growth 
  #(they define the dimensions of the array that the precalculated data are stored in)
  wt.seq = seq(0.05, 25, 0.05) #water temperature
  ra.seq = seq(0.01, 0.3, 0.001) #ration
  if(sp.idx == 2) ma.seq = c(1:100,seq(105,995,5), seq(1000,6000,10)) else ma.seq = seq(1, 500, 1) #fish mass
  if(sp.idx == 2) ration.field = "ration_ss" else ration.field = "ration"
  
  # Get data, for seg(s) 
  if(any(is.na(segs))){
    td = wq.df[wq.df$rid %in% ssn@data$rid[ssn@data$accessible == 1], c("pid", "rid", ration.field, wt.field)]
  } else {
    td = wq.df[wq.df$rid %in% segs & wq.df$rid %in% ssn@data$rid[ssn@data$accessible == 1], c("pid", "rid", ration.field, wt.field)]
  }
  row.names(td) = NULL
  
  if(nrow(td) > 0){ #ensure there are data before proceeding
    # for networks that have multiple points per seg:
    td = aggregate(td, by = list(td[,"rid"]), FUN = mean)[,3:5]

    # for each location in stream network at this time, lookup pre-calculated growth potential
    wt.idx = apply(td[,wt.field, drop = FALSE], 1, function(x) which.min(abs(wt.seq - x)))
    ra.idx = apply(td[,ration.field, drop = FALSE], 1, function(x) which.min(abs(ra.seq - x)))
    mm = which.min(abs(ma.seq - round(fweight, 2)))
    ma.idx = rep(mm, nrow(td))
    ma.idx[ma.idx < 1]<- 1
    wrm.idx = cbind(wt.idx, ra.idx, ma.idx)
    growth = apply(wrm.idx, 1, function(x) wt.growth[x[1],x[2],x[3]])
    watemp = wt.seq[wt.idx]
    
    names(growth) = td[,"rid"]
    return(growth)
    
  } else{
    cat(seg, fpid, "\n")
    browser()
    return(NA)
  }
}

# Look up best possible location to grow within current reach
fncGrowthPossInReach<- function(fweight = 1, sp.idx, ssn = parent.frame()$ssn, seg){
  
  #these need to match what was used when pre-calculating growth 
  #(they define the dimensions of the array that the precalculated data are stored in)
  wt.seq = seq(0.05, 25, 0.05) #water temperature
  ra.seq = seq(0.01, 0.3, 0.001) #ration
  if(sp.idx == 2) ma.seq = c(1:100,seq(105,995,5), seq(1000,6000,10)) else ma.seq = seq(1, 500, 1) #fish mass
  if(sp.idx == 2) ration.field = "ration_ss" else ration.field = "ration"
  
  ## Get data
  td = wq.df[wq.df$rid %in% seg & wq.df$rid %in% ssn@data$rid[ssn@data$accessible == 1], c("pid", "rid", ration.field, wt.field)]
  row.names(td) = NULL
  
  if(nrow(td) > 0){ #ensure there are data before proceeding

    # for each location at this time, lookup pre-calculated growth potential
    wt.idx = apply(td[,wt.field, drop = FALSE], 1, function(x) which.min(abs(wt.seq - x)))
    ra.idx = apply(td[,ration.field, drop = FALSE], 1, function(x) which.min(abs(ra.seq - x)))
    mm = which.min(abs(ma.seq - round(fweight, 2)))
    ma.idx = rep(mm, nrow(td))
    ma.idx[ma.idx < 1]<- 1
    wrm.idx = cbind(wt.idx, ra.idx, ma.idx)
    growth = apply(wrm.idx, 1, function(x) wt.growth[x[1],x[2],x[3]])
    watemp = wt.seq[wt.idx]
    
    names(growth) = td[,"pid"]
    return(growth)
    
  } else{
    cat(seg, fpid, "\n")
    browser()
    return(NA)
  }
  
}



# Functions below were provided by Matt Nahorniak, and adapted for our needs):

# Parameters
fncGetBioEParms <- function(spp, pred.en.dens, prey.en.dens, wt.nearest, startweights = rep(initial.mass, numFish), 
                  pvals = rep(0.5, numFish), ration = rep(0.1, numFish)){
  
  N.sites<- nrow(wt.nearest) #here, sites are fish in each reach
  N.steps<- 1 #one time step
  Species<- spp
  SimMethod<- 1 #method that predicts growth
  Pred<- pred.en.dens #predator energy density
  Oxygen<- 13560 #oxygen consumed (original parameter value from Wisconsin model)
  PFF<- 0.1 #percent indigestible prey (original parameter value from Wisconsin model)
  stab.factor<- 0.5 #stability factor for other simulation methods
  epsilon<- 0.5 #also for other simulation methods
  #startweights = rep(0.01, N.sites)
  endweights = startweights*5
  TotalConsumption = rep(100, N.sites)
  pvalues<- t(matrix(round(pvals,5)))
  sitenames <- t(matrix(wt.nearest[,"pid"]))
  temperature<- t(wt.nearest[,"WT"])
  prey.energy.density<- t(matrix(rep(prey.en.dens,N.sites))) 
  ration<- t(matrix(round(ration,5)))
  
  return(list(	
    "Species" = Species,
    "SimMethod" = SimMethod,
    "Wstart"=startweights,
    "Endweights" = endweights,
    "TotalConsumption" = TotalConsumption,
    "pp"=pvalues, "Temps"=temperature,
    "N.sites"=N.sites,
    "N.steps"=N.steps,
    "sitenames"=sitenames,
    "Pred"=Pred,
    "prey.energy.density"=prey.energy.density,
    "Oxygen"=Oxygen,
    "stab.factor"=stab.factor,
    "PFF" = PFF,
    "epsilon" = epsilon,
    "ration" = ration)
  )
  
}

# Constants, lookup from parameters file
fncReadConstants<- function(spp = "salmon"){
  fish.constants = list(parameters[spp, c("ConsEQ", "CA", "CB", "CQ", "CTO", "CTM", "CTL", "CK1", "CK4")], 
                          parameters[spp, c("RespEQ", "RA", "RB", "RQ", "RTO", "RTM", "RTL", "RK1", "RK4", "ACT", "BACT", "SDA")],
                          parameters[spp, c("ExcrEQ", "FA", "FB", "FG", "UA", "UB", "UG")])
  names(fish.constants) = c("Consumption", "Respiration", "Excretion")
  return(fish.constants)
}

# Consumption Equation 1
ConsumptionEQ1<- function(W, TEMP, PP, RATION, PREY, CA, CB, CQ) {

	CMAX = CA * (W ** CB)			#max specific feeding rate (g_prey/g_pred/d)
	
  if(sum(PP) > 0){ #if using p-values to calculate
  	CONS = CMAX * PP * exp(CQ * TEMP)	#specific consumption rate (g_prey/g_pred/d) - grams prey consumed per gram of predator mass per day
  } else if(sum(RATION) > 0){ #if using ration to calculate
    RATION[RATION >= CMAX] = CMAX[RATION >= CMAX]
    CONS = RATION * exp(CQ * TEMP)     
  }
    
	CONSj = CONS * PREY          		 #specific consumption rate (J/g_pred/d) - Joules consumed for each gram of predator for each day
	
	return(list("CMAX" = CMAX, "CONS" = CONS, "CONSj" = CONSj))
}

# Consumption Equation 2
ConsumptionEQ2<- function(W, TEMP, PP, RATION, PREY, CA, CB, CTM, CTO, CQ) {

	Y = log(CQ) * (CTM - CTO + 2)
	Z = log(CQ) * (CTM - CTO)
	X = (Z^2 * (1 + (1 + 40 / Y)^.5)^2) / 400
	V = (CTM - TEMP) / (CTM - CTO)
	CMAX = CA * (W ** CB)	
	
  if(sum(PP) > 0){ #if using p-values to calculate
  	CONS = CMAX * PP * (V ** X) * exp(X * (1 - V))
  } else if(sum(RATION) >0){ #if using ration to calculate
    RATION[RATION >= CMAX] = CMAX[RATION >= CMAX]
    CONS = RATION * (V ** X) * exp(X * (1 - V))     
  }
	
	CONSj = CONS * PREY             

	return(list("CMAX" = CMAX, "CONS" = CONS, "CONSj" = CONSj))
}

# Consumption Equation #3: Temperature Dependence for cool-cold water species
ConsumptionEQ3<- function(W, TEMP, PP, RATION, PREY, CA, CB, CK1, CTO, CQ, CK4, CTL, CTM) {
  
	G1 = (1 / (CTO - CQ)) * (log((0.98 * (1 - CK1)) / (CK1 * 0.02)))
	L1 = exp(G1 * (TEMP - CQ))
	KA = (CK1 * L1) / (1 + CK1 * (L1 - 1))
	G2 = (1 / (CTL - CTM)) * (log((0.98 * (1 - CK4)) / (CK4 * 0.02)))
	L2 = exp(G2 * (CTL - TEMP))
	KB = (CK4 * L2) / (1 + CK4 * (L2 - 1))
	CMAX = CA * (W ** CB)		#max specific feeding rate (g_prey/g_pred/d)
	
  if(sum(PP) > 0){ #if using p-values to calculate
	  CONS = CMAX * PP * KA * KB		#specific consumption rate (g_prey/g_pred/d) - grams prey consumed per gram of predator mass per day
  } else if(sum(RATION) > 0){ #if using ration to calculate
    RATION[RATION >= CMAX] = CMAX[RATION >= CMAX]
    CONS = RATION * KA * KB     
  }
	
	CONSj = CONS * PREY             #specific consumption rate (J/g_pred/d) - Joules consumed for each gram of predator for each day

	return(list("CMAX" = CMAX, "CONS" = CONS, "CONSj" = CONSj))
	}

# Excretion Equation 1
ExcretionEQ1<- function(CONS, CONSj, FA, UA) {

	EG = FA * CONS				# egestion (fecal waste) in g_waste/g_pred/d
	U = UA * (CONS - EG)	 			# excretion (nitrogenous waste) in g_waste/g_pred/d

	EGj = FA * CONSj				# egestion in J/g/d
	Uj = UA * (CONSj - EGj)			# excretion in J/g/d

	return(list("EG" = EG, "EGj" = EGj, "U" = U, "Uj" = Uj))
	}	

# Excretion Equation 2
ExcretionEQ2<- function(CONS, CONSj, TEMP, PP, RATION, CMAX, FA, UA, FB, FG, UB, UG) {

  if(sum(RATION) > 0){ #if using ration to calculate
    RATION[RATION >= CMAX] = CMAX[RATION >= CMAX]
    PP = RATION / CMAX #calculating p-value based on ration and CMax inputs
  }

	EG = FA * TEMP^FB * exp(FG * PP) * CONS			# egestion (fecal waste) in g_waste/g_pred/d
	U = UA * TEMP^UB * exp(UG * PP) * (CONS - EG)			# excretion (nitrogenous waste) in g_waste/g_pred/d

	EGj = EG * CONSj / CONS					# egestion in J/g/d
	Uj = U * CONSj / CONS					# excretion in J/g/d

	return(list("EG" = EG, "EGj" = EGj, "U" = U, "Uj" = Uj))
	}

# Excretion Equation 3 (W/ correction for indigestible prey as per Stewart 1983)
ExcretionEQ3<- function(CONS, CONSj, TEMP, PP, RATION, CMAX, FA, UA, FB, FG, UB, UG, PFF) {

	#Note: In R, "F" means "FALSE", here we use EG as the variable name for egestion instead of F (as in the FishBioE 3.0 manual)
	#Note:  PFF = 0 assumes prey are entirely digestible, making this essentially the same as Equation 2

  if(sum(RATION) > 0){ #if using ration to calculate
    RATION[RATION >= CMAX] = CMAX[RATION >= CMAX]
    PP = RATION / CMAX #calculating p-value based on ration and CMax inputs
  }
  
	PE = FA * (TEMP ** FB) * exp(FG * PP)

	PF = ((PE - 0.1) / 0.9) * (1 - PFF) + PFF

	EG = PF * CONS					# egestion (fecal waste) in g_waste/g_pred/d
	U = UA * (TEMP ** UB) * (exp(UG * PP)) * (CONS - EG)	# excretion (nitrogenous waste) in g_waste/g_pred/d

	EGj = PF * CONSj				# egestion in J/g/d
	Uj = UA * (TEMP ** UB) * (exp(UG * PP)) * (CONSj - EGj)	# excretion in J/g/d

	return(list("EG" = EG, "EGj" = EGj, "U" = U, "Uj" = Uj))
	}	

# Respiration Equation 1
RespirationEQ1<- function(W, TEMP, CONS, EG, PREY, OXYGEN, RA, RB, ACT, SDA, RQ, RTO, RK1, RK4, RTL, BACT){

	VEL = (RK1 * W^RK4) * (TEMP > RTL) +  ACT * W^RK4 * exp(BACT * TEMP) * (1 - 1 * (TEMP > RTL))
	ACTIVITY = exp(RTO * VEL)
	S = SDA * (CONS - EG)					# proportion of assimilated energy lost to SDA in g/g/d (SDA is unitless)
	Sj = S * PREY						# proportion of assimilated energy lost to SDA in J/g/d - Joules lost to digestion per gram of predator mass per day
	R = RA * (W ** RB) * ACTIVITY * exp(RQ * TEMP)       	# energy lost to respiration (metabolism) in g/g/d
	Rj = R * OXYGEN           				# energy lost to respiration (metabolism) in J/g/d - Joules per gram of predator mass per day
	return(list("R" = R, "Rj" = Rj, "S" = S, "Sj" = Sj))
	}

# Respiration Equation 2 (Temp dependent w/ ACT multiplier)
RespirationEQ2<- function(W, TEMP, CONS, EG, PREY, OXYGEN, RA, RB, ACT, SDA, RTM, RTO, RQ) {

	V = (RTM - TEMP) / (RTM - RTO)
	Z = (log(RQ)) * (RTM - RTO)
	Y = (log(RQ)) * (RTM - RTO + 2)
	X = ((Z ** 2) * (1 + (1 + 40 / Y) ** 0.5) ** 2) / 400
	S = SDA * (CONS - EG)					# proportion of assimilated energy lost to SDA in g/g/d (SDA is unitless)
	Sj = S * PREY						# proportion of assimilated energy lost to SDA in J/g/d - Joules lost to digestion per gram of predator mass per day
	R = RA * (W ** RB) * ACT * ((V ** X) * (exp(X * (1 - V) )))     	# energy lost to respiration (metabolism) in g/g/d
	Rj = R * OXYGEN           				# energy lost to respiration (metabolism) in J/g/d - Joules per gram of predator mass per day
	return(list("R" = R, "Rj" = Rj, "S" = S, "Sj" = Sj))
	}

# Calculate Growth
CalculateGrowth <- function(Constants, Input) {
  
  #attach(Input)
  #attach(Constants)
  
  pred=Input$Pred
  Oxygen=Input$Oxygen
  
  # Initialize Fish Weights
  W = array(rep(0, (Input$N.sites * (Input$N.steps + 1))), c(Input$N.steps + 1, Input$N.sites))
  W[1,] = as.numeric(Input$Wstart[1:ncol(W)])
  Growth =array(rep(0, Input$N.sites*Input$N.steps), c(Input$N.steps, Input$N.sites))
  Growth_j = Growth
  Consumpt = Growth
  Consumpt_j = Growth
  Excret = Growth
  Excret_j = Growth
  Egest = Growth
  Egest_j = Growth   
  Respirat = Growth 
  Respirat_j = Growth 
  S.resp = Growth 
  Sj.resp = Growth 
  Gg_WinBioE = Growth
  Gg_ELR = Growth
  TotalC = rep(0, Input$N.sites)
  
  
  ##Start Looping Through Time - for Known Consumption, solving for Weight
  
  t=1
  for (t in 1:(Input$N.steps)) {
    
    ### Consumption 
    if(Constants$Consumption$ConsEQ == 1) {
      Cons = with(Constants$Consumption, ConsumptionEQ1(W[t,],Input$Temps[t,],Input$pp[t,],Input$ration[t,],Input$prey.energy.density[t,], Constants$Consumption$CA, Constants$Consumption$CB, Constants$Consumption$CQ))
    } else if(Constants$Consumption$ConsEQ == 2){
      Cons = with(Constants$Consumption, ConsumptionEQ2(W[t,],Input$Temps[t,],Input$pp[t,],Input$ration[t,],Input$prey.energy.density[t,], Constants$Consumption$CA, Constants$Consumption$CB, Constants$Consumption$CTM, Constants$Consumption$CTO, Constants$Consumption$CQ))
    } else if(Constants$Consumption$ConsEQ == 3){
      Cons = with(Constants$Consumption, ConsumptionEQ3(W[t,],Input$Temps[t,],Input$pp[t,],Input$ration[t,],Input$prey.energy.density[t,], Constants$Consumption$CA, Constants$Consumption$CB, Constants$Consumption$CK1, Constants$Consumption$CTO, Constants$Consumption$CQ, Constants$Consumption$CK4, Constants$Consumption$CTL, Constants$Consumption$CTM))
    }
    
    TotalC = TotalC + Cons$CONS * W[t,]
    
    # store daily consumption 
    Consumpt[t,] = as.numeric(Cons$CONS)
    Consumpt_j[t,] = as.numeric(Cons$CONSj)
    
    
    ### Excretion / Egestion
    if(Constants$Excretion$ExcrEQ == 1) {
      ExcEgest<- with(Constants$Excretion, ExcretionEQ1(Cons$CONS, Cons$CONSj, Constants$Excretion$FA, Constants$Excretion$UA))
    } else if (Constants$Excretion$ExcrEQ == 2) {
      ExcEgest<- with(Constants$Excretion, ExcretionEQ2(Cons$CONS, Cons$CONSj, Input$Temps[t,], Input$pp[t,], Input$ration[t,], Cons$CMAX, Constants$Excretion$FA, Constants$Excretion$UA, Constants$Excretion$FB, Constants$Excretion$FG, Constants$Excretion$UB, Constants$Excretion$UG ))
    } else if (Constants$Excretion$ExcrEQ == 3) {
      ExcEgest<- with(Constants$Excretion, ExcretionEQ3(Cons$CONS, Cons$CONSj, Input$Temps[t,], Input$pp[t,], Input$ration[t,], Cons$CMAX, Constants$Excretion$FA, Constants$Excretion$UA, Constants$Excretion$FB, Constants$Excretion$FG, Constants$Excretion$UB, Constants$Excretion$UG, Input$PFF) )
    }
    
    # store daily excretion and egestion
    Excret[t,] = as.numeric(ExcEgest$U)
    Excret_j[t,] = as.numeric(ExcEgest$Uj)
    Egest[t,] = as.numeric(ExcEgest$EG)
    Egest_j[t,] = as.numeric(ExcEgest$EGj)
    
    
    ### Respiration
    if(Constants$Respiration$RespEQ == 1) {
      Resp<- with(Constants$Respiration, RespirationEQ1(W[t,], Input$Temps[t,], Cons$CONS, ExcEgest$EG, Input$prey.energy.density[t,], Input$Oxygen, Constants$Respiration$RA, Constants$Respiration$RB, Constants$Respiration$ACT, Constants$Respiration$SDA, Constants$Respiration$RQ,Constants$Respiration$RTO, Constants$Respiration$RK1, Constants$Respiration$RK4, Constants$Respiration$RTL, Constants$Respiration$BACT))
    } else if (Constants$Respiration$RespEQ == 2) {
      Resp<- with(Constants$Respiration, RespirationEQ2(W[t,], Input$Temps[t,], Cons$CONS, ExcEgest$EG, Input$prey.energy.density[t,], Input$Oxygen, Constants$Respiration$RA, Constants$Respiration$RB, Constants$Respiration$ACT, Constants$Respiration$SDA, Constants$Respiration$RTM, Constants$Respiration$RTO, Constants$Respiration$RQ))
    }
    
    #store daily respiration results
    Respirat[t,] = as.numeric(Resp$R)
    Respirat_j[t,] = as.numeric(Resp$Rj)
    S.resp[t,] = as.numeric(Resp$S)
    Sj.resp[t,] = as.numeric(Resp$Sj)
    
    
    ### Now calculate Growth
    
    # growth in J/g/d - Joules allocated to growth for each gram of predator on each day
    Gj = Cons$CONSj - Resp$Rj - ExcEgest$EGj - ExcEgest$Uj - Resp$Sj	
    G = Cons$CONS - Resp$R - ExcEgest$EG - ExcEgest$U - Resp$S
    # growth in g/d - Grams of predator growth each day
    Growth[t,] = as.numeric(Gj * W[t,]) / pred
    Growth_j[t,] = as.numeric(Gj)
    
    # growth in g/g/d (DailyWeightIncrement divided by fishWeight)
    Gg_WinBioE[t,] = as.numeric(Growth[t,] / W[t,])			
    
    # growth in g/g/d (DailyWeightIncrement divided by average of fish start end weights)
    Gg_ELR[t,] = Growth[t,] / ((as.numeric(Input$Wstart[1:ncol(W)]) + W[t,]) / 2)		
    
    # Calculate absolute weight at time t+1
    W[t + 1,] = W[t,] + Growth[t,]
    
    
  } # End of cycles through time
  
  # Return Results	
  return(list(
    "TotalC" = TotalC,
    "W" = W, 
    "Growth" = Growth, 
    "Gg_WinBioE" = Gg_WinBioE, 
    "Gg_ELR" = Gg_ELR,
    "Growth_j" = Growth_j,
    "Consumption" = Consumpt,
    "Consumption_j" = Consumpt_j,
    "Excretion" = Excret,
    "Excretion_j" = Excret_j,
    "Egestion" = Egest,
    "Egestion_j" = Egest_j,
    "Respiration" = Respirat, 
    "Respiration_j" = Respirat_j, 
    "S.resp" = S.resp,
    "Sj.resp" = Sj.resp,
    "CMAX" = Cons$CMAX
  ))
}

# The main bioenergetics function that calls other functions
BioE<- function (Input, Constants) {

# for simulation method =1 (we have p-vals, and want to solve for weights
if (Input$SimMethod == 1) {
  
	# Method 1: Calculate Growth from p-values and Temperatures
	Results = CalculateGrowth(Constants, Input)
	W = Results$W
	Growth = Results$Growth
	
} else{

# We don't know p-values, but need to iteratively solve for them
# Method 2 or 3:  Calculte P-values from Total Growth or Consumption
# Need to assume p-values are constant with time
	
	# initialize first guess p-values of .5
	Input$pp = array(rep(0.5, Input$N.sites * Input$N.step), 
		 c(Input$N.step, Input$N.sites))

	# set error at high value, iterate until it's small
	Error = rep(99, Input$N.sites)
	iteration = 0

### Interate until error is less than .1
	while(max(abs(Error)) > Input$epsilon) 
{
	iteration = iteration + 1
		Results = CalculateGrowth(Constants, Input, W)
		W = Results$W
    TConsumption = Results$TotalC

 # Find error (depending on which thing on which we're converging), and
 # and come up with new estimate for average p-value
		if (Input$SimMethod == 2)  {
		      Error = (W[Input$N.step + 1,] - Input$Endweights)
  # Delta is the amount by which we'll change the p-value (prior to
  # scaling by the stability factor
		  	Delta = (Input$Endweights) / (W[Input$N.step + 1,]) *
      	         Input$pp[1,] - Input$pp[1,]

		Pnew = as.vector(Input$p[1,] + Input$stab.factor * Delta)
	# Guard against negative p-values (maybe I shouldn't for convergence' sake)
	 for (i in 1:length(Pnew)) {Pnew[i] = max(0, Pnew[i])}
		} else {
	  Error = (TConsumption-Input$TotalConsumption) / Input$TotalConsumption
	  Delta = Input$TotalConsumption/TConsumption * Input$p[1,] - Input$p[1,]
		Pnew = as.vector(Input$p[1,] + Input$stab.factor * Delta)
	}
	# Update for user, show p-values and error, see if we're converging
		for (i in 1:Input$N.step) {Input$p[i,] = as.numeric(Pnew)}
			print(paste("Pnew = ",Pnew, "  Error = ", Error))
	            print(paste("iteration = ", iteration))	
		} # end of while statement 
	} 

if(sum(Input$pp) > 0){Results$pp = Input$pp
} else{ 
  Results$pp = Input$ration / Results$CMAX
  Results$pp[Results$pp > Results$CMAX] = 1 #can't be higher than 1; ration is capped at CMax in consumption equations.
}

return(Results)
}

#=== MORTALITY FUNCTIONS =================================================================
  #fish die based on their size (cumulative experience) and instantaneous growth rate (measure of current conditions)
    #smaller fish with lower growth rates more likely to die (sampled probabilistically)

fncSurvive<- function(df, minprob = 0.98, maxprob = 1, b = 1){
  #df: data frame of fish table with only the survivors, e.g., fish[fish$survive == 1, c("weight","growth")]
  #minprob: smallest probability any fish can have of dying in any time step (0.99*0.99 = 0.98 per day because 2 time steps)
  #b: beta, which changes the shape of the curve (b<1 flattens, b>1 sharpens the bend)

  #Weights
  w = df$weight #weight of fish that are alive at this time step
  
  #Function: 
  v = minprob + (maxprob - minprob) * (1 - 1 / exp(b * w))
  
  #Growth during this time step that reflects recent conditions (i.e., a hungry/stressed fish may behave in ways that make it more vulnerable to predation, etc.)
  g = df$growth 
  if(length(g) > 1) g = fncRescale(g, to = c(-0.001, 0.001)) else g = 0
  
  #Probability of survival
  prb.srv = v + g
  prb.srv[prb.srv > 1] = 1 #set upper bound at 1
  
  #Sample from binomial distribution with probabilities of prb.srv to determine which fish survive this time step
  survivors = rbinom(n = nrow(df), size = 1, prob = prb.srv)
  
  return(survivors)
}


# Fish are predators if:
# * they are alive and above a critical temperature threshold
# Fish are potential prey if:
# * they are alive and emerged
# * they are in the same segment as a predator
# * they are within a predator's movement distance
# * their weight <= (weight threshold) * (predator weight)
#   i.e. if weight threshold = 0.5,
#        predators can only eat prey half their size or smaller
# Whether potential prey is eaten is based on random sampling from a binomial 
# distribution, with a chance of success of predation.probability
# The maximum probability of predation per fish is set by max.pred.prob, and
# probability of predation increases linearly as temperature increases above pred.temp.crit

# Input:
#   pred.survive   = integer, predator survival (1 is alive, 0 is dead)
#   pred.temp      = numeric, predator water temperature
#   pred.seg       = integer, predator segment identificer
#   pred.dist      = numeric, distance of predator from base of segment
#   pred.move      = numeric, distance a predator can move to detect prey
#   pred.weight    = numeric, predator weight
#   pred.temp.crit = numeric, minimum temperature threshold for predators to be active
#   pred.mass.crit = numeric, ratio of predator to prey weights
#   max.pred.prob  = numeric, maximum probability of predator success
#   prey.id        = integer vector, prey identifier
#   prey.survive   = integer vector, prey survival (1 is alive, 0 is dead)
#   prey.emerge    = integer vector, prey emergence (1 is emerged, 0 is not)
#   prey.seg       = integer vector, prey segment identifier
#   prey.weight    = numeric vector, prey weight

# Output (list of four elements)
#   number.potential prey = number of potential prey for this predator
#     NA represents predation is not possible for this predator (dead or too cold)
#     0 represents eligible predator, but no potential prey
#     integer represents number of potential prey
#   predation.probability = numeric, probability of successful predation for this predator
#     NA represents predation is not possible for this predator (dead or too cold)
#     0 to 1 represents likelihood of successful predation (higher is more likely to be successful)
#   number.prey.eaten = number of potential prey that are eaten by this predator
#     NA represents predation is not possible for this predator (dead or too cold)
#     integer represents number of potential prey eaten by this predator
#   prey.eaten.id = integer or integer vector, length is number.prey.eaten
#     NA represents predation is not possible for this predator (dead or too cold)
#     integer(s) represent identifiers for prey eaten by this predator

fncOnePredator <- function(pred.survive, pred.temp, pred.seg, pred.dist,
                           pred.move, pred.weight, pred.temp.crit, 
                           pred.mass.crit, max.pred.prob, prey.id, prey.survive,
                           prey.emerge, prey.seg, prey.dist, prey.weight) {
  
  # initialize output
  output <- list(number.potential.prey = NA,
                 predation.probability = NA,
                 number.prey.eaten = NA,
                 prey.eaten.id = NA)
  
  # check if predator is alive and warm enough to eat
  if ((pred.survive == 1) & (pred.temp >= pred.temp.crit)) {
    
    # calculate how far predator can move within segment
    dist.max = pred.dist + pred.move
    dist.min = pred.dist - pred.move
    
    # check if there are potential prey nearby
    potential.prey.index = which((prey.survive == 1) & # prey are alive
                                   (prey.emerge == 1) & # prey are emerged
                                   (prey.seg == pred.seg) & # prey are in same segment
                                   (prey.dist >= dist.min) & (prey.dist <= dist.max) &  # prey are within movement distance
                                   (prey.weight <= pred.mass.crit * pred.weight)) # prey are smaller than predator by some ratio
    number.potential.prey <- length(potential.prey.index)
    # update output
    output$number.potential.prey <- number.potential.prey
    
    if (number.potential.prey > 0) {
      # calculate probability
      predation.probability <- max.pred.prob * ((pred.temp - pred.temp.crit) / pred.temp.crit)
      # set upper bound at max.pred.prob
      if (predation.probability > max.pred.prob) predation.probability <- max.pred.prob
      # update output
      output$predation.probability <- predation.probability
      
      # sample from binomial distribution to determine which potential prey are eaten
      prey.eaten <- rbinom(n = number.potential.prey, size = 1, prob = predation.probability)
      # update output
      output$number.prey.eaten <- number.prey.eaten <- sum(prey.eaten)
      if (number.prey.eaten > 0) {
        prey.eaten.index <- potential.prey.index[which(prey.eaten == 1)]
        output$prey.eaten.id <- prey.id[prey.eaten.index]
      }
    }
  }
  return(output)
}


#=== PLOTTING FUNCTIONS ========================================================

#plot.SpatialStreamNetwork modified to allow a subset of the predpoints to be passed in
plotSSN.mod<- function (x, VariableName = NULL, color.palette = NULL, nclasses = NULL, 
    breaktype = "quantile", brks = NULL, PredPointsID = NULL, add = FALSE, addWithLegend = FALSE, 
    lwdLineCol = NULL, lwdLineEx = 1, lineCol = "black", myvar = NULL,...) {
  if (missing(lwdLineEx)) 
    lwdLineEx = 1
  if (missing(lwdLineCol)) {
    x@data$lineWidth = rep(1, nrow(x@data))
    lwdLineCol = "lineWidth"
  }
  if (is.null(as.list(match.call()[-1])$pch)) {
    plch = 19
  }
  else plch = as.list(match.call()[-1])$pch
  if (is.null(as.list(match.call()[-1])$cex)) {
    chex = 1
  }
  else chex = as.list(match.call()[-1])$cex
  if (is.null(as.list(match.call()[-1])$col)) {
    colr = "black"
  }
  else colr = as.list(match.call()[-1])$col
  par.orig = par(no.readonly = TRUE)
  if (!is.null(PredPointsID)) {
    for (i in 1:length(x@predpoints@ID)) {
      if (x@predpoints@ID[i] == PredPointsID) {
        if (add == FALSE & addWithLegend == FALSE) {
          plot(x@bbox[1, ], x@bbox[2, ], type = "n", ...)
          for (j in 1:length(x@lines)) for (k in 1:length(x@lines[[j]])) if (is.null(lwdLineCol)) 
            lines((x@lines[[j]]@Lines[[k]]@coords), col = lineCol, ...)
          else lines(x@lines[[j]]@Lines[[k]]@coords, lwd = lwdLineEx * x@data[i, lwdLineCol], col = lineCol, ...)
        }
        if (add == TRUE) {
          par(new = TRUE)
          plot(x@bbox[1, ], x@bbox[2, ], type = "n", bty = "n", xlab = "", ylab = "", ...)
        }
        if (addWithLegend == TRUE) {
          par(new = TRUE)
          layout(matrix(1:2, nrow = 1), widths = c(4, 1))
          par(mar = c(5, 5, 3, 0))
          par(mfg = c(1, 1))
          plot(x@bbox[1, ], x@bbox[2, ], type = "n", bty = "n", xlab = "", ylab = "", ...)
        }
        if(!is.null(myvar)){
          points(x@predpoints@SSNPoints[[2]]@point.coords[myvar,1],
                 x@predpoints@SSNPoints[[2]]@point.coords[myvar,2],
                 pch = plch, cex = chex, col = colr)
        } else{
          points(x@predpoints@SSNPoints[[i]]@point.coords, 
                 pch = plch, cex = chex, col = colr)
        }
      }
    }
    par(par.orig)
  }
  else if (is.null(VariableName)) {
    plot(x@bbox[1, ], x@bbox[2, ], type = "n", ...)
    for (i in 1:length(x@lines)) for (j in 1:length(x@lines[[i]])) if (is.null(lwdLineCol)) 
      lines((x@lines[[i]]@Lines[[j]]@coords), col = lineCol, ...)
    else lines(x@lines[[i]]@Lines[[j]]@coords, lwd = lwdLineEx * x@data[i, lwdLineCol], col = lineCol, ...)
    points(x@obspoints@SSNPoints[[1]]@point.coords, pch = plch, cex = chex, col = colr)
    par(par.orig)
  }
  else {
    layout(matrix(1:2, nrow = 1), widths = c(4, 1))
    par(mar = c(5, 5, 3, 0))
    plot(x@bbox[1, ], x@bbox[2, ], type = "n", ...)
    for (i in 1:length(x@lines)) for (j in 1:length(x@lines[[i]])) if (is.null(lwdLineCol)) 
      lines((x@lines[[i]]@Lines[[j]]@coords), col = lineCol, ...)
    else lines(x@lines[[i]]@Lines[[j]]@coords, lwd = lwdLineEx * x@data[i, lwdLineCol], col = lineCol, ...)
    data = x@obspoints@SSNPoints[[1]]@point.data
    if (is.null(nclasses)) 
      nclasses = 10
    lower.breaks = matrix(0, nrow = nclasses, ncol = 1)
    upper.breaks = matrix(0, nrow = nclasses, ncol = 1)
    if (breaktype == "quantile") {
      brks = quantile(data[, VariableName], probs = (1:(nclasses - 1))/nclasses, na.rm = T)
      lower.breaks = c(min(data[, VariableName], na.rm = T), brks)
      upper.breaks = c(brks, max(data[, VariableName], na.rm = T))
    }
    if (breaktype == "even") {
      brks = min(data[, VariableName]) + (max(data[, VariableName]) - min(data[, VariableName])) * (1:(nclasses - 1))/nclasses
      lower.breaks = c(min(data[, VariableName], na.rm = T), brks)
      upper.breaks = c(brks, max(data[, VariableName], na.rm = T))
    }
    if (breaktype == "user") {
      if (is.null(brks)) 
        return("Must specify brks if breaktype = user")
      minD = min(data[, VariableName], na.rm = TRUE)
      maxD = max(data[, VariableName], na.rm = TRUE)
      brks = as.vector(unlist(brks))
      if (minD < min(brks)) 
        brks = c(brks, minD)
      if (maxD > max(brks)) 
        brks = c(brks, maxD)
      brks = sort(unique(unlist(brks)))
      nclasses = length(brks) - 1
      lower.breaks = brks[1:nclasses]
      upper.breaks = brks[2:(nclasses + 1)]
    }
    if (length(color.palette) == 0) 
      color.palette = rainbow(nclasses, start = 0.66, end = 0.99)
    for (j in 1:nclasses) {
      jmax = upper.breaks[j]
      jmin = lower.breaks[j]
      indj = data[, VariableName] >= jmin & data[, VariableName] <= jmax
      points(x@obspoints@SSNPoints[[1]]@point.coords[indj, , drop = F], col = color.palette[j], pch = plch, cex = chex)
    }
    dec.dig = 2
    left = as.character(as.numeric(as.integer(lower.breaks * 10^dec.dig))/10^dec.dig)
    rght = as.character(as.numeric(as.integer(upper.breaks * 10^dec.dig))/10^dec.dig)
    leglabs = paste(left, "to", rght)
    par(mar = c(0, 0, 0, 0))
    plot(c(0, 0), c(1, 1), type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n")
    legend(x = -1, y = 1.1, legend = leglabs, bty = "n", 
           pch = rep(plch, times = length(leglabs)), col = color.palette, cex = 0.8)
    par(par.orig)
    return(invisible(data.frame(lower.breaks = lower.breaks, upper.breaks = upper.breaks)))
  }
}

#Hack to allow up to 14 colors in 'Spectral' palette instead of 11
fncColBrewPlus<- function(n, paint = FALSE){
  cb = switch(n -2, 
               rgb(c(252, 255, 153), c(141, 255, 213), c(89,191, 148), maxColorValue = 255), 
               rgb(c(215, 253, 171, 43), c(25, 174, 221, 131), c(28, 97, 164, 186), maxColorValue = 255), 
               rgb(c(215, 253, 255, 171, 43),c(25, 174, 255, 221, 131), c(28, 97, 191, 164, 186), maxColorValue = 255), 
               rgb(c(213, 252, 254, 230, 153, 50), c(62, 141, 224, 245, 213, 136), c(79, 89, 139,152, 148, 189), maxColorValue = 255), 
               rgb(c(213, 252, 254, 255, 230, 153, 50), c(62, 141, 224, 255, 245, 213, 136), c(79, 89, 139, 191, 152, 148, 189), maxColorValue = 255), 
               rgb(c(213, 244, 253, 254, 230, 171, 102, 50), c(62, 109, 174, 224, 245, 221, 194,136), c(79, 67, 97, 139, 152, 164, 165, 189), maxColorValue = 255),
               rgb(c(213, 244, 253, 254, 255, 230, 171, 102, 50),c(62, 109, 174, 224, 255, 245, 221, 194, 136),c(79, 67, 97, 139, 191, 152, 164, 165, 189), maxColorValue = 255),
               rgb(c(158, 213, 244, 253, 254, 230, 171, 102, 50, 94), c(1, 62, 109, 174, 224, 245, 221, 194, 136, 79), c(66, 79, 67, 97,139, 152, 164, 165, 189, 162), maxColorValue = 255),
               rgb(c(158, 213, 244, 253, 254, 255, 230, 171, 102, 50, 94), c(1, 62, 109, 174, 224, 255, 245, 221, 194, 136, 79), c(66, 79, 67, 97, 139, 191, 152, 164, 165, 189, 162), maxColorValue = 255),
               rgb(c(158, 213, 244, 253, 252, 252, 222, 173, 118, 102, 50, 94), c(1, 62, 109, 174, 223, 247, 245, 247, 227, 194, 136, 79), c(66, 79, 67, 97, 91, 91, 100, 111, 131, 165, 189, 162), maxColorValue = 255),
               rgb(c(158, 213, 244, 253, 252, 252, 222, 173, 118, 102, 50, 94, 72), c(1, 62, 109, 174, 223, 247, 245, 247, 227, 194, 136, 79, 22), c(66, 79, 67, 97, 91, 91, 100, 111, 131, 165, 189, 162, 138), maxColorValue = 255),
               rgb(c(158, 213, 244, 253, 252, 252, 222, 173, 118, 102, 55, 50, 94, 72), c(1, 62, 109, 174, 223, 247, 245, 247, 227, 194, 172, 136, 79, 22), c(66, 79, 67, 97, 91, 91, 100, 111, 131, 165, 204, 189, 162, 138), maxColorValue = 255)
  )
  cb = cb[n:1] 
  if(paint == TRUE) image(1:n, 1, as.matrix(1:n), col = cb, axes = FALSE, ylab = "", xlab = "")
  return(cb)
}

#Title to print at top of image that conveys date and time
fncGetTitle<- function(d, t){
  #only 2 daily time steps in data (6am and 6pm)
  if(t == 6 | t == 11){thetitle = paste0(format(d, "%d %B %Y"),": 06:00")}
  if(t == 18 | t == 35){thetitle = paste0(format(d, "%d %B %Y"),": 18:00")}
  return(thetitle)
} 

#=== SENSITIVITY ANALYSIS FUNCTIONS ============================================

# Parameter Set for Global Sensisivity Analyses:
fncGetParametersSA<- function(vary_nFish = FALSE){
  set.seed(iter)
  
if(SecondSpecies == FALSE){ # salmon alone scenario
  # vary these parameters by drawing from a normal distribution
  salmon.parm.list = c("nSpawnDays", "spawn.date.shape","ATU.crit", "mvdist.shape","mvmt.scalar",
                       "ration.hi","max.density","prey.en.dens","pred.en.dens","om.mass","survival.shape")
  for(parm in salmon.parm.list){parameters["salmon",parm] = rnorm(n = 1, mean = parameters["salmon",parm], sd = 0.1 * parameters["salmon",parm])}
  
  # special cases that need a specific range to be feasible (survival.min needs to have lower minimum and survival.max needs a narrower distribution)
  parameters["salmon","survival.min"] = rnorm(n = 1, mean = parameters["salmon","survival.min"] * 0.975, sd = 0.01 * parameters["salmon","survival.min"])*1.01
  if(parameters["salmon","survival.min"] > 1) parameters["salmon","survival.min"] = 1 #can't be >1
  parameters["salmon","survival.max"] = rnorm(n = 1, mean = parameters["salmon","survival.max"], sd = 0.003 * parameters["salmon","survival.max"])
  if(parameters["salmon","survival.max"] > 1) parameters["salmon","survival.max"] = 1 #can't be >1
  # ensure that min < max!
  while(parameters["salmon","survival.min"] > parameters["salmon","survival.max"]){
    parameters["salmon","survival.min"] = rnorm(n = 1, mean = parameters["salmon","survival.min"] * 0.975, sd = 0.01 * parameters["salmon","survival.min"])*1.01
    if(parameters["salmon","survival.min"] > 1) parameters["salmon","survival.min"] = 1 #can't be >1
    parameters["salmon","survival.max"] = rnorm(n = 1, mean = parameters["salmon","survival.max"], sd = 0.003 * parameters["salmon","survival.max"])
    if(parameters["salmon","survival.max"] > 1) parameters["salmon","survival.max"] = 1 #can't be >1
  }
  
  # round to whole days
  parameters[,"nSpawnDays"] = round(parameters[,"nSpawnDays"],0) 
  
  # dates, gives a range of about a month
  # special handling for 'spawn.date.begin' because the data time series begins on September 1 so can't go earlier
  for(parm in c("spawn.date.begin","om.date.taper","om.date.end")){
    date.numeric = as.numeric(as.Date(paste0(climate.year-1, "-", parameters["salmon",parm])))
    if(parm == "spawn.date.begin") date.numeric = date.numeric + 10
    thedate = as.Date(round(rnorm(n = 1, mean = date.numeric, sd = 0.00025 * date.numeric),0), origin = "1970-01-01")
    if(parm == "spawn.date.begin" & thedate < as.Date(paste0(climate.year-1, "-", "09-01"))){
      parameters["salmon",parm] = "09-01"
    } else {
    parameters["salmon",parm] = substr(as.character(thedate),6,10)
    }
  }
  
  # vary initial number of fish from pre-determined set of 5 sizes (fish shapefiles had to be altered prior to this)
  if(vary_nFish == TRUE){
    parameters["salmon","nFish"] = floor(parameters["salmon","nFish"] * sample(c(0.5, 0.8, 1, 1.2, 1.5), 1)) #sample(c(2241, 1120, 1792, 2689, 3361), 1)
    parameters["salmon", "fish.shp"] = paste0("chin_spawn_", parameters["salmon","nFish"])
  }
  
} else { # second species scenario
  
    # Parameters directly affecting salmon (omitting 7 that were insensitive in 'alone' scenario to save processing time):
  
    # vary these parameters by drawing from a normal distribution
    salmon.parm.list = c("ATU.crit", "ration.hi","prey.en.dens","pred.en.dens","om.mass")
    for(parm in salmon.parm.list){parameters["salmon",parm] = rnorm(n = 1, mean = parameters["salmon",parm], sd = 0.1 * parameters["salmon",parm])}
    
    # special cases that need a specific range to be feasible (survival.min needs to have lower minimum and survival.max needs a narrower distribution)
    parameters["salmon","survival.min"] = rnorm(n = 1, mean = parameters["salmon","survival.min"] * 0.975, sd = 0.01 * parameters["salmon","survival.min"])*1.01
    if(parameters["salmon","survival.min"] > 1) parameters["salmon","survival.min"] = 1 #can't be >1
    parameters["salmon","survival.max"] = rnorm(n = 1, mean = parameters["salmon","survival.max"], sd = 0.003 * parameters["salmon","survival.max"])
    if(parameters["salmon","survival.max"] > 1) parameters["salmon","survival.max"] = 1 #can't be >1
    # ensure that min < max!
    while(parameters["salmon","survival.min"] > parameters["salmon","survival.max"]){
      parameters["salmon","survival.min"] = rnorm(n = 1, mean = parameters["salmon","survival.min"] * 0.975, sd = 0.01 * parameters["salmon","survival.min"])*1.01
      if(parameters["salmon","survival.min"] > 1) parameters["salmon","survival.min"] = 1 #can't be >1
      parameters["salmon","survival.max"] = rnorm(n = 1, mean = parameters["salmon","survival.max"], sd = 0.003 * parameters["salmon","survival.max"])
      if(parameters["salmon","survival.max"] > 1) parameters["salmon","survival.max"] = 1 #can't be >1
    }

    # dates, gives a range of about a month
    # special handling for 'spawn.date.begin' because the data time series begins on September 1 so can't go earlier
    for(parm in c("spawn.date.begin","om.date.end")){
      date.numeric = as.numeric(as.Date(paste0(climate.year-1, "-", parameters["salmon",parm])))
      if(parm == "spawn.date.begin") date.numeric = date.numeric + 10
      thedate = as.Date(round(rnorm(n = 1, mean = date.numeric, sd = 0.00025 * date.numeric),0), origin = "1970-01-01")
      if(parm == "spawn.date.begin" & thedate < as.Date(paste0(climate.year-1, "-", "09-01"))){
        parameters["salmon",parm] = "09-01"
      } else {
        parameters["salmon",parm] = substr(as.character(thedate),6,10)
      }
    }

    # vary initial number of fish from pre-determined set of 5 sizes (fish shapefiles had to be altered prior to this)
    if(vary_nFish == TRUE){
      parameters["salmon","nFish"] = floor(parameters["salmon","nFish"] * sample(c(0.5, 0.8, 1, 1.2, 1.5), 1)) #sample(c(2241, 1120, 1792, 2689, 3361), 1)
      parameters["salmon", "fish.shp"] = paste0("chin_spawn_", parameters["salmon","nFish"])
    }
    
    # Parameters directly affecting fish_other:
    
    # vary these parameters by drawing from a normal distribution
    fish_other.parm.list = c("max.initial.mass","mass.shape","mvdist.shape","mvmt.scalar",
                             "ration.hi","max.density","prey.en.dens","pred.en.dens","survival.shape",
                             "pred.mass.crit","max.pred.prob","pred.temp.crit", "pred.move")
    for(parm in fish_other.parm.list){parameters["fish_other",parm] = rnorm(n = 1, mean = parameters["fish_other",parm], sd = 0.1 * parameters["fish_other",parm])}
    
    # special cases that need a specific range to be feasible (survival.min needs to have lower minimum)
    parameters["fish_other","survival.min"] = rnorm(n = 1, mean = parameters["fish_other","survival.min"] * 0.975, sd = 0.01 * parameters["fish_other","survival.min"])*1.01
    if(parameters["fish_other","survival.min"] > 1) parameters["fish_other","survival.min"] = 1 #can't be >1
    
    # ensure 'max.pred.prob' is constrained between 0 and 1
    parameters["fish_other","max.pred.prob"][parameters["fish_other","max.pred.prob"] > 1] = 1
    parameters["fish_other","max.pred.prob"][parameters["fish_other","max.pred.prob"] < 0] = 0
    
    # vary initial number of fish from pre-determined set of 5 sizes (fish shapefiles had to be altered prior to this)
    if(vary_nFish == TRUE){
      parameters["fish_other","nFish"] = floor(parameters["fish_other","nFish"] * sample(c(0.5, 0.8, 1, 1.2, 1.5), 1)) #sample(c(500, 250, 400, 600, 750), 1)
      parameters["fish_other", "fish.shp"] = paste0("lmb_prevalent_", parameters["fish_other","nFish"])
    }
  }
  
  return(parameters)
}

#=== END OF FILE ===============================================================
## Global Sensitivity Analysis

# Steps
# 1. Run the model with the global variable 'SA' set to TRUE and iter.list set to '1:300'
# 2. Use this script to:
#   a) Collate parameters sets into one matrix
#   b) Collate simulated data output metrics
#   c) Use random forest regression to evaluate parameter sensitivity
#   d) Produce plots to visualize sensitivity analysis results, and
#   e) Replicate Figure 8 and all figures in Appendix 2

#----- SETUP -------------------------------------------------------------------

library(abind)
library(randomForest)
library(randomForestExplainer)

species = "salmon"
climate.year = 2014
iter.list = 1:300
missing = c() # include any runs that were omitted
iter.list = setdiff(iter.list, missing)
nFish.list = c(2241, 1120, 1792, 2689, 3361) # the number of salmon per simulation

#----- COLLATE PARAMETERS ------------------------------------------------------

# Collate parameters from 'run.info.txt' files for Baseline Scenario.

# specify directories
dataDir = "data.out/A.2014"
outputDir = paste0(dataDir, "/SA")
if (! dir.exists(outputDir)) {dir.create(outputDir)}

# specify salmon parameters
parameters = c("runtime (h)", "seed", "nFish", "mvdist.shape", "ration.lo", "survival.shape", "spawn.date.begin",
               "om.date.end", "nSpawnDays", "spawn.date.shape", "ATU.crit", "mvmt.scalar", "om.mass","om.date.taper",
               "pred.en.dens", "prey.en.dens", "max.density", "ration.hi", "survival.min", "survival.max")
parms.date = c("spawn.date.begin", "om.date.taper", "om.date.end")
parms.num = setdiff(parameters, c("runtime (h)", parms.date))

# read parameters from run.info text file from first replicate
td = read.csv(paste0(dataDir, "/run.info.", climate.year, ".", 1, ".txt"), sep = ":", stringsAsFactors = F)
td = td[which(td[ , 1] %in% parameters), ]
run.summary = td[order(td[ , 1]), ]

# read parameters from run.info text files from remaining replicates
for(iter in iter.list[-1]){
  td = read.csv(paste0(dataDir, "/run.info.", climate.year, ".", iter, ".txt"), sep = ":", stringsAsFactors = F)
  td = td[which(td[ , 1] %in% parameters),]
  run.summary = cbind(run.summary, td[order(td[ , 1]), 2], stringsAsFactors = F)
}

# reshape parameters stored in run.summary
cnams = run.summary[ , 1]
run.summary = as.data.frame(t(run.summary[ , -1]), stringsAsFactors = F)
colnames(run.summary) = cnams
rownames(run.summary) = NULL
for(i in parms.num) run.summary[ , i] = as.numeric(run.summary[,i])
for(i in parms.date){ run.summary[ , i] = as.Date(paste0(climate.year, "-", run.summary[ , i]), origin = "1970")}
colnames(run.summary)[colnames(run.summary) == "runtime (h)"] = "runtime"

# sort parameters in run.summary based on number of fish to match the simulated data
run.summary = run.summary[order(run.summary$nFish, decreasing = T), ]

# save parameters in run.summary
write.csv(run.summary, paste0(outputDir, "/run.info.SA_all.csv"), row.names = FALSE)


# Collate parameters from 'run.info.txt' files for Predator Scenario.

# specify directories
dataDir = "data.out/P.2014"
outputDir = paste0(dataDir, "/SA")
if (! dir.exists(outputDir)) {dir.create(outputDir)}

# specify salmon and bass parameters
parameters = c("runtime (h)", "seed", "nFish", "spawn.date.begin", "om.date.end", "ATU.crit", "om.mass",
               "pred.en.dens", "prey.en.dens", "ration.hi", "survival.min", "survival.max")
parameters.bass = c("nFish", "max.initial.mass", "mass.shape", "mvdist.shape", "mvmt.scalar", 
                    "ration.hi", "max.density", "prey.en.dens", "pred.en.dens", "survival.min", 
                    "survival.shape", "max.pred.prob", "pred.temp.crit", "pred.mass.crit", "pred.move")
parameters.bass = paste0(parameters.bass, ".b")
parameters = c(parameters, parameters.bass)

parms.date = c("spawn.date.begin", "om.date.end", "om.date.taper")
parms.num = setdiff(parameters, c("runtime (h)",parms.date))

# read parameters from run.info text file from first replicate
td = read.csv(paste0(dataDir, "/run.info.", climate.year, ".", 1, ".txt"), sep = ":", stringsAsFactors = F)
td[45:75, 1] = paste0(td[45:75, 1], ".b")
td = td[which(td[ , 1] %in% parameters), ]
run.summary = td
row.names(run.summary) = NULL

# read parameters from run.info text files from remaining replicates
for(iter in iter.list[-1]){
  td = read.csv(paste0(dataDir, "/run.info.", climate.year, ".", iter, ".txt"), sep = ":", stringsAsFactors = F)
  td[45:75, 1] = paste0(td[45:75, 1], ".b")
  td = td[which(td[ , 1] %in% parameters),]
  run.summary = cbind(run.summary,td[ , 2], stringsAsFactors = F)
}

# reshape parameters stored in run.summary
cnams = run.summary[ , 1]
run.summary = as.data.frame(t(run.summary[ , -1]), stringsAsFactors = F)
colnames(run.summary) = cnams
rownames(run.summary) = NULL
for(i in parms.num) run.summary[ , i] = as.numeric(run.summary[ , i])
for(i in parms.date){ run.summary[ , i] = as.Date(paste0(climate.year, "-", run.summary[ , i]), origin = "1970")}
colnames(run.summary)[colnames(run.summary) == "runtime (h)"] = "runtime"
run.summary$runtime = as.numeric(run.summary$runtime)

# sort parameters in run.summary based on number of fish to match the simulated data
run.summary = run.summary[order(run.summary$nFish, decreasing = T),]

# save parameters in run.summary
write.csv(run.summary, paste0(outputDir, "/run.info.SA_all.csv"), row.names = FALSE)


#----- COLLATE SIMULATED DATA --------------------------------------------------

# Collate simulated data from 'salmon.finalstep.201x.1.RData' files ----
for(tag in c("A", "P")){ # "A" for non-Predator scenarios, "P" for Predator scenarios
  
  # set up directories
  dataDir = paste0("data.out/", tag, ".", climate.year)
  outputDir = paste0(dataDir, "/SA")
  
  # collate arrays with the same number of fish
  # have to do these separately due to different sized arrays
  for(n in nFish.list) {
    for(iter in iter.list) {
      load(paste0(dataDir, "/", species, ".finalstep.", climate.year, ".", iter, ".RData")) 
      if(nrow(salmon.finalstep) == n & exists("fish.finalstep") == FALSE){
        # if this is the first file, initialize array by copying this file
        fish.finalstep = salmon.finalstep
      } else if(nrow(salmon.finalstep) == n & exists("fish.finalstep") == TRUE) {
        # else if this is not the first file, append the new file to the existing array
        fish.finalstep = abind(fish.finalstep, salmon.finalstep, along = 3)
      }
    }
    # rename the array
    salmon.finalstep = fish.finalstep
    # save the array
    save("salmon.finalstep", file = paste0(outputDir, "/", species, ".finalstep.", tag, ".", climate.year, ".SA", n, ".RData"))
    # remove temporary variables
    rm(fish.finalstep, salmon.finalstep)
  }

  # merge ragged arrays
  for(n in nFish.list){
    load(paste0(outputDir,"/", species, ".finalstep.", tag, ".", climate.year, ".SA", n, ".RData"))
    assign(paste0(species, ".finalstep.", n), salmon.finalstep)
    rm(salmon.finalstep)
  }
  
  # save ragged arrays
  fish.merged = array(NA, dim = c(max(nFish.list), 11, length(iter.list)))
  fish.merged[,,1:dim(salmon.finalstep.3361)[3]] = salmon.finalstep.3361
  fish.merged[1:dim(salmon.finalstep.2689)[1], , dim(salmon.finalstep.3361)[3] + 1:dim(salmon.finalstep.2689)[3]] = salmon.finalstep.2689
  fish.merged[1:dim(salmon.finalstep.2241)[1], , dim(salmon.finalstep.3361)[3] + dim(salmon.finalstep.2689)[3] + 1:dim(salmon.finalstep.2241)[3]] = salmon.finalstep.2241
  fish.merged[1:dim(salmon.finalstep.1792)[1], , dim(salmon.finalstep.3361)[3] + dim(salmon.finalstep.2689)[3] + dim(salmon.finalstep.2241)[3] + 1:dim(salmon.finalstep.1792)[3]] = salmon.finalstep.1792
  fish.merged[1:dim(salmon.finalstep.1120)[1], , dim(salmon.finalstep.3361)[3] + dim(salmon.finalstep.2689)[3] + dim(salmon.finalstep.2241)[3] + dim(salmon.finalstep.1792)[3] + 1:dim(salmon.finalstep.1120)[3]] = salmon.finalstep.1120
  colnames(fish.merged) = colnames(salmon.finalstep.3361)
  salmon.finalstep = fish.merged
  salmon.finalstep = salmon.finalstep[ , , 1:nrow(run.summary)]
  save("salmon.finalstep", file = paste0(outputDir, "/", species, ".finalstep.", tag, ".", climate.year, ".SA.all.RData"))
  rm(salmon.finalstep.1120, salmon.finalstep.1792, salmon.finalstep.2241, salmon.finalstep.2689, salmon.finalstep.3361, fish.merged)

}



# Load necessary data
tag = "A" # "A" for non-Predator scenarios, "P" for Predator scenarios

dataDir = paste0("data.out/", tag, ".", climate.year)
outputDir = paste0(dataDir, "/SA")
if (! dir.exists(paste0(outputDir, "/SA.RF"))) {dir.create(paste0(outputDir, "/SA.RF"))}
run.summary = read.csv(paste0(outputDir, "/run.info.SA_all.csv"), header = TRUE, stringsAsFactors = FALSE)
for(i in c("spawn.date.begin", "om.date.end", "om.date.taper")) {
  run.summary[ , i] = as.Date(run.summary[ , i], origin = "1970-01-01")
}
load(paste0(outputDir,"/", species, ".finalstep.", tag,".", climate.year, ".SA.all.RData"))

# Set up data filters
alive = salmon.finalstep[,"survive",]
subyearlings = alive; subyearlings[subyearlings != 2] = NA; subyearlings[subyearlings == 2] = 1
yearlings = alive; yearlings[yearlings != 1] = NA

# Plot distribution of parameters
# Figure S2-1: Baseline scneario
if(tag == "A") {
  png(paste0(outputDir,"/SA.plots/parameter_distributions_A.2014.png"),units="in",height=9,width=12,res=600)
    par(mfrow=c(4,4), cex = 0.7, mar = c(4,2,1,1), oma=c(1,1,1,1), las = 1)
    hist(as.numeric(run.summary$spawn.date.begin), main = "", xlab = "spawn.date.begin", ylab = "")
    hist(run.summary$nSpawnDays, main = "", xlab = "nSpawnDays", ylab = "")
    hist(run.summary$spawn.date.shape, main = "", xlab = "spawn.date.shape", ylab = "")
    hist(run.summary$ATU.crit, main = "", xlab = "ATU.crit", ylab = "")
    hist(run.summary$mvmt.scalar, main = "", xlab = "mvmt.scalar", ylab = "")
    hist(run.summary$mvdist.shape, main = "", xlab = "mvdist.shape", ylab = "")
    hist(run.summary$ration.hi, main = "", xlab = "ration.hi", ylab = "")
    hist(run.summary$pred.en.dens, main = "", xlab = "pred.en.dens", ylab = "")
    hist(run.summary$prey.en.dens, main = "", xlab = "prey.en.dens", ylab = "")
    hist(run.summary$max.density, main = "", xlab = "max.density", ylab = "")
    hist(run.summary$survival.min, main = "", xlab = "survival.min", ylab = "")
    hist(run.summary$survival.max, main = "", xlab = "survival.max", ylab = "")
    hist(run.summary$survival.shape, main = "", xlab = "survival.shape", ylab = "")
    hist(run.summary$om.mass, main = "", xlab = "om.mass", ylab = "")
    hist(as.numeric(run.summary$om.date.taper), main = "", xlab = "om.date.taper", ylab = "")
    hist(as.numeric(run.summary$om.date.end), main = "", xlab = "om.date.end", ylab = "")
  dev.off()
}

# Figure S2-2: Predator scenario
if(tag == "P"){
  png(paste0(outputDir,"/SA.plots/parameter_distributions_P.2014.png"),units="in",height=9,width=12,res=600)
    par(mfrow=c(5,5), cex = 0.7, mar = c(4,2,1,1), oma=c(1,1,1,1), las = 1)
    hist(as.numeric(run.summary$spawn.date.begin), main = "", xlab = "spawn.date.begin", ylab = "")
    hist(run.summary$ATU.crit, main = "", xlab = "ATU.crit", ylab = "")
    hist(run.summary$ration.hi, main = "", xlab = "ration.hi", ylab = "")
    hist(run.summary$pred.en.dens, main = "", xlab = "pred.en.dens", ylab = "")
    hist(run.summary$prey.en.dens, main = "", xlab = "prey.en.dens", ylab = "")
    hist(run.summary$survival.min, main = "", xlab = "survival.min", ylab = "")
    hist(run.summary$survival.max, main = "", xlab = "survival.max", ylab = "")
    hist(run.summary$om.mass, main = "", xlab = "om.mass", ylab = "")
    hist(as.numeric(run.summary$om.date.end), main = "", xlab = "om.date.end", ylab = "")
    plot(1:10,1:10,type = 'n', xlab="", ylab="", axes=F)
    hist(run.summary$max.initial.mass.b, main = "", xlab = "max.initial.mass.b", ylab = "")
    hist(run.summary$mass.shape.b, main = "", xlab = "mass.shape.b", ylab = "")
    hist(run.summary$mvdist.shape.b, main = "", xlab = "mvdist.shape.b", ylab = "")
    hist(run.summary$mvmt.scalar.b, main = "", xlab = "mvmt.scalar.b", ylab = "")
    hist(run.summary$ration.hi.b, main = "", xlab = "ration.hi.b", ylab = "")
    hist(run.summary$max.density.b, main = "", xlab = "max.density.b", ylab = "")
    hist(run.summary$pred.en.dens.b, main = "", xlab = "pred.en.dens.b", ylab = "")
    hist(run.summary$prey.en.dens.b, main = "", xlab = "prey.en.dens.b", ylab = "")
    hist(run.summary$survival.min.b, main = "", xlab = "survival.min.b", ylab = "")
    hist(run.summary$max.pred.prob.b, main = "", xlab = "max.pred.prob.b", ylab = "")
    hist(run.summary$pred.temp.crit.b, main = "", xlab = "pred.temp.crit.b", ylab = "")
    hist(run.summary$pred.mass.crit.b, main = "", xlab = "pred.mass.crit.b", ylab = "")
    hist(run.summary$pred.move.b, main = "", xlab = "pred.move.b", ylab = "")
  dev.off()
}


#---- RANDOM FOREST REGRESSION -------------------------------------------------

# Run random forest regressions to evaluate parameter sensitivity

parms = colnames(run.summary)[!colnames(run.summary) %in% c("runtime", "seed", "ration.lo")]
x = run.summary[ , parms]
colnames(salmon.finalstep)
metric.list = c("weight.yearlings", "weight.subyearlings", "dateOm", "survival", "prop.yearlings", "dateEm")

for (m in metric.list) {
  # filter data depending on metrics
  if (m == "weight.subyearlings") {
    y = apply(salmon.finalstep[ ,"weight", ] * subyearlings, 2, median, na.rm=T)
  } else if (m == "weight.yearlings"){
    y = apply(salmon.finalstep[ ,"weight", ] * yearlings, 2, median, na.rm=T)
  } else if (m == "survival") {
    y = survival
  } else if (m == "prop.yearlings") {
    y = prop.yearlings
  } else {
    y = apply(salmon.finalstep[,m,], 2, median, na.rm=T)
  }
  
  mydat = cbind(y, x)
  na.idx = which(is.na(mydat[,1]))
  if (length(na.idx) > 0) mydat = mydat[-na.idx,]
  set.seed(1)
  num.trees = 5000
  (rf1 = randomForest(y ~ ., data = mydat, mtry = 10, ntree = num.trees, importance = TRUE, importanceSD = TRUE))
  # alternative formulations:
    #(rf1 = randomForest(mydat[,c(2:ncol(mydat))], mydat[,1], mtry = 10, ntree = num.trees, importance = TRUE, importanceSD = TRUE))
    #(rf2 = randomForest(mydat[,c(2:ncol(mydat))], mydat[,1], mtry = 10, ntree = num.trees, localImp = TRUE))
  imp = importance(rf1) #type:	either 1 or 2, specifying the type of importance measure (1 = mean decrease in accuracy, 2 = mean decrease in node impurity).
  rf1$importanceSD
  imp = cbind(imp, "impSD" = rf1$importanceSD)
  imp[order(imp[ , 1], decreasing = TRUE), ]
  
  # plot variable importance using default from the package
  png(paste0(outputDir, "/SA.plots/", m, "_RF_default_plot.png"), units = "in", height = 6, width = 8, res = 150)
  varImpPlot(rf1)
  dev.off()
  rf1$rsq[num.trees] * 100 # variance explained

  # Visualize with RF Explainer package
  
  # min depth of variable in trees
  mindepth_frame = min_depth_distribution(rf1)
  save(mindepth_frame, file = paste0(outputDir, "/SA.RF/", m, ".mindepth.RData"))
  load(paste0(outputDir, "/SA.RF/", m, ".mindepth.RData"))
  
  # plot min depth distribution
  png(paste0(outputDir, "/SA.plots/", m, "_RF_mindepth.png"), units = "in", height = 6, width = 6, res = 150)
  plot_min_depth_distribution(mindepth_frame, mean_sample = "top_trees", min_no_of_trees = 100, mean_scale = TRUE, k = 14) #, main = var)
  dev.off()
  
  # create multi-way importance frame (has multiple evaluation metrics)
  importance_frame = measure_importance(rf1)
  save(importance_frame, file = paste0(outputDir, "/SA.RF/", m,".importance_frame.RData"))
  load(paste0(outputDir, "/SA.RF/", m,".importance_frame.RData"))
  importance_frame
  
  # Figures S2-3 and S2-4: multi-way importance frame
  # panels from Baseline (A.2014) and Predator (P.2014) scenarios were gathered and
  # annotated using GIMP image editing software
  png(paste0(outputDir, "/SA.plots/", m, "_RF_multiway_importance0.png"), units = "in", height = 4.5, width = 5.5, res = 200)
  plot_multi_way_importance(importance_frame, x_measure = "mean_min_depth", y_measure = "mse_increase", size_measure = "p_value", no_of_labels = 5)
  dev.off()
}


# Figure 8, S2-5 and S2-6: summary gridplot
# Figure 8 used "mse.increase" as importance metric
# Figures S2-5 and S2-6 used reamining importance metrics
# panels from Baseline (A.2014) and Predator (P.2014) scenarios were gathered and
# annotated using GIMP image editing software

num.parms = length(parameters)-2

rf.metrics = c("mse_increase", "mean_min_depth", "no_of_nodes", "node_purity_increase", "no_of_trees", "times_a_root", "p_value")
for(importance.metric in rf.metrics){
  m = metric.list[1]
  load(paste0(outputDir, "/SA.RF/", m,".importance_frame.RData"))
  imp.dat = importance_frame[,c("variable", importance.metric)]
  idx = which(levels(imp.dat$variable) == "ration.lo") #remove 'ration.lo' if it exists, because we didn't test it
  if(length(idx) > 0) imp.dat = imp.dat[-idx,]
  for(m in metric.list[2:5]){
    load(paste0(outputDir, "/SA.RF/", m, ".importance_frame.RData"))
    if(length(idx) > 0) importance_frame = importance_frame[-idx,]
    imp.dat = cbind.data.frame(imp.dat, importance_frame[,importance.metric])
  }
  colnames(imp.dat) = c("parameter", metric.list[1:5])
  imp.dat$dateEm = NA
  m = "dateEm"
  load(paste0(outputDir, "/SA.RF/", m, ".importance_frame.RData"))
  if(tag == "A") em.list = c("ATU.crit", "nSpawnDays", "spawn.date.begin", "spawn.date.shape") else em.list = c("ATU.crit", "spawn.date.begin")
  for(v in em.list){
    imp.dat[imp.dat[,"parameter"] == v, 7] = importance_frame[importance_frame[,"variable"] == v, importance.metric]
  }
  
  row.names(imp.dat) = imp.dat[ , 1]
  imp.dat = imp.dat[ , -1]
  imp.mat = t(as.matrix(scale(imp.dat)))
  num.parms = ncol(imp.mat)
  num.mets = nrow(imp.mat)
  
  # create plot
  png(paste0(outputDir,"/SA.plots/gridplot.", importance.metric, ".png"), units = "in", height = 6, width = 6, res = 300)
  par(las = 1, cex = 0.8, oma = c(0, 4, 0, 0))
  if(importance.metric %in% c("mean_min_depth", "p_value")){mycolor = gray.colors(num.parms)} else {mycolor = gray.colors(num.parms)[num.parms:1]}
    image(1:num.mets, 1:num.parms, imp.mat[c(num.mets,1:(num.mets - 1)), num.parms:1], axes= F, col = mycolor, xlab = "", ylab = "")
  axis(1, at = seq(1:num.mets), labels = c("dateEm", "mass.y", "mass.s", "dateOm", "survival", "p.yrlg"))
  axis(2, at = seq(1:num.parms), labels = row.names(imp.dat)[num.parms:1])
  box()
  dev.off()
}

# Evaluate and visualize interactions

# Based on above results, manually choose the top 2 parameters to investigate:
m = "dateEm"
p1 = "ATU.crit"
p2 = "spawn.date.begin"

# Make bi-colored interaction plot (e.g., panels in Figures S2-7 and S2-8, which were gathered and annotated using GIMP image editing software)
png(paste0(outputDir,"/SA.plots/", m, "_RF_interactions_top2_grid.png"), units = "in", height = 4.5, width = 5.5, res = 200)
plot_predict_interaction(forest = rf1, data = mydat, variable1 = p1, variable2 = p2) #have to use formula version in RF call, above
dev.off()

# Optional, deeper analysis - takes a long time to run
(vars <- important_variables(importance_frame, k = 5, measures = c("mean_min_depth", "no_of_trees")))
interactions_frame = min_depth_interactions(rf1, vars)
save(interactions_frame, file = paste0(outputDir, "/SA.RF/", m, ".interactions_frame.RData"))
load(paste0(outputDir, "/SA.RF/", m, ".interactions_frame.RData"))
head(interactions_frame[order(interactions_frame$occurrences, decreasing = TRUE), ])
png(paste0(outputDir,"/SA.plots/", m, "_RF_interactions.png"), units = "in", height = 6, width = 10, res = 150)
plot_min_depth_interactions(interactions_frame)
dev.off()




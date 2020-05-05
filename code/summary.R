
# Create quick summary statistics file for each run.

# Set up. ----------------------------------------------------------------------

# load package
library(dplyr)

# load number of fish from model input
nFish = parameters["salmon","nFish"]

# create fish output data frame from model output
fo = as.data.frame(salmon.finalstep)

# convert dates
fo$dateSp = as.POSIXct(fo$dateSp, origin="1970-01-01")
fo$dateEm = as.POSIXct(fo$dateEm, origin="1970-01-01")
fo$dateOm = as.POSIXct(fo$dateOm, origin="1970-01-01")
fo$dateDi = as.POSIXct(fo$dateDi, origin="1970-01-01")
fo$datePr = as.POSIXct(fo$datePr, origin="1970-01-01")
fo$dateMo = fo$dateDi; fo$dateMo[fo$survive != 0] = NA

# Calculate and store summary statistics. --------------------------------------

spnd = fo %>% group_by(as.Date(dateSp)) %>% summarise(Count = length(pid)); spnd = spnd[-nrow(spnd),]; colnames(spnd)[1] = "Date"
emgd = fo %>% group_by(as.Date(dateEm)) %>% summarise(Count = length(pid)); emgd = emgd[-nrow(emgd),]; colnames(emgd)[1] = "Date"
smolts = fo %>% group_by(as.Date(dateOm)) %>% summarise(Count = length(pid)); smolts = smolts[-nrow(smolts),]; colnames(smolts)[1] = "Date"
preds = fo %>% group_by(as.Date(datePr)) %>% summarise(Count = length(pid)); preds = preds[-nrow(preds),]; colnames(preds)[1] = "Date"
morts = fo %>% group_by(as.Date(dateDi)) %>% summarise(Count = length(pid)); morts = morts[-nrow(morts),]; colnames(morts)[1] = "Date"
ssmort = fo %>% group_by(as.Date(dateMo)) %>% summarise(Count = length(pid)); ssmort = ssmort[-nrow(ssmort),]; colnames(ssmort)[1] = "Date"

# spawning dates
spawnDates = substr(summary(fo$dateSp[fo$dateSp > as.POSIXct(paste0(climate.year-1,"-09-01"))]), 1, 10)
# emergence dates
emrgDates = substr(summary(fo$dateEm[fo$dateEm > as.POSIXct(paste0(climate.year-1,"-09-01")) & fo$emrg ==1]), 1, 10)
# outmigration dates
smoltDates = substr(summary(fo$dateOm[fo$dateOm > as.POSIXct(paste0(climate.year-1,"-09-01")) & fo$survive ==2]), 1, 10)
# predation dates
predDates = substr(summary(fo$datePr[fo$datePr > as.POSIXct(paste0(climate.year-1,"-09-01"))]), 1, 10)
# all mortality dates
mortDates = substr(summary(fo$dateDi[fo$dateDi > as.POSIXct(paste0(climate.year-1,"-09-01"))]), 1, 10)
# stochastic mortality dates
ssMortDates = substr(summary(fo$dateDi[fo$dateDi > as.POSIXct(paste0(climate.year-1,"-09-01")) & fo$survive == 0]), 1, 10)

# subyearling mass
subYMass = round(quantile(fo$weight[fo$survive ==2]), 2)
# yearling mass
ylgMass = round(quantile(fo$weight[fo$survive ==1]), 2)

# number of surviving subyearlings
subyearlings = length(fo$weight[fo$survive ==2])
# number of surviving yearlings
yearlings = length(fo$weight[fo$survive ==1])

# proportion of yearlings of total survivors
prop.yearlings = round(yearlings/(subyearlings + yearlings), 3)
# proportion of survivors over total simulated fish
survival = round((subyearlings + yearlings) / nFish, 3)


# Write summary statistics to text file. ---------------------------------------

summaryDir = paste0(outputDir, "/quick.summary.", climate.year, ".", iter, ".txt")
file.create(summaryDir)
summary.info<- c(paste0("survival: ", survival),
                 paste0("pYrlngs: ", prop.yearlings), 
                 paste0("nSubyearlings: ", subyearlings),
                 paste0("nYearlings: ", yearlings),
                 paste0("subYMass: ", subYMass),
                 paste0("ylgMass: ", ylgMass),
                 paste0("spawnDates: ", spawnDates),
                 paste0("emrgDates: ", emrgDates),
                 paste0("ssMortDates: ", ssMortDates),
                 paste0("smoltDates:", smoltDates))
write(c("Basic Summary:", summary.info), file = summaryDir)

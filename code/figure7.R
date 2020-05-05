
#----- PREPARE DATA ------------------------------------------------------------

# load packages
library(SSN)
library(RColorBrewer)
library(raster)
library(rgdal)

# save colors
light<- rgb(252,217,156,150,NULL,255) #hex #fcd99c96 orange
medium<- rgb(234,173,68,200,NULL,255) #hex #eaad44c8 orange
dark<- rgb(244,155,2,255,NULL,255) #hex #f49b02ff orange

# choose scenario of data
spp = "salmon"
climate.year = 2014 # 2014 for Cool-Predator, 2015 for Warm-Predator
tag = "P" # "P" for predator
nFish = 2241
iter = 1

# set up directories
mydir = getwd()
outputDir = dataDir = paste0("data.out/", tag, ".", climate.year)
plotDir= "data.out/Figures"
if (! dir.exists(file.path(paste0(mydir, "/", plotDir)))) {
  dir.create(file.path(paste0(mydir, "/", plotDir)))
}

# load salmon array
load(paste0(outputDir, "/",spp,".array.",climate.year,".",iter,".RData")) 

# load snoqualmie ssn object
sno.ssn = importSSN("data.in/sno.ssn", predpts='preds')
obs.df = getSSNdata.frame(sno.ssn,"preds")
if(is.factor(obs.df$rid)) obs.df$rid = as.numeric(levels(obs.df$rid)[as.numeric(obs.df$rid)])

# filter predation for final time step
# save survive column where -1 indicates death by predator
predation.final<- salmon.array[,"survive",730]
# if salmon not killed by predator, change to NA
predation.final[predation.final != -1]<- NA
# if salmon killed by predator, change to 1
predation.final[predation.final == -1]<- 1

# filter predation for all time steps
# save survive column where -1 indicates death by predator for each fish (rows) over time (columns)
daily.predation<- salmon.array[,"survive",]
# if salmon not killed by predator, change to NA
daily.predation[daily.predation != -1]<- NA
# if salmon killed by predator, change to 1
daily.predation[daily.predation == -1]<- 1

# find when & where predation occurred
segs = salmon.array[,"seg",]
td = daily.predation * segs
foo = td; foo[!is.na(foo)] = 1; foo[is.na(foo)] = 0
when.eaten = apply(foo, 1, which.max); when.eaten[when.eaten == 1] <- NA
seg.eaten = NULL
for(f in 1:nFish){
  tmp = td[f, when.eaten]
  seg.eaten[f] = tmp[!is.na(tmp)][1]
}
when.where = cbind(when.eaten, seg.eaten); colnames(when.where)=c("when.eaten","rid")
when.summary = when.where[,"when.eaten"]/2 - 122 # minus 122 because the simulation started on 1 Sep, not 1 Jan.


#---- TEMPORAL: WHEN SALMON EATEN BY PREATOR -----------------------------------

quantile(when.summary, na.rm = T)
hist(when.summary, xlab="Julian Date", main="")
  
png(paste0(getwd(),"/",plotDir,"/Figure7Temporal.png"), width = 6, height = 6, units = "in", res = 150)
  par(las = 1, cex = 1.5)
  mylabels = c("1-Jan", "19-Feb", "10-Apr", "30-May", "19-Jul")
  ats = c(0, 50, 100, 150, 200)
  y.data = when.summary[!is.na(when.summary)]
  h = hist(y.data, xlab = "Date", main = "", ylab = "Density", col = light, axes = F, freq = F, ylim = c(0,0.03))
  x.data = seq(min(y.data, na.rm = T), max(y.data, na.rm = T), length.out = 100)
  y.data2 = dnorm(x.data,mean(y.data), sd(y.data))
  lines(x.data, y.data2, col = dark, lwd = 4)
  axis(1, at = ats, labels = mylabels, cex.axis = 0.8)
  ats2 = c(0,0.01,0.02,0.03); mylabels2 = ats2
  axis(2, at = ats2, labels = mylabels2, cex.axis = 0.8)
  box()
dev.off()


#---- SPATIAL: WHERE SALMON EATEN BY PREATOR -----------------------------------

# count of how many times a reach had a predation event
where.summary = when.where[,"rid"]
where.summary = where.summary[!is.na(where.summary)]
where.summary = as.data.frame(table(where.summary)); colnames(where.summary) = c("rid","Pred.Freq")
if(is.factor(where.summary$rid)) where.summary$rid= as.numeric(levels(where.summary$rid)[as.numeric(where.summary$rid)])

# map
sno.ssn@data = cbind.data.frame(sno.ssn@data,"sort.order" = as.integer(rownames(sno.ssn@data)),"col.class" = 1, "Predation" = 0)
reaches = sno.ssn@data
if(is.factor(reaches$rid)) reaches$rid= as.numeric(levels(reaches$rid)[as.numeric(reaches$rid)])
row.names(reaches) = NULL
reaches$sort.order = as.integer(rownames(reaches))
reaches = reaches[order(reaches$rid),]
where.summary = where.summary[order(where.summary$rid),]
rid = seq(0:887)-1
where.summary2 = merge(as.data.frame(rid), where.summary, by = "rid", all.x = T)
reaches$Predation = where.summary2$Pred.Freq
reaches = reaches[order(reaches$sort.order),]
sno.ssn@data = reaches
basin <- readOGR(paste0(getwd(),"/data.in/shapefiles"),"Basin_snq")

# set plotting extent so that data overlay the background layer
ex <- extent(basin)
ex@xmin <- sno.ssn@bbox[1]
ex@xmax <- sno.ssn@bbox[3]
ex@ymin <- sno.ssn@bbox[2]
ex@ymax <- sno.ssn@bbox[4]

png(paste0(getwd(),"/",plotDir,"/Figure7Spatial.png"), width = 8, height = 6, units = "in", res = 300)
par(mar = c(1, 0, 4, 0))

# plot background
plot(basin, col="gray60", border=NA)

# assign colors based on predation data
mytitle <- "Predation events"
color.length = 7
cb <- c("gray40","#ffc04d","#ffb733","#ffae1a","#ffa500","#e69500","#cc8400")
colseq <- c(0, round(seq(min(sno.ssn@data[,"Predation"], na.rm = T), max(sno.ssn@data[,"Predation"], na.rm = T), length.out = color.length), 0))
left <- colseq[1:color.length]; rght<- colseq[2:(color.length+1)]
for(n in 1:length(cb)) {
  sno.ssn@data$col.class[sno.ssn@data[,"Predation"] >= left[n] & sno.ssn@data[,"Predation"] <= rght[n]] = n
}

# plot gray stream lines for areas without predation
for (i in 1:length(sno.ssn@lines)) for (j in 1:length(sno.ssn@lines[[i]])) lines(sno.ssn@lines[[i]]@Lines[[j]]@coords,col="gray40",lwd = 5*sno.ssn@data[i,"afvFlow"]+0.4)

# plot predation data as colored stream lines with line thickness based on cumulative drainage area
for (i in 1:length(sno.ssn@lines)) {
  for (j in 1:length(sno.ssn@lines[[i]])) {
    lines(sno.ssn@lines[[i]]@Lines[[j]]@coords,col=cb[sno.ssn@data[i,"col.class"]],lwd = 5*(sno.ssn@data[i,"afvFlow"]+0.4))
  }
}

# add legend
leglabs <- paste(round(left,1), "to", round(rght,1))
legend("right", legend = leglabs, title=mytitle,bty = "n", pch = 19, col = cb, cex = 0.8)

# add scale bar
xshift=yshift=0
rect(xleft=ex[1]+5000-xshift,ybottom=ex[3]+3000-yshift,xright=ex[1]+7500-xshift,ytop=ex[3]+3500-yshift)
rect(xleft=ex[1]+7500-xshift,ybottom=ex[3]+3000-yshift,xright=ex[1]+10000-xshift,ytop=ex[3]+3500-yshift,col=1)
rect(xleft=ex[1]+10000-xshift,ybottom=ex[3]+3000-yshift,xright=ex[1]+12500-xshift,ytop=ex[3]+3500-yshift)
rect(xleft=ex[1]+12500-xshift,ybottom=ex[3]+3000-yshift,xright=ex[1]+15000-xshift,ytop=ex[3]+3500-yshift,col=1)
segments(x0=ex[1]+5000-xshift,y0=ex[3]+3000-yshift,x1=ex[1]+5000-xshift,y1=ex[3]+2500-yshift)
segments(x0=ex[1]+10000-xshift,y0=ex[3]+3000-yshift,x1=ex[1]+10000-xshift,y1=ex[3]+2500-yshift)
segments(x0=ex[1]+15000-xshift,y0=ex[3]+3000-yshift,x1=ex[1]+15000-xshift,y1=ex[3]+2500-yshift)
text(x=ex[1]+5000-xshift,y=ex[3]+1500-yshift,"0",cex=0.8); text(x=ex[1]+10000-xshift,y=ex[3]+1500-yshift,"5",cex=0.8)
text(x=ex[1]+15000-xshift,y=ex[3]+1500-yshift,"10",cex=0.8); text(x=ex[1]+18000-xshift,y=ex[3]+1500-yshift,"km",cex=0.8) 

# add north arrow
arrows(ex[1]+2000-xshift,ex[3]+2300-yshift,ex[1]+2000-xshift,ex[3]+4000-yshift,length=0.1,lwd=5)
text(ex[1]+2000-xshift,ex[3]+1000-yshift,"N")

dev.off()

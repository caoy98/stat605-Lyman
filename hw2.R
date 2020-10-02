#Name: Yuan Cao, Email: cao234@wisc.edu

#remove all values and load package
rm(list=ls())
require("FITSio")

#exponential smoothing to remove noise
expsmooth = function(y, w){
  len=length(y)
  p=rep(0,len)
  p[1]=y[1]
  for (i in 2:len) {
    p[i]=w*y[i]+(1-w)*p[i-1]
  }
  return(p)
}
#function to find possible trough point
localmin = function(y){
  q = quantile(y, probs = 0.25)
  len = length(y)
  diff.y = diff(y)
  min_index = rep(0,len)
  for (i in 2:(len-2)) {
    if (diff.y[i]<0 & diff.y[i-1]<0 & diff.y[i+1]>0 & y[i+1]<q) {
      min_index[i+1]=i+1
    }
  }
  min_index=min_index[which(min_index!=0)]
  return(min_index)
}

#read CB58 data, standardize and smooth it
cb58=readFrameFromFITS("cB58_Lyman_break.fit")
cb58_standard = scale(cb58$FLUX)
cb58length = length(cb58$FLUX)
expcb58 = expsmooth(cb58_standard, 0.05)
cb58min_index = which.min(expcb58)

#read 100 files and do computation
filename = dir("data/")
distance = c(); spectrumID = c(); index = c()
for (name in filename) {
  spec = readFrameFromFITS(paste("data", name, sep = "/"))
  
  #exclude bad observations
  goodobs = length(which(spec$and_mask == 0))
  if (goodobs > cb58length) {
    spec$flux[which(spec$and_mask != 0)] = NA
    flux = na.omit(spec$flux)
  }
  else {
    flux = spec$flux
  }
  #standardize data and smooth it
  spec_standard = scale(flux)
  expspec=expsmooth(spec_standard, 0.05)
  
  #find possible trough point
  m = localmin(expspec)
  speclength = length(expspec)
  m = m[which(m>cb58min_index)]
  m = m[which(m+cb58length<speclength)]
  
  #for each possible trough, compute correlation and find largest one
  bestcor = -1
  for (tp in m) {
    start_index = tp-cb58min_index
    y = expspec[start_index:(start_index+cb58length-1)]
    corr = cor(expcb58, y, method = "spearman")
    if (corr > bestcor) {
      bestcor = corr
      begin_index = start_index
    }
  }
  distance = c(distance, bestcor)
  spectrumID = c(spectrumID, name)
  index = c(index, begin_index)
}

#output result to csv file
result = data.frame(distance = 1-distance, spectrumID = spectrumID, i = index)
result = result[order(result[,1],decreasing=F),]
write.csv(result, file = "hw2.csv", row.names = F)

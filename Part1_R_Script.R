# Loading XCMS package functions and also ggplot2 and gridExtra packages for plotting
library(xcms)
library(ggplot2)
library(gridExtra)

PeakPicking.function <- function(data.path, snthresh_val1, snthresh_val2, noise_val1, noise_val2, ppm_val1, ppm_val2, 
                          peakwidth_val1_low, peakwidth_val1_top, peakwidth_val2_low, peakwidth_val2_top,
                          prefilter_val1_low, prefilter_val1_top, prefilter_val2_low, prefilter_val2_top,
                          feature_mz, feature_ppm)
{
  # loading the data by separating the no-collision (msLevel. = 1) and high-collision energy (msLevel. = 2) scans 
  xcmsF1 <- readMSData(data.path, msLevel. = 1, mode="onDisk")
  xcmsF2 <- readMSData(data.path, msLevel. = 2, mode="onDisk")
  # changing the label to ensure XCMS functions work as intended
  xcmsF2@featureData@data$msLevel <- 1
  # defining the parameters using CentWave, the peak-picking algorithm
  cwp1 <- CentWaveParam(snthresh = snthresh_val1, noise = noise_val1, ppm = ppm_val1, 
                        peakwidth = c(peakwidth_val1_low, peakwidth_val1_top), 
                        prefilter = c(prefilter_val1_low, prefilter_val1_top))
  cwp2 <- CentWaveParam(snthresh = snthresh_val2, noise = noise_val2, ppm = ppm_val2, 
                        peakwidth = c(peakwidth_val2_low, peakwidth_val2_top), 
                        prefilter = c(prefilter_val2_low, prefilter_val2_top))
  # peak-picking with the parameters from CentWave
  peaksF1 <- findChromPeaks(xcmsF1, param = cwp1)
  peaksF2 <- findChromPeaks(xcmsF2, param = cwp2)
  peaksMSLevel1 <- peaksF1@msFeatureData$chromPeaks
  peaksMSLevel2 <- peaksF2@msFeatureData$chromPeaks
  # saving the filtered peaks as .CSV
  #write.csv(peaksMSLevel1,'peaksMSLevel1.csv')
  #write.csv(peaksMSLevel2,'peaksMSLevel2.csv')
  # finding the peaks of the feature of interest in the peak-picked data
  features <- chromPeaks(peaksF1, mz = feature_mz, ppm = feature_ppm) 
  
  if (length(features) == 0){
    stop("Nothing has been found. Change the peak-picking parameters!")
  }
  
  # saving the peaks of the feature of interest in the peak-picked data as .CSV
  write.csv(features,'features.csv')
  # evaluating all peaks that have the same m/z within ppm tolerance
  for(i in 1:nrow(features))
  {
    print(i)
    fmz <- features[i,"mz"]
    frt <- features[i,"rt"]
    frtmin <- features[i,"rtmin"]
    frtmax <- features[i,"rtmax"]
    eic <- chromatogram(xcmsF1, mz = c(fmz-0.01, fmz+0.01), rt = c(frtmin-1, frtmax+1))
    hc <- chromPeaks(peaksF2)
    sel <- hc[which(hc[,"rtmin"] < frt & hc[,"rtmax"] > frt),]
    eic2 <- lapply(1:nrow(sel),function(x) chromatogram(xcmsF2, mz = c(sel[x,"mz"]-0.01,sel[x,"mz"]+0.01), 
                                                        rt = c(frtmin-1, frtmax+1)))
    # selecting the fragments with the same "elution" profile as the feature of interest
    c <- lapply(1:nrow(sel),function(x) cor(intensity(eic[1,1]), intensity(eic2[[x]][1,1])[-1],
                                            use = "pairwise.complete.obs"))
    # selecting the fragments that have a correlation coefficient to the feature of interest > 0.7.
    aif <- sel[which(c > 0.7),]
    
    if (length(aif) > 11) {
      aif <- aif[which(aif[,"mz"] <= fmz),]
    } else next
    
    pspec <- as.data.frame(aif[,c("mz", "into")])
    # saving the fragmentation spectrum (spectra) of the feature of interest as .CSV
    fname <- paste('features', i, '.csv', sep='')
    write.csv(pspec, fname)
    print(paste("results writen to", fname))
  }
}

"Calling the function and supplying: the data path to the file containing the data; signal/noise ratio (snthresh) for msLevel. = 1;
signal/noise ratio (snthresh) for msLevel. = 2; noise for msLevel. = 1; noise for msLevel. = 2; ppm for msLevel. = 1; 
ppm for msLevel. = 2; peakwidth lower boundary for msLevel. = 1; peakwidth top boundary for msLevel. = 1;
peakwidth lower boundary for msLevel. = 2; peakwidth top boundary for msLevel. = 2; prefilter lower boundary for msLevel. = 1;
prefilter top boundary for msLevel. = 1; prefilter lower boundary for msLevel. = 2; prefilter top boundary for msLevel. = 2.
the m/z of the feature of interest; tolerance ppm"


PeakPicking.function("MESA_COMBI_BIO_P2_LP_iQC150.mzML", 5, 5, 100, 10, 20, 20, 2, 20, 2, 20, 3, 100, 3, 5, 286.1438, 1)


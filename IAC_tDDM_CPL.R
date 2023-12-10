pacman::p_load(
    tidyverse,
    ggplot2,
    DEoptim,
    Rcpp,
    plyr,
    parallel,
    pbmcapply
)

# https://stackoverflow.com/questions/47044068/get-the-path-of-current-script
# get path of source file
getCurrentFileLocation <-  function()
{
    this_file <- commandArgs() %>% 
        tibble::enframe(name = NULL) %>%
        tidyr::separate(col=value, into=c("key", "value"), sep="=", fill='right') %>%
        dplyr::filter(key == "--file") %>%
        dplyr::pull(value)
    if (length(this_file)==0)
    {
        this_file <- rstudioapi::getSourceEditorContext()$path
    }
    return(dirname(this_file))
}
# sets path on source file location
script_path <- getCurrentFileLocation()
setwd(script_path)

## Fit the tDDM in parallel for multiple subjects


RcppParallel::setThreadOptions(numThreads = 1) #this is critical for running on Mac OS for some reason.

### --- how to run the tDDM:
# load the c++ file containing the functions to simulate the time-varying DDM
sourceCpp("tSSM_Rcpp.cpp")

# read in behavioral data from HEALTH cue only
# this csv contains health and taste difference computed as healthy option - unhealthy options
# when both options were rated equally they were removed
# healthy option is the one with highest health rating
raw <- read_csv("filtered_data_for_DDM_with_health_as_reference.csv") %>% 
    dplyr::rename(
        td = dTaste,
        hd = dHealth,
        rt = RT
    ) %>% 
    mutate(
        label = paste(ID, context, sep = "_"),
        subID = as.numeric(as.factor(label)),
        rt = rt / 1000 # use second instead of ms
    )

# data pre-proc

# step 1: 10 equally spaced bins for health and taste ratings
# step 2: z-score those bins
# step 3: signed difference from z-scores (healthy - unhealthy)
dataBeh <-
    raw %>% 
    mutate(
        bin10_healthIzq = as.numeric(cut(healthIzq, breaks = 10, labels = 1:10)) %>% scale(),
        bin10_healthDer = as.numeric(cut(healthDer, breaks = 10, labels = 1:10)) %>% scale(),
        bin10_tasteIzq  = as.numeric(cut(tasteIzq, breaks = 10, labels = 1:10)) %>% scale(),
        bin10_tasteDer  = as.numeric(cut(tasteDer, breaks = 10, labels = 1:10)) %>% scale(),
        bin_hd = if_else(bin10_healthIzq > bin10_healthDer, bin10_healthIzq - bin10_healthDer, bin10_healthDer - bin10_healthIzq),
        bin_td = if_else(bin10_tasteIzq > bin10_tasteDer, bin10_tasteIzq - bin10_tasteDer, bin10_tasteDer - bin10_tasteIzq),
        unidentifiable = if_else(healthIzq == healthDer | tasteIzq == tasteDer, 1, 0)
    )

dataBeh$abstd = abs(dataBeh$td)
dataBeh$abshd = abs(dataBeh$hd)
dataBeh$logRT = log(dataBeh$rt)

# assign negative RTs to unhealthy option
dataBeh <- dataBeh %>% 
    mutate(
        choseL = if_else(sub_healthy_choice == "healthy", 1, 0),
        # if unhealthy is choosen then rt is negative
        RTddm = if_else(choseL == 0, rt * -1, rt)
    )
ntrials = length(dataBeh$rt)

# bin the prob. density space of RTs:
# based on xpos, we can see the size of the time bins
# if you want to adjust the bins, change the sequence length (-5,5) or 
# length.out, which tells you how many divisions are put into the sequence

# just to make sure, plot the RT distributions (with positive and negative sign)
hist(dataBeh$RTddm)  # ~-4 to 4 with significant number of observations
range(dataBeh$RTddm) # -4.5 to 4.6

xpos = seq(-5,5,length.out=1024) # this is then OK
dt = xpos[2] - xpos[1]
dataBeh$RTddm_pos = 0

# gets the probability within the pdf
for (i in 1:ntrials) {
	dataBeh$RTddm_pos[i] = which.min(abs(xpos - dataBeh$RTddm[i]))
}


## define fitting functions ---- 

ll_ddm2 <- function(x, dataBeh2, vd, hd) {
  
  d_v = x[1] # weighting for attribute 1 (here: taste)
  d_h = x[2] # weighting for attribute 2 (here: health)
  thres = x[3] # threshold
  nDT = x[4] # non-decision time
  tHin = x[5] # time Health in (relative start time of health vs. taste)
  bias = x[6] # starting point bias
  
  probs = NULL
  for (i in 1:length(vd)) {
    
    rts = ddm2_parallel(d_v,d_h,thres,nDT,tHin,bias,vd[i],hd[i],3000) # 3000 is the number of simulated RTs
    rts = rts[rts!=0]
    xdens = density(rts, from=-5, to=5, n=1024, bw=0.11) #matches the prob. space of RTs from above
    idx = which(dataBeh2$bin_td==vd[i] & dataBeh2$bin_hd==hd[i])
    # this is a vector of probabilities for all trial with same taste and health difference
    probs = c(probs, dt*xdens$y[dataBeh2$RTddm_pos[idx]])
  }
  
  probs[probs==0] = 1e-100
  return (-sum(log(probs)))
  
}

# fit individually for each subject

fitSub <- function(s, dataBeh) {
  fits=matrix(0,1,9)
  idx = which(dataBeh$subID==s)
  dataBeh2 = dataBeh[idx,]
  label <- unique(dataBeh2$label)
  cat(NULL,file=paste0(label, ".csv"))
  
  #data1 = ddply(dataBeh2, .(td, hd), summarize, acc= mean(choseL))
  # td and hd were changes to 10 equally spaced bins
  # this should reduce computation time
  data1 = ddply(dataBeh2, .(bin_td, bin_hd), summarize, acc= mean(choseL))
  
  #vd = data1$td # value difference for attribute 1 (here: taste)
  #hd = data1$hd # value difference for attribute 2 (here: health)
  
  # binned differences for taste and health
  vd = data1$bin_td
  hd = data1$bin_hd
  
  # boundaries for parameters
  # test increasing boundaries 
  # lower <- c(-2,-2,0.6,0.01,-1,-1)
  # upper <- c(2,2,3,1,1,1)
  
  lower <- c(-4,  -4, 0.6, 0.01, -4, -1)
  upper <- c( 4,   4, 3,   1,     4,  1)
  
  # force the optimizer for a given number of iterations
  # get loglikelihood of every iteration per subject
  fit_s = DEoptim(ll_ddm2, lower, upper, DEoptim.control(itermax = 150), dataBeh2=dataBeh2, vd=vd, hd=hd)
  fits[1,1:6] = fit_s$optim$bestmem #fitted parameters
  fits[1,7] = fit_s$optim$bestval #LL
  fits[1,8] = 2*fit_s$optim$bestval + length(lower)*log(length(vd)) #BIC
  fits[1,9] = 2*fit_s$optim$bestval + 2*length(lower) # AIC
  # write to csv subject_context
  f = as.data.frame(matrix(unlist(fits),ncol=9,byrow=TRUE))
  names(f)<-c("d_t", "d_h", "thres", "nDT", "timeHin", "bias", "LL", "BIC", "AIC")
  write_csv(f, paste0(label, ".csv"))
  saveRDS(fits_s, paste0(label, "_diagnostics", ".rds"))
  return(fits)
}

inputs = 1:length(unique(dataBeh$subID))
numCores <- detectCores()
fits = pbmcapply::pbmclapply(inputs, fitSub, mc.cores = numCores, dataBeh=dataBeh)
fitSub(1, dataBeh)

#save now in case unlist fails
fitsF=fits
write.csv(fitsF, file = "fits_IAC_example.csv")

fits = as.data.frame(matrix(unlist(fits),ncol=9,byrow=TRUE)) %>% 
    bind_cols(or_id)
names(fits)<-c("d_t", "d_h", "thres", "nDT", "timeHin", "bias", "LL", "BIC", "AIC", "ID", "fit_ID")
#will overwrite previous save, but that is intended
fitsF=fits
write.csv(fitsF, file = "fits_ddm_nbo.csv")

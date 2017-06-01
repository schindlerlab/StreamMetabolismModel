## KathiJo Jankowski
###### December 2016
## This code is written to run BaMM via a Mac computer that has ADMB installed. 
## Here: 'system()' calls the command line prompt through R and executes the commands
## To run command line via R on a PC, use the "shell" command - e.g., shell("BaMM")

##### Notes on results: 
## To look at maximum likelihood estimates for model run, see "BaMM.rep" file created in current directory
## R2admb package has functions to input results back into R for analysis, review, but the .rep file is 
## not currently formatted correctly to use that function

##########################################################################
#####  Compile .tpl file into executable file: 
# First, check that ADMB is properly installed and responding:
system("admb")

# Compile model tpl file into executable file "BaMM.exe"
# Important: Mac compiled file only runs on Mac, PC compiled only runs on PC
system('admb BaMM.tpl')

###########################################################################
###### Set working directory to where the BaMM.exe, BaMM.cfg, BaMM.pin, and BaMM.dat files are located
# must all be in the same directory
setwd("/Users/kj/Documents/MyDocuments/GradSchoolUW/Dissertation/BaMMUpdates/2010Streams/7th2source")
setwd("/Users/kj/Documents/MyDocuments/GradSchoolUW/Dissertation/BaMMUpdates/TwoSourceR/EmilysStreams/Smith")
setwd("/Users/kj/Documents/MyDocuments/MBLPostdoc/Data/Streams/Metabolism/Metabolism Oct 2016/Cascavel")

# Check which files are in that directory
system('ls')

# Run BaMM without mcmc
system('./BaMM')

# Run BaMM with 200,000 total iterations (draws from the posterior), save every 100th draw
system('./BaMM -mcmc 200000 -mcsave 100')

# Puts results (parameter value, likelihood value) from each of the saved draws into "BaMM_mcmc.dat" file 
system('./BaMM -mceval')

###########################################################################
### Code for trying to use R2admb to run BaMM directly
### Ended up running ADMB through "system" call in R, not using R2admb
## Setup ADMB environment (not totally sure I have to do this. on lab computer?  everytime?)
setup_admb(ADMB_HOME)

## Compile model .tpl file into executable file (.exe)
compile_admb("BaMM", verbose=TRUE, admb_errors="warn")

## Run model 
run_admb("BaMM", verbose=TRUE, admb_errors="warn")

# run with mcmc 
setwd("C:/Jankowski_BaMM/BaMMTests/BaMM.A3")
modMC1 <- run_admb("BaMM", verbose=TRUE, admb_errors="warn", mcmc=TRUE,mcmc.opts=mcmc.control(mcmcpars=c("k20","Eb", "Rref"), mcmc=100000, mcsave=1000))
# to set mcmc options use mcmc.opts=mcmc.control()
# need to specify at least one parameter in mcmcpars

results1<-read_admb("BaMM", mcmc=TRUE, mcmc.opts=mcmc.control(mcmcpars=c("k20","Eb", "Rref"), mcmc=100000, mcsave=1000)) 
"Error in scan(file, what, nmax, sep, dec, quote, skip, nlines, na.strings,  : 
scan() expected 'a real', got 'objective'
In addition: Warning message:
In read_rep(fn) : report file in non-standard format"

"Error in scan(file, what, nmax, sep, dec, quote, skip, nlines, na.strings,  : 
scan() expected 'a real', got 'objective'
In addition: Warning messages:
1: closing unused connection 4 (BaMM.mcinfo) 
2: closing unused connection 3 (BaMM.mcinfo) 
3: In read_rep(fn) : report file in non-standard format"

## Trying different datasets with no mcmc
# Dataset D: Rref =240, Eb=0.65
setwd("C:/Jankowski_BaMM/BaMMTests/BaMM.D")

modD <- run_admb("BaMM", verbose=TRUE, admb_errors="warn")

# Dataset E: Rref = 120, Eb = 1.0
setwd("C:/Jankowski_BaMM/BaMMTests/BaMM.E")
modE <- run_admb("BaMM", verbose=TRUE, admb_errors="warn")

# Dataset F: Rref=240, Eb=1.0
setwd("C:/Jankowski_BaMM/BaMMTests/BaMM.F")
modF <- run_admb("BaMM", verbose=TRUE, admb_errors="warn")



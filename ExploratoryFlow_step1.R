#Summary: This script analyze flow cytometry data. The problem was the batch effect. In the analysis I merged two groups of files (in the script the directories are called fcs_1 and fcs_2)
# I created two different directories because in fcs_1 there are 2 additional columns in the matrix (additional flow markers) that were not used in the analysis. We deleted those columns
# and then merget the two flowset in a single one called fs. Then the problem was the batch effect and I used warpSet from fdaNorm to try to get the rid of the technical error. This was my
# very first try of doing an analysis in R 2 years and a half ago so it may need to be reviewed a bit. 

#clean workspace
rm(list = ls())

#required packages
library(rstudioapi)
library(flowCore)
library(FlowSOM)
library(ggplot2)
library(dplyr)
library(stringr)
library(flowDensity) 
library(flowStats) 
library(flowVS)

# Set PrimaryDirectory where this script is located
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
PrimaryDirectory <- getwd()
PrimaryDirectory

# Define fcs_1 directory (Hanna's file)
fcs_1 <- "fcs_1"
FCS1Directory <- paste(PrimaryDirectory, fcs_1, sep = "/")
dir.create(FCS1Directory)

# Define fcs_2 directory (Rupa's files)
fcs_2 <- "fcs_2"
FCS2Directory <- paste(PrimaryDirectory, fcs_2, sep = "/")
dir.create(FCS2Directory)

# Define workingDirectory and create the directory
wdName <- "Working_Directory"
workingDirectory <- paste(PrimaryDirectory, wdName, sep = "/")
dir.create(workingDirectory)

#List FCSfiles_1
FCSfiles_1 <- list.files(FCS1Directory, pattern = ".fcs$", full = FALSE)

## flowSet_1 
fs_1 <- read.flowSet(files = FCSfiles_1, path = FCS1Directory, truncate_max_range = FALSE)
fs_1 <- fs_1[, -c(13,14)]

#List FCSfiles_2
FCSfiles_2 <- list.files(FCS2Directory, pattern = ".fcs$", full = FALSE)

## flowSet_2 
fs_2 <- read.flowSet(files = FCSfiles_2, path = FCS2Directory, truncate_max_range = FALSE)

## flowSet merged 
fs <- rbind2(fs_1, fs_2)

#Arcsinh transformation (this is necessary when you have flow citometry data in fcs format instead of csv)
# It is well explained here https://wiki.centenary.org.au/display/SPECTRE/Data+transformation. When I started I did not know these data could be 
#exported as csv

#choose cofactors for each channel
cfs <- c(750, 2200, 2300, 800, 1000, 3200, 1500, 700, 1100, 1100, 600)

#Rename the colnames with the markers corresponding to the fluorochrome (some of this channels, like FSC, SSC,.. are not of interest)
#they do not need to be renamed
colnames(fs)[c(7:9, 11:18)] <- c("KLRG1", "CD45RA","CD27","CD160","PD1","CD28","CD56",
                                 "CD57","CD8","CD3","CCR7")

channels <- colnames(fs)[c(7:9, 11:18)]

#Transform data
fs_t <- transFlowVS(fs, channels = channels, cofactor = cfs)
flowViz.par.set(theme =  trellis.par.get(), reset = TRUE)

#Look at the density plots for each channel after transformation
densityplot(~KLRG1, fs_t)
densityplot(~CD45RA, fs_t)
densityplot(~CD27, fs_t)
densityplot(~CD160, fs_t)
densityplot(~PD1, fs_t)
densityplot(~CD28, fs_t)
densityplot(~CD56, fs_t)
densityplot(~CD57, fs_t)
densityplot(~CD8, fs_t)
densityplot(~CD3, fs_t)
densityplot(~CCR7, fs_t)

## Try to correct the signal in the channels with the biggest technical issues
#Warpset from fdaNorm
fs_fda <- warpSet(fs_t, stains = channels[-c(3,4)])
densityplot(~KLRG1, fs_fda)
densityplot(~CD45RA, fs_fda)
densityplot(~CD27, fs_fda)
densityplot(~CD160, fs_fda)
densityplot(~PD1, fs_fda)
densityplot(~CD28, fs_fda)
densityplot(~CD56, fs_fda)
densityplot(~CD57, fs_fda)
densityplot(~CD8, fs_fda)
densityplot(~CD3, fs_fda)
densityplot(~CCR7, fs_fda)

#Create a directory with all the files transformed and normalized
if(!dir.exists('fcs_t')){dir.create("fcs_t", showWarnings = FALSE)}
setwd("fcs_t")

#Save flowframes wihtin flowset as fcs files using the flowCore package
write.flowSet(fs_fda, outdir='fcs_t', filename = sampleNames(fs_fda))


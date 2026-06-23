#########################################################################################################################
# This code is designed to generate hurdle models to fit population equivalent density and imidacloprid concentration   #
# data described in Liu et al. (sub judice), Sources and Mitigation of Nationwide Parasiticide Water Pollution.         #
# This code assumes that mean concentrations are being fit. It is straightforward to change 'mean' to 'median' in the   #
# code below. This code includes plotting scripts used to produce figures in the paper and its supporting information.  #
#                                                                                                                       #
# contact will.pearse@imperial.ac.uk or gareth.roberts@imperial.ac.uk for further information                           #
#                                                                                                                       #
# June 2026                                                                                                             #
#########################################################################################################################

#############################
# Headers                   #
#############################
# Packages
library(dplyr)  
library(ggplot2)
library(ggrepel)
library(rstanarm)
library(brms)

# Log conversion
convert.logs <- function(x) return(x/log(10))

#############################
# Data loading and prep     #
#############################
# Load pre-processed data and subset
WWpopulation <- readRDS("./WWpopulation.rds")
allSamples_2019OW <- readRDS("./allSamples1909RBD.rds")
pe <- WWpopulation |> select(Sample_Site_ID, PerKM2PopServed)
names(pe)[2] <- "population"

# Imidacloprid in Running Waters 
standing_waters_names <- unique(allSamples_2019OW$Sample_Site_ID[allSamples_2019OW$SMC_DESC == "POND / LAKE / RESERVOIR WATER"])
running_waters_samples_since2019 <- allSamples_2019OW |> filter(! Sample_Site_ID %in% standing_waters_names)

#################################################################################
# Preparing averaged concentration data.                                        #
#                                                                               #
# Uncomment the following blocks of code to enable the hurdle models to be      # 
# run for each of the approaches used to substitute or remove concentrations    #
# < LOD prior to running the models. Or you can run just one model, like        #
# the uncommented example below that assumes measured concentrations < LOD      #
# are replaced with LOD = 0.001 mg /L.                                          #
#################################################################################

# remove < LOD prior to averaging 
# g.data <- left_join(pe, running_waters_samples_since2019, by = "Sample_Site_ID")
# g.data <- g.data |> filter(Sample_Site_ID != "50022")
# g.data <- g.data[!is.na(g.data$Concentration),]
# stats_rm <- g.data |>
#   group_by(Sample_Site_ID) |>
#   summarise(
#     count = n(),
#     count.pos = sum(Detection),
#     detperc = (count.pos/count)*100,
#     mean_conc = mean(Concentration),
#     se_conc = sd(Concentration)/sqrt(length(Concentration)),
#     median_conc = median(Concentration),
#     mad_conc = mad(Concentration),
#     max_conc = max(Concentration)
#   ) |> 
#   filter(Sample_Site_ID %in% pe$Sample_Site_ID) |>
#   select(Sample_Site_ID, mean_conc, se_conc, median_conc, mad_conc)
# srm.data <- left_join(pe, stats_rm, by = "Sample_Site_ID")
# srm.data <- srm.data |> filter(Sample_Site_ID != "50022")
# srm.data <- srm.data[!is.na(srm.data$mean_conc),]
# srm.data$log.pop <- log10(srm.data$population)
# srm.data$log.pop[srm.data$population==0] <- 0
# srm.data$b.pop <- as.numeric(srm.data$log.pop>0)
# # fit model
# trm.data <- srm.data[srm.data$mean_conc > 0,]
# model_rm <- brm(log10(mean_conc) ~ log.pop, data=trm.data)

# replace < LOD with LOD.
stats_lod <- running_waters_samples_since2019 |>
  mutate(Concentration = ifelse(is.na(Concentration), 0.001, Concentration)) |>
  group_by(Sample_Site_ID) |>
  summarise(
    count = n(),
    count.pos = sum(Detection),
    detperc = (count.pos/count)*100,
    mean_conc = mean(Concentration),
    se_conc = sd(Concentration)/sqrt(length(Concentration)),
    median_conc = median(Concentration),
    mad_conc = mad(Concentration),
    max_conc = max(Concentration)
  ) |>
  filter(Sample_Site_ID %in% pe$Sample_Site_ID) |>
  select(Sample_Site_ID, mean_conc, se_conc, median_conc, mad_conc)
slod.data <- left_join(pe, stats_lod, by = "Sample_Site_ID")
slod.data <- slod.data |> filter(Sample_Site_ID != "50022")
slod.data$log.pop <- log10(slod.data$population)
slod.data$log.pop[slod.data$population==0] <- 0
slod.data$b.pop <- as.numeric(slod.data$log.pop>0)
# fit model
tlod.data <- slod.data[slod.data$mean_conc > 0,]
model_repLOD <- brm(log10(mean_conc) ~ log.pop, data=tlod.data)
 
# replace < LOD with 1/2 LOD.
# stats_05lod <- running_waters_samples_since2019 |>
#   mutate(Concentration = ifelse(is.na(Concentration), 0.0005, Concentration)) |>
#   group_by(Sample_Site_ID) |>
#   summarise(
#     count = n(),
#     count.pos = sum(Detection),
#     detperc = (count.pos/count)*100,
#     mean_conc = mean(Concentration),
#     se_conc = sd(Concentration)/sqrt(length(Concentration)),
#     median_conc = median(Concentration),
#     mad_conc = mad(Concentration),
#     max_conc = max(Concentration)
#   ) |>
#   filter(Sample_Site_ID %in% pe$Sample_Site_ID) |>
#   select(Sample_Site_ID, mean_conc, se_conc, median_conc, mad_conc)
# s05.data <- left_join(pe, stats_05lod, by = "Sample_Site_ID")
# s05.data <- s05.data |> filter(Sample_Site_ID != "50022")
# s05.data <- s05.data[!is.na(s05.data$mean_conc),]
# s05.data$log.pop <- log10(s05.data$population)
# s05.data$log.pop[s05.data$population==0] <- 0
# s05.data$b.pop <- as.numeric(s05.data$log.pop>0)
# # fit model
# t05.data <- s05.data[s05.data$mean_conc > 0,]
# model_rep05LOD <- brm(log10(mean_conc) ~ log.pop, data=t05.data)

# # replace < LOD with 0.
# stats_0 <- running_waters_samples_since2019 |>
#   mutate(Concentration = ifelse(is.na(Concentration), 0, Concentration)) |>
#   group_by(Sample_Site_ID) |>
#   summarise(
#     count = n(),
#     count.pos = sum(Detection),
#     detperc = (count.pos/count)*100,
#     mean_conc = mean(Concentration),
#     se_conc = sd(Concentration)/sqrt(length(Concentration)),
#     median_conc = median(Concentration),
#     mad_conc = mad(Concentration),
#     max_conc = max(Concentration)
#   ) |>
#   filter(Sample_Site_ID %in% pe$Sample_Site_ID) |>
#   select(Sample_Site_ID, mean_conc, se_conc, median_conc, mad_conc)
# s0.data <- left_join(pe, stats_0, by = "Sample_Site_ID")
# s0.data <- s0.data |> filter(Sample_Site_ID != "50022")
# s0.data$log.pop <- log10(s0.data$population)
# s0.data$log.pop[s0.data$population==0] <- 0
# s0.data$b.pop <- as.numeric(s0.data$log.pop>0)
# # fit model
# t0.data <- s0.data[s0.data$mean_conc > 0,]
# model_rep0 <- brm(log10(mean_conc) ~ log.pop, data=t0.data)

#############################
# Generate plots ############
#############################


pdf(file="./subs-plot.pdf")
# Prepare a 'pretty' dataset for plotting (p.data)
# ...making sure the 0-0 points are well-spaced for ease of viewing and are at the right point on the axes
# ...and working with log-10 concentrations of imid
prm.data <- srm.data
plod.data <- slod.data
p05.data <- s05.data
p0.data <- s0.data

no.imid <- -4
no.pop <- .2

prm.data$log.pop[prm.data$log.pop==0] <- no.pop# + seq(-.2, .2, length.out=sum(prm.data$log.pop==0))
prm.data$mean_conc <- log10(prm.data$mean_conc)
prm.data$mean_conc[!is.finite(prm.data$mean_conc)] <- no.imid# + seq(-.2, .2, length.out=sum(!is.finite(prm.data$mean_conc)))

plod.data$log.pop[plod.data$log.pop==0] <- no.pop# + seq(-.2, .2, length.out=sum(plod.data$log.pop==0))
plod.data$mean_conc <- log10(plod.data$mean_conc)
plod.data$mean_conc[!is.finite(plod.data$mean_conc)] <- no.imid# + seq(-.2, .2, length.out=sum(!is.finite(plod.data$mean_conc)))

p05.data$log.pop[p05.data$log.pop==0] <- no.pop# + seq(-.2, .2, length.out=sum(p05.data$log.pop==0))
p05.data$mean_conc <- log10(p05.data$mean_conc)
p05.data$mean_conc[!is.finite(p05.data$mean_conc)] <- no.imid# + seq(-.2, .2, length.out=sum(!is.finite(p05.data$mean_conc)))

p0.data$log.pop[p0.data$log.pop==0] <- no.pop# + seq(-.2, .2, length.out=sum(p0.data$log.pop==0))
p0.data$mean_conc <- log10(p0.data$mean_conc)
p0.data$mean_conc[!is.finite(p0.data$mean_conc)] <- no.imid# + seq(-.2, .2, length.out=sum(!is.finite(p0.data$mean_conc)))

# Numbers for legend
sum(s0.data$mean_conc==0)#27
sum(s0.data$population==0)#18
sum(s0.data$population==0 & s.data$mean_conc==0)#10

# Setup plot with site-averages with axes to demark the zeroes
with(plod.data, plot(mean_conc ~ log.pop, axes=FALSE, pch=20,
                  ylab=expression(Log[10](Imidacloprid~concentration,~mu*g/L)), xlab=expression(Log[10](population~equivalent~density)),
                  xlim=c(no.pop,3.6), ylim=c(no.imid,-0.5), type="n"))
axis(1, at=no.pop, labels="none")
axis(1, at=c(.5,1:3,3.5))
axis(2, at=no.imid, labels="no\ndetection")
axis(2, at=-c(3.5, 3:1, 0.5))

# Add threshold labels
abline(h = log10(0.2), col = "red", lty=2, lwd=1.5) #acute
abline(h = log10(0.035), col = "orange", lty=2, lwd=1.5) #chronic
abline(h = log10(0.0068), col = "#FFCC33", lty=2, lwd=1.5) #PNEC
text(0.7,log10(0.2)+0.07, "Acute toxicity", col="red")
text(0.7,log10(0.035)+0.07, "Chronic toxicity", col="orange")
text(0.7,log10(0.0068)+0.07, "Lowest PNEC", col="#FFCC33")

# Add the site-averaged trend in red
line <- fixef(model_repLOD)["Intercept","Estimate"] + fixef(model_repLOD)["log.pop","Estimate"]*c(1,3.6)
lines(c(1,3.6), line, col="red", lwd=2)

line <- fixef(model_rep05LOD)["Intercept","Estimate"] + fixef(model_rep05LOD)["log.pop","Estimate"]*c(1,3.6)
lines(c(1,3.6), line, col="blue", lwd=2)

line <- fixef(model_rep0)["Intercept","Estimate"] + fixef(model_rep0)["log.pop","Estimate"]*c(1,3.6)
lines(c(1,3.6), line, col="darkgreen", lwd=2)

line <- fixef(model_rm)["Intercept","Estimate"] + fixef(model_rm)["log.pop","Estimate"]*c(1,3.6)
lines(c(1,3.6), line, col="black", lwd=2)

# Add the data in so the lines don't obscure them
with(plod.data, points(mean_conc ~ log.pop, pch=20, col="pink"))
with(p05.data, points(mean_conc ~ log.pop, pch=20, col="lightblue"))
with(p0.data, points(mean_conc ~ log.pop, pch=20, col="lightgreen"))
with(prm.data, points(mean_conc ~ log.pop, pch=20, col="darkgrey"))

# Add values at pop = 0
points(no.pop, fixef(model_repLOD)["Intercept","Estimate"], pch=20, cex=1.5, col="red")
points(no.pop, fixef(model_rep05LOD)["Intercept","Estimate"], pch=20, cex=1.5, col="blue")
points(no.pop, fixef(model_rep0)["Intercept","Estimate"], pch=20, cex=1.5, col="darkgreen")
points(no.pop, fixef(model_rm)["Intercept","Estimate"], pch=20, cex=1.5, col="black") 

dev.off()

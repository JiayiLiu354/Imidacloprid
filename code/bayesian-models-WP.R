#############################
# Headers ###################
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
# Data loading and prep #####
#############################
# Load pre-processed data and subset
WWpopulation <- readRDS("./WWpopulation.rds")
allSamples_2019OW <- readRDS("./allSamples1909RBD.rds")
pe <- WWpopulation |> select(Sample_Site_ID, PerKM2PopServed)
names(pe)[2] <- "population"

# Imidacloprid in Running Waters 
standing_waters_names <- unique(allSamples_2019OW$Sample_Site_ID[allSamples_2019OW$SMC_DESC == "POND / LAKE / RESERVOIR WATER"])
running_waters_samples_since2019 <- allSamples_2019OW |>
  filter(! Sample_Site_ID %in% standing_waters_names)

# Prepare data without averages across sites
a.data <- left_join(pe, running_waters_samples_since2019, by = "Sample_Site_ID")
a.data <- a.data |> filter(Sample_Site_ID != "50022")

# Preparing averaged data
stats_0 <- running_waters_samples_since2019 |>
  mutate(Concentration = ifelse(is.na(Concentration), 0, Concentration)) |>
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

# Merge data and make site-averaged dataset (s.data)
s.data <- left_join(pe, stats_0, by = "Sample_Site_ID")
s.data <- s.data |> filter(Sample_Site_ID != "50022")
s.data$log.pop <- log10(s.data$population)
s.data$log.pop[s.data$population==0] <- 0
s.data$b.pop <- as.numeric(s.data$log.pop>0)

# Prepare data with 'all' site data (a.data)
a.data$log.pop <- log10(a.data$population)
a.data$log.pop[!is.finite(a.data$log.pop)] <- 0
a.data$b.pop <- as.numeric(a.data$log.pop>0)
a.data$Concentration[is.na(a.data$Concentration)] <- 0
a.data <- a.data[,c("Concentration","log.pop","b.pop","population","year","Sample_datetime","Sample_Site_ID")]
a.data$day <- as.Date(a.data$Sample_datetime, format="%m/%d/%Y")
a.data$day <- as.numeric(a.data$day - as.Date("2019-01-02"))

#############################
# Fit models ################
#############################
# Remember that the hurdle model is for whether something is 0, so it's a model of whether we *don't* detect it
model <- brm(bf(mean_conc ~ b.pop + log.pop, hu ~ b.pop + log.pop), data = s.data, family = hurdle_lognormal())
car.model <- brm(bf(Concentration ~ log.pop + b.pop + (1|Sample_Site_ID), autocor=~cor_car(day|Sample_Site_ID),
                    hu ~ log.pop + b.pop + (1|Sample_Site_ID)), data = a.data, family = hurdle_lognormal(), iter=4000)
#saveRDS(model, "model.RDS")
#saveRDS(car.model, "car.model.RDS")

library(modelsummary)
modelsummary(list(model, car.model), output="models-for-formatting.docx")

#############################
# Generate plots ############
#############################
# Load the RDS files if needed
#model <- readRDS("model.RDS")
#car.model <- readRDS("car.model.RDS")

pdf("averaged-plot.pdf")
# Prepare a 'pretty' dataset for plotting (p.data)
# ...making sure the 0-0 points are well-spaced for ease of viewing and are at the right point on the axes
# ...and working with log-10 concentrations of imid
p.data <- s.data
no.imid <- -4
no.pop <- .2
p.data$log.pop[p.data$log.pop==0] <- no.pop# + seq(-.2, .2, length.out=sum(p.data$log.pop==0))
p.data$mean_conc <- log10(p.data$mean_conc)
p.data$mean_conc[!is.finite(p.data$mean_conc)] <- no.imid# + seq(-.2, .2, length.out=sum(!is.finite(p.data$mean_conc)))

# Numbers for legend
sum(s.data$mean_conc==0)#27
sum(s.data$population==0)#18
sum(s.data$population==0 & s.data$mean_conc==0)#10

# Setup plot with site-averages with axes to demark the zeroes
with(p.data, plot(mean_conc ~ log.pop, axes=FALSE, pch=20,
                  ylab=expression(Log[10](Imidacloprid~concentration~(mu*g/L))), xlab=expression(Log[10](population~equivalent~density)),
                  xlim=c(no.pop,3.5), ylim=c(no.imid,-1), type="n"))
axis(1, at=no.pop, labels="none")
axis(1, at=c(.5,1:3,3.5))
axis(2, at=no.imid, labels="no\ndetection")
axis(2, at=-c(3.5, 3:1))
with(p.data, text(no.pop, no.imid+.1, paste0("(n=",sum(mean_conc==no.imid & log.pop==no.pop),")"), xpd=TRUE))

# Add the site-averaged trend in red (converting the logarithm base)
line <- fixef(model)["Intercept","Estimate"]+fixef(model)["b.pop","Estimate"] + fixef(model)["log.pop","Estimate"]*c(1,3.5)
lines(c(1,3.5), convert.logs(line), col="red", lwd=2)
points(no.pop, convert.logs(fixef(model)["Intercept","Estimate"]), pch=20, cex=1.5, col="red")
#lines(c(.4,.6), rep(convert.logs(fixef(model)["Intercept","Estimate"],2)), col="red", lwd=2)

# Add the CAR model in blue (converting the logarithm base)
line <- fixef(car.model)["Intercept","Estimate"]+fixef(car.model)["b.pop","Estimate"] + fixef(car.model)["log.pop","Estimate"]*c(1,3.5)
lines(c(1,3.5), convert.logs(line), col="blue", lwd=2)
points(no.pop, convert.logs(fixef(car.model)["Intercept","Estimate"]), pch=20, cex=1.5, col="blue")
#lines(c(.4,.6), rep(convert.logs(fixef(car.model)["Intercept","Estimate"],2)), col="blue", lwd=2)

# Add the data in so the lines don't obscure them
with(p.data, points(mean_conc ~ log.pop, pch=20))
dev.off()

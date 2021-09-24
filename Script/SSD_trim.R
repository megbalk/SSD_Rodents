# Trend of sexual dimorphism at the species level across Rodentia
# Meghan A. Balk
# meghan.balk@gmail.com

##load packages----
require(tidyverse)
require(nlme)
require(dplyr)
require(ggplot2)
require(reshape2)
require(plyr)
require(stringr)
require(utils)
require(taxize)
library(picante)
library(ape)
library(adephylo)
library(ade4)
library(phylobase)
library(geiger)
library(phytools)
library(AICcmodavg)
library(caper)
require(lmodel2)
require(visreg)
require(car)
require(ggtree)

##load data----

options(stringsAsFactors = FALSE)

data <- read.csv("log.normal.flagged.data.csv", header = TRUE)
tree <- read.nexus("https://data.cyverse.org/dav-anon/iplant/home/rwalls/FuTRES_data/Projects/SFritz.tre")

##manipulate data----

status <- c("juvenile.quant", "juvenile.sd", "outlier.sd", "outlier.quant", "juvenile.log.sd", "outlier.log.sd", "outlier")

df <- data[data$lifeStage == "adult" & 
          (data$sex == "male" | data$sex == "female") & 
          (data$measurementType == "{body mass}" | 
           data$measurementType == "{body length}" |
           data$measurementType == "{tail length}") &
          !(data$measurementStatus %in% status),]
df <- df %>%
  select(-(X.4, X.3, X.2, X.1, X, traits, normality, logMeasurementValue, lowerLimit, upperLImit, lowerLimitMethod, upperLimitMethod, meanValue, sdValue, meanMethod, sdMethod))

data.short <- df %>%
  pivot_wider(names_from = measurementType, values_from = measurementValue) %>%
  as.data.frame()
colnames(data.short)

data.short$`{head body length}` <- data.short$`{body length}` - data.short$`{tail length}`

df <- df %>%
  drop_na(measurementValue)

write.csv(df, "ssd.df.csv")

sp.table <- df %>%
  dplyr::group_by(scientificName, sex, measurementType) %>%
  tally() %>%
  filter(n >= 10) %>%
  mutate(sex.measType = paste(sex, ".", measurementType, sep = "")) %>%
  dplyr::select(-measurementType) %>%
  spread(key = sex.measType, value = n) %>%
  dplyr::group_by(scientificName) %>%
  dplyr::select(-sex) %>%
  dplyr::summarise(female.mass = sum(`female.{body mass}`, na.rm = TRUE),
                   male.mass = sum(`male.{body mass}`, na.rm = TRUE),
                   female.body.length = sum(`female.{body length}`, na.rm = TRUE),
                   male.body.length = sum(`male.{body length}`, na.rm = TRUE),
                   tots = isTRUE(female.mass >= 10 & male.mass >= 10) + isTRUE(female.body.length >= 10 & male.body.length >= 10)) %>%
  filter(tots > 0) %>%
  as.data.frame()
#if a species doesn't have mass for one sex, it tends not to have total length either

write.csv(sp.table, "sp.table.csv")

sp <- sp.table$scientificName

#order

sp_out <- tax_name(sp, get = 'order', db = 'ncbi')

write.csv(sp_out, "sp.taxonomy.csv")


#2-tailed t test to see if different and get direction
#add columns "significant" (T/F) and "direction" (M/F)

rodentia <- sp_out %>%
  filter(order == "Rodentia")

rodent.sp <- rodentia$query

df.rodent <- df[df$scientificName %in% rodent.sp,]

write.csv(df.rodent, "rodents.csv")

rodent.table <- df.rodent %>%
  dplyr::group_by(scientificName, sex, measurementType) %>%
  tally() %>%
  filter(n >= 10) %>%
  mutate(sex.measType = paste(sex, ".", measurementType, sep = "")) %>%
  dplyr::select(-measurementType) %>%
  spread(key = sex.measType, value = n) %>%
  dplyr::group_by(scientificName) %>%
  dplyr::select(-sex) %>%
  dplyr::summarise(female.mass = sum(`female.{body mass}`, na.rm = TRUE),
                   male.mass = sum(`male.{body mass}`, na.rm = TRUE),
                   female.total.length = sum(`female.{body length}`, na.rm = TRUE),
                   male.total.length = sum(`male.{body length}`, na.rm = TRUE),
                   tots = isTRUE(female.mass >= 10 & male.mass >= 10) + isTRUE(female.total.length >= 10 & male.total.length >= 10)) %>%
  filter(tots > 0) %>%
  as.data.frame()

write.csv(rodent.table, "rodent.table.csv")

##PeMa test----
#two sample test; two sided
pema <- df.rodent[df.rodent$scientificName == "Peromyscus maniculatus",]
tt <- t.test(pema$measurementValue[pema$measurementType == "mass" & pema$sex == "male"],
             pema$measurementValue[pema$measurementType == "mass" & pema$sex == "female"],
             alternative = "two.sided", paired = FALSE)
#are males greater than females?
tg <- t.test(pema$measurementValue[pema$measurementType == "mass" & pema$sex == "male"],
             pema$measurementValue[pema$measurementType == "mass" & pema$sex == "female"],
             alternative = "greater", paired = FALSE)
#are males smaller than females?
tl <- t.test(pema$measurementValue[pema$measurementType == "mass" & pema$sex == "male"],
             pema$measurementValue[pema$measurementType == "mass" & pema$sex == "female"],
             alternative = "less", paired = FALSE)
#p.value
tt$p.value #if TRUE tt$p.value <= 0.05, put "TRUE" in sig
#male mean
tt$estimate[[1]]
#female mean
tt$estimate[[2]]

#direction: only go into this if tt$p.value <= 0.05
tg$p.value #if TRUE tg$p.value <= 0.05, put "M" in direction
tl$p.value #if TRUE tl$p.value <= 0.05, put "F" in direction

##try on larger subset
columns <- c("sp.name",
            "n.male",
            "n.female",
            "tt.p.value",
            "df",
            "male.mean",
            "female.mean",
            "mean.diff",
            "sig",
            "tg.p.value",
            "tl.p.value",
            "direction")
test.sp <- rodent.sp[1:3]
test.t.test <- data.frame(test.sp, matrix(ncol = length(columns))) #should have 3 rows 
colnames(test.t.test) = columns

for(i in 1:length(test.sp)){
  test.t.test$sp.name[i] <- test.sp[i]
  test.t.test$n.male[i] <- nrow(df.rodent[df.rodent$scientificName == test.sp[i] & 
                                          df.rodent$measurementType == "mass" & 
                                          df.rodent$sex == "male",])
  test.t.test$n.female[i] <- nrow(df.rodent[df.rodent$scientificName == test.sp[i] & 
                                            df.rodent$measurementType == "mass" & 
                                            df.rodent$sex == "female",])
  mm <-df.rodent$measurementValue[df.rodent$scientificName == test.sp[i] & 
                                    df.rodent$measurementType == "mass" & 
                                    df.rodent$sex == "male"]
  ff <- df.rodent$measurementValue[df.rodent$scientificName == test.sp[i] & 
                                     df.rodent$measurementType == "mass" & 
                                     df.rodent$sex == "female"]
  tt <- t.test(mm, ff, alternative = "two.sided", paired = FALSE)
  test.t.test$tt.p.value[i] = tt$p.value
  test.t.test$df[i] <- tt$parameter[[1]]
  test.t.test$male.mean[i] = tt$estimate[[1]]
  test.t.test$female.mean[i] = tt$estimate[[2]]
  test.t.test$mean.diff[i] = tt$estimate[[1]] - tt$estimate[[2]]
  if(isTRUE(tt$p.value <= 0.05)){
    test.t.test$sig[i] = "TRUE"
  }
  else{
    test.t.test$sig[i] = "FALSE"
  }
  tg <- t.test(mm,ff, alternative = "greater", paired = FALSE)
  test.t.test$tg.p.value[i] = tg$p.value
  tl <- t.test(mm, ff, alternative = "less", paired = FALSE)
  test.t.test$tl.p.value[i] = tl$p.value
  if(isTRUE(tg$p.value <= 0.05)){
    test.t.test$direction[i] = "M"
  }
  else if(isTRUE(tl$p.value <= 0.05)){
    test.t.test$direction[i] = "F"
  }
  else{
    test.t.test$direction[i] = ""
  }
}

##t-test for all Rodentia----
#two sample test; two sided
#capture sample size for M and F
#extract sig level
#extract direction
#extract df
#means for M and F
#difference in means
columns <-c("sp.name",
            "n.male",
            "n.female",
            "tt.p.value",
            "df",
            "male.mean",
            "female.mean",
            "mean.diff",
            "sig",
            "tg.p.value",
            "tl.p.value",
            "direction")

rodent.t.test.mass <- data.frame(rodent.sp, matrix(ncol = length(columns))) #should have 188 rows 
colnames(rodent.t.test.mass) = columns

for(i in 1:length(rodent.sp)){
  rodent.t.test.mass$sp.name[i] <- rodent.sp[i]
  rodent.t.test.mass$n.male[i] <- nrow(df.rodent[df.rodent$scientificName == rodent.sp[i] & df.rodent$measurementType == "mass" & df.rodent$sex == "male",])
  rodent.t.test.mass$n.female[i] <- nrow(df.rodent[df.rodent$scientificName == rodent.sp[i] & df.rodent$measurementType == "mass" & df.rodent$sex == "female",])
  mm <- df.rodent$measurementValue[df.rodent$scientificName == rodent.sp[i] & df.rodent$measurementType == "mass" & df.rodent$sex == "male"]
  ff <- df.rodent$measurementValue[df.rodent$scientificName == rodent.sp[i] & df.rodent$measurementType == "mass" & df.rodent$sex == "female"]
  if(isTRUE(length(mm) < 10 & length(mm) < 10)){ #need to see if less that 10 of non NA 
    next
  }
  tt <- t.test(mm, ff, alternative = "two.sided", paired = FALSE)
  rodent.t.test.mass$tt.p.value[i] = tt$p.value
  rodent.t.test.mass$df[i] <- tt$parameter[[1]]
  rodent.t.test.mass$male.mean[i] = tt$estimate[[1]]
  rodent.t.test.mass$female.mean[i] = tt$estimate[[2]]
  rodent.t.test.mass$mean.diff[i] = tt$estimate[[1]] - tt$estimate[[2]]
  if(isTRUE(tt$p.value <= 0.05)){
    rodent.t.test.mass$sig[i] = "TRUE"
  }
  else{
    rodent.t.test.mass$sig[i] = "FALSE"
  }
  tg <- t.test(df.rodent$measurementValue[df.rodent$scientificName == rodent.sp[i] & df.rodent$measurementType == "mass" & df.rodent$sex == "male"],
               df.rodent$measurementValue[df.rodent$scientificName == rodent.sp[i] & df.rodent$measurementType == "mass" & df.rodent$sex == "female"],
               alternative = "greater", paired = FALSE)
  rodent.t.test.mass$tg.p.value[i] = tg$p.value
  tl <- t.test(df.rodent$measurementValue[df.rodent$scientificName == rodent.sp[i] & df.rodent$measurementType == "mass" & df.rodent$sex == "male"],
               df.rodent$measurementValue[df.rodent$scientificName == rodent.sp[i] & df.rodent$measurementType == "mass" & df.rodent$sex == "female"],
               alternative = "less", paired = FALSE)
  rodent.t.test.mass$tl.p.value[i] = tl$p.value
  if(isTRUE(tg$p.value <= 0.05)){
    rodent.t.test.mass$direction[i] = "M"
  }
  else if(isTRUE(tl$p.value <= 0.05)){
    rodent.t.test.mass$direction[i] = "F"
  }
  else{
    rodent.t.test.mass$direction[i] = ""
  }
}

write.csv(rodent.t.test.mass, "rodent.t.test.mass.csv")
rodent.t.test.mass <- read.csv("https://data.cyverse.org/dav-anon/iplant/home/rwalls/FuTRES_data/Projects/SSD_Rodents/rodent.t.test.mass.csv", header = TRUE)

rodent.t.test.length <- data.frame(rodent.sp, matrix(ncol = length(columns))) #should have 188 rows 
colnames(rodent.t.test.length) = columns

for(i in 1:length(rodent.sp)){
  rodent.t.test.length$sp.name[i] <- rodent.sp[i]
  rodent.t.test.length$n.male[i] <- nrow(df.rodent[df.rodent$scientificName == rodent.sp[i] & df.rodent$measurementType == "total.length" & df.rodent$sex == "male",])
  rodent.t.test.length$n.female[i] <- nrow(df.rodent[df.rodent$scientificName == rodent.sp[i] & df.rodent$measurementType == "total.length" & df.rodent$sex == "female",])
  mm <- df.rodent$measurementValue[df.rodent$scientificName == rodent.sp[i] & df.rodent$measurementType == "total.length" & df.rodent$sex == "male"]
  ff <- df.rodent$measurementValue[df.rodent$scientificName == rodent.sp[i] & df.rodent$measurementType == "total.length" & df.rodent$sex == "female"]
  if(isTRUE(length(mm) < 10 & length(mm) < 10)){ #need to see if less that 10 of non NA 
    next
  }
  tt <- t.test(mm, ff, alternative = "two.sided", paired = FALSE)
  rodent.t.test.length$tt.p.value[i] = tt$p.value
  rodent.t.test.length$df[i] <- tt$parameter[[1]]
  rodent.t.test.length$male.mean[i] = tt$estimate[[1]]
  rodent.t.test.length$female.mean[i] = tt$estimate[[2]]
  rodent.t.test.length$mean.diff[i] = tt$estimate[[1]] - tt$estimate[[2]]
  if(isTRUE(tt$p.value <= 0.05)){
    rodent.t.test.length$sig[i] = "TRUE"
  }
  else{
    rodent.t.test.length$sig[i] = "FALSE"
  }
  tg <- t.test(df.rodent$measurementValue[df.rodent$scientificName == rodent.sp[i] & df.rodent$measurementType == "total.length" & df.rodent$sex == "male"],
               df.rodent$measurementValue[df.rodent$scientificName == rodent.sp[i] & df.rodent$measurementType == "total.length" & df.rodent$sex == "female"],
               alternative = "greater", paired = FALSE)
  rodent.t.test.length$tg.p.value[i] = tg$p.value
  tl <- t.test(df.rodent$measurementValue[df.rodent$scientificName == rodent.sp[i] & df.rodent$measurementType == "total.length" & df.rodent$sex == "male"],
               df.rodent$measurementValue[df.rodent$scientificName == rodent.sp[i] & df.rodent$measurementType == "total.length" & df.rodent$sex == "female"],
               alternative = "less", paired = FALSE)
  rodent.t.test.length$tl.p.value[i] = tl$p.value
  if(isTRUE(tg$p.value <= 0.05)){
    rodent.t.test.length$direction[i] = "M"
  }
  else if(isTRUE(tl$p.value <= 0.05)){
    rodent.t.test.length$direction[i] = "F"
  }
  else{
    rodent.t.test.length$direction[i] = ""
  }
}

write.csv(rodent.t.test.length, "rodent.t.test.length.csv")

rodent.t.test.length <- read.csv("https://data.cyverse.org/dav-anon/iplant/home/rwalls/FuTRES_data/Projects/SSD_Rodents/rodent.t.test.length.csv", header = TRUE)

##plot species-level relationships----

##mass
p <- ggplot(data = rodent.t.test.mass) +
     geom_point(aes(x = log10(female.mean), y = log10(male.mean))) +
     geom_smooth(aes(x = log10(female.mean), y = log10(male.mean))) +
     ggtitle("Rodentia") +
     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
           panel.background = element_blank(), axis.line = element_line(colour = "black"),
           legend.position = "none") +
     scale_y_continuous(name = expression(log[10]~Male~Body~Mass~(g))) +
     scale_x_continuous(name = expression(log[10]~Female~Body~Mass~(g)))
ggsave(p, file=paste0("plot_log_Rodentia_mass.png"), width = 14, height = 10, units = "cm")

mass.model <- summary(lm(log10(rodent.t.test.mass$male.mean) ~ log10(rodent.t.test.mass$female.mean)))

mass.stats <- data.frame(scientificName = "Rodentia",
                         comparison = ("male mass/female mass"),
                         intercept = mass.model$coefficients[[1]],
                         slope = mass.model$coefficients[[2]],
                         resid.std.err = mass.model$sigma,
                         df = max(mass.model$df),
                         std.err.slope =  mass.model$coefficients[4],
                         std.err.intercept = mass.model$coefficients[3],
                         r.squared = mass.model$r.squared,
                         p.value = mass.model$coefficients[,4][[2]],
                         sample.size = nrow(rodent.t.test.mass),
                         male.sample.size = sum(rodent.t.test.mass$n.male),
                         female.sample.size = sum(rodent.t.test.mass$n.female))

##length
p <- ggplot(data = rodent.t.test.length) +
  geom_point(aes(x = log10(female.mean), y = log10(male.mean))) +
  geom_smooth(aes(x = log10(female.mean), y = log10(male.mean))) +
  ggtitle("Rodentia") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none") +
  scale_y_continuous(name = expression(log[10]~Male~Total~Length~(mm))) +
  scale_x_continuous(name = expression(log[10]~Female~Total~Length~(mm)))
ggsave(p, file=paste0("plot_log_Rodentia_length.png"), width = 14, height = 10, units = "cm")

length.model <- summary(lm(log10(rodent.t.test.length$male.mean) ~ log10(rodent.t.test.length$female.mean)))

length.stats <- data.frame(scientificName = "Rodentia",
                           comparison = ("male length/female length"),
                           intercept = length.model$coefficients[[1]],
                           slope = length.model$coefficients[[2]],
                           resid.std.err = length.model$sigma,
                           df = max(length.model$df),
                           std.err.slope =  length.model$coefficients[4],
                           std.err.intercept = length.model$coefficients[3],
                           r.squared = length.model$r.squared,
                           p.value = length.model$coefficients[,4][[2]],
                           sample.size = nrow(rodent.t.test.length),
                           male.sample.size = sum(rodent.t.test.length$n.male),
                           female.sample.size = sum(rodent.t.test.length$n.female))

rodentia.stats <- rbind(mass.stats, length.stats)

write.csv(rodentia.stats, "rodentia.model.stats.csv")

##PGLS----

##mass
# make trees match
#replace spaces with underscores
rodent.mass <- rodent.t.test.mass %>%
  drop_na(male.mean, female.mean)
rodent.mass$sp.name <- gsub('([[:punct:]])|\\s+', '_', rodent.mass$sp.name)
rownames(rodent.mass) <- rodent.mass$sp.name
# to drop the species from the tree that are not in the dataset
mamm_speTree <- drop.tip(tree, tree$tip.label[!is.element(tree$tip.label, rownames(rodent.mass))])
# to match the species in the tree and in the dataset
rodent.mass <- rodent.mass[match(mamm_speTree$tip.label, rodent.mass$sp.name),]
names(rodent.mass)
length(rodent.mass$sp.name)
length(mamm_speTree$tip.label)

mamm_ssd_pgls <- gls(log10(male.mean) ~ log10(female.mean),
                     data = rodent.mass,
                     correlation = corPagel(0, phy = mamm_speTree, fixed = FALSE),
                     method = "ML") ### What does REML mean??


##figure of tree with t-test results

#col_color <- c("gray", "deepsky blue4", "lightsalmon4")
rodent.mass$color <- "gray"
rodent.mass$color[rodent.mass$direct == "F"] <- "deepskyblue4"
rodent.mass$color[rodent.mass$direct == "M"] <- "lightsalmon4"

rodent.direction <- rodent.mass %>%
  dplyr::select(sp.name, direction)

p <- ggtree(mamm_speTree, 
            aes(color = direction),
            branch.length='none', 
            layout='circular') %<+% 
  rodent.direction +  

p %<+% rodent.direction +
  geom_(aes(color = direction,
            color = "black", #label font
            geom = "label", #labels not text
            label.padding = unit(0.15, "lines"), #amount of padding around labels
            label.size = 0)) + 
  scale_color_manual(values=c("gray", "deepskyblue4","lightsalmon4"))
                
                
                
+ #size of label border
  theme(legend.position = c(0.5, 0.2),
        legend.title = element_text("Direction of SSD"),
        lengend.key = element_blank())




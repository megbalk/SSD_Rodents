# Meghan A. Balk
# balkm@email.arizona.edu

##Load packages----
require(tidyverse)
require(nlme)
require(dplyr)
require(ggplot2)
require(reshape2)
require(plyr)
require(stringr)
require(OutlierDetection)
require(utils)
require(taxize)

##Load data----

options(stringsAsFactors = FALSE)

data <- read.csv("https://de.cyverse.org/dl/d/7B6941EE-1FED-4AAE-9D8A-8AC76DB0AC98/data.all.csv", header = TRUE, stringsAsFactors = FALSE)

data$sex[data$sex == "Male" | data$sex == "M" | data$sex == "males" | data$sex == "males males" | data$sex == "male male" | data$sex == "Male "] <- "male"
data$sex[data$sex == "Female" | data$sex == "F" | data$sex == "females" | data$sex == "females females" | data$sex == "female female"] <- "female"

df <- data[data$lifeStage == "Adult" & (data$sex == "male" | data$sex == "female") & (data$measurementType == "mass" | data$measurementType == "total.length"),]

df <- df %>%
  drop_na(measurementValue)

write.csv(df, "ssd.df.csv")

df <- read.csv("https://data.cyverse.org/dav-anon/iplant/home/rwalls/FuTRES_data/SSD/ssd.df.csv", header = TRUE, stringsAsFactors = FALSE)

sp.table <- df %>%
  dplyr::group_by(scientificName, sex, measurementType) %>%
  tally() %>%
  filter(n >= 10) %>%
  mutate(sex.measType = paste(sex, ".", measurementType, sep = "")) %>%
  select(-measurementType) %>%
  spread(key = sex.measType, value = n) %>%
  dplyr::group_by(scientificName) %>%
  select(-sex) %>%
  dplyr::summarise(female.mass = sum(female.mass, na.rm = TRUE),
                   male.mass = sum(male.mass, na.rm = TRUE),
                   female.total.length = sum(female.total.length, na.rm = TRUE),
                   male.total.length = sum(male.total.length, na.rm = TRUE),
                   tots = isTRUE(female.mass >= 10 & male.mass >= 10) + isTRUE(female.total.length >= 10 & male.total.length >= 10)) %>%
  filter(tots > 0) %>%
  as.data.frame()
#if a species doesn't have mass for one sex, it tends not to have total length either

write.csv(sp.table, "sp.table.csv")

sp.table <- read.csv("https://data.cyverse.org/dav-anon/iplant/home/rwalls/FuTRES_data/SSD/sp.table.csv", header = TRUE)

sp <- sp.table$scientificName

#order

sp_out <- tax_name(sp, get = 'order', db = 'ncbi')

write.csv(sp_out, "sp.taxonomy.csv")

sp_out <- read.csv("https://data.cyverse.org/dav-anon/iplant/home/rwalls/FuTRES_data/SSD/sp.taxonomy.csv", header = TRUE, stringsAsFactors = FALSE)

#2-tailed t test to see if different and get direction
#add columns "significant" (T/F) and "direction" (M/F)

rodentia <- sp_out %>%
  filter(order == "Rodentia")

rodent.sp <- rodentia$query

df.rodent <- df[df$scientificName %in% rodent.sp,]

rodent.table <- df.rodent %>%
  dplyr::group_by(scientificName, sex, measurementType) %>%
  tally() %>%
  filter(n >= 10) %>%
  mutate(sex.measType = paste(sex, ".", measurementType, sep = "")) %>%
  select(-measurementType) %>%
  spread(key = sex.measType, value = n) %>%
  dplyr::group_by(scientificName) %>%
  select(-sex) %>%
  dplyr::summarise(female.mass = sum(female.mass, na.rm = TRUE),
                   male.mass = sum(male.mass, na.rm = TRUE),
                   female.total.length = sum(female.total.length, na.rm = TRUE),
                   male.total.length = sum(male.total.length, na.rm = TRUE),
                   tots = isTRUE(female.mass >= 10 & male.mass >= 10) + isTRUE(female.total.length >= 10 & male.total.length >= 10)) %>%
  filter(tots > 0) %>%
  as.data.frame()

write.csv(df.rodent, "rodents.csv")

df.rodent <- read.csv("https://data.cyverse.org/dav-anon/iplant/home/rwalls/FuTRES_data/SSD/rodents.csv", header = TRUE, stringsAsFactors = FALSE)

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

#try on larger subset
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



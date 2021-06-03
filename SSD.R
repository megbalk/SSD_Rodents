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

itis_hierarchy(sp, what = "full")

get_ids(sp, db = "itis")

sp.id <- c()
for (i in 1:length(sp)){
  x <- get_ids(sp[i], db = "itis")
  sp.id <- append(sp.id, x$itis[[1]])
}

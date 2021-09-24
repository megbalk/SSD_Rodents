##PeMa test----
#two sample test; two sided
pema <- df.rodent[df.rodent$scientificName == "Peromyscus maniculatus",]
tt <- t.test(pema$measurementValue[pema$measurementType == "{body mass}" & pema$sex == "male"],
             pema$measurementValue[pema$measurementType == "{body mass}" & pema$sex == "female"],
             alternative = "two.sided", paired = FALSE)
#are males greater than females?
tg <- t.test(pema$measurementValue[pema$measurementType == "{body mass}" & pema$sex == "male"],
             pema$measurementValue[pema$measurementType == "{body mass}" & pema$sex == "female"],
             alternative = "greater", paired = FALSE)
#are males smaller than females?
tl <- t.test(pema$measurementValue[pema$measurementType == "{body mass}" & pema$sex == "male"],
             pema$measurementValue[pema$measurementType == "{body mass}" & pema$sex == "female"],
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
                                            df.rodent$measurementType == "{body mass}" & 
                                            df.rodent$sex == "male",])
  test.t.test$n.female[i] <- nrow(df.rodent[df.rodent$scientificName == test.sp[i] & 
                                              df.rodent$measurementType == "{body mass}" & 
                                              df.rodent$sex == "female",])
  mm <-df.rodent$measurementValue[df.rodent$scientificName == test.sp[i] & 
                                    df.rodent$measurementType == "{body mass}" & 
                                    df.rodent$sex == "male"]
  ff <- df.rodent$measurementValue[df.rodent$scientificName == test.sp[i] & 
                                     df.rodent$measurementType == "{body mass}" & 
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
  rodent.t.test.mass$n.male[i] <- nrow(df.rodent[df.rodent$scientificName == rodent.sp[i] & df.rodent$measurementType == "{body mass}" & df.rodent$sex == "male",])
  rodent.t.test.mass$n.female[i] <- nrow(df.rodent[df.rodent$scientificName == rodent.sp[i] & df.rodent$measurementType == "{body mass}" & df.rodent$sex == "female",])
  mm <- df.rodent$measurementValue[df.rodent$scientificName == rodent.sp[i] & df.rodent$measurementType == "{body mass}" & df.rodent$sex == "male"]
  ff <- df.rodent$measurementValue[df.rodent$scientificName == rodent.sp[i] & df.rodent$measurementType == "{body mass}" & df.rodent$sex == "female"]
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
  tg <- t.test(df.rodent$measurementValue[df.rodent$scientificName == rodent.sp[i] & df.rodent$measurementType == "{body mass}" & df.rodent$sex == "male"],
               df.rodent$measurementValue[df.rodent$scientificName == rodent.sp[i] & df.rodent$measurementType == "{body mass}" & df.rodent$sex == "female"],
               alternative = "greater", paired = FALSE)
  rodent.t.test.mass$tg.p.value[i] = tg$p.value
  tl <- t.test(df.rodent$measurementValue[df.rodent$scientificName == rodent.sp[i] & df.rodent$measurementType == "{body mass}" & df.rodent$sex == "male"],
               df.rodent$measurementValue[df.rodent$scientificName == rodent.sp[i] & df.rodent$measurementType == "{body mass}" & df.rodent$sex == "female"],
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
  rodent.t.test.length$n.male[i] <- nrow(df.rodent[df.rodent$scientificName == rodent.sp[i] & df.rodent$measurementType == "{body length}" & df.rodent$sex == "male",])
  rodent.t.test.length$n.female[i] <- nrow(df.rodent[df.rodent$scientificName == rodent.sp[i] & df.rodent$measurementType == "{body length}" & df.rodent$sex == "female",])
  mm <- df.rodent$measurementValue[df.rodent$scientificName == rodent.sp[i] & df.rodent$measurementType == "{body length}" & df.rodent$sex == "male"]
  ff <- df.rodent$measurementValue[df.rodent$scientificName == rodent.sp[i] & df.rodent$measurementType == "{body length}" & df.rodent$sex == "female"]
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
  tg <- t.test(df.rodent$measurementValue[df.rodent$scientificName == rodent.sp[i] & df.rodent$measurementType == "{body length}" & df.rodent$sex == "male"],
               df.rodent$measurementValue[df.rodent$scientificName == rodent.sp[i] & df.rodent$measurementType == "{body length}" & df.rodent$sex == "female"],
               alternative = "greater", paired = FALSE)
  rodent.t.test.length$tg.p.value[i] = tg$p.value
  tl <- t.test(df.rodent$measurementValue[df.rodent$scientificName == rodent.sp[i] & df.rodent$measurementType == "{body length}" & df.rodent$sex == "male"],
               df.rodent$measurementValue[df.rodent$scientificName == rodent.sp[i] & df.rodent$measurementType == "{body length}" & df.rodent$sex == "female"],
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


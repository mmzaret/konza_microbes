#Analyzing soil bacterial communities from Konza Prairie
#Study involves Bison grazing and prescribed burning
#Max Zaret, Kansas State University
#https://github.com/mmzaret/konza_microbes.git

rm(list=ls())

library(tidyverse)
library(vegan)
library(phyloseq)
library(microViz)
library(nlme)
library(lsmeans)
library(patchwork)
library(geodist)
library(picante)
library(ggpubr)
library(broom)

standard_error <- function(x) sd(x) / sqrt(length(x))

#colors for figures#
cbPalette <- c("goldenrod1", "#56B4E9", "darkorange3", "#009E73")

getwd()
setwd("/Users/mzaret/Documents/Research/KSU/Project 1/analysis") #change to your wd

#Read in data####
ps <- readRDS("cleaned_konza.ps") #files can be found on github

ps
levels(sample_data(ps)$Burn)
levels(sample_data(ps)$Graze)
levels(sample_data(ps)$Treatment)

#Look at sampling across years#
sample_data(ps) %>%
  group_by(Treatment) %>%
  dplyr::count()

#Look at sampling across years#
sample_data(ps) %>%
  group_by(Watershed) %>%
  dplyr::count(Year) 

#taxa table sorted by abundance#
View(tax_table(tax_sort(ps, by=mean)))

#pull out otu table in order to do normalizing transformation for compositional analyses
tab <- otu_table(ps)
class(tab) <- "matrix" # as.matrix() will do nothing
## you get a warning here but its nothing to worry about

#total number of sequences
sum(tab)

#Beta diversity####
#do a hellinger transformation for bray-curtis
ps.hel <- ps
otu_table(ps.hel) <- otu_table(decostand(tab, method="hellinger"), taxa_are_rows = FALSE)
ps.hel
#Calculate Bray-Curtis Distance among samples
bray_dist <- phyloseq::distance(ps.hel, method="bray")

permanova <- adonis2(bray_dist ~ sample_data(ps.hel)$Graze * sample_data(ps.hel)$Burn * sample_data(ps.hel)$Year, 
                     permutations = how(blocks=sample_data(ps.hel)$Plot, nperm=999))
permanova

#dispersion# # we don't report this in the manuscript#
#grazing effect on dispersion
mod <- betadisper(bray_dist, sample_data(ps)$Graze)
permutest(mod)

#burning effect on dispersion
mod <- betadisper(bray_dist, sample_data(ps)$Burn)
permutest(mod)

#year effect on dispersion
mod <- betadisper(bray_dist, sample_data(ps)$Year)
permutest(mod)

plot(mod)
boxplot(mod)
mod.HSD <- TukeyHSD(mod)
mod.HSD
plot(mod.HSD)

#Ordination####
#PCoA#####
ordination <- ordinate(ps.hel, method = "PCoA", distance = bray_dist)
ordination #9.1% and 7.7% variance explained by first two axes of PCOA

str(ordination$vectors)
ordination$values
scores <- as.data.frame(ordination$vectors)
scores
scores <- scores %>%
  select(Axis.1, Axis.2)

pcoa <- merge(scores, sample_data(ps.hel), by='row.names')

#calculate centroids for treatments and or years for visualization#
year.centroids <- pcoa %>%
  group_by(Year) %>%
  summarize(meanAxis1 = mean(Axis.1),
            meanAxis2 = mean(Axis.2),
            errorAxis1 = standard_error(Axis.1),
            errorAxis2 = standard_error(Axis.2))

treat.centroids <- pcoa %>%
  group_by(Treatment) %>%
  summarize(meanAxis1 = mean(Axis.1),
            meanAxis2 = mean(Axis.2),
            errorAxis1 = standard_error(Axis.1),
            errorAxis2 = standard_error(Axis.2))

year.treat.centroids <- pcoa %>%
  group_by(Year, Treatment) %>%
  summarize(meanAxis1 = mean(Axis.1),
            meanAxis2 = mean(Axis.2),
            errorAxis1 = standard_error(Axis.1),
            errorAxis2 = standard_error(Axis.2))

ordination$values$Eigenvalues
evals <- ordination$values$Eigenvalues #for normalizing plot axis rations to the relevant eigenvalue ratios

all_plot <- ggplot(pcoa, aes(Axis.1, Axis.2, color=Treatment)) +
  geom_point(alpha=0.3) +
  #stat_ellipse(type="t") +
  scale_color_manual(values=cbPalette) +
  #scale_color_brewer(palette = "Dark2") +
  coord_fixed(sqrt(evals[2] / evals[1])) +
  theme_bw()

all_plot

all_plot <- all_plot + 
  geom_point(data=treat.centroids, aes(meanAxis1, meanAxis2, color=Treatment), size=7, shape=18) +
  labs(x="Axis 1 [9.1%]", y="Axis 2 [7.7%]") +
  theme(legend.direction = "horizontal") +
  theme(legend.position= "bottom")

all_plot

#year variation#
year_plot <- ggplot(pcoa, aes(Axis.1, Axis.2, color=Year)) +
  geom_point(alpha=0.3) +
  stat_ellipse() +
  #scale_color_brewer(palette = "Dark2") +
  theme_bw() +
  coord_fixed(sqrt(evals[2] / evals[1])) 

year_plot

year_plot <- year_plot + 
  geom_point(data=year.centroids, aes(meanAxis1, meanAxis2, color=Year), size=7, shape=18) +
  labs(x="Axis 1 [9.1%]", y="Axis 2 [7.7%]")

year_plot

#treatments across years
year_treat_plot <- ggplot(pcoa, aes(Axis.1, Axis.2, color=Treatment)) +
  geom_point(alpha=0.3) +
  #stat_ellipse() +
  scale_color_manual(values=cbPalette) +
  coord_fixed(sqrt(evals[2] / evals[1])) +
  facet_wrap(~Year) +
  theme_bw() 

year_treat_plot

year_treat_plot2 <- year_treat_plot + 
  geom_point(data=year.treat.centroids, aes(meanAxis1, meanAxis2, color=Treatment), size=5, shape=18, position=position_jitter(h=0.01,w=0.01)) +
  labs(x="Axis 1 [9.1%]", y="Axis 2 [7.7%]") +
  theme(legend.position = c(0.75,0.3))

year_treat_plot2

#Effect sizes####
#Calculating variation in effect sizes of bison and fire treatments on microbial composition
#do this by measuring euclidean distance between centroids (treatments) for each year

#grazing differences
year.graze.centroids <- pcoa %>%
  group_by(Year, Graze) %>%
  summarize(meanAxis1 = mean(Axis.1),
            meanAxis2 = mean(Axis.2),
            seAxis1 = standard_error(Axis.1),
            seAxis2 = standard_error(Axis.2))

year.graze.centroids

graze2016 <- year.graze.centroids %>%
  ungroup() %>%
  filter(Year == "2016") %>%
  column_to_rownames(., var="Graze") %>%
  select(meanAxis1, meanAxis2) %>%
  dist(., method="euclidean", diag=FALSE, upper=FALSE) %>%
  as.matrix(.) %>% 
  as.data.frame(.) 

graze2018 <- year.graze.centroids %>%
  ungroup() %>%
  filter(Year == "2018") %>%
  column_to_rownames(., var="Graze") %>%
  select(meanAxis1, meanAxis2) %>%
  dist(., method="euclidean", diag=FALSE, upper=FALSE) %>%
  as.matrix(.) %>% 
  as.data.frame(.) 

graze2020 <- year.graze.centroids %>%
  ungroup() %>%
  filter(Year == "2020") %>%
  column_to_rownames(., var="Graze") %>%
  select(meanAxis1, meanAxis2) %>%
  dist(., method="euclidean", diag=FALSE, upper=FALSE) %>%
  as.matrix(.) %>% 
  as.data.frame(.) 

graze2022 <- year.graze.centroids %>%
  ungroup() %>%
  filter(Year == "2022") %>%
  column_to_rownames(., var="Graze") %>%
  select(meanAxis1, meanAxis2) %>%
  dist(., method="euclidean", diag=FALSE, upper=FALSE) %>%
  as.matrix(.) %>% 
  as.data.frame(.) 

graze2023 <- year.graze.centroids %>%
  ungroup() %>%
  filter(Year == "2023") %>%
  column_to_rownames(., var="Graze") %>%
  select(meanAxis1, meanAxis2) %>%
  dist(., method="euclidean", diag=FALSE, upper=FALSE) %>%
  as.matrix(.) %>% 
  as.data.frame(.) 

graze_effects <- c(graze2016[1,2], graze2018[1,2], graze2020[1,2], graze2022[1,2], graze2023[1,2]) #select the corresponding euclidean distance value from each year's distance matrix
year <- c(2016,2018,2020,2022,2023)

graze_effects_df <- data.frame(year,graze_effects)

graze_effects_df <- graze_effects_df %>%
  mutate(year = as.factor(year))

grazeplot <- ggplot(graze_effects_df, aes(x=year, y=graze_effects)) +
  geom_point(size=3) +
  geom_hline(yintercept=0, lty=2) +
  #scale_y_continuous(limits=c(0,0.3)) +
  labs(y="Grazing effect", x="Year") + 
  theme_bw()

grazeplot

#burning differences
year.burn.centroids <- pcoa %>%
  group_by(Year, Burn) %>%
  summarize(meanAxis1 = mean(Axis.1),
            meanAxis2 = mean(Axis.2))

burn2016 <- year.burn.centroids %>%
  ungroup() %>%
  filter(Year == "2016") %>%
  column_to_rownames(., var="Burn") %>%
  select(meanAxis1, meanAxis2) %>%
  dist(., method="euclidean", diag=FALSE, upper=FALSE) %>%
  as.matrix(.) %>% 
  as.data.frame(.) 

burn2018 <- year.burn.centroids %>%
  ungroup() %>%
  filter(Year == "2018") %>%
  column_to_rownames(., var="Burn") %>%
  select(meanAxis1, meanAxis2) %>%
  dist(., method="euclidean", diag=FALSE, upper=FALSE) %>%
  as.matrix(.) %>% 
  as.data.frame(.) 

burn2020 <- year.burn.centroids %>%
  ungroup() %>%
  filter(Year == "2020") %>%
  column_to_rownames(., var="Burn") %>%
  select(meanAxis1, meanAxis2) %>%
  dist(., method="euclidean", diag=FALSE, upper=FALSE) %>%
  as.matrix(.) %>% 
  as.data.frame(.) 

burn2022 <- year.burn.centroids %>%
  ungroup() %>%
  filter(Year == "2022") %>%
  column_to_rownames(., var="Burn") %>%
  select(meanAxis1, meanAxis2) %>%
  dist(., method="euclidean", diag=FALSE, upper=FALSE) %>%
  as.matrix(.) %>% 
  as.data.frame(.) 

burn2023 <- year.burn.centroids %>%
  ungroup() %>%
  filter(Year == "2023") %>%
  column_to_rownames(., var="Burn") %>%
  select(meanAxis1, meanAxis2) %>%
  dist(., method="euclidean", diag=FALSE, upper=FALSE) %>%
  as.matrix(.) %>% 
  as.data.frame(.) 

burn_effects <- c(burn2016[1,2], burn2018[1,2], burn2020[1,2], burn2022[1,2], burn2023[1,2])
year <- c(2016,2018,2020,2022,2023)

burn_effects_df <- data.frame(year,burn_effects)

burn_effects_df <- burn_effects_df %>%
  mutate(year = as.factor(year))

burnplot <- ggplot(burn_effects_df, aes(x=year, y=burn_effects)) +
  geom_point(size=3) +
  geom_hline(yintercept=0, lty=2) +
  #scale_y_continuous(limits=c(0,0.25)) +
  labs(y="Burning effect", x="Year") +
  theme_bw() 
burnplot

#combine the bison and fire effect size plots
effectsplot <- grazeplot + burnplot
effectsplot

trt_effects_df <- merge(graze_effects_df, burn_effects_df, by="year")

#bring in climate for comparison #need climate code for this
climate <- read.csv("konza_climatesummary.csv")


trt_effects_df <- merge(trt_effects_df, climate, by="year")
trt_effects_df
str(trt_effects_df)

#compare to precip and growing season precip
cor.test(trt_effects_df$graze_effects, trt_effects_df$total_precip)
cor.test(trt_effects_df$graze_effects, trt_effects_df$gs_precip)
cor.test(trt_effects_df$burn_effects, trt_effects_df$total_precip)
cor.test(trt_effects_df$burn_effects, trt_effects_df$gs_precip)
#only report total precip in manuscript given no difference with growing season precipitation

library(ggrepel)
graze_precip_plot <- ggplot(trt_effects_df, aes(total_precip, graze_effects)) +
  geom_point() +
  geom_text_repel(label=year) +
  labs(x="Annual Precipitation (mm)", y="Grazing effect") +
  coord_cartesian(xlim=c(500,1100), ylim=c(0.00,0.1)) +
  theme_bw()

burn_precip_plot <- ggplot(trt_effects_df, aes(total_precip, burn_effects)) +
  geom_point() +
  geom_text_repel(label=year) +
  labs(x="Annual Precipitation (mm)", y="Burning effect") +
  coord_cartesian(xlim=c(500,1100), ylim=c(0.01,0.18)) +
  theme_bw()

#visualize effects and precip

grazeplot + burnplot + graze_precip_plot + burn_precip_plot

#Distance dissimilarity####
#For each treatment in a given year - 5 years x 4 treatments = 20 df where we will calculate community distances
#then combine dfs to run an ancova of how community distance varies by geographic distance, bison, fire, and year of sampling

ps.ra <- ps.hel #just changing name of phyloseq object because I am too lazy to rename ps.ra to ps.hel over and over after originally writing this code

#UG20
ps2016ug20 <- ps.ra %>%
  ps_filter(Treatment == "UG20" & Year == "2016") %>%
  ps_select(LongDD, LatDD)

bray <- phyloseq::distance(ps2016ug20, method="bray")

latlong <- as.matrix(sample_data(ps2016ug20))

dist_matrix <- geodist(latlong, measure="cheap")

dist_matrix[!lower.tri(dist_matrix)] <- NA
dist_matrix[is.na(dist_matrix)] <- " "

bray_matrix <- as.matrix(bray)

bray_matrix[!lower.tri(bray_matrix)] <- NA
bray_matrix[is.na(bray_matrix)] <- " "

mantel(bray_matrix, dist_matrix)

df_2016ug20 <- data.frame(
  bray = as.vector(bray_matrix),
  dist = as.vector(dist_matrix)) %>%
  mutate(Year = 2016) %>%
  mutate(Treatment = "UG20") %>%
  filter(bray != " ")

ps2018ug20 <- ps.ra %>%
  ps_filter(Treatment == "UG20" & Year == "2018") %>%
  ps_select(LongDD, LatDD)

bray <- phyloseq::distance(ps2018ug20, method="bray")

latlong <- as.matrix(sample_data(ps2018ug20))

dist_matrix <- geodist(latlong, measure="cheap")

dist_matrix[!lower.tri(dist_matrix)] <- NA
dist_matrix[is.na(dist_matrix)] <- " "

bray_matrix <- as.matrix(bray)

bray_matrix[!lower.tri(bray_matrix)] <- NA
bray_matrix[is.na(bray_matrix)] <- " "

mantel(bray_matrix, dist_matrix)

df_2018ug20 <- data.frame(
  bray = as.vector(bray_matrix),
  dist = as.vector(dist_matrix)) %>%
  mutate(Year = 2018) %>%
  mutate(Treatment = "UG20") %>%
  filter(bray != " ")

ps2020ug20 <- ps.ra %>%
  ps_filter(Treatment == "UG20" & Year == "2020") %>%
  ps_select(LongDD, LatDD)

bray <- phyloseq::distance(ps2020ug20, method="bray")

latlong <- as.matrix(sample_data(ps2020ug20))

dist_matrix <- geodist(latlong, measure="cheap")

dist_matrix[!lower.tri(dist_matrix)] <- NA
dist_matrix[is.na(dist_matrix)] <- " "

bray_matrix <- as.matrix(bray)

bray_matrix[!lower.tri(bray_matrix)] <- NA
bray_matrix[is.na(bray_matrix)] <- " "

mantel(bray_matrix, dist_matrix)

df_2020ug20 <- data.frame(
  bray = as.vector(bray_matrix),
  dist = as.vector(dist_matrix)) %>%
  mutate(Year = 2020) %>%
  mutate(Treatment = "UG20") %>%
  filter(bray != " ")

ps2022ug20 <- ps.ra %>%
  ps_filter(Treatment == "UG20" & Year == "2022") %>%
  ps_select(LongDD, LatDD)

bray <- phyloseq::distance(ps2022ug20, method="bray")

latlong <- as.matrix(sample_data(ps2022ug20))

dist_matrix <- geodist(latlong, measure="cheap")

dist_matrix[!lower.tri(dist_matrix)] <- NA
dist_matrix[is.na(dist_matrix)] <- " "

bray_matrix <- as.matrix(bray)

bray_matrix[!lower.tri(bray_matrix)] <- NA
bray_matrix[is.na(bray_matrix)] <- " "

mantel(bray_matrix, dist_matrix)

df_2022ug20 <- data.frame(
  bray = as.vector(bray_matrix),
  dist = as.vector(dist_matrix)) %>%
  mutate(Year = 2022) %>%
  mutate(Treatment = "UG20") %>%
  filter(bray != " ")

ps2023ug20 <- ps.ra %>%
  ps_filter(Treatment == "UG20" & Year == "2023") %>%
  ps_select(LongDD, LatDD)

bray <- phyloseq::distance(ps2023ug20, method="bray")

latlong <- as.matrix(sample_data(ps2023ug20))

dist_matrix <- geodist(latlong, measure="cheap")

dist_matrix[!lower.tri(dist_matrix)] <- NA
dist_matrix[is.na(dist_matrix)] <- " "

bray_matrix <- as.matrix(bray)

bray_matrix[!lower.tri(bray_matrix)] <- NA
bray_matrix[is.na(bray_matrix)] <- " "

mantel(bray_matrix, dist_matrix)

df_2023ug20 <- data.frame(
  bray = as.vector(bray_matrix),
  dist = as.vector(dist_matrix)) %>%
  mutate(Year = 2023) %>%
  mutate(Treatment = "UG20") %>%
  filter(bray != " ")

##UG1##
ps2016ug1 <- ps.ra %>%
  ps_filter(Treatment == "UG1" & Year == "2016") %>%
  ps_select(LongDD, LatDD)

bray <- phyloseq::distance(ps2016ug1, method="bray")

latlong <- as.matrix(sample_data(ps2016ug1))

dist_matrix <- geodist(latlong, measure="cheap")

dist_matrix[!lower.tri(dist_matrix)] <- NA
dist_matrix[is.na(dist_matrix)] <- " "

bray_matrix <- as.matrix(bray)

bray_matrix[!lower.tri(bray_matrix)] <- NA
bray_matrix[is.na(bray_matrix)] <- " "

mantel(bray_matrix, dist_matrix)

df_2016ug1 <- data.frame(
  bray = as.vector(bray_matrix),
  dist = as.vector(dist_matrix)) %>%
  mutate(Year = 2016) %>%
  mutate(Treatment = "UG1") %>%
  filter(bray != " ")

ps2018ug1 <- ps.ra %>%
  ps_filter(Treatment == "UG1" & Year == "2018") %>%
  ps_select(LongDD, LatDD)

bray <- phyloseq::distance(ps2018ug1, method="bray")

latlong <- as.matrix(sample_data(ps2018ug1))

dist_matrix <- geodist(latlong, measure="cheap")

dist_matrix[!lower.tri(dist_matrix)] <- NA
dist_matrix[is.na(dist_matrix)] <- " "

bray_matrix <- as.matrix(bray)

bray_matrix[!lower.tri(bray_matrix)] <- NA
bray_matrix[is.na(bray_matrix)] <- " "

mantel(bray_matrix, dist_matrix)

df_2018ug1 <- data.frame(
  bray = as.vector(bray_matrix),
  dist = as.vector(dist_matrix)) %>%
  mutate(Year = 2018) %>%
  mutate(Treatment = "UG1") %>%
  filter(bray != " ")

ps2020ug1 <- ps.ra %>%
  ps_filter(Treatment == "UG1" & Year == "2020") %>%
  ps_select(LongDD, LatDD)

bray <- phyloseq::distance(ps2020ug1, method="bray")

latlong <- as.matrix(sample_data(ps2020ug1))

dist_matrix <- geodist(latlong, measure="cheap")

dist_matrix[!lower.tri(dist_matrix)] <- NA
dist_matrix[is.na(dist_matrix)] <- " "

bray_matrix <- as.matrix(bray)

bray_matrix[!lower.tri(bray_matrix)] <- NA
bray_matrix[is.na(bray_matrix)] <- " "

mantel(bray_matrix, dist_matrix)

df_2020ug1 <- data.frame(
  bray = as.vector(bray_matrix),
  dist = as.vector(dist_matrix)) %>%
  mutate(Year = 2020) %>%
  mutate(Treatment = "UG1") %>%
  filter(bray != " ")

ps2022ug1 <- ps.ra %>%
  ps_filter(Treatment == "UG1" & Year == "2022") %>%
  ps_select(LongDD, LatDD)

bray <- phyloseq::distance(ps2022ug1, method="bray")

latlong <- as.matrix(sample_data(ps2022ug1))

dist_matrix <- geodist(latlong, measure="cheap")

dist_matrix[!lower.tri(dist_matrix)] <- NA
dist_matrix[is.na(dist_matrix)] <- " "

bray_matrix <- as.matrix(bray)

bray_matrix[!lower.tri(bray_matrix)] <- NA
bray_matrix[is.na(bray_matrix)] <- " "

mantel(bray_matrix, dist_matrix)

df_2022ug1 <- data.frame(
  bray = as.vector(bray_matrix),
  dist = as.vector(dist_matrix)) %>%
  mutate(Year = 2022) %>%
  mutate(Treatment = "UG1") %>%
  filter(bray != " ")

ps2023ug1 <- ps.ra %>%
  ps_filter(Treatment == "UG1" & Year == "2023") %>%
  ps_select(LongDD, LatDD)

bray <- phyloseq::distance(ps2023ug1, method="bray")

latlong <- as.matrix(sample_data(ps2023ug1))

dist_matrix <- geodist(latlong, measure="cheap")

dist_matrix[!lower.tri(dist_matrix)] <- NA
dist_matrix[is.na(dist_matrix)] <- " "

bray_matrix <- as.matrix(bray)


bray_matrix[!lower.tri(bray_matrix)] <- NA
bray_matrix[is.na(bray_matrix)] <- " "


mantel(bray_matrix, dist_matrix)

df_2023ug1 <- data.frame(
  bray = as.vector(bray_matrix),
  dist = as.vector(dist_matrix)) %>%
  mutate(Year = 2023) %>%
  mutate(Treatment = "UG1") %>%
  filter(bray != " ")

##G20
ps2016g20 <- ps.ra %>%
  ps_filter(Treatment == "G20" & Year == "2016") %>%
  ps_select(LongDD, LatDD)

bray <- phyloseq::distance(ps2016g20, method="bray")

latlong <- as.matrix(sample_data(ps2016g20))

dist_matrix <- geodist(latlong, measure="cheap")

dist_matrix[!lower.tri(dist_matrix)] <- NA
dist_matrix[is.na(dist_matrix)] <- " "

bray_matrix <- as.matrix(bray)

bray_matrix[!lower.tri(bray_matrix)] <- NA
bray_matrix[is.na(bray_matrix)] <- " "

mantel(bray_matrix, dist_matrix)

df_2016g20 <- data.frame(
  bray = as.vector(bray_matrix),
  dist = as.vector(dist_matrix)) %>%
  mutate(Year = 2016) %>%
  mutate(Treatment = "G20") %>%
  filter(bray != " ")

ps2018g20 <- ps.ra %>%
  ps_filter(Treatment == "G20" & Year == "2018") %>%
  ps_select(LongDD, LatDD)

bray <- phyloseq::distance(ps2018g20, method="bray")

latlong <- as.matrix(sample_data(ps2018g20))

dist_matrix <- geodist(latlong, measure="cheap")

dist_matrix[!lower.tri(dist_matrix)] <- NA
dist_matrix[is.na(dist_matrix)] <- " "

bray_matrix <- as.matrix(bray)

bray_matrix[!lower.tri(bray_matrix)] <- NA
bray_matrix[is.na(bray_matrix)] <- " "

mantel(bray_matrix, dist_matrix)

df_2018g20 <- data.frame(
  bray = as.vector(bray_matrix),
  dist = as.vector(dist_matrix)) %>%
  mutate(Year = 2018) %>%
  mutate(Treatment = "G20") %>%
  filter(bray != " ")

ps2020g20 <- ps.ra %>%
  ps_filter(Treatment == "G20" & Year == "2020") %>%
  ps_select(LongDD, LatDD)

bray <- phyloseq::distance(ps2020g20, method="bray")

latlong <- as.matrix(sample_data(ps2020g20))

dist_matrix <- geodist(latlong, measure="cheap")

dist_matrix[!lower.tri(dist_matrix)] <- NA
dist_matrix[is.na(dist_matrix)] <- " "

bray_matrix <- as.matrix(bray)

bray_matrix[!lower.tri(bray_matrix)] <- NA
bray_matrix[is.na(bray_matrix)] <- " "

mantel(bray_matrix, dist_matrix)

df_2020g20 <- data.frame(
  bray = as.vector(bray_matrix),
  dist = as.vector(dist_matrix)) %>%
  mutate(Year = 2020) %>%
  mutate(Treatment = "G20") %>%
  filter(bray != " ")

ps2022g20 <- ps.ra %>%
  ps_filter(Treatment == "G20" & Year == "2022") %>%
  ps_select(LongDD, LatDD)

bray <- phyloseq::distance(ps2022g20, method="bray")

latlong <- as.matrix(sample_data(ps2022g20))

dist_matrix <- geodist(latlong, measure="cheap")

dist_matrix[!lower.tri(dist_matrix)] <- NA
dist_matrix[is.na(dist_matrix)] <- " "

bray_matrix <- as.matrix(bray)

bray_matrix[!lower.tri(bray_matrix)] <- NA
bray_matrix[is.na(bray_matrix)] <- " "

mantel(bray_matrix, dist_matrix)

df_2022g20 <- data.frame(
  bray = as.vector(bray_matrix),
  dist = as.vector(dist_matrix)) %>%
  mutate(Year = 2022) %>%
  mutate(Treatment = "G20") %>%
  filter(bray != " ")

ps2023g20 <- ps.ra %>%
  ps_filter(Treatment == "G20" & Year == "2023") %>%
  ps_select(LongDD, LatDD)

bray <- phyloseq::distance(ps2023g20, method="bray")

latlong <- as.matrix(sample_data(ps2023g20))

dist_matrix <- geodist(latlong, measure="cheap")

dist_matrix[!lower.tri(dist_matrix)] <- NA
dist_matrix[is.na(dist_matrix)] <- " "

bray_matrix <- as.matrix(bray)

bray_matrix[!lower.tri(bray_matrix)] <- NA
bray_matrix[is.na(bray_matrix)] <- " "

mantel(bray_matrix, dist_matrix)

df_2023g20 <- data.frame(
  bray = as.vector(bray_matrix),
  dist = as.vector(dist_matrix)) %>%
  mutate(Year = 2023) %>%
  mutate(Treatment = "G20") %>%
  filter(bray != " ")

##G1##
ps2016g1 <- ps.ra %>%
  ps_filter(Treatment == "G1" & Year == "2016") %>%
  ps_select(LongDD, LatDD)

bray <- phyloseq::distance(ps2016g1, method="bray")

latlong <- as.matrix(sample_data(ps2016g1))

dist_matrix <- geodist(latlong, measure="cheap")

dist_matrix[!lower.tri(dist_matrix)] <- NA
dist_matrix[is.na(dist_matrix)] <- " "

bray_matrix <- as.matrix(bray)

bray_matrix[!lower.tri(bray_matrix)] <- NA
bray_matrix[is.na(bray_matrix)] <- " "

mantel(bray_matrix, dist_matrix)

df_2016g1 <- data.frame(
  bray = as.vector(bray_matrix),
  dist = as.vector(dist_matrix)) %>%
  mutate(Year = 2016) %>%
  mutate(Treatment = "G1") %>%
  filter(bray != " ")

ps2018g1 <- ps.ra %>%
  ps_filter(Treatment == "G1" & Year == "2018") %>%
  ps_select(LongDD, LatDD)

bray <- phyloseq::distance(ps2018g1, method="bray")

latlong <- as.matrix(sample_data(ps2018g1))

dist_matrix <- geodist(latlong, measure="cheap")

dist_matrix[!lower.tri(dist_matrix)] <- NA
dist_matrix[is.na(dist_matrix)] <- " "

bray_matrix <- as.matrix(bray)

bray_matrix[!lower.tri(bray_matrix)] <- NA
bray_matrix[is.na(bray_matrix)] <- " "

mantel(bray_matrix, dist_matrix)

df_2018g1 <- data.frame(
  bray = as.vector(bray_matrix),
  dist = as.vector(dist_matrix)) %>%
  mutate(Year = 2018) %>%
  mutate(Treatment = "G1") %>%
  filter(bray != " ")

ps2020g1 <- ps.ra %>%
  ps_filter(Treatment == "G1" & Year == "2020") %>%
  ps_select(LongDD, LatDD)

bray <- phyloseq::distance(ps2020g1, method="bray")

latlong <- as.matrix(sample_data(ps2020g1))

dist_matrix <- geodist(latlong, measure="cheap")

dist_matrix[!lower.tri(dist_matrix)] <- NA
dist_matrix[is.na(dist_matrix)] <- " "

bray_matrix <- as.matrix(bray)

bray_matrix[!lower.tri(bray_matrix)] <- NA
bray_matrix[is.na(bray_matrix)] <- " "

mantel(bray_matrix, dist_matrix)

df_2020g1 <- data.frame(
  bray = as.vector(bray_matrix),
  dist = as.vector(dist_matrix)) %>%
  mutate(Year = 2020) %>%
  mutate(Treatment = "G1") %>%
  filter(bray != " ")

ps2022g1 <- ps.ra %>%
  ps_filter(Treatment == "G1" & Year == "2022") %>%
  ps_select(LongDD, LatDD)

bray <- phyloseq::distance(ps2022g1, method="bray")

latlong <- as.matrix(sample_data(ps2022g1))

dist_matrix <- geodist(latlong, measure="cheap")

dist_matrix[!lower.tri(dist_matrix)] <- NA
dist_matrix[is.na(dist_matrix)] <- " "

bray_matrix <- as.matrix(bray)

bray_matrix[!lower.tri(bray_matrix)] <- NA
bray_matrix[is.na(bray_matrix)] <- " "

mantel(bray_matrix, dist_matrix)

df_2022g1 <- data.frame(
  bray = as.vector(bray_matrix),
  dist = as.vector(dist_matrix)) %>%
  mutate(Year = 2022) %>%
  mutate(Treatment = "G1") %>%
  filter(bray != " ")

ps2023g1 <- ps.ra %>%
  ps_filter(Treatment == "G1" & Year == "2023") %>%
  ps_select(LongDD, LatDD)

bray <- phyloseq::distance(ps2023g1, method="bray")

latlong <- as.matrix(sample_data(ps2023g1))

dist_matrix <- geodist(latlong, measure="cheap")

dist_matrix[!lower.tri(dist_matrix)] <- NA
dist_matrix[is.na(dist_matrix)] <- " "

bray_matrix <- as.matrix(bray)

bray_matrix[!lower.tri(bray_matrix)] <- NA
bray_matrix[is.na(bray_matrix)] <- " "

mantel(bray_matrix, dist_matrix)

df_2023g1 <- data.frame(
  bray = as.vector(bray_matrix),
  dist = as.vector(dist_matrix)) %>%
  mutate(Year = 2023) %>%
  mutate(Treatment = "G1") %>%
  filter(bray != " ")

#combine dfs#
dd_df <- df_2016ug20 %>%
  rbind(df_2016g20) %>%
  rbind(df_2016ug1) %>%
  rbind(df_2016g1) %>%
  rbind(df_2018ug20) %>%
  rbind(df_2018ug1) %>%
  rbind(df_2018g20) %>%
  rbind(df_2018g1) %>%
  rbind(df_2020ug20) %>%
  rbind(df_2020ug1) %>%
  rbind(df_2020g20) %>%
  rbind(df_2020g1) %>%
  rbind(df_2022ug20) %>%
  rbind(df_2022ug1) %>%
  rbind(df_2022g20) %>%
  rbind(df_2022g1) %>%
  rbind(df_2023ug20) %>%
  rbind(df_2023g20) %>%
  rbind(df_2023ug1) %>%
  rbind(df_2023g1)

#have to rename and relevel factors
dd_df <- dd_df %>%
  mutate(Treatment = as.factor(Treatment),
         Graze = case_when(Treatment == "UG20" ~ "No Grazers",
                           Treatment == "UG1" ~ "No Grazers",
                           Treatment == "G20" ~ "Bison Present",
                           Treatment == "G1" ~ "Bison Present"),
         Burn = case_when(Treatment == "UG20" | Treatment == "G20" ~ "Infrequent Burn",
                          Treatment == "UG1" | Treatment == "G1" ~ "Frequent Burn"))

dd_df <- dd_df %>%
  mutate(Year = as.factor(as.character(Year)),
         Graze = as.factor(Graze),
         Burn = as.factor(Burn),
         bray = as.numeric(bray),
         dist = as.numeric(dist))

levels(dd_df$Year)

levels(dd_df$Treatment)
dd_df$Treatment <- relevel(dd_df$Treatment, ref='G1')
dd_df$Treatment <- relevel(dd_df$Treatment, ref='UG1')
dd_df$Treatment <- relevel(dd_df$Treatment, ref='G20')
dd_df$Treatment <- relevel(dd_df$Treatment, ref='UG20')
levels(dd_df$Treatment)

levels(dd_df$Graze)
dd_df$Graze <- relevel(dd_df$Graze, ref='Bison Present')
dd_df$Graze <- relevel(dd_df$Graze, ref='No Grazers')
levels(dd_df$Graze)

levels(dd_df$Burn)
dd_df$Burn <- relevel(dd_df$Burn, ref='Frequent Burn')
dd_df$Burn <- relevel(dd_df$Burn, ref='Infrequent Burn')
levels(dd_df$Burn)

#now run ancova and visualize#
#we log transform geographic distance to match sampling scheme which was done on a log scale
model <- (lm((bray) ~ log(dist) * Graze * Burn * Year, data=dd_df))
summary(model)
ancova <- anova(model)
ancova

#do distance dissimilarity slopes differ across treatments?
#Run t tests that determine if the slope between community distance and geographic distance varies between experimental treatments
#Bison
m.lst <- lstrends(model, ~Graze, var="log(dist)")
m.lst
pairs(m.lst)
#Higher slope with bison present

#Burn
m.lst <- lstrends(model, ~Burn, var="log(dist)")
m.lst
pairs(m.lst)


#interaction
m.lst <- lstrends(model, ~Graze & Burn, var="log(dist)")
m.lst
pairs(m.lst)
#this is meaningful comparison that is accounting for both
#use these estimates and p values for the subset inside the the DDR figure

#graze and year interaction
m.lst <- lstrends(model, ~Graze & Year, var="log(dist)")
m.lst
pairs(m.lst)

#burn and year interaction
m.lst <- lstrends(model, ~Burn & Year, var="log(dist)")
m.lst
pairs(m.lst)

#graze and burn and year interaction
m.lst <- lstrends(model, ~Graze & Burn & Year, var="log(dist)")
m.lst
pairs(m.lst)
#this last one spits out every slope from every treatment and year if someone wants that...

#visuals#
ggplot(dd_df, aes(x=log(dist), y=(bray))) +
  geom_point(alpha=0.1) +
  stat_smooth(method=lm, se=TRUE, fullrange=FALSE) +
  scale_color_manual(values=cbPalette) +
  stat_cor() +
  #scale_x_continuous(limits=c(-2.5,8)) +
  labs(x="Distance (log meters)", y="Bray-Curtis Dissimilarity")

#All Treatments visualized#
dd_df <- dd_df %>%
  mutate(dist_labels = as.factor(dist),
         logdist = log(dist))

ddplot <- ggplot(dd_df, aes(x=logdist, y=(bray), color=Treatment)) +
  geom_point(alpha=0.1) +
  stat_smooth(method=lm, se=TRUE, fullrange=TRUE) +
  scale_color_manual(values=cbPalette) +
  scale_y_continuous(limits=c(0.3,0.6)) +
  scale_x_continuous(breaks =c(-2.1970281, 0.0001964805, 2.303711, 4.612754, 7.058103), labels=c("0.1", "1", "10", "100", "1000")) +
  labs(x="Geographic Distance (Meters)", y="Bray-Curtis Dissimilarity") +
  theme_bw()
ddplot

#treatment by year
ggplot(dd_df, aes(x=log(dist), y=(bray), color=Treatment)) +
  geom_point(alpha=0.1) +
  stat_smooth(method=lm, se=FALSE, fullrange=TRUE) +
  scale_color_manual(values=cbPalette) +
  #stat_cor() +
  labs(x="Distance (log meters)", y="Bray-Curtis Dissimilarity") +
  facet_wrap(~Year) +
  theme_bw() +
  theme(legend.position = "bottom")
#year
ggplot(dd_df, aes(x=log(dist), y=(bray), color=Year)) +
  geom_point(alpha=0.1) +
  stat_smooth(method=lm, se=TRUE, fullrange=TRUE) +
  #scale_color_brewer(palette="Dark2") +
  stat_cor() +
  labs(x="Distance (log meters)", y="Bray-Curtis Dissimilarity")

#graze
ggplot(dd_df, aes(x=log(dist), y=log(bray), color=Graze)) +
  geom_point(alpha=0.1) +
  stat_smooth(method=lm, se=FALSE, fullrange=TRUE) +
  scale_color_brewer(palette="Dark2") +
  #scale_x_continuous(limits=c(-2.5,8)) +
  stat_cor() +
  labs(x="Distance (meters)", y="Bray-Curtis Dissimilarity")

#burn
ggplot(dd_df, aes(x=log(dist), y=log(bray), color=Burn)) +
  geom_point(alpha=0.1) +
  stat_smooth(method=lm, se=FALSE, fullrange=TRUE) +
  scale_color_brewer(palette="Dark2") +
  #scale_x_continuous(limits=c(-2.5,8)) +
  stat_cor() +
  labs(x="Distance (meters)", y="Bray-Curtis Dissimilarity")

#make a figure that shows slope estimates with CL and add significance groups
#interaction
m.lst <- lstrends(model, ~Graze & Burn, var="log(dist)")
m.lst
pairs(m.lst)
slopes <- as.data.frame(m.lst)
str(slopes)
slopes <- slopes %>%
  mutate(Treatment = case_when(Graze == "No Grazers" & Burn == "Infrequent Burn" ~ "UG20",
                               Graze == "No Grazers" & Burn == "Frequent Burn" ~ "UG1",
                               Graze == "Bison Present" & Burn == "Infrequent Burn" ~ "G20",
                               Graze == "Bison Present" & Burn == "Frequent Burn" ~ "G1",)) %>%
  mutate(Treatment = as.factor(Treatment),
         slope = `log(dist).trend`)
levels(slopes$Treatment)
slopes$Treatment <- relevel(slopes$Treatment, ref='G1')
slopes$Treatment <- relevel(slopes$Treatment, ref='UG1')
slopes$Treatment <- relevel(slopes$Treatment, ref='G20')
slopes$Treatment <- relevel(slopes$Treatment, ref='UG20')

ddslopes <- ggplot(slopes, aes(Treatment, slope)) +
  geom_point() +
  geom_errorbar(aes(ymin=slope-SE, ymax=slope+SE),width=0.2) +
  annotate("text", x=1, y=0.008, label = "A", size=5) +
  annotate("text", x=3, y=0.008, label = "A", size=5) +
  annotate("text", x=2, y=0.008, label = "B", size=5) +
  annotate("text", x=4, y=0.008, label = "A", size=5) +
  labs(y="Slope", x="") +
  coord_cartesian(ylim=c(0.004,0.011)) +
  theme_bw()
ddslopes

ddplot + inset_element(ddslopes, 0.01,0.6,0.4,0.99, align_to = "panel")

ggsave(filename="ddplot.pdf", 
       path="/Users/mzaret/Documents/Research/KSU/Project 1/analysis/figures",
       width=8, height = 6)

#Alpha diversity####
#read in rarefied asv table (rarefaction with 1000 iterations)
ps.rare <- readRDS("konza_asv_rare7222024.ps")
otu_table(ps.rare) <- round(otu_table(ps.rare)) #need to round to integers for function to work

#phyloseq estimate_richness approach to alpha diversity
alphadiv <- estimate_richness(ps.rare)
alphadiv <- merge(alphadiv, sample_data(ps.rare), by="row.names")

#calculate phylogenetic diversity (Faith's D)
#need to add tree to rarified dataset#
taxa_names(ps.rare) <- paste("ASV_", 1:ntaxa(ps.rare), sep="")
phy_tree(ps.rare) <- phy_tree(ps)

phylodiv <- pd(otu_table(ps.rare), phy_tree(ps.rare), include.root=FALSE)

phylodiv <- phylodiv %>%
  dplyr::select(PD) %>%
  mutate(SampleID = row.names(phylodiv))

#merge phylogenetic diversity with other alpha diversity measures
alphadiv <- merge(alphadiv, phylodiv, by="SampleID")

cor.test(alphadiv$Observed, alphadiv$PD)

#remove this outlier sample
#alphadiv <- alphadiv %>%
#  filter(Row.names != "Samp131_2023")

levels(alphadiv$Treatment)
alphadiv$Treatment <- relevel(alphadiv$Treatment, ref='G1')
alphadiv$Treatment <- relevel(alphadiv$Treatment, ref='UG1')
alphadiv$Treatment <- relevel(alphadiv$Treatment, ref='G20')
alphadiv$Treatment <- relevel(alphadiv$Treatment, ref='UG20')
levels(alphadiv$Treatment)
levels(alphadiv$Graze)
levels(alphadiv$Burn)
#now stats and visualize
#Richness
hist((alphadiv$Observed))

lme1 <- lme(fixed = (Observed) ~ Graze * Burn * Year,
            random = ~1 | Plot,
            data=alphadiv)

plot(lme1)
summary(lme1)
anova(lme1)

lsmeans(lme(fixed = (Observed) ~ Graze * Burn * Year,
            random = ~1 | Plot,
            data=alphadiv),
        pairwise ~ Graze * Burn|Year, type="Tukey")

#in 2023:
#12% reduction in observed asv richness

#Phylogenetic diversity
hist((alphadiv$PD))

lme1 <- lme(fixed = (PD) ~ Graze * Burn * Year,
            random = ~1 | Plot,
            data=alphadiv)

plot(lme1)
summary(lme1)
anova(lme1)

lsmeans(lme(fixed = PD ~ Graze * Burn * Year,
            random = ~1 | Plot,
            data=alphadiv),
        pairwise ~ Graze * Burn|Year, type="Tukey")

#visualize alpha diversity
richnessplot <- alphadiv %>%
  ggplot(., aes(Treatment, Observed, color=Treatment)) +
  scale_color_manual(values=cbPalette) +
  geom_jitter(alpha=0.2) +
  stat_summary(fun.data="mean_cl_boot", size=1) +
  theme_bw() +
  labs(y="Microbial Richness (ASV)", x="") +
  #theme(axis.ticks.x=element_blank(),
  #      axis.text.x = element_blank(),
  #      axis.title.x = element_blank()) +
  coord_cartesian(ylim=c(2910,4750)) +
  theme(legend.position = "none")

richnessplot

ggsave(filename="richnessplot.pdf", 
       path="/Users/mzaret/Documents/Research/KSU/Project 1/analysis/figures",
       width=6, height = 6)

alphadiv %>%
  ggplot(., aes(Treatment, Observed, color=Treatment)) +
  scale_color_manual(values=cbPalette) +
  geom_jitter(alpha=0.2) +
  stat_summary(fun.data="mean_cl_boot", size=1) +
  theme_bw() +
  labs(y="Microbial Richness (ASV)", x="") +
  theme(axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank()) +
  theme(legend.position = "bottom") +
  facet_wrap(~Year, scales="free_y")

#figure of 2023 for talk given strong grazing effects that year
alphadiv %>%
  filter(Year == 2023) %>%
  ggplot(., aes(Treatment, Observed, color=Treatment)) +
  scale_color_manual(values=cbPalette) +
  geom_jitter(alpha=0.2) +
  stat_summary(fun.data="mean_cl_boot", size=1) +
  theme_bw() +
  labs(y="Microbial Richness (ASV)", x="") +
  #theme(axis.ticks.x=element_blank(),
  #      axis.text.x = element_blank(),
  #      axis.title.x = element_blank()) +
  coord_cartesian(ylim=c(2750,4500)) +
  theme(legend.position = "none")




pdplot <- alphadiv %>%
  ggplot(., aes(Treatment, PD, color=Treatment)) +
  scale_color_manual(values=cbPalette) +
  geom_jitter(alpha=0.2) +
  stat_summary(fun.data="mean_cl_boot", size=1) +
  labs(y="Phylogenetic Diversity (Faith's D)", x="") +
  theme_bw() +
  theme(legend.position = "none") +
  coord_cartesian(ylim=c(250,370))

pdplot

ggsave(filename="pdplot.pdf", 
       path="/Users/mzaret/Documents/Research/KSU/Project 1/analysis/figures",
       width=6, height = 6)

alphadiv %>%
  ggplot(., aes(Treatment, PD, color=Treatment)) +
  scale_color_manual(values=cbPalette) +
  geom_jitter(alpha=0.2) +
  stat_summary(fun.data="mean_cl_boot", size=1) +
  theme_bw() +
  labs(y="PD", x="") +
  theme(axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank()) +
  theme(legend.position = "bottom") +
  facet_wrap(~Year)

alphadivplot <- richnessplot + pdplot + plot_layout(ncol=1, nrow=2)
alphadivplot

#merge PCOA with alpha div plot (code below)
all_plot + alphadivplot + plot_layout(ncol=2, nrow=1, widths=c(0.9,0.6))

#this is our figure 1 ^^^^^

#Merged samples for visualizing composition averaged for treatments or years####
#By Treatment
mergedps = merge_samples(ps, "Treatment")
mergedps
(sample_data(mergedps))

mergedps <- mergedps %>%
  ps_mutate(Treatment = case_when(Treatment == 1 ~ "UG20",
                                  Treatment == 2 ~ "G20",
                                  Treatment == 3 ~ "UG1",
                                  Treatment == 4 ~ "G1"))
mergedps.ra <- transform_sample_counts(mergedps, function(x) x / sum(x))

glom <- tax_glom(mergedps.ra, taxrank = 'Phylum')
glom # should list # taxa as # phyla
data_glom<- psmelt(glom) # create dataframe from phyloseq object
data_glom$Phylum <- as.character(data_glom$Phylum) #convert to character

#simple way to rename phyla with < 1% abundance
data_glom$Phylum[data_glom$Abundance < 0.01] <- "< 1% abund."

ggplot(data=data_glom, aes(x=Sample, y=100*Abundance, fill=Phylum)) +
  geom_bar(aes(), stat="identity", position="stack") +
  labs(x="Treatment", y="Percent Relative Abundance") +
  scale_fill_manual(values=c("green", "purple", "orange", "red", "blue", "yellow", "grey", "black",
                                     "pink", "brown", "green4", "purple4", "orange4")) +
  theme_bw()

#By Year                                         
mergedps = merge_samples(ps, "Year")
mergedps
(sample_data(mergedps))

mergedps <- mergedps %>%
  ps_mutate(Treatment = case_when(Year == 1 ~ "2016",
                                  Year == 2 ~ "2018",
                                  Year == 3 ~ "2020",
                                  Year == 4 ~ "2022",
                                  Year == 5 ~ "2023"))
mergedps.ra <- transform_sample_counts(mergedps, function(x) x / sum(x))

glom <- tax_glom(mergedps.ra, taxrank = 'Phylum')
glom # should list # taxa as # phyla
data_glom<- psmelt(glom) # create dataframe from phyloseq object
data_glom$Phylum <- as.character(data_glom$Phylum) #convert to character

#simple way to rename phyla with < 1% abundance
data_glom$Phylum[data_glom$Abundance < 0.01] <- "< 1% abund."

ggplot(data=data_glom, aes(x=Sample, y=100*Abundance, fill=Phylum)) +
  geom_bar(aes(), stat="identity", position="stack") +
  labs(x="Year", y="Percent Relative Abundance")

#RDA####
#Redundancy analysis
#filter data so we only have samples where we have soilchem data
ps2 <- ps_drop_incomplete(ps.hel)
ps2
ps2 <- ps2 %>%
  ps_filter(GWC != "NA") %>%
  ps_filter(C != ".")

ps2

#environmental data
env <- as.data.frame(as.matrix(sample_data(ps2)))
env
env <- env[,13:16]
str(env)
env <- env %>%
  mutate(C = as.numeric(C),
         N = as.numeric(N),
         pH = as.numeric(pH),
         GWC = as.numeric(GWC))

#need to convert C and N to moles - then probably calculate C:N to remove collinearity in rda
env <- env %>%
  mutate(C = C*0.084,
         N = N*0.07139,
         C_N = C/N) %>%
  select(pH, GWC, C, N, C_N)

asv_hell <- otu_table(ps2)

env <- scale(env, center=TRUE, scale=TRUE) #standardize env data

env <- as.data.frame(env)

dca <- decorana(asv_hell)
dca
#if axis length of DCA1 is <3 RDA, if >4, CCA - between go with whatever you want


#dimensions match thanks to our filtering in making ps2 object
dim(env)
dim(asv_hell)

PC1=rda(asv_hell~., env)
summary(PC1)
#anova.cca(PC1, permutations=1000, by="term")

#visualize rda
#extract sites
x <- summary(PC1)
sites <- as.data.frame(x$sites) %>%
  select(RDA1, RDA2) %>%
  mutate(SampleID = row.names(.))

#extract vectors
PC1$CCA$biplot
vscores <- data.frame(PC1$CCA$biplot) %>%
  select(RDA1, RDA2)

ps2 <- ps2 %>%
  ps_join(., sites, by="SampleID", type="left")

df <- as.data.frame(as.matrix(sample_data(ps2)))
str(df)

df <- df %>%
  mutate(RDA1 = as.numeric(RDA1),
         RDA2 = as.numeric(RDA2),
         C = as.numeric(C),
         N = as.numeric(N),
         pH = as.numeric(pH),
         GWC = as.numeric(GWC),
         Year = as.factor(Year),
         Treatment = as.factor(Treatment),
         Graze = as.factor(Graze),
         Burn = as.factor(Burn),
         Plot = as.factor(Plot)) %>%
  mutate(C = C*0.084,
         N = N*0.07139,
         C_N = C/N)

levels(df$Treatment)
df$Treatment <- relevel(df$Treatment, ref='G1')
df$Treatment <- relevel(df$Treatment, ref='UG1')
df$Treatment <- relevel(df$Treatment, ref='G20')
df$Treatment <- relevel(df$Treatment, ref='UG20')
levels(df$Treatment)

ggplot(df, aes(x=RDA1, y=RDA2, color=Treatment)) +
  geom_point() +
  geom_text(data = vscores, aes(x = RDA1, y = RDA2, label = rownames(vscores)), col = 'black') +
  geom_segment(data = vscores, aes(x = 0, y = 0, xend = RDA1, yend = RDA2), arrow=arrow(length=unit(0.2,"cm")),
               alpha = 0.75, color = 'black')+
  scale_color_manual(values=cbPalette) +
  theme_bw() +
  theme(strip.text.y = element_text(angle = 0))

#variance partitioning
#need to bring in lat long info as geo data
geo <- as.data.frame(as.matrix(sample_data(ps2)))
geo
geo <- geo[,11:12]
geo

#add year#
year <- as.data.frame(as.matrix(sample_data(ps2)))
year
year <- year[,2]
year <- as.factor(year)

spe.part.all <- varpart(asv_hell, env, geo, year)
spe.part.all

plot(spe.part.all,
     Xnames = c("Env", "Geo", "Year"), # name the partitions
     bg = c("seagreen3", "mediumpurple", "red"), alpha = 80, # colour the circles
     digits = 2, # only show 2 digits
     cex = 1.5)

#We can make a nicer one than that
#custom vpa plot
library(ggforce)

df.venn <- data.frame(x = c(3, 1, 2),y = c(1, 1,2.8),labels = c('Year', 'Env',"Geo"))
p <- ggplot(df.venn, aes(x0 = x, y0 = y, r = 1.5, fill = df.venn$labels)) +geom_circle(alpha = .5, size = 2, colour = 'black',show.legend = FALSE ) +
  coord_fixed() +
  scale_fill_manual(values=c("white","purple", "black")) #annotate("text", x = df.venn$x , y = df.venn$y,label=df.venn$labels ,size = 5)
p

p <- p+
  annotate("text", x=2, y=4.5, label="Geographic Distance", size=6) +
  annotate("text", x=1, y=-0.75, label="Environment", size=6) +
  annotate("text", x=3, y=-0.75, label="Year", size=6) +
  annotate("text", x=3.25, y=1, label="8.1%", size=6) +
  annotate("text", x=0.6, y=1, label="2.4%", size=6) +
  annotate("text", x=2, y=3, label="7.4%", size=6) +
  annotate("text", x = 2 , y =1,label="6.4%" ,size = 4)+
  annotate("text", x = 1.35 , y =2,label="2.3%" ,size = 4)+
  annotate("text", x = 2.7 , y =2,label="1.0%" ,size = 4) +
  theme_void()
p

ggsave(filename="vpaplot.pdf", 
       path="/Users/mzaret/Documents/Research/KSU/Project 1/analysis/figures",
       width=6, height = 6)

envplot <- phplot + cnplot + plot_layout(ncol=1, nrow=2)
figs14 <- p + envplot
figs14

#partial mantel tests#####
#looking at correlations between community dissimilarity and geographic distance between samples while accounting for environmental differences (RDA axes)
#function to remove 0 values when comparing matrices
select_all_but_diag <- function(x) matrix(x[lower.tri(x, diag = F) | upper.tri(x, diag = F)], nrow = nrow(x) - 1, ncol = ncol(x))

#UG20
ug20<- ps2 %>%
  ps_filter(Treatment == "UG20")

#bray curtis distance
bray <- phyloseq::distance(ug20, method="bray")
bray_matrix <- as.matrix(bray)

bray_matrix <- select_all_but_diag(bray_matrix)


#geographic distance
latlong <- df %>%
  filter(Treatment == "UG20") %>%
  select(LatDD, LongDD) %>%
  as.matrix(.)
dist_matrix <- geodist(latlong, measure="cheap")
dist_matrix <- log(dist_matrix+1)
dist_matrix <- select_all_but_diag(dist_matrix)

#RDA1 distance
RDA1_matrix <- df %>%
  filter(Treatment == "UG20") %>%
  select(RDA1) %>%
  as.matrix(.) %>%
  dist(., method="euclidean") %>%
  as.matrix(.)
RDA1_matrix <- select_all_but_diag(RDA1_matrix)


#RDA2 distance
RDA2_matrix <- df %>%
  filter(Treatment == "UG20") %>%
  select(RDA2) %>%
  as.matrix(.) %>%
  dist(., method="euclidean") %>%
  as.matrix(.)
RDA2_matrix <- select_all_but_diag(RDA2_matrix)


mantel(bray_matrix, dist_matrix)
mantel.partial(bray_matrix, dist_matrix, RDA1_matrix)
mantel.partial(bray_matrix, dist_matrix, RDA2_matrix)

#UG1
ug1<- ps2 %>%
  ps_filter(Treatment == "UG1")

#bray curtis distance
bray <- phyloseq::distance(ug1, method="bray")
bray_matrix <- as.matrix(bray)
bray_matrix <- select_all_but_diag(bray_matrix)


#geographic distance
latlong <- df %>%
  filter(Treatment == "UG1") %>%
  select(LatDD, LongDD) %>%
  as.matrix(.)
dist_matrix <- geodist(latlong, measure="cheap")
dist_matrix <- log(dist_matrix+1)
dist_matrix <- select_all_but_diag(dist_matrix)

#RDA1 distance
RDA1_matrix <- df %>%
  filter(Treatment == "UG1") %>%
  select(RDA1) %>%
  as.matrix(.) %>%
  dist(., method="euclidean") %>%
  as.matrix(.)
RDA1_matrix <- select_all_but_diag(RDA1_matrix)

#RDA2 distance
RDA2_matrix <- df %>%
  filter(Treatment == "UG1") %>%
  select(RDA2) %>%
  as.matrix(.) %>%
  dist(., method="euclidean") %>%
  as.matrix(.)
RDA2_matrix <- select_all_but_diag(RDA2_matrix)

mantel(bray_matrix, dist_matrix)
mantel.partial(bray_matrix, dist_matrix, RDA1_matrix)
mantel.partial(bray_matrix, dist_matrix, RDA2_matrix)


#G20
g20<- ps2 %>%
  ps_filter(Treatment == "G20")

#bray curtis distance
bray <- phyloseq::distance(g20, method="bray")
bray_matrix <- as.matrix(bray)
bray_matrix <- select_all_but_diag(bray_matrix)

#geographic distance
latlong <- df %>%
  filter(Treatment == "G20") %>%
  select(LatDD, LongDD) %>%
  as.matrix(.)
dist_matrix <- geodist(latlong, measure="cheap")
dist_matrix <- log(dist_matrix+1)
dist_matrix <- select_all_but_diag(dist_matrix)

#RDA1 distance
RDA1_matrix <- df %>%
  filter(Treatment == "G20") %>%
  select(RDA1) %>%
  as.matrix(.) %>%
  dist(., method="euclidean") %>%
  as.matrix(.)
RDA1_matrix <- select_all_but_diag(RDA1_matrix)

#RDA2 distance
RDA2_matrix <- df %>%
  filter(Treatment == "G20") %>%
  select(RDA2) %>%
  as.matrix(.) %>%
  dist(., method="euclidean") %>%
  as.matrix(.)
RDA2_matrix <- select_all_but_diag(RDA2_matrix)

mantel(bray_matrix, dist_matrix)
mantel.partial(bray_matrix, dist_matrix, RDA1_matrix)
mantel.partial(bray_matrix, dist_matrix, RDA2_matrix)


#G1
g1<- ps2 %>%
  ps_filter(Treatment == "G1")

#bray curtis distance
bray <- phyloseq::distance(g1, method="bray")
bray_matrix <- as.matrix(bray)
bray_matrix <- select_all_but_diag(bray_matrix)

#geographic distance
latlong <- df %>%
  filter(Treatment == "G1") %>%
  select(LatDD, LongDD) %>%
  as.matrix(.)
dist_matrix <- geodist(latlong, measure="cheap")
dist_matrix <- log(dist_matrix+1)
dist_matrix <- select_all_but_diag(dist_matrix)

#RDA1 distance
RDA1_matrix <- df %>%
  filter(Treatment == "G1") %>%
  select(RDA1) %>%
  as.matrix(.) %>%
  dist(., method="euclidean") %>%
  as.matrix(.)
RDA1_matrix <- select_all_but_diag(RDA1_matrix)


#RDA2 distance
RDA2_matrix <- df %>%
  filter(Treatment == "G1") %>%
  select(RDA2) %>%
  as.matrix(.) %>%
  dist(., method="euclidean") %>%
  as.matrix(.)
RDA2_matrix <- select_all_but_diag(RDA2_matrix)

mantel(bray_matrix, dist_matrix)
mantel.partial(bray_matrix, dist_matrix, RDA1_matrix)
mantel.partial(bray_matrix, dist_matrix, RDA2_matrix)


#Soil Chemistry results####
df <- as.data.frame(as.matrix(sample_data(ps2)))
str(df)
#have to transform and relevel variables again..same patterns that we did for the RDA
df <- df %>%
  mutate(C = as.numeric(C),
         N = as.numeric(N),
         pH = as.numeric(pH),
         GWC = as.numeric(GWC),
         Year = as.factor(Year),
         Treatment = as.factor(Treatment),
         Graze = as.factor(Graze),
         Burn = as.factor(Burn),
         Plot = as.factor(Plot)) %>%
  mutate(C = C*0.084,
         N = N*0.07139,
         C_N = C/N)

#reorder factors for analysis (again)
levels(df$Year)

levels(df$Treatment)
df$Treatment <- relevel(df$Treatment, ref='G1')
df$Treatment <- relevel(df$Treatment, ref='UG1')
df$Treatment <- relevel(df$Treatment, ref='G20')
df$Treatment <- relevel(df$Treatment, ref='UG20')
levels(df$Treatment)


levels(df$Graze)
df$Graze <- relevel(df$Graze, ref='Bison Present')
df$Graze <- relevel(df$Graze, ref='No Grazers')
levels(df$Graze)

levels(df$Burn)
df$Burn <- relevel(df$Burn, ref='Frequent Burn')
df$Burn <- relevel(df$Burn, ref='Infrequent Burn')
levels(df$Burn)

#models
#pH
hist(df$pH)
lme1 <- lme(fixed = pH ~ Graze * Burn * Year, 
            random = ~1 | Plot, data=df)
plot(lme1)
summary(lme1)
anova(lme1)

lsmeans(lme(fixed = pH ~ Graze * Burn * Year, 
            random = ~1 | Plot, data=df),
        pairwise ~ Graze|Year, type="Tukey")

#gwc
hist(log(df$GWC))

lme1 <- lme(fixed = (GWC) ~ Graze * Burn * Year, 
            random = ~1 | Plot, data=df)
plot(lme1)
summary(lme1)
anova(lme1)

lsmeans(lme(fixed = GWC ~ Graze * Burn * Year, 
            random = ~1 | Plot, data=df),
        pairwise ~ Graze|Year, type="Tukey")
#C_N
hist((df$C_N))

lme1 <- lme(fixed = C_N ~ Graze * Burn * Year, 
            random = ~1 | Plot, data=df)
plot(lme1)
summary(lme1)
anova(lme1)

lsmeans(lme(fixed = C_N ~ Graze * Burn * Year, 
            random = ~1 | Plot, data=df),
        pairwise ~ Graze*Burn|Year, type="Tukey")

#C
hist((df$C))

lme1 <- lme(fixed = (C) ~ Graze * Burn * Year, 
            random = ~1 | Plot, data=df)
plot(lme1)
summary(lme1)
anova(lme1)

lsmeans(lme(fixed = C ~ Graze * Burn * Year, 
            random = ~1 | Plot, data=df),
        pairwise ~ Graze|Year, type="Tukey")

#N
hist(log(df$N))

lme1 <- lme(fixed = log(N) ~ Graze * Burn * Year, 
            random = ~1 | Plot, data=df)
plot(lme1)
summary(lme1)
anova(lme1)

lsmeans(lme(fixed = N ~ Graze * Burn * Year, 
            random = ~1 | Plot, data=df),
        pairwise ~ Graze|Year, type="Tukey")

##visualize soil chemistry trends
phplot <- df %>%
  ggplot(., aes(Treatment, pH, color=Treatment)) +
  scale_color_manual(values=cbPalette) +
  geom_jitter(alpha=0.2) +
  stat_summary(fun.data="mean_cl_boot", size=1) +
  labs(y="pH", x="") +
  theme_bw() +
  #theme(axis.ticks.x=element_blank(),
  #      axis.text.x = element_blank(),
  #      axis.title.x = element_blank()) +
  theme(legend.position = "none")

phplot

df %>%
  ggplot(., aes(Treatment, pH, color=Treatment)) +
  scale_color_manual(values=cbPalette) +
  geom_jitter(alpha=0.2) +
  stat_summary(fun.data="mean_cl_boot", size=0.75) +
  labs(y="soil pH", x="") +
  theme_bw() +
  theme(legend.position = "bottom") +
  facet_wrap(~Year)

gwcplot <- df %>%
  ggplot(., aes(Treatment, GWC, color=Treatment)) +
  scale_color_manual(values=cbPalette) +
  geom_jitter(alpha=0.2) +
  stat_summary(fun.data="mean_cl_boot", size=1) +
  labs(y="Water Content (Percent)", x="") +
  #theme(axis.ticks.x=element_blank(),
  #      axis.text.x = element_blank(),
  #      axis.title.x = element_blank()) +
  theme(legend.position = "none")

gwcplot

df %>%
  ggplot(., aes(Treatment, GWC, color=Treatment)) +
  scale_color_manual(values=cbPalette) +
  geom_jitter(alpha=0.2) +
  stat_summary(fun.data="mean_cl_boot", size=0.75) +
  labs(y=~paste("soil field water content ", "(", "g ", "g"^-1,")"), x="") +
  theme_bw() +
  theme(legend.position = "bottom") +
  facet_wrap(~Year)


cnplot <- df %>%
  ggplot(., aes(Treatment, C_N, color=Treatment)) +
  scale_color_manual(values=cbPalette) +
  geom_jitter(alpha=0.2) +
  stat_summary(fun.data="mean_cl_boot", size=1) +
  labs(y="Soil C:N Ratio", x="") +
  theme_bw() +
  #theme(axis.ticks.x=element_blank(),
  #      axis.text.x = element_blank(),
  #      axis.title.x = element_blank()) +
  theme(legend.position = "none") 

cnplot

df %>%
  ggplot(., aes(Treatment, C_N, color=Treatment, fill=Treatment)) +
  scale_color_manual(values=cbPalette) +
  geom_jitter(alpha=0.2) +
  stat_summary(fun.data="mean_cl_boot", size=0.75) +
  labs(y="Soil total C:N", x="") +
  theme_bw() +
  theme(legend.position = "bottom") +
  facet_wrap(~Year)

cplot <- df %>%
  ggplot(., aes(Treatment, C, color=Treatment)) +
  scale_color_manual(values=cbPalette) +
  geom_jitter(alpha=0.2) +
  stat_summary(fun.data="mean_cl_boot", size=1) +
  labs(y="Soil Total Carbon", x="") +
  #theme(axis.ticks.x=element_blank(),
  #      axis.text.x = element_blank(),
  #      axis.title.x = element_blank()) +
  theme(legend.position = "none")

cplot

df %>%
  ggplot(., aes(Treatment, C, color=Treatment, fill=Treatment)) +
  scale_color_manual(values=cbPalette) +
  geom_jitter(alpha=0.2) +
  stat_summary(fun.data="mean_cl_boot", size=0.75) +
  labs(y="Soil %C", x="") +
  theme_bw() +
  theme(legend.position = "bottom") +
  facet_wrap(~Year)

nplot <- df %>%
  ggplot(., aes(Treatment, N, color=Treatment)) +
  scale_color_manual(values=cbPalette) +
  geom_jitter(alpha=0.2) +
  stat_summary(fun.data="mean_cl_boot", size=1) +
  labs(y="Soil Total N", x="")
  #theme(axis.ticks.x=element_blank(),
  #      axis.text.x = element_blank(),
  #      axis.title.x = element_blank()) +
  theme(legend.position = "none")

nplot

df %>%
  ggplot(., aes(Treatment, N, color=Treatment, fill=Treatment)) +
  scale_color_manual(values=cbPalette) +
  geom_jitter(alpha=0.2) +
  stat_summary(fun.data="mean_cl_boot", size=0.75) +
  labs(y="Soil %N", x="") +
  theme_bw() +
  theme(legend.position = "bottom") +
  facet_wrap(~Year)

phplot + cnplot

ggsave(filename="phcnplot.pdf", 
       path="/Users/mzaret/Documents/Research/KSU/Project 1/analysis/figures",
       width=8, height = 5)

gwcplot + cplot + nplot

#last is to combine VPA with some of our metrics for figure 3
envplot <- phplot + cnplot + plot_layout(ncol=1, nrow=2)
envplot
fig3 <- p + envplot
fig3

#differential abundance####
#look at what specific taxa are responding to treatments at ASV and phylum level
#i toggled this code back and forth between using ASV (ps) object or phylum aggregated tables (ps.phylum below)
library(DESeq2)

#combine asvs into phylum level
ps.phylum <- tax_glom(ps, taxrank="Phylum")
ps.phylum

#sort tax table by abundance
View(tax_table(tax_sort(ps.phylum, by=mean)))

#grazing impact
ds <- phyloseq_to_deseq2(ps, ~Graze) #toggle this with ps.phylum depending on focus
ds <- DESeq(ds)

alpha=0.001
res <- results(ds, alpha=alpha) #contrast=c("Graze", "Bison Present", "No Grazers"))
res = res[order(res$padj, na.last=NA), ]
res_sig = res[(res$padj < alpha), ]
res_sig
results <- (as.data.frame(res_sig))
results <- merge(results, tax_table(ps), by="row.names")
results %>%
  summarize(mean = mean(log2FoldChange))

reduced <- results %>%
  filter(log2FoldChange <0)
increased <- results %>%
  filter(log2FoldChange >0)

grazing_asv_distribution <- ggplot(results, aes(log2FoldChange, fill=Phylum)) +
  geom_histogram(bins=100) +
  geom_vline(xintercept=0, lty=2) +
  geom_vline(xintercept=-0.253, lty=1, col="red") +
  theme_bw() +
  labs(x="Change in ASV Abundance (log2)", y="Count") +
  scale_y_continuous(expand = c(0,0), limits=c(0,210), breaks=c(0,50,100,150,200))
grazing_asv_distribution       

res_sig = cbind(as(res_sig, "data.frame"), as(tax_table(ps)[rownames(res_sig), ], "matrix"))

#rerun deseq analysis with ps phylum for these results
results <- results %>%
  arrange(log2FoldChange) %>%    # First sort by val. This sort the dataframe but NOT the factor levels
  mutate(Phylum=factor(Phylum, levels=Phylum))   # This trick update the factor levels
  
gdplot <- ggplot(results, aes(x=(Phylum), y=log2FoldChange)) +
  geom_point(size=3) +
  geom_hline(yintercept=0, lty=2) +
  scale_y_continuous(limits=c(-2.5,2.5)) +
  ggtitle("Grazing effect") +
  theme_bw() +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, face="bold"))

gdplot

#(2^log2fc -1) *100 to calculate percent change
#Crenarchaeota
(2^1.36532824 - 1) * 100
#Firmicutes
(2^0.87862550 - 1) * 100

#fire treatment analysis now
ds <- phyloseq_to_deseq2(ps, ~Burn)
ds <- DESeq(ds)

alpha=0.001
res <- results(ds, alpha=alpha) #contrast=c("Graze", "Bison Present", "No Grazers"))
res = res[order(res$padj, na.last=NA), ]
res_sig = res[(res$padj < alpha), ]
res_sig
results <- (as.data.frame(res_sig))
results <- merge(results, tax_table(ps), by="row.names")
results %>%
  summarize(mean = mean(log2FoldChange))

burn_asv_distribution <- ggplot(results, aes(log2FoldChange)) +
  geom_histogram(bins=100) +
  geom_vline(xintercept=0, lty=2) +
  geom_vline(xintercept=0.006, lty=1, col="red") +
  labs(x="Change in ASV Abundance (log2)", y="Count") +
  theme_bw() +
  scale_y_continuous(expand=c(0,0), limits=c(0,210), breaks=c(0,50,100,150,200))

burn_asv_distribution  


results <- results %>%
  arrange(log2FoldChange) %>%    # First sort by val. This sort the dataframe but NOT the factor levels
  mutate(Phylum=factor(Phylum, levels=Phylum))   # This trick update the factor levels

results

fdplot <- ggplot(results, aes(x=(Phylum), y=log2FoldChange)) +
  geom_point(size=3) +
  geom_hline(yintercept=0, lty=2) +
  scale_y_continuous(limits=c(-3,2.5)) +
  ggtitle("Fire effect") +
  theme_bw() +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, face="bold"))

fdplot

#(2^log2fc -1) *100 to calculate percent change
#Proteobacteria
(2^0.1761 - 1) * 100
#Bacteroidota
(2^0.27941823 - 1) * 100

gdplot + fdplot

grazing_asv_distribution + labs(title="Grazing effect") + burn_asv_distribution +
  labs(title="Burn Effect")
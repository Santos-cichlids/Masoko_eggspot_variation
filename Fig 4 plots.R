library(tidyverse)
library (ggbeeswarm)
library(ggpubr)

####settings ----------

#For customising graphs themes
theme <- theme(
  legend.position = "none",
  axis.text = element_text(color = "black"),
  axis.line = element_line(color = "black",linewidth = 0.8),
  axis.ticks = element_line(color = "black"))


#colours for plots
sibling <- "#528F77"
albino <- "#FFB701"


#for size and appearance of plot geoms
violinwidth = 0.8
boxwidth = 0.1
jitterwidth = 0.07

violinalpha = 0.3
boxalpha = 0.6
jitteralpha = 0.2
pointsize <- 2
linesize <- 0.4

fourshapes <- c(16, 17, 15, 18)
fiveshapes <- (c(20, 16, 17, 15, 18))
nineshapes <- c(15, 16, 17, 18, 15, 16, 17, 18, 15)
scale = "count"


###common data ----------
#import standard length data
SL_import <- read_csv("mutant clutches SL grouped days.csv")

#make dates workable
SL_import$Date <- dmy(SL_import$Date)

#create unique identifiers for each row matching between datasets
SL_import$Indiv_date <- paste(SL_import$Individual, SL_import$Date, sep="_")

###iridophore coverage data -------
#import iridophore coverage meaurements
irid_cov <- read_csv("Iridophore coverage 2024 .csv")

#make dates workable
irid_cov$Date <- dmy(irid_cov$Date)

#create unique identifiers for each row matching between datasets
irid_cov$Indiv_date <- paste(irid_cov$Individual, irid_cov$Date, sep="_")

#join datasets by unique identifiers 
irid_cov_SL <- irid_cov %>% left_join(SL_import, by = "Indiv_date")

##remove NA values from growth column
irid_cov_SL <- irid_cov_SL[!is.na(irid_cov_SL$Clutch), ]

### iridophore coverage plot -------------------------------

ggplot(irid_cov_SL, mapping = aes(x= Dpf_estimated.x, y=Percent_area, 
                                               group_by = Individual.x, 
                                             colour = Pheno.x, shape = Clutch))+
  #add data points
  geom_point(size=3, alpha=0.6)+
  #custom scales
  scale_colour_manual(values = c(albino, sibling))+
  scale_shape_manual(values=fiveshapes)+
  #labels and theme
  ylab("Percentage of cropped flank area above value threshold")+
  xlab("dpf")+
  xlim(10, 25)+
  theme_classic(base_size = 16)+
  theme

#save high res
ggsave(filename="irid coverage clutch by dpf vertical .pdf", path="irid_plots", width=2.7, height=4)



## early aggregation timing data--------------

#import event dated timing data
events <- read_csv("early ag events.csv")

events$Event <- factor(events$Event, levels=c("stellate associations 3-5", 
                                              "overlaying 3-5",
                                              "iridophores 4-5"))

#create unique identifiers for each row matching between datasets
events$Indiv_date <- paste(events$Individual, events$Date, sep="_")

#remove duplicate columns from SL dataframe
SL = subset(SL_import, select = c(Standard_length_mm, Indiv_date))

#join datasets by unique identifiers 
events <- events %>% left_join(SL, by = "Indiv_date")


### early aggregation event timing plot------------------------

#violin plots for event timing between amelanistics and siblings
ggplot(data = events, mapping = aes(x=Day, 
                                    y=Event,
                                    colour=Pheno, fill=Pheno))+
  #add violins
  geom_violin(alpha=violinalpha, width=violinwidth, trim=FALSE, scale = scale)+
  
  # Boxplot overlay (narrow)
  geom_boxplot(
    aes(fill = Pheno),
    width = boxwidth,
    alpha = boxalpha,
    outlier.shape = NA,
    color = "black",
    linewidth = linesize
  )+
  
  #add data points, jittered by quasirandom
  geom_quasirandom(method = 'quasirandom',
                   alpha = jitteralpha, width = 0.04, 
                   size=2,stroke=0.3,
                   mapping=aes(colour=Pheno), dodge.width = 0.9)+
  #custom scales
  scale_color_manual(values = c(albino, sibling))+
  scale_fill_manual(values=c(albino, sibling))+
  
  #theme and labels
  theme_classic(base_size = 16)+
  theme+
  scale_y_discrete(labels = function(x) lapply(strwrap(x, width = 10, simplify = FALSE), paste, collapse="\n"))

#save high res
ggsave(filename="early ag events with box .pdf", path="early_ag_events", width=6, height=4)

### early aggregation event timing stats-----------
#subsets to separate data for each event
irids <- subset (events, Event == "iridophores 4-5")
overlay <- subset (events, Event == "overlaying 3-5")
assoc <- subset (events, Event == "stellate associations 3-5")

#test distribution of each to determine parametric or non-parametric test  
shapiro.test(irids$Day)
#normally distributed 
shapiro.test(overlay$Day)
#normally distributed 
shapiro.test(assoc$Day)
#normally distributed 

#iridophore appearance test with ANOVA
res.aov2.irids <- aov(Day ~ Pheno * Sex, data = irids)
anova(res.aov2.irids)
aov2_irids <- anova(res.aov2.irids)
write.csv(aov2_irids, "aov2_agevents_irids .csv")

#overlaying test with ANOVA
res.aov2.overlay <- aov(Day ~ Pheno * Sex, data = overlay)
anova(res.aov2.overlay)
aov2_overlay <- anova(res.aov2.overlay)
write.csv(aov2_overlay, "aov2_agevents_overlay .csv")

#stellate association test with ANOVA
res.aov2.assoc <- aov(Day ~ Pheno * Sex, data = assoc)
anova(res.aov2.assoc)
aov2_assoc <- anova(res.aov2.assoc)
write.csv(aov2_assoc, "aov2_agevents_assoc .csv")

##early aggregation size comparisons --------
### early aggregation size comparisons data -------------------------

#import grouped by days measurements
SL_import <- read_csv("mutant clutches SL grouped days.csv")


#create a column for SL growth
SL_import$Growth <- (SL_import$Standard_length_mm)- (SL_import$Initial_SL)

#make dates workable
SL_import$Date <- dmy(SL_import$Date)


#remove NA values from SL column
SL_skip <- SL_import[!is.na(SL_import$Standard_length_mm), ]

#remove NA values from growth column
Growth_skip <- SL_import[!is.na(SL_import$Growth), ]


### subsets for initial aggregation plots
early_first <- subset(SL_import, 
                      Day == "40" | 
                        Day == "41" |
                        Day == "35")

early_first$Dpf_extrap <- early_first$Day + early_first$Start_dpf


###early aggregation calculations

#create column for area that includes whole aggregation
early_first$Area <- (early_first$Ag_length)*(early_first$Ag_width)
early_first$Area_relative <- (early_first$Ag_length)*(early_first$Ag_relative_width)

#create columns for length and area normalised by SL
early_first$nArea <- (early_first$Area)/(early_first$Standard_length_mm)
early_first$nAg_length <- (early_first$Ag_length)/(early_first$Standard_length_mm)
early_first$nArea_relative <- (early_first$Area_relative)/(early_first$Standard_length_mm)

#early aggregations size plots--------


#check SL not sig diff
ggboxplot(early_first, x = "Clutch", y = "Standard_length_mm", color = "Pheno",
          palette = c(albino, sibling))
ggdensity(early_first$Standard_length_mm)
ggqqplot(early_first$Standard_length_mm)
shapiro.test(early_first$Standard_length_mm)
#normal distributin, parametric test ok 

res.aov3.earlyag.SL <- aov(Standard_length_mm ~ Pheno * Day * Sex, data = early_first)
aov3_earlyag_SL <- anova(res.aov3.earlyag.SL)
aov3_earlyag_SL
###no sig diff in SL so no need to SL-normalise this section


##early ag length
 ggplot(early_first, mapping = aes(x=Pheno, y=Ag_length,
                     colour = Pheno, fill =Pheno,
                     shape=Clutch))+
   #violin
  geom_violin(alpha=violinalpha, width = violinwidth, trim=F, scale = "area", 
              mapping=aes(group=Pheno))+
   #boxplot
   geom_boxplot(
     aes(fill = Pheno, group=Pheno),
     width = boxwidth,
     alpha = boxalpha,
     outlier.shape = NA,
     color = "black",
     linewidth = linesize)+
   #data points jittered by quasirandom
  geom_quasirandom(method = 'quasirandom',alpha = jitteralpha, width = 0.1,size=3,stroke=0.5)+
   #custom scales
  scale_colour_manual(values = c(albino, sibling))+
  scale_fill_manual(values = c(albino, sibling))+
  scale_shape_manual(values=c(fourshapes))+
   #labels and theme
  ylab("Aggregation length (mm)")+
  theme_classic(base_size = 16)+
  theme

 #save
ggsave(filename="early ag length with box.pdf", path="early_ag_plots", width=5.25, height=5)

##early aggregation length stats ------------
#test normality
ggdensity(early_first$Ag_length)
ggqqplot(early_first$Ag_length)
shapiro.test(early_first$Ag_length)
#normal

#ANOVA
res.aov3.earlyag.length <- aov(Ag_length ~ Pheno * Day * Sex, data = early_first)
aov3_earlyag_length <- anova(res.aov3.earlyag.length)
aov3_earlyag_length
write.csv(aov3_earlyag_length, "aov3_earlyag_length .csv")


### spot number data ------------

#spot number calculations and wrangling

#create a column for SL growth
SL_import$Growth <- (SL_import$Standard_length_mm)- (SL_import$Initial_SL)


#remove NA values from SL column
SL_skip <- SL_import[!is.na(SL_import$Standard_length_mm), ]

#remove NA values from growth and sex columns
Growth_skip <- SL_import[!is.na(SL_import$Growth), ]
SL <- SL_import[!is.na(SL_import$Sex), ]

#change egg-spot number format to tidy table, leaving other columns unchanged
SL <- gather(SL, maturity, number, -Individual, -Date, 
             -Day, -Grouped_day, -Pheno, -Geno, -Growth, -Standard_length_mm, 
             -Sex, -Clutch, -Start_date, -Initial_SL, -Start_date,
             -Stage_CFR, -Dpf_estimated, -Start_dpf,
             -Isolation, -Line, -Ag_length, -Ag_width, -Ag_actual_width,
             -Ag_normaliser_width, -Ag_relative_width, -Irids_present, -SL_type,
             -Indiv_date)

#created number normalised by SL value
SL$nNumber <- SL$number/SL$Standard_length_mm

#calculate dpf
SL$Dpf_extrap <- SL$Day + SL$Start_dpf

#remove unmatched datapoints
SL <- subset (SL, Day < 170 | Day > 172)

###spot number subsets to isolate final spot numbers in males

SL_spots <- subset (SL,
                    maturity == "Spots")

SL_spotsM <- subset (SL_spots, Sex == "M")

SL_late_spotsM <- subset (SL_spotsM, Day > 180)


###spot number plots -----------------

#plot of day 200 male spot number
ggplot(SL_late_spotsM, mapping = aes(x=number, y=Pheno,
                                     colour = Pheno, fill =Pheno))+
  #violin
  geom_violin(alpha=violinalpha, width=violinwidth, trim=FALSE, scale = "count")+
  #boxplot
  geom_boxplot(
    aes(fill = Pheno, group=Pheno),
    width = boxwidth,
    alpha = boxalpha,
    outlier.shape = NA,
    color = "black",
    linewidth = linesize)+
  #data points jittered by quasirandom
  geom_quasirandom(size=pointsize, alpha=jitteralpha, 
                   mapping = aes ( shape = Clutch), width=0.15)+
  #custom scales
  scale_colour_manual(values = c(albino, sibling))+
  scale_fill_manual(values = c(albino, sibling))+
  scale_shape_manual(values=fourshapes)+
  #theme and labels
  theme+
  theme_classic(base_size = 16)+
  theme(legend.position="none",
        axis.title.x = element_blank())+
  xlim(0,5)+
  coord_flip()

#save
ggsave(filename="d200 M with box.pdf", path="spot_plots", width=5.25, height=5)



###proximal edge measurements data-------------------------

#import grouped by days measurements
prox_all <- read_csv("mutant clutches prox spot.csv")

#create a column for SL growth
prox_all$Growth <- (prox_all$Standard_length_mm)- (prox_all$Initial_SL)

#make dates workable
prox_all$Date <- dmy(prox_all$Date)

#remove NA values 
prox_all <- prox_all[!is.na(prox_all$Prox_edge), ]

#create a column for proportion (centre)
prox_all$Prox_proportion <- ((prox_all$Prox_calc_centre) / (prox_all$Para_full_length)) * 100

#create a column for proportion (edge)
prox_all$Prox_proportion_edge <- ((prox_all$Prox_edge) / (prox_all$Para_full_length)) * 100

#subset by timepoint
prox_200 <- subset(prox_all, Day > 190)
prox <- subset(prox_all, Day < 170)


### proximal measurements plots ------------

#plot measurments of the distances between the proximal 
#edge of the fin and the proximal edge of the first aggregation

ggplot(prox_200, mapping = aes(x=Pheno, y=Prox_edge,
                               colour = Pheno, fill = Pheno))+
  #violin
  geom_violin(alpha=violinalpha, width=violinwidth, trim=FALSE, scale = scale)+
  #boxplot
  geom_boxplot(
    aes(fill = Pheno, group=Pheno),
    width = boxwidth,
    alpha = boxalpha,
    outlier.shape = NA,
    color = "black",
    linewidth = linesize)+
  #data points
  geom_point(size=pointsize, alpha=jitteralpha, mapping = aes ( shape = Clutch))+
  #custom scales
  scale_colour_manual(values = c(albino, sibling))+
  scale_fill_manual(values = c(albino, sibling))+
  #labels and themes
  ylab("Distance (mm)")+
  ylim(5,-0.1)+
  scale_shape_manual(values = c(fourshapes))+
  theme+
  theme_classic(base_size = 16)+
  theme(legend.position="none",
        axis.title.x = element_blank())+
  #stats
geom_signif(comparisons = list(c("albino", "sibling")),   
            map_signif_level=TRUE, colour="black", y_position=-0.1,
            tip_length=0)+

#save
ggsave(filename="200_prox_edge with box.pdf", path="prox_plots", width=5.25, height=5)




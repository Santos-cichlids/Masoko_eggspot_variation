library(tidyverse)
library (ggbeeswarm)
library(ggpubr)
library(RColorBrewer)

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

female<-"#8D8A91"
male<-"#BF675E"

D <- "palegreen3"
G <-"skyblue3"
H <-"tan2"
I <-"dodgerblue4"

brewer.pal(n=7,"Set2")
segpal <- c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494",
            "#66C2D9", "#FC8C22", "#8DB0FB", "#E78AF3", "#A6B011", "#FFB92A", "#E5A494")

segpalrev <- c("#E5A494", "#FFB92A", "#A6B011", "#E78AF3", "#8DB0FB", "#FC8C22", "#66C2D9",
               "#E5C494", "#FFD92F", "#A6D854", "#E78AC3", "#8DA0CB", "#FC8D62", "#66C2A5")

segpalrev10 <-  c("#8DB0FB", "#FC8C22", "#66C2D9",
                  "#E5C494", "#FFD92F", "#A6D854", "#E78AC3", "#8DA0CB", "#FC8D62", "#66C2A5")

segpalrev12 <- c("#A6B011", "#E78AF3", "#8DB0FB", "#FC8C22", "#66C2D9",
                 "#E5C494", "#FFD92F", "#A6D854", "#E78AC3", "#8DA0CB", "#FC8D62", "#66C2A5")


#for size and appearance of plot geoms
violinwidth = 0.8
boxwidth = 0.1
jitterwidth = 0.07

alpha=0.7
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




#import growth measurements

lengths <- read_csv("albino_growth.csv")
prop <- read_csv("albino_propgrowth.csv")


#make dates workable
lengths$Date <- dmy(lengths$Date)


#change format to tidy table
growth <- gather(lengths, measurement, value, 
                 -Individual, -Date, -Day, -Grouped_day, -Sex, -Geno,
                 -Pheno, -Clutch)

#remove NA values from value column
growth <- growth[!is.na(growth$value), ]

#remove NA values from lengths df FR6 col
lengths <- lengths[!is.na(lengths$Ray6_length), ]

growth$Indivmeas <-  paste(growth$Individual, growth$measurement, sep="_")

growth <- growth %>%
  group_by(Indivmeas) %>%
  mutate(recentgrowth = value - lag(value, default = first(value), order_by = Date))



whole_lengths <- subset(growth, 
                        measurement == "Body_fin_connect_length" | 
                          measurement == "Ray6_length" |
                          measurement == "Standard_length_mm")



####how the albino_propgrowths.csv was made 
#propgrowth <- whole_lengths
#propgrowth <- propgrowth %>% subset(Grouped_day == "a_start" |
#                                     Grouped_day == "h_end")
#propgrowth <- propgrowth[!is.na(propgrowth$Sex), ]
#write_csv(propgrowth, "albino_propgrowth.csv")
####In excel sort separate into two columns: start and end (quicker than troubleshooting pivot wider)


###SL through time data -----------


#create a column for SL growth
SL_import$Growth <- (SL_import$Standard_length_mm)- (SL_import$Initial_SL)


#remove NA values from SL column
SL_skip <- SL_import[!is.na(SL_import$Standard_length_mm), ]


#remove NA values from growth column
Growth_skip <- SL_import[!is.na(SL_import$Growth), ]

### SL through time ----------------

#growth v Date by Clutch
#SL by day by clutch
ggplot(SL_skip, mapping = aes(x=Day, 
                            y=Standard_length_mm, 
                             group_by = Individual, 
                             colour = Clutch))+
  #lines
  geom_line(linewidth=linesize, alpha=alpha)+
  #custom scales
  scale_colour_manual(values = c(D,G,H,I))+
  #theme and labels
  theme_classic(base_size = 16)+ 
  theme

ggsave(filename="growth_clutch_day_sq s7.pdf", path="sl_plots", units=c("mm"), width=84.1, height=57.6)


#SL  v Day by phenotype
ggplot(Growth_skip, 
                   mapping = aes(x=Day, 
                          y=Standard_length_mm, 
                          group_by = Individual, 
                          colour = Pheno))+
  #line
  geom_line(linewidth=linesize, alpha=alpha)+
  #custom scale
  scale_colour_manual(values = c(albino, sibling))+
  #theme
  theme_classic(base_size = 16)+ 
  theme

ggsave(filename="growth_pheno_day_sq s7.pdf", path="sl_plots", units=c("mm"), width=84.1, height=57.6)


#SL v Day by sex
ggplot(Growth_skip, mapping = aes(x=Day, 
                                 y=Standard_length_mm, 
                                group_by = Individual,
                                colour = Sex))+
  #add lines
  geom_line(linewidth=linesize, alpha=alpha)+
  #custom scale
  scale_colour_manual(values = c(female, male))+
  #theme
  theme_classic(base_size = 16)+ 
  theme

ggsave(filename="growth_sex_day_sq s7.pdf", path="sl_plots", units=c("mm"), width=84.1, height=57.6)

#growth v Date by Clutch
ggplot(Growth_skip, mapping = aes(x=Date, y=Growth, 
                         group_by = Individual, 
                         colour = Clutch))+
  #lines
  geom_line(linewidth=linesize, alpha=alpha)+
  #custom scale
  scale_colour_manual(values = c(D,G,H,I))+
  #theme
  theme_classic(base_size = 16)+ 
  theme

ggsave(filename="growth_clutch_date_sq s7.pdf", path="sl_plots", units=c("mm"), width=84.1, height=57.6)



### fin growth data -------------------------------

#calculate proportional growth
prop$prop <- (prop$end_value / prop$start_value)
prop$actual <- (prop$end_value - prop$start_value)
prop_gathered <- gather(prop, growthcalcu, growthvalue, 
                        -Individual, -Date, -Day, -Grouped_day, -Sex, -Geno,
                        -Pheno, -Clutch, -start_value, -end_value, -measurement, 
                        -Indivmeas, -recentgrowth)

#male separate prop and actual dfs
prop_red <- prop_gathered %>% select(Individual, Date, Day, Grouped_day,
                                     Sex, Pheno, Geno, Clutch,
                                     growthcalcu, growthvalue, measurement )
prop_wide <- pivot_wider(prop_red, names_from = measurement, values_from = growthvalue)
prop_only <- subset(prop_wide, 
                    growthcalcu == "prop")
actual_only <- subset(prop_wide, 
                      growthcalcu == "actual")



###### fin growth plots ------------------

###AP

#Body_fin_connect_length
APstart <- subset(growth, 
                  measurement == "Body_fin_connect_length" &
                    Day == "0" )


APend <- subset(growth, 
                measurement == "Body_fin_connect_length" &
                  Day >140 & Day < 160)

APstart$start <- APstart$value
APend$end <- APend$value

AP <- APstart %>% left_join(APend, by = "Individual", suffix = c("", ".y")) %>%
  select(-ends_with(".y"))

AP$totalgrowth <- AP$end - AP$start


##anova total growth AP
aov3.AP.total.growth <- aov(totalgrowth ~ Pheno * Clutch * Sex, data=AP)
anova(aov3.AP.total.growth)
APtotalgrowth <- anova(aov3.AP.total.growth)
write.csv(APtotalgrowth, "aov3 AP total growth .csv")


##violin plots total growth AP

##pheno
ggplot(AP, mapping = aes(y=totalgrowth, x=Pheno, colour=Pheno, fill=Pheno,
                         shape = Clutch))+
  #violin
  geom_violin(alpha=0.3, size=0.5, trim=T, scale = scale, mapping=aes(group=Pheno))+
  #boxplot
  geom_boxplot(
    aes(fill = Pheno, group=Pheno),
    width = boxwidth,
    alpha = boxalpha,
    outlier.shape = NA,
    color = "black",
    linewidth = linesize)+
  #data points jittered by quasirandom
  geom_quasirandom(method = 'quasirandom',alpha = jitteralpha*2, width = 0.05, size=3,stroke=0.5,
                   mapping=aes(colour=Sex))+
  #custom scales
  scale_colour_manual(values = c(albino, female, male, sibling))+
  scale_fill_manual(values = c(albino, sibling))+
  scale_shape_manual(values=c(fourshapes))+
  #theme and labels
  theme_classic(base_size = 16)+ 
  theme+
  coord_flip()

ggsave(filename="AP total growth Pheno s7.pdf", path="albsib_fin_growth_plots",
       width=5.25, height=5)


##clutch
ggplot(AP, mapping = aes(y=totalgrowth, x=Clutch, colour=Clutch, fill=Clutch,
                         shape = Clutch))+
  #violin
  geom_violin(alpha=0.3, size=0.5, trim=T, scale = scale, mapping=aes(group=Clutch))+
  #boxplot
  geom_boxplot(
    aes(fill = Clutch, group=Clutch),
    width = boxwidth,
    alpha = boxalpha,
    outlier.shape = NA,
    color = "black",
    linewidth = linesize)+
  #datapoints jittered by quasirandom
  geom_quasirandom(method = 'quasirandom',alpha = jitteralpha*2, width = 0.05, size=3,stroke=0.5,
                   mapping=aes(colour=Pheno))+
  #custom scales
  scale_colour_manual(values = c(albino, D, G, H, I, sibling))+
  scale_fill_manual(values = c(D, G, H, I))+
  scale_shape_manual(values=c(fourshapes))+
  #themes and labels
  ylab("Total increase in length of base of fin (mm)")+
  theme_classic(base_size = 16)+
  theme+
  coord_flip()

ggsave(filename="AP total growth Clutch (pheno) s7.pdf", path="albsib_fin_growth_plots",
       width=5.25, height=5)


##sex
ggplot(AP, mapping = aes(y=totalgrowth, x=Sex, colour=Sex, fill=Sex, 
                         shape = Clutch))+
  #violin
  geom_violin(alpha=0.3, size=0.5, trim=T, scale = scale, mapping=aes(group=Sex))+
  #boxplot
  geom_boxplot(
    aes(fill = Sex, group=Sex),
    width = boxwidth,
    alpha = boxalpha,
    outlier.shape = NA,
    color = "black",
    linewidth = linesize)+
  #datapoints jittered by quasirandom
  geom_quasirandom(method = 'quasirandom',alpha = jitteralpha*2, width = 0.05, size=3,stroke=0.5, 
                   mapping=aes(colour=Pheno))+
  #custom scales
  scale_colour_manual(values = c(albino, female, male, sibling))+
  scale_fill_manual(values = c(female, male))+
  scale_shape_manual(values=c(fourshapes))+
  #theme and labels
  ylab("Total increase in length of base of fin (mm)")+
  theme_classic(base_size = 16)+
  theme+
  coord_flip()

ggsave(filename="AP total growth Sex s7.pdf", path="albsib_fin_growth_plots", 
       width=5.25, height=5)


###PD

#Ray6length

PDstart <- subset(growth, 
                  measurement == "Ray6_length" &
                    Day == "0" )


PDend <- subset(growth, 
                measurement == "Ray6_length" &
                  Day >140 & Day < 160)

PDstart$start <- PDstart$value
PDend$end <- PDend$value

PD <- PDstart %>% left_join(PDend, by = "Individual", suffix = c("", ".y")) %>%
  select(-ends_with(".y"))

PD$totalgrowth <- PD$end - PD$start

##anova total growth PD

aov3.PD.total.growth <- aov(totalgrowth ~ Pheno * Clutch * Sex, data=PD)
anova(aov3.PD.total.growth)
PDtotalgrowth <- anova(aov3.PD.total.growth)
write.csv(PDtotalgrowth, "aov3 PD total growth .csv")

##violin plots total growth PD

##Pheno
ggplot(PD, mapping = aes(y=totalgrowth, x=Pheno, colour=Pheno, fill=Pheno,
                         shape = Clutch))+
  #violin
  geom_violin(alpha=0.3, size=0.5, trim=T, scale = scale, mapping=aes(group=Pheno))+
  #boxplot
  geom_boxplot(
    aes(fill = Pheno, group=Pheno),
    width = boxwidth,
    alpha = boxalpha,
    outlier.shape = NA,
    color = "black",
    linewidth = linesize)+
  #data
  geom_quasirandom(method = 'quasirandom',alpha = jitteralpha*2, width = 0.05, size=3,stroke=0.5,
                   mapping=aes(colour=Sex))+
  #scales
  scale_colour_manual(values = c(albino, female, male, sibling))+
  scale_fill_manual(values = c(albino, sibling))+
  scale_shape_manual(values=c(fourshapes))+
  #theme
  ylab("Total increase in length of fin ray 6 (mm)")+
  theme+
  theme_classic(base_size = 16)+
  theme(legend.position="none",
        axis.title.x = element_blank())

ggsave(filename="PD total growth Pheno s7.pdf", path="albsib_fin_growth_plots", 
       width=5.25, height=5)
       

##Clutch
ggplot(PD, mapping = aes(y=totalgrowth, x=Clutch, colour=Clutch, fill=Clutch,
                         shape = Clutch))+
  #violin
  geom_violin(alpha=0.3, size=0.5, trim=T, scale = scale, mapping=aes(group=Clutch))+
  geom_quasirandom(method = 'quasirandom',alpha = jitteralpha*2, width = 0.05, size=3,stroke=0.5,
                   mapping=aes(colour=Pheno))+
  #boxplot
  geom_boxplot(
    aes(fill = Clutch, group=Clutch),
    width = boxwidth,
    alpha = boxalpha,
    outlier.shape = NA,
    color = "black",
    linewidth = linesize)+
  #custom scale
  scale_colour_manual(values = c(albino, D, G, H, I, sibling))+
  scale_fill_manual(values = c(D, G, H, I))+
  scale_shape_manual(values=c(fourshapes))+
  #label and theme
  ylab("Total increase in length of fin ray 6 (mm)")+
  theme+
  theme_classic(base_size = 16)+
    theme(legend.position="none",
        axis.title.x = element_blank())

ggsave(filename="PD total growth Clutch (Pheno) s7.pdf", path="albsib_fin_growth_plots", 
       width=5.25, height=5)


## Sex
ggplot(PD, mapping = aes(y=totalgrowth, x=Sex, colour=Sex, fill=Sex, 
                         shape = Clutch))+
  #violin
  geom_violin(alpha=0.3, size=0.5, trim=T, scale = scale, mapping=aes(group=Sex))+
  #boxplot
  geom_boxplot(
    aes(fill = Sex, group=Sex),
    width = boxwidth,
    alpha = boxalpha,
    outlier.shape = NA,
    color = "black",
    linewidth = linesize)+
  #data
  geom_quasirandom(method = 'quasirandom',alpha = jitteralpha*2, width = 0.05, size=3,stroke=0.5, 
                   mapping=aes(colour=Pheno))+
  #scale
  scale_colour_manual(values = c(albino, female, male, sibling))+
  scale_fill_manual(values = c(female, male))+
  scale_shape_manual(values=c(fourshapes))+
  #theme and laebls
  ylab("Total increase in length of fin ray 6 (mm)")+
  theme+
  theme_classic(base_size = 16)+
  theme(legend.position="none",
        axis.title.x = element_blank())

ggsave(filename="PD total growth Sex s7.pdf", path="albsib_fin_growth_plots", 
       width=5.25, height=5)


### segment growth data ----------------



seg_lengths <- subset(growth, 
                      measurement == "Ray6_seg1_length" | 
                        measurement == "Ray6_seg2_length" |
                        measurement == "Ray6_seg3_length" | 
                        measurement == "Ray6_seg4_length" |
                        measurement == "Ray6_seg5_length" | 
                        measurement == "Ray6_seg6_length" |
                        measurement == "Ray6_seg7_length" | 
                        measurement == "Ray6_seg8_length" |
                        measurement == "Ray6_seg9_length"|
                        measurement == "Ray6_seg10_length" |
                        measurement == "Ray6_seg11_length" | 
                        measurement == "Ray6_seg12_length" |
                        measurement == "Ray6_seg13_length"|
                        measurement == "Ray6_seg14_length" 
)

###segment growth plots --------------------


seg_lengths_G <- subset(seg_lengths, Clutch == "G")
seg_lengths_G$measurement <- factor(seg_lengths_G$measurement, 
                                    levels=c("Ray6_seg1_length", "Ray6_seg2_length",
                                             "Ray6_seg3_length", "Ray6_seg4_length",
                                             "Ray6_seg5_length", "Ray6_seg6_length",
                                             "Ray6_seg7_length", "Ray6_seg8_length",
                                             "Ray6_seg9_length", "Ray6_seg10_length",
                                             "Ray6_seg11_length", "Ray6_seg12_length",
                                             "Ray6_seg13_length", "Ray6_seg14_length"))


 ggplot(seg_lengths_G, mapping = aes(x=Day, y=value, colour=measurement)) +
  geom_jitter(mapping = aes(x=Day, y=value, colour=measurement), alpha = 0.5, size = 2)+
  geom_smooth(method="lm")+
  scale_colour_manual(values=c(segpal))+
  facet_grid(vars(Pheno))+
  xlim(0,150)+
  theme

ggsave(filename="seg lengths G facet wrap pheno s7.pdf", path="albsib_fin_growth_plots", width=5.88, height=5.34)



ggplot(seg_lengths_G, mapping = aes(x=Day, y=value, colour=measurement)) +
  geom_point(mapping = aes(x=Day, y=value, colour=measurement, shape=Individual), 
             alpha = 0.5)+
  geom_smooth(method="lm")+
  scale_colour_manual(values=c(segpal))+
  facet_grid(rows=vars(measurement), cols=vars(Pheno))+
  scale_shape_manual(values=c(7,8,9,10,11,12,13,14))+
  stat_regline_equation(mapping=aes(colour=measurement))+
  xlim(0,150)+
  theme


ggsave(filename="seg lengths G facet grid pheno s7.pdf", path="albsib_fin_growth_plots", width=5.88, height=16)


##representative individual illustrations

seg_lengths_G$measurement <- factor(seg_lengths_G$measurement, 
                                    levels=c("Ray6_seg12_length",
                                             "Ray6_seg11_length",
                                             "Ray6_seg10_length",
                                             "Ray6_seg9_length",
                                             "Ray6_seg8_length",
                                             "Ray6_seg7_length", 
                                             "Ray6_seg6_length",
                                             "Ray6_seg5_length",
                                             "Ray6_seg4_length",
                                             "Ray6_seg3_length", 
                                             "Ray6_seg2_length",
                                             "Ray6_seg1_length"))

seg_lengths_Ga2 <- subset(seg_lengths_G, 
                          Individual == "Ga2_8") 

seg_lengths_Gs3 <- subset(seg_lengths_G, 
                          Individual == "Gs3_22") 


ggplot(seg_lengths_Ga2, mapping = aes(x=Day, y=value, colour=measurement, fill=measurement)) + 
  geom_bar(position="stack", stat="identity")+
  scale_colour_manual(values=c(segpalrev12))+
  scale_fill_manual(values=c(segpalrev12))+
  ggtitle("Representative albino (M, G2)")+
  ylim(6,0)+
  theme

ggsave(filename="rep albino Ga2 fin growth 12 s7 .pdf", path="albsib_fin_growth_plots", width=5.88, height=3.34)

ggplot(seg_lengths_Gs3, mapping = aes(x=Day, y=value, colour=measurement, fill=measurement)) + 
  geom_bar(position="stack", stat="identity")+
  scale_colour_manual(values=c(segpalrev12))+
  scale_fill_manual(values=c(segpalrev12))+
  ylim(6,0)+
  ggtitle("Representative sibling (M, G3)")+
  theme

ggsave(filename="rep sib Gs3 fin growth 12 s7  .pdf", path="albsib_fin_growth_plots", width=5.88, height=3.34)


### SL 200 data ------------

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


#remove unmatched datapoints
SL <- subset (SL, Day < 170 | Day > 172)

###spot number subsets to isolate final spot numbers in males

SL_spots <- subset (SL,
                    maturity == "Spots")

SL_spotsM <- subset (SL_spots, Sex == "M")

SL_late_spotsM <- subset (SL_spotsM, Day > 180)


### SL day 200 plots--------------- 



ggboxplot(SL_late_spotsM, x = "Clutch", y = "Standard_length_mm", color = "Pheno",
          palette = c(albino, sibling))

ggsave(filename="200_SLM_clutches_sq.pdf", path="sl_plots", units=c("mm"), width=84.1, height=57.6)


ggplot(SL_late_spotsM, mapping = aes(x=Standard_length_mm, y=Pheno,
                                     colour = Pheno, fill =Pheno))+
  #add violin
  geom_violin(alpha=0.3, size=0.6, trim=FALSE, scale = "count")+
  #add data points jittered by quasirandom
  geom_quasirandom(size=2, alpha=jitteralpha, mapping = aes ( shape = Clutch), width=0.15)+
  #boxplot
  geom_boxplot(
    aes(fill = Pheno, group=Pheno),
    width = boxwidth,
    alpha = boxalpha,
    outlier.shape = NA,
    color = "black",
    linewidth = linesize)+
  #custom scales
  scale_colour_manual(values = c(albino, sibling))+
  scale_fill_manual(values = c(albino, sibling))+
  scale_shape_manual(values=fourshapes)+
  #theme and labels
  theme_classic(base_size = 16)+ 
  theme+
  coord_flip()

ggsave(filename="200_SLM_sqsq.pdf", path="sl_plots", units=c("mm"), width=88.8, height=55.6)


#sl
ggdensity(SL_late_spotsM$Standard_length_mm)
ggqqplot(SL_late_spotsM$Standard_length_mm)
shapiro.test(SL_late_spotsM$Standard_length_mm)
#normal
res.aov2.SLM <- aov(Standard_length_mm ~ Pheno * Clutch, data = SL_late_spotsM)
aov2_SLM <- anova(res.aov2.SLM)
write.csv(aov2_SLM, "aov2_SL_day200_M .csv")
anova(res.aov2.SLM)



### aggregation area data ------------


#import size measurements

sizes <- read_csv("mutant clutches size and spacing.csv")



### size calculations


#create columns for areas that include each whole aggregation 
sizes$Area1 <- (sizes$ag1length)*(sizes$ag1calcwidth)
sizes$Area2 <- (sizes$ag2length)*(sizes$ag2calcwidth)
sizes$Area3 <- (sizes$ag3length)*(sizes$ag3calcwidth)

#convert NAs to zeros so they can be added for total area
sizes <- sizes %>% mutate(Area1 = ifelse(is.na(Area1), 0, Area1))
sizes <- sizes %>% mutate(Area2 = ifelse(is.na(Area2), 0, Area2))
sizes <- sizes %>% mutate(Area3 = ifelse(is.na(Area3), 0, Area3))

#add for total area within which all aggregations included
sizes$Areatotal <- (sizes$Area1)+(sizes$Area2)+(sizes$Area3)

#divide by number of aggregations in each individual
sizes$Areatotal_div <- sizes$Areatotal/sizes$Ags

#subset by days and sex
sizes_mid <- sizes %>% subset(Day < 140)
View(sizes_mid)


### aggregation area plots -----------

##AVERAGE AREA VS PHENO

#average area vs pheno
ggplot(sizes_mid, mapping = aes(x=Areatotal_div, y=Pheno,
                                colour = Pheno, fill =Pheno,
                                shape=Clutch))+
  #violin
  geom_violin(alpha=0.3, size=0.5, trim=T, scale = "area", mapping=aes(group=Pheno))+
  #boxplot
  geom_boxplot(
    aes(fill = Pheno, group=Pheno),
    width = boxwidth,
    alpha = boxalpha,
    outlier.shape = NA,
    color = "black",
    linewidth = linesize)+
  #datapoints jittered by quasirandom
  geom_quasirandom(method = 'quasirandom',alpha = jitteralpha, width = 0.1,size=3,stroke=0.5)+
  #scales
  scale_colour_manual(values = c(albino, sibling))+
  scale_fill_manual(values = c(albino, sibling))+
  scale_shape_manual(values=c(16, 17, 15, 18))+
  #labels and theme
  xlab("Average aggregation area (total area/number of ags)")+
  geom_signif(comparisons = list(c("albino", "sibling")),   
              map_signif_level=TRUE, colour="black", y_position=2.7,
              tip_length=0)+
  theme_classic(base_size = 16)+ 
  theme

ggsave(filename="average area s7.pdf", path="sizes",
       width=5.25, height=5)

aov3.avarea <- aov(Areatotal_div ~ Pheno * Clutch * Sex, data=sizes_mid)
anova(aov3.avarea)

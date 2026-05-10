library(ggplot2)
library(dplyr)
library(ggpmisc)

dat <- read.csv("Masoko_meta_master_20230331.csv", header = FALSE)
dat <- dat[-1, ]
colnames(dat) <- dat[1, ]
dat <- dat[-1, ]
names(dat) <- make.names(names(dat), unique = TRUE)

dat <- dat %>%
  mutate(
    SL                 = as.numeric(SL),
    EggspotNumberTheis = as.numeric(EggspotNumberTheis),
    spot.hue           = as.numeric(spot.hue),
    black.L            = as.numeric(black.L),
    Group = case_when(
      PCA_Group == "High"   ~ "Benthic",
      PCA_Group == "Low"    ~ "Littoral",
      PCA_Group == "Middle" ~ "Intermediate"
    )
  )

group_cols <- c("Benthic" = "steelblue2", "Littoral" = "goldenrod2", "Intermediate" = "grey")
facet_labs <- c("High" = "Benthic", "Low" = "Littoral", "Middle" = "Intermediate")

d_num <- dat %>% filter(!is.na(SL), !is.na(EggspotNumberTheis), !is.na(Group))
d_hue <- dat %>% filter(!is.na(SL), !is.na(spot.hue),           !is.na(Group))
d_blk <- dat %>% filter(!is.na(SL), !is.na(black.L),            !is.na(Group))

# ── Shared layers: single black regression line + stats ─────────────────────
stat_layers <- list(
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.9),
  stat_poly_eq(aes(label = ..eq.label..),       formula = y ~ x,
               eq.with.lhs = "italic(y)~`=`~", parse = TRUE,
               size = 3.5, color = "black", label.x = 0.05, label.y = 0.98),
  stat_poly_eq(aes(label = ..rr.label..),        formula = y ~ x,
               parse = TRUE, size = 3.5, color = "black",
               label.x = 0.05, label.y = 0.88),
  stat_poly_eq(aes(label = ..p.value.label..),   formula = y ~ x,
               parse = TRUE, size = 3.5, color = "black",
               label.x = 0.05, label.y = 0.78)
)

base_theme <- list(
  theme_classic(base_size = 14),
  theme(
    strip.text      = element_text(size = 13, face = "bold"),
    strip.background = element_rect(color = "grey80", fill = "grey95"),
    axis.title      = element_text(face = "bold"),
    legend.position = "bottom",
    legend.title    = element_text(size = 13),
    legend.text     = element_text(size = 12)
  )
)

# ── FACETED ─────────────────────────────────────────────────────────────────
facet_layer <- facet_wrap(~ PCA_Group, labeller = labeller(PCA_Group = facet_labs))

pA_facet <- ggplot(d_num, aes(x = SL, y = EggspotNumberTheis)) +
  geom_point(aes(color = Group), size = 2.5, alpha = 0.8) +
  stat_layers + facet_layer + base_theme +
  scale_color_manual(values = group_cols) +
  labs(x = "Standard Length (mm)", y = "Eggspot Number", color = "Group", tag = "A")

pB_facet <- ggplot(d_hue, aes(x = SL, y = spot.hue)) +
  geom_point(aes(color = Group), size = 2.5, alpha = 0.8) +
  stat_layers + facet_layer + base_theme +
  scale_color_manual(values = group_cols) +
  labs(x = "Standard Length (mm)", y = "Eggspot Hue", color = "Group", tag = "B")

pC_facet <- ggplot(d_blk, aes(x = SL, y = black.L)) +
  geom_point(aes(color = Group), size = 2.5, alpha = 0.8) +
  stat_layers + facet_layer + base_theme +
  scale_color_manual(values = group_cols) +
  labs(x = "Standard Length (mm)", y = "Black Anal Fin Luminance (L*)", color = "Group", tag = "C")

# ── COMBINED (all groups, single regression line) ────────────────────────────
pA_combined <- ggplot(d_num, aes(x = SL, y = EggspotNumberTheis)) +
  geom_point(aes(color = Group), size = 2.5, alpha = 0.8) +
  stat_layers + base_theme +
  scale_color_manual(values = group_cols) +
  labs(x = "Standard Length (mm)", y = "Eggspot Number", color = "Group", tag = "A")

pB_combined <- ggplot(d_hue, aes(x = SL, y = spot.hue)) +
  geom_point(aes(color = Group), size = 2.5, alpha = 0.8) +
  stat_layers + base_theme +
  scale_color_manual(values = group_cols) +
  labs(x = "Standard Length (mm)", y = "Eggspot Hue", color = "Group", tag = "B")

pC_combined <- ggplot(d_blk, aes(x = SL, y = black.L)) +
  geom_point(aes(color = Group), size = 2.5, alpha = 0.8) +
  stat_layers + base_theme +
  scale_color_manual(values = group_cols) +
  labs(x = "Standard Length (mm)", y = "Black Anal Fin Luminance (L*)", color = "Group", tag = "C")

# ── Print ────────────────────────────────────────────────────────────────────
pA_facet;    pB_facet;    pC_facet
pA_combined; pB_combined; pC_combined




d_egg_blk <- dat %>% filter(!is.na(EggspotNumberTheis), !is.na(black.L), !is.na(Group))

pD_combined <- ggplot(d_egg_blk, aes(x = EggspotNumberTheis, y = black.L)) +
  geom_point(aes(color = Group), size = 2.5, alpha = 0.8) +
  stat_layers + base_theme +
  scale_color_manual(values = group_cols) +
  labs(x = "Eggspot Number", y = "Black Anal Fin Luminance (L*)", color = "Group", tag = "D")

pD_facet <- ggplot(d_egg_blk, aes(x = EggspotNumberTheis, y = black.L)) +
  geom_point(aes(color = Group), size = 2.5, alpha = 0.8) +
  stat_layers + facet_layer + base_theme +
  scale_color_manual(values = group_cols) +
  labs(x = "Eggspot Number", y = "Black Anal Fin Luminance (L*)", color = "Group", tag = "D")

pD_combined
pD_facet






#####other comb
# Hue vs Number
d_hue_num <- dat %>% filter(!is.na(spot.hue), !is.na(EggspotNumberTheis), !is.na(Group))

pE_combined <- ggplot(d_hue_num, aes(x = EggspotNumberTheis, y = spot.hue)) +
  geom_point(aes(color = Group), size = 2.5, alpha = 0.8) +
  stat_layers + base_theme +
  scale_color_manual(values = group_cols) +
  labs(x = "Eggspot Number", y = "Eggspot Hue", color = "Group", tag = "E")

pE_facet <- ggplot(d_hue_num, aes(x = EggspotNumberTheis, y = spot.hue)) +
  geom_point(aes(color = Group), size = 2.5, alpha = 0.8) +
  stat_layers + facet_layer + base_theme +
  scale_color_manual(values = group_cols) +
  labs(x = "Eggspot Number", y = "Eggspot Hue", color = "Group", tag = "E")

# Hue vs black.L
d_hue_blk <- dat %>% filter(!is.na(spot.hue), !is.na(black.L), !is.na(Group))

pF_combined <- ggplot(d_hue_blk, aes(x = black.L, y = spot.hue)) +
  geom_point(aes(color = Group), size = 2.5, alpha = 0.8) +
  stat_layers + base_theme +
  scale_color_manual(values = group_cols) +
  labs(x = "Black Anal Fin Luminance (L*)", y = "Eggspot Hue", color = "Group", tag = "F")

pF_facet <- ggplot(d_hue_blk, aes(x = black.L, y = spot.hue)) +
  geom_point(aes(color = Group), size = 2.5, alpha = 0.8) +
  stat_layers + facet_layer + base_theme +
  scale_color_manual(values = group_cols) +
  labs(x = "Black Anal Fin Luminance (L*)", y = "Eggspot Hue", color = "Group", tag = "F")

# Print
pE_combined; pE_facet
pF_combined; pF_facet

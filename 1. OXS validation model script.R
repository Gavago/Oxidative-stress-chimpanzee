library(tidyverse)
library(magrittr)
library(lmerTest) # glmer, lmer
library(ggcorrplot) # cor_pmat

# i. load functions----
source("functions - biomarks stats and dfs for regression.R")
source("functions - corr & glmer info extract.R")
z. <- function(x) scale(x)
select <- dplyr::select

# ii.load data -----
load("oxs datasets for pub.Rdata", verbose = T)

# 1. Biomarker correlations ---------

#time adjusted biomarkers
cmta <- oxs_final_all_pub %>%
  select(OHDG_SG, Iso_SG, TBARS_SG_TAR, Neo_SG_TAR, TAC_SG_TAR) %>%
  rename(`8-OHdG` = OHDG_SG, Isop = Iso_SG, `MDA-TBARS` = TBARS_SG_TAR,
         Neo = Neo_SG_TAR, TAC = TAC_SG_TAR) %>%
  cor(., method = "spearman" , use = "pairwise.complete.obs")


cmta_pvalue <- oxs_final_all_pub %>%
  select(OHDG_SG, Iso_SG, TBARS_SG_TAR, Neo_SG_TAR, TAC_SG_TAR) %>%
  rename(`8-OHdG` = OHDG_SG, Isop = Iso_SG, `MDA-TBARS` = TBARS_SG_TAR,
         Neo = Neo_SG_TAR, TAC = TAC_SG_TAR) %>%
  cor_pmat(., method = "spearman" , use = "pairwise.complete.obs")

cmta
cmta_pvalue


# 2. Before during and after epidemic ----


ohdge <- ohdg_bdae_df %>%
  glmer(OHDG_SG ~ C(bdae_timing_ext, base = 1) + sex + z.(Age) + (1|ChimpID), family = Gamma(link = "log"), data = .)
summary(ohdge)

isoe <- iso_bdae_df %>%
  glmer(Iso_SG ~ C(bdae_timing_ext, base = 1) + sex + z.(Age) + (1|ChimpID), family = Gamma(link = "log"), data = .)
summary(isoe)

tbarse <- tbars_bdae_df_tar %>%
  lmer(TBARS_SG_TAR ~ C(bdae_timing_ext, base = 1) + sex + z.(Age) + (1|ChimpID), data = .)
summary(tbarse)

neoe <- neo_bdae_df_tar %>%
  glmer(neo_pos ~ C(bdae_timing_ext, base = 1) + sex + z.(Age) + (1|ChimpID), family = Gamma(link = "log"), data = .)
summary(neoe)

tace <- tac_bdae_df_tar %>%
  lmer(TAC_SG_TAR ~ C(bdae_timing_ext, base = 1) + sex + z.(Age) + (1|ChimpID), data = .)
summary(tace)

# 3. Before during after injury ------
ohdgi <- oxs_final_avg_pub %>%
  glmer(OHDG_SG ~ C(bdai_timing_ext, base = 1) + sex + z.(Age) + (1|ChimpID), family = Gamma(link = "log"), data = .)
summary(ohdgi)

isoi <-  oxs_final_avg_pub %>%
  glmer(Iso_SG ~ C(bdai_timing_ext, base = 1) + sex + z.(Age) + (1|ChimpID), family = Gamma(link = "log"), data = .)
summary(isoi)

tbarsi <- oxs_final_avg_pub %>%  
  glmer(tbars_pos ~ C(bdai_timing_ext, base = 1) + sex + z.(Age) + (1|ChimpID), family = Gamma(link = "log"), data = .)
summary(tbarsi)

neoi <- oxs_final_avg_pub %>% 
  glmer(neo_pos ~ C(bdai_timing_ext, base = 1) + sex + z.(Age) + (1|ChimpID), family = Gamma(link = "log"), data = .)
summary(neoi)

taci <- oxs_final_avg_pub %>%  
  glmer(tac_pos ~ C(bdai_timing_ext, base = 1) + sex + z.(Age) + (1|ChimpID), family = Gamma(link = "log"), data = .)
summary(taci)

# 4. Mixed-longitudinal -----

ohdg_x <- cross_pub %>%
  glmer(OHDG_SG ~ sex + z.(Age) + (1|ChimpID), family = Gamma(link = "log"), data = .) #sig neg age
ohdg_x_i <- cross_pub %>%
  glmer(OHDG_SG ~ sex + z.(Age) + sex*z.(Age) + (1|ChimpID), family = Gamma(link = "log"), data = .) #nada

summary(ohdg_x)
summary(ohdg_x_i)


iso_x <- cross_pub %>%
  glmer(Iso_SG ~ sex + z.(Age) + (1|ChimpID), family = Gamma(link = "log"), data = .)
iso_x_i <- cross_pub %>%
  glmer(Iso_SG ~ sex + z.(Age) + sex*z.(Age) + (1|ChimpID), family = Gamma(link = "log"), data = .)

summary(iso_x)
summary(iso_x_i)

tbars_x <- cross_pub %>%
  glmer(tbars_pos ~ sex + z.(Age) + (1|ChimpID), family = Gamma(link = "log"), data = .)
tbars_x_i <- cross_pub %>%
  glmer(tbars_pos ~ sex + z.(Age) + sex*z.(Age) + (1|ChimpID), family = Gamma(link = "log"), data = .)

summary(tbars_x)
summary(tbars_x_i)

neo_x <- cross_pub %>%
  glmer(neo_pos ~ sex + z.(Age) + (1|ChimpID), family = Gamma(link = "log"), data = .)
neo_x_i <- cross_pub %>%
  glmer(neo_pos ~ sex + z.(Age) + sex*z.(Age) + (1|ChimpID), family = Gamma(link = "log"), data = .)

summary(neo_x)
summary(neo_x_i)

tac_x <- cross_pub %>%
  lmer(TAC_SG_TAR ~ sex + z.(Age) + (1|ChimpID), data = .)
tac_x_i <- cross_pub %>%
  lmer(TAC_SG_TAR ~ sex + z.(Age) + sex*z.(Age) + (1|ChimpID), data = .)

summary(tac_x)
summary(tac_x_i)


# 4i. Mixed-longitudinal past-prime ----

ohdg_x_f <- cross_frail_pub %>%
  glmer(OHDG_SG ~ sex + z.(Age) + (1|ChimpID), family = Gamma(link = "log"), data = .) #sig neg age
ohdg_x_f_i <- cross_frail_pub %>%
  glmer(OHDG_SG ~ sex + z.(Age) + sex*z.(Age) + (1|ChimpID), family = Gamma(link = "log"), data = .) #nada

summary(ohdg_x_f)
summary(ohdg_x_f_i)

iso_x_f <- cross_frail_pub %>%
  glmer(Iso_SG ~ sex + z.(Age) + (1|ChimpID), family = Gamma(link = "log"), data = .)
iso_x_f_i <- cross_frail_pub %>%
  glmer(Iso_SG ~ sex + z.(Age) + sex*z.(Age) + (1|ChimpID), family = Gamma(link = "log"), data = .)

summary(iso_x_f)
summary(iso_x_f_i)

tbars_x_f <- cross_frail_pub %>%
  glmer(tbars_pos ~ sex + z.(Age) + (1|ChimpID), family = Gamma(link = "log"), data = .)
tbars_x_f_i <- cross_frail_pub %>%
  glmer(tbars_pos ~ sex + z.(Age) + sex*z.(Age) + (1|ChimpID), family = Gamma(link = "log"), data = .)

summary(tbars_x_f)
summary(tbars_x_f_i)

neo_x_f <- cross_frail_pub %>%
  glmer(neo_pos ~ sex + z.(Age) + (1|ChimpID), family = Gamma(link = "log"), data = .)
neo_x_f_i <- cross_frail_pub %>%
  glmer(neo_pos ~ sex + z.(Age) + sex*z.(Age) + (1|ChimpID), family = Gamma(link = "log"), data = .)


summary(neo_x_f)
summary(neo_x_f_i)


tac_x_f <- cross_frail_pub %>%
  lmer(TAC_SG_TAR ~ sex + z.(Age) + (1|ChimpID), data = .)
tac_x_f_i <- cross_frail_pub %>%
  lmer(TAC_SG_TAR ~ sex + z.(Age) + sex*z.(Age) + (1|ChimpID), data = .)


summary(tac_x_f)
summary(tac_x_f_i)

# 5. Lead up to death -----

lutd_ohdg <- oxs_final_avg_pub %>% filter(!is.na(lutd)) %>%
  glmer(OHDG_SG ~  z.(lutd_years_until) + sex + (1|ChimpID), data = ., family = Gamma(link = "log"))
summary(lutd_ohdg)

lutd_iso <- oxs_final_avg_pub %>% filter(!is.na(lutd)) %>%
  glmer(Iso_SG ~ z.(lutd_years_until) + sex + (1|ChimpID), family = Gamma(link = "log"), data = .)
summary(lutd_iso)

lutd_tbars <- oxs_final_avg_pub %>% filter(!is.na(lutd)) %>%
  glmer(tbars_pos ~ z.(lutd_years_until) + sex + (1|ChimpID), data = ., family = Gamma(link = "log"))
summary(lutd_tbars)

lutd_neo <- oxs_final_avg_pub %>% filter(!is.na(lutd)) %>%
  lmer(Neo_SG_TAR ~  z.(lutd_years_until) + sex + (1|ChimpID), data = .)
summary(lutd_neo)

lutd_tac <- oxs_final_avg_pub %>% filter(!is.na(lutd)) %>%
  lmer(TAC_SG_TAR~ z.(lutd_years_until) + sex + (1|ChimpID), data = .) 
summary(lutd_tac)

# 5i.Lead up to death: no young individual (KK) -----

lutd_ohdgnk <- oxs_final_avg_pub %>% filter(!is.na(lutd) & ChimpID != "KK") %>%
  glmer(OHDG_SG ~  z.(lutd_years_until) + sex + (1|ChimpID), data = ., family = Gamma(link = "log"))
summary(lutd_ohdgnk)

lutd_isonk <- oxs_final_avg_pub %>% filter(!is.na(lutd) & ChimpID != "KK") %>%
  glmer(Iso_SG ~ z.(lutd_years_until) + sex + (1|ChimpID), family = Gamma(link = "log"), data = .)
summary(lutd_isonk)

lutd_tbarsnk <- oxs_final_avg_pub %>% filter(!is.na(lutd) & ChimpID != "KK") %>%
  glmer(tbars_pos ~ z.(lutd_years_until) + sex + (1|ChimpID), data = ., family = Gamma(link = "log"))
summary(lutd_tbarsnk)

lutd_neonk <- oxs_final_avg_pub %>% filter(!is.na(lutd) & ChimpID != "KK") %>%
  lmer(Neo_SG_TAR ~  z.(lutd_years_until) + sex + (1|ChimpID), data = .)
summary(lutd_neonk)

lutd_tacnk <- oxs_final_avg_pub %>% filter(!is.na(lutd) & ChimpID != "KK") %>%
  lmer(TAC_SG_TAR~ z.(lutd_years_until) +  sex + (1|ChimpID), data = .) 
summary(lutd_tacnk)


# 6. Rerun age-related analysis of MDA-TBARS: -----
# checking possible storage time effect using samples stored <= 3 years (m-l) or within 3 of death (lutd)

# mixed longitudinal (m-l)
load("df of samples w less than 3 yrs storage time per biomarker for x-sec test.Rdata", verbose = T)

cross_recent <- recent %>%
  filter(ChimpID != "BT" & (!bdae_timing_ext %in% c("dur","aft") | is.na(bdae_timing_ext)) &
           (!bdai_timing_ext %in% c("dur","aft") | is.na(bdai_timing_ext)))

tbars_xr <- cross_recent %>%
  glmer(tbars_pos ~ sex + z.(Age) + (1|ChimpID), family = Gamma(link = "log"), data = .)
summary(tbars_xr)

load("df samples collected within 3 years of individual death.Rdata", verbose = T)

# lead up to death (lutd)
lutd_tbars3 <- lutd_3yrs %>% filter(!is.na(lutd)) %>%
  glmer(tbars_pos ~  z.(lutd_years_until) + sex + (1|ChimpID), data = ., family = Gamma(link = "log"))
summary(lutd_tbars3)

lutd_tbars3k <- lutd_3yrs %>% filter(!is.na(lutd) & ChimpID != "KK") %>%
  glmer(tbars_pos ~  z.(lutd_years_until) + sex + (1|ChimpID), data = ., family = Gamma(link = "log"))
summary(lutd_tbars3k)

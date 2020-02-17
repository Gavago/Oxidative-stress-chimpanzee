library(tidyverse)
library(ggcorrplot)

#colorblind palette
cbPalette <- c( "#CC79A7","#56B4E9","#009E73", "#0072B2", "#F0E442", "#D55E00", "#999999","#E69F00")



# i.  load functions ----
fac_bda <- function(data) {
  fac_data <- data %>%
    mutate(bdae_timing_ext = case_when( #replace bdae and bdai ext levels
      bdae_timing_ext == "bef" ~ "before",
      bdae_timing_ext == "dur" ~ "during",
      bdae_timing_ext == "aft" ~ "after",
    ),
    bdai_timing_ext = case_when(
      bdai_timing_ext == "bef" ~ "before",
      bdai_timing_ext == "dur" ~ "during",
      bdai_timing_ext == "aft" ~ "after",
    )) %>%
    mutate(bdae_timing_ext = factor(bdae_timing_ext, levels = c("before","during","after")), #order bda levels
           bdai_timing_ext = factor(bdai_timing_ext, levels = c("before","during","after"))) 
  return(fac_data)
}

# ii. load data ----
load("oxs datasets for pub.Rdata", verbose = T)

# 1.  Biomarker correlations ----
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

# For homemade heat map
rho <- cmta %>%
  data.frame() %>%
  rownames_to_column() %>%
  rename(var1 =  rowname, `8-OHdG` = X8.OHdG, `MDA-TBARS` = MDA.TBARS) %>%
  gather(var2, Rho, -var1)

rho_p_val <- cmta_pvalue %>%
  data.frame() %>%
  rownames_to_column() %>%
  rename(var1 =  rowname, `8-OHdG` = X8.OHdG, `MDA-TBARS` = MDA.TBARS) %>%
  gather(var2, p_val, -var1)

long_corr <- merge(rho, rho_p_val, by = c("var1", "var2")) %>%
  distinct(Rho, p_val, .keep_all = T) %>%
  filter(var1 != var2) %>%
  mutate(Rho = round(Rho, 2)) %>%
  mutate(sig = case_when(
    p_val < 0.0001 ~ "****",
    p_val < 0.001 ~ "***",
    p_val < 0.01 ~ "**",
    p_val < 0.05 ~ "*",
    p_val < 0.1 ~ "†",
    TRUE ~ ""
  )) %>%
  unite("Rho_sig", Rho, sig, sep = " ", remove = F)

ggplot(long_corr, aes(var1, var2, fill = Rho)) +
  geom_tile() +
  scale_fill_gradient2(low="blue", mid="white", high="red") +
  geom_text(aes(label = Rho_sig)) +
  labs(x = "", y = "") +
  coord_flip() +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12)) 

# 2.  Before during and after epidemic (BDAE) -----
trail_o_the_dead <- c("LR", "KK", "ST", "PG") 

# 2a. BDAE OHDG -----
ohdg_bdae_dead <- ohdg_bdae_df %>%
  filter(ChimpID %in% trail_o_the_dead) %>%
  fac_bda()

ohdg_bdae_vis <- ohdg_bdae_df %>%
  filter(!ChimpID %in% trail_o_the_dead) %>% 
  fac_bda() %>%
  ggplot(., aes(x = bdae_timing_ext, y = OHDG_SG)) + 
  geom_boxplot(outlier.shape = NA, color = "darkblue") + 
  geom_point(position = "jitter", shape = 1, color = "darkblue") +
  geom_jitter(data = ohdg_bdae_dead, aes(x = bdae_timing_ext, y = OHDG_SG), color = "black", shape = 18, size = 3.5) +
  ggtitle("8-OHdG") + 
  labs(x = "period of epidemic", y = "8-OHdG ng/mL") +
  geom_segment(aes(x = "before", y = 60, xend = "during", yend = 60), linetype = 2 )+
  geom_segment(aes(x = "before", y = 70, xend = "after", yend = 70), linetype = 2) +
  annotate("text", x = 0.6, y = 80, label = "A)", size = 10) +
  annotate("text", x = 1.5, y = 61, label = "*", size = 12) +
  annotate("text", x = 2, y = 71, label = "*", size = 12) +
  theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5, size = 18)) 

ohdg_bdae_vis

# 2b. BDAE Isoprostanes ------

iso_bdae_vis <- iso_bdae_df %>%
  fac_bda() %>%
  ggplot(., aes(x = bdae_timing_ext, y = Iso_SG)) + 
  geom_boxplot(outlier.shape = NA, color = "forestgreen") + 
  geom_point(position = "jitter", shape = 1, color = "forestgreen") +
  ggtitle("Isoprostanes") +
  labs(x = "period of epidemic", y = "Isoprostanes ng/mL") +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5, size = 18)) +
  annotate("text", x = 0.6, y = 14, label = "B)", size = 10)

iso_bdae_vis

# 2c. BDAE TBARS -------

tbars_bdae_dead <- tbars_bdae_df_tar %>%
  filter(ChimpID %in% trail_o_the_dead) %>%
  fac_bda()

tbars_bdae_vis <- tbars_bdae_df_tar %>%
  fac_bda() %>%
  filter(!ChimpID %in% trail_o_the_dead) %>% 
  ggplot(., aes(x = bdae_timing_ext, y = TBARS_SG_TAR)) + 
  geom_boxplot(outlier.shape = NA, color = "firebrick3") + 
  geom_point(position = "jitter", shape = 1, color = "firebrick3") +
  geom_jitter(data = tbars_bdae_dead, aes(x = bdae_timing_ext, y = TBARS_SG_TAR), color = "black", shape = 18, size = 3.5) +
  ggtitle("MDA-TBARS") + 
  labs(x = "period of epidemic", y = "MDA-TBARS µM (time adjusted residual)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 18)) +
  geom_segment(aes(x = "during", y = 20, xend = "after", yend = 20), linetype = 2)  +
  annotate("text", x = 0.6, y = 35, label = "C)", size = 10) +
  annotate("text", x = 2.5, y = 20.5, label = "*", size = 12) 

tbars_bdae_vis

# 2d. BDAE Neopterin -----
neo_bdae_dead <- neo_bdae_df_tar %>%
  filter(ChimpID %in% trail_o_the_dead) %>%
  fac_bda()

neo_bdae_vis <- neo_bdae_df_tar %>%
  fac_bda() %>%
  filter(!ChimpID %in% trail_o_the_dead) %>% 
  ggplot(., aes(x = bdae_timing_ext, y = Neo_SG_TAR)) + 
  geom_boxplot(outlier.shape = NA, color = "darkorchid4") + 
  geom_point(position = "jitter", shape = 1, color = "darkorchid4") +
  geom_jitter(data = neo_bdae_dead, aes(x = bdae_timing_ext, y = Neo_SG_TAR), color = "black", shape = 18, size = 3.5) +
  ggtitle("Neopterin") + 
  labs(x = "period of epidemic", y = "Neopterin nM (time adjusted residual)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 18)) +
  geom_segment(aes(x = "before", y = 5800, xend = "after", yend = 5800), linetype = 2 ) +
  geom_segment(aes(x = "during", y = 4800, xend = "after", yend = 4800), linetype = 2 ) +
  annotate("text", x = 0.6, y = 7000, label = "D)", size = 10) +
  annotate("text", x = 2, y = 5900, label = "*", size = 12) +
  annotate("text", x = 2.5, y = 4900, label = "*", size = 12)
neo_bdae_vis

# 2e. BDAE TAC -----
tac_bdae_vis <- tac_bdae_df_tar %>%
  fac_bda() %>%
  ggplot(., aes(x = bdae_timing_ext, y = TAC_SG_TAR)) + 
  geom_boxplot(outlier.shape = NA, color = "deeppink3") + 
  geom_point(position = "jitter", shape = 1, color = "deeppink3") +
  ggtitle("TAC") + 
  labs(x = "period of epidemic", y = "TAC mM (time adjusted residual)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 18)) +
  annotate("text", x = 0.6, y = 18, label = "E)", size = 10)

tac_bdae_vis

# 3.  Before during and afer injury (BDAI) ----

# 3a. BDAI OHDG----

average_ohdg <- oxs_final_avg_pub %>%
  fac_bda() %>%
  filter(!is.na(bdai_timing_ext)) %>%
  group_by(ChimpID, bdai_timing_ext) %>%
  summarise(avg_ohdg = mean(OHDG_SG, na.rm = T))

bdai_ohdg <- oxs_final_avg_pub %>%
  fac_bda() %>%
  filter(!is.na(bdai_timing_ext)) %>%
  ggplot() + 
  geom_line(data = average_ohdg, aes(x = bdai_timing_ext, y = avg_ohdg, group = ChimpID, color = ChimpID)) +
  geom_jitter(aes(x = bdai_timing_ext, y = OHDG_SG, color = ChimpID), width = .1) +
  ggtitle("8-OHdG") +
  labs(x = "period of injury", y = "8-OHdG ng/mL") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 18), legend.position = "none") +
  annotate("text", x = 0.6, y = 92, label = "A)", size = 10)
bdai_ohdg

# 3b. BDAI Isoprostanes ----

average_iso <- oxs_final_avg_pub %>%
  fac_bda() %>%
  filter(!is.na(bdai_timing_ext)) %>%
  group_by(ChimpID, bdai_timing_ext) %>%
  summarise(avg_iso = mean(Iso_SG, na.rm = T))

bdai_iso <- oxs_final_avg_pub %>%
  fac_bda() %>%
  filter(!is.na(bdai_timing_ext)) %>%
  ggplot() + 
  geom_line(data = average_iso, aes(x = bdai_timing_ext, y = avg_iso, group = ChimpID, color = ChimpID)) + 
  geom_jitter(aes(x = bdai_timing_ext, y = Iso_SG, color = ChimpID), width = .1) +
  ggtitle("Isoprostanes") +
  labs(x = "period of injury", y = "F-isoprostane ng/mL") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 18), legend.position = "none") +
  geom_segment(aes(x = "before", y = 7.8, xend = "during", yend = 7.8), linetype = 2, color = "black") +
  annotate("text", x = 0.6, y = 11.5, label = "B)", size = 10) +
  annotate("text", x = 1.55, y = 8, label = "*", size = 12)
bdai_iso

# 3c. BDAI TBARS ----

average_tbars <- oxs_final_avg_pub %>%
  fac_bda() %>%
  filter(!is.na(bdai_timing_ext)) %>%
  group_by(ChimpID, bdai_timing_ext) %>%
  summarise(avg_tbars = mean(TBARS_SG_TAR, na.rm = T))

bdai_tbars <- oxs_final_avg_pub %>%
  fac_bda() %>%
  filter(!is.na(bdai_timing_ext)) %>%
  ggplot() + 
  geom_line(data = average_tbars, aes(x = bdai_timing_ext, y = avg_tbars, group = ChimpID, color = ChimpID)) + 
  geom_jitter(aes(x = bdai_timing_ext, y = TBARS_SG_TAR, group = ChimpID, color = ChimpID), width = .1) +
  ggtitle("MDA-TBARS") +
  labs(x = "period of injury", y = "MDA-TBARS ?M (time adjusted residual)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 18), legend.position = "none") +
  annotate("text", x = 0.6, y = 25, label = "C)", size = 10)

bdai_tbars

# 3d. BDAI Neopterin  -----
average_neo <- oxs_final_avg_pub %>%
  filter(!(ChimpID == "PB" & Date == lubridate::ymd("2012-11-21"))) %>%
  fac_bda() %>%
  filter(!is.na(bdai_timing_ext)) %>%
  group_by(ChimpID, bdai_timing_ext) %>%
  summarise(avg_neo = mean(Neo_SG_TAR, na.rm = T))

bdai_neo <- oxs_final_avg_pub %>%
  filter(!(ChimpID == "PB" & Date == lubridate::ymd("2012-11-21"))) %>%
  fac_bda() %>%
  filter(!is.na(bdai_timing_ext)) %>%
  ggplot() + 
  geom_line(data =average_neo, aes(x = bdai_timing_ext, y = avg_neo, group = ChimpID, color = ChimpID)) + 
  geom_jitter(aes(x = bdai_timing_ext, y = Neo_SG_TAR, color = ChimpID), width = .1) +
  ggtitle("Neopterin") + 
  labs(color = "Individual" , x = "period of injury", y = "Neopterin nmol/L (time adjusted residual)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 18), legend.position = "none") +
  geom_segment(aes(x = "before", y = 1900, xend = "during", yend = 1900), linetype = 2, color = "black") +
  geom_segment(aes(x = "during", y = 2000, xend = "after", yend = 2000), linetype = 2, color = "black") +
  annotate("text", x = 0.6, y = 3500, label = "D)", size = 10) +
  annotate("text", x = 1.5, y = 2000, label = "*", size = 12) +
  annotate("text", x = 2.5, y = 2200, label = "*", size = 12)

bdai_neo

# 3e. BDAI TAC -----
average_tac <- oxs_final_avg_pub %>%
  fac_bda() %>%
  filter(!is.na(bdai_timing_ext)) %>%
  group_by(ChimpID, bdai_timing_ext) %>%
  summarise(avg_tac = mean(TAC_SG_TAR, na.rm = T))


bdai_tac <- oxs_final_avg_pub %>%
  fac_bda() %>%
  filter(!is.na(bdai_timing_ext)) %>%
  ggplot() + 
  geom_line(data = average_tac,  aes(x = bdai_timing_ext, y = avg_tac, group = ChimpID, color = ChimpID)) + 
  geom_jitter( aes(x = bdai_timing_ext, y = TAC_SG_TAR, color = ChimpID), width = .1) +
  ggtitle("TAC") +
  labs(x = "period of injury", y = "TAC mM (time adjusted residual)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 18), legend.position = "bottom", legend.text = element_text(size = 12),
        legend.title = element_text(size = 12)) +
  annotate("text", x = 0.6, y = 10, label = "E)", size = 10)

bdai_tac

# 4.  Mixed longitudinal -----

# 4a. m-l OHDG ----
cross_pub %>%
  ggplot(., aes(x = Age/365.25, y = OHDG_SG, group = sex, color = sex, shape = sex))  +
  geom_point() +
  geom_smooth(aes(linetype = sex), method = "lm") +
  scale_color_manual(values = cbPalette) +
  labs(x = "Age", y = "8-OHdG ng/mL", title = "8-OHdG") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 18))

# 4b. m-l Isoprostanes ----
cross_pub %>%
  ggplot(., aes(x = Age/365.25, y = Iso_SG, group = sex, color = sex))  +
  geom_point() +
  geom_smooth(method = "loess") +
  labs(x = "Age", y = "15-F2t-Isop ng/mL") +
  ggtitle("Isoprostanes Among M-L Sample (non-sig)") +
  theme(plot.title = element_text(hjust = 0.5, size = 12))

# 4c. m-l TBARS ----
cross_pub %>%
  ggplot(., aes(x = Age/365.25, y = TBARS_SG_TAR, group = sex, color = sex))  +
  geom_point() +
  geom_smooth(method = "loess") +
  labs(x = "Age", y = "TBARS mM") +
  ggtitle("MDA-TBARS Among M-L Sample (non-sig)") +
  theme(plot.title = element_text(hjust = 0.5, size = 12))

# 4d. m-l Neopterin ----
cross_pub %>%
  ggplot(., aes(x = Age/365.25, y = Neo_SG_TAR, group = sex, color = sex))  +
  geom_point() +
  geom_smooth(method = "loess") +
  labs(x = "Age", y = "Neopterin nM") +
  ggtitle("Neopterin Among M-L Sample (non-sig)") +
  theme(plot.title = element_text(hjust = 0.5, size = 12))

# 4e. m-l TAC ----
cross_pub %>%
  ggplot(., aes(x = Age/365.25, y = TAC_SG_TAR, group = sex, color = sex))  +
  geom_point() +
  geom_smooth(method = "loess") +
  labs(x = "Age", y = "TAC mM") +
  ggtitle("TAC Among M-L Sample (non-sig)") +
  theme(plot.title = element_text(hjust = 0.5, size = 12))

# 5.  Lead up to death (LUTD) ----

# 5a. LUTD ohdg ----

ohdg_pop_median <- oxs_final_avg_pub %>%
  filter((!bdae_timing_ext %in% c("dur","aft") | is.na(bdae_timing_ext)) &
           (!bdai_timing_ext %in% c("dur","aft") | is.na(bdai_timing_ext)) & 
           (lutd != "1" | is.na(lutd))) %>% pull(OHDG_SG) %>% median(., na.rm = T)

oxs_final_avg_pub %>% filter(!is.na(lutd)) %>%
  filter(OHDG_SG < 50) %>%
  ggplot(., aes(x = lutd_years_until, y = OHDG_SG))  +
  geom_hline(yintercept = ohdg_pop_median, linetype = 3, color = "black") +
  geom_smooth(method = "lm", aes(color = ChimpID)) +
  geom_point(aes(color = ChimpID)) +
  labs(x = "Years Before Death", y = "8-OHdG ng/mL") +
  ggtitle("8-OHdG") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 18), legend.position = "none") +
  annotate("text", x = -7.8, y = 43, label = "A)", size = 10) +
  annotate("text", x = -4, y = 40, label = "*", size = 15)

# 5b. LUTD Isoprostanes ----

iso_pop_median <- oxs_final_avg_pub %>%
  filter((!bdae_timing_ext %in% c("dur","aft") | is.na(bdae_timing_ext)) &
           (!bdai_timing_ext %in% c("dur","aft") | is.na(bdai_timing_ext)) & 
           (lutd != "1" | is.na(lutd))) %>% pull(Iso_SG) %>% median(., na.rm = T)

oxs_final_avg_pub %>% filter(!is.na(lutd)) %>%
  filter(Iso_SG < 10) %>%
  ggplot(., aes(x = lutd_years_until, y = Iso_SG))  +
  geom_hline(yintercept = iso_pop_median, linetype = 3, color = "black") +
  geom_smooth(method = "lm", aes(color = ChimpID)) +
  geom_point(aes(color = ChimpID)) +
  labs(x = "Years Before Death", y = "Isoprostane ng/mL") +
  ggtitle("Isoprostanes") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 18), legend.position = "none") +
  annotate("text", x = -7.8, y = 8, label = "B)", size = 10) +
  annotate("text", x = -4, y = 7.5, label = "*", size = 15)

# 5c. LUTD TBARS ----

tbars_pop_median <- oxs_final_avg_pub %>%
  filter((!bdae_timing_ext %in% c("dur","aft") | is.na(bdae_timing_ext)) &
           (!bdai_timing_ext %in% c("dur","aft") | is.na(bdai_timing_ext)) & 
           (lutd != "1" | is.na(lutd))) %>% pull(TBARS_SG_TAR) %>% median(., na.rm = T)

oxs_final_avg_pub %>% filter(!is.na(lutd)) %>%
  ggplot(., aes(x = lutd_years_until, y = TBARS_SG_TAR))  +
  geom_hline(yintercept = tbars_pop_median, linetype = 3, color = "black") +
  geom_smooth(method = "lm", aes(color = ChimpID)) +
  geom_point(aes(color = ChimpID)) +
  labs(x = "Years Before Death", y = "MDA-TBARS nM (time adjusted residual)") +
  ggtitle("MDA-TBARS") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 18), legend.position = "none") +
  annotate("text", x = -8, y = 38, label = "C)", size = 10)


# 5d. LUTD Neopterin ----
neo_pop_median <- oxs_final_avg_pub %>%
  filter((!bdae_timing_ext %in% c("dur","aft") | is.na(bdae_timing_ext)) &
           (!bdai_timing_ext %in% c("dur","aft") | is.na(bdai_timing_ext)) & 
           (lutd != "1" | is.na(lutd))) %>% pull(Neo_SG_TAR) %>% median(., na.rm = T)

oxs_final_avg_pub %>% filter(!is.na(lutd)) %>%
  filter(Neo_SG < 3000) %>%
  ggplot(., aes(x = lutd_years_until, y = Neo_SG_TAR))  +
  geom_hline(yintercept = neo_pop_median, linetype = 3, color = "black") +
  geom_smooth(method = "lm", aes(color = ChimpID)) +
  geom_point(aes(color = ChimpID)) +
  labs(x = "Years Before Death", y = "Neopterin nM (time adjusted residual)") +
  ggtitle("Neopterin") + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 18), legend.position = "none") +
  annotate("text", x = -7.8, y = 1900, label = "D)", size = 10) +
  annotate("text", x = -4, y = 1600, label = "*", size = 15) +
  annotate("text", x = -3.3, y = 1600, label = "(no KK)", size = 4)

# 5e. LUTD TAC ----
tac_pop_median <- oxs_final_avg_pub %>%
  filter((!bdae_timing_ext %in% c("dur","aft") | is.na(bdae_timing_ext)) &
           (!bdai_timing_ext %in% c("dur","aft") | is.na(bdai_timing_ext)) & 
           (lutd != "1" | is.na(lutd))) %>% pull(TAC_SG_TAR) %>% median(., na.rm = T)

oxs_final_avg_pub %>% filter(!is.na(lutd)) %>%
  ggplot(., aes(x = lutd_years_until, y = TAC_SG_TAR))  +
  geom_hline(yintercept = tac_pop_median, linetype = 3, color = "black") +
  geom_smooth(method = "lm", aes(color = ChimpID)) +
  geom_point(aes(color = ChimpID)) +
  labs(x = "Years Before Death", y = "TAC mM (time adjusted residual)") +
  ggtitle("TAC") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 18), legend.position = "bottom",
        legend.text = element_text(size = 12), legend.title = element_text(size = 12)) +
  annotate("text", x = -7.8, y = 12.5, label = "E)", size = 10)



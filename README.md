# Recombinastics
Advances with the GA recombination site

Recombinastics_analysis
================
Kenneth Matreyek and Nisha D. Kamath
initialized 6/17/2020 - last updated 6/7/2023

``` r
rm(list = ls())
library(tidyverse)
```

    ## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.2 ──
    ## ✔ ggplot2 3.4.0      ✔ purrr   0.3.5 
    ## ✔ tibble  3.1.8      ✔ dplyr   1.0.10
    ## ✔ tidyr   1.2.1      ✔ stringr 1.5.0 
    ## ✔ readr   2.1.3      ✔ forcats 0.5.2 
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()

``` r
library(ggrepel)
library(reshape)
```

    ## 
    ## Attaching package: 'reshape'
    ## 
    ## The following object is masked from 'package:dplyr':
    ## 
    ##     rename
    ## 
    ## The following objects are masked from 'package:tidyr':
    ## 
    ##     expand, smiths

``` r
library(abind)
library(gridExtra)
```

    ## 
    ## Attaching package: 'gridExtra'
    ## 
    ## The following object is masked from 'package:dplyr':
    ## 
    ##     combine

``` r
library(ggbeeswarm)
library(factoextra)
```

    ## Welcome! Want to learn more? See two factoextra-related books at https://goo.gl/ve3WBa

``` r
library(ggfortify)
library(patchwork)

theme_set(theme_bw())
theme_update(panel.grid.minor = element_blank())

R.Version()
```

    ## $platform
    ## [1] "aarch64-apple-darwin20"
    ## 
    ## $arch
    ## [1] "aarch64"
    ## 
    ## $os
    ## [1] "darwin20"
    ## 
    ## $system
    ## [1] "aarch64, darwin20"
    ## 
    ## $status
    ## [1] ""
    ## 
    ## $major
    ## [1] "4"
    ## 
    ## $minor
    ## [1] "2.2"
    ## 
    ## $year
    ## [1] "2022"
    ## 
    ## $month
    ## [1] "10"
    ## 
    ## $day
    ## [1] "31"
    ## 
    ## $`svn rev`
    ## [1] "83211"
    ## 
    ## $language
    ## [1] "R"
    ## 
    ## $version.string
    ## [1] "R version 4.2.2 (2022-10-31)"
    ## 
    ## $nickname
    ## [1] "Innocent and Trusting"

## Look at the recombination method test data. Note: this is not part of the manuscript

``` r
methods_data <- read.csv(file = "data/Recombination_method_tests.csv", header = T)

methods_data_replicates <- ncol(methods_data) - 5
methods_data$mean <- rowMeans(methods_data[,c("Recombined_Rep1","Recombined_Rep2","Recombined_Rep3","Recombined_Rep4")], na.rm = T)
methods_data$sd <- sqrt((methods_data$Recombined_Rep1 - methods_data$mean)^2 + 
                        (methods_data$Recombined_Rep2 - methods_data$mean)^2 + 
                        (methods_data$Recombined_Rep3 - methods_data$mean)^2 +
                        (methods_data$Recombined_Rep4 - methods_data$mean)^2) 
methods_data$se <- methods_data$sd / sqrt(methods_data_replicates - 1)

methods_data$upper_conf <- methods_data$mean + methods_data$se * 1.96
methods_data$lower_conf <- methods_data$mean - methods_data$se * 1.96

methods_data2 <- methods_data %>% filter(!is.na(Recombined_Rep1) & Method != "none") %>% mutate(concat = paste("Method: ",Method,"\n","Bxb1: ",Bxb1,sep = ""))
methods_data2$WellSize <- factor(methods_data2$WellSize)
methods_data2$concat <- factor(methods_data2$concat, levels = methods_data2$concat[c(1,3,2,4)])

methods_data2[methods_data2$lower_conf < 0,"lower_conf"] <- 0

Recombination_method_plot <- ggplot() + 
  theme_classic() + theme(axis.text.x = element_text(angle = -90, vjust = 0.5)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,65)) +
  ylab("% mCherry positive cells") +
  xlab(NULL) +
  geom_hline(yintercept = 0) + 
  geom_errorbar(data = methods_data2, aes(x = concat, ymin = lower_conf, ymax = upper_conf, color = WellSize), 
                width = 0.2, position=position_dodge(width=0.3), size = 0.3 ) +
  geom_point(data = methods_data2, aes(x= concat, y=mean, color = WellSize), 
             position=position_dodge(width=0.3), shape = 95, size = 2) +
  geom_jitter(data = methods_data2, aes(x= concat, y=Recombined_Rep1, color = WellSize), 
             size = 0.9, position=position_dodge(width=0.3), alpha = 0.4) +
  geom_jitter(data = methods_data2, aes(x= concat, y=Recombined_Rep2, color = WellSize), 
             size = 0.9, position=position_dodge(width=0.3), alpha = 0.4) +
  geom_jitter(data = methods_data2, aes(x= concat, y=Recombined_Rep3, color = WellSize), 
             size = 0.9, position=position_dodge(width=0.3), alpha = 0.4) +
  geom_jitter(data = methods_data2, aes(x= concat, y=Recombined_Rep4, color = WellSize), 
             size = 0.9, position=position_dodge(width=0.3), alpha = 0.4)
```

    ## Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
    ## ℹ Please use `linewidth` instead.

``` r
print(Recombination_method_plot)
```

![](Recombinastics_analysis_files/figure-gfm/Testing%20recombination%20efficiency%20in%20two%20different%20plate%20formats-1.png)<!-- -->
\## GT and GA Ortholognality \## Relevant to Figure 1 This next chunk is
for testing the orthogonality of the GT and GA Bxb1 sequences

``` r
rep_1_frame <- read.csv(file = "data/Orthog_recomb/GT_vs_GA_200206.csv", header = T, stringsAsFactors = F) %>% arrange(attp, bxb1, mix) %>% filter(bxb1 != "none")
rep_2_frame <- read.csv(file = "data/Orthog_recomb/GT_vs_GA_200213.csv", header = T, stringsAsFactors = F) %>% arrange(attp, bxb1, mix) %>% filter(bxb1 != "none")
rep_3_frame <- read.csv(file = "data/Orthog_recomb/GT_vs_GA_200219.csv", header = T, stringsAsFactors = F) %>% arrange(attp, bxb1, mix) %>% filter(bxb1 != "none")

replicates <- 2
rep_1 <- log10(rep_1_frame[,5:8])
rep_2 <- log10(rep_2_frame[,5:8])
rep_3 <- log10(rep_3_frame[,5:8])
means <- (rep_1 + rep_2 + rep_3)/3

#https://stackoverflow.com/questions/32609926/performing-element-wise-standard-deviation-in-r-with-two-matrices
m <- abind(rep_1, rep_2, rep_3, along=3)
standard_devs <- data.frame(apply(m, 1:2, sd))
standard_errors <-  standard_devs / sqrt(replicates)

upper_conf <- means + standard_errors * 1.96
lower_conf <- means - standard_errors * 1.96

label_frame <- rep_1_frame[,1:4] %>% mutate(comboname = paste(bxb1, mix, sep = "\n"))

means2 <- cbind(label_frame,means)
upper_conf2 <- cbind(label_frame,upper_conf)
lower_conf2 <- cbind(label_frame,lower_conf)
#ga_data_filtered$attp <- factor(ga_data_filtered$attp, levels = c("gt","ga"))

upper_conf2_melted <- melt(upper_conf2[,c("comboname","attp","gfp","mcherry")], id = c("comboname","attp"))
lower_conf2_melted <- melt(lower_conf2[,c("comboname","attp","gfp","mcherry")], id = c("comboname","attp"))
means2_melted <- melt(means2[,c("comboname","attp","gfp","mcherry")], id = c("comboname","attp"))


green_red_colorscale <- c(gfp = "green", mcherry = "red")

GT_vs_GA_plot <- ggplot() + 
  theme_classic() + theme(axis.text.x = element_text(angle = -90, vjust = 0.5)) +
  ylab("% cells recombined") +
  xlab(NULL) +
  scale_y_log10() +
  scale_color_manual(values = green_red_colorscale) + 
  geom_errorbar(data = upper_conf2_melted, aes(x = comboname, color = variable, ymin = 10^lower_conf2_melted$value, ymax = 10^value),
                alpha = 0.5, width = 0.2, position = position_dodge(width = 0.5)) +
  geom_point(data = means2_melted, aes(x = comboname, y = 10^value, color = variable), alpha = 0.5, position = position_dodge(width = 0.5)) +
  facet_wrap(~attp)
print(GT_vs_GA_plot)
```

![](Recombinastics_analysis_files/figure-gfm/Initial%20testing%20of%20orthogonality%20of%20GT%20and%20GA%20recombination%20site%20pairs-1.png)<!-- -->

``` r
ggsave(file = "Plots/GT_vs_GA_plot.pdf", GT_vs_GA_plot, height = 1.8, width = 3.8)
```

## Flanking recombinase sites to excise unwanted bacterial DNA sequences

## Relevant to Figure 2

``` r
## Making a data frame to keep track of how many unexcised cells are observed with G718A
flanking_flow_df <- data.frame("flow" = NA, "frac_unexcised" = NA)
flank_red_pos_cutoff <- 3e3


## F84
f84_none <- read.csv(file = "data/flow/Flanking/F84/F84_None_unselected_Sample(1).csv.gz") %>% mutate("flow" = "F84", sample = "None")
f84_g718a <- read.csv(file = "data/flow/Flanking/F84/F84_G718A_hygro_Sample(4).csv.gz") %>% mutate("flow" = "F84", sample = "Flanked")
f84_g747a <- read.csv(file = "data/flow/Flanking/F84/F84_G747A_hygro_Sample(5).csv.gz") %>% mutate("flow" = "F84", sample = "Control")

flank_expt_cell_num <- 90000
flank_red_pos_cutoff <- 5e3

f84 <- rbind(f84_none[1:flank_expt_cell_num,], f84_g718a[1:flank_expt_cell_num,], f84_g747a[1:flank_expt_cell_num,])

#ggplot() + theme_bw() + scale_x_log10(limits = c(1e1,1e5), expand = c(0,0)) + scale_y_log10(limits = c(1e1,1e5), expand = c(0,0)) + geom_vline(xintercept = flank_red_pos_cutoff, linetype = 2) + geom_point(data = f84, aes(x = YL2.A, y = BL1.A), alpha = 0.01) + facet_grid(rows = vars(sample))

f84$ratio <- f84$BL1.A / f84$YL2.A
f84_subset <- f84 %>% filter(YL2.A >= flank_red_pos_cutoff)

f84_control_95pct_interval <- c(quantile((f84_subset %>% filter(sample == "Control"))$ratio,0.05),quantile((f84_subset %>% filter(sample == "Control"))$ratio,0.975))

#ggplot() + theme_bw() + scale_x_continuous(limits = c(-0.05, 0.5)) + geom_histogram(data = f84_subset, aes(x = ratio), binwidth = 0.01) + facet_grid(rows = vars(sample)) + geom_vline(xintercept  = f84_control_95pct_interval, linetype = 2)

f84_fraction_unexcised <- sum((f84_subset %>% filter(sample == "Flanked"))$ratio > f84_control_95pct_interval[1])/nrow((f84_subset %>% filter(sample == "Flanked")))
flanking_flow_df <- rbind(flanking_flow_df, data.frame("flow" = "F84", "frac_unexcised" = f84_fraction_unexcised))

## F114
f114_g718a <- read.csv(file = "data/flow/Flanking/F114/F114_G718A_unselected_A2.csv.gz") %>% mutate("flow" = "F114", sample = "Flanked")
f114_g747a <- read.csv(file = "data/flow/Flanking/F114/F114_G747A_unselected_A3.csv.gz") %>% mutate("flow" = "F114", sample = "Control")

flank_expt_cell_num <- 45000
flank_red_pos_cutoff <- 3e3

f114 <- rbind(f114_g718a[1:flank_expt_cell_num,], f114_g747a[1:flank_expt_cell_num,])

#ggplot() + theme_bw() + scale_x_log10(limits = c(1e1,1e5), expand = c(0,0)) + scale_y_log10(limits = c(1e1,1e5), expand = c(0,0)) + geom_vline(xintercept = flank_red_pos_cutoff, linetype = 2) + geom_point(data = f114, aes(x = YL2.A, y = BL1.A), alpha = 0.1) + facet_grid(rows = vars(sample))

f114$ratio <- f114$BL1.A / f114$YL2.A
f114_subset <- f114 %>% filter(YL2.A >= flank_red_pos_cutoff)

f114_control_95pct_interval <- c(quantile((f114_subset %>% filter(sample == "Control"))$ratio,0.05),quantile((f114_subset %>% filter(sample == "Control"))$ratio,0.975))

#ggplot() + theme_bw() + scale_x_continuous(limits = c(-0.05, 0.5)) + geom_histogram(data = f114_subset, aes(x = ratio), binwidth = 0.01) + facet_grid(rows = vars(sample)) + geom_vline(xintercept  = f114_control_95pct_interval, linetype = 2)

f114_fraction_unexcised <- sum((f114_subset %>% filter(sample == "Flanked"))$ratio > f114_control_95pct_interval[1])/nrow((f114_subset %>% filter(sample == "Flanked")))
flanking_flow_df <- rbind(flanking_flow_df, data.frame("flow" = "F114", "frac_unexcised" = f114_fraction_unexcised))


## F116
f116_g718a <- read.csv(file = "data/flow/Flanking/F116/F116_G718A_unselected_H11.csv.gz") %>% mutate("flow" = "F116", sample = "Flanked")
f116_g747a <- read.csv(file = "data/flow/Flanking/F116/F116_G747A_unselected_H12.csv.gz") %>% mutate("flow" = "F116", sample = "Control")
flank_expt_cell_num <- 25000
flank_red_pos_cutoff <- 5e3
f116 <- rbind(f116_g718a[1:flank_expt_cell_num,], f116_g747a[1:flank_expt_cell_num,])

#ggplot() + theme_bw() + scale_x_log10(limits = c(1e1,1e5), expand = c(0,0)) + scale_y_log10(limits = c(1e1,1e5), expand = c(0,0)) + geom_vline(xintercept = flank_red_pos_cutoff, linetype = 2) + geom_point(data = f116, aes(x = YL2.A, y = BL1.A), alpha = 0.1) + facet_grid(rows = vars(sample))

f116$ratio <- f116$BL1.A / f116$YL2.A
f116_subset <- f116 %>% filter(YL2.A >= flank_red_pos_cutoff)

f116_control_95pct_interval <- c(quantile((f116_subset %>% filter(sample == "Control"))$ratio,0.05),quantile((f116_subset %>% filter(sample == "Control"))$ratio,0.975))

#ggplot() + theme_bw() + scale_x_continuous(limits = c(-0.05, 0.5)) + geom_histogram(data = f116_subset, aes(x = ratio), binwidth = 0.01) + facet_grid(rows = vars(sample)) + geom_vline(xintercept  = f116_control_95pct_interval, linetype = 2)

f116_fraction_unexcised <- sum((f116_subset %>% filter(sample == "Flanked"))$ratio > f116_control_95pct_interval[1])/nrow((f116_subset %>% filter(sample == "Flanked")))
flanking_flow_df <- rbind(flanking_flow_df, data.frame("flow" = "F116", "frac_unexcised" = f116_fraction_unexcised))


## F119
f119_none <- read.csv(file = "data/flow/Flanking/F119/F119_none_unselected_G11.csv.gz") %>% mutate("flow" = "F119", sample = "None")
f119_g718a <- read.csv(file = "data/flow/Flanking/F119/F119_G718A_unselected_H11.csv.gz") %>% mutate("flow" = "F119", sample = "Flanked")
f119_g747a <- read.csv(file = "data/flow/Flanking/F119/F119_G747A_unselected_H12.csv.gz") %>% mutate("flow" = "F119", sample = "Control")
flank_expt_cell_num <- 9000
flank_red_pos_cutoff <- 5e3
f119 <- rbind(f119_none[1:flank_expt_cell_num,], f119_g718a[1:flank_expt_cell_num,], f119_g747a[1:flank_expt_cell_num,]) %>% filter(!is.na(sample))

#ggplot() + theme_bw() + scale_x_log10(limits = c(1e1,1e5), expand = c(0,0)) + scale_y_log10(limits = c(1e1,1e5), expand = c(0,0)) + geom_vline(xintercept = flank_red_pos_cutoff, linetype = 2) + geom_point(data = f119, aes(x = YL2.A, y = BL1.A), alpha = 0.1) + facet_grid(rows = vars(sample))

f119$ratio <- f119$BL1.A / f119$YL2.A
f119_subset <- f119 %>% filter(YL2.A >= flank_red_pos_cutoff)

f119_control_95pct_interval <- c(quantile((f119_subset %>% filter(sample == "Control"))$ratio,0.05),quantile((f119_subset %>% filter(sample == "Control"))$ratio,0.975))

#ggplot() + theme_bw() + scale_x_continuous(limits = c(-0.05, 0.5)) + geom_histogram(data = f119_subset, aes(x = ratio), binwidth = 0.01) + facet_grid(rows = vars(sample)) + geom_vline(xintercept  = f119_control_95pct_interval, linetype = 2)

f119_fraction_unexcised <- sum((f119_subset %>% filter(sample == "Flanked"))$ratio > f119_control_95pct_interval[1])/nrow((f119_subset %>% filter(sample == "Flanked")))
flanking_flow_df <- rbind(flanking_flow_df, data.frame("flow" = "F119", "frac_unexcised" = f119_fraction_unexcised))

## F130
f130_none <- read.csv(file = "data/flow/Flanking/F130/F130_none_unselected_D10.csv.gz") %>% mutate("flow" = "F130", sample = "None")
f130_g718a <- read.csv(file = "data/flow/Flanking/F130/F130_G718A_unselected_D11.csv.gz") %>% mutate("flow" = "F130", sample = "Flanked")
f130_g747a <- read.csv(file = "data/flow/Flanking/F130/F130_G747A_unselected_D12.csv.gz") %>% mutate("flow" = "F130", sample = "Control")
flank_expt_cell_num <- 89000
flank_red_pos_cutoff <- 5e3
f130 <- rbind(f130_none[1:flank_expt_cell_num,], f130_g718a[1:flank_expt_cell_num,], f130_g747a[1:flank_expt_cell_num,]) %>% filter(!is.na(sample))

#ggplot() + theme_bw() + scale_x_log10(limits = c(1e1,1e5), expand = c(0,0)) + scale_y_log10(limits = c(1e1,1e5), expand = c(0,0)) + geom_vline(xintercept = flank_red_pos_cutoff, linetype = 2) + geom_point(data = f130, aes(x = YL2.A, y = BL1.A), alpha = 0.1) + facet_grid(rows = vars(sample))

f130$ratio <- f130$BL1.A / f130$YL2.A
f130_subset <- f130 %>% filter(YL2.A >= flank_red_pos_cutoff)

f130_control_95pct_interval <- c(quantile((f130_subset %>% filter(sample == "Control"))$ratio,0.05),quantile((f130_subset %>% filter(sample == "Control"))$ratio,0.975))

#ggplot() + theme_bw() + scale_x_continuous(limits = c(-0.05, 0.5)) + geom_histogram(data = f130_subset, aes(x = ratio), binwidth = 0.01) + facet_grid(rows = vars(sample)) + geom_vline(xintercept  = f130_control_95pct_interval, linetype = 2)

f130_fraction_unexcised <- sum((f130_subset %>% filter(sample == "Flanked"))$ratio > f130_control_95pct_interval[1])/nrow((f130_subset %>% filter(sample == "Flanked")))
flanking_flow_df <- rbind(flanking_flow_df, data.frame("flow" = "F130", "frac_unexcised" = f130_fraction_unexcised))

## F131
f131_none <- read.csv(file = "data/flow/Flanking/F131/F131_none_unselected_D4.csv.gz") %>% mutate("flow" = "F131", sample = "None")
f131_g718a <- read.csv(file = "data/flow/Flanking/F131/F131_G718A_unselected_D5.csv.gz") %>% mutate("flow" = "F131", sample = "Flanked")
f131_g747a <- read.csv(file = "data/flow/Flanking/F131/F131_G747A_unselected_D6.csv.gz") %>% mutate("flow" = "F131", sample = "Control")
flank_expt_cell_num <- 190000
flank_red_pos_cutoff <- 3e3
f131 <- rbind(f131_none[1:flank_expt_cell_num,], f131_g718a[1:flank_expt_cell_num,], f131_g747a[1:flank_expt_cell_num,]) %>% filter(!is.na(sample))

#ggplot() + theme_bw() + scale_x_log10(limits = c(1e1,1e5), expand = c(0,0)) + scale_y_log10(limits = c(1e1,1e5), expand = c(0,0)) + geom_vline(xintercept = flank_red_pos_cutoff, linetype = 2) + geom_point(data = f131, aes(x = YL2.A, y = BL1.A), alpha = 0.1) + facet_grid(rows = vars(sample))

f131$ratio <- f131$BL1.A / f131$YL2.A
f131_subset <- f131 %>% filter(YL2.A >= flank_red_pos_cutoff)

f131_control_95pct_interval <- c(quantile((f131_subset %>% filter(sample == "Control"))$ratio,0.05),quantile((f131_subset %>% filter(sample == "Control"))$ratio,0.975))

#ggplot() + theme_bw() + scale_x_continuous(limits = c(-0.05, 0.5)) + geom_histogram(data = f131_subset, aes(x = ratio), binwidth = 0.01) + facet_grid(rows = vars(sample)) + geom_vline(xintercept  = f131_control_95pct_interval, linetype = 2)

f131_fraction_unexcised <- sum((f131_subset %>% filter(sample == "Flanked"))$ratio > f131_control_95pct_interval[1])/nrow((f131_subset %>% filter(sample == "Flanked")))
flanking_flow_df <- rbind(flanking_flow_df, data.frame("flow" = "F131", "frac_unexcised" = f131_fraction_unexcised))

## F132
f132_g718a <- read.csv(file = "data/flow/Flanking/F132/F132_G718A_unselected_G12.csv.gz") %>% mutate("flow" = "F132", sample = "Flanked")
f132_g747a <- read.csv(file = "data/flow/Flanking/F132/F132_G747A_unselected_H1.csv.gz") %>% mutate("flow" = "F132", sample = "Control")
flank_expt_cell_num <- 180000
flank_red_pos_cutoff <- 3e3
f132 <- rbind(f132_g718a[1:flank_expt_cell_num,], f132_g747a[1:flank_expt_cell_num,])

#ggplot() + theme_bw() + scale_x_log10(limits = c(1e1,1e5), expand = c(0,0)) + scale_y_log10(limits = c(1e1,1e5), expand = c(0,0)) + geom_vline(xintercept = flank_red_pos_cutoff, linetype = 2) + geom_point(data = f132, aes(x = YL2.A, y = BL1.A), alpha = 0.1) + facet_grid(rows = vars(sample))

f132$ratio <- f132$BL1.A / f132$YL2.A
f132_subset <- f132 %>% filter(YL2.A >= flank_red_pos_cutoff)

f132_control_95pct_interval <- c(quantile((f132_subset %>% filter(sample == "Control"))$ratio,0.05),quantile((f132_subset %>% filter(sample == "Control"))$ratio,0.975))

#ggplot() + theme_bw() + scale_x_continuous(limits = c(-0.05, 0.5)) + geom_histogram(data = f132_subset, aes(x = ratio), binwidth = 0.01) + facet_grid(rows = vars(sample)) + geom_vline(xintercept  = f132_control_95pct_interval, linetype = 2)

f132_fraction_unexcised <- sum((f132_subset %>% filter(sample == "Flanked"))$ratio > f132_control_95pct_interval[1])/nrow((f132_subset %>% filter(sample == "Flanked")))
flanking_flow_df <- rbind(flanking_flow_df, data.frame("flow" = "F132", "frac_unexcised" = f132_fraction_unexcised))

## f280
f280_g718a <- read.csv(file = "data/flow/Flanking/F280/F280_G718A_unselected.csv.gz") %>% mutate("flow" = "F280", sample = "Flanked")
f280_g747a <- read.csv(file = "data/flow/Flanking/F280/F280_G747A_unselected.csv.gz") %>% mutate("flow" = "F280", sample = "Control")
flank_expt_cell_num <- 180000
flank_red_pos_cutoff <- 3e3
f280 <- rbind(f280_g718a[1:flank_expt_cell_num,], f280_g747a[1:flank_expt_cell_num,])

#ggplot() + theme_bw() + scale_x_log10(limits = c(1e1,1e5), expand = c(0,0)) + scale_y_log10(limits = c(1e1,1e5), expand = c(0,0)) + geom_vline(xintercept = flank_red_pos_cutoff, linetype = 2) + geom_point(data = f280, aes(x = YL2.A, y = BL1.A), alpha = 0.1) + facet_grid(rows = vars(sample))

f280$ratio <- f280$BL1.A / f280$YL2.A
f280_subset <- f280 %>% filter(YL2.A >= flank_red_pos_cutoff)

f280_control_95pct_interval <- c(quantile((f280_subset %>% filter(sample == "Control"))$ratio,0.05),quantile((f280_subset %>% filter(sample == "Control"))$ratio,0.975))

#ggplot() + theme_bw() + scale_x_continuous(limits = c(-0.05, 0.5)) + geom_histogram(data = f280_subset, aes(x = ratio), binwidth = 0.01) + facet_grid(rows = vars(sample)) + geom_vline(xintercept  = f280_control_95pct_interval, linetype = 2)

f280_fraction_unexcised <- sum((f280_subset %>% filter(sample == "Flanked"))$ratio > f280_control_95pct_interval[1])/nrow((f280_subset %>% filter(sample == "Flanked")))
flanking_flow_df <- rbind(flanking_flow_df, data.frame("flow" = "F280", "frac_unexcised" = f280_fraction_unexcised))


## f281
f281_g718a <- read.csv(file = "data/flow/Flanking/F281/F281_G718A_unselected_A1.csv.gz") %>% mutate("flow" = "F281", sample = "Flanked")
f281_g747a <- read.csv(file = "data/flow/Flanking/F281/F281_G747A_unselected_A2.csv.gz") %>% mutate("flow" = "F281", sample = "Control")
flank_expt_cell_num <- 310000
flank_red_pos_cutoff <- 3e3
f281 <- rbind(f281_g718a[1:flank_expt_cell_num,], f281_g747a[1:flank_expt_cell_num,])

#ggplot() + theme_bw() + scale_x_log10(limits = c(1e1,1e5), expand = c(0,0)) + scale_y_log10(limits = c(1e1,1e5), expand = c(0,0)) + geom_vline(xintercept = flank_red_pos_cutoff, linetype = 2) + geom_point(data = f281, aes(x = YL2.A, y = BL1.A), alpha = 0.1) + facet_grid(rows = vars(sample))

f281$ratio <- f281$BL1.A / f281$YL2.A
f281_subset <- f281 %>% filter(YL2.A >= flank_red_pos_cutoff)

f281_control_95pct_interval <- c(quantile((f281_subset %>% filter(sample == "Control"))$ratio,0.05),quantile((f281_subset %>% filter(sample == "Control"))$ratio,0.975))

#ggplot() + theme_bw() + scale_x_continuous(limits = c(-0.05, 0.5)) + geom_histogram(data = f281_subset, aes(x = ratio), binwidth = 0.01) + facet_grid(rows = vars(sample)) + geom_vline(xintercept  = f281_control_95pct_interval, linetype = 2)

f281_fraction_unexcised <- sum((f281_subset %>% filter(sample == "Flanked"))$ratio > f281_control_95pct_interval[1])/nrow((f281_subset %>% filter(sample == "Flanked")))
flanking_flow_df <- rbind(flanking_flow_df, data.frame("flow" = "F281", "frac_unexcised" = f281_fraction_unexcised))


## f281b
f281b_g718a <- read.csv(file = "data/flow/Flanking/F281/F281b_G718A_unselected_A3.csv.gz") %>% mutate("flow" = "F281b", sample = "Flanked")
flank_expt_cell_num <- 230000
flank_red_pos_cutoff <- 3e3
f281b <- rbind(f281b_g718a[1:flank_expt_cell_num,], f281_g747a[1:flank_expt_cell_num,])

#ggplot() + theme_bw() + scale_x_log10(limits = c(1e1,1e5), expand = c(0,0)) + scale_y_log10(limits = c(1e1,1e5), expand = c(0,0)) + geom_vline(xintercept = flank_red_pos_cutoff, linetype = 2) + geom_point(data = f281b, aes(x = YL2.A, y = BL1.A), alpha = 0.1) + facet_grid(rows = vars(sample))

f281b$ratio <- f281b$BL1.A / f281b$YL2.A
f281b_subset <- f281b %>% filter(YL2.A >= flank_red_pos_cutoff)

f281b_control_95pct_interval <- c(quantile((f281b_subset %>% filter(sample == "Control"))$ratio,0.05),quantile((f281b_subset %>% filter(sample == "Control"))$ratio,0.975))

#ggplot() + theme_bw() + scale_x_continuous(limits = c(-0.05, 0.5)) + geom_histogram(data = f281b_subset, aes(x = ratio), binwidth = 0.01) + facet_grid(rows = vars(sample)) + geom_vline(xintercept  = f281b_control_95pct_interval, linetype = 2)

f281b_fraction_unexcised <- sum((f281b_subset %>% filter(sample == "Flanked"))$ratio > f281b_control_95pct_interval[1])/nrow((f281b_subset %>% filter(sample == "Flanked")))
flanking_flow_df <- rbind(flanking_flow_df, data.frame("flow" = "F281b", "frac_unexcised" = f281b_fraction_unexcised))


## f282
f282_g718a <- read.csv(file = "data/flow/Flanking/F282/F282_G718A_unselected_A1.csv.gz") %>% mutate("flow" = "F282", sample = "Flanked")
f282_g747a <- read.csv(file = "data/flow/Flanking/F282/F282_G747A_unselected_A2.csv.gz") %>% mutate("flow" = "F282", sample = "Control")
flank_expt_cell_num <- 280000
flank_red_pos_cutoff <- 3e3
f282 <- rbind(f282_g718a[1:flank_expt_cell_num,], f282_g747a[1:flank_expt_cell_num,])

#ggplot() + theme_bw() + scale_x_log10(limits = c(1e1,1e5), expand = c(0,0)) + scale_y_log10(limits = c(1e1,1e5), expand = c(0,0)) + geom_vline(xintercept = flank_red_pos_cutoff, linetype = 2) + geom_point(data = f282, aes(x = YL2.A, y = BL1.A), alpha = 0.1) + facet_grid(rows = vars(sample))

f282$ratio <- f282$BL1.A / f282$YL2.A
f282_subset <- f282 %>% filter(YL2.A >= flank_red_pos_cutoff)

f282_control_95pct_interval <- c(quantile((f282_subset %>% filter(sample == "Control"))$ratio,0.05),quantile((f282_subset %>% filter(sample == "Control"))$ratio,0.975))

#ggplot() + theme_bw() + scale_x_continuous(limits = c(-0.05, 0.5)) + geom_histogram(data = f282_subset, aes(x = ratio), binwidth = 0.01) + facet_grid(rows = vars(sample), scales="free_y") + geom_vline(xintercept  = f282_control_95pct_interval, linetype = 2)

f282_fraction_unexcised <- sum((f282_subset %>% filter(sample == "Flanked"))$ratio > f282_control_95pct_interval[1])/nrow((f282_subset %>% filter(sample == "Flanked")))
flanking_flow_df <- rbind(flanking_flow_df, data.frame("flow" = "F282", "frac_unexcised" = f282_fraction_unexcised))


## f282b
f282b_g718a <- read.csv(file = "data/flow/Flanking/F282/F282b_G718A_unselected_A3.csv.gz") %>% mutate("flow" = "F282b", sample = "Flanked")
flank_expt_cell_num <- 930000
flank_red_pos_cutoff <- 3e3
f282b <- rbind(f282b_g718a[1:flank_expt_cell_num,], f282_g747a[1:flank_expt_cell_num,])

#ggplot() + theme_bw() + scale_x_log10(limits = c(1e1,1e5), expand = c(0,0)) + scale_y_log10(limits = c(1e1,1e5), expand = c(0,0)) + geom_vline(xintercept = flank_red_pos_cutoff, linetype = 2) + geom_point(data = f282b, aes(x = YL2.A, y = BL1.A), alpha = 0.1) + facet_grid(rows = vars(sample))

f282b$ratio <- f282b$BL1.A / f282b$YL2.A
f282b_subset <- f282b %>% filter(YL2.A >= flank_red_pos_cutoff)

f282b_control_95pct_interval <- c(quantile((f282b_subset %>% filter(sample == "Control"))$ratio,0.05),quantile((f282b_subset %>% filter(sample == "Control"))$ratio,0.975))

#ggplot() + theme_bw() + scale_x_continuous(limits = c(-0.05, 0.5)) + geom_histogram(data = f282b_subset, aes(x = ratio), binwidth = 0.01) + facet_grid(rows = vars(sample), scales="free_y") + geom_vline(xintercept  = f282b_control_95pct_interval, linetype = 2)

f282b_fraction_unexcised <- sum((f282b_subset %>% filter(sample == "Flanked"))$ratio > f282b_control_95pct_interval[1])/nrow((f282b_subset %>% filter(sample == "Flanked")))
flanking_flow_df <- rbind(flanking_flow_df, data.frame("flow" = "F282b", "frac_unexcised" = f282b_fraction_unexcised))


## f283
f283_g718a <- read.csv(file = "data/flow/Flanking/F283/F283_G718A_unselected_A1.csv.gz") %>% mutate("flow" = "F283", sample = "Flanked")
f283_g747a <- read.csv(file = "data/flow/Flanking/F283/F283_G747A_unselected_A2.csv.gz") %>% mutate("flow" = "F283", sample = "Control")
flank_expt_cell_num <- 305000
flank_red_pos_cutoff <- 3e3
f283 <- rbind(f283_g718a[1:flank_expt_cell_num,], f283_g747a[1:flank_expt_cell_num,])

#ggplot() + theme_bw() + scale_x_log10(limits = c(1e1,1e5), expand = c(0,0)) + scale_y_log10(limits = c(1e1,1e5), expand = c(0,0)) + geom_vline(xintercept = flank_red_pos_cutoff, linetype = 2) +geom_point(data = f283, aes(x = YL2.A, y = BL1.A), alpha = 0.1) +facet_grid(rows = vars(sample))

f283$ratio <- f283$BL1.A / f283$YL2.A
f283_subset <- f283 %>% filter(YL2.A >= flank_red_pos_cutoff)

f283_control_95pct_interval <- c(quantile((f283_subset %>% filter(sample == "Control"))$ratio,0.05),quantile((f283_subset %>% filter(sample == "Control"))$ratio,0.975))

#ggplot() + theme_bw() + scale_x_continuous(limits = c(-0.05, 0.5)) + geom_histogram(data = f283_subset, aes(x = ratio), binwidth = 0.01) + facet_grid(rows = vars(sample), scales="free_y") + geom_vline(xintercept  = f283_control_95pct_interval, linetype = 2)

f283_fraction_unexcised <- sum((f283_subset %>% filter(sample == "Flanked"))$ratio > f283_control_95pct_interval[1])/nrow((f283_subset %>% filter(sample == "Flanked")))
flanking_flow_df <- rbind(flanking_flow_df, data.frame("flow" = "F283", "frac_unexcised" = f283_fraction_unexcised))


## f283b
f283b_g718a <- read.csv(file = "data/flow/Flanking/F283/F283b_G718A_unselected_A3.csv.gz") %>% mutate("flow" = "F283b", sample = "Flanked")
flank_expt_cell_num <- 80000
flank_red_pos_cutoff <- 3e3
f283b <- rbind(f283b_g718a[1:flank_expt_cell_num,], f283_g747a[1:flank_expt_cell_num,])

#ggplot() + theme_bw() + scale_x_log10(limits = c(1e1,1e5), expand = c(0,0)) + scale_y_log10(limits = c(1e1,1e5), expand = c(0,0)) + geom_vline(xintercept = flank_red_pos_cutoff, linetype = 2) + geom_point(data = f283b, aes(x = YL2.A, y = BL1.A), alpha = 0.1) + facet_grid(rows = vars(sample))

f283b$ratio <- f283b$BL1.A / f283b$YL2.A
f283b_subset <- f283b %>% filter(YL2.A >= flank_red_pos_cutoff)

f283b_control_95pct_interval <- c(quantile((f283b_subset %>% filter(sample == "Control"))$ratio,0.05),quantile((f283b_subset %>% filter(sample == "Control"))$ratio,0.975))

#ggplot() + theme_bw() + scale_x_continuous(limits = c(-0.05, 0.5)) + geom_histogram(data = f283b_subset, aes(x = ratio), binwidth = 0.01) + facet_grid(rows = vars(sample), scales="free_y") + geom_vline(xintercept  = f283b_control_95pct_interval, linetype = 2)

f283b_fraction_unexcised <- sum((f283b_subset %>% filter(sample == "Flanked"))$ratio > f283b_control_95pct_interval[1])/nrow((f283b_subset %>% filter(sample == "Flanked")))
flanking_flow_df <- rbind(flanking_flow_df, data.frame("flow" = "F283b", "frac_unexcised" = f283b_fraction_unexcised))
```

``` r
## Some summary graphs
## An example set of samples
f131_example <- rbind(f131_none[1:20000,], f131_g718a[1:20000,], f131_g747a[1:20000,]) %>% filter(!is.na(sample))
f131_example$sample <- factor(f131_example$sample, levels = c("None", "Flanked", "Control"))
F131_example_scatterplot <- ggplot() + theme(panel.grid.minor = element_blank()) + 
  scale_x_log10(limits = c(1e1,1e5), expand = c(0,0), breaks = c(1e2,1e4)) + 
  scale_y_log10(limits = c(1e1,1e4), expand = c(0,0), breaks = c(1e2,1e3)) + 
  geom_point(data = f131_example, aes(x = YL2.A, y = BL1.A), alpha = 0.1, size = 0.5) +
  facet_grid(cols = vars(sample))
ggsave(file = "plots/Flanking_example_scatterplot.pdf", F131_example_scatterplot, height = 1.6, width = 4.5)
```

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Removed 40467 rows containing missing values (`geom_point()`).

``` r
F131_example_scatterplot
```

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Removed 40467 rows containing missing values (`geom_point()`).

![](Recombinastics_analysis_files/figure-gfm/Some%20more%20flanking%20data%20-%20Flow%20Cyometry-1.png)<!-- -->

``` r
f131_example <- rbind(f131_none[1:194000,], f131_g718a[1:194000,], f131_g747a[1:194000,]) %>% filter(!is.na(sample))
f131_example$sample <- factor(f131_example$sample, levels = c("None", "Flanked", "Control"))
f131_example$ratio <- f131_example$BL1.A / f131_example$YL2.A
f131_example_subset <- f131_example %>% filter(YL2.A >= flank_red_pos_cutoff & sample != "None")

f131_example_control_95pct_interval <- c(quantile((f131_example_subset %>% filter(sample == "Control"))$ratio,0.05))

F131_ratio_histogram <- ggplot() + theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank()) + 
  labs(x = "Green / Red ratio", y = "Number\nof cells") +
  scale_x_continuous(limits = c(-0.05, 0.5), breaks = c(0,0.2,0.4)) + scale_y_continuous(breaks = c(0,100)) +
  geom_histogram(data = f131_example_subset, aes(x = ratio), binwidth = 0.01) + facet_grid(rows = vars(sample)) + geom_vline(xintercept  = f131_example_control_95pct_interval, linetype = 2)
ggsave(file = "Plots/F131_ratio_histogram.pdf", F131_ratio_histogram, height = 1.4, width = 1.8)
```

    ## Warning: Removed 39 rows containing non-finite values (`stat_bin()`).

    ## Warning: Removed 4 rows containing missing values (`geom_bar()`).

``` r
F131_ratio_histogram
```

    ## Warning: Removed 39 rows containing non-finite values (`stat_bin()`).
    ## Removed 4 rows containing missing values (`geom_bar()`).

![](Recombinastics_analysis_files/figure-gfm/Some%20more%20flanking%20data%20-%20Flow%20Cyometry-2.png)<!-- -->

``` r
## Looking at fraction excised over time
flanking_flow_days <- read.csv(file = "data/Flanking_expt_timepoints.csv", header = T, stringsAsFactors = F)
flanking_flow_df2 <- merge(flanking_flow_df, flanking_flow_days, by = "flow")

flanking_flow_df2$replicate <- as.factor(flanking_flow_df2$replicate)

Flanking_timeplot <- ggplot() + theme(panel.grid.minor = element_blank(), legend.position = "none") + 
  scale_y_log10(limits = c(0.001,1)) + scale_x_log10() +
  labs(x = "Days after\ntransfection", y = "Fraction of\ncells unexcised") +
  geom_point(data = flanking_flow_df2 %>% filter(instrument == "Attune"), aes(x = days, y = frac_unexcised, color = replicate), color = "red")
ggsave(file = "plots/Flanking_timeplot.pdf", Flanking_timeplot, height = 1.6, width = 2)
```

    ## Warning: Transformation introduced infinite values in continuous y-axis

``` r
Flanking_timeplot
```

    ## Warning: Transformation introduced infinite values in continuous y-axis

![](Recombinastics_analysis_files/figure-gfm/Flanking%20data%20averages-1.png)<!-- -->

``` r
## Looking at the reproducibility of recombination
flank_indep_recombs <- flanking_flow_df2 %>% filter(frac_unexcised != 0) %>% group_by(replicate) %>% summarize(frac_unexcised = 10^mean(log10(frac_unexcised)))
flank_flow_geomean <- 10^mean(log10(flank_indep_recombs$frac_unexcised))
flank_flow_upper_conf <- 10^(mean(log10(flank_indep_recombs$frac_unexcised)) + (sd(log10(flank_indep_recombs$frac_unexcised))/sqrt(nrow(flank_indep_recombs)) * 1.96))
flank_flow_lower_conf <- 10^(mean(log10(flank_indep_recombs$frac_unexcised)) - (sd(log10(flank_indep_recombs$frac_unexcised))/sqrt(nrow(flank_indep_recombs)) * 1.96))

Repeated_excisions_plot <- ggplot() + theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  labs(x = NULL, y = "Fraction of cells unexcised") + 
  scale_y_log10(limits = c(0.001,1)) + #scale_x_continuous(limits = c(-0.1,0.1)) +
  geom_point(data = flank_indep_recombs, aes(x = 0, y = frac_unexcised), alpha = 0.5, color = "Red") +
  geom_point(aes(x = 0, y = flank_flow_geomean), shape = 95, size = 10, color = "Black") +
  geom_errorbar(aes(x = 0, ymax = flank_flow_upper_conf, ymin = flank_flow_lower_conf), width = 0.01, alpha = 0.4)
ggsave(file = "plots/Repeated_excisions_plot.pdf", Repeated_excisions_plot, height = 1.1, width = 1)
Repeated_excisions_plot
```

![](Recombinastics_analysis_files/figure-gfm/Flanking%20data%20averages-2.png)<!-- -->

``` r
## f280
f280_g718a <- read.csv(file = "data/flow/Flanking/F280/F280_G718A_unselected.csv.gz") %>% mutate("flow" = "F280", sample = "Flanked", treatment = "none")
f280_g747a <- read.csv(file = "data/flow/Flanking/F280/F280_G747A_unselected.csv.gz") %>% mutate("flow" = "F280", sample = "Control", treatment = "none")
f280_g718a_ap1903 <- read.csv(file = "data/flow/Flanking/F280/F280_G718A_AP1903.csv.gz") %>% mutate("flow" = "F280", sample = "Flanked", treatment = "AP1903")
f280_g747a_ap1903 <- read.csv(file = "data/flow/Flanking/F280/F280_G747A_AP1903.csv.gz") %>% mutate("flow" = "F280", sample = "Control", treatment = "AP1903")
flank_expt_cell_num <- 374000
flank_red_pos_cutoff <- 3e3
f280 <- rbind(f280_g718a[1:flank_expt_cell_num,], f280_g747a[1:flank_expt_cell_num,], f280_g718a_ap1903[1:flank_expt_cell_num,], f280_g747a_ap1903[1:flank_expt_cell_num,])
#ggplot() + theme_bw() + scale_x_log10(limits = c(1e1,1e5), expand = c(0,0)) + scale_y_log10(limits = c(1e1,1e5), expand = c(0,0)) + geom_vline(xintercept = flank_red_pos_cutoff, linetype = 2) + geom_point(data = f280, aes(x = YL2.A, y = BL1.A), alpha = 0.1) + facet_grid(rows = vars(sample))
f280$ratio <- f280$BL1.A / f280$YL2.A
f280_subset <- f280 %>% filter(YL2.A >= flank_red_pos_cutoff)
f280_control_95pct_interval <- c(quantile((f280_subset %>% filter(sample == "Control"))$ratio,0.05),quantile((f280_subset %>% filter(sample == "Control"))$ratio,0.975))
f280_subset$label <- paste0(f280_subset$sample,"_",f280_subset$treatment)
#ggplot() + theme_bw() + scale_x_continuous(limits = c(-0.05, 0.5)) + geom_histogram(data = f280_subset, aes(x = ratio), binwidth = 0.01) + facet_wrap(~label, scales = "free_y") + geom_vline(xintercept  = f280_control_95pct_interval, linetype = 2)
f280_subset_summary1 <- f280_subset %>% group_by(label) %>% summarize(unexcised = sum(ratio > f280_control_95pct_interval))
```

    ## Warning in ratio > f280_control_95pct_interval: longer object length is not a
    ## multiple of shorter object length

    ## Warning in ratio > f280_control_95pct_interval: longer object length is not a
    ## multiple of shorter object length

``` r
f280_subset_summary2 <- f280_subset %>% group_by(label) %>% count()
f280_subset_summary <- merge(f280_subset_summary1, f280_subset_summary2, by = "label")
f280_subset_summary$frac_unexcised <- f280_subset_summary$unexcised / f280_subset_summary$n
f280_subset_summary$rep <- "F280"

## F281
f281_g718a <- read.csv(file = "data/flow/Flanking/F281/F281_G718A_unselected_A1.csv.gz") %>% mutate("flow" = "F281", sample = "Flanked", treatment = "none")
f281_g747a <- read.csv(file = "data/flow/Flanking/F281/F281_G747A_unselected_A2.csv.gz") %>% mutate("flow" = "F281", sample = "Control", treatment = "none")
f281_g718a_ap1903 <- read.csv(file = "data/flow/Flanking/F281/F281_G718A_AP1903_B1.csv.gz") %>% mutate("flow" = "F281", sample = "Flanked", treatment = "AP1903")
f281_g747a_ap1903 <- read.csv(file = "data/flow/Flanking/F281/F281_G747A_AP1903_B2.csv.gz") %>% mutate("flow" = "F281", sample = "Control", treatment = "AP1903")
flank_expt_cell_num <- 287000
flank_red_pos_cutoff <- 3e3
f281 <- rbind(f281_g718a[1:flank_expt_cell_num,], f281_g747a[1:flank_expt_cell_num,], f281_g718a_ap1903[1:flank_expt_cell_num,], f281_g747a_ap1903[1:flank_expt_cell_num,])
f281$ratio <- f281$BL1.A / f281$YL2.A
f281_subset <- f281 %>% filter(YL2.A >= flank_red_pos_cutoff)
f281_control_95pct_interval <- c(quantile((f281_subset %>% filter(sample == "Control"))$ratio,0.05),quantile((f281_subset %>% filter(sample == "Control"))$ratio,0.975))
f281_subset$label <- paste0(f281_subset$sample,"_",f281_subset$treatment)
#ggplot() + theme_bw() + scale_x_continuous(limits = c(-0.05, 0.5)) + geom_histogram(data = f281_subset, aes(x = ratio), binwidth = 0.01) + facet_wrap(~label, scales = "free_y") + geom_vline(xintercept  = f281_control_95pct_interval, linetype = 2)
f281_subset_summary1 <- f281_subset %>% group_by(label) %>% summarize(unexcised = sum(ratio > f281_control_95pct_interval))
```

    ## Warning in ratio > f281_control_95pct_interval: longer object length is not a
    ## multiple of shorter object length

``` r
f281_subset_summary2 <- f281_subset %>% group_by(label) %>% count()
f281_subset_summary <- merge(f281_subset_summary1, f281_subset_summary2, by = "label")
f281_subset_summary$frac_unexcised <- f281_subset_summary$unexcised / f281_subset_summary$n
f281_subset_summary$rep <- "F281"

## F281b
f281b_g718a <- read.csv(file = "data/flow/Flanking/F281/F281b_G718A_unselected_A3.csv.gz") %>% mutate("flow" = "F281b", sample = "Flanked", treatment = "none")
f281b_g747a <- read.csv(file = "data/flow/Flanking/F281/F281b_G747A_unselected_A4.csv.gz") %>% mutate("flow" = "F281b", sample = "Control", treatment = "none")
f281b_g718a_ap1903 <- read.csv(file = "data/flow/Flanking/F281/F281b_G718A_AP1903_B3.csv.gz") %>% mutate("flow" = "F281b", sample = "Flanked", treatment = "AP1903")
f281b_g747a_ap1903 <- read.csv(file = "data/flow/Flanking/F281/F281b_G747A_AP1903_B4.csv.gz") %>% mutate("flow" = "F281b", sample = "Control", treatment = "AP1903")
flank_expt_cell_num <- 157000
flank_red_pos_cutoff <- 3e3
f281b <- rbind(f281b_g718a[1:flank_expt_cell_num,], f281b_g747a[1:flank_expt_cell_num,], f281b_g718a_ap1903[1:flank_expt_cell_num,], f281b_g747a_ap1903[1:flank_expt_cell_num,])
f281b$ratio <- f281b$BL1.A / f281b$YL2.A
f281b_subset <- f281b %>% filter(YL2.A >= flank_red_pos_cutoff)
f281b_control_95pct_interval <- c(quantile((f281b_subset %>% filter(sample == "Control"))$ratio,0.05),quantile((f281b_subset %>% filter(sample == "Control"))$ratio,0.975))
f281b_subset$label <- paste0(f281b_subset$sample,"_",f281b_subset$treatment)
#ggplot() + theme_bw() + scale_x_continuous(limits = c(-0.05, 0.5)) + geom_histogram(data = f281b_subset, aes(x = ratio), binwidth = 0.01) + facet_wrap(~label, scales = "free_y") + geom_vline(xintercept  = f281b_control_95pct_interval, linetype = 2)
f281b_subset_summary1 <- f281b_subset %>% group_by(label) %>% summarize(unexcised = sum(ratio > f281b_control_95pct_interval))
```

    ## Warning in ratio > f281b_control_95pct_interval: longer object length is not a
    ## multiple of shorter object length

    ## Warning in ratio > f281b_control_95pct_interval: longer object length is not a
    ## multiple of shorter object length

    ## Warning in ratio > f281b_control_95pct_interval: longer object length is not a
    ## multiple of shorter object length

``` r
f281b_subset_summary2 <- f281b_subset %>% group_by(label) %>% count()
f281b_subset_summary <- merge(f281b_subset_summary1, f281b_subset_summary2, by = "label")
f281b_subset_summary$frac_unexcised <- f281b_subset_summary$unexcised / f281b_subset_summary$n
f281b_subset_summary$rep <- "F281b"


## F282
f282_g718a <- read.csv(file = "data/flow/Flanking/F282/F282_G718A_unselected_A1.csv.gz") %>% mutate("flow" = "F282", sample = "Flanked", treatment = "none")
f282_g747a <- read.csv(file = "data/flow/Flanking/F282/F282_G747A_unselected_A2.csv.gz") %>% mutate("flow" = "F282", sample = "Control", treatment = "none")
f282_g718a_ap1903 <- read.csv(file = "data/flow/Flanking/F282/F282_G718A_AP1903_B1.csv.gz") %>% mutate("flow" = "F282", sample = "Flanked", treatment = "AP1903")
f282_g747a_ap1903 <- read.csv(file = "data/flow/Flanking/F282/F282_G747A_AP1903_B2.csv.gz") %>% mutate("flow" = "F282", sample = "Control", treatment = "AP1903")
flank_expt_cell_num <- 190000
flank_red_pos_cutoff <- 3e3
f282 <- rbind(f282_g718a[1:flank_expt_cell_num,], f282_g747a[1:flank_expt_cell_num,], f282_g718a_ap1903[1:flank_expt_cell_num,], f282_g747a_ap1903[1:flank_expt_cell_num,])
f282$ratio <- f282$BL1.A / f282$YL2.A
f282_subset <- f282 %>% filter(YL2.A >= flank_red_pos_cutoff)
f282_control_95pct_interval <- c(quantile((f282_subset %>% filter(sample == "Control"))$ratio,0.05),quantile((f282_subset %>% filter(sample == "Control"))$ratio,0.975))
f282_subset$label <- paste0(f282_subset$sample,"_",f282_subset$treatment)
#ggplot() + theme_bw() + scale_x_continuous(limits = c(-0.05, 0.5)) + geom_histogram(data = f282_subset, aes(x = ratio), binwidth = 0.01) + facet_wrap(~label, scales = "free_y") + geom_vline(xintercept  = f282_control_95pct_interval, linetype = 2)
f282_subset_summary1 <- f282_subset %>% group_by(label) %>% summarize(unexcised = sum(ratio > f282_control_95pct_interval))
```

    ## Warning in ratio > f282_control_95pct_interval: longer object length is not a
    ## multiple of shorter object length

    ## Warning in ratio > f282_control_95pct_interval: longer object length is not a
    ## multiple of shorter object length

    ## Warning in ratio > f282_control_95pct_interval: longer object length is not a
    ## multiple of shorter object length

    ## Warning in ratio > f282_control_95pct_interval: longer object length is not a
    ## multiple of shorter object length

``` r
f282_subset_summary2 <- f282_subset %>% group_by(label) %>% count()
f282_subset_summary <- merge(f282_subset_summary1, f282_subset_summary2, by = "label")
f282_subset_summary$frac_unexcised <- f282_subset_summary$unexcised / f282_subset_summary$n
f282_subset_summary$rep <- "F282"

## F282b
f282b_g718a <- read.csv(file = "data/flow/Flanking/F282/F282b_G718A_unselected_A3.csv.gz") %>% mutate("flow" = "F282b", sample = "Flanked", treatment = "none")
f282b_g747a <- read.csv(file = "data/flow/Flanking/F282/F282b_G747A_unselected_A4.csv.gz") %>% mutate("flow" = "F282b", sample = "Control", treatment = "none")
f282b_g718a_ap1903 <- read.csv(file = "data/flow/Flanking/F282/F282b_G718A_AP1903_B3.csv.gz") %>% mutate("flow" = "F282b", sample = "Flanked", treatment = "AP1903")
f282b_g747a_ap1903 <- read.csv(file = "data/flow/Flanking/F282/F282b_G747A_AP1903_B4.csv.gz") %>% mutate("flow" = "F282b", sample = "Control", treatment = "AP1903")
flank_expt_cell_num <- 93000
flank_red_pos_cutoff <- 3e3
f282b <- rbind(f282b_g718a[1:flank_expt_cell_num,], f282b_g747a[1:flank_expt_cell_num,], f282b_g718a_ap1903[1:flank_expt_cell_num,], f282b_g747a_ap1903[1:flank_expt_cell_num,])
f282b$ratio <- f282b$BL1.A / f282b$YL2.A
f282b_subset <- f282b %>% filter(YL2.A >= flank_red_pos_cutoff)
f282b_control_95pct_interval <- c(quantile((f282b_subset %>% filter(sample == "Control"))$ratio,0.05),quantile((f282b_subset %>% filter(sample == "Control"))$ratio,0.975))
f282b_subset$label <- paste0(f282b_subset$sample,"_",f282b_subset$treatment)
#ggplot() + theme_bw() + scale_x_continuous(limits = c(-0.05, 0.5)) + geom_histogram(data = f282b_subset, aes(x = ratio), binwidth = 0.01) + facet_wrap(~label, scales = "free_y") + geom_vline(xintercept  = f282b_control_95pct_interval, linetype = 2)
f282b_subset_summary1 <- f282b_subset %>% group_by(label) %>% summarize(unexcised = sum(ratio > f282b_control_95pct_interval))
```

    ## Warning in ratio > f282b_control_95pct_interval: longer object length is not a
    ## multiple of shorter object length

    ## Warning in ratio > f282b_control_95pct_interval: longer object length is not a
    ## multiple of shorter object length

``` r
f282b_subset_summary2 <- f282b_subset %>% group_by(label) %>% count()
f282b_subset_summary <- merge(f282b_subset_summary1, f282b_subset_summary2, by = "label")
f282b_subset_summary$frac_unexcised <- f282b_subset_summary$unexcised / f282b_subset_summary$n
f282b_subset_summary$rep <- "F282b"


## F283
f283_g718a <- read.csv(file = "data/flow/Flanking/F283/F283_G718A_unselected_A1.csv.gz") %>% mutate("flow" = "F283", sample = "Flanked", treatment = "none")
f283_g747a <- read.csv(file = "data/flow/Flanking/F283/F283_G747A_unselected_A2.csv.gz") %>% mutate("flow" = "F283", sample = "Control", treatment = "none")
f283_g718a_ap1903 <- read.csv(file = "data/flow/Flanking/F283/F283_G718A_AP1903_B1.csv.gz") %>% mutate("flow" = "F283", sample = "Flanked", treatment = "AP1903")
f283_g747a_ap1903 <- read.csv(file = "data/flow/Flanking/F283/F283_G747A_AP1903_B2.csv.gz") %>% mutate("flow" = "F283", sample = "Control", treatment = "AP1903")
f283_g718a_hygro <- read.csv(file = "data/flow/Flanking/F283/F283_G718A_Hygro_C1.csv.gz") %>% mutate("flow" = "F283b", sample = "Flanked", treatment = "Hygro")
f283_g747a_hygro <- read.csv(file = "data/flow/Flanking/F283/F283_G747A_Hygro_C2.csv.gz") %>% mutate("flow" = "F283b", sample = "Control", treatment = "Hygro")
flank_expt_cell_num <- 12000
flank_red_pos_cutoff <- 3e3
f283 <- rbind(f283_g718a[1:flank_expt_cell_num,], f283_g747a[1:flank_expt_cell_num,], f283_g718a_ap1903[1:flank_expt_cell_num,], f283_g747a_ap1903[1:flank_expt_cell_num,], f283_g718a_hygro[1:flank_expt_cell_num,], f283_g747a_hygro[1:flank_expt_cell_num,])
f283$ratio <- f283$BL1.A / f283$YL2.A
f283_subset <- f283 %>% filter(YL2.A >= flank_red_pos_cutoff)
f283_control_95pct_interval <- c(quantile((f283_subset %>% filter(sample == "Control"))$ratio,0.05),quantile((f283_subset %>% filter(sample == "Control"))$ratio,0.975))
f283_subset$label <- paste0(f283_subset$sample,"_",f283_subset$treatment)
#ggplot() + theme_bw() + scale_x_continuous(limits = c(-0.05, 0.5)) + geom_histogram(data = f283_subset, aes(x = ratio), binwidth = 0.01) + facet_wrap(~label, scales = "free_y") + geom_vline(xintercept  = f283_control_95pct_interval, linetype = 2)
f283_subset_summary1 <- f283_subset %>% group_by(label) %>% summarize(unexcised = sum(ratio > f283_control_95pct_interval))
```

    ## Warning in ratio > f283_control_95pct_interval: longer object length is not a
    ## multiple of shorter object length

    ## Warning in ratio > f283_control_95pct_interval: longer object length is not a
    ## multiple of shorter object length

    ## Warning in ratio > f283_control_95pct_interval: longer object length is not a
    ## multiple of shorter object length

    ## Warning in ratio > f283_control_95pct_interval: longer object length is not a
    ## multiple of shorter object length

``` r
f283_subset_summary2 <- f283_subset %>% group_by(label) %>% count()
f283_subset_summary <- merge(f283_subset_summary1, f283_subset_summary2, by = "label")
f283_subset_summary$frac_unexcised <- f283_subset_summary$unexcised / f283_subset_summary$n
f283_subset_summary$rep <- "F283"

## F283b
f283b_g718a <- read.csv(file = "data/flow/Flanking/F283/F283b_G718A_unselected_A3.csv.gz") %>% mutate("flow" = "F283b", sample = "Flanked", treatment = "none")
f283b_g747a <- read.csv(file = "data/flow/Flanking/F283/F283b_G747A_unselected_A4.csv.gz") %>% mutate("flow" = "F283b", sample = "Control", treatment = "none")
f283b_g718a_ap1903 <- read.csv(file = "data/flow/Flanking/F283/F283b_G718A_AP1903_B3.csv.gz") %>% mutate("flow" = "F283b", sample = "Flanked", treatment = "AP1903")
f283b_g747a_ap1903 <- read.csv(file = "data/flow/Flanking/F283/F283b_G747A_AP1903_B4.csv.gz") %>% mutate("flow" = "F283b", sample = "Control", treatment = "AP1903")
flank_expt_cell_num <- 83000
flank_red_pos_cutoff <- 3e3
f283b <- rbind(f283b_g718a[1:flank_expt_cell_num,], f283b_g747a[1:flank_expt_cell_num,], f283b_g718a_ap1903[1:flank_expt_cell_num,], f283b_g747a_ap1903[1:flank_expt_cell_num,])
f283b$ratio <- f283b$BL1.A / f283b$YL2.A
f283b_subset <- f283b %>% filter(YL2.A >= flank_red_pos_cutoff)
f283b_control_95pct_interval <- c(quantile((f283b_subset %>% filter(sample == "Control"))$ratio,0.05),quantile((f283b_subset %>% filter(sample == "Control"))$ratio,0.975))
f283b_subset$label <- paste0(f283b_subset$sample,"_",f283b_subset$treatment)
#ggplot() + theme_bw() + scale_x_continuous(limits = c(-0.05, 0.5)) + geom_histogram(data = f283b_subset, aes(x = ratio), binwidth = 0.01) + facet_wrap(~label, scales = "free_y") + geom_vline(xintercept  = f283b_control_95pct_interval, linetype = 2)
f283b_subset_summary1 <- f283b_subset %>% group_by(label) %>% summarize(unexcised = sum(ratio > f283b_control_95pct_interval))
```

    ## Warning in ratio > f283b_control_95pct_interval: longer object length is not a
    ## multiple of shorter object length

    ## Warning in ratio > f283b_control_95pct_interval: longer object length is not a
    ## multiple of shorter object length

``` r
f283b_subset_summary2 <- f283b_subset %>% group_by(label) %>% count()
f283b_subset_summary <- merge(f283b_subset_summary1, f283b_subset_summary2, by = "label")
f283b_subset_summary$frac_unexcised <- f283b_subset_summary$unexcised / f283b_subset_summary$n
f283b_subset_summary$rep <- "F283b"

## F284
f284_g718a_ap1903_hygro <- read.csv(file = "data/flow/Flanking/F284/F284_G718A_AP1903_Hygro_A1.csv.gz") %>% mutate("flow" = "F284", sample = "Flanked", treatment = "AP1903_Hygro")
f284_g747a_ap1903_hygro <- read.csv(file = "data/flow/Flanking/F284/F284_G747A_AP1903_Hygro_A2.csv.gz") %>% mutate("flow" = "F284", sample = "Control", treatment = "AP1903_Hygro")
f284_g718a_hygro <- read.csv(file = "data/flow/Flanking/F284/F284_G718A_Hygro_B1.csv.gz") %>% mutate("flow" = "F284", sample = "Flanked", treatment = "Hygro")
f284_g747a_hygro <- read.csv(file = "data/flow/Flanking/F284/F284_G747A_Hygro_B2.csv.gz") %>% mutate("flow" = "F284", sample = "Control", treatment = "Hygro")
flank_expt_cell_num <- 51000
flank_red_pos_cutoff <- 3e3
f284 <- rbind(f284_g718a_ap1903_hygro[1:flank_expt_cell_num,], f284_g747a_ap1903_hygro[1:flank_expt_cell_num,], f284_g718a_hygro[1:flank_expt_cell_num,], f284_g747a_hygro[1:flank_expt_cell_num,])
f284$ratio <- f284$BL1.A / f284$YL2.A
f284_subset <- f284 %>% filter(YL2.A >= flank_red_pos_cutoff)
f284_control_95pct_interval <- c(quantile((f284_subset %>% filter(sample == "Control"))$ratio,0.05),quantile((f284_subset %>% filter(sample == "Control"))$ratio,0.975))
f284_subset$label <- paste0(f284_subset$sample,"_",f284_subset$treatment)
#ggplot() + theme_bw() + scale_x_continuous(limits = c(-0.05, 0.5)) + geom_histogram(data = f284_subset, aes(x = ratio), binwidth = 0.01) + facet_wrap(~label, scales = "free_y") + geom_vline(xintercept  = f284_control_95pct_interval, linetype = 2)
f284_subset_summary1 <- f284_subset %>% group_by(label) %>% summarize(unexcised = sum(ratio > f284_control_95pct_interval))
```

    ## Warning in ratio > f284_control_95pct_interval: longer object length is not a
    ## multiple of shorter object length

``` r
f284_subset_summary2 <- f284_subset %>% group_by(label) %>% count()
f284_subset_summary <- merge(f284_subset_summary1, f284_subset_summary2, by = "label")
f284_subset_summary$frac_unexcised <- f284_subset_summary$unexcised / f284_subset_summary$n
f284_subset_summary$rep <- "F284"

## F284b
f284b_g718a_hygro <- read.csv(file = "data/flow/Flanking/F284/F284b_G718A_Hygro_B3.csv.gz") %>% mutate("flow" = "F284b", sample = "Flanked", treatment = "Hygro")
f284b_g747a_hygro <- read.csv(file = "data/flow/Flanking/F284/F284b_G747A_Hygro_B4.csv.gz") %>% mutate("flow" = "F284b", sample = "Control", treatment = "Hygro")
f284b_g718a_ap1903_hygro <- read.csv(file = "data/flow/Flanking/F284/F284b_G718A_AP1903_Hygro_A3.csv.gz") %>% mutate("flow" = "F284b", sample = "Flanked", treatment = "AP1903_Hygro")
f284b_g747a_ap1903_hygro <- read.csv(file = "data/flow/Flanking/F284/F284b_G747A_AP1903_Hygro_A4.csv.gz") %>% mutate("flow" = "F284b", sample = "Control", treatment = "AP1903_Hygro")
flank_expt_cell_num <- 35000
flank_red_pos_cutoff <- 3e3
f284b <- rbind(f284b_g718a_hygro[1:flank_expt_cell_num,], f284b_g747a_hygro[1:flank_expt_cell_num,], f284b_g718a_ap1903_hygro[1:flank_expt_cell_num,], f284b_g747a_ap1903_hygro[1:flank_expt_cell_num,])
f284b$ratio <- f284b$BL1.A / f284b$YL2.A
f284b_subset <- f284b %>% filter(YL2.A >= flank_red_pos_cutoff)
f284b_control_95pct_interval <- c(quantile((f284b_subset %>% filter(sample == "Control"))$ratio,0.05),quantile((f284b_subset %>% filter(sample == "Control"))$ratio,0.975))
f284b_subset$label <- paste0(f284b_subset$sample,"_",f284b_subset$treatment)
#ggplot() + theme_bw() + scale_x_continuous(limits = c(-0.05, 0.5)) + geom_histogram(data = f284b_subset, aes(x = ratio), binwidth = 0.01) + facet_wrap(~label, scales = "free_y") + geom_vline(xintercept  = f284b_control_95pct_interval, linetype = 2)
f284b_subset_summary1 <- f284b_subset %>% group_by(label) %>% summarize(unexcised = sum(ratio > f284b_control_95pct_interval))
```

    ## Warning in ratio > f284b_control_95pct_interval: longer object length is not a
    ## multiple of shorter object length

    ## Warning in ratio > f284b_control_95pct_interval: longer object length is not a
    ## multiple of shorter object length

``` r
f284b_subset_summary2 <- f284b_subset %>% group_by(label) %>% count()
f284b_subset_summary <- merge(f284b_subset_summary1, f284b_subset_summary2, by = "label")
f284b_subset_summary$frac_unexcised <- f284b_subset_summary$unexcised / f284b_subset_summary$n
f284b_subset_summary$rep <- "F284b"


## F285
f285_g718a_hygro <- read.csv(file = "data/flow/Flanking/F285/F285_G718A_Hygro_A3.csv.gz") %>% mutate("flow" = "F285", sample = "Flanked", treatment = "Hygro")
f285_g747a_hygro <- read.csv(file = "data/flow/Flanking/F285/F285_G747A_Hygro_A4.csv.gz") %>% mutate("flow" = "F285", sample = "Control", treatment = "Hygro")
flank_expt_cell_num <- 235000
flank_red_pos_cutoff <- 3e3
f285 <- rbind(f285_g718a_hygro[1:flank_expt_cell_num,], f285_g747a_hygro[1:flank_expt_cell_num,])
f285$ratio <- f285$BL1.A / f285$YL2.A
f285_subset <- f285 %>% filter(YL2.A >= flank_red_pos_cutoff)
f285_control_95pct_interval <- c(quantile((f285_subset %>% filter(sample == "Control"))$ratio,0.05),quantile((f285_subset %>% filter(sample == "Control"))$ratio,0.975))
f285_subset$label <- paste0(f285_subset$sample,"_",f285_subset$treatment)
#ggplot() + theme_bw() + scale_x_continuous(limits = c(-0.05, 0.5)) + geom_histogram(data = f285_subset, aes(x = ratio), binwidth = 0.01) + facet_wrap(~label, scales = "free_y") + geom_vline(xintercept  = f285_control_95pct_interval, linetype = 2)
f285_subset_summary1 <- f285_subset %>% group_by(label) %>% summarize(unexcised = sum(ratio > f285_control_95pct_interval))
```

    ## Warning in ratio > f285_control_95pct_interval: longer object length is not a
    ## multiple of shorter object length

    ## Warning in ratio > f285_control_95pct_interval: longer object length is not a
    ## multiple of shorter object length

``` r
f285_subset_summary2 <- f285_subset %>% group_by(label) %>% count()
f285_subset_summary <- merge(f285_subset_summary1, f285_subset_summary2, by = "label")
f285_subset_summary$frac_unexcised <- f285_subset_summary$unexcised / f285_subset_summary$n
f285_subset_summary$rep <- "F285"

## F286
f286_g718a_hygro <- read.csv(file = "data/flow/Flanking/F286/F286_G718A_Hygro_B1.csv.gz") %>% mutate("flow" = "F286", sample = "Flanked", treatment = "Hygro")
f286_g747a_hygro <- read.csv(file = "data/flow/Flanking/F286/F286_G747A_Hygro_B2.csv.gz") %>% mutate("flow" = "F286", sample = "Control", treatment = "Hygro")
flank_expt_cell_num <- 88000
flank_red_pos_cutoff <- 3e3
f286 <- rbind(f286_g718a_hygro[1:flank_expt_cell_num,], f286_g747a_hygro[1:flank_expt_cell_num,])
f286$ratio <- f286$BL1.A / f286$YL2.A
f286_subset <- f286 %>% filter(YL2.A >= flank_red_pos_cutoff)
f286_control_95pct_interval <- c(quantile((f286_subset %>% filter(sample == "Control"))$ratio,0.05),quantile((f286_subset %>% filter(sample == "Control"))$ratio,0.975))
f286_subset$label <- paste0(f286_subset$sample,"_",f286_subset$treatment)
#ggplot() + theme_bw() + scale_x_continuous(limits = c(-0.05, 0.5)) + geom_histogram(data = f286_subset, aes(x = ratio), binwidth = 0.01) + facet_wrap(~label, scales = "free_y") + geom_vline(xintercept  = f286_control_95pct_interval, linetype = 2)
f286_subset_summary1 <- f286_subset %>% group_by(label) %>% summarize(unexcised = sum(ratio > f286_control_95pct_interval))
```

    ## Warning in ratio > f286_control_95pct_interval: longer object length is not a
    ## multiple of shorter object length

    ## Warning in ratio > f286_control_95pct_interval: longer object length is not a
    ## multiple of shorter object length

``` r
f286_subset_summary2 <- f286_subset %>% group_by(label) %>% count()
f286_subset_summary <- merge(f286_subset_summary1, f286_subset_summary2, by = "label")
f286_subset_summary$frac_unexcised <- f286_subset_summary$unexcised / f286_subset_summary$n
f286_subset_summary$rep <- "F286"

## F286b
f286b_g718a_hygro <- read.csv(file = "data/flow/Flanking/F286/F286b_G718A_Hygro_B3.csv.gz") %>% mutate("flow" = "F286b", sample = "Flanked", treatment = "Hygro")
f286b_g747a_hygro <- read.csv(file = "data/flow/Flanking/F286/F286b_G747A_Hygro_B4.csv.gz") %>% mutate("flow" = "F286b", sample = "Control", treatment = "Hygro")
f286b_g718a_ap1903_hygro <- read.csv(file = "data/flow/Flanking/F286/F286b_G718A_AP1903_Hygro_A1.csv.gz") %>% mutate("flow" = "F286b", sample = "Flanked", treatment = "AP1903_Hygro")
f286b_g747a_ap1903_hygro <- read.csv(file = "data/flow/Flanking/F286/F286b_G747A_AP1903_Hygro_A2.csv.gz") %>% mutate("flow" = "F286b", sample = "Control", treatment = "AP1903_Hygro")
flank_expt_cell_num <- 51000
flank_red_pos_cutoff <- 3e3
f286b <- rbind(f286b_g718a_hygro[1:flank_expt_cell_num,], f286b_g747a_hygro[1:flank_expt_cell_num,], f286b_g718a_ap1903_hygro[1:flank_expt_cell_num,], f286b_g747a_ap1903_hygro[1:flank_expt_cell_num,])
f286b$ratio <- f286b$BL1.A / f286b$YL2.A
f286b_subset <- f286b %>% filter(YL2.A >= flank_red_pos_cutoff)
f286b_control_95pct_interval <- c(quantile((f286b_subset %>% filter(sample == "Control"))$ratio,0.05),quantile((f286b_subset %>% filter(sample == "Control"))$ratio,0.975))
f286b_subset$label <- paste0(f286b_subset$sample,"_",f286b_subset$treatment)
#ggplot() + theme_bw() + scale_x_continuous(limits = c(-0.05, 0.5)) + geom_histogram(data = f286b_subset, aes(x = ratio), binwidth = 0.01) + facet_wrap(~label, scales = "free_y") + geom_vline(xintercept  = f286b_control_95pct_interval, linetype = 2)
f286b_subset_summary1 <- f286b_subset %>% group_by(label) %>% summarize(unexcised = sum(ratio > f286b_control_95pct_interval))
```

    ## Warning in ratio > f286b_control_95pct_interval: longer object length is not a
    ## multiple of shorter object length

    ## Warning in ratio > f286b_control_95pct_interval: longer object length is not a
    ## multiple of shorter object length

    ## Warning in ratio > f286b_control_95pct_interval: longer object length is not a
    ## multiple of shorter object length

``` r
f286b_subset_summary2 <- f286b_subset %>% group_by(label) %>% count()
f286b_subset_summary <- merge(f286b_subset_summary1, f286b_subset_summary2, by = "label")
f286b_subset_summary$frac_unexcised <- f286b_subset_summary$unexcised / f286b_subset_summary$n
f286b_subset_summary$rep <- "F286b"


## Now combining everything
flank_ap1903_summary <- rbind(f280_subset_summary, f281_subset_summary, f281b_subset_summary, f282_subset_summary, f282b_subset_summary, f283_subset_summary, f283b_subset_summary, f284_subset_summary, f284b_subset_summary, f285_subset_summary, f286_subset_summary, f286b_subset_summary)

flank_ap1903_summary2 <- flank_ap1903_summary %>% group_by(label) %>% summarize(geomean = 10^mean(log10(frac_unexcised)))

flank_ap1903_summary$label <- factor(flank_ap1903_summary$label, levels = c("Control_none", "Control_AP1903", "Control_Hygro", "Control_AP1903_Hygro", "Flanked_none", "Flanked_AP1903", "Flanked_Hygro", "Flanked_AP1903_Hygro"))
flank_ap1903_summary2$label <- factor(flank_ap1903_summary2$label, levels = c("Control_none", "Control_AP1903", "Control_Hygro", "Control_AP1903_Hygro", "Flanked_none", "Flanked_AP1903", "Flanked_Hygro", "Flanked_AP1903_Hygro"))

Flanking_AP1903_plot <- ggplot() + theme(panel.grid.major.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  scale_y_log10(breaks = c(0.01,0.1,1)) + 
  labs(x = NULL, y = "Fraction of cells\nunexcised") + 
  geom_beeswarm(data = flank_ap1903_summary, aes(x = label, y = frac_unexcised), color = "red", alpha = 0.5, size = 1) +
  geom_point(data = flank_ap1903_summary2, aes(x = label, y = geomean), shape = 95, size = 8, color = "Black", alpha = 0.4) +
  NULL
ggsave(file = "plots/Flanking_AP1903_plot.pdf", Flanking_AP1903_plot, height = 2.25, width = 2.25)
Flanking_AP1903_plot
```

![](Recombinastics_analysis_files/figure-gfm/Effect%20of%20adding%20AP1903%20to%20flanking%20recombinations-1.png)<!-- -->

``` r
flanked_seq <- read.csv(file = "data/Flanked_sequencing.csv", header = T, stringsAsFactors = F)

flanked_seq_geomean <- flanked_seq %>% group_by(sample, result) %>% summarize(geomean_fraction = 10^(mean(log10(fraction))))
```

    ## `summarise()` has grouped output by 'sample'. You can override using the
    ## `.groups` argument.

``` r
Flanking_sequencing <- ggplot() + theme(panel.grid.major.x = element_blank(), legend.position = "top") + 
  labs(x = NULL, y = "Fraction of reads") +
  geom_point(data = flanked_seq, aes(x = sample, y = fraction, color = result), position = position_dodge(width = 0.5), alpha = 0.4) +
  geom_point(data = flanked_seq_geomean, aes(x = sample, y = geomean_fraction, color = result), position = position_dodge(width = 0.5), shape = 95, size = 10)
ggsave(file = "plots/Flanking_sequencing.pdf", Flanking_sequencing, height = 2, width = 2.2)
Flanking_sequencing
```

![](Recombinastics_analysis_files/figure-gfm/Flanking%20sequencing%20data-1.png)<!-- -->

## Creating and initially testing the double landing pad

## This is relevant to Figure 3

``` r
c2_none_raw <- read.csv(file = "data/flow/201002_F69_G783A_Clone_Comparisons/Clone2_none.csv.gz")
c2_none <- c2_none_raw[,c("FSC.A","SSC.A","VL1.A","BL1.A","YL2.A","RL1.A")]; colnames(c2_none) <- c("fsc","ssc","blu","grn","red","nir")
c2_none$clone <- "c2"; c2_none$treatment <- "none"
c2_recomb_raw <- read.csv(file = "data/flow/201002_F69_G783A_Clone_Comparisons/Clone2_recomb.csv.gz")
c2_recomb <- c2_recomb_raw[,c("FSC.A","SSC.A","VL1.A","BL1.A","YL2.A","RL1.A")]; colnames(c2_recomb) <- c("fsc","ssc","blu","grn","red","nir")
c2_recomb$clone <- "c2"; c2_recomb$treatment <- "recomb"
c2_ap1903_raw <- read.csv(file = "data/flow/201002_F69_G783A_Clone_Comparisons/Clone2_ap1903.csv.gz")
c2_ap1903 <- c2_ap1903_raw[,c("FSC.A","SSC.A","VL1.A","BL1.A","YL2.A","RL1.A")]; colnames(c2_ap1903) <- c("fsc","ssc","blu","grn","red","nir")
c2_ap1903$clone <- "c2"; c2_ap1903$treatment <- "ap1903"


cell_number <- 2000
combined_data <- rbind(c2_none[1:cell_number,], c2_recomb[1:cell_number,], c2_ap1903[1:cell_number,])
combined_data$treatment <- factor(combined_data$treatment, levels = c("recomb","none","ap1903"))

plot_alpha <- 0.2
axis_limits <- c(10,1e6)

custom_color_scale <- c("none" = "magenta", "recomb" = "black", "ap1903" = "cyan")

c2_bn_plot <- ggplot() + theme_bw() + theme(panel.grid.minor = element_blank(), legend.position = "none") + 
  scale_color_manual(values = custom_color_scale) +
  scale_x_log10(limits = axis_limits) + scale_y_log10(limits = axis_limits) +
  geom_point(data = combined_data, aes(x = blu, y = nir, color = treatment), alpha = plot_alpha)

c2_bg_plot <- ggplot() + theme_bw() + theme(panel.grid.minor = element_blank(), legend.position = "none") + 
  scale_color_manual(values = custom_color_scale) +
  scale_x_log10(limits = axis_limits) + scale_y_log10(limits = axis_limits) +
  geom_point(data = combined_data, aes(x = blu, y = grn, color = treatment), alpha = plot_alpha)

c2_br_plot <- ggplot() + theme_bw() + theme(panel.grid.minor = element_blank(), legend.position = "none") + 
  scale_color_manual(values = custom_color_scale) +
  scale_x_log10(limits = axis_limits) + scale_y_log10(limits = axis_limits) +
  geom_point(data = combined_data, aes(x = blu, y = red, color = treatment), alpha = plot_alpha)

c2_ng_plot <- ggplot() + theme_bw() + theme(panel.grid.minor = element_blank(), legend.position = "none") + 
  scale_color_manual(values = custom_color_scale) +
  scale_x_log10(limits = axis_limits) + scale_y_log10(limits = axis_limits) +
  geom_point(data = combined_data, aes(x = nir, y = grn, color = treatment), alpha = plot_alpha)

c2_nr_plot <- ggplot() + theme_bw() + theme(panel.grid.minor = element_blank(), legend.position = "none") + 
  scale_color_manual(values = custom_color_scale) +
  scale_x_log10(limits = axis_limits) + scale_y_log10(limits = axis_limits) +
  geom_point(data = combined_data, aes(x = nir, y = red, color = treatment), alpha = plot_alpha)

c2_gr_plot <- ggplot() + theme_bw() + theme(panel.grid.minor = element_blank(), legend.position = "none") + 
  scale_color_manual(values = custom_color_scale) +
  scale_x_log10(limits = axis_limits) + scale_y_log10(limits = axis_limits) +
  geom_point(data = combined_data, aes(x = grn, y = red, color = treatment), alpha = plot_alpha) +
  geom_rect(mapping = aes(xmin = 3e4, xmax = 1e6, ymin = 3e4, ymax = 1e6), alpha = 0, color = "black", linetype = 2)
c2_gr_plot
```

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Removed 2745 rows containing missing values (`geom_point()`).

![](Recombinastics_analysis_files/figure-gfm/Initial%20validation%20data%20for%20G542Ac3%20and%20G783Ac2%20cells-1.png)<!-- -->

``` r
c2_arranged_plot <- grid.arrange(c2_bn_plot, c2_bg_plot, c2_br_plot, c2_ng_plot, c2_nr_plot, c2_gr_plot, ncol=6, nrow=1)
```

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Removed 242 rows containing missing values (`geom_point()`).

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Removed 1956 rows containing missing values (`geom_point()`).

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Removed 1941 rows containing missing values (`geom_point()`).

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Removed 1981 rows containing missing values (`geom_point()`).

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Removed 1942 rows containing missing values (`geom_point()`).

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Removed 2745 rows containing missing values (`geom_point()`).

![](Recombinastics_analysis_files/figure-gfm/Initial%20validation%20data%20for%20G542Ac3%20and%20G783Ac2%20cells-2.png)<!-- -->

``` r
ggsave(file = "plots/201002/c2_arranged_plot.png", c2_arranged_plot, height = 3, width = 18)

paste("Percent double positive before selection:", round(sum(c2_recomb$grn >= 3e4 & c2_recomb$red >= 3e4) / nrow(c2_recomb) * 100,2))
```

    ## [1] "Percent double positive before selection: 3.11"

``` r
paste("Percent double positive after selection:", round(sum(c2_ap1903$grn >= 3e4 & c2_ap1903$red >= 3e4) / nrow(c2_ap1903) * 100,2))
```

    ## [1] "Percent double positive after selection: 82.98"

``` r
c6_none_raw <- read.csv(file = "data/flow/201002_F69_G783A_Clone_Comparisons/Clone6_none.csv.gz")
c6_none <- c6_none_raw[,c("FSC.A","SSC.A","VL1.A","BL1.A","YL2.A","RL1.A")]; colnames(c6_none) <- c("fsc","ssc","blu","grn","red","nir")
c6_none$clone <- "c6"; c6_none$treatment <- "none"
c6_recomb_raw <- read.csv(file = "data/flow/201002_F69_G783A_Clone_Comparisons/Clone6_recomb.csv.gz")
c6_recomb <- c6_recomb_raw[,c("FSC.A","SSC.A","VL1.A","BL1.A","YL2.A","RL1.A")]; colnames(c6_recomb) <- c("fsc","ssc","blu","grn","red","nir")
c6_recomb$clone <- "c6"; c6_recomb$treatment <- "recomb"
c6_ap1903_raw <- read.csv(file = "data/flow/201002_F69_G783A_Clone_Comparisons/Clone6_ap1903.csv.gz")
c6_ap1903 <- c6_ap1903_raw[,c("FSC.A","SSC.A","VL1.A","BL1.A","YL2.A","RL1.A")]; colnames(c6_ap1903) <- c("fsc","ssc","blu","grn","red","nir")
c6_ap1903$clone <- "c6"; c6_ap1903$treatment <- "ap1903"


cell_number <- 2000
combined_data <- rbind(c6_none[1:cell_number,], c6_recomb[1:cell_number,], c6_ap1903[1:cell_number,])
combined_data$treatment <- factor(combined_data$treatment, levels = c("recomb","none","ap1903"))

plot_alpha <- 0.2
axis_limits <- c(10,1e6)

custom_color_scale <- c("none" = "magenta", "recomb" = "black", "ap1903" = "cyan")

c6_bn_plot <- ggplot() + theme_bw() + theme(panel.grid.minor = element_blank(), legend.position = "none") + 
  scale_color_manual(values = custom_color_scale) +
  scale_x_log10(limits = axis_limits) + scale_y_log10(limits = axis_limits) +
  geom_point(data = combined_data, aes(x = blu, y = nir, color = treatment), alpha = plot_alpha)

c6_bg_plot <- ggplot() + theme_bw() + theme(panel.grid.minor = element_blank(), legend.position = "none") + 
  scale_color_manual(values = custom_color_scale) +
  scale_x_log10(limits = axis_limits) + scale_y_log10(limits = axis_limits) +
  geom_point(data = combined_data, aes(x = blu, y = grn, color = treatment), alpha = plot_alpha)

c6_br_plot <- ggplot() + theme_bw() + theme(panel.grid.minor = element_blank(), legend.position = "none") + 
  scale_color_manual(values = custom_color_scale) +
  scale_x_log10(limits = axis_limits) + scale_y_log10(limits = axis_limits) +
  geom_point(data = combined_data, aes(x = blu, y = red, color = treatment), alpha = plot_alpha)

c6_ng_plot <- ggplot() + theme_bw() + theme(panel.grid.minor = element_blank(), legend.position = "none") + 
  scale_color_manual(values = custom_color_scale) +
  scale_x_log10(limits = axis_limits) + scale_y_log10(limits = axis_limits) +
  geom_point(data = combined_data, aes(x = nir, y = grn, color = treatment), alpha = plot_alpha)

c6_nr_plot <- ggplot() + theme_bw() + theme(panel.grid.minor = element_blank(), legend.position = "none") + 
  scale_color_manual(values = custom_color_scale) +
  scale_x_log10(limits = axis_limits) + scale_y_log10(limits = axis_limits) +
  geom_point(data = combined_data, aes(x = nir, y = red, color = treatment), alpha = plot_alpha)

c6_gr_plot <- ggplot() + theme_bw() + theme(panel.grid.minor = element_blank(), legend.position = "none") + 
  scale_color_manual(values = custom_color_scale) +
  scale_x_log10(limits = axis_limits) + scale_y_log10(limits = axis_limits) +
  geom_point(data = combined_data, aes(x = grn, y = red, color = treatment), alpha = plot_alpha) +
  geom_rect(mapping = aes(xmin = 3e4, xmax = 1e6, ymin = 3e4, ymax = 1e6), alpha = 0, color = "black", linetype = 2)
c6_gr_plot
```

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Removed 3131 rows containing missing values (`geom_point()`).

![](Recombinastics_analysis_files/figure-gfm/Initial%20validation%20data%20for%20G542Ac3%20and%20G783Ac6%20cells-1.png)<!-- -->

``` r
c6_arranged_plot <- grid.arrange(c6_bn_plot, c6_bg_plot, c6_br_plot, c6_ng_plot, c6_nr_plot, c6_gr_plot, ncol=6, nrow=1)
```

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Removed 182 rows containing missing values (`geom_point()`).

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Removed 2135 rows containing missing values (`geom_point()`).

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Removed 2349 rows containing missing values (`geom_point()`).

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Removed 2132 rows containing missing values (`geom_point()`).

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Removed 2349 rows containing missing values (`geom_point()`).

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Removed 3131 rows containing missing values (`geom_point()`).

![](Recombinastics_analysis_files/figure-gfm/Initial%20validation%20data%20for%20G542Ac3%20and%20G783Ac6%20cells-2.png)<!-- -->

``` r
ggsave(file = "plots/201002/c6_arranged_plot.png", c6_arranged_plot, height = 3, width = 18)

paste("Percent double positive before selection:", round(sum(c6_recomb$grn >= 3e4 & c6_recomb$red >= 3e4) / nrow(c6_recomb) * 100,2))
```

    ## [1] "Percent double positive before selection: 3.9"

``` r
paste("Percent double positive after selection:", round(sum(c6_ap1903$grn >= 3e4 & c6_ap1903$red >= 3e4) / nrow(c6_ap1903) * 100,2))
```

    ## [1] "Percent double positive after selection: 60.28"

``` r
c11_none_raw <- read.csv(file = "data/flow/201002_F69_G783A_Clone_Comparisons/Clone11_none.csv.gz")
c11_none <- c11_none_raw[,c("FSC.A","SSC.A","VL1.A","BL1.A","YL2.A","RL1.A")]; colnames(c11_none) <- c("fsc","ssc","blu","grn","red","nir")
c11_none$clone <- "c11"; c11_none$treatment <- "none"
c11_recomb_raw <- read.csv(file = "data/flow/201002_F69_G783A_Clone_Comparisons/Clone11_recomb.csv.gz")
c11_recomb <- c11_recomb_raw[,c("FSC.A","SSC.A","VL1.A","BL1.A","YL2.A","RL1.A")]; colnames(c11_recomb) <- c("fsc","ssc","blu","grn","red","nir")
c11_recomb$clone <- "c11"; c11_recomb$treatment <- "recomb"
c11_ap1903_raw <- read.csv(file = "data/flow/201002_F69_G783A_Clone_Comparisons/Clone11_ap1903.csv.gz")
c11_ap1903 <- c11_ap1903_raw[,c("FSC.A","SSC.A","VL1.A","BL1.A","YL2.A","RL1.A")]; colnames(c11_ap1903) <- c("fsc","ssc","blu","grn","red","nir")
c11_ap1903$clone <- "c11"; c11_ap1903$treatment <- "ap1903"


cell_number <- 2000
combined_data <- rbind(c11_none[1:cell_number,], c11_recomb[1:cell_number,], c11_ap1903[1:cell_number,])
combined_data$treatment <- factor(combined_data$treatment, levels = c("recomb","none","ap1903"))

plot_alpha <- 0.2
axis_limits <- c(10,1e6)

custom_color_scale <- c("none" = "magenta", "recomb" = "black", "ap1903" = "cyan")

c11_bn_plot <- ggplot() + theme_bw() + theme(panel.grid.minor = element_blank(), legend.position = "none") + 
  scale_color_manual(values = custom_color_scale) +
  scale_x_log10(limits = axis_limits) + scale_y_log10(limits = axis_limits) +
  geom_point(data = combined_data, aes(x = blu, y = nir, color = treatment), alpha = plot_alpha)

c11_bg_plot <- ggplot() + theme_bw() + theme(panel.grid.minor = element_blank(), legend.position = "none") + 
  scale_color_manual(values = custom_color_scale) +
  scale_x_log10(limits = axis_limits) + scale_y_log10(limits = axis_limits) +
  geom_point(data = combined_data, aes(x = blu, y = grn, color = treatment), alpha = plot_alpha)

c11_br_plot <- ggplot() + theme_bw() + theme(panel.grid.minor = element_blank(), legend.position = "none") + 
  scale_color_manual(values = custom_color_scale) +
  scale_x_log10(limits = axis_limits) + scale_y_log10(limits = axis_limits) +
  geom_point(data = combined_data, aes(x = blu, y = red, color = treatment), alpha = plot_alpha)

c11_ng_plot <- ggplot() + theme_bw() + theme(panel.grid.minor = element_blank(), legend.position = "none") + 
  scale_color_manual(values = custom_color_scale) +
  scale_x_log10(limits = axis_limits) + scale_y_log10(limits = axis_limits) +
  geom_point(data = combined_data, aes(x = nir, y = grn, color = treatment), alpha = plot_alpha)

c11_nr_plot <- ggplot() + theme_bw() + theme(panel.grid.minor = element_blank(), legend.position = "none") + 
  scale_color_manual(values = custom_color_scale) +
  scale_x_log10(limits = axis_limits) + scale_y_log10(limits = axis_limits) +
  geom_point(data = combined_data, aes(x = nir, y = red, color = treatment), alpha = plot_alpha)

c11_gr_plot <- ggplot() + theme_bw() + theme(panel.grid.minor = element_blank(), legend.position = "none") + 
  scale_color_manual(values = custom_color_scale) +
  scale_x_log10(limits = axis_limits) + scale_y_log10(limits = axis_limits) +
  geom_point(data = combined_data, aes(x = grn, y = red, color = treatment), alpha = plot_alpha) +
  geom_rect(mapping = aes(xmin = 3e4, xmax = 1e6, ymin = 3e4, ymax = 1e6), alpha = 0, color = "black", linetype = 2)
c11_gr_plot
```

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Removed 2484 rows containing missing values (`geom_point()`).

![](Recombinastics_analysis_files/figure-gfm/Initial%20validation%20data%20for%20G542Ac3%20and%20G783Ac11%20cells-1.png)<!-- -->

``` r
c11_arranged_plot <- grid.arrange(c11_bn_plot, c11_bg_plot, c11_br_plot, c11_ng_plot, c11_nr_plot, c11_gr_plot, ncol=6, nrow=1)
```

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Removed 243 rows containing missing values (`geom_point()`).

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Removed 1648 rows containing missing values (`geom_point()`).

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Removed 1654 rows containing missing values (`geom_point()`).

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Removed 1708 rows containing missing values (`geom_point()`).

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Removed 1678 rows containing missing values (`geom_point()`).

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Removed 2484 rows containing missing values (`geom_point()`).

![](Recombinastics_analysis_files/figure-gfm/Initial%20validation%20data%20for%20G542Ac3%20and%20G783Ac11%20cells-2.png)<!-- -->

``` r
ggsave(file = "plots/201002/c11_arranged_plot.png", c11_arranged_plot, height = 3, width = 18)

paste("Percent double positive before selection:", round(sum(c11_recomb$grn >= 3e4 & c11_recomb$red >= 3e4) / nrow(c11_recomb) * 100,2))
```

    ## [1] "Percent double positive before selection: 2.97"

``` r
paste("Percent double positive after selection:", round(sum(c11_ap1903$grn >= 3e4 & c11_ap1903$red >= 3e4) / nrow(c11_ap1903) * 100,2))
```

    ## [1] "Percent double positive after selection: 85.21"

Below is now repeat data for testing these three clones to confirm which
is best

``` r
c2_none_raw <- read.csv(file = "data/flow/201016_F74_G783A_Clone_Comparisons/Clone2_none.csv.gz")
c2_none <- c2_none_raw[,c("FSC.A","SSC.A","VL1.A","BL1.A","YL2.A","RL1.A")]; colnames(c2_none) <- c("fsc","ssc","blu","grn","red","nir")
c2_none$clone <- "c2"; c2_none$treatment <- "none"
c2_recomb_raw <- read.csv(file = "data/flow/201016_F74_G783A_Clone_Comparisons/Clone2_recomb.csv.gz")
c2_recomb <- c2_recomb_raw[,c("FSC.A","SSC.A","VL1.A","BL1.A","YL2.A","RL1.A")]; colnames(c2_recomb) <- c("fsc","ssc","blu","grn","red","nir")
c2_recomb$clone <- "c2"; c2_recomb$treatment <- "recomb"
c2_ap1903_raw <- read.csv(file = "data/flow/201016_F74_G783A_Clone_Comparisons/Clone2_ap1903.csv.gz")
c2_ap1903 <- c2_ap1903_raw[,c("FSC.A","SSC.A","VL1.A","BL1.A","YL2.A","RL1.A")]; colnames(c2_ap1903) <- c("fsc","ssc","blu","grn","red","nir")
c2_ap1903$clone <- "c2"; c2_ap1903$treatment <- "ap1903"


cell_number <- 2000
combined_data <- rbind(c2_none[1:cell_number,], c2_recomb[1:cell_number,], c2_ap1903[1:cell_number,])
combined_data$treatment <- factor(combined_data$treatment, levels = c("recomb","none","ap1903"))

plot_alpha <- 0.2
axis_limits <- c(10,1e6)

custom_color_scale <- c("none" = "magenta", "recomb" = "black", "ap1903" = "cyan")

c2_bn_plot <- ggplot() + theme_bw() + theme(panel.grid.minor = element_blank(), legend.position = "none") + 
  scale_color_manual(values = custom_color_scale) +
  scale_x_log10(limits = axis_limits) + scale_y_log10(limits = axis_limits) +
  geom_point(data = combined_data, aes(x = blu, y = nir, color = treatment), alpha = plot_alpha)

c2_bg_plot <- ggplot() + theme_bw() + theme(panel.grid.minor = element_blank(), legend.position = "none") + 
  scale_color_manual(values = custom_color_scale) +
  scale_x_log10(limits = axis_limits) + scale_y_log10(limits = axis_limits) +
  geom_point(data = combined_data, aes(x = blu, y = grn, color = treatment), alpha = plot_alpha)

c2_br_plot <- ggplot() + theme_bw() + theme(panel.grid.minor = element_blank(), legend.position = "none") + 
  scale_color_manual(values = custom_color_scale) +
  scale_x_log10(limits = axis_limits) + scale_y_log10(limits = axis_limits) +
  geom_point(data = combined_data, aes(x = blu, y = red, color = treatment), alpha = plot_alpha)

c2_ng_plot <- ggplot() + theme_bw() + theme(panel.grid.minor = element_blank(), legend.position = "none") + 
  scale_color_manual(values = custom_color_scale) +
  scale_x_log10(limits = axis_limits) + scale_y_log10(limits = axis_limits) +
  geom_point(data = combined_data, aes(x = nir, y = grn, color = treatment), alpha = plot_alpha)

c2_nr_plot <- ggplot() + theme_bw() + theme(panel.grid.minor = element_blank(), legend.position = "none") + 
  scale_color_manual(values = custom_color_scale) +
  scale_x_log10(limits = axis_limits) + scale_y_log10(limits = axis_limits) +
  geom_point(data = combined_data, aes(x = nir, y = red, color = treatment), alpha = plot_alpha)

c2_gr_plot <- ggplot() + theme_bw() + theme(panel.grid.minor = element_blank(), legend.position = "none") + 
  scale_color_manual(values = custom_color_scale) +
  scale_x_log10(limits = axis_limits) + scale_y_log10(limits = axis_limits) +
  geom_point(data = combined_data, aes(x = grn, y = red, color = treatment), alpha = plot_alpha) +
  geom_rect(mapping = aes(xmin = 3e4, xmax = 1e6, ymin = 3e4, ymax = 1e6), alpha = 0, color = "black", linetype = 2)
c2_gr_plot
```

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Removed 3250 rows containing missing values (`geom_point()`).

![](Recombinastics_analysis_files/figure-gfm/Initial%20validation%20data%20for%20G542Ac3%20and%20G783Ac2%20cells%20repeated-1.png)<!-- -->

``` r
c2_arranged_plot <- grid.arrange(c2_bn_plot, c2_bg_plot, c2_br_plot, c2_ng_plot, c2_nr_plot, c2_gr_plot, ncol=6, nrow=1)
```

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Removed 248 rows containing missing values (`geom_point()`).

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Removed 2704 rows containing missing values (`geom_point()`).

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Removed 2317 rows containing missing values (`geom_point()`).

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Removed 2712 rows containing missing values (`geom_point()`).

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Removed 2309 rows containing missing values (`geom_point()`).

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Removed 3250 rows containing missing values (`geom_point()`).

![](Recombinastics_analysis_files/figure-gfm/Initial%20validation%20data%20for%20G542Ac3%20and%20G783Ac2%20cells%20repeated-2.png)<!-- -->

``` r
ggsave(file = "plots/201016/c2_arranged_plot.png", c2_arranged_plot, height = 3, width = 18)

paste("Percent double positive before selection:", round(sum(c2_recomb$grn >= 3e4 & c2_recomb$red >= 3e4) / nrow(c2_recomb) * 100,2))
```

    ## [1] "Percent double positive before selection: 5.01"

``` r
paste("Percent double positive after selection:", round(sum(c2_ap1903$grn >= 3e4 & c2_ap1903$red >= 3e4) / nrow(c2_ap1903) * 100,2))
```

    ## [1] "Percent double positive after selection: 78.43"

``` r
c6_none_raw <- read.csv(file = "data/flow/201016_F74_G783A_Clone_Comparisons/Clone6_none.csv.gz")
c6_none <- c6_none_raw[,c("FSC.A","SSC.A","VL1.A","BL1.A","YL2.A","RL1.A")]; colnames(c6_none) <- c("fsc","ssc","blu","grn","red","nir")
c6_none$clone <- "c6"; c6_none$treatment <- "none"
c6_recomb_raw <- read.csv(file = "data/flow/201016_F74_G783A_Clone_Comparisons/Clone6_recomb.csv.gz")
c6_recomb <- c6_recomb_raw[,c("FSC.A","SSC.A","VL1.A","BL1.A","YL2.A","RL1.A")]; colnames(c6_recomb) <- c("fsc","ssc","blu","grn","red","nir")
c6_recomb$clone <- "c6"; c6_recomb$treatment <- "recomb"
c6_ap1903_raw <- read.csv(file = "data/flow/201016_F74_G783A_Clone_Comparisons/Clone6_ap1903.csv.gz")
c6_ap1903 <- c6_ap1903_raw[,c("FSC.A","SSC.A","VL1.A","BL1.A","YL2.A","RL1.A")]; colnames(c6_ap1903) <- c("fsc","ssc","blu","grn","red","nir")
c6_ap1903$clone <- "c6"; c6_ap1903$treatment <- "ap1903"


cell_number <- 2000
combined_data <- rbind(c6_none[1:cell_number,], c6_recomb[1:cell_number,], c6_ap1903[1:cell_number,])
combined_data$treatment <- factor(combined_data$treatment, levels = c("recomb","none","ap1903"))

plot_alpha <- 0.2
axis_limits <- c(10,1e6)

custom_color_scale <- c("none" = "magenta", "recomb" = "black", "ap1903" = "cyan")

c6_bn_plot <- ggplot() + theme_bw() + theme(panel.grid.minor = element_blank(), legend.position = "none") + 
  scale_color_manual(values = custom_color_scale) +
  scale_x_log10(limits = axis_limits) + scale_y_log10(limits = axis_limits) +
  geom_point(data = combined_data, aes(x = blu, y = nir, color = treatment), alpha = plot_alpha)

c6_bg_plot <- ggplot() + theme_bw() + theme(panel.grid.minor = element_blank(), legend.position = "none") + 
  scale_color_manual(values = custom_color_scale) +
  scale_x_log10(limits = axis_limits) + scale_y_log10(limits = axis_limits) +
  geom_point(data = combined_data, aes(x = blu, y = grn, color = treatment), alpha = plot_alpha)

c6_br_plot <- ggplot() + theme_bw() + theme(panel.grid.minor = element_blank(), legend.position = "none") + 
  scale_color_manual(values = custom_color_scale) +
  scale_x_log10(limits = axis_limits) + scale_y_log10(limits = axis_limits) +
  geom_point(data = combined_data, aes(x = blu, y = red, color = treatment), alpha = plot_alpha)

c6_ng_plot <- ggplot() + theme_bw() + theme(panel.grid.minor = element_blank(), legend.position = "none") + 
  scale_color_manual(values = custom_color_scale) +
  scale_x_log10(limits = axis_limits) + scale_y_log10(limits = axis_limits) +
  geom_point(data = combined_data, aes(x = nir, y = grn, color = treatment), alpha = plot_alpha)

c6_nr_plot <- ggplot() + theme_bw() + theme(panel.grid.minor = element_blank(), legend.position = "none") + 
  scale_color_manual(values = custom_color_scale) +
  scale_x_log10(limits = axis_limits) + scale_y_log10(limits = axis_limits) +
  geom_point(data = combined_data, aes(x = nir, y = red, color = treatment), alpha = plot_alpha)

c6_gr_plot <- ggplot() + theme_bw() + theme(panel.grid.minor = element_blank(), legend.position = "none") + 
  scale_color_manual(values = custom_color_scale) +
  scale_x_log10(limits = axis_limits) + scale_y_log10(limits = axis_limits) +
  geom_point(data = combined_data, aes(x = grn, y = red, color = treatment), alpha = plot_alpha) +
  geom_rect(mapping = aes(xmin = 3e4, xmax = 1e6, ymin = 3e4, ymax = 1e6), alpha = 0, color = "black", linetype = 2)
c6_gr_plot
```

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Removed 3455 rows containing missing values (`geom_point()`).

![](Recombinastics_analysis_files/figure-gfm/Initial%20validation%20data%20for%20G542Ac3%20and%20G783Ac6%20cells%20repeated-1.png)<!-- -->

``` r
c6_arranged_plot <- grid.arrange(c6_bn_plot, c6_bg_plot, c6_br_plot, c6_ng_plot, c6_nr_plot, c6_gr_plot, ncol=6, nrow=1)
```

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Removed 235 rows containing missing values (`geom_point()`).

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Removed 2395 rows containing missing values (`geom_point()`).

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Removed 2666 rows containing missing values (`geom_point()`).

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Removed 2436 rows containing missing values (`geom_point()`).

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Removed 2708 rows containing missing values (`geom_point()`).

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Removed 3455 rows containing missing values (`geom_point()`).

![](Recombinastics_analysis_files/figure-gfm/Initial%20validation%20data%20for%20G542Ac3%20and%20G783Ac6%20cells%20repeated-2.png)<!-- -->

``` r
ggsave(file = "plots/201016/c6_arranged_plot.png", c6_arranged_plot, height = 3, width = 18)

paste("Percent double positive before selection:", round(sum(c6_recomb$grn >= 3e4 & c6_recomb$red >= 3e4) / nrow(c6_recomb) * 100,2))
```

    ## [1] "Percent double positive before selection: 5.75"

``` r
paste("Percent double positive after selection:", round(sum(c6_ap1903$grn >= 3e4 & c6_ap1903$red >= 3e4) / nrow(c6_ap1903) * 100,2))
```

    ## [1] "Percent double positive after selection: 52.41"

``` r
c11_none_raw <- read.csv(file = "data/flow/201016_F74_G783A_Clone_Comparisons/Clone11_none.csv.gz")
c11_none <- c11_none_raw[,c("FSC.A","SSC.A","VL1.A","BL1.A","YL2.A","RL1.A")]; colnames(c11_none) <- c("fsc","ssc","blu","grn","red","nir")
c11_none$clone <- "c11"; c11_none$treatment <- "none"
c11_recomb_raw <- read.csv(file = "data/flow/201016_F74_G783A_Clone_Comparisons/Clone11_recomb.csv.gz")
c11_recomb <- c11_recomb_raw[,c("FSC.A","SSC.A","VL1.A","BL1.A","YL2.A","RL1.A")]; colnames(c11_recomb) <- c("fsc","ssc","blu","grn","red","nir")
c11_recomb$clone <- "c11"; c11_recomb$treatment <- "recomb"
c11_ap1903_raw <- read.csv(file = "data/flow/201016_F74_G783A_Clone_Comparisons/Clone11_ap1903.csv.gz")
c11_ap1903 <- c11_ap1903_raw[,c("FSC.A","SSC.A","VL1.A","BL1.A","YL2.A","RL1.A")]; colnames(c11_ap1903) <- c("fsc","ssc","blu","grn","red","nir")
c11_ap1903$clone <- "c11"; c11_ap1903$treatment <- "ap1903"


cell_number <- 2000
combined_data <- rbind(c11_none[1:cell_number,], c11_recomb[1:cell_number,], c11_ap1903[1:cell_number,])
combined_data$treatment <- factor(combined_data$treatment, levels = c("recomb","none","ap1903"))

plot_alpha <- 0.2
axis_limits <- c(10,1e6)

custom_color_scale <- c("none" = "magenta", "recomb" = "black", "ap1903" = "cyan")

c11_bn_plot <- ggplot() + theme_bw() + theme(panel.grid.minor = element_blank(), legend.position = "none") + 
  scale_color_manual(values = custom_color_scale) +
  scale_x_log10(limits = axis_limits) + scale_y_log10(limits = axis_limits) +
  geom_point(data = combined_data, aes(x = blu, y = nir, color = treatment), alpha = plot_alpha)

c11_bg_plot <- ggplot() + theme_bw() + theme(panel.grid.minor = element_blank(), legend.position = "none") + 
  scale_color_manual(values = custom_color_scale) +
  scale_x_log10(limits = axis_limits) + scale_y_log10(limits = axis_limits) +
  geom_point(data = combined_data, aes(x = blu, y = grn, color = treatment), alpha = plot_alpha)

c11_br_plot <- ggplot() + theme_bw() + theme(panel.grid.minor = element_blank(), legend.position = "none") + 
  scale_color_manual(values = custom_color_scale) +
  scale_x_log10(limits = axis_limits) + scale_y_log10(limits = axis_limits) +
  geom_point(data = combined_data, aes(x = blu, y = red, color = treatment), alpha = plot_alpha)

c11_ng_plot <- ggplot() + theme_bw() + theme(panel.grid.minor = element_blank(), legend.position = "none") + 
  scale_color_manual(values = custom_color_scale) +
  scale_x_log10(limits = axis_limits) + scale_y_log10(limits = axis_limits) +
  geom_point(data = combined_data, aes(x = nir, y = grn, color = treatment), alpha = plot_alpha)

c11_nr_plot <- ggplot() + theme_bw() + theme(panel.grid.minor = element_blank(), legend.position = "none") + 
  scale_color_manual(values = custom_color_scale) +
  scale_x_log10(limits = axis_limits) + scale_y_log10(limits = axis_limits) +
  geom_point(data = combined_data, aes(x = nir, y = red, color = treatment), alpha = plot_alpha)

c11_gr_plot <- ggplot() + theme_bw() + theme(panel.grid.minor = element_blank(), legend.position = "none") + 
  scale_color_manual(values = custom_color_scale) +
  scale_x_log10(limits = axis_limits) + scale_y_log10(limits = axis_limits) +
  geom_point(data = combined_data, aes(x = grn, y = red, color = treatment), alpha = plot_alpha) +
  geom_rect(mapping = aes(xmin = 3e4, xmax = 1e6, ymin = 3e4, ymax = 1e6), alpha = 0, color = "black", linetype = 2)
c11_gr_plot
```

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Removed 2731 rows containing missing values (`geom_point()`).

![](Recombinastics_analysis_files/figure-gfm/Initial%20validation%20data%20for%20G542Ac3%20and%20G783Ac11%20cells%20repeated-1.png)<!-- -->

``` r
c11_arranged_plot <- grid.arrange(c11_bn_plot, c11_bg_plot, c11_br_plot, c11_ng_plot, c11_nr_plot, c11_gr_plot, ncol=6, nrow=1)
```

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Removed 253 rows containing missing values (`geom_point()`).

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Removed 1829 rows containing missing values (`geom_point()`).

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Removed 1920 rows containing missing values (`geom_point()`).

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Removed 1857 rows containing missing values (`geom_point()`).

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Removed 1920 rows containing missing values (`geom_point()`).

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Removed 2731 rows containing missing values (`geom_point()`).

![](Recombinastics_analysis_files/figure-gfm/Initial%20validation%20data%20for%20G542Ac3%20and%20G783Ac11%20cells%20repeated-2.png)<!-- -->

``` r
ggsave(file = "plots/201016/c11_arranged_plot.png", c11_arranged_plot, height = 3, width = 18)

paste("Percent double positive before selection:", round(sum(c11_recomb$grn >= 3e4 & c11_recomb$red >= 3e4) / nrow(c11_recomb) * 100,2))
```

    ## [1] "Percent double positive before selection: 3.83"

``` r
paste("Percent double positive after selection:", round(sum(c11_ap1903$grn >= 3e4 & c11_ap1903$red >= 3e4) / nrow(c11_ap1903) * 100,2))
```

    ## [1] "Percent double positive after selection: 73.16"

## We’re moving forward with clone 11, so show that data here

``` r
plot_alpha <- 0.2
axis_limits <- c(10,1e6)

c11_none_bn_plot <- ggplot() + theme_bw() + theme(panel.grid.minor = element_blank(), legend.position = "none") + 
  scale_x_log10(limits = axis_limits) + scale_y_log10(limits = axis_limits) +
  geom_point(data = combined_data %>% filter(treatment == "none"), aes(x = blu, y = nir), alpha = plot_alpha)
c11_none_bn_plot
```

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Removed 63 rows containing missing values (`geom_point()`).

![](Recombinastics_analysis_files/figure-gfm/G542Ac3%20and%20G783Ac11%20cell%20plots%20for%20manuscript-1.png)<!-- -->

``` r
c11_recomb_bn_plot <- ggplot() + theme_bw() + theme(panel.grid.minor = element_blank(), legend.position = "recomb") + 
  scale_x_log10(limits = axis_limits) + scale_y_log10(limits = axis_limits) +
  geom_point(data = combined_data %>% filter(treatment == "recomb"), aes(x = blu, y = nir), alpha = plot_alpha)
c11_recomb_bn_plot
```

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Removed 43 rows containing missing values (`geom_point()`).

![](Recombinastics_analysis_files/figure-gfm/G542Ac3%20and%20G783Ac11%20cell%20plots%20for%20manuscript-2.png)<!-- -->

``` r
c11_ap1903_bn_plot <- ggplot() + theme_bw() + theme(panel.grid.minor = element_blank(), legend.position = "ap1903") + 
  scale_x_log10(limits = axis_limits) + scale_y_log10(limits = axis_limits) +
  geom_point(data = combined_data %>% filter(treatment == "ap1903"), aes(x = blu, y = nir), alpha = plot_alpha)
c11_ap1903_bn_plot
```

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Removed 147 rows containing missing values (`geom_point()`).

![](Recombinastics_analysis_files/figure-gfm/G542Ac3%20and%20G783Ac11%20cell%20plots%20for%20manuscript-3.png)<!-- -->

``` r
c11_none_gr_plot <- ggplot() + theme_bw() + theme(panel.grid.minor = element_blank(), legend.position = "none") + 
  scale_x_log10(limits = axis_limits) + scale_y_log10(limits = axis_limits) +
  geom_point(data = combined_data %>% filter(treatment == "none"), aes(x = grn, y = red), alpha = plot_alpha)
c11_none_gr_plot
```

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Removed 1142 rows containing missing values (`geom_point()`).

![](Recombinastics_analysis_files/figure-gfm/G542Ac3%20and%20G783Ac11%20cell%20plots%20for%20manuscript-4.png)<!-- -->

``` r
c11_recomb_gr_plot <- ggplot() + theme_bw() + theme(panel.grid.minor = element_blank(), legend.position = "recomb") + 
  scale_x_log10(limits = axis_limits) + scale_y_log10(limits = axis_limits) +
  geom_point(data = combined_data %>% filter(treatment == "recomb"), aes(x = grn, y = red), alpha = plot_alpha)
c11_recomb_gr_plot
```

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Removed 1281 rows containing missing values (`geom_point()`).

![](Recombinastics_analysis_files/figure-gfm/G542Ac3%20and%20G783Ac11%20cell%20plots%20for%20manuscript-5.png)<!-- -->

``` r
c11_ap1903_gr_plot <- ggplot() + theme_bw() + theme(panel.grid.minor = element_blank(), legend.position = "ap1903") + 
  scale_x_log10(limits = axis_limits) + scale_y_log10(limits = axis_limits) +
  geom_point(data = combined_data %>% filter(treatment == "ap1903"), aes(x = grn, y = red), alpha = plot_alpha)
c11_ap1903_gr_plot
```

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Removed 308 rows containing missing values (`geom_point()`).

![](Recombinastics_analysis_files/figure-gfm/G542Ac3%20and%20G783Ac11%20cell%20plots%20for%20manuscript-6.png)<!-- -->

``` r
c11_arranged_plot <- grid.arrange(c11_bn_plot, c11_bg_plot, c11_br_plot, c11_ng_plot, c11_nr_plot, c11_gr_plot, ncol=6, nrow=1)
```

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Removed 253 rows containing missing values (`geom_point()`).

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Removed 1829 rows containing missing values (`geom_point()`).

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Removed 1920 rows containing missing values (`geom_point()`).

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Removed 1857 rows containing missing values (`geom_point()`).

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Removed 1920 rows containing missing values (`geom_point()`).

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Removed 2731 rows containing missing values (`geom_point()`).

![](Recombinastics_analysis_files/figure-gfm/G542Ac3%20and%20G783Ac11%20cell%20plots%20for%20manuscript-7.png)<!-- -->

``` r
ggsave(file = "plots/201016/c11_arranged_plot.png", c11_arranged_plot, height = 3, width = 18)

paste("Percent double positive before selection:", round(sum(c11_recomb$grn >= 3e4 & c11_recomb$red >= 3e4) / nrow(c11_recomb) * 100,2))
```

    ## [1] "Percent double positive before selection: 3.83"

``` r
paste("Percent double positive after selection:", round(sum(c11_ap1903$grn >= 3e4 & c11_ap1903$red >= 3e4) / nrow(c11_ap1903) * 100,2))
```

    ## [1] "Percent double positive after selection: 73.16"

Plotting the double landing pad PCA data for all relevant colors

``` r
# http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/118-principal-component-analysis-in-r-prcomp-vs-princomp/

combined_data <- rbind(c11_none[1:2000,], c11_recomb[1:2000,], c11_ap1903[1:2000,])

ggplot() + theme_bw() + scale_x_log10() + scale_y_log10() + geom_point(data = combined_data, aes(x = red, y = nir))
```

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Removed 1808 rows containing missing values (`geom_point()`).

![](Recombinastics_analysis_files/figure-gfm/Principal%20component%20analysis-1.png)<!-- -->

``` r
combined_data$zir <- combined_data$nir - combined_data$red * 0.01
ggplot() + theme_bw() + scale_x_log10() + scale_y_log10() + geom_point(data = combined_data, aes(x = red, y = zir))
```

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Removed 3467 rows containing missing values (`geom_point()`).

![](Recombinastics_analysis_files/figure-gfm/Principal%20component%20analysis-2.png)<!-- -->

``` r
combined_for_pca <- combined_data[,c("blu","grn","red","zir")]

clone2.pca <- prcomp(combined_for_pca, center = TRUE,scale. = TRUE)
fviz_eig(clone2.pca)
```

![](Recombinastics_analysis_files/figure-gfm/Principal%20component%20analysis-3.png)<!-- -->

``` r
autoplot(clone2.pca, alpha = 0.2) +
  scale_x_continuous(limits = c(-0.3,0.3)) + scale_y_continuous(limits = c(-0.3,0.3))
```

![](Recombinastics_analysis_files/figure-gfm/Principal%20component%20analysis-4.png)<!-- -->

``` r
fviz_pca_var(clone2.pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
             )
```

![](Recombinastics_analysis_files/figure-gfm/Principal%20component%20analysis-5.png)<!-- -->

``` r
clone2.pca$rotation
```

    ##            PC1        PC2        PC3        PC4
    ## blu  0.4694194 -0.6135853  0.5941415 -0.2239517
    ## grn -0.5248392 -0.4194442 -0.2773776 -0.6867839
    ## red -0.5118940 -0.5278099  0.1113592  0.6685659
    ## zir  0.4920912 -0.4110917 -0.7467641  0.1766158

``` r
pca2 <- data.frame("pc1" = clone2.pca$x[,1], "pc2" = clone2.pca$x[,2])

pca2$treatment <- combined_data[,"treatment"]

pca2$treatment <- factor(pca2$treatment, levels = c("none","recomb","ap1903"))
custom_color_scale <- c("none" = "magenta", "ap1903" = "cyan", "recomb" = "green")

ggplot() + theme_bw() + theme(panel.grid.minor = element_blank()) + 
  scale_color_manual(values = custom_color_scale) +
  scale_x_continuous(limits = c(-5,5)) +
  scale_y_continuous(limits = c(-3,2)) +
  geom_point(data = pca2, aes(x = pc1, y = pc2, color = treatment), alpha = 0.05)
```

    ## Warning: Removed 18 rows containing missing values (`geom_point()`).

![](Recombinastics_analysis_files/figure-gfm/Principal%20component%20analysis-6.png)<!-- -->

``` r
ggplot() + theme_bw() + theme(panel.grid.minor = element_blank()) + 
  scale_color_manual(values = custom_color_scale) +
  scale_x_continuous(limits = c(-5,5)) +
  scale_y_continuous(limits = c(-3,2)) +
  geom_point(data = pca2, aes(x = pc1, y = pc2), alpha = 0.05) +
  facet_grid(rows = vars(treatment))
```

    ## Warning: Removed 18 rows containing missing values (`geom_point()`).

![](Recombinastics_analysis_files/figure-gfm/Principal%20component%20analysis-7.png)<!-- -->

``` r
pca_details <- cbind(combined_for_pca,pca2)

pca_details_melted <- melt(pca_details[,c("treatment","blu","grn","red","zir")], id. = "treatment")

ggplot() + theme_bw() +
  scale_x_log10() + 
  geom_histogram(data = pca_details_melted, aes(x = value)) +
  facet_grid(cols = vars(variable))
```

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

    ## Warning: Removed 5287 rows containing non-finite values (`stat_bin()`).

![](Recombinastics_analysis_files/figure-gfm/Principal%20component%20analysis-8.png)<!-- -->

``` r
blu_cutoff <- 1e3
grn_cutoff <- 1e4
red_cutoff <- 1e4
zir_cutoff <- 1e3

ggplot() + theme_bw() + scale_x_log10() + geom_histogram(data = pca_details, aes(x = blu)) + geom_vline(xintercept = blu_cutoff)
```

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning in self$trans$transform(x): Transformation introduced infinite values in
    ## continuous x-axis

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

    ## Warning: Removed 99 rows containing non-finite values (`stat_bin()`).

![](Recombinastics_analysis_files/figure-gfm/Principal%20component%20analysis-9.png)<!-- -->

``` r
ggplot() + theme_bw() + scale_x_log10() + geom_histogram(data = pca_details, aes(x = grn)) + geom_vline(xintercept = grn_cutoff)
```

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning in self$trans$transform(x): Transformation introduced infinite values in
    ## continuous x-axis

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

    ## Warning: Removed 1651 rows containing non-finite values (`stat_bin()`).

![](Recombinastics_analysis_files/figure-gfm/Principal%20component%20analysis-10.png)<!-- -->

``` r
ggplot() + theme_bw() + scale_x_log10() + geom_histogram(data = pca_details, aes(x = red)) + geom_vline(xintercept = red_cutoff)
```

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning in self$trans$transform(x): Transformation introduced infinite values in
    ## continuous x-axis

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

    ## Warning: Removed 1743 rows containing non-finite values (`stat_bin()`).

![](Recombinastics_analysis_files/figure-gfm/Principal%20component%20analysis-11.png)<!-- -->

``` r
ggplot() + theme_bw() + scale_x_log10() + geom_histogram(data = pca_details, aes(x = zir)) + geom_vline(xintercept = zir_cutoff)
```

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning in self$trans$transform(x): Transformation introduced infinite values in
    ## continuous x-axis

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

    ## Warning: Removed 1794 rows containing non-finite values (`stat_bin()`).

![](Recombinastics_analysis_files/figure-gfm/Principal%20component%20analysis-12.png)<!-- -->

``` r
pca_details$blu2 <- 0
pca_details$grn2 <- 0
pca_details$red2 <- 0
pca_details$zir2 <- 0

for(x in 1:nrow(pca_details)){
  if(pca_details$blu[x] > blu_cutoff){pca_details$blu2[x] <- 1}
  if(pca_details$grn[x] > grn_cutoff){pca_details$grn2[x] <- 1}
  if(pca_details$red[x] > red_cutoff){pca_details$red2[x] <- 1}
  if(pca_details$zir[x] > zir_cutoff){pca_details$zir2[x] <- 1}
}

pca_details$bgrn <- paste(pca_details$blu2,pca_details$grn2,pca_details$red2,pca_details$zir2,sep="")

bgrn_table <- data.frame(table(pca_details$bgrn)) %>% arrange(desc(Freq))

pca_details$label <- "Unannotated"

for(x in 1:nrow(pca_details)){
  if(pca_details$bgrn[x] == "1001"){pca_details$label[x] <- "Unrecombined"}
  if(pca_details$bgrn[x] == "0110"){pca_details$label[x] <- "Doubly_recombined"}
  if(pca_details$bgrn[x] == "1010"){pca_details$label[x] <- "GA_only"}
  if(pca_details$bgrn[x] == "0101"){pca_details$label[x] <- "GT_only"}
}

custom_color_scale2 <- c("Unannotated" = "black", "Unrecombined" = "blue", "Doubly_recombined" = "purple", "GA_only" = "red", "GT_only" = "green")

Clone11_PCA_plot <- ggplot() + theme_bw() + theme(panel.grid.minor = element_blank()) + 
  scale_color_manual(values = custom_color_scale2) +
  scale_x_continuous(limits = c(-5,5)) +
  scale_y_continuous(limits = c(-3,2)) +
  geom_point(data = pca_details, aes(x = pc1, y = pc2, color = label), alpha = 0.1, size = 0.5) +
  facet_grid(rows = vars(treatment))
Clone11_PCA_plot
```

    ## Warning: Removed 18 rows containing missing values (`geom_point()`).

![](Recombinastics_analysis_files/figure-gfm/Principal%20component%20analysis-13.png)<!-- -->

``` r
ggsave(file = "plots/Clone11_PCA_plot.pdf", Clone11_PCA_plot, height = 3.4, width = 4.5)
```

    ## Warning: Removed 18 rows containing missing values (`geom_point()`).

``` r
plot_alpha = 0.02
point_size = 0.25

c11_none_bz_plot <- ggplot() + theme_bw() + theme(panel.grid.minor = element_blank(), legend.position = "none") + 
  scale_x_log10(limits = axis_limits) + scale_y_log10(limits = axis_limits) +
  geom_point(data = pca_details %>% filter(treatment == "none"), aes(x = blu, y = zir), alpha = plot_alpha, size = point_size)
c11_none_bz_plot
```

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Removed 64 rows containing missing values (`geom_point()`).

![](Recombinastics_analysis_files/figure-gfm/Individual%20flow%20plots-1.png)<!-- -->

``` r
c11_recomb_bz_plot <- ggplot() + theme_bw() + theme(panel.grid.minor = element_blank(), legend.position = "none") + 
  scale_x_log10(limits = axis_limits) + scale_y_log10(limits = axis_limits) +
  geom_point(data = pca_details %>% filter(treatment == "recomb"), aes(x = blu, y = zir), alpha = plot_alpha, size = point_size)
c11_recomb_bz_plot
```

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Removed 309 rows containing missing values (`geom_point()`).

![](Recombinastics_analysis_files/figure-gfm/Individual%20flow%20plots-2.png)<!-- -->

``` r
c11_ap1903_bz_plot <- ggplot() + theme_bw() + theme(panel.grid.minor = element_blank(), legend.position = "none") + 
  scale_x_log10(limits = axis_limits) + scale_y_log10(limits = axis_limits) +
  geom_point(data = pca_details %>% filter(treatment == "ap1903"), aes(x = blu, y = zir), alpha = plot_alpha, size = point_size)
c11_ap1903_bz_plot
```

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Removed 1544 rows containing missing values (`geom_point()`).

![](Recombinastics_analysis_files/figure-gfm/Individual%20flow%20plots-3.png)<!-- -->

``` r
c11_none_gr_plot <- ggplot() + theme_bw() + theme(panel.grid.minor = element_blank(), legend.position = "none") + 
  scale_x_log10(limits = axis_limits) + scale_y_log10(limits = axis_limits) +
  geom_point(data = pca_details %>% filter(treatment == "none"), aes(x = grn, y = red), alpha = plot_alpha, size = point_size)
c11_none_gr_plot
```

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Removed 1142 rows containing missing values (`geom_point()`).

![](Recombinastics_analysis_files/figure-gfm/Individual%20flow%20plots-4.png)<!-- -->

``` r
c11_recomb_gr_plot <- ggplot() + theme_bw() + theme(panel.grid.minor = element_blank(), legend.position = "recomb") + 
  scale_x_log10(limits = axis_limits) + scale_y_log10(limits = axis_limits) +
  geom_point(data = pca_details %>% filter(treatment == "recomb"), aes(x = grn, y = red), alpha = plot_alpha, size = point_size)
c11_recomb_gr_plot
```

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Removed 1281 rows containing missing values (`geom_point()`).

![](Recombinastics_analysis_files/figure-gfm/Individual%20flow%20plots-5.png)<!-- -->

``` r
c11_ap1903_gr_plot <- ggplot() + theme_bw() + theme(panel.grid.minor = element_blank(), legend.position = "ap1903") + 
  scale_x_log10(limits = axis_limits) + scale_y_log10(limits = axis_limits) +
  geom_point(data = pca_details %>% filter(treatment == "ap1903"), aes(x = grn, y = red), alpha = plot_alpha, size = point_size)
c11_ap1903_gr_plot
```

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Removed 308 rows containing missing values (`geom_point()`).

![](Recombinastics_analysis_files/figure-gfm/Individual%20flow%20plots-6.png)<!-- -->

``` r
c11_arranged_plot <- grid.arrange(c11_none_bz_plot, c11_recomb_bz_plot, c11_ap1903_bz_plot, c11_none_gr_plot, c11_recomb_gr_plot, c11_ap1903_gr_plot, ncol=3, nrow=2)
```

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Removed 64 rows containing missing values (`geom_point()`).

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Removed 309 rows containing missing values (`geom_point()`).

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Removed 1544 rows containing missing values (`geom_point()`).

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Removed 1142 rows containing missing values (`geom_point()`).

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Removed 1281 rows containing missing values (`geom_point()`).

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Removed 308 rows containing missing values (`geom_point()`).

![](Recombinastics_analysis_files/figure-gfm/Individual%20flow%20plots-7.png)<!-- -->

``` r
ggsave(file = "plots/c11_arranged_plot.pdf", c11_arranged_plot, height = 2.2, width = 4.1)

#paste("Percent double positive before selection:", round(sum(c11_recomb$grn >= 3e4 & c11_recomb$red >= 3e4) / nrow(c11_recomb) * 100,2))
#paste("Percent double positive after selection:", round(sum(c11_ap1903$grn >= 3e4 & c11_ap1903$red >= 3e4) / nrow(c11_ap1903) * 100,2))

pca_vectorplot <- fviz_pca_var(clone2.pca,
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
             )
pca_vectorplot
```

![](Recombinastics_analysis_files/figure-gfm/Individual%20flow%20plots-8.png)<!-- -->

``` r
ggsave(file = "plots/pca_vectorplot.pdf", pca_vectorplot, height = 2, width = 2)
```

## Double landing pad cells for biosensors

## GCaMP and STIM1 analysis relevant to figure 4

``` r
rep1 = read.csv("Data/jgcamp_data/rep1_signal.csv")
rep1 = rep1[-1]
rep2 = read.csv("Data/jgcamp_data/rep2_signal.csv")
rep2 = rep2[-1]

## Shifting times so addition of drug is time zero
rep1$ND.T <- rep1$ND.T - 30
rep2$ND.T <- rep2$ND.T - 30

#separate based on stimulus and variant
rep1_mch_cch = rep1 %>% filter(Name == "mCherry_cch")
rep1_005_cch = rep1 %>% filter(Name == "pNK005_cch")
rep1_013_cch = rep1 %>% filter(Name == "pNK013_cch")

rep1_mch_iono = rep1 %>% filter(Name == "mCherry_ionomycin")
rep1_005_iono = rep1 %>% filter(Name == "pNK005_ionomycin")
rep1_013_iono = rep1 %>% filter(Name == "pNK013_ionomycin")

rep2_mch_cch = rep2 %>% filter(Name == "mCherry_cch")
rep2_005_cch = rep2 %>% filter(Name == "pNK005_cch")
rep2_013_cch = rep2 %>% filter(Name == "pNK013_cch")

rep2_mch_iono = rep2 %>% filter(Name == "mCherry_ionomycin")
rep2_005_iono = rep2 %>% filter(Name == "pNK005_ionomycin")
rep2_013_iono = rep2 %>% filter(Name == "pNK013_ionomycin")
```

Figure 4B - Representative plot of MFI vs time for fluorescence
microscopy measurements - add carbachol and use neg control

``` r
rep_4b = ggplot() +
  scale_x_continuous(limits = c(-30,150), breaks = c(-30,0,30,60,90,120,150,180)) + 
  scale_y_continuous(limits = c(0,2500))+
  geom_point(data = rep1_mch_cch, aes(x = ND.T, y = signal, color = "100uM Carbachol"), size =3, alpha = 0.2, color = "black")+
  geom_point(data = rep1_mch_iono, aes(x = ND.T, y = signal, color = "2uM Ionomycin"), size =3, alpha = 0.2, color = "orange")+
  labs(x = "Time (seconds)", y = "Green MFI", color = "Condition") +
  theme(legend.position = "top", panel.grid.major = element_blank())
rep_4b
```

    ## Warning: Removed 1 rows containing missing values (`geom_point()`).

    ## Warning: Removed 2 rows containing missing values (`geom_point()`).

![](Recombinastics_analysis_files/figure-gfm/Representative%20GcaMP%20plot%20with%20negative%20control-1.png)<!-- -->

``` r
ggsave("Plots/rep_4b.pdf", rep_4b, width = 3.5, height = 1.75, useDingbats = F)
```

    ## Warning: Removed 1 rows containing missing values (`geom_point()`).
    ## Removed 2 rows containing missing values (`geom_point()`).

Supplementary Figure 3 - un-normalized fluorescence measurement curves
for variants in carbachol and ionomycin

``` r
#fluorescence measurements for cch of each replicate
cch_rep1 = bind_rows(rep1_mch_cch, rep1_005_cch, rep1_013_cch) 
cch_rep1 = cch_rep1 %>% mutate(rep = '1')
cch_rep2 = bind_rows(rep2_mch_cch, rep2_005_cch, rep2_013_cch)
cch_rep2 = cch_rep2 %>% mutate(rep = '2')
cch_reps = bind_rows(cch_rep1, cch_rep2)
cch_reps[cch_reps == "mCherry_cch"] = "Negative control"
cch_reps[cch_reps == "pNK005_cch"] = "WT STIM1"
cch_reps[cch_reps == "pNK013_cch"] = "R429C STIM1"

cch_reps$Name <- factor(cch_reps$Name, levels = c("Negative control", "WT STIM1", "R429C STIM1"))

b = ggplot(data = cch_reps, aes(x = ND.T, y = signal, color = Name))+ theme(panel.grid.major = element_blank(), legend.position = "top") +
  scale_y_continuous(limits = c(0,2500))+ 
  scale_x_continuous(limits= c(-30,150), breaks = c(0,60,120))+
  labs(title = "100uM Carbachol", x = "Time (s)", y = "Green MFI", color = "Variant")+ 
  geom_point(size = 1, alpha = 0.2) +
  facet_wrap(~rep)
b
```

    ## Warning: Removed 7 rows containing missing values (`geom_point()`).

![](Recombinastics_analysis_files/figure-gfm/Raw%20GCaMP%20curves%20for%20supplement-1.png)<!-- -->

``` r
ggsave("Plots/cch_reps.pdf", b, width = 7, height = 4)
```

    ## Warning: Removed 7 rows containing missing values (`geom_point()`).

``` r
#fluorescence measurements for ionomycin of each replicate
iono_rep1 = bind_rows(rep1_mch_iono, rep1_005_iono, rep1_013_iono) 
iono_rep1 = iono_rep1 %>% mutate(rep = '1')
iono_rep2 = bind_rows(rep2_mch_iono, rep2_005_iono, rep2_013_iono)
iono_rep2 = iono_rep2 %>% mutate(rep = '2')
iono_reps = bind_rows(iono_rep1, iono_rep2)
iono_reps[iono_reps == "mCherry_ionomycin"] = "Negative control"
iono_reps[iono_reps == "pNK005_ionomycin"] = "WT STIM1"
iono_reps[iono_reps == "pNK013_ionomycin"] = "R429C STIM1"

iono_reps$Name <- factor(iono_reps$Name, levels = c("Negative control", "WT STIM1", "R429C STIM1"))

c = ggplot(data = iono_reps, aes(x = ND.T, y = signal, color = Name))+ theme(panel.grid.major = element_blank(), legend.position = "top") +
  scale_y_continuous(limits = c(0,6000))+ scale_x_continuous(limits= c(-30,150), breaks = c(0,60,120))+ 
  labs(title = "2uM Ionomycin", x = "Time (s)", y = "Green MFI")+ 
  geom_point(size = 1, alpha = 0.2) +
  facet_wrap(~rep)
c
```

    ## Warning: Removed 9 rows containing missing values (`geom_point()`).

![](Recombinastics_analysis_files/figure-gfm/Raw%20GCaMP%20curves%20for%20supplement-2.png)<!-- -->

``` r
ggsave("Plots/iono_reps.pdf", c, width = 7, height = 4)
```

    ## Warning: Removed 9 rows containing missing values (`geom_point()`).

``` r
raw_gcamp_plots <- b | c

raw_gcamp_plots
```

    ## Warning: Removed 7 rows containing missing values (`geom_point()`).

    ## Warning: Removed 9 rows containing missing values (`geom_point()`).

![](Recombinastics_analysis_files/figure-gfm/Raw%20GCaMP%20curves%20for%20supplement-3.png)<!-- -->

``` r
ggsave(file = "Plots/Raw_gcamp_plots.pdf", raw_gcamp_plots, height = 3, width = 5)
```

    ## Warning: Removed 7 rows containing missing values (`geom_point()`).
    ## Removed 9 rows containing missing values (`geom_point()`).

Figure 4C - Differences in initial vs peak fluorescence for each variant
for carbachol and ionomycin

``` r
#100uM carbachol 
#replicate 1
comp_rep1_cch = data.frame(matrix(ncol = 3, nrow = 3))
colnames(comp_rep1_cch) = c("variant", "Start value", "Peak value")
comp_rep1_cch[1,1] = "Negative control"
comp_rep1_cch[2,1] = "WT STIM1"
comp_rep1_cch[3,1] = "R429C STIM1"
comp_rep1_cch[1,2] = mean(rep1_mch_cch$signal[1:30]) #min(rep1_mch_cch$signal)
comp_rep1_cch[2,2] = mean(rep1_005_cch$signal[1:30]) #min(rep1_005_cch$signal)
comp_rep1_cch[3,2] = mean(rep1_013_cch$signal[1:30]) #min(rep1_013_cch$signal)
comp_rep1_cch[1,3] = quantile(rep1_mch_cch$signal, 0.95) #max(rep1_mch_cch$signal)
comp_rep1_cch[2,3] = quantile(rep1_005_cch$signal, 0.95) #max(rep1_005_cch$signal)
comp_rep1_cch[3,3] = quantile(rep1_013_cch$signal, 0.95) #max(rep1_013_cch$signal)
comp_rep1_cchlg = comp_rep1_cch %>% melt() %>% mutate(replicate = '1')
```

    ## Using variant as id variables

``` r
#replicate 2
comp_rep2_cch = data.frame(matrix(ncol = 3, nrow = 3))
colnames(comp_rep2_cch) = c("variant", "Start value", "Peak value")
comp_rep2_cch[1,1] = "Negative control"
comp_rep2_cch[2,1] = "WT STIM1"
comp_rep2_cch[3,1] = "R429C STIM1"
comp_rep2_cch[1,2] = mean(rep2_mch_cch$signal[1:30]) #min(rep2_mch_cch$signal)
comp_rep2_cch[2,2] = mean(rep2_005_cch$signal[1:30]) #min(rep2_005_cch$signal)
comp_rep2_cch[3,2] = mean(rep2_013_cch$signal[1:30]) #min(rep2_013_cch$signal)
comp_rep2_cch[1,3] = quantile(rep2_mch_cch$signal, 0.95) #max(rep2_mch_cch$signal)
comp_rep2_cch[2,3] = quantile(rep2_005_cch$signal, 0.95) #max(rep2_005_cch$signal)
comp_rep2_cch[3,3] = quantile(rep2_013_cch$signal, 0.95) #max(rep2_013_cch$signal)
comp_rep2_cchlg = comp_rep2_cch %>% melt() %>% mutate(replicate = '2')
```

    ## Using variant as id variables

``` r
comp_cch = bind_rows(comp_rep1_cchlg, comp_rep2_cchlg)

cch_4c = ggplot(data = comp_cch, aes(x =variable, y= value, color = variant))+ 
  geom_point(position=position_dodge(width = 0.25)) + 
  stat_summary(geom = "point", shape = 8,  position = position_dodge(width = 0.25))+
  scale_x_discrete(labels = c("Initial", "Peak"))+ 
  labs(x = "Time point", y = "Green MFI", color = "Variant")+ 
  theme(text=element_text(size=20), legend.position = "none")
cch_4c
```

    ## No summary function supplied, defaulting to `mean_se()`

![](Recombinastics_analysis_files/figure-gfm/GCaMP%20initial%20and%20peak%20comparisons-1.png)<!-- -->

``` r
ggsave("Plots/cch_4c.pdf", cch_4c, width = 4, height = 3)
```

    ## No summary function supplied, defaulting to `mean_se()`

``` r
#ionomycin
#replicate 1
comp_rep1_iono = data.frame(matrix(ncol = 3, nrow = 3))
colnames(comp_rep1_iono) = c("variant", "Start value", "Peak value")
comp_rep1_iono[1,1] = "Negative control"
comp_rep1_iono[2,1] = "WT STIM1"
comp_rep1_iono[3,1] = "R429C STIM1"
comp_rep1_iono[1,2] = mean(rep1_mch_iono$signal[1:30]) #min(rep1_mch_iono$signal)
comp_rep1_iono[2,2] = mean(rep1_005_iono$signal[1:30]) #min(rep1_005_iono$signal)
comp_rep1_iono[3,2] = mean(rep1_013_iono$signal[1:30]) #min(rep1_013_iono$signal)
comp_rep1_iono[1,3] =  #max(rep1_mch_iono$signal)
comp_rep1_iono[2,3] = quantile(rep1_005_iono$signal, 0.95) #max(rep1_005_iono$signal)
comp_rep1_iono[3,3] = quantile(rep1_013_iono$signal, 0.95) #max(rep1_013_iono$signal)
comp_rep1_ionolg = comp_rep1_iono %>% melt() %>% mutate(replicate = '1')
```

    ## Using variant as id variables

``` r
#replicate 2
comp_rep2_iono = data.frame(matrix(ncol = 3, nrow = 3))
colnames(comp_rep2_iono) = c("variant", "Start value", "Peak value")
comp_rep2_iono[1,1] = "Negative control"
comp_rep2_iono[2,1] = "WT STIM1"
comp_rep2_iono[3,1] = "R429C STIM1"
comp_rep2_iono[1,2] = mean(rep2_mch_iono$signal[1:30]) #min(rep2_mch_iono$signal)
comp_rep2_iono[2,2] = mean(rep2_005_iono$signal[1:30]) #min(rep2_005_iono$signal)
comp_rep2_iono[3,2] = mean(rep2_013_iono$signal[1:30]) #min(rep2_013_iono$signal)
comp_rep2_iono[1,3] = quantile(rep2_mch_iono$signal, 0.95) #max(rep2_mch_iono$signal)
comp_rep2_iono[2,3] = quantile(rep2_005_iono$signal, 0.95) #max(rep2_005_iono$signal)
comp_rep2_iono[3,3] = quantile(rep2_013_iono$signal, 0.95) #max(rep2_013_iono$signal)
comp_rep2_ionolg = comp_rep2_iono %>% melt() %>% mutate(replicate = '2')
```

    ## Using variant as id variables

``` r
comp_iono = bind_rows(comp_rep1_ionolg, comp_rep2_ionolg)

iono_4c = ggplot(data = comp_iono, aes(x =variable, y= value, color = variant))+ scale_y_continuous(limits = c(0,6000))+
 geom_point(position=position_dodge(width = 0.25)) + 
  stat_summary(geom = "point", shape = 8,  position = position_dodge(width = 0.25))+
  scale_x_discrete(labels = c("Initial", "Peak"))+ 
  labs(x = "Time point", y = "Green MFI", color = "Variant")+ 
  theme(text=element_text(size=20), legend.position = "none")
iono_4c
```

    ## No summary function supplied, defaulting to `mean_se()`

![](Recombinastics_analysis_files/figure-gfm/GCaMP%20initial%20and%20peak%20comparisons-2.png)<!-- -->

``` r
ggsave("Plots/iono_4c.pdf", iono_4c, width = 4, height = 3)
```

    ## No summary function supplied, defaulting to `mean_se()`

``` r
## Another version of the plots, from the way I do it (KAM)

comp_cch$variant <- factor(comp_cch$variant, levels = c("Negative control", "WT STIM1", "R429C STIM1"))

comp_iono$variant <- factor(comp_iono$variant, levels = c("Negative control", "WT STIM1", "R429C STIM1"))

cch4c_version2 <- ggplot() + theme(legend.position = "top", panel.grid.major.x = element_blank(), axis.text.x.bottom = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  labs(title = "Carbachol", x = NULL, y = "Mean fluorescent intensity") +
  scale_y_continuous(limits = c(-200,6000)) +
  geom_hline(yintercept = 0, alpha = 0.5) + 
  geom_point(data = comp_cch, aes(x = variable, y = value, color = variant), position = position_dodge(width = 0.6)) +
  geom_point(data = comp_cch %>% group_by(variant, variable) %>% summarize(value = mean(value)), aes(x = variable, y = value, color = variant), position = position_dodge(width = 0.6), shape = 95, size = 5)
```

    ## `summarise()` has grouped output by 'variant'. You can override using the
    ## `.groups` argument.

``` r
cch4c_version2
```

![](Recombinastics_analysis_files/figure-gfm/GCaMP%20initial%20and%20peak%20comparisons-3.png)<!-- -->

``` r
iono4c_version2 <- ggplot() + theme(legend.position = "top", panel.grid.major.x = element_blank(), axis.text.x.bottom = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  labs(title = "Ionomycin", x = NULL, y = "Mean fluorescent intensity") +
  scale_y_continuous(limits = c(-200,6000)) +
  geom_hline(yintercept = 0, alpha = 0.5) + 
  geom_point(data = comp_iono, aes(x = variable, y = value, color = variant), position = position_dodge(width = 0.6)) +
  geom_point(data = comp_iono %>% group_by(variant, variable) %>% summarize(value = mean(value)), aes(x = variable, y = value, color = variant), position = position_dodge(width = 0.6), shape = 95, size = 5)
```

    ## `summarise()` has grouped output by 'variant'. You can override using the
    ## `.groups` argument.

``` r
#iono4c_version2

peak_figure <- cch4c_version2 | iono4c_version2
peak_figure
```

![](Recombinastics_analysis_files/figure-gfm/GCaMP%20initial%20and%20peak%20comparisons-4.png)<!-- -->

``` r
ggsave(file = "Plots/Peak_value.pdf", peak_figure, height = 2.5, width = 3)
```

Figure 4D - normalized curves to peak values

``` r
#carbachol
#rep 1
rep1_cch_mch_decay = rep1_mch_cch %>% mutate(norm_peak = signal/quantile(signal, 0.95))
rep1_cch_005_decay = rep1_005_cch %>% mutate(norm_peak = signal/quantile(signal, 0.95))
rep1_cch_013_decay = rep1_013_cch %>% mutate(norm_peak = signal/quantile(signal, 0.95))
rep1_cch_decay_normpeak = bind_rows(rep1_cch_mch_decay, rep1_cch_005_decay, rep1_cch_013_decay)

#rep 2
rep2_cch_mch_decay = rep2_mch_cch %>% mutate(norm_peak = signal/quantile(signal, 0.95))
rep2_cch_005_decay = rep2_005_cch %>% mutate(norm_peak = signal/quantile(signal, 0.95))
rep2_cch_013_decay = rep2_013_cch %>% mutate(norm_peak = signal/quantile(signal, 0.95))
rep2_cch_decay_normpeak = bind_rows(rep2_cch_mch_decay, rep2_cch_005_decay, rep2_cch_013_decay)

rep_normpeak_decay = bind_rows(rep1_cch_decay_normpeak, rep2_cch_decay_normpeak)
rep_normpeak_decay[rep_normpeak_decay == "mCherry_cch"] = "Negative control"
rep_normpeak_decay[rep_normpeak_decay == "pNK005_cch"] = "WT STIM1"
rep_normpeak_decay[rep_normpeak_decay == "pNK013_cch"] = "R429C STIM1"

cch_normpeak_4d = ggplot(data = rep_normpeak_decay, aes(x = ND.T, y = norm_peak, color = Name)) +
  theme_bw() + 
  scale_y_continuous(breaks = c(0,0.5,1)) + 
  scale_x_continuous(breaks = c(-30,0,30,60,90,120,150))+ 
  labs(x = "Time (s)", y = "Green MFI", color = "Variant")+ theme(text = element_text(size = 20)) + 
  geom_line(size = 1) +
  facet_wrap(~rep)
cch_normpeak_4d
```

![](Recombinastics_analysis_files/figure-gfm/Normalized%20GCaMP%20data%20for%20plot-1.png)<!-- -->

``` r
## Kenny mods

rep_normpeak_decay2 <- rep_normpeak_decay %>% group_by(Name, ND.T) %>% summarize(mean_normpeak_decay = mean(norm_peak))
```

    ## `summarise()` has grouped output by 'Name'. You can override using the
    ## `.groups` argument.

``` r
rep_normpeak_decay2$Name <- factor(rep_normpeak_decay2$Name, levels = c("Negative control", "WT STIM1", "R429C STIM1"))

cch_normpeak_average_decay <- ggplot() + theme(legend.position = "top") +
  labs(x = "Time (seconds)", y = "Normalized fluorescence", title = "Carbachol") +
  scale_x_continuous(breaks = c(-30,0,30,60,90,120,150))+ 
  geom_hline(yintercept = 0, alpha = 0.5) + 
  geom_point(data = rep_normpeak_decay2, aes(x = ND.T, y = mean_normpeak_decay, color = Name), alpha = 0.2)
cch_normpeak_average_decay
```

![](Recombinastics_analysis_files/figure-gfm/Normalized%20GCaMP%20data%20for%20plot-2.png)<!-- -->

``` r
#ionomycin
#rep 1
rep1_iono_mch_decay = rep1_mch_iono %>% mutate(norm_peak = signal/quantile(signal, 0.95))
rep1_iono_005_decay = rep1_005_iono %>% mutate(norm_peak = signal/quantile(signal, 0.95))
rep1_iono_013_decay = rep1_013_iono %>% mutate(norm_peak = signal/quantile(signal, 0.95))
rep1_iono_decay_normpeak = bind_rows(rep1_iono_mch_decay, rep1_iono_005_decay, rep1_iono_013_decay)

#rep 2
rep2_iono_mch_decay = rep2_mch_iono %>% mutate(norm_peak = signal/quantile(signal, 0.95))
rep2_iono_005_decay = rep2_005_iono %>% mutate(norm_peak = signal/quantile(signal, 0.95))
rep2_iono_013_decay = rep2_013_iono %>% mutate(norm_peak = signal/quantile(signal, 0.95))
rep2_iono_decay_normpeak = bind_rows(rep2_iono_mch_decay, rep2_iono_005_decay, rep2_iono_013_decay)

rep_normpeak_decay_iono = bind_rows(rep1_iono_decay_normpeak, rep2_iono_decay_normpeak)
rep_normpeak_decay_iono[rep_normpeak_decay_iono$Name == "mCherry_ionomycin","Name"] = "Negative control"
rep_normpeak_decay_iono[rep_normpeak_decay_iono$Name == "pNK005_ionomycin","Name"] = "WT STIM1"
rep_normpeak_decay_iono[rep_normpeak_decay_iono$Name == "pNK013_ionomycin","Name"] = "R429C STIM1"

iono_normpeak_4d = ggplot(data = rep_normpeak_decay_iono, aes(x = ND.T, y = norm_peak, color = Name))+theme_bw()+ scale_y_continuous(breaks = c(0,0.5,1))+ 
  scale_x_continuous(breaks = c(-30,0,30,60,90,120,150))+ 
  labs(x = "Time (s)", y = "Green MFI", color = "Variant")+ theme(text = element_text(size = 20)) + 
  geom_line(size = 1) +
  facet_wrap(~rep)
iono_normpeak_4d
```

![](Recombinastics_analysis_files/figure-gfm/Normalized%20GCaMP%20data%20for%20plot-3.png)<!-- -->

``` r
## Kenny mods
rep_normpeak_decay_iono2 <- rep_normpeak_decay_iono %>% group_by(Name, ND.T) %>% summarize(mean_normpeak_decay = mean(norm_peak))
```

    ## `summarise()` has grouped output by 'Name'. You can override using the
    ## `.groups` argument.

``` r
rep_normpeak_decay_iono2$Name <- factor(rep_normpeak_decay_iono2$Name, levels = c("Negative control", "WT STIM1", "R429C STIM1"))

iono_normpeak_average_decay <- ggplot() + theme(legend.position = "top") +
  labs(x = "Time (seconds)", y = "Normalized fluorescence", title = "Ionomycin") +
  scale_x_continuous(breaks = c(-30,0,30,60,90,120,150))+ 
  geom_hline(yintercept = 0, alpha = 0.5) + 
  geom_point(data = rep_normpeak_decay_iono2, aes(x = ND.T, y = mean_normpeak_decay, color = Name), alpha = 0.2)
iono_normpeak_average_decay
```

![](Recombinastics_analysis_files/figure-gfm/Normalized%20GCaMP%20data%20for%20plot-4.png)<!-- -->

``` r
# ggsave("cch_normpeak_4d.pdf", cch_normpeak_4d, width = 7, height = 4)
# ggsave("iono_normpeak_4d.pdf", iono_normpeak_4d, width = 7, height = 4)

grid.arrange(cch_normpeak_4d, iono_normpeak_4d)
```

![](Recombinastics_analysis_files/figure-gfm/Normalized%20GCaMP%20data%20for%20plot-5.png)<!-- -->

``` r
Decay_curves <- cch_normpeak_average_decay | iono_normpeak_average_decay
Decay_curves
```

![](Recombinastics_analysis_files/figure-gfm/Normalized%20GCaMP%20data%20for%20plot-6.png)<!-- -->

``` r
ggsave(file = "Plots/Decay_curves.pdf", Decay_curves, height = 2.5, width = 5)
```

## Double landing pad with protein interaction pairs

## Relevant to Figure 5 Vif \<-\> Apobec experiments

``` r
vif_apobec <- read.csv(file = "data/flow/Vif_Apobec/Orientations.csv")
vif_apobec_fusions <- vif_apobec %>% filter(gene != "None")
vif_apobec_egfp <- vif_apobec %>% filter(gene == "None")


Apobec_fusions_plot <- ggplot() + theme_bw() + theme(panel.grid.minor = element_blank()) +
  labs(x = "Geometric mean of MFI (Apobec only)", y = "Geometric mean of MFI\n(HIV-1 vif coexpressed)") +
  scale_x_log10() + scale_y_log10() +
  geom_abline(slope = 1, alpha = 0.2, linetype = 2, alpha = 0.5, size = 1.5) +
  geom_point(data = vif_apobec_fusions, aes(x = stop, y = wt, color = gene, shape = orientation), size = 3, alpha = 0.7) +
  geom_point(data = vif_apobec_egfp, aes(x = stop, y = wt), size = 3, shape = 15, alpha = 0.7)
```

    ## Warning: Duplicated aesthetics after name standardisation: alpha

``` r
Apobec_fusions_plot
```

![](Recombinastics_analysis_files/figure-gfm/Import%20the%20vif%20orientations%20data-1.png)<!-- -->

``` r
ggsave(file = "plots/Apobec_fusions_plot.pdf", Apobec_fusions_plot, height = 2.2, width = 4)
```

``` r
va_rep1 <- read.csv(file = "Data/Vif_apobec/Vif_apobec_rep1.csv", header = T)
for(x in 1:nrow(va_rep1)){
  for(y in 2:ncol(va_rep1)){
    va_rep1[x,y] <- va_rep1[x,y] / va_rep1[6,y]
  }
}
va_rep1_melt <- melt(va_rep1) %>% filter(variable != "EGFP.A3F")
```

    ## Using vif as id variables

``` r
colnames(va_rep1_melt) <- c("vif","apobec","normalized")
va_rep1_melt$rep <- 1

va_rep2 <- read.csv(file = "Data/Vif_apobec/Vif_apobec_rep2.csv", header = T)
for(x in 1:nrow(va_rep2)){
  for(y in 2:ncol(va_rep2)){
    va_rep2[x,y] <- va_rep2[x,y] / va_rep2[8,y]
  }
}
va_rep2_melt <- melt(va_rep2)
```

    ## Using vif as id variables

``` r
colnames(va_rep2_melt) <- c("vif","apobec","normalized")
va_rep2_melt$rep <- 2

va_rep3 <- read.csv(file = "Data/Vif_apobec/Vif_apobec_rep3.csv", header = T)
for(x in 1:nrow(va_rep3)){
  for(y in 2:ncol(va_rep3)){
    va_rep3[x,y] <- va_rep3[x,y] / va_rep3[8,y]
  }
}
va_rep3_melt <- melt(va_rep3)
```

    ## Using vif as id variables

``` r
colnames(va_rep3_melt) <- c("vif","apobec","normalized")
va_rep3_melt$rep <- 2

va_rep_combined <- rbind(va_rep1_melt, va_rep2_melt, va_rep3_melt)

va_rep_combined_summary <- va_rep_combined %>% group_by(vif, apobec) %>% summarize(mean = mean(normalized), sd = sd(normalized)) %>% mutate(mean_log2 = log2(mean), sd_log2 = log2(sd))
```

    ## `summarise()` has grouped output by 'vif'. You can override using the `.groups`
    ## argument.

``` r
va_rep_combined_summary$apobec <- factor(va_rep_combined_summary$apobec, levels = c("EGFP","A3G.EGFP","EGFP.A3D","EGFP.A3H2"))
va_rep_combined_summary$vif <- factor(va_rep_combined_summary$vif, levels = rev(c("WT","STOP","W89R","A149Y","W11A","H43E","F39A","G84K")))

vif_apobec_heatmap1 <- ggplot() + theme(axis.text.x.bottom = element_text(angle = 45, hjust = 1, vjust = 1)) +
  labs(x = "EGFP or Apobec fusion", y = "HIV-1 Vif variant") +
  #scale_fill_gradient2(low = "white", mid = "darkgreen", high = "red", na.value = "brown4", limits = c(0.25,1.3), midpoint = 1) + 
  scale_fill_gradient2(low = "white", mid = "darkgreen", high = "red", na.value = "green4", limits = c(log2(0.25),log2(1.3)), midpoint = 1) + 
  geom_tile(data = va_rep1_melt, aes(x = apobec, y = vif, fill = log2(normalized)))
vif_apobec_heatmap1
```

![](Recombinastics_analysis_files/figure-gfm/Vif%20Apobec%20heatmap-1.png)<!-- -->

``` r
ggsave(file = "Plots/Vif_apobec_heatmap1.pdf", vif_apobec_heatmap1, height = 2, width = 3.5)

vif_apobec_heatmap2 <- ggplot() + theme(axis.text.x.bottom = element_text(angle = 45, hjust = 1, vjust = 1)) +
  labs(x = "EGFP or Apobec fusion", y = "HIV-1 Vif variant") +
  #scale_fill_gradient2(low = "white", mid = "darkgreen", high = "red", na.value = "brown4", limits = c(0.25,1.3), midpoint = 1) + 
  scale_fill_gradient2(low = "white", mid = "darkgreen", high = "red", na.value = "green4", limits = c(log2(0.18),log2(1.32)), midpoint = 1) + 
  geom_tile(data = va_rep2_melt, aes(x = apobec, y = vif, fill = log2(normalized)))
vif_apobec_heatmap2
```

![](Recombinastics_analysis_files/figure-gfm/Vif%20Apobec%20heatmap-2.png)<!-- -->

``` r
ggsave(file = "Plots/Vif_apobec_heatmap2.pdf", vif_apobec_heatmap2, height = 2, width = 3.5)

vif_apobec_heatmap3 <- ggplot() + theme(axis.text.x.bottom = element_text(angle = 45, hjust = 1, vjust = 1)) +
  labs(x = "EGFP or Apobec fusion", y = "HIV-1 Vif variant") +
  #scale_fill_gradient2(low = "white", mid = "darkgreen", high = "red", na.value = "brown4", limits = c(0.25,1.3), midpoint = 1) + 
  scale_fill_gradient2(low = "white", mid = "darkgreen", high = "red", na.value = "green4", limits = c(log2(0.2),log2(1.62)), midpoint = 1) + 
  geom_tile(data = va_rep3_melt, aes(x = apobec, y = vif, fill = log2(normalized)))
vif_apobec_heatmap3
```

![](Recombinastics_analysis_files/figure-gfm/Vif%20Apobec%20heatmap-3.png)<!-- -->

``` r
ggsave(file = "Plots/Vif_apobec_heatmap3.pdf", vif_apobec_heatmap3, height = 2, width = 3.5)

vif_apobec_heatmap_combined <- ggplot() + theme(axis.text.x.bottom = element_text(angle = 45, hjust = 1, vjust = 1)) +
  labs(x = "EGFP or Apobec fusion", y = "HIV-1 Vif variant") +
  #scale_fill_gradient2(low = "white", mid = "darkgreen", high = "red", na.value = "brown4", limits = c(0.25,1.3), midpoint = 1) + 
  scale_fill_gradient2(low = "white", mid = "darkgreen", high = "red", na.value = "green4", limits = c(-2.15,0.4), midpoint = 1) + 
  geom_tile(data = va_rep_combined_summary, aes(x = apobec, y = vif, fill = mean_log2))
vif_apobec_heatmap_combined
```

![](Recombinastics_analysis_files/figure-gfm/Vif%20Apobec%20heatmap-4.png)<!-- -->

``` r
ggsave(file = "Plots/vif_apobec_heatmap_combined.pdf", vif_apobec_heatmap_combined, height = 2, width = 3)
```

``` r
va_rep_combined$vif <- factor(va_rep_combined$vif, levels = c("WT","STOP","W89R","A149Y","W11A","H43E","F39A","G84K"))
va_rep_combined_summary$vif <- factor(va_rep_combined_summary$vif, levels = c("WT","STOP","W89R","A149Y","W11A","H43E","F39A","G84K"))

EGFP_plot <- ggplot() + theme(panel.grid.major.x = element_blank(), plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  scale_y_log10(limits = c(0.18,1.62), breaks = c(0.25,0.5,1)) + labs(title = "EGFP", x = NULL) +
  geom_hline(yintercept = 1, alpha = 0.4) +
  geom_point(data = va_rep_combined %>% filter(apobec == "EGFP"), aes(x = vif, y = normalized), alpha = 0.4, color = "red") +
  geom_point(data = va_rep_combined_summary %>% filter(apobec == "EGFP"), aes(x = vif, y = mean), shape = 95, size = 6)
EGFP_plot
```

![](Recombinastics_analysis_files/figure-gfm/Individual%20Vif-Apobec%20data%20for%20supplementary%20figure-1.png)<!-- -->

``` r
A3G.EGFP_plot <- ggplot() + theme(panel.grid.major.x = element_blank(), plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  scale_y_log10(limits = c(0.18,1.62), breaks = c(0.25,0.5,1)) + labs(title = "A3G-EGFP", x = NULL) +
  geom_hline(yintercept = 1, alpha = 0.4) +
  geom_point(data = va_rep_combined %>% filter(apobec == "A3G.EGFP"), aes(x = vif, y = normalized), alpha = 0.4, color = "red") +
  geom_point(data = va_rep_combined_summary %>% filter(apobec == "A3G.EGFP"), aes(x = vif, y = mean), shape = 95, size = 6)
A3G.EGFP_plot
```

![](Recombinastics_analysis_files/figure-gfm/Individual%20Vif-Apobec%20data%20for%20supplementary%20figure-2.png)<!-- -->

``` r
EGFP.A3D_plot <- ggplot() + theme(panel.grid.major.x = element_blank(), plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  scale_y_log10(limits = c(0.18,1.62), breaks = c(0.25,0.5,1)) + labs(title = "EGFP-A3D", x = NULL) +
  geom_hline(yintercept = 1, alpha = 0.4) +
  geom_point(data = va_rep_combined %>% filter(apobec == "EGFP.A3D"), aes(x = vif, y = normalized), alpha = 0.4, color = "red") +
  geom_point(data = va_rep_combined_summary %>% filter(apobec == "EGFP.A3D"), aes(x = vif, y = mean), shape = 95, size = 6)
EGFP.A3D_plot
```

![](Recombinastics_analysis_files/figure-gfm/Individual%20Vif-Apobec%20data%20for%20supplementary%20figure-3.png)<!-- -->

``` r
EGFP.A3H2_plot <- ggplot() + theme(panel.grid.major.x = element_blank(), plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  scale_y_log10(limits = c(0.18,1.62), breaks = c(0.25,0.5,1)) + labs(title = "EGFP-A3H2", x = NULL) +
  geom_hline(yintercept = 1, alpha = 0.4) +
  geom_point(data = va_rep_combined %>% filter(apobec == "EGFP.A3H2"), aes(x = vif, y = normalized), alpha = 0.4, color = "red") +
  geom_point(data = va_rep_combined_summary %>% filter(apobec == "EGFP.A3H2"), aes(x = vif, y = mean), shape = 95, size = 6)
EGFP.A3H2_plot
```

![](Recombinastics_analysis_files/figure-gfm/Individual%20Vif-Apobec%20data%20for%20supplementary%20figure-4.png)<!-- -->

``` r
Apobec_vif_supp_fig <- EGFP_plot | A3G.EGFP_plot | EGFP.A3D_plot | EGFP.A3H2_plot
Apobec_vif_supp_fig
```

![](Recombinastics_analysis_files/figure-gfm/Individual%20Vif-Apobec%20data%20for%20supplementary%20figure-5.png)<!-- -->

``` r
ggsave(file = "Plots/Apobec_vif_supp_fig.pdf", Apobec_vif_supp_fig, height = 2, width = 8)
```

## Testing additional GN dinucleotide pairs in combination

## Relevant to Figure 6

``` r
attB_plasmid_mixture <- read.delim(file = "data/GN_recombinations/attBmix_outfile.tsv", header = T, stringsAsFactors = F)
attB_plasmid_mixture_table <- data.frame(table(attB_plasmid_mixture$attB))
attB_plasmid_mixture_table$Freq2 <- attB_plasmid_mixture_table$Freq / sum(attB_plasmid_mixture_table$Freq)
colnames(attB_plasmid_mixture_table) <- c("attB","start_count","start_freq")

## Rep 1
rep1_ga <- read.delim(file = "data/GN_recombinations/Rep1_GA_outfile.tsv", header = T, stringsAsFactors = F)
rep1_ga2 <- data.frame(table(rep1_ga$attB)) %>% mutate(attP = "GA")
rep1_gt <- read.delim(file = "data/GN_recombinations/Rep1_GT_outfile.tsv", header = T, stringsAsFactors = F)
rep1_gt2 <- data.frame(table(rep1_gt$attB)) %>% mutate(attP = "GT")
rep1_gc <- read.delim(file = "data/GN_recombinations/Rep1_GC_outfile.tsv", header = T, stringsAsFactors = F)
rep1_gc2 <- data.frame(table(rep1_gc$attB)) %>% mutate(attP = "GC")
rep1_gg <- read.delim(file = "data/GN_recombinations/Rep1_GG_outfile.tsv", header = T, stringsAsFactors = F)
rep1_gg2 <- data.frame(table(rep1_gg$attB)) %>% mutate(attP = "GG")

rep1_comb <- rbind(rep1_ga2, rep1_gt2, rep1_gc2, rep1_gg2)
colnames(rep1_comb) <- c("attB","Freq","attP")

gn_rep1_results_frame <- data.frame("attB" = NA, "attP" = NA, "final_freq" = NA, "start_freq" = NA)
for( x in c("GT","GA","GC","GG")){
  temp_subset <- subset(rep1_comb, attP == x)
  temp_subset$final_freq <- temp_subset$Freq / sum(temp_subset$Freq)
  temp_subset2 <- merge(temp_subset[,c("attP","attB","final_freq")], attB_plasmid_mixture_table[,c("attB","start_freq")])
  gn_rep1_results_frame <- rbind(gn_rep1_results_frame, temp_subset2)
}
gn_rep1_results_frame <- gn_rep1_results_frame %>% filter(!is.na(attB)) %>% mutate(rep = 1)

## Rep 2
gn_rep2 <- read.delim(file = "data/GN_recombinations/Rep2_outfile.tsv", header = T, stringsAsFactors = F)
gn_rep2_table <- data.frame(with(gn_rep2,table(attP,attB)))

gn_rep2_results_frame <- data.frame("attB" = NA, "attP" = NA, "final_freq" = NA, "start_freq" = NA)
for( x in c("GT","GA","GC","GG")){
  temp_subset <- subset(gn_rep2_table, attP == x)
  temp_subset$final_freq <- temp_subset$Freq / sum(temp_subset$Freq)
  temp_subset2 <- merge(temp_subset[,c("attP","attB","final_freq")], attB_plasmid_mixture_table[,c("attB","start_freq")])
  gn_rep2_results_frame <- rbind(gn_rep2_results_frame, temp_subset2)
}
gn_rep2_results_frame <- gn_rep2_results_frame %>% filter(!is.na(attB)) %>% mutate(rep = 2)

## Rep 3
gn_rep3 <- read.delim(file = "data/GN_recombinations/Rep3_outfile.tsv", header = T, stringsAsFactors = F)
gn_rep3_table <- data.frame(with(gn_rep3,table(attP,attB)))

gn_rep3_results_frame <- data.frame("attB" = NA, "attP" = NA, "final_freq" = NA, "start_freq" = NA)
for( x in c("GT","GA","GC","GG")){
  temp_subset <- subset(gn_rep3_table, attP == x)
  temp_subset$final_freq <- temp_subset$Freq / sum(temp_subset$Freq)
  temp_subset2 <- merge(temp_subset[,c("attP","attB","final_freq")], attB_plasmid_mixture_table[,c("attB","start_freq")])
  gn_rep3_results_frame <- rbind(gn_rep3_results_frame, temp_subset2)
}
gn_rep3_results_frame <- gn_rep3_results_frame %>% filter(!is.na(attB)) %>% mutate(rep = 3)

## Combine the results of all three replicates

gn_combined <- rbind(gn_rep1_results_frame, gn_rep2_results_frame, gn_rep3_results_frame)

gn_combined$enrichment <- gn_combined$final_freq/gn_combined$start_freq

write.csv(file = "data/GN_recombination_results.csv", gn_combined, row.names = FALSE)

gn_combined_summary <- gn_combined %>% group_by(attP, attB) %>% summarize(mean = mean(final_freq), sd = sd(final_freq))
```

    ## `summarise()` has grouped output by 'attP'. You can override using the
    ## `.groups` argument.

``` r
gn_combined_summary2 <- rbind(gn_combined_summary, data.frame("attP" = "Input", "attB" = attB_plasmid_mixture_table$attB, "mean" = attB_plasmid_mixture_table$start_freq, "sd" = NA))
```

``` r
gn_combined$attP <- factor(gn_combined$attP, levels = c("GT","GA","GC","GG"))
gn_combined_summary2$attP <- factor(gn_combined_summary2$attP, levels = c("Input", "GT","GA","GC","GG"))

Orthogonal_plasmid_recombinations <- ggplot() + scale_y_log10() + theme(panel.grid.major.x = element_blank()) +
  labs(x = "AttP plasmid (individually tested)", y = "Fraction of\nrecombinants", color = "AttB\nplasmid\n(mixed)") + 
  geom_point(data = gn_combined_summary2, aes(x = attP, y = mean, color = attB), position = position_dodge(width = 0.6), shape = 95, size = 4) +
  geom_point(data = gn_combined, aes(x = attP, y = final_freq, color = attB), position = position_dodge(width = 0.6), alpha = 0.4) +
  #geom_errorbar(data = gn_combined_summary2, aes(x = attP, ymin = mean - sd/sqrt(2)*1.96, ymax = mean + sd/sqrt(3)*1.96, color = attB), position = position_dodge(width = 0.6), width = 0.5, alpha = 0.5) +
  geom_point(data = gn_combined_summary2, aes(x = attP, y = mean, color = attB), position = position_dodge(width = 0.6), shape = 95, size = 4) +
  NULL
ggsave(file = "plots/Orthogonal_plasmid_recombinations.pdf", Orthogonal_plasmid_recombinations, height = 2, width = 4.2)
Orthogonal_plasmid_recombinations
```

![](Recombinastics_analysis_files/figure-gfm/Recombination%20mixture%20experiment%20plot-1.png)<!-- -->

``` r
recomb_mixture_test <- gn_combined_summary[1:16,c("attP","attB")]
colnames(recomb_mixture_test) <- c("lp1","lp2")

simulation_test_size = 50000
recomb_mixture_test$on_target <- NA

for(x in 1:nrow(recomb_mixture_test)){
  recomb_mix_lp1 <- recomb_mixture_test$lp1[x]
  recomb_mix_lp2 <- recomb_mixture_test$lp2[x]
  if(recomb_mix_lp1 == recomb_mix_lp2){
    lp1_sampling <- sample(x = c(1,0), prob = c(0.5,0.5) , size = simulation_test_size, replace = T)
    lp2_sampling <- sample(x = c(1,0), prob = c(0.5,0.5) , size = simulation_test_size, replace = T)
    results_frame <- data.frame("lp1_on_target" = lp1_sampling == 1, "lp2_on_target" = lp2_sampling == 1)
    results_frame$sum <- results_frame$lp1_on_target + results_frame$lp2_on_target
    recomb_mixture_test$on_target[x] <- sum(results_frame$sum == 1)
  }
  if(recomb_mix_lp1 != recomb_mix_lp2){
    lp1_probs <- gn_combined_summary %>% filter(attP == recomb_mix_lp1 & attB %in% c(recomb_mix_lp1, recomb_mix_lp2))
    lp2_probs <- gn_combined_summary %>% filter(attP == recomb_mix_lp2 & attB %in% c(recomb_mix_lp1, recomb_mix_lp2))
    lp1_sampling <- sample(x = lp1_probs$attB, prob = lp1_probs$mean , size = simulation_test_size, replace = T)
    lp2_sampling <- sample(x = lp2_probs$attB, prob = lp2_probs$mean , size = simulation_test_size, replace = T)
    results_frame <- data.frame("lp1_on_target" = lp1_sampling == recomb_mix_lp1, "lp2_on_target" = lp2_sampling == recomb_mix_lp2)
    results_frame$sum <- results_frame$lp1_on_target + results_frame$lp2_on_target
    recomb_mixture_test$on_target[x] <- sum(results_frame$sum == 2)
  }
}
recomb_mixture_test$frac_on_target <- recomb_mixture_test$on_target / simulation_test_size

Double_LP_heatmap <- ggplot() + labs(x = "First landing pad", y = "Second landing pad", fill = "Fraction\nwith unique\nrecombinant\nplasmid pair") + 
  scale_fill_gradient2(low = "white", high = "red", midpoint = 0.8) + 
  geom_tile(data = recomb_mixture_test, aes(x = lp1, y = lp2, fill = frac_on_target)) +
  geom_text(data = recomb_mixture_test, aes(x = lp1, y = lp2, label = round(frac_on_target,2)), size = 3)
ggsave(file = "plots/Double_LP_heatmap.pdf", Double_LP_heatmap, height = 2, width = 3.3)
Double_LP_heatmap
```

![](Recombinastics_analysis_files/figure-gfm/Modeling%20predicted%20pairwise%20specificities%20based%20on%20the%20GN%20data-1.png)<!-- -->

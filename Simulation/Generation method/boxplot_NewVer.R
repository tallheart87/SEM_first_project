# read all replication results
path <- "REP/REP3"
files <- list.files(path, pattern = "\\.rds$", full.names = TRUE)
criteria_list <- lapply(files, readRDS)
names(criteria_list) <- paste0("rep", c(1:length(files)))

# simulation time
path2 <- "REP/REP3/time"
files2 <- list.files(path2, pattern = "\\.rds$", full.names = TRUE)
criteria_list2 <- lapply(files2, readRDS)

library(dplyr)
library(purrr)
library(ggplot2)
# 1. Tidy the nested list into a long data frame
## criteria_list: list of replications
#methods <- c("SUM", "SEMReg", "SEMBar", 
 #                 "SAM", "SGCCA", "rESEM", 
  #                "SEMRule", "PLS", "GLM", "ELA")
methods <- c("SumScore", 
             "SEM_Reg", "SEM_Bar",
             "SAM","SAM_Reg",
             "SGCCA","rESEM", 
             "SEM_BASED", "PLS", 
             "GLM", "elastic")
df_mae <- imap_dfr(criteria_list, function(rep_obj, rep_id) {
  #rep_obj: list("M:weak" = ..., "M:strong" = ...)
  map_dfr(names(rep_obj), function(m_name) {
    s_level <- rep_obj[[m_name]]          #list("S:weak" = ..., "S:strong" = ...)
    map_dfr(names(s_level), function(s_name) {
      size_level <- s_level[[s_name]]     #list("100" = ..., "200" = ..., ...)
      map_dfr(names(size_level), function(n_name) {
        metrics <- size_level[[n_name]]   #list with MAE, RMSE, etc.
        mae_vec <- as.numeric(metrics$MAE) #MAE
        tibble(
          replication = rep_id,
          M = m_name,
          S = s_name,
          condition = paste(m_name, s_name, sep = " / "),
          n = as.integer(n_name),
          method = methods,
          MAE = mae_vec
        )
      })
    })
  })
})
df_mae <- df_mae %>% 
  mutate(method = factor(method, levels = methods))

# convergence rate
df_mae %>%
  group_by(method) %>%
  summarise(NA_rate = 1 - mean(is.na(MAE)))

#2. Boxplots of MAE by method, faceted by condition and sample size
p <- ggplot(df_mae, aes(x = method, y = MAE, fill = method)) +
  geom_boxplot(outlier.shape = NA) +                     # box per method (33 values)
  geom_jitter(width = 0.15, alpha = 0.5, size = 1) +     # show the 33 rep points
  facet_grid(condition ~ n, scales = "free_y") +         # one panel per condition × n
  labs(title = paste0("MAE across " , length(files), " replications"),
       x = "Method", y = "MAE") +
  scale_y_continuous(limits = c(0, 2)) + #scale the y-axis
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))
print(p)

# no SAM
## delete form the matrix
df_mae_noSAM <- df_mae %>% 
  filter(method != c("SAM")) %>% 
  filter(method != c("SAM_Reg"))
## plot the boxplot
p_noSAM <- ggplot(df_mae_noSAM, aes(x = method, y = MAE, fill = method)) +
  geom_boxplot(outlier.shape = NA) +                     # box per method (33 values)
  geom_jitter(width = 0.15, alpha = 0.5, size = 1) +     # show the 33 rep points
  facet_grid(condition ~ n, scales = "free_y") +         # one panel per condition × n
  labs(title = paste0("MAE across" , length(files), "replications"),
       x = "Method", y = "MAE") +
  #scale_y_continuous(limits = c(0, 1)) + #scale the y-axis
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))
print(p_noSAM)

# sperate into M:weak and M:strong
## Plot for M:weak
p_weak <- df_mae_noSAM %>%
  filter(M == "M:weak") %>%
  ggplot(aes(x = method, y = MAE, fill = method)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1) +
  facet_grid(S ~ n, scales = "free_y") +
  labs(title = paste0("MAE (M:weak) across ", length(files), " replications"),
       x = "Method", y = "MAE") +
  scale_y_continuous(limits = c(0, 1)) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

## Plot for M:strong
p_strong <- df_mae_noSAM %>%
  filter(M == "M:strong") %>%
  ggplot(aes(x = method, y = MAE, fill = method)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1) +
  facet_grid(S ~ n, scales = "free_y") +
  labs(title = paste0("MAE (M:strong) across ", length(files), " replications"),
       x = "Method", y = "MAE") +
  scale_y_continuous(limits = c(0, 1)) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

## Print both
print(p_weak)
print(p_strong)

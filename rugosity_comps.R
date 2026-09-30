library(tidyverse)

folder_path <- "/Users/ellag/Desktop/PhD/academic_projects/3d_reconstruction/data/rugosities/"

files <- list.files(folder_path, pattern = "\\.csv$", full.names = TRUE)

combined_df <- map_dfr(files, function(file) {
  
  df <- read_csv(file)
  
  filename <- basename(file)
  name <- str_remove(filename, "\\.csv$")
  
  parts <- str_split(name, "_")[[1]]
  
  # Extract radius
  radius <- tail(parts, 1)
  
  # Find "rugosity"
  rugosity_index <- which(parts == "rugosity")
  
  # Extract trial ID
  trial_id <- str_c(parts[1:(rugosity_index - 1)], collapse = "_")
  
  df %>%
    mutate(
      trial_id = trial_id,
      radius = radius
    )
})

filtered <- combined_df %>%
  filter(trial_id == "LG_EH_062523_7_L_3fps_opt-ahojnnes_settings" | trial_id == "LG_EH_062523_7_L_3fps_opt")

combined_df1 <- filtered %>%
  mutate(
    trial=str_extract(trial_id, "^[^_]+_[^_]+"),
    fps = str_extract(trial_id, "\\d+fps"),
    optimizations = str_extract(trial_id, "(?<=fps_).*")
  ) %>%
  mutate(optimizations=case_when(optimizations=="noopt1" ~ "noopt",
                                 optimizations=="opt2" ~ "opt", 
                                 TRUE ~ optimizations)) %>%
  mutate(fps_opt = str_c(fps, optimizations, sep = "_"))


# Include radius as a thing
combined_wide<- combined_df1 %>%
  dplyr::select(trial, frame, fps_opt, radius, local_rugosity) %>%pivot_wider(
    names_from = fps_opt,
    values_from = local_rugosity
  )

### Graph Making Below. Many options! 
### Test for linear fit: x

combined_wide %>%
  ggplot(aes(color=trial, y = `1fps_noopt`-`3fps_noopt` ,x=`3fps_noopt`)) + 
  facet_grid(trial~radius)+
  geom_point(size=.5, alpha=.5)+
  geom_smooth(method = "lm", color = "black", linewidth = .5)+
  #geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed", linewidth = 1)+ # Add 1:1 line
  theme_bw()

combined_wide %>%
  ggplot(aes(color=trial, y = `1fps_opt`-`3fps_noopt` ,x=`3fps_noopt`)) + 
  facet_grid(trial~radius)+
  geom_point(size=.5, alpha=.5)+
  geom_smooth(method = "lm", color = "black", linewidth = .5)+
  #geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed", linewidth = 1)+ # Add 1:1 line
  theme_bw()

combined_wide %>%
  ggplot(aes(color=trial, y = `3fps_opt`-`3fps_noopt` ,x=`3fps_noopt`)) + 
  facet_grid(trial~radius)+
  geom_point(size=.5, alpha=.5)+
  geom_smooth(method = "lm", color = "black", linewidth = .5)+
  #geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed", linewidth = 1)+ # Add 1:1 line
  theme_bw()

combined_wide %>%
  ggplot(aes(color=trial, y = `1fps_opt`-`3fps_opt` ,x=`3fps_opt`)) + 
  facet_grid(trial~radius)+
  geom_point(size=.5, alpha=.5)+
  geom_smooth(method = "lm", color = "black", linewidth = .5)+
  #geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed", linewidth = 1)+ # Add 1:1 line
  theme_bw()

combined_wide %>%
  ggplot(aes(color=trial, y = `1fps_opt`-`1fps_noopt` ,x=`1fps_noopt`)) + 
  facet_grid(trial~radius)+
  geom_point(size=.5, alpha=.5)+
  geom_smooth(method = "lm", color = "black", linewidth = .5)+
  #geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed", linewidth = 1)+ # Add 1:1 line
  theme_bw()

combined_wide %>%
  ggplot(aes(color=trial, y = `1fps_opt`-`3fps_opt` ,x=`3fps_opt`)) + 
  facet_grid(trial~radius)+
  geom_point(size=.5, alpha=.5)+
  geom_smooth(method = "lm", color = "black", linewidth = .5)+
  #geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed", linewidth = 1)+ # Add 1:1 line
  theme_bw()

### Next, plot actual data to compare settings ###
# Plot 1 (upper left)
combined_wide %>%
  filter(radius=="0.5m")%>%
  ggplot(aes(color=trial, y = `3fps_opt-ahojnnes_settings` ,x=`3fps_opt`)) + 
  geom_point(size=.5, alpha=.5)+
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed", linewidth = 1)+ # Add 1:1 line
  theme_classic()

combined_wide %>%
  filter(radius=="1m")%>%
  ggplot(aes(color=trial, y = `3fps_opt` ,x=`3fps_noopt`)) + 
  geom_point(size=.5, alpha=.5)+
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed", linewidth = 1)+ # Add 1:1 line
  theme_classic()


# Plot 2 (bottom left)
combined_wide %>%
  filter(radius=="3m")%>%
  ggplot(aes(color=trial, y = `1fps_noopt` ,x=`3fps_noopt`)) + 
  geom_point(size=.5, alpha=.5)+
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed", linewidth = 1)+ # Add 1:1 line
  theme_classic()

# Plot 3 (top right)
combined_wide %>%
  filter(radius=="3m")%>%
  ggplot(aes(color=trial, y = `1fps_opt` ,x=`1fps_noopt`)) + 
  geom_point(size=.5, alpha=.5)+
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed", linewidth = 1)+ # Add 1:1 line
  theme_classic()

# Plot 4 (bottom right)
combined_wide %>%
  filter(radius=="3m")%>%
  ggplot(aes(color=trial, y = `1fps_opt` ,x=`3fps_opt`)) + 
  geom_point(size=.5, alpha=.5)+
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed", linewidth = 1)+ # Add 1:1 line
  theme_classic()

# Plot 5 (next page)
combined_wide %>%
  filter(radius=="3m")%>%
  ggplot(aes(color=trial, y = `1fps_opt` ,x=`3fps_noopt`)) + 
  geom_point(size=.5, alpha=.5)+
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed", linewidth = 1)+ # Add 1:1 line
  theme_classic()





### Testing correlations ###
combined_wide_3m <- combined_wide %>% filter(radius == "3m")
combined_wide_1m <- combined_wide %>% filter(radius == "1m")
combined_wide_half <- combined_wide %>% filter(radius == "0.5m")

cor.test(combined_wide_3m$"3fps_noopt", combined_wide_3m$"3fps_opt", na.action = na.exclude)
cor.test(combined_wide_3m$"3fps_noopt", combined_wide_3m$"1fps_noopt", na.action = na.exclude)
cor.test(combined_wide_3m$"1fps_noopt", combined_wide_3m$"1fps_opt", na.action = na.exclude)
cor.test(combined_wide_3m$"3fps_opt", combined_wide_3m$"1fps_opt", na.action = na.exclude)
cor.test(combined_wide_3m$`3fps_noopt`, combined_wide_3m$`1fps_opt`,na.action = na.exclude)

cor.test(combined_wide_1m$"3fps_noopt", combined_wide_1m$"3fps_opt", na.action = na.exclude)
cor.test(combined_wide_1m$"3fps_noopt", combined_wide_1m$"1fps_noopt", na.action = na.exclude)
cor.test(combined_wide_1m$"1fps_noopt", combined_wide_1m$"1fps_opt", na.action = na.exclude)
cor.test(combined_wide_half$"3fps_opt", combined_wide_half$"3fps_opt-ahojnnes_settings", na.action = na.exclude)
cor.test(combined_wide_1m$`3fps_noopt`, combined_wide_1m$`1fps_opt`,na.action = na.exclude)

cor.test(combined_wide_half$"3fps_noopt", combined_wide_half$"3fps_opt", na.action = na.exclude)
cor.test(combined_wide_half$"3fps_noopt", combined_wide_half$"1fps_noopt", na.action = na.exclude)
cor.test(combined_wide_half$"1fps_noopt", combined_wide_half$"1fps_opt", na.action = na.exclude)
cor.test(combined_wide_half$"3fps_opt", combined_wide_half$"1fps_opt", na.action = na.exclude)
cor.test(combined_wide_half$`3fps_noopt`, combined_wide_half$`1fps_opt`,na.action = na.exclude)

#### Then do a regression to see if it differs from a slope of one
library(lme4)
library(car)
library(emmeans)

## 3m radius models
m1_3m <- lm(`3fps_opt-ahojnnes_settings` ~ `3fps_opt`, data = combined_wide_half)
#plot(m1_3m)
Anova(m1_3m)
summary(m1_3m) # Use the estimate from this for slideshow

m2_3m <- lm(`3fps_noopt` ~ `1fps_noopt`, data = combined_wide_3m)
#plot(m1_3m)
Anova(m2_3m)
summary(m2_3m) # Use the estimate from this for slideshow

m3_3m <- lm(`1fps_noopt` ~ `1fps_opt`, data = combined_wide_3m)
#plot(m1_3m)
Anova(m3_3m)
summary(m3_3m) # Use the estimate from this for slideshow

m4_3m <- lm(`3fps_opt` ~ `1fps_opt`, data = combined_wide_3m)
#plot(m1_3m)
Anova(m4_3m)
summary(m4_3m) # Use the estimate from this for slideshow

## 1m radius models
m1_1m <- lm(`3fps_noopt` ~ `3fps_opt`, data = combined_wide_1m)
#plot(m1_3m)
Anova(m1_1m)
summary(m1_1m) # Use the estimate from this for slideshow

m2_1m <- lm(`3fps_noopt` ~ `1fps_noopt`, data = combined_wide_1m)
#plot(m1_3m)
Anova(m2_1m)
summary(m2_1m) # Use the estimate from this for slideshow

m3_1m <- lm(`1fps_noopt` ~ `1fps_opt`, data = combined_wide_1m)
#plot(m1_3m)
Anova(m3_1m)
summary(m3_1m) # Use the estimate from this for slideshow

m4_1m <- lm(`3fps_opt` ~ `1fps_opt`, data = combined_wide_1m)
#plot(m1_3m)
Anova(m4_1m)
summary(m4_1m) # Use the estimate from this for slideshow


## 0.5m radius models
m1_half <- lm(`3fps_noopt` ~ `3fps_opt`, data = combined_wide_half)
#plot(m1_3m)
Anova(m1_half)
summary(m1_half) # Use the estimate from this for slideshow

m2_half <- lm(`3fps_noopt` ~ `1fps_noopt`, data = combined_wide_half)
#plot(m1_3m)
Anova(m2_half)
summary(m2_half) # Use the estimate from this for slideshow

m3_half <- lm(`1fps_noopt` ~ `1fps_opt`, data = combined_wide_half)
#plot(m1_3m)
Anova(m3_half)
summary(m3_half) # Use the estimate from this for slideshow

m4_half <- lm(`3fps_opt` ~ `1fps_opt`, data = combined_wide_half)
#plot(m1_3m)
Anova(m4_half)
summary(m4_half) # Use the estimate from this for slideshow



### Optional: Also test trial-specific correlations 
# Then do a regression to see if it differs from a slope of one

m1.test <- lm(`3fps_noopt` ~ trial*`3fps_opt`*radius, data = combined_all)
plot(m1.test)

summary(m1.test)

Anova(m1.test)


# Get slopes
emtrends(m1.test, ~ trial*radius, var = '3fps_opt')



m2.test <- lm(`3fps_noopt` ~ trial*`1fps_opt`*radius, data = combined_all)
plot(m2.test)

summary(m2.test)
Anova(m2.test)

# Get slopes
emtrends(m2.test, ~ trial*radius, var = '1fps_opt')

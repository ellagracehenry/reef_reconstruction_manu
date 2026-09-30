library(ggplot2);library(reshape2)
library(dplyr)
library(tidyr)

JM151 <- read.csv("/Users/ellag/Desktop/PhD/academic_projects/3d_reconstruction/data/raycast_data/coral_sand_results_JM151_20260701_184840.csv")
ER23 <- read.csv("/Users/ellag/Desktop/PhD/academic_projects/3d_reconstruction/data/raycast_data/coral_sand_results_ER23_20260701_142828.csv")

hist(JM151$HRI)
x <- data.frame(JM151$HRI, JM151$HRI_01, JM151$HRI_1, JM151$HRI_sqrt05)
data <- melt(x)
ggplot(data,aes(x=value, fill=variable)) + geom_density(alpha=0.25)
ggplot(data,aes(x=value, fill=variable)) + geom_density(alpha=0.25)
ggplot(data,aes(x=variable, y=value, fill=variable)) + geom_boxplot(0.25)




plot(JM151$avg_distance, JM151$H, 
     col= penguins$species,
     pch=16,
     xlab="Bill Length (mm)", #label for x-axis
     ylab="Bill Depth (mm)", #label for y-axis
     main="A: Bill Length by Bill Depth ") #plot title


df_long <- JM151 %>%
pivot_longer(
  cols = starts_with("HRI"),
  names_to = "Method",
  values_to = "HRI"
)


#HISTOGRAMS
par(mfrow = c(2, 2))
min(JM151$HRI)
max(JM151$HRI)
hist(JM151$HRI)

min(JM151$HRI_01)
max(JM151$HRI_01)
hist(JM151$HRI_01)

min(JM151$HRI_1)
max(JM151$HRI_1)
hist(JM151$HRI_1)

min(JM151$HRI_sqrt05)
max(JM151$HRI_sqrt05)
hist(JM151$HRI_sqrt05)


#ALL HRI TOGETHER WITH DISTANCE
ggplot(df_long, aes(x = avg_distance, y = HRI, colour = Method)) +
  #geom_line() +
  geom_point(alpha = 0.4) +
  theme_minimal() +
  ylim(0,2.5)


par(mfrow = c(2, 2))
#HRI INDIVIDUALS WITH DISTANCE
plot(JM151$avg_distance, JM151$HRI)

plot(JM151$avg_distance, JM151$HRI_01)

plot(JM151$avg_distance, JM151$HRI_1)

plot(JM151$avg_distance, JM151$HRI_sqrt05)

par(mfrow = c(2, 2))
#Ray histograms
pos_44 <- JM151_rays %>%
  filter(frame == 44)
hist(pos_44$norm_height, xlim = c(0, 0.8), breaks = 10)
hist(pos_44$distance, xlim = c(0.3, 2), breaks = 10)

pos_1044 <- JM151_rays %>%
  filter(frame == 1044)
hist(pos_1044$norm_height, xlim = c(0, 0.8), breaks = 10)
hist(pos_1044$distance, xlim = c(0.3, 2), breaks = 10)

nrow(pos_44)
nrow(pos_1044)

comb_pos <- rbind(pos_1044, pos_44)

comb_pos$frame <- as.factor(comb_pos$frame)

comb_pos$frame <- factor(comb_pos$frame, levels = c("1044", "44"))

distance_plot <- ggplot(comb_pos) +
  geom_histogram(aes(distance, fill = frame), alpha = 0.7, bins = 10) +
  facet_wrap(~frame) +
  scale_x_continuous(limits = c(0, 2)) +
  theme_classic() +
  theme(plot.title = element_blank(),
        strip.text.x = element_blank()) +
  guides(fill="none") +
  labs(x = "Ray hit distance", y = "Count")
  
  

height_plot <- ggplot(comb_pos) +
  geom_histogram(aes(norm_height, fill = frame), alpha = 0.7, bins = 10) +
  facet_wrap(~frame) +
  scale_x_continuous(limits = c(0, 0.8)) +
  theme_classic() +
  theme(plot.title = element_blank(),
        strip.text.x = element_blank()) +
  guides(fill="none") +
  labs(x = "Ray hit height")


grid.arrange(distance_plot, height_plot, ncol=1)


ggplot(JM151, aes(x = bin, y = proportion)) +
  geom_col() +
  labs(x = "Bin", y = "Proportion") +
  theme_minimal()


# Convert distance bins to long format
distance_long <- JM151 %>%
  pivot_longer(
    cols = starts_with("distance_bin_"),
    names_to = "bin",
    values_to = "Proportion"
  ) %>%
  mutate(
    type = "Distance",
    midpoint = sapply(bin, function(x) {
      nums <- gsub("distance_bin_|m", "", x)
      bounds <- as.numeric(strsplit(nums, "_")[[1]])
      mean(bounds)
    })
  )

# Convert height bins to long format
height_long <- JM151 %>%
  pivot_longer(
    cols = starts_with("height_bin_"),
    names_to = "bin",
    values_to = "Proportion"
  ) %>%
  mutate(
    type = "Height",
    midpoint = sapply(bin, function(x) {
      nums <- gsub("height_bin_|m", "", x)
      bounds <- as.numeric(strsplit(nums, "_")[[1]])
      mean(bounds)
    })
  )

# Combine
JM151_long <- bind_rows(distance_long, height_long)

JM151_long_distance <- JM151_long %>%
  filter(frame %in% c(1044,44)) %>%
  filter(type == "Distance")

JM151_long_height <- JM151_long %>%
  filter(frame %in% c(1044,44)) %>%
  filter(type == "Height")

JM151_long_distance$frame <- factor(JM151_long_distance$frame, levels = c("1044", "44"))
JM151_long_height$frame <- factor(JM151_long_height$frame, levels = c("1044", "44"))


plot1 <- ggplot(JM151_long_distance, aes(x=midpoint, y=Proportion, fill = frame)) +
  geom_col(alpha = 0.7,
           width = 0.5) +
  facet_wrap(~frame) +
  scale_x_continuous(limits = c(0, 10)) +
  scale_fill_manual(values = c("#005EA0","#C0FF00")) +
  theme_classic() +
  theme(plot.title = element_blank(),
        strip.text.x = element_blank()) +
  guides(fill="none") +
  labs(x = "Ray hit distance", y = "Proportion")

plot2 <- ggplot(JM151_long_height, aes(x=midpoint, y=Proportion, fill = frame)) +
  geom_col(alpha = 0.7,
           width = 0.2) +
  facet_wrap(~frame) +
  scale_x_continuous(limits = c(0, 1.5)) +
  scale_fill_manual(values = c("#005EA0","#C0FF00")) +
  theme_classic() +
  theme(plot.title = element_blank(),
        strip.text.x = element_blank()) +
  guides(fill="none") +
  labs(x = "Ray hit height", y = "Proportion")

#10 for axis text, 12 for actual labels. 

p <- grid.arrange(plot1, plot2, ncol=1)

ggsave(
  "figure.png",
  plot = p,
  width = 4.25,
  height = 2.85,
  units = "in",
  dpi = 300
)



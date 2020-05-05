
#----- PREPARE DATA ------------------------------------------------------------

# load packages
library(lubridate)
library(ggplot2)
library(tidyr)
library(dplyr)

# would you like to save these figures?
save.figures <- TRUE
if (save.figures) {
  # specify directory for plot output
  plot.directory <- paste0(getwd(), "/data.out/Figures")
  if (!dir.exists(plot.directory)) {
    dir.create(plot.directory)
  }
}

# load four scenarios of data
load(file = "data.out/A.2014/salmon.finalstep.2014.1.RData")
baseline.data <- as.data.frame(salmon.finalstep)
rm(salmon.finalstep)
load(file = "data.out/A.2015/salmon.finalstep.2015.1.RData")
warm.data <- as.data.frame(salmon.finalstep)
rm(salmon.finalstep)
load(file = "data.out/P.2014/salmon.finalstep.2014.1.RData")
predator.data <- as.data.frame(salmon.finalstep)
rm(salmon.finalstep)
load(file = "data.out/P.2015/salmon.finalstep.2015.1.RData")
warm.predator.data <- as.data.frame(salmon.finalstep)
rm(salmon.finalstep)

# combine four scenarios into one dataframe
baseline.data <- cbind(scenario = "Baseline", baseline.data)
warm.data <- cbind(scenario = "Warm", warm.data)
predator.data <- cbind(scenario = "Predator", predator.data)
warm.predator.data <- cbind(scenario = "Warm-Predator", warm.predator.data)
all.scenarios.data <- rbind(baseline.data, warm.data, 
                            predator.data, warm.predator.data)
rm(baseline.data, warm.data, predator.data, warm.predator.data)

# prep data for ggplot2
all.scenarios.data <- all.scenarios.data %>% 
  select(scenario, survive, weight, dateSp, dateEm, dateOm) %>% 
  transmute(Scenario = as.factor(scenario),
            FinalState = as.factor(survive),
            Weight = as.numeric(weight),
            DateSpawn = date(as_datetime(dateSp, origin = "1970-01-01")),
            DateEmerge = date(as_datetime(dateEm, origin = "1970-01-01")),
            DateOutmigrate = date(as_datetime(dateOm, origin = "1970-01-01")))
levels(all.scenarios.data$FinalState) <- c("Predation", "Stochastic", "Yearling", "Subyearling")

#--- FIGURE 6: FINAL STATE -----------------------------------------------------

# Final state barplot
all.scenarios.data %>% 
  # filter data to surviving fish
  filter(FinalState == "Subyearling" | FinalState == "Yearling") %>%
  # plot final state vs. scenario
  ggplot(aes(x = Scenario, color = FinalState, fill = FinalState)) + 
  # add barplot
  geom_bar(alpha = 0.5) +
  # set theme
  theme_classic() +
  # remove legend title
  theme(legend.title = element_blank()) +
  # adjust y-axis label position
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
  # manually set color
  scale_color_manual(values=c("#008000", "#808080")) +
  # manually set fill
  scale_fill_manual(values=c("#008000", "#808080")) +
  # adjust y-axis label text
  labs(y = "Simulated\nsalmon\n(count)")
if (save.figures) {
  ggsave(path = plot.directory, filename = "Figure6.png", plot = last_plot(), 
         width = 5.5, height = 4, units = "in", dpi = 300)  
}

#--- FIGURE 5: MASS ------------------------------------------------------------

# Subyearling weight histogram
all.scenarios.data %>% 
  # filter data to surviving yearlings
  filter(FinalState == "Subyearling") %>%
  # plot weight
  ggplot(aes(x = Weight, color = FinalState, fill = FinalState)) + 
  # add histogram plot
  geom_histogram(alpha = 0.5) +
  # split plot by life history stage
  facet_wrap( ~ Scenario, nrow = 4, ncol = 1) +
  # set theme
  theme_classic() +
  # remove legend title
  theme(legend.position = "none") +
  # remove bottom axis line and ticks
  theme(axis.line.x.bottom = element_blank(), axis.ticks.x.bottom = element_blank()) +
  # remove subplot label background
  theme(strip.background = element_blank()) +
  # adjust y-axis label position
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
  # add axis and title text
  labs(x = "Subyearling mass (g)", y = "Simulated\nsalmon\n(count)") +
  # manually set color
  scale_color_manual(values="#808080") +
  # manually set fill
  scale_fill_manual(values="#808080") +
  # add x-axis to each plot
  geom_hline(yintercept = 0)
if (save.figures) {
  ggsave(path = plot.directory, filename = "Figure5Subyearling.png", plot = last_plot(), 
         width = 4, height = 6, units = "in", dpi = 300)
}

# Yearling weight histogram
all.scenarios.data %>% 
  # filter data to surviving yearlings
  filter(FinalState == "Yearling") %>%
  # plot weight
  ggplot(aes(x = Weight, color = FinalState, fill = FinalState)) + 
  # add histogram plot
  geom_histogram(alpha = 0.5) +
  # split plot by life history stage
  facet_wrap( ~ Scenario, nrow = 4, ncol = 1) +
  # set theme
  theme_classic() +
  # remove legend
  theme(legend.position = "none") +
  # remove bottom axis line and ticks
  theme(axis.line.x.bottom = element_blank(), axis.ticks.x.bottom = element_blank()) +
  # remove subplot label background
  theme(strip.background = element_blank()) +
  # adjust y-axis label position
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
  # add axis and title text
  labs(x = "Yearling mass (g)", y = "Simulated\nsalmon\n(count)") +
  # manually set color
  scale_color_manual(values="#008000") +
  # manually set fill
  scale_fill_manual(values="#008000") +
  # add x-axis to each plot
  geom_hline(yintercept = 0)
if (save.figures) {
  ggsave(path = plot.directory, filename = "Figure5Yearling.png", plot = last_plot(), 
         width = 4, height = 6, units = "in", dpi = 300)
}


#--- FIGURE 4: PHENOLOGY -------------------------------------------------------

# change year in warm scenarios to match year in cool scenarios,
# so that day-month is the comparable across scenarios.
phenology.data <- all.scenarios.data
scenario.column <- all.scenarios.data$Scenario
spawn.column <- all.scenarios.data$DateSpawn
emergence.column <- all.scenarios.data$DateEmerge
outmigration.column <- all.scenarios.data$DateOutmigrate
# this is slow, but troubleshooting a vectorized way was wasting time
for (i in 1:length(scenario.column)) {
  if (scenario.column[i] == "Warm" | scenario.column[i] == "Warm-Predator") {
    date.spawn <- spawn.column[i]
    year(date.spawn) <- (year(date.spawn) - 1)
    spawn.column[i] <- date.spawn
    
    date.emerge <- emergence.column[i]
    year(date.emerge) <- (year(date.emerge) - 1)
    emergence.column[i] <- date.emerge
    
    date.outmigrate <- outmigration.column[i]
    year(date.outmigrate) <- (year(date.outmigrate) - 1)
    outmigration.column[i] <- date.outmigrate
  }
}

# plot emergence and outmigration across scenarios
phenology.data %>% 
  mutate(Emergence = emergence.column, Outmigration = outmigration.column) %>% 
  # filter data to surviving fish
  filter(FinalState == "Subyearling" | FinalState == "Yearling") %>% 
  # select relevant columns
  select(Scenario, Emergence, Outmigration) %>% 
  # combine emergence and outmigration into a single column
  gather(key = "Event", value = "Date", c(Emergence, Outmigration)) %>% 
  # filter data to fish that experienced each event
  filter(!is.na(Date)) %>% 
  # plot event vs. day-month
  ggplot(aes(x = Date, fill = Event, color = Event)) + 
  # add histogram
  geom_histogram(alpha = 0.5, position = "identity") + 
  # split plot by scenario
  facet_wrap(~ Scenario, nrow = 4, ncol = 1) +
  # set theme
  theme_classic() +
  # adjust y-axis label position
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
  # remove legend title
  theme(legend.title = element_blank()) +
  # remove bottom axis line and ticks
  theme(axis.line.x.bottom = element_blank(), axis.ticks.x.bottom = element_blank()) +
  # remove subplot label background
  theme(strip.background = element_blank()) +
  # manually set color
  scale_color_manual(values=c("#00008B", "#FF8C00")) +
  # manually set fill
  scale_fill_manual(values=c("#00008B", "#FF8C00")) +
  # set x-axis labels to day-month
  scale_x_date(date_labels = "%d-%b") +
  # adjust y-axis label text
  labs(y = "Simulated\nsalmon\n(count)") +
  # add x-axis to each plot
  geom_hline(yintercept = 0)
if (save.figures) {
  ggsave(path = plot.directory, filename = "Figure4.png", plot = last_plot(), 
         width = 6.5, height = 5.5, units = "in", dpi = 300)
}


#' Estimating the rate at which individuals drop an antibody dilution titre level
#' Using data from a point source outbreak (Lupidi et al) 
#' Collected data over five time points (a-e) over 5-year period from 18 individuals
#' Three different serovars
#' E Rees
#' 

library(tidyverse)
library(lme4)

## Read in data
dat <- read_csv("data/lupidi_full_data.csv")

## process data to find the max between a and b and the lowest titre time (individuals undetectable at different times)

dat <- dat %>%
  rowwise() %>% 
  mutate(highest_titre = max(a,b)) %>% 
  mutate(highest_titre_time = ifelse(b > a, 9,1))

dat <- dat %>% 
  mutate(lowest_titre = ifelse(!is.na(e), e, 
                               ifelse(!is.na(d), d, 
                                      ifelse(!is.na(c),c,b)))) %>% 
  mutate(lowest_titre_time = ifelse(!is.na(e), 54, 
                                    ifelse(!is.na(d), 36, 
                                           ifelse(!is.na(c),18,9))))

## Look at three different serovars independently
## 1. Bratislava

dat_bratislava_all <- dat %>% 
  filter(Serovar == "Bratislava") %>% 
  select(1:6)

dat_bratislava_all <- dat_bratislava_all %>%
  gather(key = "time", value = "logTitre", 2:6) %>% 
  mutate(time = ifelse(time == "a", 1, 
                       ifelse(time == "b", 9,
                              ifelse(time == "c" , 18, 
                                     ifelse(time == "d", 36, 54)))))


m1 <- lmer(logTitre ~ time + (1|ID), data = dat_bratislava_all)
summary(m1)
plot(m1)
confint(m1)

fig1 <- ggplot(data=dat_bratislava_all, aes(x=time, y=logTitre, group=ID, colour=ID)) +
  geom_line() +
  geom_abline(aes(slope=-0.043546,
                  intercept=3.151250), size=.8)+
  theme(legend.position="none") +
  labs(x="time (months)", y="Scaled Titer (log10)")

fig1

## calculating the time it takes to drop from 1:100 to 1:50 (log base 10) - since these are present in both datasets
calculateAbDrop <- function(slope,intercept,y1=2,y2=1.69897){
  x1 = (y1-intercept)/slope
  x2 =  (y2-intercept)/slope
  x1 - x2
}

calculateAbDrop(slope=-0.043546,intercept = 3.151250)

# upper CI
intercept = 3.48584385
slope = -0.03325409

# lower CI
intercept = 2.81748360
slope = -0.05343897


## 1. Australis

dat_australis_all <- dat %>% 
  filter(Serovar == "Australis") %>% 
  select(1:6)

dat_australis_all <- dat_australis_all %>%
  gather(key = "time", value = "logTitre", 2:6) %>% 
  mutate(time = ifelse(time == "a", 1, 
                       ifelse(time == "b", 9,
                              ifelse(time == "c" , 18, 
                                     ifelse(time == "d", 36, 54)))))


m2 <- lmer(logTitre ~ time + (1|ID), data = dat_australis_all)
summary(m2)
plot(m2)
confint(m2)


fig2 <- ggplot(data=dat_australis_all, aes(x=time, y=logTitre, group=ID, colour=ID)) +
  geom_line() +
  geom_abline(aes(slope=-0.032357,
                  intercept=2.584991), size=.8)+
  theme(legend.position="none") +
  labs(x="time (months)", y="Scaled Titer (log10)")

fig2

calculateAbDrop(slope=-0.032357,intercept = 2.584991)

# upper
intercept = 2.90335580
slope = -0.01996802

# lower
intercept = 2.26952936
slope = -0.04372622

## 3. Iora

dat_iora_all <- dat %>% 
  filter(Serovar == "Iora") %>% 
  select(1:6)

dat_iora_all <- dat_iora_all %>%
  gather(key = "time", value = "logTitre", 2:6) %>% 
  mutate(time = ifelse(time == "a", 1, 
                       ifelse(time == "b", 9,
                              ifelse(time == "c" , 18, 
                                     ifelse(time == "d", 36, 54)))))


m3 <- lmer(logTitre ~ time + (1|ID), data = dat_iora_all)
summary(m3)
plot(m3)
confint(m3)

# Very similar model results, could look at AIC to see which is better

fig3 <- ggplot(data=dat_iora_all, aes(x=time, y=logTitre, group=ID, colour=ID)) +
  geom_line() +
  geom_abline(aes(slope=-0.040110,
                  intercept=3.367196), size=.8)+
  theme(legend.position="none") +
  labs(x="time (months)", y="Scaled Titer (log10)")

fig3

calculateAbDrop(slope=-0.040110,intercept = 3.367196)

# upper
intercept = 3.65692042
slope = -0.03301864

# lower
intercept = 3.0779675
slope = -0.0471598

##Antibody titre drops 

# Bratislava 6.937 (5.633 - 9.052)
# Australis 9.303 (6.884 - 15.076)
# Iora 7.505 (6.383191- 9.116972)

mean(c(6.937,9.303,7.505))
mean(c(5.633,6.884,6.383191))
mean(c(9.052,15.076,9.116972))


library(tidyverse)
install.packages("ggpubr")
library(ggpubr)
install.packages("dplyr")
library(dplyr)
library(tidylog)
library(ggplot2)
library(aod)
library(caret)
library(pROC)
library(MuMIn)


## model evaluation functions
evaluation = function(fit_model,df){
  
  prbs = predict(fit_model,type='response')
  df$predicted_class <- ifelse(prbs > 0.5, 1, 0)
  df$predicted_pbs <- prbs
  
  
  df$event <- factor(df$event)
  print(df$event)
  
  df$predicted_class <- factor(df$predicted_class, levels = levels(df$event))
  conf_matrix <- confusionMatrix(df$predicted_class, df$event)
  
  roc_curve <- roc(df$event, df$predicted_pbs)
  plot(roc_curve)
  print(auc(roc_curve))
  return (conf_matrix)
}

get_predict = function(fit_model,df){
  
  prbs = predict(fit_model,type='response')
  df$predicted_class <- ifelse(prbs > 0.5, 1, 0)
  df$predicted_pbs <- prbs
  return(df)
}


goodft<- function(tst){
  return (1-pchisq(tst$deviance,tst$df.residual))
}

## data preprocessing
cd = read.csv("final_cardiac_data.csv", header = TRUE)
cd = mutate(cd, ethnic2 = ifelse(ethnic1 == 1, "Other", 
                                 ifelse(ethnic1 == 2, "Other",
                                        ifelse(ethnic1 == 3, "White",
                                               ifelse(ethnic1 == 4, "Black",
                                                      ifelse(ethnic1 == 5, "Other", NA)))))) %>%
mutate(diabetes = ifelse(diabetes == 9, NA, 
                         ifelse(diabetes == 7, NA, diabetes))) %>%
mutate(smoker = ifelse(smoker == 9, NA,
                       ifelse(smoker == 7, NA, smoker))) %>%
mutate(diab.new = ifelse(diabetes == 2, 0,
                         ifelse(diabetes == 3, 1, 
                                ifelse(diabetes == 1, 2, NA)))) %>%
mutate(educ = ifelse(educ == 9, NA,
                     ifelse(educ == 7, NA, educ))) %>%
#new coding for educ.new: 0 = <9 + 9-11th grade, 1 = hs/ged, 2= some college, 3 = college graduate)
mutate(educ.new = ifelse(educ == 1, 0, 
                         ifelse(educ == 2, 0, 
                                ifelse(educ == 3, 1, 
                                       ifelse(educ == 4, 2,
                                              ifelse(educ == 5, 3, NA))))))




cd.nn = subset(cd, is.na(bmi) == FALSE & is.na(smoker) == FALSE & 
                 is.na(diabetes) == FALSE & is.na(educ) == FALSE & is.na(sleep.hrs) == FALSE)
mean = mean(cd.nn$bmi)
sdev = sd(cd.nn$bmi)
cd.nn = subset(cd.nn, bmi <= mean+3.5*sdev & bmi >= mean-3.5*sdev)
count_tb = table(event=cd.nn$event,col=cd.nn[['educ']])
print(count_tb)
print(prop.table(count_tb,margin = 2))
count_tb_df <- as.data.frame(count_tb)
 


#create subset of data that's the tossed out subjects w/ NA BMI (for descriptive analysis)
cd.na = subset(cd, is.na(bmi) == TRUE | is.na(smoker) == FALSE |
                 is.na(diabetes) == FALSE | is.na(educ) == FALSE | bmi > mean + 3.5*sdev)

cd.nn$educ.new <- factor(cd.nn$educ.new, levels = c(0,1,2,3),labels = c("<=-11th Grade",
                                                                        "High school grad/GED","Some college or AA degree","College graduate or above"))

cd.nn$age_group <- cut(cd.nn$age, breaks = c(0, 20, 30, 40, 50,60,70,80), 
                       labels = c("0-20", "21-30", "31-40", "41-50",'51-60','61-70','71-80'), 
                       include.lowest = TRUE)

cd.nn$bmi_group <- cut(cd.nn$bmi, breaks = c(10, 20, 30, 40, 50,60), 
                       labels = c("10-20", "21-30", "31-40", "41-50",'51-60'), 
                       include.lowest = TRUE)

cd.nn$sleep.hrs_group <- cut(cd.nn$sleep.hrs, breaks = c(0,3, 6, 9, 14), 
                       labels = c( '0-3',"3-6", "6-9",'9+'), 
                       include.lowest = TRUE)

# Use ggplot2 to create a bar plot
ggplot(count_tb_df, aes(x = count_tb_df$col, y = count_tb_df$Freq, fill = count_tb_df$event)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Educ", x = "Education Level", y = "Count") +
  scale_fill_manual(values = c("lightblue", "lightgreen"), name = "Event Case") +
  theme_minimal()


#create subset of data that tosses out all subjects w/ NAs for BMI
ggplot(cd.nn, aes(x = "", y = bmi)) +
  geom_boxplot(fill = "lightblue", color = "blue") +
  geom_point(size=4) +
  labs(title = "BMI Box Plot", x = "BMI") +
  theme_minimal()

ggplot(cd.nn, aes(x = "", y = sleep.hrs)) +
  geom_boxplot(fill = "lightblue", color = "blue") +
  labs(title = "AGE Box Plot", x = "age") +
  theme_minimal()


cate_data <-c("gender"  ,  "age"    ,   "ethnic1"  , "educ"      ,"sleep.hrs" ,"diabetes" , "smoker"  ,  "bmi"     ,  "ethnic2"  , "diab.new" ,"educ.new"  ,"age_group" ,
              "bmi_group")
get_mode<-function(colname){
  
  frequency_table <- table(cd.nn[[colname]])
  
  modes <-   names(frequency_table)[frequency_table == max(frequency_table)] 
 
  return(modes)
}
 

## slice polot for sleep hr,
gp_count <- cd.nn %>% group_by(sleep.hrs_group,event) %>% count()
gp_count %>% pivot_wider(names_from = event, values_from = n, names_prefix = "event_", values_fill = 0)

widen <- gp_count %>% pivot_wider(names_from = event, values_from = n, names_prefix = "event_", values_fill = 0)
widen$sum_ct <- widen$event_0 + widen$event_1
widen$ratio <- widen$event_1/(widen$event_0 + widen$event_1)
widen$logodds <- log((widen$ratio+0.01)/(1-widen$ratio+0.01))
ggplot(widen, aes(x = sleep.hrs_group, y = logodds)) +
  geom_point(size=5) +
  labs(title = "Plot of log odds by sleep_hrs group", x = "sleep_hrs Group", y = "logit")


## make sliced plots
gp_count <- cd.nn %>% group_by(age_group,event) %>% count()
gp_count %>% pivot_wider(names_from = event, values_from = n, names_prefix = "event_", values_fill = 0)

widen <- gp_count %>% pivot_wider(names_from = event, values_from = n, names_prefix = "event_", values_fill = 0)
widen$sum_ct <- widen$event_0 + widen$event_1
widen$ratio <- widen$event_1/(widen$event_0 + widen$event_1)
widen$logodds <- log((widen$ratio+0.01)/(1-widen$ratio+0.01))
library(ggplot2)
# Create a scatter plot with points

ggplot(widen, aes(x = age_group, y = logodds)) +
  geom_point(size=5) +
  labs(title = "Plot of log odds by age_group", x = "Age Group", y = "logit")

gp_count <- cd.nn %>% group_by(bmi_group,event) %>% count()
gp_count %>% pivot_wider(names_from = event, values_from = n, names_prefix = "event_", values_fill = 0)

widen <- gp_count %>% pivot_wider(names_from = event, values_from = n, names_prefix = "event_", values_fill = 0)
widen$sum_ct <- widen$event_0 + widen$event_1
widen$ratio <- widen$event_1/(widen$event_0 + widen$event_1)
widen$logodds <- log((widen$ratio+0.01)/(1-widen$ratio+0.01))
library(ggplot2)
# Create a scatter plot with points
ggplot(widen, aes(x = bmi_group, y = logodds)) +
  geom_point(size=5) +
  labs(title = "Plot of log odds by bmi_group", x = "BMI Group", y = "logit")



## model selection
first <- glm(event~gender+age+ethnic2+educ.new+diab.new+smoker+sleep.hrs+bmi,data=cd.nn,family = binomial)
anova(first)
summary(first)

test.f = glm(event ~ age+bmi+diab.new+ethnic2,data = cd.nn, family = binomial)
anova(test.f)
summary(test.f)

tst1 <- glm(event~age+diab.new+ethnic2+bmi+gender,data=cd.nn,family = binomial)
anova(tst1)
summary(tst1)

tst2 <- glm(event~age+diab.new+ethnic2+bmi+gender+educ.new,data=cd.nn,family = binomial)
anova(tst2)
summary(tst2)

tst3 <- glm(event~age+diab.new+ethnic2+bmi+gender+educ.new+smoker,data=cd.nn,family = binomial)
anova(tst3)
summary(tst3)

tst4 <- glm(event~age+diab.new+ethnic2+bmi+gender+educ.new+smoker+sleep.hrs,data=cd.nn,family = binomial)
anova(tst4)
summary(tst4)

## call autonation fitting function
set
pdredge(final)

#gender educ.new bmi ethnic2 smoker sleep.hrs diab.new age
final = glm(event ~ gender*educ.new + gender*bmi + gender*ethnic2 + age*ethnic2 + educ.new*smoker + ethnic2*smoker + sleep.hrs + diab.new,data = cd.nn, family = binomial)
summary(final)

final.test = glm(event ~ gender + age + bmi + ethnic2 + diab.new + educ.new + smoker, data = cd.nn, family = binomial)
## call evaluation function

evaluation(test.f,cd.nn)
evaluation(tst1,cd.nn)
evaluation(tst2,cd.nn)
evaluation(tst3,cd.nn)
evaluation(tst4,cd.nn)
evaluation(final,cd.nn)
evaluation(final.test,cd.nn)

goodft(tst1)
goodft(tst2)
goodft(tst3)
goodft(tst4)
goodft(final)


predicted_df <- get_predict(final,cd.nn)
plot(x=predicted_df$age,y=predicted_df$predicted_pbs)





final <- glm(event ~ gender*educ.new + gender*bmi + gender*ethnic2 + age*ethnic2 + educ.new*smoker + ethnic2*smoker + sleep.hrs + diab.new, data = cd.nn, family = binomial)

get_mode<-function(colname){
  
  frequency_table <- table(cd.nn[[colname]])
  
  modes <- names(frequency_table)[frequency_table == max(frequency_table)] 
  print(modes)
  return(modes)
}
get_mode('smoker')


# Create new data for prediction
new_data <- data.frame(
  gender = c(rep(get_mode('gender'),70)),   
  educ.new = c(rep(get_mode('educ.new'),70)),   
  bmi =  c(rep(median(cd.nn$bmi),70)),      
  ethnic2 = c(rep(get_mode('ethnic2'),70)),   
  age = c(rep(median(cd.nn$age),70)),   
  smoker = c(rep(get_mode('smoker'),70)), 
  sleep.hrs = c(rep(median(cd.nn$sleep.hrs),70)),    
  diab.new =  c(rep(get_mode('diab.new'),70)), 
)


# Assuming you have x values defined
x <- seq(10, 80, length.out = 1000)
logit_exprdia <- function(x, gender) {
  exp(
    intercept +
      beta_gender * gender +
      beta_age * x +
      beta_bmi * 29.1 +
      beta_ethnic2Other * 1 +
      beta_ethnic2White * 0 +
      beta_diab_new * 0+
      beta_educ_new * 2 +
      beta_smoker ) / (1 +   exp(
        intercept +
          beta_gender * gender +
          beta_age * x +
          beta_bmi * 29.1 +
          beta_ethnic2Other * 1 +
          beta_ethnic2White * 0 +
          beta_diab_new * 0+
          beta_educ_new * 2 +
          beta_smoker ))
}
# Plotting the first curve
curve(logit_exprdia(x, 0), from = 10, to = 80, n = 1000, col = "blue", ylab = "Probability", xlab = "Age")

# Adding the second curve
lines(from = 10, to = 80,logit_exprdia(x, 1), col = "red")



curve(logit_exprgen(x, 0), from = 10, to = 80, n = 1000, col = "blue", ylab = "Probability", xlab = "Age")
# Adding the second curve
lines(x, logit_exprgen(x, 1), col = "red")
legend("topright", legend = c("female", "male"), col = c("red", "blue"), lty = 1)
 



intercept <- -4.728681
beta_gender <- -0.110615
beta_age <- 0.065176
beta_bmi <- 0.062702
beta_ethnic2Other <- -0.598345
beta_ethnic2White <- -0.567687
beta_diab_new <- 0.414654
beta_educ_new <- 0.017423
beta_smoker <- -0.304745

logit_exprgen <- function(x, gen) {
  exp(
    intercept +
      beta_gender * gen +
      beta_age * x +
      beta_bmi * 29.1 +
      beta_ethnic2Other * 1 +
      beta_ethnic2White * 0 +
      beta_diab_new * 0+
      beta_educ_new * 2 +
      beta_smoker )/ (1 +   exp(
    intercept +
      beta_gender * gen +
      beta_age * x +
      beta_bmi * 29.1 +
      beta_ethnic2Other * 1 +
      beta_ethnic2White * 0 +
      beta_diab_new * 0+
      beta_educ_new * 2 +
      beta_smoker ))
}
df <- data.frame(
  x = seq(10, 80, length.out = 1000),
  probability_female = logit_exprgen(df$x, 0),
  probability_male = logit_exprgen(df$x, 1)
)

# Reshape the data to long format
df_long <- tidyr::gather(df, key = "gender", value = "probability", -x)

# Plotting with ggplot2
ggplot(df_long, aes(x = x, y = probability, color = gender)) +
  geom_line() +
  labs(
    title = "Probability Subplot with Age by Gender",
    x = "Age",
    y = "Probability"
  ) +
  scale_color_manual(values = c("red", "blue")) 




logit_exprdia <- function(x, diabetes) {
  exp(
    intercept +
      beta_gender * 1 +
      beta_age * x +
      beta_bmi * 29.1 +
      beta_ethnic2Other * 1 +
      beta_ethnic2White * 0 +
      beta_diab_new * diabetes+
      beta_educ_new * 2 +
      beta_smoker ) / (1 +   exp(
    intercept +
      beta_gender * 1 +
      beta_age * x +
      beta_bmi * 29.1 +
      beta_ethnic2Other * 1 +
      beta_ethnic2White * 0 +
      beta_diab_new * diabetes+
      beta_educ_new * 2 +
      beta_smoker ))
}

df <- data.frame(
  x = seq(10, 80, length.out = 1000),
  no_diaetes = logit_exprdia(seq(10, 80, length.out = 1000), 0),
  borderline = logit_exprdia(seq(10, 80, length.out = 1000), 1),
  wdiaetes = logit_exprdia(seq(10, 80, length.out = 1000), 2)
)

# Convert data to long format for ggplot2
df_long <- tidyr::gather(df, key = "diabetes_status", value = "probability", -x)

# Plot using ggplot2
ggplot(df_long, aes(x = x, y = probability, color = diabetes_status)) +
  geom_line() +
  labs(
    title = "Probability Subplot with Age by Diabetes Status",
    x = "Age",
    y = "Probability"
  ) +
  scale_color_manual(values = c("blue", "red", "green"))
 

 
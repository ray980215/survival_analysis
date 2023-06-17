library(survival)
library(survminer)
library(psych)
library(ggplot2)

setwd("C:/Users/ray98/Desktop/class/survival analysis/final project")
data <- read.csv("METABRIC_RNA_Mutation.csv")
data <- data[, c(1, 2, 7, 16, 24, 25, 27, 29, 30, 31)]
data <- data[complete.cases(data),]
data <- data[-(data$death_from_cancer==""),]
data$cs <- ifelse(data$death_from_cancer == "Died of Disease", 1, 0)


#EDA#######################
#KM
fitkm <- survfit(Surv(overall_survival_months, cs) ~ 1, data = data)
ggsurvplot(fitkm, data = data,  xlab = "Months")

#X summary
p1 <- ggplot(data, aes(x=overall_survival_months)) + geom_histogram(binwidth=20, color="black", fill="white")
p2 <- ggplot(data, aes(y=overall_survival_months)) + geom_boxplot()
p3 <- ggplot(data, aes(x=cs)) + geom_histogram(binwidth=1, color="black", fill="white") + xlab("censoring status")
p4 <- ggplot(data, aes(x=age_at_diagnosis)) + geom_histogram(binwidth=5, color="black", fill="white")
p5 <- ggplot(data, aes(y=age_at_diagnosis)) +  geom_boxplot()
p6 <- ggplot(data, aes(x=tumor_stage)) + geom_histogram(binwidth=1, color="black", fill="white")
p7 <- ggplot(data, aes(x=chemotherapy)) + geom_histogram(binwidth=1, color="black", fill="white")
p8 <- ggplot(data, aes(x=hormone_therapy)) + geom_histogram(binwidth=1, color="black", fill="white")
p9 <- ggplot(data, aes(x=radio_therapy)) + geom_histogram(binwidth=1, color="black", fill="white")
p10 <- ggplot(data, aes(x=tumor_size)) + geom_histogram(binwidth=10, color="black", fill="white")
ggarrange(p1, p2)
ggarrange(p3, p4)
ggarrange(p5, p6)
ggarrange(p7, p8)
ggarrange(p9, p10)
#MH test&Gehan's test
myfun <- function(data, i){
  name <- colnames(data)[i]
  fmla <- as.formula(paste("Surv(overall_survival_months, cs) ~ ", name))
  surv_diff <- survdiff(fmla, data = data)
  print(surv_diff)
  fit1 <- surv_fit(fmla, data = data)
  print(ggsurvplot(fit1, pval=T, data = data))
  surv_diff <- survdiff(fmla, data = data, rho = 1)
  print(surv_diff)
}
myfun(data, 3)
myfun(data, 4)
myfun(data, 7)
myfun(data, 9)
data1 <- data[data$tumor_stage == 0|data$tumor_stage == 1,]
data2 <- data[data$tumor_stage == 2,]
data3 <- data[data$tumor_stage == 3|data$tumor_stage == 4,]
myfun(data1, 3)
myfun(data2, 3)
myfun(data3, 3)
myfun(data1, 4)
myfun(data2, 4)
myfun(data3, 4)
myfun(data1, 7)
myfun(data2, 7)
myfun(data3, 7)
myfun1 <- function(i){
  splots <- list()
  name <- colnames(data)[i]
  fmla <- as.formula(paste("Surv(overall_survival_months, cs) ~ ", name))
  fit1 <- surv_fit(fmla, data = data1)
  splots[[1]] <- ggsurvplot(fit1, pval=T, data = data1)
  fit1 <- surv_fit(fmla, data = data2)
  splots[[2]] <- ggsurvplot(fit1, pval=T, data = data2)
  fit1 <- surv_fit(fmla, data = data3)
  splots[[3]] <- ggsurvplot(fit1, pval=T, data = data3)
  arrange_ggsurvplots(splots, print = TRUE,
                      ncol = 3, nrow = 1)
}
myfun1(3)
myfun1(4)
myfun1(7)











#Data Analysis##############
#Cox

fitcox <- coxph(Surv(overall_survival_months, cs) ~  hormone_therapy + chemotherapy + radio_therapy + tumor_stage + tumor_size + age_at_diagnosis , data = data)
summary(fitcox)

csr <- data$cs - fitcox$residuals
csmodel <- survfit(Surv(csr, data$cs) ~ 1)
ggplot() +
  geom_point(aes(x = csmodel$time, 
                 y = -log(csmodel$surv))) +
  labs(title = "Cox-snell residuals versus its cumulative hazard function",
       x = expression(r[Ci] ),
       y = expression(-log(S[r[C]](r[Ci])) )) +
  geom_abline(intercept = 0, slope = 1)
cox.zph(fitcox, transform="km", global=TRUE)


#AFT######################################
#lognormal
fit_ln <- survreg(Surv(overall_survival_months, cs) ~ hormone_therapy + chemotherapy + radio_therapy + tumor_stage + tumor_size + age_at_diagnosis, data = data, dist = "lognormal")
summary(fit_ln)
#loglogistic
fit_ll <- survreg(Surv(overall_survival_months, cs) ~ hormone_therapy + chemotherapy + radio_therapy + tumor_stage + tumor_size + age_at_diagnosis, data = data, dist = "loglogistic")
summary(fit_ll)

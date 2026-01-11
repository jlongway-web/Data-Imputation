library(car)
library(ggplot2)
library(GGally)
library(dplyr)
library(modeest)
library(olsrr)
library(kernlab)
library(pROC)

set.seed(9)
url <- "https://archive.ics.uci.edu/ml/machine-learning-databases/breast-cancer-wisconsin/breast-cancer-wisconsin.data"
cancer <- read.table(url, sep = ",", header = FALSE, stringsAsFactors = FALSE)

head(cancer)
#predict V11
sapply(cancer, class)
#according to data source, V1 is an ID number, Class is Binary, and v7 is integer
cancer <- cancer %>%
  mutate(across(V11, ~  . /2))
cancer <- cancer %>%
  mutate(across(V11, ~  . -1))
#Class 0 benign, 1 malignant for logistic regression if needed

cor(cancer$V1, cancer$V11)
#We can see the V1 does not have have any correlation value as an integer
#since character/integer correlation does not work, we will run logistic
#regression models with and without V1 for comparison
#these models
cancer$V1 <- as.character(cancer$V1)
model <- glm(V11 ~., family="binomial"(link="logit"), data=cancer)
summary(model)

#Null deviance: 9.0053e+02  on 698  degrees of freedom
#Residual deviance: 4.1015e-09  on  38  degrees of freedom
#AIC: 1322 - Number of Fisher Scoring iterations: 25

cancer <- cancer[,-1]

model <- glm(V11 ~., family="binomial"(link="logit"), data=cancer)
summary(model)

#regression model ran to compare with and without V1 as character to validate
#non-usefulness of V1, AIC is substantially better 1322 -> 139.5
#this model still has V7 as character but is sufficient for our purpose

#not worth using outside of labeling data points, as such, we will not use it as
#an integral part of any models


#splitting up dataset into missing and non-missing



missing <- subset(cancer, V7 == "?")
notm <- subset(cancer, V7 != "?")
#there are 16 observations with missing values, approx 2%, which is
#well below the recommended 5% threshold

notm <- type.convert(notm, as.is = TRUE)
sapply(notm, class)

mean <- mean(notm$V7)
mode <- mfv(notm$V7)[1]

cor(notm$V7, notm$V11)
#after V7 is converted to an integer, it shows a 0.8226959 correlation with 
#V11 with the missing values excluded and is useful for modeling

#creating training and testing dataset with 20/80 split, missing data was
#excluded and readded to ensure it was part of the training set, and
#split values were adjusted to account for the extra values.

# Define descriptive column names
column_names <- c("Clump_Thickness", "Uniformity_Cell_Size", 
                  "Uniformity_Cell_Shape", "Marginal_Adhesion", 
                  "Single_Epithelial_Cell_Size", "Bare_Nuclei", 
                  "Bland_Chromatin", "Normal_Nucleoli", "Mitoses", "Class")

# Apply names to the dataframe
colnames(cancer) <- column_names
colnames(notm) <- column_names
colnames(missing) <- column_names




ind <- sample(2, nrow(cancer), replace = TRUE, prob = c(.2,.80))
#13 ? in training, 3 in testing, close to 20/80 split
ind1 <- sample(2, nrow(notm), replace = TRUE, prob = c(.2,.80))

#-----------------------------------input model using regression--------------  

inputmodel <- lm(formula = Bare_Nuclei ~ Clump_Thickness + Uniformity_Cell_Size + Uniformity_Cell_Shape + Marginal_Adhesion + Single_Epithelial_Cell_Size + Bland_Chromatin + Normal_Nucleoli + Mitoses, data = notm)
summary(inputmodel) # 0.6104 

inputfor <- ols_step_forward_aic(inputmodel)
inputfor #AIC =  3067.279 R2 = 0.611
plot(inputfor)

inputfinal <- lm(formula = Bare_Nuclei ~ Clump_Thickness + Uniformity_Cell_Shape + Marginal_Adhesion + Bland_Chromatin, data = notm)
summary(inputfinal)
inputfinal

reginput <- cancer
for(x in 1:nrow(reginput))
  {
  if(reginput$Bare_Nuclei[x] == "?")
    {
    reginput$Bare_Nuclei[x] = round(-0.53601
    + reginput$Clump_Thickness[x]*0.22617
    + reginput$Uniformity_Cell_Shape[x]*0.31729
    + reginput$Marginal_Adhesion[x]*0.33227
    + reginput$Bland_Chromatin[x]*0.32378
    )}
}


reginput$Bare_Nuclei <- as.integer(reginput$Bare_Nuclei)
cor(reginput$Bare_Nuclei,reginput$Class)
#After regression input, Bare_Nuclei correlation (notm) changed from 0.8226959 to 0.8193201

regtrain <- reginput[ind==2,]
regtest <- reginput[ind==1,]

#building a second test data set with an additional column to test models
#that add an additional column for Bare_Nuclei


#-----------------------------------input model using perturbation-------------  

perinput <- missing
head(perinput)
for(x in 1:nrow(perinput))
{
  perinput$Bare_Nuclei[x] = (-0.53601
                           + perinput$Clump_Thickness[x]*0.22617
                           + perinput$Uniformity_Cell_Shape[x]*0.31729
                           + perinput$Marginal_Adhesion[x]*0.33227
                           + perinput$Bland_Chromatin[x]*0.32378)
}

perinput$Bare_Nuclei <- as.numeric(perinput$Bare_Nuclei)
cor(perinput$Bare_Nuclei,perinput$Class)
#for imputed values, correlation to Class is significantly lower at 0.5346644
sd <- sd(notm$Bare_Nuclei)
plot(density(notm$Bare_Nuclei))
perturb <- as.data.frame(rnorm(nrow(perinput),0, sd))
perturb <- data.frame(perturb)
perturb <- data.frame(round(perturb + perinput$Bare_Nuclei))
#Values do not exceed 10, but some are negative and will be adjusted to 0 to
#match value range
perturb[perturb<0] <- 0
colnames(perturb)[1] <- "Bare_Nuclei"
perturb <- unlist(perturb)

perinput$Bare_Nuclei <- data.frame(perturb)
perinput$Bare_Nuclei <- unlist(perinput$Bare_Nuclei)

fullper <- data.frame(rbind(perinput, notm))

pertrain <- fullper[ind==2,]
pertest <- fullper[ind==1,]


#-----------------------------------input model using mean--------------  

meaninput <- cancer
for(x in 1:nrow(meaninput))
{
  if(meaninput$Bare_Nuclei[x] == "?")
  {
    meaninput$Bare_Nuclei[x] = mean}
}

meaninput$Bare_Nuclei <- as.integer(meaninput$Bare_Nuclei)
sapply(meaninput, class)
cor(meaninput$Bare_Nuclei,meaninput$Class)
#cor  0.8226959 to 0.8174423

meantrain <- meaninput[ind==2,]
meantest <- meaninput[ind==1,]


#-----------------------------------input model using mode--------------  

modeinput <- cancer
for(x in 1:nrow(modeinput))
{
  if(modeinput$Bare_Nuclei[x] == "?")
  {
    modeinput$Bare_Nuclei[x] = mode}
}

modeinput$Bare_Nuclei <- as.integer(modeinput$Bare_Nuclei)
sapply(modeinput, class)
cor(modeinput$Bare_Nuclei,modeinput$Class)
#cor  0.8226959 to  0.8189679

modetrain <- modeinput[ind==2,]
modetest <- modeinput[ind==1,]

#-----------------------------------input model using binary value-------------  

bininput <- cancer
bininput$Bare_Nucleia <- 0
for(x in 1:nrow(bininput))
{
  if(bininput$Bare_Nuclei[x] == "?")
  {
    bininput$Bare_Nuclei[x] = 0
    bininput$Bare_Nucleia[x] = 1}
}

sapply(bininput, class)
bininput$Bare_Nuclei <- as.integer(bininput$Bare_Nuclei)
bininput <- bininput[, c("Clump_Thickness", "Uniformity_Cell_Size", "Uniformity_Cell_Shape", "Marginal_Adhesion", "Single_Epithelial_Cell_Size", "Bare_Nuclei", "Bare_Nucleia",
                         "Bland_Chromatin", "Normal_Nucleoli", "Mitoses", "Class")]

bintrain <- bininput[ind==2,]
bintest <- bininput[ind==1,]

#------------not missing input-------------------------


notmtrain <- notm[ind1==2,]
notmtest <- notm[ind1==1,]


#-----------------------SVM----------------------------------------------


#regression input--------------------------------

svmreginput <- ksvm(as.matrix(regtrain[,1:9]),as.factor(regtrain[,10]),
                    type="C-svc",kernel="vanilladot",C=0.1,scaled=TRUE)
a <- colSums(svmreginput@xmatrix[[1]] * svmreginput@coef[[1]])
print(svmreginput@coef[[1]])
a
a0 <- -svmreginput@b
a0
pred <- predict(svmreginput,regtrain[,1:9])
pred
sum(pred == regtrain[,10]) / nrow(regtrain)
#0.9711712
pred <- predict(svmreginput,regtest[,1:9])
sum(pred == regtest[,10]) / nrow(regtest)
# 0.9583333

#mode input--------------------------------

svmmodeinput <- ksvm(as.matrix(modetrain[1:9]),as.factor(modetrain[,10]),
                    type="C-svc",kernel="vanilladot",C=0.1,scaled=TRUE)
a <- colSums(svmmodeinput@xmatrix[[1]] * svmmodeinput@coef[[1]])
a0 <- -svmmodeinput@b
pred <- predict(svmmodeinput,modetrain[,1:9])
sum(pred == modetrain[,10]) / nrow(modetrain)
#0.972973
pred <- predict(svmmodeinput,modetest[,1:9])
sum(pred == modetest[,10]) / nrow(modetest)
# 0.9583333

#mean input--------------------------------

svmmeaninput <- ksvm(as.matrix(meantrain[,1:9]),as.factor(meantrain[,10]),
                     type="C-svc",kernel="vanilladot",C=0.1,scaled=TRUE)
a <- colSums(svmmeaninput@xmatrix[[1]] * svmmeaninput@coef[[1]])
a0 <- -svmmeaninput@b
pred <- predict(svmmeaninput,meantrain[,1:9])
sum(pred == meantrain[,10]) / nrow(meantrain)
# 0.972973
pred <- predict(svmmeaninput,meantest[,1:9])
sum(pred == meantest[,10]) / nrow(meantest)
# 0.9583333

#binary input--------------------------------

svmbininput <- ksvm(as.matrix(bintrain[,1:10]),as.factor(bintrain[,11]),
                     type="C-svc",kernel="vanilladot",C=0.1,scaled=TRUE)
a <- colSums(svmbininput@xmatrix[[1]] * svmbininput@coef[[1]])
a0 <- -svmbininput@b
pred <- predict(svmbininput,bintrain[,1:10])
sum(pred == bintrain[,11]) / nrow(bintrain)
# 0.9747748
pred <- predict(svmbininput,bintest[,1:10])
sum(pred == bintest[,11]) / nrow(bintest)
#  0.9583333

#remove missing value rows input--------------------------------
svmnotminput <- ksvm(as.matrix(notmtrain[,1:9]),as.factor(notmtrain[,10]),
                     type="C-svc",kernel="vanilladot",C=0.1,scaled=TRUE)
a <- colSums(svmnotminput@xmatrix[[1]] * svmnotminput@coef[[1]])
a0 <- -svmnotminput@b
pred <- predict(svmnotminput,notmtrain[,1:9])
sum(pred == notmtrain[,10]) / nrow(notmtrain)
# 0.975
pred <- predict(svmnotminput,notmtest[,1:9])
sum(pred == notmtest[,10]) / nrow(notmtest)
# 0.9593496

#perturb input--------------------------------

svmperinput <- ksvm(as.matrix(pertrain[,1:9]),as.factor(pertrain[,10]),
                     type="C-svc",kernel="vanilladot",C=0.1,scaled=TRUE)
a <- colSums(svmperinput@xmatrix[[1]] * svmperinput@coef[[1]])
a0 <- -svmperinput@b
pred <- predict(svmperinput,pertrain[,1:9])
sum(pred == pertrain[,10]) / nrow(pertrain)
#  0.9693694
pred <- predict(svmperinput ,pertest[,1:9])
sum(pred == pertest[,10]) / nrow(pertest)
# 0.9583333

#------------------------------------Logistic Regression------------------------
#regression-------------

regmodel <- lm(Class ~., data = regtrain)
summary(regmodel) # 0.845

reg <- ols_step_forward_aic(regmodel)
reg
#Clump_Thickness, Uniformity_Cell_Size, Uniformity_Cell_Shape, Single_Epithelial_Cell_Size, Bare_Nuclei, Bland_Chromatin, Normal_Nucleoli AIC =  -279.164, R2 = 0.845  
plot(reg)

lmregmod <- glm(formula = Class ~ Clump_Thickness + Uniformity_Cell_Size + Uniformity_Cell_Shape + Bare_Nuclei + Bland_Chromatin + Normal_Nucleoli,
                family="binomial"(link="logit"), data = regtrain)

summary(lmregmod) #Adjusted R-squared: 99.28

yhat <- predict(lmregmod, regtrain, type = "response")
roc(regtrain$Class, round(yhat))
threshold <- 0.5
ythresh <- as.integer(yhat > threshold)
cmatrix <- as.matrix(table(ythresh, regtrain$Class))
cmatrix
(cmatrix[1] + cmatrix[2,2])/(cmatrix[1]+cmatrix[1,2]+cmatrix[2,1] + cmatrix[2,2])
#  0.972973

yhat <- predict(lmregmod, regtest, type = "response")
ythresh <- as.integer(yhat > .5)
cmatrix <- as.matrix(table(ythresh, regtest$Class))
(cmatrix[1] + cmatrix[2,2])/(cmatrix[1]+cmatrix[1,2]+cmatrix[2,1] + cmatrix[2,2])
#   0.9583333

#Mode-------------

lmmodemod <- glm(formula = Class ~ Clump_Thickness + Uniformity_Cell_Size + Uniformity_Cell_Shape + Bare_Nuclei + Bland_Chromatin + Normal_Nucleoli,
                family="binomial"(link="logit"), data = modetrain)

summary(lmmodemod) #AIC  96.079
yhat <- predict(lmmodemod, modetrain, type = "response")
ythresh <- as.integer(yhat > threshold)
cmatrix <- as.matrix(table(ythresh, modetrain$Class))
cmatrix
(cmatrix[1] + cmatrix[2,2])/(cmatrix[1]+cmatrix[1,2]+cmatrix[2,1] + cmatrix[2,2])
#   0.9782609
yhat <- predict(lmmodemod, modetest, type = "response")
ythresh <- as.integer(yhat > .5)
cmatrix <- as.matrix(table(ythresh, modetest$Class))
(cmatrix[1] + cmatrix[2,2])/(cmatrix[1]+cmatrix[1,2]+cmatrix[2,1] + cmatrix[2,2])
#  0.9583333

#Mean-------------
lmmeanmod <- glm(formula = Class ~ Clump_Thickness + Uniformity_Cell_Size + Uniformity_Cell_Shape + Bare_Nuclei + Bland_Chromatin + Normal_Nucleoli, 
                 family="binomial"(link="logit"),
                 data = meantrain)

summary(lmmeanmod) #AIC: 97.835
yhat <- predict(lmmeanmod, meantrain, type = "response")
ythresh <- as.integer(yhat > threshold)
cmatrix <- as.matrix(table(ythresh, meantrain$Class))
cmatrix
(cmatrix[1] + cmatrix[2,2])/(cmatrix[1]+cmatrix[1,2]+cmatrix[2,1] + cmatrix[2,2])
#   0.972973
yhat <- predict(lmmeanmod, meantest, type = "response")
ythresh <- as.integer(yhat > .5)
cmatrix <- as.matrix(table(ythresh, meantest$Class))
cmatrix
(cmatrix[1] + cmatrix[2,2])/(cmatrix[1]+cmatrix[1,2]+cmatrix[2,1] + cmatrix[2,2])
# 0.9583333

#not missing-------------
notmodel <- lm(Class ~., data = notmtrain)
summary(notmodel) # 0.8326 

reg <- ols_step_forward_aic(notmodel)
reg #Clump_Thickness, Uniformity_Cell_Size, Uniformity_Cell_Shape, Marginal_Adhesion, Single_Epithelial_Cell_Size, Bare_Nuclei, Bland_Chromatin, Normal_Nucleoli AIC = -231.868 , R2 = 0.833 
plot(inputfor)

lmnotmod <- glm(formula = Class ~ Clump_Thickness + Uniformity_Cell_Size + Uniformity_Cell_Shape + Marginal_Adhesion + Single_Epithelial_Cell_Size + Bare_Nuclei + Bland_Chromatin + Normal_Nucleoli,
               family="binomial"(link="logit"), data = notmtrain)
summary(lmnotmod) # AIC 116.49

yhat <- predict(lmnotmod, notmtrain, type = "response")
ythresh <- as.integer(yhat > threshold)
cmatrix <- as.matrix(table(ythresh, notmtrain$Class))
cmatrix
(cmatrix[1] + cmatrix[2,2])/(cmatrix[1]+cmatrix[1,2]+cmatrix[2,1] + cmatrix[2,2])
# 0.9714286
yhat <- predict(lmnotmod, notmtest, type = "response")
ythresh <- as.integer(yhat > .5)
cmatrix <- as.matrix(table(ythresh, notmtest$Class))
cmatrix
(cmatrix[1] + cmatrix[2,2])/(cmatrix[1]+cmatrix[1,2]+cmatrix[2,1] + cmatrix[2,2])
# 0.9593496


#binary values------------------
binmodel <- lm(Class ~., data = bintrain)
summary(binmodel) #0.8452

reg <- ols_step_forward_aic(binmodel)
reg #Clump_Thickness, Uniformity_Cell_Size, Uniformity_Cell_Shape, Marginal_Adhesion, Single_Epithelial_Cell_Size, Bare_Nuclei, Bland_Chromatin, Normal_Nucleoli AIC = -231.868 , R2 = 0.833 
#model did not select Bare_Nucleia, we will add it since the goal is to test this
#input, further adjustments needed to refine the model with Bare_Nucleia
lmbinmod <- glm(formula = Class ~ Clump_Thickness + Uniformity_Cell_Size + Uniformity_Cell_Shape + Bare_Nuclei + Bare_Nucleia + Bland_Chromatin,
                family="binomial"(link="logit"), data = bintrain)
summary(lmbinmod) # AIC 95.946

yhat <- predict(lmbinmod, bintrain, type = "response")
ythresh <- as.integer(yhat > threshold)
cmatrix <- as.matrix(table(ythresh, bintrain$Class))
cmatrix
(cmatrix[1] + cmatrix[2,2])/(cmatrix[1]+cmatrix[1,2]+cmatrix[2,1] + cmatrix[2,2])
#     0.9765766
yhat <- predict(lmbinmod, bintest, type = "response")
ythresh <- as.integer(yhat > .5)
cmatrix <- as.matrix(table(ythresh, bintest$Class))
cmatrix
(cmatrix[1] + cmatrix[2,2])/(cmatrix[1]+cmatrix[1,2]+cmatrix[2,1] + cmatrix[2,2])
#  0.9513889

#Perturbation------------------
permodel <- lm(Class ~., data = pertrain)
summary(permodel) #0.8452

reg <- ols_step_forward_aic(permodel)
reg #Clump_Thickness, Uniformity_Cell_Size, Uniformity_Cell_Shape, Marginal_Adhesion, Single_Epithelial_Cell_Size, Bare_Nuclei, Bland_Chromatin, Normal_Nucleoli AIC = -231.868 , R2 = 0.833 
#model did not select Bare_Nucleia, we will add it since the goal is to test this
#input, further adjustments needed to refine the model with Bare_Nucleia
lmpermod <- glm(formula = Class ~ Clump_Thickness + Uniformity_Cell_Size + Single_Epithelial_Cell_Size + Bare_Nuclei + Bland_Chromatin + Normal_Nucleoli,
                family="binomial"(link="logit"), data = bintrain)
summary(lmpermod) # AIC  96.293

yhat <- predict(lmpermod, pertrain, type = "response")
ythresh <- as.integer(yhat > threshold)
cmatrix <- as.matrix(table(ythresh, pertrain$Class))
cmatrix
(cmatrix[1] + cmatrix[2,2])/(cmatrix[1]+cmatrix[1,2]+cmatrix[2,1] + cmatrix[2,2])
#     0.9765766
yhat <- predict(lmpermod, pertest, type = "response")
ythresh <- as.integer(yhat > .5)
cmatrix <- as.matrix(table(ythresh, pertest$Class))
cmatrix
(cmatrix[1] + cmatrix[2,2])/(cmatrix[1]+cmatrix[1,2]+cmatrix[2,1] + cmatrix[2,2])
#  0.9722222

performance_data <- data.frame(
  Method = c("Mean", "Mode", "Regression", "Perturbation", "Binary Indicator", "Remove Rows"),
  Accuracy = c(0.9583333, 0.9583333, 0.9583333, 0.9583333, 0.9583333, 0.9593496) # Replace with your actual results
)

# Generate the plot
library(ggplot2)
ggplot(performance_data, aes(x = reorder(Method, -Accuracy), y = Accuracy, fill = Method)) +
  geom_bar(stat = "identity") +
  coord_cartesian(ylim = c(0.9, 1.0)) + # Zoom in to see the differences
  labs(title = "Impact of Imputation Method on Model Accuracy",
       x = "Imputation Technique",
       y = "Test Set Accuracy") +
  theme_minimal() +
  guides(fill = "none")


performance_data <- data.frame(
  Method = c("Mean", "Mode", "Regression", "Perturbation", "Binary Indicator", "Remove Rows"),
  Accuracy = c(0.972973, 0.9583333, 0.9583333, 0.9722222, 0.9513889, 0.9593496) # Replace with your actual results
)

# Generate the plot
library(ggplot2)
ggplot(performance_data, aes(x = reorder(Method, -Accuracy), y = Accuracy, fill = Method)) +
  geom_bar(stat = "identity") +
  coord_cartesian(ylim = c(0.9, 1.0)) + # Zoom in to see the differences
  labs(title = "Impact of Imputation Method on Model Accuracy",
       x = "Imputation Technique",
       y = "Test Set Accuracy") +
  theme_minimal() +
  guides(fill = "none")

library(ggplot2)
library(tidyr)
library(dplyr)

# 1. Create a "Tidy" dataframe containing both model results
performance_data <- data.frame(
  Method = c("Mean", "Mode", "Regression", "Perturbation", "Binary Indicator", "Remove Rows"),
  Logistic_Regression = c(0.973, 0.958, 0.958, 0.972, 0.951, 0.959),
  SVM = c(0.958, 0.958, 0.958, 0.958, 0.958, 0.959)
)

# 2. Reshape data for ggplot (Long Format)
plot_data <- performance_data %>%
  pivot_longer(cols = c(Logistic_Regression, SVM), 
               names_to = "Model", 
               values_to = "Accuracy")

ggplot(plot_data, aes(x = reorder(Method, -Accuracy), y = Accuracy, fill = Model)) +
  # Use position_dodge to put bars side-by-side
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  # Add text labels on top of bars for immediate clarity
  geom_text(aes(label = round(Accuracy, 3)), 
            position = position_dodge(width = 0.8), 
            vjust = -0.5, size = 3, fontface = "bold") +
  # Focus the Y-axis to highlight the differences
  coord_cartesian(ylim = c(0.94, 0.98)) + 
  # Professional styling
  scale_fill_manual(values = c("#2c3e50", "#e74c3c"), labels = c("Logistic Regression", "SVM")) +
  labs(title = "Model Accuracy by Imputation Strategy",
       subtitle = "Breast Cancer Wisconsin Diagnostic Dataset",
       x = "Imputation Technique",
       y = "Accuracy (Test Set)",
       fill = "Model Type") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    legend.position = "top",
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )



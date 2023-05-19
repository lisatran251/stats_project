# Load necessary libraries
library(dplyr)
library(tidyr)
library(devtools)
library(VariantAnnotation)
library(vcfR)
library(genetics)
library(pROC)
library(Matrix)
library(glmnet)
library(rpart) 
library(randomForest)
library(caret)
library(rpart)
library(pROC)
library(pracma)

# Read in VCF files
asip_afr <- read.vcfR('ASIP_AFR.vcf.gz')
asip_eas <- read.vcfR('ASIP_EAS.vcf.gz')

# Create dataframes for ASIP_AFR and ASIP_EAS
df_asip_afr <- cbind(as.data.frame(getFIX(asip_afr)), INFO2df(asip_afr))
df_asip_eas <- cbind(as.data.frame(getFIX(asip_eas)), INFO2df(asip_eas))

# Filter out rows where AF < 0.001 and QUAL < 20 
df_filtered_asip_afr <- df_asip_afr %>% 
         select_if(~ any(!is.na(.))) %>% 
         filter(AF < 0.001 & QUAL < 20)
df_filtered_asip_eas <- df_asip_eas %>% 
         select_if(~ any(!is.na(.))) %>% 
         filter(AF < 0.001 & QUAL < 20)

## Check minimum DP
# min(df_filtered_asip_afr$DP)
# min(df_filtered_asip_eas$DP)

# Calculate HWE
df_filtered_asip_eas$AF <- as.numeric(df_filtered_asip_eas$AF)
df_filtered_asip_eas$HWE <- c((df_filtered_asip_eas$AF + 2*(df_filtered_asip_eas$AF)*(1-df_filtered_asip_eas$AF)+ (1-df_filtered_asip_eas$AF)))
df_filtered_asip_afr$AF <- as.numeric(df_filtered_asip_afr$AF)
df_filtered_asip_afr$HWE <- c((df_filtered_asip_afr$AF + 2*(df_filtered_asip_afr$AF)*(1-df_filtered_asip_afr$AF)+ (1-df_filtered_asip_afr$AF)))

# Remove all rows with NA values 
df_filtered_asip_eas <- df_filtered_asip_eas[!is.na(df_filtered_asip_eas$HWE),]
df_filtered_asip_afr <- df_filtered_asip_afr[!is.na(df_filtered_asip_afr$HWE),]

# Filter out all the values that have HWE greater than 1.009
df_filtered_asip_eas$HWE[df_filtered_asip_eas$HWE >= 1.009]
df_filtered_asip_afr$HWE[df_filtered_asip_afr$HWE >= 1.009]

# Filter rows where VT is not 'SNP'
final_asip_afr <- df_filtered_asip_afr[df_filtered_asip_afr$VT == 'SNP', ]
final_asip_eas <- df_filtered_asip_eas[df_filtered_asip_eas$VT == 'SNP', ]

# Select required columns from final_asip_afr and final_asip_eas
df_asip <- bind_rows(
         dplyr::select(final_asip_afr, POS, REF, ALT, AFR_AF),
         dplyr::select(final_asip_eas, POS, REF, ALT, EAS_AF)
)
# Create a new column 'SUPERPOP_AF' based on 'afr_AF' and 'eas_AF' values
df_asip <- df_asip %>%
         mutate(SUPERPOP_AF = ifelse(is.na(AFR_AF), EAS_AF, AFR_AF),
                ETH = ifelse(is.na(AFR_AF), 'EAS', 'AFR')) %>%
         dplyr::select(-c(EAS_AF, AFR_AF))

# Convert columns to appropriate types
df_asip$SUPERPOP_AF <- as.numeric(df_asip$SUPERPOP_AF)
df_asip$POS <- as.numeric(df_asip$POS)
df_asip$ETH <- as.factor(df_asip$ETH)
df_asip$REF <- as.factor(df_asip$REF)
df_asip$ALT <- as.factor(df_asip$ALT)

# Attach data frame to search for column names more easily
attach(df_asip)

# Logistic Regression
glm.fit = glm(ETH~SUPERPOP_AF+ POS +ALT + REF, data = df_asip, family='binomial')
summary(glm.fit)
coef(glm.fit)
glm.probs=predict(glm.fit, type='response')
glm.probs[1:10]

contrasts(ETH)

# Predict ETH using the logistic regression model and compare with actual ETH
glm.pred = rep("EAS", 3092)
glm.pred[glm.probs > .5] = "AFR"
table(glm.pred, df_asip$ETH)
(360+1034)/3092
mean(glm.pred==df_asip$ETH)

# ROC Curve
roc.curve = roc(df_asip$ETH, glm.probs, plot = TRUE, print.auc = TRUE, legacy.axes = TRUE, col = "blue",main='Figure 1. ROC Curve of ASIP')

# Residual Plot
glm.pred = ifelse(glm.probs > 0.5, "AFR", "EAS")
plot(resid(glm.fit), predict(glm.fit), pch = 20, col = ifelse(df_asip$ETH == "AFR", "blue", "red"), main = "Figure 3. Residual Plot of ASIP", xlab = "Residuals", ylab = "Predicted Values")
legend("bottomleft", legend = c("AFR", "EAS"), pch = 20, col = c("blue", "red"), title = "Ethnicity")


# Calculate the middle row
mid_row <- nrow(final_asip_afr) %/% 2

# Print the middle row of the column
print(final_asip_afr[mid_row, "POS"])

# Subset data and train a logistic regression model on the subset
train=(POS<34226249)
df_asip.34226249=df_asip[!train,]
dim(df_asip.34226249)
ETH.34226249=ETH[!train]
glm.fit=glm(ETH~SUPERPOP_AF+ POS +ALT + REF, data = df_asip, family='binomial',subset=train)

# Predict ETH using the logistic regression model on the subset and compare with actual ETH
glm.probs=predict(glm.fit, df_asip.34226249, type='response')
#This number is from the dim above
glm.pred=rep("EAS", 1548)
glm.pred[glm.probs>.5]="AFR"
table(glm.pred, ETH.34226249)
mean(glm.pred==ETH.34226249)
mean(glm.pred!=ETH.34226249)

# Linear Discriminant Analysis
lda.fit = lda(ETH ~ SUPERPOP_AF + POS + ALT + REF, subset = train)
lda.pred = predict(lda.fit, newdata = df_asip.34226249)
lda.class = lda.pred$class
table(lda.class, ETH.34226249)

# Create confusion matrix
cm <- confusionMatrix(lda.class, ETH.34226249)

# Extract confusion matrix values
cm_df <- as.data.frame.matrix(cm$table)
cm_df$actual <- rownames(cm_df)
cm_df <- gather(cm_df, predicted, value, -actual)

# Create confusion matrix plot
ggplot(cm_df, aes(x = predicted, y = actual, fill = value)) +
         geom_tile(color = "white") +
         scale_fill_gradient(low = "white", high = "blue") +
         labs(title = "LDA Confusion Matrix",
              x = "Predicted",
              y = "Actual",
              fill = "Count")

mean(lda.class==ETH.34226249)
sum(lda.pred$posterior[,1] >= .5)
sum(lda.pred$posterior[,1] <.5)
lda.pred$posterior[1:20,1]
lda.class[1:20]
sum(lda.pred$posterior[,1]>.9)

# K-Nearest Neighbor
library(class)
train.X=cbind(POS, REF, ALT, SUPERPOP_AF)[train,]
test.X=cbind(POS, REF, ALT, SUPERPOP_AF)[!train,]
train.ETH=ETH[train]
set.seed(10)
# KNN algorithm
knn.pred <- knn(train.X, test.X, train.ETH, k = 1)
# Display a confusion matrix to evaluate the accuracy
table(knn.pred, ETH.34226249)
# Calculate the accuracy
(356 + 380) / 1548

knn.pred <- knn(train.X, test.X, train.ETH, k = 3)
table(knn.pred, ETH.34226249)
(383 + 374) / 1548

# Classification Tree using the tree package
library(tree)
set.seed(345)

trainIndex <- sample(1:nrow(df_asip), 0.7 * nrow(df_asip))
train <- df_asip[trainIndex,]
test <- df_asip[-trainIndex]
# Create the model
model <- tree(ETH ~ SUPERPOP_AF + POS + ALT + REF, data = train)

# Plot the decision tree
plot(model)
text(model)

# Use the model to make predictions on the test data
predictions <- predict(model, test, type = 'class')

# Calculate the accuracy of the predictions
accuracy <- sum(predictions == test$ETH) / nrow(test)
print(paste("Accuracy: ", round(accuracy, 2)))


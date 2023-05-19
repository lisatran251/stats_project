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
library(class)
library(pROC)
library(pracma)
library(plot.matrix)

# Read in VCF files
oca2_afr <- read.vcfR('OCA2_AFR.vcf.gz')
oca2_eas <- read.vcfR('OCA2_EAS.vcf.gz')

# Create dataframes for OCA2_AFR and OCA2_EAS
df_oca2_afr <- cbind(as.data.frame(getFIX(oca2_afr)), INFO2df(oca2_afr))
df_oca2_eas <- cbind(as.data.frame(getFIX(oca2_eas)), INFO2df(oca2_eas))

# Filter out rows where AF < 0.001 and QUAL < 20 
df_filtered_oca2_afr <- df_oca2_afr %>% 
         select_if(~ any(!is.na(.))) %>% 
         filter(AF < 0.001 & QUAL < 20)
df_filtered_oca2_eas <- df_oca2_eas %>% 
         select_if(~ any(!is.na(.))) %>% 
         filter(AF < 0.001 & QUAL < 20)

## Check minimum DP
# min(final_oca2_eas$DP)
# min(final_oca2_afr$DP)

# Remove all rows with double values
df_filtered_oca2_afr <- df_filtered_oca2_afr[!grepl(",", df_filtered_oca2_afr$AF), ]
df_filtered_oca2_eas <- df_filtered_oca2_eas[!grepl(",", df_filtered_oca2_eas$AF), ]

# Calculate HWE
df_filtered_oca2_eas$AF <- as.numeric(df_filtered_oca2_eas$AF)
df_filtered_oca2_eas$HWE <- c((df_filtered_oca2_eas$AF + 2*(df_filtered_oca2_eas$AF)*(1-df_filtered_oca2_eas$AF)+ (1-df_filtered_oca2_eas$AF)))
df_filtered_oca2_afr$AF <- as.numeric(df_filtered_oca2_afr$AF)
df_filtered_oca2_afr$HWE <- c((df_filtered_oca2_afr$AF + 2*(df_filtered_oca2_afr$AF)*(1-df_filtered_oca2_afr$AF)+ (1-df_filtered_oca2_afr$AF)))


# Remove all rows with NA values 
df_filtered_oca2_eas <- df_filtered_oca2_eas[!is.na(df_filtered_oca2_eas$HWE),]
df_filtered_oca2_afr <- df_filtered_oca2_afr[!is.na(df_filtered_oca2_afr$HWE),]

# Filter out all the values that have HWE greater than 1.009
df_filtered_oca2_eas$HWE[df_filtered_oca2_eas$HWE >= 1.009]
df_filtered_oca2_afr$HWE[df_filtered_oca2_afr$HWE >= 1.009]

# Filter rows where VT is not 'SNP'
final_oca2_afr <- df_filtered_oca2_afr[df_filtered_oca2_afr$VT == 'SNP', ]
final_oca2_eas <- df_filtered_oca2_eas[df_filtered_oca2_eas$VT == 'SNP', ]

set.seed(331)
sample_idx <- sample(1:nrow(final_oca2_eas), 1565)
sample_idx2 <- sample(1:nrow(final_oca2_afr), 1565)

final_oca2_eas <- final_oca2_eas[sample_idx, ]
final_oca2_afr <- final_oca2_afr[sample_idx2, ]


# Select required columns from final_oca2_afr and final_oca2_eas
df_oca2 <- bind_rows(
         dplyr::select(final_oca2_afr, POS, REF, ALT, AFR_AF),
         dplyr::select(final_oca2_eas, POS, REF, ALT, EAS_AF)
)
# Create a new column 'SUPERPOP_AF' based on 'afr_AF' and 'eas_AF' values
df_oca2 <- df_oca2 %>%
         mutate(SUPERPOP_AF = ifelse(is.na(AFR_AF), EAS_AF, AFR_AF),
                ETH = ifelse(is.na(AFR_AF), 'EAS', 'AFR')) %>%
         dplyr::select(-c(EAS_AF, AFR_AF))

# Convert columns to appropriate types
df_oca2$SUPERPOP_AF <- as.numeric(df_oca2$SUPERPOP_AF)
df_oca2$POS <- as.numeric(df_oca2$POS)
df_oca2$ETH <- as.factor(df_oca2$ETH)
df_oca2$REF <- as.factor(df_oca2$REF)
df_oca2$ALT <- as.factor(df_oca2$ALT)

attach(df_oca2)

# Logistic Regression
glm.fit = glm(ETH~SUPERPOP_AF+ POS +ALT + REF, data = df_oca2, family='binomial')
summary(glm.fit)
coef(glm.fit)
glm.probs=predict(glm.fit, type='response')
glm.probs[1:10]
contrasts(ETH)

# Predict ETH using the logistic regression model and compare with actual ETH
glm.pred = rep("EAS", 3092)
glm.pred[glm.probs > .5] = "AFR"
table(glm.pred, df_oca2$ETH)
(989+341)/3092

# Fit logistic regression model
glm.fit = glm(ETH ~ SUPERPOP_AF + POS + ALT + REF, data = df_oca2, family = 'binomial')
summary(glm.fit)
coef(glm.fit)
glm.probs = predict(glm.fit, type = 'response')

# ROC Curve
roc.curve = roc(df_oca2$ETH, glm.probs, plot = TRUE, print.auc = TRUE, legacy.axes = TRUE, col = "blue", main='Figure 2. ROC Curve of OCA2')


# Residual Plot
glm.pred = ifelse(glm.probs > 0.5, "AFR", "EAS")
plot(resid(glm.fit), predict(glm.fit), pch = 20, col = ifelse(df_oca2$ETH == "AFR", "blue", "red"), main = "Figure 4. Residual Plot of OCA2", xlab = "Residuals", ylab = "Predicted Values")
legend("bottomleft", legend = c("AFR", "EAS"), pch = 20, col = c("blue", "red"), title = "Ethnicity")



# Calculate the middle row
mid_row <- nrow(final_oca2_afr) %/% 2

# Print the middle row of the column
print(final_oca2_afr[mid_row, "POS"])

# Subset data and train a logistic regression model on the subset
train=(POS<27947360)
df_oca2.27947360=df_oca2[!train,]
dim(df_oca2.27947360)
ETH.27947360=ETH[!train]
glm.fit=glm(ETH~SUPERPOP_AF+ POS +ALT + REF, data = df_oca2, family='binomial',subset=train)

# Predict ETH using the logistic regression model on the subset and compare with actual ETH
glm.probs=predict(glm.fit, df_oca2.27947360, type='response')
#This number is from the dim above
glm.pred=rep("EAS", 1236)
glm.pred[glm.probs>.5]="AFR"
table(glm.pred, ETH.27947360)
(446+113)/1236
mean(glm.pred==ETH.27947360)
mean(glm.pred!=ETH.27947360)


# Linear Discriminant Analysis
lda.fit = lda(ETH ~ SUPERPOP_AF + POS + ALT + REF, subset = train)
lda.pred = predict(lda.fit, newdata = df_oca2.27947360)
lda.class = lda.pred$class
table(lda.class, ETH.27947360)

# Create confusion matrix
cm <- confusionMatrix(lda.class, ETH.27947360)

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

mean(lda.class==ETH.27947360)
sum(lda.pred$posterior[,1] >= .5)
sum(lda.pred$posterior[,1] <.5)
lda.pred$posterior[1:20,1]
lda.class[1:20]
sum(lda.pred$posterior[,1]>.9)

# K-Nearest Neighbor
train.X=cbind(POS, REF, ALT, SUPERPOP_AF)[train,]
test.X=cbind(POS, REF, ALT, SUPERPOP_AF)[!train,]
train.ETH=ETH[train]
set.seed(10)
# KNN algorithm
knn.pred <- knn(train.X, test.X, train.ETH, k = 1)
# Display a confusion matrix to evaluate the accuracy
table(knn.pred, ETH.27947360)
# Calculate the accuracy
(604 + 0) / 1236

knn.pred <- knn(train.X, test.X, train.ETH, k = 3)
table(knn.pred, ETH.27947360)
(632 + 0) / 1236


# Classification Tree using the tree library
library(tree)
set.seed(345)

trainIndex <- sample(1:nrow(df_oca2), 0.7 * nrow(df_oca2))
train <- df_oca2[trainIndex,]
test <- df_oca2[-trainIndex]
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


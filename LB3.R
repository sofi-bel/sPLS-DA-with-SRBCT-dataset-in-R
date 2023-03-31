library(mixOmics) # import the mixOmics library

set.seed(5249) # for reproducibility, remove for normal use
data(srbct) # extract the small round bull cell tumour data
X <- srbct$gene # use the gene expression data as the X matrix
Y <- srbct$class # use the class data as the Y matrix

dim(X) # check the dimensions of the X dataframe
summary(Y) # check the distribution of class labels

# run pca method on data
pca.srbct = pca(X, ncomp = 10, center = TRUE, scale = TRUE) 
plot(pca.srbct)  # barplot of the eigenvalues (explained variance per component)

plotIndiv(pca.srbct, group = srbct$class, ind.names = FALSE, # plot the samples projected
          legend = TRUE, title = 'PCA on SRBCT, comp 1 - 2') # onto the PCA subspace
srbct.splsda <- splsda(X, Y, ncomp = 10)  # set ncomp to 10 for performance assessment later

# plot the samples projected onto the first two components of the PLS-DA subspace
plotIndiv(srbct.splsda , comp = 1:2, 
          group = srbct$class, ind.names = FALSE,  # colour points by class
          ellipse = TRUE, # include 95% confidence ellipse for each class
          legend = TRUE, title = '(a) PLSDA with confidence ellipses')

# use the max.dist measure to form decision boundaries between classes based on PLS-DA data
background = background.predict(srbct.splsda, comp.predicted=2, dist = "max.dist")

# plot the samples projected onto the first two components of the PLS-DA subspace
plotIndiv(srbct.splsda, comp = 1:2,
          group = srbct$class, ind.names = FALSE, # colour points by class
          background = background, # include prediction background for each class
          legend = TRUE, title = " (b) PLSDA with prediction background")

# undergo performance evaluation in order to tune the number of components to use
perf.splsda.srbct <- perf(srbct.splsda, validation = "Mfold", 
                          folds = 5, nrepeat = 10, # use repeated cross-validation
                          progressBar = FALSE, auc = TRUE) # include AUC values

# plot the outcome of performance evaluation across all ten components
plot(perf.splsda.srbct, col = color.mixo(5:7), sd = TRUE,
     legend.position = "horizontal")

perf.splsda.srbct$choice.ncomp # what is the optimal value of components according to perf()

# grid of possible keepX values that will be tested for each component
list.keepX <- c(1:10,  seq(20, 300, 10))

# undergo the tuning process to determine the optimal number of variables
tune.splsda.srbct <- tune.splsda(X, Y, ncomp = 4, # calculate for first 4 components
                                 validation = 'Mfold',
                                 folds = 5, nrepeat = 10, # use repeated cross-validation
                                 dist = 'max.dist', # use max.dist measure
                                 measure = "BER", # use balanced error rate of dist measure
                                 test.keepX = list.keepX,
                                 cpus = 2) # allow for paralleliation to decrease runtime

plot(tune.splsda.srbct, col = color.jet(4)) # plot output of variable number tuning

tune.splsda.srbct$choice.ncomp$ncomp # what is the optimal value of components according to tune.splsda()
tune.splsda.srbct$choice.keepX # what are the optimal values of variables according to tune.splsda()
optimal.ncomp <- tune.splsda.srbct$choice.ncomp$ncomp
optimal.keepX <- tune.splsda.srbct$choice.keepX[1:optimal.ncomp]
# form final model with optimised values for component and variable count
final.splsda <- splsda(X, Y, 
                       ncomp = optimal.ncomp, 
                       keepX = optimal.keepX)
plotIndiv(final.splsda, comp = c(1,2), # plot samples from final model
          group = srbct$class, ind.names = FALSE, # colour by class label
          ellipse = TRUE, legend = TRUE, # include 95% confidence ellipse
          title = ' (a) sPLS-DA on SRBCT, comp 1 & 2')

plotIndiv(final.splsda, comp = c(1,3), # plot samples from final model
          group = srbct$class, ind.names = FALSE,  # colour by class label
          ellipse = TRUE, legend = TRUE, # include 95% confidence ellipse
          title = '(b) sPLS-DA on SRBCT, comp 1 & 3')
# set the styling of the legend to be homogeneous with previous plots
legend=list(legend = levels(Y), # set of classes
            col = unique(color.mixo(Y)), # set of colours
            title = "Tumour Type", # legend title
            cex = 0.7) # legend size

# generate the CIM, using the legend and colouring rows by each sample's class
cim <- cim(final.splsda, row.sideColors = color.mixo(Y), 
           legend = legend)
# form new perf() object which utilises the final model
perf.splsda.srbct <- perf(final.splsda, 
                          folds = 5, nrepeat = 10, # use repeated cross-validation
                          validation = "Mfold", dist = "max.dist",  # use max.dist measure
                          progressBar = FALSE)

# plot the stability of each feature for the first three components, 'h' type refers to histogram
par(mfrow=c(1,3))
plot(perf.splsda.srbct$features$stable[[1]], type = 'h', 
     ylab = 'Stability', 
     xlab = 'Features', 
     main = '(a) Comp 1', las =2)
plot(perf.splsda.srbct$features$stable[[2]], type = 'h', 
     ylab = 'Stability', 
     xlab = 'Features', 
     main = '(b) Comp 2', las =2)
plot(perf.splsda.srbct$features$stable[[3]], type = 'h', 
     ylab = 'Stability', 
     xlab = 'Features',
     main = '(c) Comp 3', las =2)
var.name.short <- substr(srbct$gene.name[, 2], 1, 10) # form simplified gene names

plotVar(final.splsda, comp = c(1,2), var.names = list(var.name.short), cex = 3) # generate correlation circle plot

train <- sample(1:nrow(X), 50) # randomly select 50 samples in training
test <- setdiff(1:nrow(X), train) # rest is part of the test set

# store matrices into training and test set:
X.train <- X[train, ]
X.test <- X[test,]
Y.train <- Y[train]
Y.test <- Y[test]

# train the model
train.splsda.srbct <- splsda(X.train, Y.train, ncomp = optimal.ncomp, keepX = optimal.keepX)

# use the model on the Xtest set
predict.splsda.srbct <- predict(train.splsda.srbct, X.test, 
                                dist = "mahalanobis.dist")
# evaluate the prediction accuracy for the first two components
predict.comp2 <- predict.splsda.srbct$class$mahalanobis.dist[,2]
table(factor(predict.comp2, levels = levels(Y)), Y.test)

# evaluate the prediction accuracy for the first three components
predict.comp3 <- predict.splsda.srbct$class$mahalanobis.dist[,3]
table(factor(predict.comp3, levels = levels(Y)), Y.test)

auc.splsda = auroc(final.splsda, roc.comp = 1, print = FALSE) # AUROC for the first component

auc.splsda = auroc(final.splsda, roc.comp = 3, print = FALSE) # AUROC for all three components

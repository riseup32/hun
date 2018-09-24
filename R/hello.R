# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'


##########################################################################
edaplot <- function(x,y){
  par(mfrow=c(2,2))
  if(is.numeric(x)) {hist(x,prob=TRUE)}
  if(is.numeric(x)) {lines(sort(x),dnorm(sort(x),mean(x),sd(x)),col='red')}
  if(is.character(x) | is.factor(x)) {plot(x)}

  hist(y,prob=TRUE)
  lines(sort(y),dnorm(sort(y),mean(y),sd(y)),col='red')

  if(is.numeric(x)) {plot(x,y,main='linear regression')}
  if(is.numeric(x)) {abline(lm(y~x),col='red')}
  if(is.character(x) | is.factor(x)) {plot(y~x)}

  if(is.numeric(x)) {plot(x,y,main='lowess')}
  if(is.numeric(x)) {lines(lowess(x,y),col='red')}
  gmm <- 0
  for(i in 1:length(levels(x))){
    gmm <- gmm+sum(x==levels(x)[i])/length(x)*dnorm(min(y):max(y),mean(y[x==levels(x)[i]]),sd(y[x==levels(x)[i]]))
  }
  if(is.character(x) | is.factor(x)) {hist(y,breaks=25,prob=T,main='GMM')}
  if(is.character(x) | is.factor(x)) {lines(min(y):max(y),gmm,col='red')}
}
##########################################################################
rmse <- function(y, yhat){
  sqrt(mean((y - yhat)^2))
}
##########################################################################
k.fold.cv <- function(model,data,y,k,mtry=round(sqrt(ncol(data))),ntree=500,kernel='radial',cost=1,scale=TRUE){
  n <- round((nrow(data)/k)-0.5)
  sum.mse <- 0
  for(i in 1:k){
    m <- (n*(i-1)+1):(n*i)
    x1 <- data[,-grep(y,colnames(data))]
    y1 <- data[,grep(y,colnames(data))]
    cv.test_x <- x1[m,]
    cv.train_x <- x1[-m,]
    cv.test_y <- y1[m]
    cv.train_y <- y1[-m]
    for(j in model){
      if(j=='lr') {cv.model <- glm(cv.train_y~.,data=cv.train_x,family='binomial')}
      else if(j=='mn') {cv.model <- multinom(cv.train_y~.,data=cv.train_x)}
      else if(j=='dt') {cv.model <- ctree(cv.train_y~.,data=cv.train_x)}
      else if(j=='rf') {cv.model <- randomForest(cv.train_y~.,data=cv.train_x,mtry=mtry,ntree=ntree)}
      else if(j=='nb') {cv.model <- naiveBayes(cv.train_y~.,data=cv.train_x)}
      else if(j=='lm') {cv.model <- lm(cv.train_y~.,data=cv.train_x)}
      else if(j=='svm') {cv.model <- svm(cv.train_y~.,data=cv.train_x,kernel=kernel,cost=cost,scale=scale)}
    }
    predict.value <- predict(cv.model,cv.test_x)

    if(model=='lr') {predict.value <- predict.value <- round(1/(1+exp(-predict.value)))}
    if(model=='lr') {predict.value <- as.factor(predict.value)}

    mse <- mean((as.numeric(predict.value)-as.numeric(cv.test_y))^2)
    sum.mse <- sum.mse+mse
  }
  return(sum.mse)
}
##########################################################################
machine <- function(model,data,test,y,outlier=TRUE,method='normal',mtry=round(sqrt(ncol(data))),ntree=500,kernel='radial',cost=1,scale=TRUE){
  if(length(levels(data[,grep(y,colnames(data))])) > length(levels(test[,grep(y,colnames(test))]))) {levels(test[,grep(y,colnames(test))]) <- levels(data[,grep(y,colnames(data))])}
  if(length(levels(data[,grep(y,colnames(data))])) < length(levels(test[,grep(y,colnames(test))]))) {levels(data[,grep(y,colnames(data))]) <- levels(test[,grep(y,colnames(test))])}

  for(i in 1:ncol(data)){
    if(is.character(data[,i])) {data[,i] <- as.factor(data[,i])}
  }

  for(i in 1:ncol(test)){
    if(is.character(test[,i])) {test[,i] <- as.factor(test[,i])}
  }

  A <- NULL
  for(i in 1:ncol(data)){
    if(is.factor(data[,i]) & length(levels(data[,i]))>15) {A <- c(A,FALSE)}
    else {A <- c(A,TRUE)}
  }

  data <- data[,A]

  for(i in 1:ncol(data)){
    if(sum(is.na(data[,i]))>nrow(data)/4) {data[,i] <- NULL}
  }

  test <- test[,colnames(data)]

  A <- NULL
  for(i in 1:ncol(data)){
    A <- c(A,is.numeric(data[,i]))
  }
  B <- cor(data[apply(is.na(data),1,sum)==0,A])

  for(i in 1:ncol(data)){
    for(j in 1:nrow(data)){
      if(is.numeric(data[,i])) {C <- colnames(B)[abs(B[,grep(colnames(data)[i],colnames(B))])==sort(abs(B[,grep(colnames(data)[i],colnames(B))]),decreasing=T)[2]]}
      if(is.numeric(data[,i])) {lm <- lm(data[,i][!is.na(data[,grep(C,colnames(data))])]~data[,grep(C,colnames(data))][!is.na(data[,grep(C,colnames(data))])],data=data)}
      if(is.numeric(data[j,i]) & is.na(data[j,i])) {data[j,i] <- predict(lm,list=data[,grep(C,colnames(data))])[j]}
      else if(is.factor(data[j,i]) & is.na(data[j,i])) {data[j,i] <- levels(data[,i])[which(table(data[,i])==max(table(data[,i])))]}
    }
  }

  for(i in 1:ncol(test)){
    for(j in 1:nrow(test)){
      if(is.numeric(test[,i])) {C <- colnames(B)[abs(B[,grep(colnames(data)[i],colnames(B))])==sort(abs(B[,grep(colnames(data)[i],colnames(B))]),decreasing=T)[2]]}
      if(is.numeric(test[,i])) {lm <- lm(data[,i][!is.na(data[,grep(C,colnames(data))])]~data[,grep(C,colnames(data))][!is.na(data[,grep(C,colnames(data))])],data=data)}
      if(is.numeric(test[j,i]) & is.na(test[j,i])) {test[j,i] <- predict(lm,list=test[,grep(C,colnames(test))])[j]}
      else if(is.factor(test[j,i]) & is.na(test[j,i])) {test[j,i] <- levels(test[,i])[which(table(test[,i])==max(table(test[,i])))]}
    }
  }

  for(i in 1:ncol(data)){
    if(outlier==T & is.numeric(data[,i])) {data <- data[data[,i]>=quantile(data[,i],0.005) & data[,i]<=quantile(data[,i],0.995),]}
  }

  for(i in 1:ncol(data)){
    if(is.factor(data[,i]) & (length(levels(data[,i]))>length(levels(test[,i])))) {levels(test[,i]) <- levels(data[,i])}
    if(is.factor(data[,i]) & (length(levels(data[,i]))<length(levels(test[,i])))) {levels(data[,i]) <- levels(test[,i])}
  }

  num <- NULL
  for(i in 1:ncol(data)){
    if(is.numeric(data[,i])) {num <- c(num,TRUE)}
    else {num <- c(num,FALSE)}
  }
  if(scale==TRUE) {data[,num] <- scale(data[,num])}

  n_pca <- which(summary(prcomp(data[,num]))$importance[3,]>=0.85)[1]
  if(method=='pca') {data <- data.frame(prcomp(data[,num])$x[,1:n_pca],data[,!num])}

  num2 <- NULL
  for(i in 1:ncol(test)){
    if(is.numeric(test[,i])) {num2 <- c(num2,TRUE)}
    else {num2 <- c(num2,FALSE)}
  }
  if(scale==TRUE) {test[,num2] <- scale(test[,num2])}

  if(method=='pca') {test <- data.frame(prcomp(test[,num2])$x[,1:n_pca],test[,!num2])}

  x1 <- data[,-grep(y,colnames(data))]
  y1 <- data[,grep(y,colnames(data))]
  y2 <- test[,grep(y,colnames(test))]

  k.fold.cv <- function(model,k,data){
    n <- round((nrow(data)/k)-0.5)
    sum.mse <- 0
    for(i in 1:k){
      m <- (n*(i-1)+1):(n*i)
      cv.test_x <- x1[m,]
      cv.train_x <- x1[-m,]
      cv.test_y <- y1[m]
      cv.train_y <- y1[-m]
      for(j in model){
        if(j=='lr') {cv.model <- glm(cv.train_y~.,data=cv.train_x,family='binomial')}
        else if(j=='mn') {cv.model <- multinom(cv.train_y~.,data=cv.train_x)}
        else if(j=='dt') {cv.model <- ctree(cv.train_y~.,data=cv.train_x)}
        else if(j=='rf') {cv.model <- randomForest(cv.train_y~.,data=cv.train_x,mtry=mtry,ntree=ntree)}
        else if(j=='nb') {cv.model <- naiveBayes(cv.train_y~.,data=cv.train_x)}
        else if(j=='lm') {cv.model <- lm(cv.train_y~.,data=cv.train_x)}
        else if(j=='svm') {cv.model <- svm(cv.train_y~.,data=cv.train_x,kernel=kernel,cost=cost,scale=scale)}
      }
      predict.value <- predict(cv.model,cv.test_x)

      if(model=='lr') {predict.value <- predict.value <- round(1/(1+exp(-predict.value)))}
      if(model=='lr') {predict.value <- as.factor(predict.value)}

      mse <- mean((as.numeric(predict.value)-as.numeric(cv.test_y))^2)
      sum.mse <- sum.mse+mse
    }
    return(sum.mse)
  }

  cv.value <- k.fold.cv(model,10,data)

  model.Nm <- model
  for(i in model.Nm){
    if(i=='lr') {modelNm <- 'logistic regression'}
    else if(i=='mn') {modelNm <- 'multinomial logistic regression'}
    else if(i=='dt') {modelNm <- 'decision tree'}
    else if(i=='rf') {modelNm <- 'random forest'}
    else if(i=='nb') {modelNm <- 'naive bayes'}
    else if(i=='lm') {modelNm <- 'linear regression'}
    else if(i=='svm') {modelNm <- 'support vector machine'}
  }

  for(i in model){
    if(i=='lr') {model <- glm(y1~.,data=x1,family='binomial')}
    else if(i=='mn') {model <- multinom(y1~.,data=x1)}
    else if(i=='dt') {model <- ctree(y1~.,data=x1)}
    else if(i=='rf') {model <- randomForest(y1~.,data=x1,mtry=mtry,ntree=ntree)}
    else if(i=='nb') {model <- naiveBayes(y1~.,data=x1)}
    else if(i=='lm') {model <- lm(y1~.,data=x1)}
    else if(i=='svm') {model <- svm(y1~.,data=x1,kernel=kernel,cost=cost,scale=scale)}
  }

  predict <- predict(model,test)

  if(model.Nm=='lr') {predict <- round(1/(1+exp(-predict)))}
  if(model.Nm=='lr') {predict <- as.factor(predict)}

  value <- mean(as.numeric(predict)==as.numeric(y2))
  prediction <- paste0(round(value*100,2),'%')
  rmse <- rmse(predict,y2) 

  y <- y2
  if(is.factor(y)) {cat('model:',modelNm,'\n',
      'cv.value',':',cv.value,'\n',
      'prediction',':',prediction,'\n')}
  if(is.factor(y)) {cat('\n','Confusion matrix','\n')}
  if(is.factor(y)) {print(table(y,predict))}
  if(is.factor(y)) {for(i in 1:length(levels(y))){
    cat(levels(y)[i],':',table(y,predict)[i,i]/sum(table(y,predict)[i,]),'\n')
  }}

  if(is.numeric(y)) {cat('model:',modelNm,'\n',
      'cv.value',':',cv.value,'\n',
      'rmse',':',rmse,'\n')}
}
##########################################################################
as.fourier <- function(data,x,period,k){
 A <- NULL
 for(i in 1:k){
  Sin <- sin(i*(2*pi/period)*data[,grep(x,colnames(data))])
  Cos <- cos(i*(2*pi/period)*data[,grep(x,colnames(data))])
  A <- cbind(A,Sin,Cos)
 }
 
 df <- as.data.frame(A)
 colName <- colnames(data)[-grep(x,colnames(data))]
 
 for(i in 1:k){
  colnames(df)[2*i-1] <- paste0(x,'_','Sin',i)
  colnames(df)[2*i] <- paste0(x,'_','Cos',i)
 }
 data <- data.frame(data[,-grep(x,colnames(data))],df)
 
 colnames(data)[1:length(colName)] <- colName
 data <- data
}
##########################################################################
as.Bspline <- function(data,x,kernel,interval){
 round.kernel <- function(val){
  return(sqrt((interval+val)*(interval-val)))
  }

 angle.kernel <- function(val){
  if(val<=interval) {value <- val+interval}
  else {value <- val-interval}
  return(value)
 }

 for(i in kernel){
  if(i=='round') {kernel <- round.kernel}
  else if(i=='angle') {kernel <- angle.kernel}
 }
 
 B <- NULL
 for(i in 0:round(nrow(data)/interval)){
  A <- vector(length=nrow(data))
  for(j in 1:nrow(data)){
   if((data[,grep(x,colnames(data))][j] >= interval*(i-1)) & (data[,grep(x,colnames(data))][j] <= interval*(i+1))) {A[j] <- kernel(j-interval*i)}
   else {A[j] <- 0}
  }
  B <- cbind(B,A)
 }

 df <- as.data.frame(B)
 colName <- colnames(data)[-grep(x,colnames(data))]

 for(i in 1:(round(nrow(data)/interval)+1)){
  colnames(df)[i] <- paste0(x,'_B',i)
 }
 data <- data.frame(data[,-grep(x,colnames(data))],df)

 colnames(data)[1:length(colName)] <- colName
 data <- data
}
##########################################################################
pipeline <- function(data,y){
 for(i in 1:ncol(data)){
  if(is.numeric(data[,i]) & !(sum(is.na(data[,i]))==0)) {data[,i][is.na(data[,i])] <- median(data[,i],na.rm=T)}
  if(is.numeric(data[,i]) & colnames(data)[i]!=colnames(wine.data)[grep(y,colnames(data))]) {data[,i] <- scale(data[,i])}
  if(is.character(data[,i])) {data[,i] <- as.factor(data[,i])}
  if(is.factor(data[,i]) & !(sum(is.na(data[,i]))==0)) {data[,i][is.na(data[,i])] <- levels(data[,i])[which(table(data[,i])==max(table(data[,i])))]}
 }
 data <- data
}
##########################################################################
precision_score <- function(y,predict){
 T <- table(y,predict)
 return(T[2,2]/(T[1,2]+T[2,2]))
}
##########################################################################
recall_score <- function(y,predict){
 T <- table(y,predict)
 return(T[2,2]/(T[2,1]+T[2,2]))
}
##########################################################################
f1_score <- function(y,predict){
 T <- table(y,predict)
 precision_score <- T[2,2]/(T[1,2]+T[2,2])
 recall_score <- T[2,2]/(T[2,1]+T[2,2])
 return(2/((1/precision_score)+(1/recall_score)))
}
##########################################################################
confusion_matrix <- function(y,predict){
 value <- mean(as.numeric(predict)==as.numeric(y))
 prediction <- paste0(round(value*100,2),'%')
 cat('prediction',':',prediction,'\n')
 print(table(y,predict))
 for(i in 1:length(levels(y))){
  cat(levels(y)[i],':',table(y,predict)[i,i]/sum(table(y,predict)[i,]),'\n')
 }
}
##########################################################################
hunhelp <- function(func){
 for(i in func){
  if(i=='machine') {cat('machine(model,data,test,y,outlier=TRUE,method=normal,mtry=round(sqrt(ncol(data))),ntree=500,kernel=radial,cost=1,scale=TRUE)','\n')}
  else if(i=='edaplot') {cat('edaplot(x,y)','\n')}
  else if(i=='rmse') {cat('rmse(y,yhat)','\n')}
  else if(i=='k.fold.cv') {cat('k.fold.cv(model,data,y,k,mtry=round(sqrt(ncol(data))),ntree=500,kernel=radial,cost=1,scale=TRUE)','\n')}
  else if(i=='as.fourier') {cat('as.fourier(data,x,period,k)','\n')}
  else if(i=='as.Bspline') {cat('as.Bspline(data,x,kernel,interval)','\n')}
 }
}

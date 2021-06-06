# loading the package
library(leaps)
library(knitr)
library(mice)
library(MASS)

# Load the data into R studio
xe <- read.table("c:\\Users\\16318\\Desktop\\AMS 578\\578_E.csv", head = T, sep = ",")
xg <- read.table("c:\\Users\\16318\\Desktop\\AMS 578\\578_G.csv", head = T, sep = ",")
y <- read.table("c:\\Users\\16318\\Desktop\\AMS 578\\578_Y.csv", head = T, sep = ",")

# Merge the independ variables into R studio
xeg <- merge(xe,xg,by="ID", all=TRUE)
dat <- merge(y,xeg,by="ID", all=TRUE)
# Remove the ID column
dat <- dat[,-1]

# checking coefficents
r=c()
for(i in 2:32){
  p=cor.test(dat[,1],dat[,i])
  r=c(r,p$p.value)
}

# the number of NA missing value is 174 in columns Y, total entries is 380 and total row is 455
sum(is.na(dat$Y))
sum(complete.cases(dat))
table(is.na(dat))

# Mean and SD in all DV/IV
sapply(xeg,mean,na.rm=TRUE)
sapply(xe,sd,na.rm=TRUE)
summary(dat)


# Imputation method1:dat1 mean value add into NA ----------------------------------------------------------
dat1 <- dat
dim1 <- sum(!complete.cases(dat))
View(dim1)
dim2 <- dim(dat)[2]
View(dim2)
temp <- rep(1, dim1)%*%t(apply(dat, 2, mean, na.rm=TRUE))
dim(temp)
dat1 <- as.data.frame(dat1)
dat1[!complete.cases(dat), ] <- data.frame(rep(1, dim1)%*%t(apply(dat, 2, mean, na.rm=TRUE)))
for(j in 8:dim(dat1)[2]){
  dat1[dat1[,j]>=0.5, j] <- 1
  dat1[dat1[,j]<0.5, j] <- 0
}
View(dat1)

#-----------------------------------------------------------

# Imputation method 2:dat2 Make up the missing value in the data set using mice

dat2 <- mice(dat,m=3, maxit=50, meth='pmm',seed=500)
dat2 <- complete(dat2, 1)
#write.csv(dat, "c:\\Users\\16318\\Desktop\\AMS 578\\dat2.csv")
#dat2 <- read.csv("c:\\Users\\16318\\Desktop\\AMS 578\\dat2.csv", row.names = NULL)

#---------------------------------------------------------------
# Imputation method 3: dat_a Make up the missing value in the data set using average of its neighbors
dat_a <- dat
for( i in 1:dim(dat_a)[1]){
  if(is.na(dat_a[i,1])){
    dat_a[i,1] <- (dat_a[i-1,1]+dat_a[i+1,1]/2)
  }
}
dim1 <- sum(!complete.cases(dat_a))
temp <- rep(1, dim1)%*%t(apply(dat_a, 2, mean, na.rm=TRUE))
dat_a <- as.data.frame(dat_a)
dat_a[!complete.cases(dat_a), ] <- data.frame(rep(1, dim1)%*%t(apply(dat_a, 2, mean, na.rm=TRUE)))
for(j in 8:dim(dat_a)[2]){
  dat_a[dat_a[,j]>=0.5, j] <- 1
  dat_a[dat_a[,j]<0.5, j] <- 0
}
View(dat_a)
#------------------------------------------------------------------
# Imputation method 4: dat_omit
dat_omit <- dat[complete.cases(dat),]
View(dat_omit)
#-----------------------------------------------------------------

# Box_Cox power trsformation plot profile loglikehood
M <- lm(Y ~(.),data=dat2)
boxcox(M, lambda = seq(-3, 3, length.out = 200))
plot(M)
QQplot(M)
plot(resid(M~fitted(M), main="NEw Residual"))
summary(M)

M1 <- lm(Y ~(.)^2,data=dat2)
summary(M1)

# m2 <- lm(Y ~(.)^3,data=dat)
# summary(m2)

# check significant E variables
E=paste0("E",1:6,collapse="+")
M_E <- lm(paste0("Y~",E),data=dat2)
summary(M_E)

# check significant G variables
G=paste0("G",1:25,collapse="+")
M_G <- lm(paste0("Y~", G), data=dat1)
summary(M_G)

#-------------------------------------------------------------------
# Selecting model
M_raw1<- lm(Y~(E1+E4+E6+G9+G5+G22)^3,data=dat2)
summary(M_raw1) #E1+E4+E6+E6:G9 R^2=0.4483 F=20.8 P=2.2e-16 T>=4.4  BIC=6326.38
M_raw2<- lm(Y~(E1+E4+E6+G9+G5+G22)^2,data=dat2)
summary(M_raw2) #E1+E4+E6:G9  R^2=0.4504 F=39.99 P=2.2e-16 T>=6.7  BIC=6204.97
M_raw3<- lm(Y~(E1+E4+E6+G9+G5+G22),data=dat2)
summary(M_raw3) #E1+E4+E6   R^2=0.4041 F=113.9 P=2.2e-16  T>=10.3  BIC=6197.49
BIC(M_raw3)
#-------------------------------------------------------------------

M_raw <- lm(Y ~ (.)^3, data=dat2)
plot(M_raw)
qqnorm(resid(M_raw))
qqline(resid(M_raw))

M_raw_2 <- lm(Y ~ (.)^2, data=dat2)
plot(M_raw_2)
qqnorm(resid(M_raw_2))
qqline(resid(M_raw_2))

# Stepwise forward
# regsubsets function to perform variable selection
M=regsubsets(model.matrix(M_raw_2)[,-1], as.matrix(dat2$Y[complete.cases(dat2)]) ,data=dat2, method='forward', intercept=TRUE,
             nbest=1, nvmax=5)
temp <- summary(M)
# extract selected variables-------------------------------------------
Var <- colnames(model.matrix(M_raw_2))
M_select <- apply(temp$which,1,function(x)paste0(Var[x], collapse="+"))
kable(data.frame(cbind(model=M_select, adjR2=temp$adjr2, BIC=temp$bic, CP=temp$cp)))
# ----------------------------------------------------------------------
temp=summary(M)
summary(M)
kable(data.frame(cbind(model=M, adjR2=temp$adjr2, BIC=temp$bic)), caption='Model Summary')
M_main <- lm(Y~(.), data=dat2)
temp <- summary(M_main)
kable(temp$coefficients[abs(temp$coefficients[,4]) <= 0.001,], caption = 'Sig Coefficients')
M_2nd <- lm(Y~(.)^2,data=dat2)
temp <- summary(M_2nd)
kable(temp$coefficients[abs(temp$coefficients[,4]) <= 0.001,], caption = '2nd Interaction')
M_2stage <- lm(Y~(.)^3,data=dat2)
temp <- summary(M_2stage)
temp$coefficients[abs(temp$coefficients[,3])>=4,]

test_1 <- lm(Y~(E1+E4+E6+E6:G9), data=dat2) 
summary(test_1)

test_2 <- lm(Y~(E1+E4+E6:G9), data=dat2) #dat2 R^2=0.4422 F=265 t>=14//dat1 R^2=0.4857 F=315.4 t>=15//
summary(test_2)   #dat_a R^2=0.05346 F=19.81 t>=0.83//dat_omit R^2=0.4847 F=195.1 t>=11.8
pairwise.t.test(dat2$E1, dat2$E4,dat2$E6,dat2$G9,p.adjust.method = "bonferroni")

plot(test_2)
QQlot(test_2)

aov(test_2)
anova(test_2)

test_3 <- lm(Y~(E1+E4+E6), data=dat2)
summary(test_3)


# data_set: dat1

M_raw_1 <- lm(Y ~ (.)^2, data=dat1)
M=regsubsets(model.matrix(M_raw_1)[,-1], as.matrix(dat1$Y[complete.cases(dat1)]) ,data=dat1, method='forward', intercept=TRUE,
             nbest=1, nvmax=5)
temp <- summary(M)
# extract selected variables-------------------------------------------
Var <- colnames(model.matrix(M_raw_1))
M_select <- apply(temp$which,1,function(x)paste0(Var[x], collapse="+"))
kable(data.frame(cbind(model=M_select, adjR2=temp$adjr2, BIC=temp$bic, CP=temp$cp)))
# --------------------------------------------------------


# data_set: dat_a
M_raw_a <- lm(Y ~ (.)^2, data=dat_a)
M=regsubsets(model.matrix(M_raw_a)[,-1], as.matrix(dat1$Y[complete.cases(dat_a)]) ,data=dat_a, method='forward', intercept=TRUE,
             nbest=1, nvmax=5)
temp <- summary(M)
# extract selected variables-------------------------------------------
Var <- colnames(model.matrix(M_raw_2))
M_select <- apply(temp$which,1,function(x)paste0(Var[x], collapse="+"))
kable(data.frame(cbind(model=M_select, adjR2=temp$adjr2, BIC=temp$bic, CP=temp$cp)))
# ----------------------------------------------------------------------

# data_set: dat_omit
M_raw_omit <- lm(Y ~ (.)^2, data=dat_omit)
M=regsubsets(model.matrix(M_raw_omit)[,-1], as.matrix(dat_omit$Y[complete.cases(dat_omit)]) ,data=dat_omit, method='forward', intercept=TRUE,
             nbest=1, nvmax=5)
temp <- summary(M)
# extract selected variables-------------------------------------------
Var <- colnames(model.matrix(M_raw_omit))
M_select <- apply(temp$which,1,function(x)paste0(Var[x], collapse="+"))
kable(data.frame(cbind(model=M_select, adjR2=temp$adjr2, BIC=temp$bic, CP=temp$cp)))
# --------------------

# Stepwise backward
# regsubsets function to perform variable selection
M_raw_2 <- lm(Y ~ (.)^2, data=dat2)
M=regsubsets(model.matrix(M_raw_2)[,-1], as.matrix(dat2$Y[complete.cases(dat2)]) ,data=dat2, method='backward', intercept=TRUE,
             nbest=1, nvmax=5)
temp <- summary(M)
Var <- colnames(model.matrix(M_raw_2))
M_select <- apply(temp$which,1,function(x)paste0(Var[x], collapse="+"))
kable(data.frame(cbind(model=M_select, adjR2=temp$adjr2, BIC=temp$bic, CP=temp$cp)))
# ----------------------------------------------------------------------
temp=summary(M)
summary(M)
kable(data.frame(cbind(model=M, adjR2=temp$adjr2, BIC=temp$bic)), caption='Model Summary')
M_main <- lm(Y~(.), data=dat2)
temp <- summary(M_main)
kable(temp$coefficients[abs(temp$coefficients[,4]) <= 0.001,], caption = 'Sig Coefficients')
M_2nd <- lm(Y~(.)^2,data=dat2)
temp <- summary(M_2nd)
kable(temp$coefficients[abs(temp$coefficients[,4]) <= 0.001,], caption = '2nd Interaction')
M_2stage <- lm(Y~(.)^3,data=dat2)
temp <- summary(M_2stage)
temp$coefficients[abs(temp$coefficients[,3])>=4,]

test_1 <- lm(Y~(E1+E4+E6+E6:G9), data=dat2) 
summary(test_1)

test_2 <- lm(Y~(E1+E4+E6:G9), data=dat2) #dat2 R^2=0.4422 F=265 t>=14//dat1 R^2=0.4857 F=315.4 t>=15//
summary(test_2)   #dat_a R^


x = c(125,123,117,123,115,112,128,118,124,111,116,109,125,120,113,123,112,118,121,118,122,115,105,118,131)

qqnorm(x)

shapiro.test(x)

x1 = c(1203,830,372,346,321,244,163,148,95,87,81,69,47,41,37,29,29,26,26,24,23,17,12,5,5,1)
x11 = log(x1)
mean(x11)
sd(x11)
qqnorm(x11)
qqline(x11)

x2 = c(2746,1698,1656,978,703,489,430,334,303,275,275,255,243,201,199,130,119,118,115,92,41,33,31,18,8,4)
x22 = log(x2)
mean(x22)
sd(x22)
qqnorm(x22)
qqline(x22)

x <- qf(c(0.025,0.975),25,25)
x
x1 = c(1203,830,372,346,321,244,163,148,95,87,81,69,47,41,37,29,29,26,26,24,23,17,12,5,5,1)
x11 = log(x1)
x2 = c(2746,1698,1656,978,703,489,430,334,303,275,275,255,243,201,199,130,119,118,115,92,41,33,31,18,8,4)
x22 = log(x2)
var.test(x11,x22,alternative = c("two.sided"), con =0.95)


x1 = c(488,478,480,426,440,410,458,460)
qqnorm(x1)
qqline(x1)

x2 = c(484,478,492,444,436,398,464,476)
qqnorm(x2)
qqline(x2)

affected <- c(488,478,480,426,440,410,458,460)
nonaffected <- c(484,478,492,444,436,398,464,476)
d <- c(affected-nonaffected)
mean(d)
sd(d)
t.test(affected,nonaffected,paired=TRUE, alternative = "two.sided", conf.level = 0.9)
my_data <- data.frame(
  group = rep(c("affected", "nonaffetced", each = 8), glaucoma = c(affected, nonaffected) )
)


summary(my_data)

alpha = 0.025
qf(alpha,25,25)
qf(1-alpha,25,25)

alpha = 0.05
qf(1-alpha,9,4)
qf(1-alpha,4,9)
a1 <- c(12.0129,12.0072,12.0064,12.0054,12.0016,11.9853,11.9949,11.9985,12.0077,12.0061)
a2 <- c(12.0318,12.0246,12.0069,12.0006,12.0075)
var.test(a1,a2,alternative = "two.sided",conf.level = 0.9)

t.text1  <- function(m1,m2,s1,s2,n1,n2)
se <- sqrt((1/n1+1/n2)*((n1-1)*s1^2+(n2-1)*s2^2)/(n1+n2-2))
t <- (m1-m2)/se
df <- n1+n2-2

library(BSDA)



alpha = 0.05
qf(alpha,20,20)
qf(1-alpha,20,20)

alpha = 0.05
qf(1-alpha,9,4)

install.packages("BSDA")
library(BSDA)
tsum.test(mean.x=4.46-4.44, n.x=10, mu = 0, s.x =0.6, conf.level = 0.95)
t.test(x)

x <- rnorm(20)
tsum.test(mean(x),sd(x),n.x=20)
t.test(x)

z.test(mean=17410-19817, sd.x=6423, sd.y=7123, conf.level = 0.95)


zsum.test(mean.x = 17410, mean.y =19817, sigma.x = 6423, sigma.y = 7123, conf.level = 0.95, n.x = 72, n.y =72)

zsum.test(mean.x = 9.709, mean.y =9.840, sigma.x = 0.325, sigma.y = 0.325, conf.level = 0.95, n.x = 72, n.y =72)


tsum.test(mean.x = 4.44, s.x = 0.73, n.x = 20,
     mean.y =4.46, s.y=0.64, n.y = 20)

tsum.test(mean.x = 4.46-4.44, n.x = 10, s.x = 0.4)

tsum.test(mean.x = 4.44-4.46, n.x = 10, s.x = 0.4)
tsum.test(mean.x = 4.46-4.44, n.x = 20, s.x = 0.4)

tsum.test(mean.x = 4.8-4.44, n.x = 20, s.x = 0.4)



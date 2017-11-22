library(sp)
library(lattice)

data(meuse)
coordinates(meuse)=~x+y
data(meuse.grid)
gridded(meuse.grid) = ~x+y

# Cross validation
loo <- krige.cv(log(zinc) ~ sqrt(dist), meuse, meuse.grid, 
              model = vgm(.149, "Sph", 700, .0674), 
              nmax = 40, nmin = 20, maxdist = 1000)
loo <- na.omit(loo)


# Does the value fall in the right PI ? (Goovaert 2001, equ. 13)
q <- seq(from=0.01, to=0.99, by=0.01)
m <- rep(NA,99)
a=0
ch <- rep(NA,length(loo$observed))

for (j in q){
  
  a=a+1
  min_q <- qnorm((1-j)/2, mean=loo$var1.pred, sd=sqrt(loo$var1.var))
  max_q <- qnorm((1+j)/2, mean=loo$var1.pred, sd=sqrt(loo$var1.var))
  for (i in 1:length(loo$observed)){
    if(loo$observed[i]>=min_q[i] & loo$observed[i]<=max_q[i]){ch[i]=1}else{ch[i]=0}
  }
  m[a] <- mean(ch, na.rm=T)
  
}

plot(q,m,xlim=c(0,1), ylim=c(0,1))
abline(0,1)

#Accuracy (Deutsch, 1997. equ. 6):
FUN <- function(x){ 
  ap =c()
  min_q <- qnorm((1-x)/2, mean=loo$var1.pred, sd=sqrt(loo$var1.var))
  max_q <- qnorm((1+x)/2, mean=loo$var1.pred, sd=sqrt(loo$var1.var))
  for (i in 1:length(loo$observed)){
    if(loo$observed[i]>=min_q[i] & loo$observed[i]<=max_q[i]){ch[i]=1}else{ch[i]=0}
  }
  m <- mean(na.omit(ch))
  ifelse(m>=x, 1, 0)
} 
Accuracy <- integrate(f=FUN, lower = 0, upper=1)$value
Accuracy 


#Precision (Deutsch, 1997. equ. 7):
FUN <- function(x){ 
  min_q <- qnorm((1-x)/2, mean=loo$var1.pred, sd=sqrt(loo$var1.var))
  max_q <- qnorm((1+x)/2, mean=loo$var1.pred, sd=sqrt(loo$var1.var))
  for (i in 1:length(loo$observed)){
    if(loo$observed[i]>=min_q[i] & loo$observed[i]<=max_q[i]){ch[i]=1}else{ch[i]=0}
  }
  m <- mean(na.omit(ch))
  
  ifelse(m >= x, 1*(m-x), 0*(m-x))
  
} 
Precision <- 1- 2*integrate(f=FUN, lower = 0, upper=1)$value
Precision 


# Goodness (Deutsch, 1997. equ. 8):
FUN <- function(x){ 
  min_q <- qnorm((1-x)/2, mean=loo$var1.pred, sd=sqrt(loo$var1.var))
  max_q <- qnorm((1+x)/2, mean=loo$var1.pred, sd=sqrt(loo$var1.var))
  for (i in 1:length(loo$observed)){
    if(loo$observed[i]>=min_q[i] & loo$observed[i]<=max_q[i]){ch[i]=1}else{ch[i]=0}
  }
  m <- mean(na.omit(ch))
  
  ifelse(m >= x, (3*1-2)*(m-x), (3*0-2)*(m-x))
  
} 
Goodness <- 1- integrate(f=FUN, lower = 0, upper=1)$value
Goodness 


# Spread of the predictions (Goovaerts, 2001. equ. 15)
q <- seq(from=0.01, to=0.99, by=0.01)
W <- rep(NA,99)
a=0
ch <- rep(NA,length(loo$observed))
t <- rep(NA,99)

N <- length(loo$observed)
for (j in q){
  
  a=a+1
  min_q <- qnorm((1-j)/2, mean=loo$var1.pred, sd=sqrt(loo$var1.var))
  max_q <- qnorm((1+j)/2, mean=loo$var1.pred, sd=sqrt(loo$var1.var))
  for (i in 1:length(loo$observed)){
    if(loo$observed[i]>=min_q[i] & loo$observed[i]<=max_q[i]){ch[i]=1}else{ch[i]=0}
  }
  t[a] <- mean(ch, na.rm=T)
  W[a] <- (1/(N*t[a]))*sum(ch*(max_q-min_q))
  
}

plot(q,W,xlim=c(0,1), ylim=c(0,2)) # change the scale of the Y-axis if needed














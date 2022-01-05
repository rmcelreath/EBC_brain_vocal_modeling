# sythnetic simulation of Pan brain/vocal data

library(rethinking)

n_individuals <- 40

# age and sex have no observed influences
age <- runif(n_individuals,1,20)
sex <- sample(1:2,size=n_individuals,replace=TRUE)

# simulate a contextual variable that can influence both B and V
X <- rbern(n_individuals)
bXB <- 0
bXV <- 0

# brain structure influenced by age and sex
# assume logistic relationship with age for each sex
b_rate <- c(0.9,1.1)
alpha <- 7
# curve( inv_logit( b_rate[1]*(x-5) ) , from=0 , to=20 , ylim=c(0,1) )
brain <- inv_logit( rnorm(n_individuals, b_rate[sex]*(age-alpha) + bXB*X , 0.5 ) )

# vocal behavior influenced by age,sex,brain
# simulate count of unique utterences
# assume for dramatic effect that brain has NO influence
v_rate <- c(0,0)
v_rate <- c(0,1)
delta <- 0.75
lambda <- brain^v_rate[sex] * age^(delta + bXV*X)
vocal <- rpois(n_individuals, lambda )

# plot pairwise associations
#pairs( ~ brain + vocal + age + sex , col=2 , lwd=2 , cex=1.5 )

# estimate

dat <- list(
    n_individuals = n_individuals,
    A = age,
    S = sex,
    logitB = logit(brain),
    V = vocal,
    X = X
)

formula <- alist(
    # V model
     V|V>=0 ~ poisson(lambda),
     lambda <- A^delta * inv_logit(logitB)^(gamma[S] + bXV*X),
    # B model
     logitB ~ normal(mu,sigma),
     mu <- beta[S]*(A-alpha) + bXB*X,
    # priors
     delta ~ normal(0.5,1),
     vector[2]:gamma ~ normal(0,1),
     alpha ~ normal(5,1),
     vector[2]:beta ~ normal(0,1),
     sigma ~ exponential(1),
     bXV ~ normal(0,0.5),
     bXB ~ normal(0,0.5)
)

m1 <- ulam(
    formula,
    data=dat, 
    constraints=list(delta="lower=0",gamma="lower=0",alpha="lower=0"),
    chains=3 , cores=3 , sample=TRUE )

# now with missing data

# RB as function of age
RB <- rbern(n_individuals,prob=inv_logit(age*0.1))
Bo <- brain
Bo[RB==1] <- NA # missing

# RV as random, 10% missing
RV <- rbern(n_individuals,prob=ifelse(RB==1,0.2,0.9))
Vo <- vocal
Vo[RV==1] <- NA # missing

dat2 <- list(
    n_individuals = n_individuals,
    A = age,
    S = sex,
    logitB = logit(Bo),
    V = ifelse( is.na(Vo) , -10 , Vo ),
    X = X
)

m2 <- ulam(
    formula,
    data=dat2, 
    constraints=list(delta="lower=0",gamma="lower=0",alpha="lower=0"),
    chains=3 , cores=3 , sample=TRUE )

#blank(w=3,h=0.8)
yoff <- (-0.3)

par(mfrow=c(1,2))
plot(precis(m1,2,pars=c("gamma","beta","delta")))
points( delta , 1+yoff , col=2 , pch=3 , lwd=3 )
points( b_rate[2] , 2+yoff , col=2 , pch=3 , lwd=3)
points( b_rate[1] , 3+yoff , col=2 , pch=3 , lwd=3)
points( v_rate[2] , 4+yoff , col=2 , pch=3 , lwd=3)
points( v_rate[1] , 5+yoff , col=2 , pch=3 , lwd=3)

plot(precis(m2,2,pars=c("gamma","beta","delta")))
points( delta , 1+yoff , col=2 , pch=3 , lwd=3 )
points( b_rate[2] , 2+yoff , col=2 , pch=3 , lwd=3)
points( b_rate[1] , 3+yoff , col=2 , pch=3 , lwd=3)
points( v_rate[2] , 4+yoff , col=2 , pch=3 , lwd=3)
points( v_rate[1] , 5+yoff , col=2 , pch=3 , lwd=3)

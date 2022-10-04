library(foreign)
lead <- read.dta("https://wps.pearsoned.com/wps/media/objects/11422/11696965/data3eu/lead_mortality_dta")
#computing correlation coefficients to to find two variables for my further analysis
#removing state and city as those are non-numeric variables to compute correlation coef. lead_1 <- subset(lead, select = -c(state, city))
round(cor(lead_1), 2)

#using MatchBalance function to check the balance of the two variables across treatment #and control groups.
mb <- MatchBalance(lead ~ age + temperature, data = lead,
print.level = 1)

#plotting the balance
mb1 <- matchit(lead ~ age + temperature, data = lead, method = "full") plot(mb1, type = "qq", which.xs = c("age", "temperature"))
		 	 
#finding the means of the outcome (the infrant mortality rate) divided by those who receive the treatment or not (lead)
means <- aggregate(lead$infrate, list(lead$lead), FUN=mean)
means

#calculating the Prima Facie treatment effect by simply subtracting the means of the outcomes of the control group from the treated group
Prima_Facie_TE <- (means[2,] - means[1,])
Prima_Facie_TE

#propensity score model
prop1 <- glm(lead ~ age + temperature, data = lead, family = "binomial")

#defining variables
X <- prop1$fitted Y <- lead$infrate Tr <- lead$lead

#matching the treated and control groups
mout <- Match(Y = Y, Tr = Tr, X=X, M=1, replace = FALSE) summary(mout)
 Checking the level of balance
mb2 <- MatchBalance(lead ~ age + temperature, data = lead, match.out = mout, nboots = 1000)

#creating a visualisation to show the balance
plot1 <- bal.plot(mout, formula = lead ~ age + temperature, data=lead, var.name = "age") 
plot2 <- bal.plot(mout, formula = lead ~ age + temperature, data=lead, var.name = "temperature") 
grid.arrange(plot1, plot2, nrow=2)

#finding means of the outcomes for the treatment and control groups
mean_tps = mean(lead$infrate[mout$index.treated]) mean_cps= mean(lead$infrate[mout$index.control])
Mean_tps

#calculating the treatment effect
te_ps= mean_tps - mean_cps te_ps

#the treatment effect I found closely resembles the Prima Facie treatment effect with a minor change from 0.22 to 0.21

#using the function psens which checks for the sensitivity, and setting Gamma values myself
library(rbounds)
psens(mout, Gamma = 2, GammaInc=.1)

#storing variables
library(rgenoud)
X <- cbind(lead$age, lead$temperature)
Y <- lead$infrate
Tr <- lead$lead

#running a genetic matching model
genout <- GenMatch(Tr = lead$lead, X = X, M=1, replace = FALSE)
 mout2 <- Match(Y = Y,X = X, Tr = Tr, Weight.matrix=genout)

#checking the level of balance
mb3 <- MatchBalance(lead ~ age + temperature, data = lead, match.out = mout2, nboots = 1000)

#creating a visualisation to show the balance
plot3 <- bal.plot(mout2, formula = lead ~ age + temperature, data=lead, var.name = "age")
plot4 <- bal.plot(mout2, formula = lead ~ age + temperature, data=lead, var.name = "temperature")
grid.arrange(plot3, plot4, nrow=2)

#running a genetic matching model with a narrow and a wide caliper as per instructions
genout1 <- GenMatch(Tr = lead$lead, X = X, estimand="ATT", M=1, replace = FALSE, caliper = c(1e-2,1e5))
 mout3 <- Match(Y = Y,X = X, Tr = Tr, estimand="ATT", Weight.matrix=genout1)

#checking the level of balance
mb3 <- MatchBalance(lead ~ age + temperature, data = lead, match.out = mout3, nboots = 1000)

#running genetic matching
genout3 <- GenMatch(Tr = lead$lead, X = X, estimand="ATT", M=1, replace = FALSE)

#instead of setting a caliper like in the previous step, I use the "exact" option
mout4 <- Match(Y = Y,X = X, Tr = Tr, estimand="ATT", exact = TRUE, Weight.matrix=genout3) # Checking the level of balance
mb3 <- MatchBalance(lead ~ age + temperature, data = lead, match.out = mout4, nboots = 1000)

#finding means of the outcomes for the treatment and control groups
mean_tgm = mean(lead$infrate[mout2$index.treated]) mean_cgm = mean(lead$infrate[mout2$index.control])
Mean_tgm

#calculating the treatment effect
te_gm= mean_tgm - mean_cgm te_gm

#the treatment effect I found doesn't resemble the previous two treatment effects and is equal to 0.036


#using the function psens which checks for the sensitivity
psens(mout2, Gamma=2, GammaInc=.1)$bounds

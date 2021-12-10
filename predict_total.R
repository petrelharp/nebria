x <- read.csv("sdm_numbers.csv")

the_lm <- glm(cbind(yosemite_occupied, yosemite_number - yosemite_occupied) ~ sdm_class,
              family='binomial', data=x)


sdm <- seq(0, 1, length.out=101)
plot(sdm, 1/(1+exp(-predict(the_lm, newdata=data.frame(sdm_class=sdm)))), type='l')

x$occupation_prop <- 1/(1 + exp(-predict(the_lm)))

sprintf("The total number predicted is %d", round(sum(x$total_number * x$occupation_prop)))

set.seed(1);
data = rpois(n=10000, lambda=5)

png("Example_Poisson_lambda5.png")
hist(data, col="grey65", xlab="Observed Expression Value", main="Poisson Distribution", breaks=0:20)
dev.off()

png("Example_Poisson_Weights.png")
hist(data, col="grey65", xlab="Observed Expression Value", main="Poisson Distribution", breaks=0:20)
hist(data[data<=2 | data >=8], col="cyan", xlab="Observed Expression Value", main="Poisson Distribution", breaks=0:20, add=T)
abline(v=2, col="blue", lwd=6)
dev.off()


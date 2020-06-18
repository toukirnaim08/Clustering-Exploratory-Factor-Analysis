residual_stats<-function(matrix) {
	residuals<-as.matrix(matrix[upper.tri(matrix)])
	large.resid<-abs(residuals) > 0.05
	numberLargeResids<-sum(large.resid)
	propLargeResid<-numberLargeResids/nrow(residuals)
	rmsr<-sqrt(mean(residuals^2))
	
	cat("Root means squared residual = ", rmsr, "\n")
	cat("Number of absolute residuals > 0.05 = ", numberLargeResids, "\n")
	cat("Proportion of absolute residuals > 0.05 = ", propLargeResid, "\n")
	#hist(residuals)
}

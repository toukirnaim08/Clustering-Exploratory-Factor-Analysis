mosthighlycorrelated <- function(mydataframe,numtoreport) 
{ 
# find the correlations 
cormatrix <- round(cor(mydataframe),4)
# set the correlations on the diagonal or lower triangle to zero, 
# so they will not be reported as the highest ones: 
diag(cormatrix) <- 0 
cormatrix[lower.tri(cormatrix)] <- 0 
# find the dimensions of the matrix, and the row names: 
numrows <- nrow(cormatrix) 
therownames <- rownames(cormatrix) 
# find the highest correlations 
sorted <- sort(abs(cormatrix),decreasing=TRUE) 
for (i in 1:numtoreport) 
{ 
corri <- sorted[i] 
# find the pair of variables with this correlation 
for (j in 1:(numrows-1)) 
{ 
for (k in (j+1):numrows) 
{ 
corrjk <- cormatrix[j,k]
if (corri == abs(corrjk)) 
{ 
rowname <- therownames[j] 
colname <- therownames[k] 
print(paste("i=",i,"variables",rowname,"and",colname,"correlation=",corrjk)) 
} 
} 
} 
} 
} 

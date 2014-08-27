# 
# Process permutations test between two paired matrices and return results
# for individual and global tests. Global test is processed on rows sum and individual tests are 
# processed for each row. The null hypothesis of this test is that there is no differences between 
# data obtained with the two paired matrices. 
#
# Input:   
#   matrix1:                             a first matrix containing data. Each column is 
#                                        paired to the corresponding column in matrix2, which is
#                                        assumed. Each row represents a sub-matrix (or a variable)
#                                        corresponding to the same nth row in matrix2, which is
#                                        assumed.
#   matrix2:                             a second matrix/matrix containing data. Each column is 
#                                        paired to the corresponding column in matrix1, which is
#                                        assumed. Each row represents a sub-matrix (or a variable)
#                                        corresponding to the same nth row in matrix1, which is
#                                        assumed.
#   confidenceLevel:                     a fraction between 0 and 1 representing the confidence level
#                                        used to test the null hypothesis. The most common value 
#                                        is 0.95, which is associated to a threshold of 0.05 (5%)
#
# Prerequisites: 
#   The matrix1 argument is a numeric matrix. If it contains NA, those will be replaced by zeros.
#   The matrix2 argument is a numeric matrix. If it contains NA, those will be replaced by zeros.
#   The dimensions of matrix1 and matrix2 must be identical.
#   The confidenceLevel must be a number between 0 and 1. 0 being the situation where every result is 
#   considered as statistically significant and 1 being the situation where every result is considered
#   as NOT statistically significant.
#
# Output: 
#   A matrix (permTest) of seven columns:   
#                           1. Row names
#                           2. Observed values of the statistic
#                           3. Values associated to the upper limit of the permutations tests
#                           4. Values associated to the observed mean with the first matrix
#                           5. Values associated to the observed mean with the second matrix
#                           6. matrix who gives significantly higher data (if there is)
#                           7. P-values assiciated to permutations tests
#
#   Warning: The global permutations test results will be the first row of the resulting matrix.
#
permutationTest <- function(matrix1, matrix2, confidenceLevel=0.95){
  
  #######################################
  # Test prerequisites
  #######################################
  
  # The matrix1 must be a matrix. If it is a vector, than transform it in a one row matrix
  if (is.vector(matrix1)) {
    matrix1=t(as.matrix(matrix1))
  }
  if (!is.matrix(matrix1)) {
    stop("The matrix1 argument must be a matrix.")
  }
  
  # The matrix2 must be a matrix. If it is a vector, than transform it in a one row matrix
  if (is.vector(matrix2)) {
    matrix2=t(as.matrix(matrix2))
  }
  if (!is.matrix(matrix2)) {
    stop("The matrix2 argument must be a matrix.")
  }
  
  # The matrix1 must contains numeric elements
  if (!is.numeric(matrix1)) {
    stop("The matrix1 argument must contains numeric elements.")
  }
  
  # The matrix2 must contains numeric elements
  if (!is.numeric(matrix2)) {
    stop("The matrix2 argument must contains numeric elements.")
  }
  
  # If there are NA in matrix1 or matrix2, then transform them into zeros
  matrix1[is.na(matrix1)] <- 0
  matrix2[is.na(matrix2)] <- 0
  
  # The dimensions of matrix1 and matrix2 must be identicals
  if ((dim(matrix1)[1] != dim(matrix2)[1]) | (dim(matrix1)[2] != dim(matrix2)[2])) {
    stop("The dimensions of matrix1 and matrix2 must be identicals.")
  }
  
  # The confidence level must be numeric
  if (!is.numeric(confidenceLevel)) {
    stop("The confidence level must be numeric.")
  }
  # The confidence level must be a fraction between 0 and 1
  if (confidenceLevel<0 | confidenceLevel>1) {
    stop("The confidence level must be a fraction between 0 and 1.")
  }
  
  ##############################################
  # Set basic parameters for permutations tests.
  ##############################################
  
  # nreps is the number of permutations that will be applyed 
  # The number of columns in matrix1 are decreased of one because the number of permutations is reduced to the unique ones 
  nreps <- 2^(dim(matrix1)[2]-1)
  
  # base is the matrix containing the two vectors (one full of '-1' and to other full of '1') that are permuted to obtain all
  # possibilities of permutation
  base <- as.data.frame(matrix(rep(c(-1,1),dim(matrix1)[2]), nr =2))
  
  # permutations is the matrix of every single possibilities of permutations
  permutations <- as.matrix(expand.grid(base))
  
  # testPer is a matrix containing half of the rows in the 'permutations' matrix,  which are the unique permutations
  testPer <- permutations[1:((dim(permutations)[1])/2),]
  
  # upperLimit gives the index of the upper limit who consists in the last value who doesn't reject the nul hypothesis 
  # according to the confidenceLevel sets at the beginning
  upperLimit <- ceiling(nreps*confidenceLevel)
  
  # Get rows sum for each column to make a global permutation test
  # Those will be added as first rows of matrix1 and matrix2
  matrixSum = apply(matrix1, 2, sum)
  matrix1=rbind(matrixSum, matrix1)
  
  matrixSum = apply(matrix2, 2, sum)
  matrix2=rbind(matrixSum, matrix2)
  
  # Set an empty matrix to contain permutations tests results
  permTest = matrix(ncol = 7, nrow=dim(matrix1)[1])
  
  # First column contains row names
  permTest[,1] = rownames(matrix1)
  
  # obs.diff is a matrix containing differences between observed paired rows (one for each matrix) 
  obs.diff <- matrix1 - matrix2
  
  # Second column contains the observed values of the statistic (mean difference between paired data).
  permTest[,2] = apply(obs.diff, 1, function(x) abs(mean(x)))
  
  # perm1 gives the upper limit of ordered statistics obtained with permutations
  perm1 = function(position){
    
    # allmeans is a matrix containing absolute values of statistics associated with each permutation processed.
    allmeans <- abs(testPer %*% position)/(dim(matrix1)[2])
    
    # Value associated to the upper limit
    up.lim = (sort(allmeans))[upperLimit]
    
    return(up.lim)
  }  
  
  # Third column contains the value associated to the upper limit
  permTest[,3] = apply(obs.diff, 1, function(x) perm1(x))
  
  # Fourth and fifth columns contain the values associated to the observed mean with the first and second matrix
  permTest[,4] = apply(matrix1, 1, mean)
  permTest[,5] = apply(matrix2, 1, mean)
  
  # perm2 gives the matrix who give higher depths
  perm2 = function(position){
    
    if (as.numeric(position[2]) > as.numeric(position[3])){
      if(as.numeric(position[4]) > as.numeric(position[5])){
        matrix.high = 1
      } else {matrix.high = 2}
    } else {matrix.high = NA}
    
    return(matrix.high)
  }
  
  # Sixth column contains the matrix who give higher depths when the difference between rows is statistically significant
  permTest[,6] = apply(permTest, 1, function(x) perm2(x))
  
  # perm3 gives the p-value associated to the permutations test
  perm3 = function(position){
    
    # allmeans is a matrix containing absolute values of statistics associated with each permutation processed
    allmeans <- abs(testPer %*% position)/(dim(matrix1)[2])
    
    p.val = length(allmeans[allmeans >= abs(mean(position))])/nreps
    
    return(p.val)
  }  
  
  # Seventh column contains the p-value associated to permutation tests, for each row
  permTest[,7] = apply(obs.diff, 1, function(x) perm3(x))

  colnames(permTest) = c("Rowname", "Obs value", "Upper limit value", "Mean matrix1", "Mean matrix2", "Higher matrix", "P-value")

  # Gives the global permutation test result
  if (is.na(permTest[1,6])){
    best="None"
  } else if (as.numeric(permTest[1,6])==1){
    best="matrix1"} else if (as.numeric(permTest[1,6])==2){
      best="matrix2"
    } 

  writeLines(paste("The global permutation gives a p-value of", permTest[1,7], ", we can see that", best,"is the matrix giving the higher mean data."))
  
  totP = dim(permTest)[1]-1
  totSignP = sum(as.numeric(permTest[-1,7])<=(1-confidenceLevel))
  
  # Gives the percentage of significant p-values
  writeLines(paste("There are", paste(round(totSignP/totP*100, digit=2),"%", sep=""),
                   "of p-values that are significants with a threshold limit of", paste((1-confidenceLevel)*100,"%.", sep="")))
  
    # Give the proportion of significant p-values associated with each case
  if(totSignP>0){writeLines(paste("From those significant p-values, ", paste(round(sum(as.numeric(permTest[-1,6]) == 1, na.rm=TRUE)/totSignP*100, 
                   digit=2),"%", sep=""), " are associated to higher data with matrix1 and ", 
                   paste(round(sum(as.numeric(permTest[-1,6]) == 2, na.rm=TRUE)/totSignP*100, digit=2),"%", sep=""), 
                   " are associated to higher data with matrix2.", sep=""))}
  
  # Return the permTest matrix containing permutation results
  return(permTest)
}  
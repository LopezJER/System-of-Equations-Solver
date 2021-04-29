#' This R script solves systems of linear equations through the Gaussian and
#' Gauss-Jordan methods of elimination. Each method is encapsulated in a function
#' so that the results of each method may be compared. On the one hand, 
#' Gaussian elimination converts the system into its row-echelon form. Hence, 
#' a phase for backward substitution is needed for the solution set. 
#' On the other hand, Gauss-Jordan elimination converts the system into its 
#' reduced row-echelon form by normalizing each row. The constants column 
#' of the reduced row-echelon matrix then contains the solution set.
#
#' @author Jose Enrique Lopez
#' @date 12-10-20 20:50


# NOTE: Uncomment install method below if rstudioapi library does not exist
# install.packages("rstudioapi") 

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
source ("AugCoeffMatrix.R")

GaussJordanElimination <- function (system){
  matrix_data <- AugCoeffMatrix(system);
  n <- length(matrix_data$variables)
  mat <- matrix_data$augcoeffmatrix

  for (i in 1:n){
    
    if (i!=n){
      pivotRow <- which(mat == max (abs(mat[i:n, i])), arr.ind = T)
      if (mat[pivotRow[1], i] == 0) break;
      mat[c(pivotRow[1], i),] <- mat[c(i, pivotRow[1]),]
    }

    mat[i, ] = mat[i, ] / mat[i, i]
    for (j in 1 : n){
      if (i == j) next
      normalizedRow = mat[j, i] * mat[i, ]
      mat[j, ] = mat[j, ] - normalizedRow
    }

  }
  solutionSet = mat[, ncol(mat)]
  
  result = list(solutionSet = solutionSet, variables = matrix_data$variables, matrix = mat)
  return(result)
}
  
GaussianElimination <-function (system){
  matrix_data <- AugCoeffMatrix(system);
  n <- length(matrix_data$variables)
  mat <- matrix_data$augcoeffmatrix
  for (i in 1:(n-1)){
      pivotRow <- which(mat == max (abs(mat[i:n, i])), arr.ind = T)
      if (mat[pivotRow[1], i] == 0) break;
      mat[c(pivotRow[1], i),] <- mat[c(i, pivotRow[1]),]
      
    for (j in (i+1) : n){
      pivotElement = mat[i, i]
      multiplier = mat [j, i]/pivotElement
      normalizedRow = multiplier * mat[i, ]
      mat[j, ] = mat[j, ] - normalizedRow
      
    }
    
    solutionSet = numeric(length(matrix_data$variables))
    solutionSet[length(solutionSet)] = mat[nrow(mat), ncol(mat)]/mat[nrow(mat), ncol(mat)-1]
    for (i in (n-1):1){
      solutionSet[i] = (mat[i, n+1] - sum(mat[i, (i+1):n] * solutionSet[(i+1): n])) / mat[i,i]
    }
  }
  result = list(solutionSet = solutionSet, variables = matrix_data$variables, matrix = mat)
  return(result)
}





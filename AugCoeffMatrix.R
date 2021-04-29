#' This R script takes a system of equations and converts it
#' into its matrix representation. The system is deparsed and each equation
#' is tokenized into its corresponding addends and factors. Each token is then
#' parsed to arrive at the constants, coefficients, and variables of the system.
#' This information is  compiled into a list containing 1) an augmented coefficient
#' matrix representing the system and 2) a vector of the variables. 
#' This list is returned once the system of equations has been successfully parsed;  
#' otherwise, a value of NA is returned.
#' 
#' @author Jose Enrique Lopez
#' @date 10-10-20 19:49


getConsts <- function (addends){

  constants <- c()
  for (i in 1:length(addends)){
    constant <- -1*as.numeric(trimws(addends[[i]][[length(addends[[i]])]]))
    constants<-append(constants, constant)
  }
  return (constants)
}

getVars <- function (addends){

  matrix_data = parseAddends(addends)
  if (length(matrix_data)){
    return(matrix_data$vars)
  } else{
    print ("Variables are not the same")
    return (matrix_data)
  }
}

getCoeffs <-function (addends){
  matrix_data = parseAddends(addends)
  return(matrix_data$coeffs)
}

parseAddends <- function (addends){
  vars <- c()
  coeffs <- c()
  
  for (i in 1: length (addends)){
    eq_vars <- c ()
    eq_coeffs <- c()
    
    for (j in 1: length(addends[[i]][-length(addends[[i]])])){
      factors = strsplit(addends[[i]][[j]], "\\*")
      
      coeff = as.numeric(trimws(factors[[1]][[1]]))
      var = trimws(factors[[1]][[2]])
      
      eq_vars <- append(eq_vars, var)
      eq_coeffs <- append (eq_coeffs, coeff)
    }
    
    if (!length(vars)){
      vars <- sort(eq_vars)
    } 
    if (!setequal(eq_vars, vars)){
      return (list())
    }
    
    sorted_eq_coeffs <-numeric (length(eq_vars))
    for (i in 1:length(eq_vars)){
      degree = strtoi(strsplit(eq_vars[i], "")[[1]][2])
      sorted_eq_coeffs[degree] <- eq_coeffs[i]
    }
    coeffs <- append (coeffs, sorted_eq_coeffs)
  }
  return (list(vars = vars, coeffs = coeffs))
}

deparseSystem <-function(system){
  equations<-c()
  
  for (i in 1:length(system)){
    eq <- deparse(system[[i]])
    eq_str <- c()
    for (i in 2: length (eq)){
      eq_str <- paste(eq_str, eq[[i]])
    }
    equations <- append(equations, eq_str) 
  }
  
  return (equations)
}

buildMatrix <-function (equations){
  addends <- strsplit (equations, "\\+")
  vars <- getVars(addends)
  if (!length(vars)) return (NA)
  
  consts <- getConsts(addends)
  coeffs <- getCoeffs(addends)
    
  mat <- matrix (coeffs, nrow = length(equations), ncol=length(vars), byrow = TRUE)
  
  rownames(mat)<-(1:length(equations))
  mat <-cbind(mat, consts)
  
  colnames(mat)<-c(vars, "RHS")
  
  return (list(variables = vars, augcoeffmatrix = mat))
}


AugCoeffMatrix <-function (system){
  equations <- deparseSystem (system)
  matrix_data <- buildMatrix(equations)
  return (matrix_data)

}






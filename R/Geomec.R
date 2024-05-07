#Transform text input into a matrix function
text2matrix <- function(input_text) {
  values <- strsplit(input_text, "\\s+")
  if (length(unlist(values)) %% 3 != 0) {
    stop("The number of values must be a multiple of 3.")
  }
  values <- as.numeric(unlist(values))
  result_matrix <- matrix(values, ncol = 3, byrow = TRUE)
  return(result_matrix)
}

#Mirror point A function
mirror_A <- function(result_matrix) {
  A <- result_matrix[1, ]
  B <- result_matrix[2, ]
  C <- result_matrix[3, ]
  D <- result_matrix[4, ]
  cross_product <- function(u, v) {
    return(c(u[2]*v[3] - u[3]*v[2], u[3]*v[1] - u[1]*v[3], u[1]*v[2] - u[2]*v[1]))
  }
  BC <- C - B
  BD <- D - B
  n <- cross_product(BC, BD)
  d <- abs(sum((A - B) * n) / sqrt(sum(n^2)))
  n_norm <- n / sqrt(sum(n^2))
  A_Prime <- A + 2 * d * n_norm
  result_matrix_A_Prime <- rbind(result_matrix, A_Prime)
  assign("A_Prime", A_Prime, envir = .GlobalEnv)
  assign("result_matrix_A_Prime", result_matrix_A_Prime, envir = .GlobalEnv)
  return(A_Prime)
}

#Create point M function
point_M <- function(result_matrix_A_Prime) {
  A <- result_matrix_A_Prime[1, ]
  B <- result_matrix_A_Prime[2, ]
  C <- result_matrix_A_Prime[3, ]
  D <- result_matrix_A_Prime[4, ]
  A_Prime <- result_matrix_A_Prime[5, ]
  AB <- A_Prime - A
  projection <- (sum((B - A) * AB) / sum(AB^2)) * AB
  M <- A + projection
  MB <- B - M
  assign("point_M_coordinates", M, envir = .GlobalEnv)
  return(list(M = M, MB = MB))
}

#Calculate euclidean distance in mm between points function - Godinho et al., 2020 protocol
CalcEstDistancesR <- function(input_text) {
  text2matrix <- function(input_text) {
    values <- strsplit(input_text, "\\s+")
    if (length(unlist(values)) %% 3 != 0) {
      stop("The number of values must be a multiple of 3.")
    }
    values <- as.numeric(unlist(values))
    result_matrix <- matrix(values, ncol = 3, byrow = TRUE)
    return(result_matrix)
  }
  est_distances_matrix <- text2matrix(input_text)
  M_coordinates <- point_M_coordinates
  est_distances_matrix <- rbind(est_distances_matrix, M_coordinates)
  distances <- c(
    "Corpus length" = sqrt(sum((est_distances_matrix[1, ] - est_distances_matrix[11, ])^2)),
    "Ramus lateral length" = sqrt(sum((est_distances_matrix[6, ] - est_distances_matrix[13, ])^2)),
    "Ramus width" = sqrt(sum((est_distances_matrix[10, ] - est_distances_matrix[21, ])^2)),
    "Dental arch breadth" = sqrt(sum((est_distances_matrix[2, ] - est_distances_matrix[22, ])^2))
  )
  assign("est_distances_matrix", est_distances_matrix, envir = .GlobalEnv)
  return(list(distances = distances))
}

#Calculate euclidean distance in mm between points function - Standard protocol
CalcEstDistances <- function(input_text) {
  text2matrix <- function(input_text) {
    values <- strsplit(input_text, "\\s+")
    if (length(unlist(values)) %% 3 != 0) {
      stop("The number of values must be a multiple of 3.")
    }
    values <- as.numeric(unlist(values))
    result_matrix <- matrix(values, ncol = 3, byrow = TRUE)
    return(result_matrix)
  }
  est_distances_matrix <- text2matrix(input_text)
  M_coordinates <- point_M_coordinates
  B_coordinates <- result_matrix[2, ]
  est_distances_matrix <- rbind(est_distances_matrix, M_coordinates, B_coordinates)
  distances <- c(
    "Corpus length" = sqrt(sum((est_distances_matrix[1, ] - est_distances_matrix[2, ])^2)),
    "Ramus lateral length" = sqrt(sum((est_distances_matrix[3, ] - est_distances_matrix[4, ])^2)),
    "Ramus width" = sqrt(sum((est_distances_matrix[5, ] - est_distances_matrix[6, ])^2)),
    "Dental arch breadth" = sqrt(sum((est_distances_matrix[7, ] - est_distances_matrix[8, ])^2))
  )
  assign("est_distances_matrix", est_distances_matrix, envir = .GlobalEnv)
  return(list(distances = distances))
}

#Extra - Calculate the euclidean distance in mm between two points function
PointDistancemm <- function(point_A, point_B) {
  distance <- sqrt(sum((point_A - point_B)^2))
  cat("The distance between point A and point B is:", distance, "\n")
}

# Extra - Transform text input into R vector function
text2vector <- function(input_text) {
  lines <- strsplit(input_text, "\n", fixed = TRUE)[[1]]
  values <- numeric(0)
  for (line in lines) {
    coords <- strsplit(line, " ")[[1]]
    coords <- as.numeric(trimws(coords))
    values <- c(values, coords)
  }
  return(values)
}

#Extra - Calculate the cross-product of two vector function
cross_product <- function(u_vector, v_vector) {
  return(c(u[2]*v[3] - u[3]*v[2], u[3]*v[1] - u[1]*v[3], u[1]*v[2] - u[2]*v[1]))
}


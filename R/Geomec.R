#Transform text input into a matrix with all the points labeled
getmatrix <- function(input_text) {
  values <- strsplit(input_text, "\\s+")
  if (length(unlist(values)) %% 3 != 0) {
    stop("The number of values must be a multiple of 3.")
  }
  values <- as.numeric(unlist(values))
  result_matrix <- matrix(values, ncol = 3, byrow = TRUE)
  row.names(result_matrix) <- c("A", "B", "C", "D", "Gonion", "CP3", "Cond_Midpoint", "Ramus_Post", "Ramus_Root")
  return(result_matrix)
}

#ALEX - Transform text input into a matrix function - Godinho et al., 2020 protocol
t2m <- function(input_text) {
  values <- strsplit(input_text, "\\s+")
  if (length(unlist(values)) %% 3 != 0) {
    stop("The number of values must be a multiple of 3")
  }
  values <- as.numeric(unlist(values))
  result_matrix <- matrix(values, ncol = 3, byrow = TRUE)
  row.names(result_matrix) <- c("A", "B", "C", "D")
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
  n_norm <- n / sqrt(sum(n^2))
  d <- sum((A - B) * n_norm)
  A_Line <- A - 2 * d * n_norm
  result_matrix <- rbind(result_matrix, A_Line)
  assign("A_Line", A_Line, envir = .GlobalEnv)
  assign("result_matrix", result_matrix, envir = .GlobalEnv)
  return(A_Line)
}

#Create point M function
point_M_Vs <- function(result_matrix) {
  A <- result_matrix[1, ]
  B <- result_matrix[2, ]
  C <- result_matrix[3, ]
  D <- result_matrix[4, ]
  A_Line <- result_matrix[10, ]
  AA_Line <- A_Line - A
  projection <- (sum((B - A) * AA_Line) / sum(AA_Line^2)) * AA_Line
  M <- A + projection
  MB <- B - M
  result_matrix <- rbind(result_matrix, M, MB, AA_Line)
  assign("point_M_coordinates", M, envir = .GlobalEnv)
  assign("result_matrix", result_matrix, envir = .GlobalEnv)
  return(list(M = M, MB = MB, AA_Line = AA_Line))
}

#ALEX - Create point M function - Godinho et al., 2020 protocol
pMV <- function(result_matrix) {
  A <- result_matrix[1, ]
  B <- result_matrix[2, ]
  C <- result_matrix[3, ]
  D <- result_matrix[4, ]
  A_Line <- result_matrix[5, ]
  AA_Line <- A_Line - A
  projection <- (sum((B - A) * AA_Line) / sum(AA_Line^2)) * AA_Line
  M <- A + projection
  MB <- B - M
  result_matrix <- rbind(result_matrix, M, MB, AA_Line)
  assign("point_M_coordinates", M, envir = .GlobalEnv)
  assign("result_matrix", result_matrix, envir = .GlobalEnv)
  return(list(M = M, MB = MB, AA_Line = AA_Line))
}

#ALEX - Calculate euclidean distance in mm between points function - Godinho et al., 2020 protocol
distR <- function(input_text) {
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
    "Mandibular superior length" = sqrt(sum((est_distances_matrix[2, ] - est_distances_matrix[22, ])^2))
  )
  Estandardization_distances <- matrix(distances, nrow = 1, ncol = 4)
  colnames(Estandardization_distances) <- c("Corpus length", "Ramus lateral length", "Ramus width", "Mandibular superior length") 
  assign("Estandardization_distances", Estandardization_distances, envir = .GlobalEnv)
  return(Estandardization_distances)
}

#Calculate euclidean distance in mm between points function - Standard protocol
CalcEstDistances <- function(input_matrix) {
  distances <- c(
    "Corpus length" = sqrt(sum((input_matrix[3, ] - input_matrix[5, ])^2)),
    "Ramus lateral length" = sqrt(sum((input_matrix[6, ] - input_matrix[7, ])^2)),
    "Ramus width" = sqrt(sum((input_matrix[8, ] - input_matrix[9, ])^2)),
    "Mandibular superior length" = sqrt(sum((input_matrix[2, ] - input_matrix[11, ])^2))
  )
  Estandardization_distances <- matrix(distances, nrow = 1, ncol = 4)
  colnames(Estandardization_distances) <- c("Corpus length", "Ramus lateral length", "Ramus width", "Mandibular superior length") 
  assign("Estandardization_distances", Estandardization_distances, envir = .GlobalEnv)
  return(Estandardization_distances)
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

#Extra - Transform text input into a matrix
text2matrix <- function(input_text) {
  values <- strsplit(input_text, "\\s+")
  if (length(unlist(values)) %% 3 != 0) {
    stop("The number of values must be a multiple of 3")
  }
  values <- as.numeric(unlist(values))
  result_matrix <- matrix(values, ncol = 3, byrow = TRUE)
  return(result_matrix)
}

#Extra - Fill NA cells of a column preserving pre-exinting mean 
fill_na_mean_pres <- function(df, column) {
  # Calculate the mean ignoring NA values
  initial_mean <- mean(df[[column]], na.rm = TRUE)
  # Count the number of NA values
  num_NA <- sum(is.na(df[[column]]))
  # Calculate the SD ignoring NA cells
  sd_value <- sd(df[[column]], na.rm = TRUE)
  #Calculate 95% CI
  n <- sum(!is.na(df[[column]]))  # Tamaño de la muestra sin NAs
  error_margin <- qt(0.975, df = n - 1) * sd_value / sqrt(n)
  conf_interval_min <- initial_mean - error_margin
  conf_interval_max <- initial_mean + error_margin
  # Calculate the sum that the values must have to maintain the mean
  required_sum <- initial_mean * length(df[[column]]) - sum(df[[column]], na.rm = TRUE)
  # Generate random values within the allowed range
  set.seed(123) # para reproducibilidad
  na_values <- runif(num_NA, min = conf_interval_min, max = conf_interval_max)
  na_values <- na_values / sum(na_values) * required_sum
  # Adjust values to be within the allowed range (existing min and max)
  while (any(na_values < conf_interval_min) || any(na_values > conf_interval_max)) {
    na_values[na_values < conf_interval_min] <- runif(sum(na_values < conf_interval_min), min = conf_interval_min, max = conf_interval_max)
    na_values[na_values > conf_interval_max] <- runif(sum(na_values > conf_interval_max), min = conf_interval_min, max = conf_interval_max)
    na_values <- na_values / sum(na_values) * required_sum
  }
  # Fill the NA values with the generated values
  df[[column]][is.na(df[[column]])] <- na_values
  # Verify the new mean
  final_mean <- mean(df[[column]])
  # Display the means for verification
  print(paste("Initial mean:", initial_mean))
  print(paste("Final mean:", final_mean))
  if (initial_mean == final_mean) {
    print("Means are equal!")
  }
  #Create the new df with NA cells filled
  assign("NA_filled_df", df, envir = .GlobalEnv)
  # Return the dataset with filled values
  return(df)
}

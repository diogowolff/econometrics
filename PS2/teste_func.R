
design_matrix_base_transpose = t(matrix(1, nrow(clean_df), ncol = 2))

fast_beta_estimator = function(point) {
  Kv = map_dbl(clean_df$age, ~ K_quartic(.x-point)/3)
  
  design_matrix_base_transpose[2, ] = t(clean_df$age-point)
  
  solve(design_matrix_base_transpose %*% t(design_matrix_base_transpose * Kv)) %*% 
    colSums((t(design_matrix_base_transpose * (Kv*clean_df$chiptime))))
}

fast_betas <- lapply(x0_grid_q10, fast_beta_estimator)


################### abaixo = rapido e funciona


left_side_for_xi = function(xi, x0) {
  K = K_quartic((xi-x0)/3)
  
  return(matrix(c(K, K*(xi-x0), K*(xi-x0), K*(xi-x0)^2), nrow=2))
}

left_side_over_all_xi = function(x0) {
  Reduce('+', map(clean_df$age, ~ left_side_for_xi(.x, x0 = x0)))
}

right_side_for_xi = function(xi, yi, x0) {
  K = K_quartic((xi-x0)/3)
  
  return(matrix(c(K*yi, K*yi*(xi-x0))))
}

right_side_over_all_xi = function(x0) {
  Reduce('+', map2(clean_df$age, clean_df$chiptime, ~ right_side_for_xi(.x, .y, x0 = x0)))
}

somewhat_fast_beta_estimator = function(x0) {
  solve(left_side_over_all_xi(x0)) %*% right_side_over_all_xi(x0)
}


beta_q10 = lapply(x0_grid_q10, somewhat_fast_beta_estimator)


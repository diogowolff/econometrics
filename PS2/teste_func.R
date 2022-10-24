
design_matrix_base_transpose = t(matrix(1, nrow(clean_df), ncol = 2))

fast_beta_estimator = function(point) {
  Kv = map_dbl(clean_df$age, ~ K_quartic(.x-point)/3)
  
  design_matrix_base_transpose[2, ] = t(clean_df$age-point)
  
  solve(design_matrix_base_transpose %*% t(design_matrix_base_transpose * Kv)) %*% 
    colSums((t(design_matrix_base_transpose * (Kv*clean_df$chiptime))))
}

fast_betas <- lapply(x0_grid_q10, fast_beta_estimator)

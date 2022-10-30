###########################################
# Econometrics - Problem Set 2            #
# Students: Luan Borelli and Diogo Wolff  #
###########################################


# Defining Nelder-mead parameters:
reflection <- 1
expansion <- 1.7
contraction <- 0.3 
shrinkage <- 0.6

# Defining tolerance variables:
max_run <- 10000 # Loop max run times.
run_count <- 0 # Loop counter.
tolerance <- 10^(-6) # Tolerance.
tol <- 999999999999 # An initial tolerance value to be updated.

# Defining the Himmelblau function, to be minimized:
f <- function(theta) { 
  (theta[1]^2 + theta[2] - 11)^2 + (theta[1]+theta[2]^2 - 7)^2
}

# Defining the algorithm function.
# This function receives a list of three 2-dimensional vectors (i.e., a list of vertices)
# E.g.: list(c(2, 3), c(7, 3), c(6,4))

nelder_mead <- function(verts) {
        vertices <- verts
        while(run_count <= max_run & tol >= tolerance) {
          
          f_at_vertices <- sapply(vertices,f)
          best_pos <- which.min(f_at_vertices)
          worst_pos <- which.max(f_at_vertices)
          new_vertices <- vertices[-worst_pos]
          f_at_new_vertices <- sapply(new_vertices, f)
          
          centroid <- 1/length(new_vertices) * (unlist(new_vertices[1]) + unlist(new_vertices[2]))
          reflection_point <- (1+reflection)*centroid - reflection*unlist(vertices[worst_pos])
            
            if(f(reflection_point) < f_at_vertices[best_pos]) {
              theta_e = (1 + reflection*expansion)*centroid - reflection*expansion*unlist(vertices[worst_pos])
              if(f(theta_e) < f(reflection_point)) {
                vertices <- c(new_vertices, list(theta_e))
              } else {
                vertices <- c(new_vertices, list(reflection_point))
              }
            } else if(f(reflection_point) >= f_at_vertices[best_pos] & f(reflection_point) < f(unlist(new_vertices[which.max(f_at_new_vertices)]))) {
                vertices <- c(new_vertices, list(reflection_point))  
            } else if(f(unlist(new_vertices[which.max(f_at_new_vertices)])) <= f(reflection_point & f(reflection_point) < f_at_vertices[worst_pos])) {
                theta_c <- (1 + reflection*contraction)*centroid - reflection*contraction*unlist(vertices[worst_pos])
                if(f(theta_c) < f(reflection_point)) {
                  vertices <- c(new_vertices, list(theta_c))  
                } else {
                  tilde_theta_1 <- unlist(vertices[best_pos])
                  tilde_theta_2 <- shrinkage*unlist(vertices[best_pos]) + (1-shrinkage)*unlist(new_vertices[which.max(f_at_new_vertices)])
                  tilde_theta_3 <- shrinkage*unlist(vertices[best_pos]) + (1-shrinkage)*unlist(vertices[worst_pos])
                  vertices <- list(tilde_theta_1, tilde_theta_2, tilde_theta_3)
                } 
            } else if(f(reflection_point) >= f_at_vertices[worst_pos]) {
                theta_c2 <- (1 - contraction)*centroid + contraction*unlist(vertices[worst_pos])
                if(f(theta_c2)<f_at_vertices[worst_pos]) {
                  vertices <- c(new_vertices, list(theta_c2))
                } else {
                  tilde_theta_1 <- unlist(vertices[best_pos])
                  tilde_theta_2 <- shrinkage*unlist(vertices[best_pos]) + (1-shrinkage)*unlist(new_vertices[which.max(f_at_new_vertices)])
                  tilde_theta_3 <- shrinkage*unlist(vertices[best_pos]) + (1-shrinkage)*unlist(vertices[worst_pos])
                  vertices <- list(tilde_theta_1, tilde_theta_2, tilde_theta_3)
                }
            }
          run_count = run_count + 1
          tol = max(sapply(vertices,f)) - min(sapply(vertices,f))
          print(c("Round", run_count))
        }
        
        print('Pontos:')
        print(vertices)
        print('Valores:')
        print(format(sapply(vertices,f), scientific=F))
}

##############
# RESULTADOS #
##############

# Essa função tem quatro mínimos locais idênticos e um máximo local. 
# Vou tentar fazer o algoritmo pegar cada um desses pontos chutando diferentes vértices iniciais.

# Minima examples 
nelder_mead(list(c(-0.2, -1), c(-0.3, -0.9), c(0,-0.5))) # Pegou (3,2)
nelder_mead(list(c(2, 3), c(7, 3), c(6,4))) # Pegou (-2.8, 3.13)
nelder_mead(list(c(1, 7), c(2, 1), c(3,5))) # Pegou (3.58, -1.84)
nelder_mead(list(c(-5, -4), c(-3, -2), c(-1,-4))) # Pegou (-3.77. -3.2)

# Maxima (NÃO CONSIGO PEGAR ESSE MÁXIMO DE JEITO NENHUM)
nelder_mead(list(c(-0.27, -0.92), c(-0.26, -0.91), c(-0.28, -0.91))) # ... (-0.27, -0.92)

# What? Esse esquema de tolerância parece problemático, faz a função pegar uns pontos nonsense.

nelder_mead(list(c(1, 3), c(1, 5), c(1,1))) # ???
nelder_mead(list(c(2, 3), c(2, 3), c(1,1))) # ???


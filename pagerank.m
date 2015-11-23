# Dilyn Fullerton

0;

function k = get_k_j(A, j)
# Parameter: A, incidence matrix
# Parameter: j, node (column position)
# Return:    k, number of outgoing paths
# Purpose:   Return the number of outgoing paths for a given node
  m = rows(A);
  k0 = 0;
  
  for i = 1:m
    if A(i, j) == -1
      k0 = k0 + 1;
    endif
  endfor
  
  k = k0;
endfunction

function B = addPaths(A)
# Parameter: A, incidence matrix
# Return:    B, modified incidence matrix
# Purpose:   For each node with no outgoing paths, add a path to each other node
  n = columns(A);
  m = rows(A);
  k = [];
  
  for j = 1:n
    k(j) = get_k_j(A, j);
    
    if k(j) == 0
      for col = 1:n
        if col != j
          A = [A; zeros(1, n)];
          A(m+1, col) = 1;
          A(m+1, j) = -1;
        endif
        m = rows(A);
      endfor
    endif
    
  endfor
  
  B = A;
endfunction

function ans = isPathFromTo(A, a, b)
# Parameter: A, incidence matrix
# Parameter: a, index of (from) node
# Parameter: b, index of (to) node
# Return:    ans, either 0 or 1 (boolean)
# Purpose:   Determine whether there is a direct path from node a to node b
  m = rows(A);
  
  ans0 = 0;
  
  for i = 1:m
    if A(i, a) == -1
      if A(i, b) == 1
        ans0 = 1;
      endif
    endif
  endfor
  
  ans = ans0;
endfunction

function P = createP(A)
# Parameter: A, incidence matrix
# Return:    P, stochastic matrix of probabilities
# Purpose:   Create the matrix P, which determines the next state
  n = columns(A);
  P0 = [];
  k = [];
  
  for j = 1:n
    k(j) = get_k_j(A, j);
    for i = 1:n
      
      if i == j
        P0(i, j) = 0;
      elseif isPathFromTo(A, j, i) == 1
        P0(i, j) = 1/k(j);
      else
        # In this case this line is redundant, as I have already modified
        # the matrix to add the extra paths. I included this here anyway
        # to make this function individually more robust.
        P0(i, j) = 1/(n-1);
      endif
      
    endfor
  endfor
  
  P = P0;
endfunction

function Q = createQ(A)
# Parameter: A, incidence matrix
# Return:    Q, stochastic matrix
# Purpose:   Return the matrix Q, which assigns even probability to
#            each other node
  n = columns(A);
  Q0 = [];
  
  for i = 1:n
    for j = 1:n
      if i == j
        Q0(i, j) = 0;
      else
        Q0(i, j) = 1/(n-1);
      endif
    endfor
  endfor
  
  Q = Q0;
endfunction

function [lambda1, v1] = largestEigenvalue(A, N=100000)
# Parameter: A, matrix to evaluate
# Parameter: N, number of iterations of power method
# Return:    lambda1, the largest eigenvalue of A (in magnitude)
# Reutrn:    v1, the associated eigenvector
# Purpose:   Compute the larges eigenvalue and assoc. eigenvector of A,
#            normalized with 2-norm

  # determine random starting vector
  x_0 = rand(rows(A), 1);
  x_1 = x_0;
  
  # apply power method
  for i = 1:N
    x_1 = A*x_1;
    x_1 = x_1/norm(x_1);
  endfor

  v1 = x_1/norm(x_1, 1);
  # determine eigenvalue
  lambda1 = x_1'*A*x_1/norm(x_1)^2;
endfunction

function [lambda, x_inf] = pageRank(A, alpha=1)
# Parameter: A, incidence matrix of pages graph
# Parameter: alpha, weight to assign to regular probability matrix P.
#                   (1-alpha) is assigned to Q.
# Return:    lambda, the eigenvalue for the state vector
# Return:    x_inf, the final state vector
# Purpose:   Compute the final state vector for a given graph of pages
  A = addPaths(A); 
  P = createP(A);
  Q = createQ(A); 
  S = (1-alpha)*Q + alpha*P;
  [lambda, x_inf] = largestEigenvalue(S);

endfunction

# My name
B = [-1   0   0   0   0   0   0   0   0   1   0   0   0   0       0       0;
      0  -1   0   1   0   0   0   0   0   0   0   0   0   0       0       0;
      0   0  -1   0   0   1   0   0   0   0   0   0   0   0       0       0;
      0   0   0  -1   0   0   1   0   0   0   0   0   0   0       0       0;
      0   0   0  -1   0   0   0   0   0   1   0   0   0   0       0       0;
      0   0   0   0  -1   0   0   0   0   0   0   1   0   0       0       0;
      0   0   0   1   0  -1   0   0   0   0   0   0   0   0       0       0;
      0   0   0   0   0  -1   1   0   0   0   0   0   0   0       0       0;
      0   0   0   1   0   0  -1   0   0   0   0   0   0   0       0       0;
      0   0   0   0   0   0  -1   0   0   0   0   0   1   0       0       0;
      0   0   0   0   0   0  -1   0   0   0   0   0   0   0       1       0;
      0   0   0   0   0   0   0  -1   0   0   0   0   0   0       1       0;
      0   0   0   0   0   0   0  -1   0   0   0   0   0   0       0       1;
      0   1   0   0   0   0   0   0  -1   0   0   0   0   0       0       0;
      0   0   0   0   0   0   0   1  -1   0   0   0   0   0       0       0;
      0   0   0   0   0   1   0   0   0  -1   0   0   0   0       0       0;
      0   0   0   0   0   0   0   0   1  -1   0   0   0   0       0       0;
      0   0   0   0   0   0   0   0   0  -1   1   0   0   0       0       0;
      0   0   0   0   0   0   0   0   1   0  -1   0   0   0       0       0;
      0   0   0   0   0   0   0   0   0   0  -1   0   0   0       1       0;
      0   0   0   0   0   0   1   0   0   0   0  -1   0   0       0       0;
      0   0   0   0   0   0   0   1   0   0   0   0  -1   0       0       0;
      1   0   0   0   0   0   0   0   0   0   0   0   0  -1       0       0;
      1   0   0   0   0   0   0   0   0   0   0   0   0   0      -1       0;
      0   0   0   0   1   0   0   0   0   0   0   0   0   0      -1       0;
      0   0   0   0   0   0   0   0   0   1   0   0   0   0      -1       0];
       

[lambda3, x_inf3] = pageRank(B, alpha=.85)



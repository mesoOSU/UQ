function f = isposdef(M)
  
    eig_mat = eig(M);
    flag = 0;
    for i = 1:rank(M)
      if eig_mat(i) <= 0 
      flag = 1;
      end
    end
    if flag == 1
       %disp('the matrx is not positive definite')
       f = 0;
    else
       %
       %disp('the matrix is positive definite')
       f = 1;
    end
end
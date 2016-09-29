function [Q,R] = givensQR(A)
  [m,n] = size(A);
  Q = eye(m);
  R = A;
  for j = 1:n
    for i = m:-1:(j+1)
      G = eye(m);
      [c,s] = givensSinCos( R(i-1,j),R(i,j) );
      G(i-1,i-1) = c;
      G(i-1, i) = -s;
      G(i, i-1) = s;
      G(i, i) = c; 

      Q = Q*G;
    end
  end
  R = Q'*A;
end
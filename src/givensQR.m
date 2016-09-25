function [Q,R] = givensQR(A)
  [m,n] = size(A);
  Q = eye(m);
  prevG = 1;
  R = A;
  count = 0;

  for j = n-1 : 1
    for i = m:-1:(j+1)
      G = eye(m);
      if (R(i,j) == 0)
        continue;
      end
      [c,s] = givensSinCos( R(i-1,j),R(i,j) );
      count ++;
      G(i-1,i-1) = c;
      G(i-1, i) = -s;
      G(i, i-1) = s;
      G(i, i) = c;

      R(i,j)
      if (prevG != 1)
        1
        prevG = cat(1, prevG, G);
      else 
        prevG = G;
      end
      prevG
    end
  end
  Q = prevG'
  R = Q' * A;
  Q * R
  R * Q
end
function E = eigenvaluesIterative(m, L, delta1, delta2, alpha, betha)
  A = discreteJacobian(m, L, delta1, delta2, alpha, betha);
  %A = [6 5 0; 5 1 4; 0 4 3];
  %A = [0 1; -1 -1];
  n = length(A);
  i = 1;
  tol = 0.01;
  while (n > 2)
		[Q,R]=qrGS(A);
		A = R*Q;
		if ( abs(A(n,n-1)) < tol*(abs(A(n-1,n-1))+abs(A(n,n))) )  %single shift
			E(i) = A(n,n);
			i=i+1; n = n-1;	
			A = A(1:n,1:n);	
		elseif ( abs(A(n-1,n-2)) < tol*(abs(A(n-1,n-1))+abs(A(n-2,n-2))) )  %double shift
			Eaux = eig2p2 (A(n-1:n,n-1:n));
			E(i) = Eaux(1); E(i+1) = Eaux(2);
			i=i+2;n = n - 2;
			A = A(1:n,1:n);
		end
		fflush(stdout);
		n
	end
	if (n==2)
		Eaux = eig2p2 (A);
		E(i) = Eaux(1); E(i+1) = Eaux(2);
	elseif (n==1)
		E(i) = A(1,1);
	end

end;

function E = eig2p2 (A)
	p = [ 1 , -A(1,1)-A(2,2) , A(1,1)*A(2,2) - A(1,2)*A(2,1) ];
	E = roots(p);
end;
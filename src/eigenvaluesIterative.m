function E = eigenvaluesIterative(m, L, delta1, delta2, alpha, betha)
  A = discreteJacobian(m, L, delta1, delta2, alpha, betha);
  A = hessenbergTransformation(A);
  n = length(A);
  i = 1;
  tol = 0.001;
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

function A = hessenbergTransformation(A)
	n = length(A);
	for k = 1 : n - 2
	    v = A(k+1:n,k);
	    alpha = -norm(v);
	    if (v(1) < 0) alpha = -alpha; end
	    v(1) = v(1) - alpha; 
	    v = v / norm(v);
	    A(k+1:n,k+1:n) = A(k+1:n,k+1:n) - 2 * v * (v.' * A(k+1:n,k+1:n));
	    A(k+1,k) = alpha;
	    A(k+2:n,k) = 0;
	    A(1:n,k+1:n) = A(1:n,k+1:n) - 2 * (A(1:n,k+1:n) * v) * v.';
	end
end;
function e = eigenvaluesIterative(m, L, delta1, delta2, alpha, betha)
  % A = discreteJacobian(m, L, delta1, delta2, alpha, betha);
  A = [6 5 0; 5 1 4; 0 4 3];
	flag = true;
	tic()
	while (flag)
		[Q, R]= givensQR(A);
		Aaux = A;
		A = R * Q;
		Error = Aaux - A;
		errors = 0;
		for i=1 : 1 : size(A)[2];
			errors += abs(Error(i,i));
		end;
		if (errors / size(A)(1) < 1e-20)
			flag = false;
		end;
		errors / size(A)(1)
		fflush(stdout);
	end;
	toc()
	e = diag(A);
end;
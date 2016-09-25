function e = eigenvaluesIterative(m, L, delta1, delta2, alpha, betha)
	A = discreteJacobian(m, L, delta1, delta2, alpha, betha);
	flag = true;
	while (flag)
		[Q,R]= qr(A);
		Aaux = A;
		A = R * Q;
		Error = Aaux - A;
		errors = 0;
		for i=1 : 1 : size(A)[2];
			errors += abs(Error(i,i));
		end;
		if (errors / size(A)(1) < 2e-1)
			flag = false;
		end;
		errors / size(A)(1)
		fflush(stdout);
	end;
	e = A;
end;
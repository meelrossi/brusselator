function A = discreteJacobian (m, L, delta1, delta2, alpha, betha)
	% m	 		1/2 of the order of the matrix
	% L	 		bifurcation parameter (L^2 divides delta1 and delta2)
	% alpha		coefficient in reaction term for x
	% betha		coefficient in reaction term for y
	% delta1	difussion coefficient for x
	% delta2	diffusion coefficient for y
	h = 1 / (m + 1);
	square_h = h * h;
	square_L = L * L;
	square_alpha = alpha * alpha;
	I = eye(m);
	r1 = delta1 / (square_h * square_L);
	r2 = delta2 / (square_h * square_L);
	T = tridiagonal (m, 1, -2, 1);
	A = [ r1 * T + (betha - 1) * I, square_alpha * I ; -betha * I, r2 * T - square_alpha * I]; 
end;

function T = tridiagonal (m, n1, n2, n3)
	T = zeros (m, m);
	for i = 1 : 1 : m
		T (i ,i) = n2;
		if (i < m)
			T (i + 1, i) = n1;
			T (i, i + 1) = n3;
		end;
	end;
end;
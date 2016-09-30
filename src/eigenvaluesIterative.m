function E = eigenvaluesIterative(m, L, delta1, delta2, alpha, betha)
    i = 1;
    tol = 0.001;
    A = discreteJacobian(m, L, delta1, delta2, alpha, betha);
    A = hessenbergTransformation(A);
    n = length(A);
    while (n > 2)
        [Q, R] = qrGramSchmidt(A);
        A = R * Q;
        a1 = abs(A(n, n-1));
        a2 = abs(A(n-1, n-1));
        a3 = abs(A(n, n));
        a4 = abs(A(n-1, n-2));
        a5 = abs(A(n-2, n-2));
        if (tol * (a2 + a3) > a1)
            E(i) = A(n, n);
            i += 1;
            n -= 1;
            A = A(1:n, 1:n);
        elseif (tol * (a2 + a5) > a4)
            M = eigvals(A(n-1:n, n-1:n));
            E(i:i+1) = [M(1), M(2)];
            i += 2;
            n -= 2;
            A = A(1:n, 1:n);
        end
    end
    if (n == 2)
        M = eigvals(A);
        E(i:i+1) = [M(1), M(2)];
    elseif (n == 1)
        E(i) = A(1, 1);
    end
end;

function R = eigvals(M)
	R = roots([1, -M(1,1) - M(2,2), M(1,1) * M(2,2) - M(1,2) * M(2,1)]);
end;

function A = hessenbergTransformation(A)
	n = length(A);
	for x = 1:(n - 2)
	    v = A(x+1:n, x);
	    alpha = -norm(v);
	    if (v(1) < 0)
            alpha = -alpha;
        end
	    v(1) -= alpha;
	    v /= norm(v);
	    A(x+1:n, x+1:n) -= 2 * v * (v.' * A(x+1:n, x+1:n));
	    A(x+1, x) = alpha;
	    A(x+2:n, x) = 0;
	    A(1:n, x+1:n) -= 2 * (A(1:n, x+1:n) * v) * v.';
	end
end;

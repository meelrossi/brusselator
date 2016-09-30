function [rta] = start(m)

	L = 0.51302;
	delta1 = 0.008;
	delta2 = 0.004;
	alpha = 2;
	beta = 5.45;

	
	rta = double(0);

	%Para ser mas precisos hacemos corrida de 100 y promediamos
	for j=1:100
    t = time(); 
	Aeig = eigenvalues(m, L, delta1, delta2, alpha, beta);
	%Aeig = eigenvalues(m, L, delta1, delta2, alpha, beta);
	%Aeig = eigenvaluesIterative(m, L, delta1, delta2, alpha, beta);	
	rta = rta + (time()-t);	

	end

	rta = rta/100;

endfunction
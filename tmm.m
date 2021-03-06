% Main function: 
% takes the list of refractive indicies, the list of thicknesses, how many layers there are, wavelength, angle of incident, and polarity
% returns the reflection, and transfer percentages 
function [r0,t0] = main(n_list, d_list, layerCount, wvl, initial_th_i, pol, n_first, n_last)
	switch nargin
		case 4
			initial_th_i = 0;
			pol = "s";
			n_first = 1.00293;
			n_last  = 1.00293;
		case 5
			initial_th_i = initial_th_i * pi/180;
			pol = "s";
			n_first = 1.00293;
			n_last  = 1.00293;
		case 6
			initial_th_i = initial_th_i * pi/180;
			n_first = 1.00293;
			n_last  = 1.00293;
		case 7
			initial_th_i = initial_th_i * pi/180;
			n_last = 1.00293;
		case 8
			initial_th_i = initial_th_i * pi/180;
		otherwise 
			error("4 inputs are needed");
	end

	n1 = n_first;
	n2 = n_list(1);
	th_i = initial_th_i;
	th_f = snells_law(n1, n2, initial_th_i);
	[r t] = rt_layer(n1, n2, initial_th_i, th_f, pol);
	T_first = (1/t) * [ 1 r; r 1 ];

	th_f_list = snells_list(n_list, initial_th_i);
	[T_list, P_list] = matrix_list(n_list, d_list, th_f_list, wvl, pol);
	T = P_list{1};
	structureLayers = size(n_list, 1) * size(n_list,2);
	for i = 2:1:structureLayers
		T = T * T_list{i};
		T = T * P_list{i};
	end

	structCount = floor(layerCount/structureLayers);
	if (structCount > 1)
		Tstruct = T_list{1} * T;
		n = structCount - 1;
		T = T * (expRec(Tstruct, n));

		%T = T * (T_list{1} * T)^(structCount - 1);;
	end
	leftOverLayers = mod(layerCount, structureLayers);
	for i = 1:1:leftOverLayers
		T = T * T_list{i};
		T = T * P_list{i};
	end
	if (leftOverLayers == 0)
		i = structureLayers;
	else 
		i = leftOverLayers;
	end

	n1 = n_list(i);
	n2 = n_last;
	th_i = th_f_list(i);
	th_f = snells_law(n1, n2, th_i);
	[r t] = rt_layer(n1, n2, th_i, th_f, pol);
	T_last = (1/t) * [ 1 r; r 1 ];
	
	T = T_first * T * T_last;
	
	r  = T(2,1)/ T(1,1);
	t  =   1   / T(1,1);
	if (isnan(r) && ( isnan(t) || t == 0))
		r = 1;
		t = 0;
	end
	r0 = (abs(r)^2);
	t0 = (abs(t)^2) * ((real(n_last*cos(th_f)))) /...
					(real((n_first*cos(initial_th_i))));
end

% Calculates each layer's T and P matrix
function [T_list, P_list] = matrix_list(n_list, d_list, th_f_list, wvl, pol);
	n1 = n_list(size(n_list,2));
	th_i = th_f_list(size(th_f_list, 2));
	t_list = cell(1,size(n_list,2));
	p_list = cell(1,size(n_list,2));
	for i=1:1:size(n_list, 2);
		d  = d_list(i);
		n2 = n_list(i);
		th_f = th_f_list(i);
		phi = 2*pi * n2 * d * cos(th_f) / wvl;
		[r_i, t_i] = rt_layer(n1, n2, th_i, th_f, pol);
		P_i = [ exp(-1j*phi) 0; 0 exp(1j*phi) ];
		T_i = (1/t_i) * [ 1 r_i; r_i 1 ];
		t_list{i} = T_i;
		p_list{i} = P_i;
		n1 = n2;
		th_i = th_f;
	end
	T_list = t_list;
	P_list = p_list;
end

% Creates a list of all refracted angles for each layer
function th_f_list = snells_list( n_list, initial_th_i);
	n1   = 1;
	n2   = n_list(1); th_f = snells_law(n1, n2, initial_th_i); temp = [th_f]; for i=2:1:size(n_list,2)
		n1 = n2;
		n2 = n_list(i);
		th_f = snells_law(n1, n2, th_f);
		temp = [temp th_f ];
	end
	th_f_list = temp;
end

function th_f = snells_law(n1, n2, th_i);
	 th_f = asin((n1/n2)*sin(th_i));
end

% Calculates the reflection and transmission of a layer using fresnel equations
function [r,t] = rt_layer(n1, n2, th_i, th_f, pol);
	if (pol == "s") 
		[r t] = fresnel_S(n1, n2, th_i, th_f);
	elseif (pol == "p")
		[r t] = fresnel_P(n1, n2, th_i, th_f);
	end
end

function [r,t] = fresnel_S(n1, n2, th_i, th_f);
	r = (n1*cos(th_i) - n2*cos(th_f)) /...
	    (n1*cos(th_i) + n2*cos(th_f));
	t = 1 + r;
end

function [r,t] = fresnel_P(n1, n2, th_i, th_f)
	r = (n2*cos(th_i) - n1*cos(th_f)) /...
	    (n2*cos(th_i) + n1*cos(th_f));
	t = (n1/n2)*(1 + r);
end

% Recursive exponential function. This makes it run a little bit faster than normal
function result = expRec(A, n)
	result = [1 0; 0 1];
	if (n == 0)
		result = [1 0; 0 1];
		return

	elseif (n == 1)
		result = A;
		return

	elseif (mod(n, 2) == 0)
		result = expRec(A * A, n / 2);
		return

	elseif (mod(n,2) == 1) 
		result = A * expRec(A * A, (n-1) / 2);
		return
	end
	return
end

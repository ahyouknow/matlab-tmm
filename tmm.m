classdef tmm
    %	Calculates the fresnel coefficients for a multilayered material using transfer matrix method.
    %   Uses the fact that electromagnetic waves passing to and from each layer have a linear relationship
	%	More information on tmm are at http://www.photonics.intec.ugent.be/download/ocs129.pdf and https://www.osapublishing.org/ao/abstract.cfm?uri=ao-41-19-3978
	%
	%%% TODO %%%
	%
	% Impelementing incoherent and partially incoherent layers
	% More explanations and comments
	%
	%%% END %%%
    
    properties
		n_matrix   % Index of refraction matrix. Columns are for each unique material and rows are at each wavelength
		d_list     % List of layer thicknesses in meters
		num_layers % Total number of layers
		wvl_list   % List of all wavelengths to simulate in meters
		th_i_list  % List of all initial incident angles to simulate inputed as degrees
		type       % can either be S (TE) or P (TM). default is S.
    end
    
    methods
        function obj = tmm(n_matrix, d_list, num_layers, type, wvl_list, th_i_list) 
			obj.n_matrix   = n_matrix;
			obj.d_list     = d_list;
			obj.num_layers = num_layers;
			obj.wvl_list   = wvl_list;
			obj.th_i_list  = th_i_list * pi/180;
			if (type == "TM" || type == "P")
				obj.type = "P";
			else
				obj.type = "S";
			end

        end

		% Goes through each wavelength and each angle.
		% returns the fresnel coefficients in matrix form.
		% each row is for every wavelength and each 
		% column is for every initial incident angle.
		function [r_results, t_results] = run(obj)
			t_matrix = [];
			r_matrix = [];
			for wvl_index = 1:1:size(obj.wvl_list, 2)
				wvl    = obj.wvl_list(wvl_index);
				n_row  = mod(wvl_index - 1, size(obj.n_matrix, 1)) + 1;
				n_layers = obj.n_matrix(n_row, :);
				t_list = [];
				r_list = [];
				for th_i_index = 1:1:size(obj.th_i_list, 2)
					th_i = obj.th_i_list(th_i_index);

					[r t] = obj.sim_coh(n_layers, wvl, th_i); 
					t_list = [ t_list t ];
				    r_list = [ r_list r ];
				end
				t_matrix = [ t_matrix; t_list ];
				r_matrix = [ r_matrix; r_list ];
			end
			t_results = t_matrix;
			r_results = r_matrix;
		end
        
        function [r0,t0] = sim_coh(obj, n_list, wvl, initial_th_i)
			unique_layers = size(n_list, 2);

			n2   = 1.000293; % index of refraction for air
			th_f = initial_th_i;
			th_f_list = [];
			for i=1:1:unique_layers
				th_i = th_f;
				n1   = n2;
				n2   = n_list(i);
				th_f = snells_law(n1, n2, th_i);
				th_f_list = [th_f_list th_f];

			end

			n1 = n_list(unique_layers);
			th_i = th_f_list(unique_layers);
			r_list = [];
			t_list = [];
			for i=1:1:unique_layers
				n2   = n_list(i);
				th_f = th_f_list(i);

				[r t] = obj.rt_layer(n1, n2, th_i, th_f);
				r_list = [ r_list r ];
				t_list = [ t_list t ];

				n1 = n2;
				th_i = th_f;
			end
			n1 = 1.000293;
			n2 = n_list(i);
			d  = obj.d_list(1);
			th_f = th_f_list(i);
			[ r1 t1 ] = obj.rt_layer(n1, n2, initial_th_i, th_f);

			delta = (2*pi * n1 * d*cos(th_i)) /...
						          wvl;

			P_1 = [ exp(i*delta)   0;...
			          0      exp(-i*delta) ];

            T_01 = (1 / t1) * [ 1  r1;...
					     		r1  1 ];
			T_0n = T_01;

			layer_index = 1;
			for layer = 2:1:(obj.num_layers - 1)
				layer_index = mod(layer - 1, unique_layers) + 1;
				last_layer_index = mod(layer - 2, unique_layers) + 1; 

				n1 = n_list(last_layer_index);
				n2 = n_list(layer_index);
				d  = obj.d_list(layer_index);

				th_i = th_f_list(last_layer_index);
				th_f = th_f_list(layer_index);

				ri = r_list(layer_index);
				ti = t_list(layer_index);

				delta = (2*pi * n2 * d*cos(th_i)) /...
						          wvl;

				P_j = [ exp(i*delta)   0;...
				          0      exp(-i*delta) ];

            	T_ij = (1 / ti) * [ 1  ri;...
						     		ri  1 ];

				T_0n = T_0n *  P_j * T_ij;

			end
			rn = r_list(layer_index);
			tn = t_list(layer_index);

			T = T_0n;
			r = T(2,1) / T(1,1);
			t = 1 / T(1,1);

			r0 = abs(r) ^ 2;

			t0 = abs(t ^ 2) * (real(n2*cos(th_f))) /...
							(real(1.00293*cos(initial_th_i)));
        end

		function [r,t] = rt_layer(obj, n1, n2, th_i, th_f)
			if (obj.type == "S")
				r = (n1*cos(th_i) - n2*cos(th_f)) /...
				    (n1*cos(th_i) + n2*cos(th_f));
				t = 1 + r;
			else
				r = (n2*cos(th_i) - n1*cos(th_f)) /...
				    (n2*cos(th_i) + n1*cos(th_f));
				t = (n1/n2)*(1 + r);
			end
		end
        
    end
end


function delta = phase_change(n2, d, th_i, wvl)

end

function th_f = snells_law(n1, n2, th_i)
	 th_f = asin((n1/n2)*sin(th_i));
end


function DmDn = transmission_matrix(t)
	r = t - 1;

	DmDn = [ 1  r;...
		     r  1 ];
end

function DmDn = transmission_matrix_rcalc(n1, n2, th_i, th_f)
	r = r_layer(n1, n2, th_i, th_f)


	DmDn = [ 1  r; ...
		     r  1 ];
end

function p = propagation_matrix(n2, d, th_i, wvl)
	delta = phase_change(n2, d, th_i, wvl)

	p = [ exp(i*delta)     0;...
		       0     exp(-i*delta)];
end


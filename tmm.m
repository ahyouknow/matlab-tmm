classdef tmmNew
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
		n_matrix   % Index of refraction matrix. Columns are for each unique material and rows are at each wavelength.
		d_list     % List of layer thicknesses in meters.
		cell_num   % Number of cell units.
		wvl_list   % List of all wavelengths to simulate in meters.
		th_i_list  % List of all initial incident angles to simulate inputed as degrees.
		polar      % Polarization: can be S (TE) or P (TM)
    end
    
    methods
        function obj = tmmNew(n_matrix, d_list, cell_num, polar, wvl_list, th_i_list) 
			obj.n_matrix   = n_matrix;
			obj.d_list     = d_list;
			obj.cell_num   = cell_num;
			obj.wvl_list   = wvl_list;
			obj.th_i_list  = th_i_list * pi/180;

			polar = lower(polar);
			if (polar == "tm" || polar == "p")
				obj.polar = "p";
			elseif (polar == "te" || polar == "s")
				obj.polar = "s";
			end


        end

		function [r_results, t_results] = run(obj)
			% Goes through each wavelength and each angle.
			% returns the fresnel coefficients in matrix form.
			% each row is for every wavelength and each 
			% column is for every initial incident angle.
			% To make a surface plot with the output it would be "surf(th_i_list, wvl_list, r)".
			% To make a 2d plot with wavelength as the x-axis "plot(wvl_list, r(:, specific_th_i))".
			% To make a 2d plot with angle of incident as the x-axis "plot(th_i_list, r(specific_wvl, :))".
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
			t_results = t_matrix; r_results = r_matrix;
		end
        
        function [r0,t0] = sim_coh(obj, n_list, wvl, initial_th_i)
			% Main tmm function
			% tmm works by calculating a transmission matrix of one layer with fresnel's equation then
			% since the wave propagate through the material there is a propagation matrix that calculates how
			% much the wave changes within the material. Then that process is repeated for each layer in a unit cell.
			% Since the only difference between a single unit cell's transfer matrix and 2 stacked unit cells are 
			% the propagation matrix through the last layer and the transmission matrix from the last layer to the first layer.
			% So the final transfer matrix is calculated by exponenting the product of the unit cell's transmission matrix 
			% and the product of the last layer's propagation matrix and the first layer's transfer matrix. 
			% 
			% n_list is the list of index of refraction for each layer in the unit cell
			% wvl is the wavelength of the light that is propagating through the stack in meters
			% initial_th_i is the incident angle that the light first makes contact with the stack in radians
			%
			% n1 is the index of the material the light is coming from
			% n2 is the index of the material the light is hitting
			% th_i is the angle of incident in radians
			% th_f is the angle of refraction in radians
			% T_first is the transfer matrix from air to the first layer
			% T_cell is the transmission matrix of the unit cell
			% r0 is the final reflection coefficient 
			% t0 is the final transmission coefficient
			%
			n1 = 1.000293;
			n2 = n_list(1);
			th_i = initial_th_i;
			th_f = obj.snells_law(n1, n2, initial_th_i);
			[r t] = obj.rt_layer(n1, n2, initial_th_i, th_f);
			T_first = (1/t) * [ 1 r; r 1 ];
			T_cell  = [ 1 0; 0 1 ];

			for layer_index = 2:1:size(n_list, 2);

				n1  = n_list(layer_index - 1);
				n2  = n_list(layer_index);
				d   = obj.d_list(layer_index - 1);
				phi = (2*pi/wvl) * n1 * d / cos(th_f);
				P_i = [ exp(phi*j) 0; 0 exp(-1*phi*j) ];
				th_i = th_f;
				th_f = obj.snells_law(n1, n2, th_i);
				[r t] = obj.rt_layer(n1, n2, th_i, th_f);
				T_ij = (1/t) * [ 1 r; r 1 ];

				T_cell = T_cell * P_i * T_ij;
			end
			if (obj.cell_num > 1)
				n1 = n_list(size(n_list, 2));
				n2 = n_list(1);
				phi = (2*pi/wvl) * n1 * obj.d_list(size(obj.d_list, 2)) / cos(th_f);
				P_n = [ exp(phi*j) 0; 0 exp(-1*phi*j) ];
				th_i = th_f;
				th_f = obj.snells_law(n1, n2, th_i);
				[r t] = obj.rt_layer(n1, n2, th_i, th_f);
				T_n1  = (1/t) * [ 1 r; r 1 ];

				T = ( T_cell * P_n * T_n1 ) ^ (obj.cell_num - 1);
				T = T * T_cell;
			else
				T = T_cell;
			end
			T = T_first * T_cell;

			r  = T(2,1)/ T(1,1);
			t  =   1   / T(1,1);

			r0 = abs(r) ^ 2;
			t0 = (abs(t) ^ 2) * (real(n2*cos(th_f))) /...
							(real(1.00293*cos(initial_th_i)));
        end

		function th_f = snells_law(obj, n1, n2, th_i)
			 th_f = asin((n1/n2)*sin(th_i));
		end

		function [r,t] = rt_layer(obj, n1, n2, th_i, th_f)
			if (obj.polar == "s") 
				[r t] = obj.fresnel_S(n1, n2, th_i, th_f);
			elseif (obj.polar == "p")
				[r t] = obj.fresnel_P(n1, n2, th_i, th_f);
			end
		end

		function [r,t] = fresnel_S(obj, n1, n2, th_i, th_f)
			r = (n1*cos(th_i) - n2*cos(th_f)) /...
			    (n1*cos(th_i) + n2*cos(th_f));
			t = 1 + r;
		end

		function [r,t] = fresnel_P(obj, n1, n2, th_i, th_f)
			r = (n2*cos(th_i) - n1*cos(th_f)) /...
			    (n2*cos(th_i) + n1*cos(th_f));
			t = (n1/n2)*(1 + r);
		end
    end
end

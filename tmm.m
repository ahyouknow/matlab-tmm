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
		n_matrix   % Index of refraction matrix. Columns are for each unique material and rows are at each wavelength.
		d_list     % List of layer thicknesses in meters.
		layerCount % Number of Layers
		wvl_list   % List of all wavelengths to simulate in meters.
		th_i_list  % List of all initial incident angles to simulate inputed as degrees.
		pol        % Polarization: can be S (TE) or P (TM)
    end
    
    methods
        function obj = tmm(n_matrix, d_list, layerCount, pol, wvl_list, th_i_list) 
			obj.n_matrix   = n_matrix;
			obj.d_list     = d_list;
			obj.layerCount = layerCount;
			obj.wvl_list   = wvl_list;
			obj.th_i_list  = th_i_list * pi/180;
			pol = lower(pol);
			if (pol == "tm" || pol == "p")
				obj.pol = "p";
			elseif (pol == "te" || pol == "s")
				obj.pol = "s";
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
			if (size(obj.wvl_list, 1) > size(obj.wvl_list,2))
				num_wvls = size(obj.wvl_list, 1);
			else
				num_wvls = size(obj.wvl_list, 2);
			end

			if (size(obj.th_i_list, 1) > size(obj.th_i_list,2))
				num_angles = size(obj.th_i_list, 1);
			else
				num_angles = size(obj.th_i_list, 2);
			end
			t_matrix = zeros(num_wvls, num_angles);
			r_matrix = zeros(num_wvls, num_angles);
			parfor wvl_index = 1:1:num_wvls
				wvl    = obj.wvl_list(wvl_index);
				n_row  = mod(wvl_index - 1, size(obj.n_matrix, 1)) + 1;
				n_layers = obj.n_matrix(n_row, :);
				for th_i_index = 1:1:num_angles
					th_i = obj.th_i_list(th_i_index);

					[r t] = obj.sim_coh(n_layers, wvl, th_i); 
					t_matrix(wvl_index, th_i_index) = t;
					r_matrix(wvl_index, th_i_index) = r;
				end
			end
			t_results = t_matrix;
		   	r_results = r_matrix;
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
			n1 = 1.00293;
			n2 = n_list(1);
			th_i = initial_th_i;
			th_f = obj.snells_law(n1, n2, initial_th_i);
			[r t] = obj.rt_layer(n1, n2, initial_th_i, th_f);
			T_first = (1/t) * [ 1 r; r 1 ];

			th_f_list = obj.snells_list(n_list, initial_th_i);
			[T_list, P_list] = obj.matrix_list(n_list, th_f_list, wvl);
			T =  T_first;
			layer_index = 1;
			for i = 2:1:obj.layerCount
				T = T * P_list{layer_index};
				layer_index = mod(i - 1, size(obj.d_list, 2)) + 1;
				T = T * T_list{layer_index};
			end


			n1 = n_list(layer_index);
			n2 = 1.00293;
			th_i = th_f_list(layer_index);
			th_f = obj.snells_law(n1, n2, th_i);
			[r t] = obj.rt_layer(n1, n2, th_i, th_f);
			P = P_list{layer_index};
			T_last = (1/t) * [ 1 r; r 1 ];
			T = T * P * T_last;

			r  = T(2,1)/ T(1,1);
			t  =   1   / T(1,1);
			if (isnan(r) && ( isnan(t) || t == 0))
				r = 1;
				t = 0;
			end
			r0 = (abs(r)^2);
			t0 = (abs(t)^2) * ((real(1.000293*cos(th_f)))) /...
							((1.000293*cos(initial_th_i)));
        end

		function [T_list, P_list] = matrix_list(obj, n_list, th_f_list, wvl)
			n1 = n_list(size(n_list,2));
			th_i = th_f_list(size(th_f_list, 2));
			t_list = cell(1,size(n_list,2));
			p_list = cell(1,size(n_list,2));
			for i=1:1:size(n_list, 2);
				d  = obj.d_list(i);
				n2 = n_list(i);
				th_f = th_f_list(i);
				phi = 2*pi * n2 * d * cos(th_f) / wvl;
				[r_i, t_i] = obj.rt_layer(n1, n2, th_i, th_f);
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

		function th_f_list = snells_list(obj, n_list, initial_th_i)
			n1   = 1;
			n2   = n_list(1);
			th_f = obj.snells_law(n1, n2, initial_th_i);
			temp = [th_f];
			for i=2:1:size(n_list,2)
				n1 = n2;
				n2 = n_list(i);
				th_f = obj.snells_law(n1, n2, th_f);
				temp = [temp th_f ];
			end
			th_f_list = temp;
		end

		function th_f = snells_law(obj, n1, n2, th_i)
			 th_f = real(asin((n1/n2)*sin(th_i)));
		end

		function [r,t] = rt_layer(obj, n1, n2, th_i, th_f)
			if (obj.pol == "s") 
				[r t] = obj.fresnel_S(n1, n2, th_i, th_f);
			elseif (obj.pol == "p")
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

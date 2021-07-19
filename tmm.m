classdef tmm
    %	Calculates the fresnel coefficients for a multilayered material using transfer matrix method.
    %   Uses the fact that electromagnetic waves passing to and from each layer have a linear relationship
	%	More information on tmm are at http://www.photonics.intec.ugent.be/download/ocs129.pdf and https://www.osapublishing.org/ao/abstract.cfm?uri=ao-41-19-3978
	%
	%%% TODO %%%
	%
	% Picking between S and P waves
	% Impelementing incoherent and partially incoherent layers
	% Make calculations faster by creating a list of refraction angles instead of calculating it each time
	% More explanations and comments
	%
	%%% END %%%
    
    properties
		n_list     % List of refraction indexes
		d_list     % List of layer thicknesses in meters
		num_layers % Total number of layers
		wvl_list   % List of all wavelengths to simulate in meters
		th_i_list  % List of all initial incident angles to simulate inputed as degrees
    end
    
    methods
        function obj = tmm(n_list, d_list, num_layers, wvl_list, th_i_list) 
			obj.n_list	   = n_list;
			obj.d_list     = d_list;
			obj.num_layers = num_layers;
			obj.wvl_list   = wvl_list;
			obj.th_i_list  = th_i_list * pi/180;
        end

		% Goes through each wavelength and each angle.
		% returns the fresnel coefficients in matrix form.
		% each row is for every wavelength and each 
		% column is for every initial incident angle.
		function [r_results, t_results] = run(obj)
			t_matrix = [];
			r_matrix = [];
			for wvl_index = 1:1:size(obj.wvl_list, 2)
				t_list = [];
				r_list = [];
				for th_i_index = 1:1:size(obj.th_i_list, 2)
					wvl  = obj.wvl_list(wvl_index);
					th_i = obj.th_i_list(th_i_index);

					[r t] = obj.sim_coh(wvl, th_i); 
					t_list = [ t_list t ];
				    r_list = [ r_list r ];
				end
				t_matrix = [ t_matrix; t_list ];
				r_matrix = [ r_matrix; r_list ];
			end
			t_results = t_matrix;
			r_results = r_matrix;
		end
        
        function [r0,t0] = sim_coh(obj, wvl, initial_th_i)
			layer_repeat = size(obj.n_list, 2);

			n1   = 1.000293; % index of refraction for air
			n2   = obj.n_list(1);
			th_f = snells_law(n1, n2, initial_th_i);
			dis  = obj.d_list(1);
			tm   = t_layer_TE(n1, n2, initial_th_i, th_f);
            T    = (1 / tm) * transmission_matrix(tm);

			for layer = 2:1:(obj.num_layers - 1)
				layer_index = mod(layer - 1, layer_repeat) + 1;

				n1  = n2;
				n2  = obj.n_list(layer_index);


				th_i = th_f;
				th_f = snells_law(n1, n2, th_i);
				tm   = t_layer_TE(n1, n2, th_i, th_f);

				Pm   = propagation_matrix(n2, dis, th_i, wvl);
				DmDn = transmission_matrix(tm);
				T    = T * (1/tm) * Pm * DmDn;

				dis = obj.d_list(layer_index);
			end
			r = T(3) / T(1);
			t = 1 / T(1);

			r0 = abs(r) ^ 2;
			t0 = abs(t) ^ 2;

			%t0 = abs(t ^ 2) * (real(n2*cos(th_f))) /...
			%						(real(1.00293*cos(initial_th_i)));
        end
        
        function [R0,T0] = sim_wavelengths(obj, wvl_lineSpace)
            R = [];
            T = [];
            for wvl = wvl_lineSpace
                obj.wvl = wvl;
                [r t] = obj.run;
                R = [R r];
                T = [T t];
            end
            R0 = R;
            T0 = T;
        end

		function [R0,T0] = sim_aoi(obj, aoi_linespace)
			R = [];
			T = [];
			for aoi = aoi_linespace
				obj.th_i = aoi * pi / 180;
				[r t] = obj.run;
				R = [R r];
				T = [T t];
			end
			R0 = R;
			T0 = T;
		end
    end
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
	delta = phase_change(n2, d, th_i, wvl);

	p = [ exp(i*delta)     0;...
		       0     exp(-i*delta)];
end

function delta = phase_change(n2, d, th_i, wvl)

	delta = 2*pi * n2 * d*cos(th_i) /...
						wvl;
end

function t = t_layer_TE(n1, n2, th_i, th_f)
	t =      (2*n1*cos(th_i)) / ...
	    (n1*cos(th_i) + n2*cos(th_f));
end	

function r = r_layer_TE(n1, n2, th_i, th_f)
	r = (n1*cos(th_i) - n2*cos(th_f)) /...
	    (n1*cos(th_i) + n2*cos(th_f));
end

function t = t_layer_TM(n1, n2, th_i, th_f)
	t =      (2*n1*cos(th_i)) / ...
	    (n1*cos(th_f) + n2*cos(th_f));
end	

function r = r_layer_TM(n1, n2, th_i, th_f)
	r = (n2*cos(th_i) - n1*cos(th_f)) /...
	    (n2*cos(th_i) + n1*cos(th_f));
end

function th_f = snells_law(n1, n2, th_i)
	 th_f = asin((n1/n2)*sin(th_i));
end

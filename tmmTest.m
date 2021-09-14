classdef tmmTest
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
		wlLinSpace
		d_list
		n_matrix
    end
    
    methods
        function obj = tmmTest(wlLinSpace)
			cd rt_test;
			obj.wlLinSpace = wlLinSpace;
			si = matDataToN("Si_index.xlsx", wlLinSpace);
			au = matDataToN("Au.csv", wlLinSpace);
			sio2 = matDataToN("SiO2.csv", wlLinSpace);
			ti = matDataToN("Ti.csv", wlLinSpace);
			al2o3 = matDataToN("al2o3.xlsx", wlLinSpace);

			obj.n_matrix = [ti.n_matrix sio2.n_matrix au.n_matrix al2o3.n_matrix si.n_matrix ]; 
			obj.d_list = [ 10e-9 50e-9 15e-9 40e-9 100e-9 ];
			cd ..;

        end
        
        function run(obj)
			[r_30 t_30] = obj.tmm_run(0, "s");
			[r_30P t_30P] = obj.tmm_run(0, "p");
			tiledlayout(3,2);
			sgtitle("Ti / SiO2 / Au / Al2O3 / Si - 3 unit cells ");
			nexttile;
			plot(obj.wlLinSpace, r_30, obj.wlLinSpace, t_30, obj.wlLinSpace, 1 - (r_30+t_30), 'LineWidth', 2.0);
			legend("R", "T", "A");
			title("angle of incident 0 degrees S-polarized")

			nexttile;
			plot(obj.wlLinSpace, r_30P, obj.wlLinSpace, t_30P, obj.wlLinSpace, 1 - (r_30P+t_30P), 'LineWidth', 2.0);
			legend("R", "T", "A");
			title("angle of incident 0 degrees P-polarized")

			[r_30 t_30] = obj.tmm_run(30, "s");
			[r_30P t_30P] = obj.tmm_run(30, "p");
			nexttile;
			plot(obj.wlLinSpace, r_30, obj.wlLinSpace, t_30, obj.wlLinSpace, 1 - (r_30+t_30), 'LineWidth', 2.0);
			legend("R", "T", "A");
			title("angle of incident 30 degrees S-polarized")

			nexttile;
			plot(obj.wlLinSpace, r_30P, obj.wlLinSpace, t_30P, obj.wlLinSpace, 1 - (r_30P+t_30P), 'LineWidth', 2.0);
			legend("R", "T", "A");
			title("angle of incident 30 degrees P-polarized")

			[r_30 t_30] = obj.tmm_run(60, "s");
			[r_30P t_30P] = obj.tmm_run(60, "p");
			nexttile;
			plot(obj.wlLinSpace, r_30, obj.wlLinSpace, t_30, obj.wlLinSpace, 1 - (r_30+t_30), 'LineWidth', 2.0);
			legend("R", "T", "A");
			title("angle of incident 60 degrees S-polarized")

			nexttile;
			plot(obj.wlLinSpace, r_30P, obj.wlLinSpace, t_30P, obj.wlLinSpace, 1 - (r_30P+t_30P), 'LineWidth', 2.0);
			legend("R", "T", "A");
			title("angle of incident 60 degrees P-polarized")
        end

		function [r t] = tmm_run(obj, th_i, pol)
			r_1 = [];
			t_1 = [];
			for i = 1:1:(size(obj.wlLinSpace,1) * size(obj.wlLinSpace, 2))
				[r_0 t_0] = tmm(obj.n_matrix(i,:), obj.d_list, 15, obj.wlLinSpace(i), th_i, pol);
				r_1 = [r_1 r_0];
				t_1 = [t_1 t_0];
			end
			r = r_1;
			t = t_1;
		end

    end
end


classdef matDataToN
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
		n_matrix;
    end
    
    methods
        function obj = matDataToN(filename, wlLinSpace)
			startwl = wlLinSpace(1);
			endwl   = wlLinSpace(size(wlLinSpace,2)*size(wlLinSpace,1));
			T = readtable(filename, 'NumHeaderLines', 1);
			orig_wls = T{:,1};			
			if (orig_wls(1) < 1)
				orig_wls = orig_wls * 1e-6;
			else
				orig_wls = orig_wls * 1e-9;
			end
			start = find(startwl < orig_wls, 1) - 1;
			finish = find(endwl < orig_wls, 1) - 1;
			orig_n = T{:,2};
			orig_k = T{:,3};
			
			n = orig_n(start:finish);
			k = orig_k(start:finish);

			n_est = [];
			k_est = [];
			for i =1:1:((size(wlLinSpace,1)*size(wlLinSpace,2)))
				wl = wlLinSpace(i);
				x1 = find(wl >= orig_wls, 1, 'last');
				x2 = x1 + 1;
				nest = linearEstimate(orig_wls(x1), orig_wls(x2), orig_n(x1), orig_n(x2), wl);
				kest = linearEstimate(orig_wls(x1), orig_wls(x2), orig_k(x1), orig_k(x2), wl);
				n_est = [n_est; nest ];
				k_est = [k_est; kest ];
			end

			obj.n_matrix = n_est + 1j*k_est;
        end
        
    end
end

function est = linearEstimate(x1, x2, y1, y2, point)
	if (y2 == y1)
		est = y1;
		return
	end
	m = (y2-y1)/(x2-x1);
	est = m*(point - x1) + y1;
end

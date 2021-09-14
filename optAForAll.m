function [maxD l maxA] = main(wls, showgraph, targetA, startingL, estimate)
	switch nargin
		case 1
			showgraph = true;
			targetA = 0.9;
			startingL = 1;
			estimate  = true;
		case 2
			targetA = 0.9;
			startingL = 1;
			estimate  = true;
		case 3
			startingL = 1;
			estimate  = true;
		case 4
			estimate  = true;
	end
	l = startingL - 2;
	n_list = matDataToN("MoS2_monolayer_nk.xlsx", wls);
	n_matrix = [ n_list, 1.000293*ones(size(n_list, 1), 1) ];
	d_mono = 0.7e-9;
	maxA = 0;
	maxD = 0;
	maxAList = [];
	lList = [];
	while ((maxA - targetA) < 0)
		if (estimate)
			maxAList = [maxAList; maxA];
			lList = [lList; l];
			if (size(lList, 1) > 4)
				warning('');
				f = fit(maxAList, lList, 'poly2', 'Normalize', 'on', 'Robust', 'Bisquare');
				[msgstr, msgid] = lastwarn;
				if (msgid)
					estimate = false;
					break;
				else
					newL = floor(f(targetA))-2;
				end
				if (mod(newL, 2) == 0)
					newL = newL - 1;
				end
				if (newL > l)
					l = newL;
				end
				maxAList = [];
				lList = [];
			end
		end
		l = l + 2
		maxA = 0;
		maxD = 0;
		for d = 1e-9:1e-9:1000e-9
			aList = [];
			for i = 1:1:size(wls,1) * size(wls,2)
				n_list = n_matrix(i, :);
				wl = wls(i);
				[r t] = tmm(n_list, [d_mono, d], l, wl);
				a = 1-(r+t);
				aList = [aList, a];
			end
			a = mean(aList);
			if ( maxA < a )
				maxA = a;
				maxD = d;
			end
		end
	end
	if (showgraph)
		createGraph(maxD, l, d_mono, wls);
	end
	return;
end

function createGraph(maxD, l, d_mono, wls)
	tiledlayout(1,1);
	wlLinSpace = 350e-9:1e-9:1150e-9;
	n_list = matDataToN("MoS2_monolayer_nk.xlsx", wlLinSpace);
	n_matrix = [ n_list, 1.000293*ones(size(n_list, 1), 1) ];

	nexttile;
	a_list = [];
	for i = 1:1:size(wlLinSpace,1)*size(wlLinSpace, 2)
		[r t] = tmm(n_matrix(i, :), [d_mono, maxD], l, wlLinSpace(i));
		a_list = [ a_list, 1-(r+t) ];
	end
	title("MoS2 monolayer, thickness 0.7 nm with air in between");
	name = "h = " + maxD / 1e-9 + " nm, " + (l+1)/2 + " monolayers. Optimized for range " + wls(1) / 1e-9 + " - " + wls(size(wls,2)) / 1e-9 + " nm";
	plot(wlLinSpace / 1e-9, a_list, 'DisplayName', name, 'LineWidth', 2.0)
	xlabel("Wavelength (nm)");
	ylabel("Absorption");
	set(gca,'FontSize',14)
	h = legend;
end

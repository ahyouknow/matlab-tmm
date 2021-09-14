function [maxD l maxR] = main(wls, showgraph, targetR, startingL, estimate)
	switch nargin
		case 1
			showgraph = true;
			targetR = 0.9;
			startingL = 1;
			estimate  = false;
		case 2
			targetR = 0.9;
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
	maxR = 0;
	maxD = 0;
	maxRList = [];
	lList = [];
	while ((maxR - targetR) < 0)
		if (estimate)
			maxRList = [maxRList; maxR];
			lList = [lList; l];
			if (size(lList, 1) > 4)
				warning('');
				f = fit(maxRList, lList, 'poly3');
				[msgstr, msgid] = lastwarn;
				if (msgid)
					estimate = false;
					break;
				else
					l = floor(f(targetR))-2;
				end
				if (mod(l, 2) == 0)
					l = l - 1;
				end
				if (l < 0)
					l = -1*l;
				end
				maxRList = [];
				lList = [];
			end
		end
		l = l + 2
		for i = 1:1:size(wls,1) * size(wls,2)
			maxR = 0;
			maxD = 0;
			wl = wls(i);
			n_list = real(n_matrix(i, :));
			for d = 1e-9:0.1e-9:1001e-9
				[r t] = tmm(n_list, [d_mono, d], l, wl);
				if ( maxR < r )
					maxR = r;
					maxD = d;
				end
			end
		end
		maxR = maxR
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
	r_list = [];
	for i = 1:1:size(wlLinSpace,1)*size(wlLinSpace, 2)
		[r t] = tmm(n_matrix(i, :), [d_mono, maxD], l, wlLinSpace(i));
		r_list = [ r_list, r ];
	end
	title("MoS2 monolayer, thickness 0.7 nm with air in between");
	name = "h = " + maxD / 1e-9 + " nm, " + (l+1)/2 + " monolayers. Optimized for range " + wls(1) / 1e-9 + " - " + wls(size(wls,2)) / 1e-9 + " nm";
	plot(wlLinSpace / 1e-9, r_list, 'DisplayName', name, 'LineWidth', 2.0)
	xlabel("Wavelength (nm)");
	ylabel("Absorption");
	set(gca,'FontSize',14)
	h = legend;
end

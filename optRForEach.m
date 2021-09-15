function [mmaxDandL maxR] = main(wls, showgraph, targetR, startingL, estimate)
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
	n_list = matDataToN("MoS2_monolayer_nk.xlsx", wls);
	n_matrix = [ n_list, 1.000293*ones(size(n_list, 1), 1) ];
	d_mono = 0.7e-9;
	maxDandL = [];
	for i = 1:1:size(wls,1) * size(wls,2)
		lList = [];
		maxRList = [];
		maxR = 0;
		maxD = 0;
		l = startingL-2;
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
			maxR = 0;
			maxD = 0;
			wl = wls(i);
			n_list = n_matrix(i, :);
			for d = 1e-9:0.1e-9:1001e-9
				[r t] = tmm(n_list, [d_mono, d], l, wl);
				if ( maxR < r )
					maxR = r;
					maxD = d;
				end
			end
			maxR = maxR
		end
		maxDandL = [maxDandL; maxD l];
	end
	if (showgraph)
		createGraph(maxDandL, d_mono, wls);
	end
	return;
end

function createGraph(maxDandL, d_mono, wls)
	tiledlayout(1,1);
	wlLinSpace = 350e-9:1e-9:1150e-9;
	n_list = matDataToN("MoS2_monolayer_nk.xlsx", wlLinSpace);
	n_matrix = [ n_list, 1.000293*ones(size(n_list, 1), 1) ];

	nexttile;
	hold on;
	for j = 1:1:size(wls,2)
		maxD = maxDandL(j, 1);
		l = maxDandL(j, 2);
		r_list = [];
		for i = 1:1:size(wlLinSpace,1)*size(wlLinSpace, 2)
			[r t] = tmm((n_matrix(i, :)), [d_mono, maxD], l, wlLinSpace(i));
			r_list = [ r_list, r ];
		end
		title("MoS2 monolayer, thickness 0.7 nm with air in between");
		name = "h = " + maxD / 1e-9 + " nm, " + (l+1)/2 + " monolayers. Optimized for " + wls(j) / 1e-9 + " nm";
		plot(wlLinSpace / 1e-9, r_list, 'DisplayName', name, 'LineWidth', 2.0);
	end
	nmin = min(real(n_matrix(:,1)));
	nmax = max(real(n_matrix(:,1)));
	plot(wlLinSpace / 1e-9, (real(n_matrix(:, 1)) - nmin) / nmax, 'DisplayName', 'n of material', 'LineWidth', 2.0);
	kmin = min(imag(n_matrix(:,1)));
	kmax = max(imag(n_matrix(:,1)));
	plot(wlLinSpace / 1e-9, (imag(n_matrix(:, 1)) - kmin) / kmax, 'DisplayName', 'imaginary of material', 'LineWidth', 2.0);
	hold off;
	xlabel("Wavelength (nm)");
	ylabel("Resistance");
	set(gca,'FontSize',14)
	h = legend;
end

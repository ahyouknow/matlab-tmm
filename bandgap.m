function main(wl)
	maxDandL = optRForEach(wl, false, 0.9, 115, false);

	tiledlayout(2,2)
	nexttile
	wls = 325e-9:1e-9:1100e-9;
	n_matrix = matDataToN("MoS2_monolayer_nk.xlsx", wls);
	n_matrix = [ n_matrix, ones(size(n_matrix, 1), size(n_matrix, 2)) * 1.000293 ];
	wls = 325e-9:1e-9:1100e-9;
	for i=1:1:size(wl, 1) * size(wl, 2)
		l = maxDandL(i, 2);
		for d = [maxDandL(i, 1) / 2, maxDandL(i,1), maxDandL(i, 1) * 2]
			d_list = [0.7e-9, d ];
			rlist = plottmm(wls, n_matrix, d_list, l);
			bandGapSize = calcBandGap(wl, wls, rlist);
			plot(wls * 1e9, rlist, 'DisplayName', "Air thickness: " + d * 1e9 + " nm. Bandgap Size = " + bandGapSize, 'LineWidth' , 2.0)
			hold on;
		end
	end
	xlabel("Wavelength (nm)");
	ylabel("Resistance");
	legend("location", "northwest");
	hold off;
	exportResults(gca, "thickness");

	nexttile
	for i=1:1:size(wl, 1) * size(wl, 2)
		d = [0.7e-9, maxDandL(i, 1) ];
		for l = [ ceil(maxDandL(i, 2) / 2), maxDandL(i, 2), (maxDandL(i, 2) * 2) + 1]
			disp(l);
			rlist = plottmm(wls, n_matrix, d, l);
			bandGapSize = calcBandGap(wl, wls, rlist);
			plot(wls * 1e9, rlist, 'DisplayName', "monolayers: " + (l+1)/2 + ". Bandgap Size = " + bandGapSize, 'LineWidth' , 2.0)
			hold on;
		end
	end
	xlabel("Wavelength (nm)");
	ylabel("Resistance");
	legend("location", "northwest");
	hold off;
	exportResults(gca, "layers");

	nexttile
	for i=1:1:size(wl, 1) * size(wl, 2)
		l = maxDandL(i, 2);
		d = [0.7e-9, maxDandL(i, 1) ];
		for k = [1/2, 1, 2]
			n_chg = real(n_matrix) + k*imag(n_matrix);
			rlist = plottmm(wls, n_chg, d, l);
			bandGapSize = calcBandGap(wl, wls, rlist);
			plot(wls * 1e9, rlist, 'DisplayName', "k = k*" + k + ". Bandgap Size = " + bandGapSize, 'LineWidth' , 2.0)
			hold on;
		end
	end
	xlabel("Wavelength (nm)");
	ylabel("Resistance");
	legend("location", "northwest");
	hold off;
	exportResults(gca, "imaginary");

	nexttile
	for i=1:1:size(wl, 1) * size(wl, 2)
		l = maxDandL(i, 2);
		d = [0.7e-9, maxDandL(i, 1) ];
		for n = [1/2, 1, 2]
			n_chg = n*real(n_matrix) + imag(n_matrix);
			rlist = plottmm(wls, n_chg, d, l);
			bandGapSize = calcBandGap(wl, wls, rlist);
			plot(wls * 1e9, rlist, 'DisplayName', "n = n*" + n + ". Bandgap Size = " + bandGapSize, 'LineWidth' , 2.0)
			hold on;
		end
	end
	xlabel("Wavelength (nm)");
	ylabel("Resistance");
	legend("location", "northwest");
	hold off;
	exportResults(gca, "real");

	set(gca, 'FontSize', 14);
	print -dpdf -bestfit resultbandgap.pdf 
end

function rlist = plottmm(wls, n_matrix, d_list, l)
	rlist = [];
	for j=1:1:size(wls,1) * size(wls,2)
		n_list = n_matrix(j, :);
		wvl = wls(j);
		[r t] = tmm(n_list, d_list, l, wvl);
		rlist = [rlist, r];
	end
	return;
end
 
function bandgapSize = calcBandGap(wl, wls, rlist)
	index = findPeak(wl, wls, rlist);
	last = index;
	startLeft = index-1;

	while(rlist(startLeft)-rlist(last) > 0)
		last = startLeft;
		startLeft = startLeft-1;
	end

	last = index;
	startRight = index+1;

	while(rlist(startRight)-rlist(last) < 0)
		last = startRight;
		startRight = startRight+1;
	end

	bandgapSize = startRight - startLeft;
end

function peak = findPeak(wl, wls, rlist)
	index = find(wls == wl, 1);
	last = rlist(index);
	next = rlist(index+1);

	if (next - last > 0)
		step = +1;
	else
		step = -1;
	end

	last = index;
	next = index+step;
	while(rlist(next) - rlist(last) > 0)
		last = next;
		next = last+step;
	end
	peak = last;
end

function exportResults(grph, name)
	set(grph, 'FontSize', 14);
	set(grph, 'FontSize', 14);
	exportgraphics(grph, "resultsBandGap/" + name +  ".png");
end




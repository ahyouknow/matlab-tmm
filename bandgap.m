function main(wl, startingL)
	switch (nargin);
		case 0
			wl = 1000e-9;
			startingL = 1;
		case 1
			startingL = 1;
	end
	maxDandL = optRForEach(wl, false, 0.7, startingL, false);

	tiledlayout(2,2);
	nexttile
	wls = 325e-9:1e-9:1100e-9;
	n_matrix = matDataToN("MoS2_monolayer_nk.xlsx", wls);
	n_matrix = [ n_matrix, ones(size(n_matrix, 1), size(n_matrix, 2)) * 1.000293 ];
	wls = 325e-9:1e-9:1100e-9;
	hold on;
	for i=1:1:size(wl, 1) * size(wl, 2)
		l = maxDandL(i, 2);
		for d = [maxDandL(i, 1) / 2, maxDandL(i,1), maxDandL(i, 1) * 2]
			d_list = [0.7e-9, d ];
			rlist = plottmm(wls, n_matrix, d_list, l);
			bandGapSize = calcBandGap(wl, wls, rlist);
			plot(wls * 1e9, rlist, 'DisplayName', "Air thickness: " + d * 1e9 + " nm. Bandgap Size = " + bandGapSize, 'LineWidth' , 2.0)
		end
	end
	title("thickness");
	xlabel("Wavelength (nm)");
	ylabel("Resistance");
	legend("location", "northwest");
	hold off;
	exportResults(gca, "thickness");

	nexttile
	for i=1:1:size(wl, 1) * size(wl, 2)
		d = [0.7e-9, maxDandL(i, 1) ];
		for l = [ ceil(maxDandL(i, 2) / 2), maxDandL(i, 2), (maxDandL(i, 2) * 2) + 1]
			rlist = plottmm(wls, n_matrix, d, l);
			bandGapSize = calcBandGap(wl, wls, rlist);
			plot(wls * 1e9, rlist, 'DisplayName', "monolayers: " + (l+1)/2 + ". Bandgap Size = " + bandGapSize, 'LineWidth' , 2.0)
			hold on;
		end
	end
	title("layers");
	xlabel("Wavelength (nm)");
	ylabel("Resistance");
	legend("location", "northwest");
	hold off;
	exportResults(gca, "layers");

	nexttile
	for i=1:1:size(wl, 1) * size(wl, 2)
		l = maxDandL(i, 2);
		d = [0.7e-9, maxDandL(i, 1) ];
		for k = [1/100, 1, 2]
			n_chg = real(n_matrix) + k*imag(n_matrix)*j;
			rlist = plottmm(wls, n_chg, d, l);
			bandGapSize = calcBandGap(wl, wls, rlist);
			plot(wls * 1e9, rlist, 'DisplayName', "k = k*" + k + ". Bandgap Size = " + bandGapSize, 'LineWidth' , 2.0)
			hold on;
		end
	end
	title("imaginary")
	xlabel("Wavelength (nm)");
	ylabel("Resistance");
	legend("location", "northwest");
	hold off;
	exportResults(gca, "imaginary");

	nexttile
	for i=1:1:size(wl, 1) * size(wl, 2)
		l = maxDandL(i, 2);
		d = [2.4e-9, maxDandL(i, 1) ];
		rlist = plottmm(wls, n_chg, d - (2.4e-9 - 0.7e-9), 3);
		plot(wls * 1e9, rlist) 
	end
	title("0.7nm to 100 nm");
	xlabel("Wavelength (nm)");
	ylabel("Resistance");
	legend("location", "northwest");
	sgtitle("Bandgap tests");
	exportResults(gca, "real");

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
	index = find(abs(wls - wl) < 0.1e-9, 1)
	startLeft = index;

	disp(rlist(startLeft))
	while(rlist(startLeft) > 0.05)
		startLeft = startLeft-1
	end

	startRight = index;

	while(rlist(startRight) > 0.05)
		startRight = startRight+1
	end

	bandgapSize = (wls(startRight) - wls(startLeft)) * 1e9
	if (isempty(bandgapSize))
		bandgapSize = -1;
		return;
	end

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
	exportgraphics(grph, "resultsBandGap/" + name +  ".png");
end

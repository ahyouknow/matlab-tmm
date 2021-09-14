function main(wlsAvg)
	switch nargin
		case 0
			wlsAvg = 325e-9:1e-9:660e-9;
	end
	targets  = [ 0.997, 0.998, 0.99803 ];
	tCount  = size(targets, 2);
	layers = [];
	ds     = [];
	means  = [];
	l = 1;
	for i = 1:1:tCount;
		target = targets(i);
		[d l a]  = optAForAll(wlsAvg, false, target, l+2, true);
		means  = [means a];
		layers = [layers, l];
		ds     = [ds, d];
	end

	tiledlayout(1,1)
	nexttile;
	title("Approaching 1 absorption with MoS2 at 610 nm wavelength, thickness 0.7 nm, with air inbetween");
	wls = 325e-9:1e-9:1100e-9;
	n_matrix = matDataToN("MoS2_monolayer_nk.xlsx", wls);
	n_matrix = [ n_matrix, ones(size(n_matrix, 1), size(n_matrix, 2)) * 1.000293 ];
	hold on;
	for i = 1:1:tCount;
		l = layers(i);
		d = ds(i);
		a = means(i);
		aRes = [];
		for j = 1:1:size(wls,1) * size(wls,2)
			n_list = n_matrix(j, :);
			wl = wls(j);
			[r t] = tmm(n_list, [ 0.7e-9, d ], l, wl);
			aRes = [ aRes, 1-(r+t) ];
		end
		plotName = "h = " + d * 1e9 + " nm, " + (l+1)/2 + " monolayers. Mean absorption ≈ " + a;
		plot(wls * 1e9, aRes, "DisplayName", plotName, "LineWidth", 2.0);
	end
	xlabel("Wavelength (nm)");
	ylabel("Absorption");
	set(gca,'FontSize',14);
	h = legend;
	hold off;
end

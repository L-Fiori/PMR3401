function plotA(Az, title_def, M, N, dx, isE2, dt, num)

	if ~exist('dt', 'var')
		dt = 0;
	end

	if ~exist('num', 'var')
		num = 0;
	end
	
	if isE2 == false
		figure

		h = heatmap(Az, 'xlabel', 'x (cm)', 'ylabel', 'y (cm)', 'Colormap', flipud(hot));

		h.GridVisible = 'off';
		title([title_def, num2str(dx)]);

		XLabels = linspace(0, 22, N);
		CustomXLabels = string(XLabels);
		CustomXLabels(mod(XLabels, 2) ~= 0) = " ";
		h.XDisplayLabels = CustomXLabels;

		YLabels = linspace(10, -10, 2*M-1);
		CustomYLabels = string(YLabels);
		CustomYLabels(mod(YLabels, 2) ~= 0) = " ";
		h.YDisplayLabels = CustomYLabels;

	else
		figure

		h = heatmap(Az, 'xlabel', 'x (cm)', 'ylabel', 'y (cm)', 'Colormap', flipud(hot));

		h.GridVisible = 'off';
		title([title_def, num2str(dx), 'cm, instante ', num2str(num*dt), 's']);

		XLabels = linspace(0, 22, N);
		CustomXLabels = string(XLabels);
		CustomXLabels(mod(XLabels, 2) ~= 0) = " ";
		h.XDisplayLabels = CustomXLabels;

		YLabels = linspace(10, -10, 2*M-1);
		CustomYLabels = string(YLabels);
		CustomYLabels(mod(YLabels, 2) ~= 0) = " ";
		h.YDisplayLabels = CustomYLabels;
    end
        
end
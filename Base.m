function fout = Base_Schoof(x,parameters)
%elevation of base below sea level; needs to be consistent with that used
%in evolution problem for h

fout=nan.*ones(size(x));

% b_icedivide = parameters.icedivide;                  %bed elevation at ice divide
% slope = parameters.bedslope;                         %slope of bed
% sill_min = parameters.sill_min;
% sill_max = parameters.sill_max;
% slope_sill = parameters.sill_slope;

% fout = b_icedivide + (slope*x);

% fout(x < sill_min) = b_icedivide + (slope*x(x < sill_min));
% fout(x >= sill_min & x <= sill_max) = b_icedivide + (slope*sill_min) +...
%     (slope_sill*(x(x >= sill_min & x <= sill_max)-sill_min));
% fout(x > sill_max) = b_icedivide + (slope*sill_min) +...
%     (slope_sill*(sill_max-sill_min)) + (slope.*(x(x > sill_max)-sill_max));

fout = (729 - (2184.8.*(x/750e3).^2) + (1031.72.*(x/750e3).^4) -(151.72.*(x/750e3).^6));

end
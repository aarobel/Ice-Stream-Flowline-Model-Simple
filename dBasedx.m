function fout = dBasedx_Schoof(x,parameters)
%derivative of elevation of base below sea level; needs to be consistent with that used
%in evolution problem for h

fout=nan.*ones(length(x),1);

% slope = parameters.bedslope;                         %slope of bed
% sill_min = parameters.sill_min;
% sill_max = parameters.sill_max;
% slope_sill = parameters.sill_slope;
% 
% fout(x<sill_min) = slope;
% fout(x>=sill_min & x<=sill_max) = slope_sill;
% fout(x>sill_max) = slope;

fout = (-(2/750e3).*(2184.8.*(x/750e3).^1) + ((4/750e3).*1031.72.*(x/750e3).^3) -((6/750e3).*151.72.*(x/750e3).^5));
% fout = (-(2/750e3).*(parameters.bedparam1.*(x/750e3).^1) + ((4/750e3).*1031.72.*(x/750e3).^3) -((6/750e3).*151.72.*(x/750e3).^5));
end
function plotcorr(x,y,xname,yname,titlestrm)

% function plotcorr(x,y,xname,yname,titlestr)
%
% plots scatter plot of x against y and calculates r and p values
% xname,yname are optional arguments for axis labels
    
[r,p] = corrcoef(cat(2, x, y), 'Rows','pairwise');

idxy = (isnan(x)+isnan(y)) > 0;

linearCoef = polyfit(x(~idxy),y(~idxy),1);
linearFit = polyval(linearCoef,x(~idxy));
% plot(x(~idxy),y(~idxy), 'bo'); 
hold on
plot(x(~idxy),linearFit,'r--')

statstr = sprintf('r = %.3f, p = %.3f',r(2, 1),p(2, 1));
if exist('xname','var')
    xlabel({xname; statstr})
else
    xlabel(statstr)
end

if exist('titlestr','var'), title(titlestr), end
if exist('yname','var'), ylabel(yname), end


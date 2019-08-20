function [r, p, statstr, slope] = plotpartialcorr(x, y, z, xname, yname, titlestrm, linelabel, color)

% function plotcorr(x,y,xname,yname,titlestr)
%
% plots scatter plot of x against y and calculates r and p values
% xname,yname are optional arguments for axis labels
    
[r, p] = partialcorr(x, y, z, 'Rows','pairwise');

idxy = (isnan(x)+isnan(y)+isnan(z)) > 0;

linearCoef = polyfit(x(~idxy),y(~idxy),1);
linearFit = polyval(linearCoef,x(~idxy));
% plot(x(~idxy),y(~idxy), 'bo'); 
slope = linearCoef(1);
% intercept =  linearCoef(2);

if slope >= 0
    linetype = '-';
else
    linetype = ':';
end

plot(x(~idxy), linearFit, 'LineStyle', linetype, 'Color', color)

statstr = sprintf('r = %.3f, p = %.3f', r, p);

% if exist('xname','var')
%     xlabel({xname; statstr})
% else
%     xlabel(statstr)
% end
xlabel(xname)
if exist('titlestr','var'), title(titlestr), end
if exist('yname','var'), ylabel(yname), end

% xlim([min(x(~idxy))-1 max(x(~idxy))+8])
    
format short g

if exist('linelabel')

%     text(x(find(linearFit == linearFit(end)))+0.5, linearFit(find(linearFit == linearFit(end))), [linelabel ', b = ' num2str(slope) ', ' statstr], 'Color', color, 'FontSize', 8, 'HorizontalAlignment', 'left')
    text(x(find(linearFit == linearFit(end)))+.2, linearFit(find(linearFit == linearFit(end))), linelabel, 'Color', color, 'FontSize', 8, 'HorizontalAlignment', 'left')
%     text(x(find(linearFit == linearFit(end)))+5, linearFit(find(linearFit == linearFit(end))+5), statstr, 'Color', color, 'FontSize', 8, 'HorizontalAlignment', 'left')

end

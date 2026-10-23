function SimPlot(dr,figh,fs,nsi,Nth,nXstore,tlines,tplv,xpih,t0)
%SIMPLOT Draw allele histories, the dominance inset, and diversity over time.
%   xpih columns: [x, d_f, d_m, parent allele index]. Each Nth entry stores
%   [allele counts; generations]. nXstore columns are [n, ne, xr].
set(0,'CurrentFigure',figh); clf(figh);set(figh,'Color','w')

Nth_2 = Nth(~cellfun(@isempty, Nth));

tmin = t0;
tmax = t0 + max(cellfun(@(x) max(x(2,:),[],2),Nth_2));
Nmax =      max(cellfun(@(x) max(x(1,:),[],2),Nth_2));

% --- Fix y-ticks --- %
ytl = [];c=0;
while length(ytl)<3
    c=c+1;
    ytl = sort(floor(log10(tmax)*c):-1:ceil(log10(tmin)*c))/c;
    
    if c<=4
        ytl = log10(round(10.^ytl,1,'significant')); % Round to one significant digit
    end
    ytl = ytl(ytl<=log10(tmax)); % Remove ticks outside the plotted interval after rounding
end
if c==1
    % Round the final tick downward to one or two significant digits.
    ndig = 1 + (mod(log10(tmax),1)<log10(2));
    scale = 10^(ndig - 1 - floor(log10(tmax)));
    tmax_tick = floor(tmax*scale)/scale;
    ytl = unique([ytl log10(tlines(tlines>0)) log10(tmax_tick)]);
end
yti  = 10.^ytl;

if c==1
    ytlc = cell(size(ytl));
    
    ind = ~~mod(ytl,1) & 10.^mod(ytl,1)>2.9;

    ytlc( ~mod(ytl,1)) = cellfun(@(x)([                       '10^{' num2str(      x ) '}']),num2cell(ytl( ~mod(ytl,1))), 'UniformOutput', false);
    ytlc(         ind) = cellfun(@(x)([num2str(10.^mod(x,1)) '×10^{' num2str(floor(x)) '}']),num2cell(ytl(         ind)), 'UniformOutput', false);
    
    ytl = ytlc;
else
    ytl = round(10.^ytl);
end

subplot(1,2,1);
AllelePlot(xpih,Nth,dr,  yti,ytl, tmin,tmax,Nmax, t0, tlines,1/2,true,0.1,fs);

subplot(1,2,2)
nPlot(nsi,tlines,nXstore,tplv,fs,yti, tmin,tmax,t0)

end

% ============================ %
% --- Local helper functions --- %
% ============================ %

%   % ----------------------- %
%   % --- Function: nPlot --- %
%   % ----------------------- %
function nPlot(nsi,tlines,nXstore,tplv,fs ,yti, tmin,tmax,t0)
%NPLOT Plot n and ne against generation on logarithmic axes.
%   nXstore columns: [n, ne, xr]; tplv contains their sampling generations.

necol = hsl2rgb(0.08+0.25,0.6,0.50);

hold off

ti = tplv(1:nsi)+t0;

plot(    nXstore(1:nsi,2),ti,'Color',necol*0.9,'LineStyle','-' ,'LineWidth',1.8);hold on % n_e
plot(    nXstore(1:nsi,1),ti,'LineStyle','none','Marker','o','MarkerSize',2,'MarkerFaceColor','k','MarkerEdgeColor','k')  % n

ylim([tmin tmax])

xr = [0.9 max(nXstore(:,[1,2]),[],'all')*1.1];

xr(2) = max(20, xr(2)); % To standardize between plots, set n range to at least 20

xlim(xr)
xlabel('Allele count')

ax = gca;
set(ax,'YScale','log')
set(ax,'XScale','log')%,'XAxisLocation','top')
if xr(2)<20
    xticv = [1:6, 8, 12, 16];
elseif xr(2)<24
    xticv = [1:6, 8, 12, 20];
elseif xr(2)<64
    xticv = [1:4, 6, 8, 12, 20, 32, 48];
else
    xticv = [1:4, 2.^(3:10)];
end

set(ax,'XTick',xticv,'XAxisLocation','bottom','XTickLabelRotation',0, 'XMinorTick','off')

pbaspect(ax,[1 1.4 1])

set(ax,'ytick'     ,yti,'YMinorTick','on')
set(ax,'yticklabel',[])
set(ax,'FontSize',fs,'TickLength',[0.025 0.025])

for tl = tlines
    lw = 1.5;lst='--';
    plot([0.1 100], [tl tl],'Color',0.4*[1 1 1],'LineWidth',lw,'LineStyle',lst)
end

leglab = {'n_{e}','n'};

lgd = legend(ax,leglab,'FontSize',fs-1,'Location','northwest','Box','off');
lgd.ItemTokenSize = 10;
lgd.BackgroundAlpha = 0.8;

box on
% Ensure right y-axis is active and set its color to black (not the default orange)
yyaxis(ax,'right')
ax.YAxis(2).Color = 'k';
set(gca, 'YTick', [], 'YTickLabel', []);
yyaxis left

% Add xtick labels at the top as well as the bottom
xt = get(ax,'XTick');
xtl = get(ax,'XTickLabel');
% Ensure top ticks match bottom ticks and labels
set(ax,'XAxisLocation','bottom') % keep bottom
% Create a secondary (top) axes for labels without affecting the main axes
axTop = axes('Position',ax.Position, 'YTick', [], 'XTick', xt, 'XTickLabel', xtl, 'Box', 'off', 'Color','none', 'XAxisLocation','top','XScale', get(ax,'XScale'), 'FontSize', fs,'XTickLabelRotation',0, 'XMinorTick','off');pbaspect(axTop,[1 1.4 1]);axis(axis(ax))

linkaxes([ax axTop],'x');

end

%   % ------------------------- %
%   % --- Function: hsl2rgb --- %
%   % ------------------------- %
    function rgb = hsl2rgb(h,s,l)
        h = max(0,min(1,h(:))); 
        s = max(0,min(1,s(:))); 
        l = max(0,min(1,l(:))); 
        % --- Applying hue (h) --- %
        r = max(0,min(1, -1 + 6*abs(h - 12/24)));
        g = max(0,min(1,  2 - 6*abs(h -  8/24)));
        b = max(0,min(1,  2 - 6*abs(h - 16/24)));

        l = l./(1 + min([r,g,b],[],2));

        % --- Applying saturation (s) --- %
        r = 1/2 - s/2 + s.*r;
        g = 1/2 - s/2 + s.*g;
        b = 1/2 - s/2 + s.*b;

        % --- Applying luminance (l) --- %

        r(l<1/2) = r(l<1/2).*l(l<1/2)*2; r(l>=1/2) = r(l>=1/2) + (1-r(l>=1/2)).*(2*l(l>=1/2)-1);
        g(l<1/2) = g(l<1/2).*l(l<1/2)*2; g(l>=1/2) = g(l>=1/2) + (1-g(l>=1/2)).*(2*l(l>=1/2)-1);
        b(l<1/2) = b(l<1/2).*l(l<1/2)*2; b(l>=1/2) = b(l>=1/2) + (1-b(l>=1/2)).*(2*l(l>=1/2)-1);

        rgb = [r,g,b];
    end

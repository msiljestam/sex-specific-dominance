function [ax, dom_ax] = AllelePlot(Xpih, Nth,dr, yti,ytl, tmin,tmax,Nmax, t0, tlines,Np,logt,mindiff,fs)
%ALLELEPLOT Draw the allele tree plot colored by the two regulators.
%   Xpih columns: [x, d_f, d_m, parent allele index]. Nth allele histories contain
%   [allele counts; generations]. Nmax is the plot's count normalization,
%   not the census parameter in Main. dr gives the two log10 affinity ranges.

warningState = warning('off','MATLAB:polyshape:repairedBySimplify');
restoreWarning = onCleanup(@() warning(warningState)); % Restore the caller's setting

% --- Defaults & Constants --- %
if nargin < 10, tlines  = []   ; end % Vertical lines
if nargin < 11, Np      = 1/2  ; end % Branch width power, 0.5->sqrt
if nargin < 12, logt    = true ; end % Use logarithmic time
if nargin < 13, mindiff = 0.1  ; end % Minimum relative change before a new branch vertex
if nargin < 14, fs      = 7    ; end % Font size

% -- Dominance color map -- %
d_colmap = @(x,y) min(1,max(0,[1 + 0.96*y, 1 - 0.52*x - 0.48*y, 1 + 1.00*x]/2)); % cyan-orange color map

tubw = 0.12*sqrt(1/2)*sqrt(Np); % Branch width
wh = 1.4; % Height-to-width ratio

% Prepare data indices
ind     = find(cellfun(@(x) size(x,2) > 1, Nth));

% Determine ranges
x_vec = Xpih(:,1:end-1) ; % Values to plot
maxr  = diff([-0.5,0.5]); % Max range
x0 = x_vec(1,1);

if nargin<3, dr = 2*max(abs(x_vec(:,2:3)),[], 1); end
ar = max(abs(dr)); % Common log10 affinity range for both axes and colors
if ar == 0, ar = 4; end

hold on;

yr = [0 maxr*wh];
ylift = yr(2)*3e-3;
yr(2) = yr(2) + ylift;

x1r2 = [-1 1]*max(0.6,tubw*1.5 + max(abs(x_vec(:,1))));

hold on

ax = gca;
logtfunc = @(t)(log10(t)-log10(tmin))/(log10(tmax)-log10(tmin))*maxr*wh;

yti = logtfunc(yti);

hold on
set(gca,'FontSize',fs,'FontName', 'Arial')

if logt; t = cellfun(@(x)(logtfunc(x(2,:) + t0)             ),Nth(ind),'UniformOutput',false);
else   ; t = cellfun(@(x)((        x(2,:) + t0)/tmax*maxr*wh),Nth(ind),'UniformOutput',false);
end
for i = length(ind):-1:1
    N{i} = Nth{ind(i)}(1,:)/Nmax;
    temp = false;
    while sum(~temp)
        temp = [true abs(diff(N{i}.^Np)./min(N{i}(1:end-1).^Np,N{i}(2:end).^Np)) >= mindiff];
        temp(1:2:end) = true;
        temp(end) = true;
        t{i} = t{i}(temp);
        N{i} = N{i}(temp);
    end
end

cmi = zeros(length(N),3);

maxN = cellfun(@max,N);

[~,iv] = sort(maxN,'descend'); % Alleles ordered by decreasing maximum abundance

lw = 1.5;lst='--';
for tl = tlines
    if tl<=tmax
        plot([-3 3], logtfunc([tl tl]),'Color',0.4*[1 1 1],'LineWidth',lw,'LineStyle',lst)
    end
    lw = lw*0.66;
    lst = ':';
end

for i = iv
    ci1 = x_vec(ind(i),2)*2/ar;
    ci2 = x_vec(ind(i),3)*2/ar;

    cmi(i,:) = d_colmap(ci1,ci2);

    X = N{i}.^Np*tubw*maxr;
    X = [X, -fliplr(X)] + x_vec(ind(i),1);
    Y = [t{i}, fliplr(t{i})] + ylift;

    pol = polyshape(X,Y);

    plot(pol,'FaceColor',cmi(i,:),'EdgeColor',cmi(i,:)*0.85,'FaceAlpha',1,'EdgeAlpha',1);

end

ylim(yr)

xlim(x1r2)

% --- Fix y-ticks --- %
if logt
    set(ax,'ytick'     ,yti)
    if yti(2)==1 %c==1
        ytminor = log10(2:9)' + (log10(tmin):ceil(log10(tmax)));
        ytminori  = logtfunc(10.^ytminor(:));

        set(ax,'YMinorTick','on')
        ax.YAxis.MinorTickValues = ytminori;
    end

    set(ax,'yticklabel',ytl,'TickLabelInterpreter', 'tex','FontSize',fs-1)
else
    set(ax,'ytick',yti)
    set(ax,'yticklabel',ytl)
end

box on
% Ensure right y-axis is active and set its color to black (not the default orange)
yyaxis right
ax.YAxis(2).Color = 'k';
set(gca, 'YTick', [], 'YTickLabel', []);
yyaxis left

pbaspect(ax,[1 wh 1])

xt = fix(xlim*2)/2;
set(gca,'XTick',xt(1):0.5:xt(2),'TickLabelInterpreter','tex');
xtickangle(gca, 0)
xlabel('Allelic value (x)')
ylabel('Generation')

hold off
set(gca,'fontsize',fs,'fontname','arial','TickLength',[0.025 0.025])
% --- Add dominance color map --- %
    % Define normalized position [left bottom width height]

    if -0.25<x0 && x0<0
        insetPos  = [0.14  0.31  0.10  0.10];
    else
        insetPos  = [0.34  0.31  0.10  0.10];
    end

    % Create inset axes in that corner
    dom_ax = axes('Units','normalized', 'Position', insetPos, 'Box','on');

    % Draw the RGB dominance map in the inset
    n=50;
    x = linspace(-1,1,n);           % x-coordinates
    y = linspace(-1,1,n);           % y-coordinates
    [X,Y] = meshgrid(x,y);
    Z = zeros(n,n,3);
    Z(:) = d_colmap(Y(:),X(:));

    imagesc(x, y, Z, 'Parent', dom_ax);

    aTL = ["10^{-" + ar/2 + "}", "1", "10^{" + ar/2 + "}"];
    set(dom_ax, 'YDir','normal',...
        'XTick',[-1 0 1], 'XTickLabel',aTL,...
        'YTick',[-1 0 1], 'YTickLabel',aTL,...
        'FontSize',fs-3.0);
    xtickangle(dom_ax, 0)
    axis(dom_ax, 'equal', 'tight');

    % -- Plot all dominance markers -- %
    hold on
    lims = axis;
    for i = iv
        if max(t{i})==wh % Allele history reaches the final plotted generation
            scatter(dom_ax,x_vec(ind(i),3)*2/ar,x_vec(ind(i),2)*2/ar,N{i}(end)*38,'filled','MarkerFaceColor',cmi(i,:),'MarkerEdgeColor','k','LineWidth',0.5);
        end
    end
    axis(lims)

    xl = xlabel('α_♂','FontSize',fs-1.5,'FontName','DejaVu Sans');
    yl = ylabel('α_♀','FontSize',fs-1.5,'FontName','DejaVu Sans','Rotation',0);

    pos = get(xl,'Position');
    set(xl,'Position', pos + [0.5 0.3 0])

    pos = get(yl,'Position');
    set(yl,'Position', pos + [0.28 0.18 0])

    pos = get(dom_ax,'Position');
    set(dom_ax,'Position', pos + [0.024 -0.000 0 0])

end

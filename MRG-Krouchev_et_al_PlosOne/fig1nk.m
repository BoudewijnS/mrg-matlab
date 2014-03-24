% fig1nk = draws the figure.

%  =========================================================

% ---------------------------
clr = 'kr';

figNo = 1;  
figure(figNo); clf; hold on; zoom on

%  ==========================
%# Hide the y axis:
set( gca, 'YTick', [] );
set(gca, 'FontSize', 20 ); 
set(gca, 'Color', 'none', 'Box', 'off' ); 
cParentFig = get(get(gca,'parent'), 'Color' ); set(gca, 'yColor', cParentFig);

% -----------------------

% ---------------------------------------------------------

%  ==========================
% scale-bars:

Yrn = Y(:,i_nodes); 
yL = min( Yrn(:) ); yU = max( Yrn(:) ); yRange = yU - yL; 

xLim = [0, dur];
  
% xBar = 10^round( log10(dur/5) );  %ms
xBar = 1e3*dur/7; %us
xBarTag = num2str( xBar, '%3.0f [us]');

% -----------------------

yOffset = yRange/(N_nodes-1);  

yBar = 50; %% 10^round( log10(yRange) ); 

if yRange < yBar, yBar = round(yRange/2); %mV
end

yBarTag = num2str( yBar, '%3d [mV]');

yLim = [yL-N_nodes*yOffset, yU];

% -----------------------

set(gca, 'xLim', xLim, 'yLim', yLim );

%  ==========================

xx = xLim(1)+0.01*diff(xLim) + [0, 1e-3*xBar];
% xx = xLim(2) + [-1e-3*xBar 0]; 
yy = yLim(2) + [-yBar 0];

plot(xx, yy(1)+[0 0], 'k', 'lineWidth', 2 ); 
text(xx(2),yy(1),xBarTag, 'FontSize', 20, 'Horiz', 'right', 'Vertic', 'bot');

plot(xx(1)+[0 0], yy, 'k', 'lineWidth', 2 ); 
text(xx(1),yy(2),yBarTag, 'FontSize', 20, 'rotation', 90, 'Horiz', 'right', 'Vertic', 'top');

%  ==========================

for j = 1:N_nodes

if ismember(j, kStim), ic = 2; else, ic = 1; end

plot(t, Yrn(:,j) - yOffset*(j-1), clr(ic), 'lineWidth', 1 );
end %j

% txt = sprintf( 'N_{nodes} = %d', N_nodes ); title( txt );
%
xlabel( 't [ms]' ); 

%  =========================================================

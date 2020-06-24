% two cycles glued along a path
tt = 60; % number of vertices in each cycle
m = 5; % number of shared vertices
n = 2*tt - m; % total number of vertices
P = zeros(n,n); % probability transition matrix

addpath('util')

% fill first cycle
for ii = 1:(tt-1)
    P(ii,ii+1) = 1;
end

% last vertex in the first cycle splits
P(tt,1) = 0.5; 
P(tt,tt+1) = 0.5; 

% second cycle
for ii = (tt+1):(n-1)
    P(ii,ii+1) = 1;
end

% last vertex of second cycle joins back to the tt-m+1 vertex
P(n,tt-m+1) = 1; 

%P = P';

R=extEffRes(P); % effective resistance distance
R=sqrt(R);

d12 = -log(get_Ahp(P)); % beta=1/2 hitting probabilities distance
for i=1:n; d12(i,i)=0; end

d1 = -log(get_Ahp(P,1)); % beta=1 hitting probabilities distance
for i=1:n; d1(i,i)=0; end

assert(max(max((d1 - d1')./d1)) < .001); % ensure symmetry

r = round((tt-m)/2); % where to put the bidirectional edge

%plot(1:n,d12(r,:)/max(d12(r,:)),1:n,d1(r,:)/max(d1(r,:)),1:n,R(r,:)/max(R(r,:)));
%xlabel('Node');
%ylabel('Distance');
%legend('$d^{\frac12}$','$d^1$','Effective resistance metric','Location','SouthEast','Interpreter','latex');
%ylim([-.05,1.05]);

G=digraph(P); % make P a digraph for layout
ddd={d12,d1,R};
titl={'d12','d1','R'};
tiledlayout(3,1);   % If using a pre-2019 matlab, comment this out and change nexttile to figure
for i=1:3
	nexttile;
    %figure;
	p=plot(G);
	p.MarkerSize = 10;
	dd=ddd{i};
	tmp = dd(:,r);
	%tmp(r) = Inf;
	%tmp(r) = .9*min(tmp);
	p.NodeCData=tmp;
	cc=colormap(jet);
	colormap(cc(1:end-30,:));
	colorbar;
	title(titl(i));
	p.NodeLabel=[];
	highlight(p,r);
	p.ArrowPosition=.8;
	p.ArrowSize=10;
	%labelnode(p,r,num2str(r));
	
	% get color data
	c=colormap;
	cinds = p.NodeCData;
	cinds = cinds/max(cinds) * size(c,1);
	cinds = max(cinds,1);
	cinds = round(cinds);
	RGBs=c(cinds,:);

	% write coordinates
	datt=[(1:n)' p.XData' p.YData' 255*RGBs];
	datt=array2table(datt,'VariableNames',{'id','x','y','R','G','B'});
	writetable(datt,strcat(['dat/nodes' titl{i} '.dat']))
end
[u, v] = find(P);
writetable(array2table([u v],'VariableNames',{'u','v'}),'dat/edges.dat')

%xdata = zeros(1,n);
%ydata = zeros(1,n);
%for i=1:tt-m
%	xdata(i) = sin(pi*(i+m/2)/tt);
%	ydata(i) = cos(pi*(i+m/2)/tt);
%end
%for i=tt-m+1:tt
%	xdata(i) = 0;
%	ydata(i) = -1 + 2*(i-tt+m)/m;
%end
%for i=tt+1:n
%	xdata(i) = -sin(pi*(i-tt+m/2)/tt);
%	ydata(i) = cos(pi*(i-tt+m/2)/tt);
%end
%p.XData=xdata;
%p.YData=ydata;
%
%tiledlayout(3,1);
%for i=1:3
%	nexttile;
%	x=pca(ddd{i});
%	scatter(x(:,1),x(:,2),1:size(dd,1));
%	title(titl(i))
%end
%
%tiledlayout(3,1);
%for i=1:3
%	nexttile;
%	imagesc(ddd{i});
%	colorbar
%	title(titl(i));
%end

% replicate comparisons to geometric distances from the paper
clear all;
close all;
addpath util
distance_generators; % define the functions that give us distance matrices
%IND=4; % choose which distance matrix to use

%tiledlayout(2,3);
for IND=1:6
	IND
	if IND < 6 % get the adjacency matrix
		D = gd{IND};
		gamma = 1.0; f = @(D) exp(-gamma*D.^2);
		Adj=f(D);
		for ii = 1:size(Adj,1), Adj(ii,ii) = 0; end 
	else
		[D,Adj]=gd{IND}();
		Adj=full(Adj);
	end

	%figure;
	%colormap('hot')
	%imagesc(Dhp)
	%colorbar

	[Dhp,Dhpalt] = get_distances(Adj);
	%nexttile;
    	figure;
	make_figure(D,Dhp,Dhpalt)
end

function [Dhp,Dhpalt]=get_distances(Adj)

	assert(is_connected(Adj));

	d = full(sum(Adj,2));
	n=size(Adj,1);
	P = spdiags(1./d,0,n,n)*Adj;

	Dhp = get_dhp(P);
	Dhpalt = get_dhp(P,1);
end

function make_figure(D,Dhp,Dhpalt)

	global ind
	[~,ind]=sort(Dhp(1,:));
	ind = ind(2:end);
	%ind = 2:size(D,1);
	x = 2:size(D,1);

	%figure;
	scatter(x,nm(D),2,'k')
	hold on;
	scatter(x,nm(Dhp),2,'b')
	hold on;
	scatter(x,nm(Dhpalt),2,'r')
	%plot(x,nm(D),'-k',x,nm(Dhp),'--b',x,nm(Dhpalt),'-.r','LineWidth',2)
	%legend('Euclidean Distance','Hitting Prob. \beta = 1/2 Distance','Hitting Prob. \beta = 1 Distance','Location','southeast')
	set(gca,'FontSize',18,'TickLength',[.02 0])
    xlabel('Vertex Indices')
    ylabel('Distance from Vertex 1')
end

% normalization for plotting
function Y = nm(y)
	y=y(1,:);
	global ind
	Y = (y(ind) - min(y(ind)))/max(y(ind)-min(y(ind)));
end

% replicate the community detection results from the paper

Ngrid=80; % number of grid points to test in rho,delta space
nhat=3; % number of communities
reps=5; % number of times to try each method at each grid point
nmeth=3; % number of methods being tested (things break if you change this)
Nnode=300; % number of nodes in each network

addpath('util')

% setup
rho = linspace(.01,.5,Ngrid);
delta = linspace(.01,.34,Ngrid);
pu = zeros(nmeth,Ngrid,Ngrid);
directed=true;

% main for loop over rho, delta space
parfor i=1:Ngrid
disp(i)
rho_ = linspace(.01,.5,Ngrid); % linter is suggesting that sharing rho is bad
for j=1:Ngrid
	pin  = delta(i) * 2/3 + rho_(j); %in-edge density
	pout = rho_(j) - delta(i)/3; % out-edge density
	if pout >= 0 % not all delta, rho configutations are possible
		t=-1*ones(reps,nmeth); % temporarily store purity scores from the trials
		for k=1:reps
			try %#ok<TRYNC> %						
				t(k,:)=dhp_kmedoids(pin,pout,'directed',directed,'N',Nnode); % this is where all the work is done at each grid point
			end
		end
		pu(:,i,j) = max(t); % pick best of reps (5) trials
		%snr(i,j) = N*delta(i)^2/(nhat*(1-rho(j)) * rho(j));
	end
end
end
pu=max(pu,0); % negative purity is an artifact
p=p_scores(pu,Nnode,nhat,reps); % get p values compared to random clusterings


%%
% figure 1
close all;
%tiledlayout(nmeth,2); % if you get an error here, you are using an older matlab. Comment this and replace "nexttile" with "figure"

ptemp=p;
ptemp(p==min(p(:))) = .5;
for method=1:nmeth
	%nexttile;
    figure;
	plt(max(squeeze(pu(method,:,:)),.33),rho,delta,5);
	if method==1
		title('Accuracy')
	end
	colorbar;
	if method==1
		ylabel('\begin{tabular}{c} \textbf{PCA + $k$-means on $A$} \\ $\Delta$ \end{tabular}','Interpreter','latex');
	elseif method==2
		ylabel('\begin{tabular}{c} \textbf{PCA + $k$-means on $d^{\frac{1}{2}}$} \\ $\Delta$ \end{tabular}','Interpreter','latex');
	else
		ylabel('\begin{tabular}{c} \textbf{$k$-medoids $d^{\frac{1}{2}}$} \\ $\Delta$ \end{tabular}','Interpreter','latex');
	end
	%nexttile;
    figure;
	plt(squeeze(-log10(1-ptemp(method,:,:))),rho,delta,5);
	colorbar('Ticks',[1, 2, 3],'TickLabels',{'.1','.01','.001'});
	if method==1
		title('p value');
		%title('\textbf{Negative $\log$ p-value}','Interpreter','latex');
	end
end
set(gcf,'Position',[0 0 1000 1000]);

% figure 2
[res,resp,best]=best_method_matrix(pu,p,rho,delta(1:24));

% %plt(log(snr),rho,delta,'log detectability');
%subplot(2,nmeth,nmeth)
%plt(squeeze(pu(1,:,:)-pu(2,:,:)),rho,delta,'Rand index');
%subplot(2,nmeth,6)
%plt(squeeze(pu(method,:,:)),rho,delta,'Rand index');

function p=p_scores(pu,N,nhat,nrep)
	% pu -- purity scores
	% N -- node count
	% nhat -- community count
	% nrep -- how many trials you take the best of
	pus = random_pu_maxes(N,nhat,nrep);
	pus=sort(pus);
	p=zeros(size(pu));
	pusmax=max(pus);
	for i=1:numel(pu)
		if pu(i) >= pusmax
			p(i)=1;
		else
			[~,ind]=max(pu(i)<pus);
			p(i) = ind/numel(pus);
		end
	end
end

function pus = random_pu_maxes(N,nhat,nrep)
	% get max purity after nrep random clusterings
	Npart = 4000;

	gt = ones(N,1);
	n=N/nhat;
	for a=2:nhat
		gt( (1+(a-1)*n) : (a*n) ) = a;
	end

	parts = randi(nhat,Npart*nrep,N);
	pus = zeros(Npart,1);
	for i=1:Npart
		t=zeros(nrep,1);
		for r=1:nrep
			t(r)=purityMeas(parts(nrep*(i-1)+r,:),gt);
		end
		pus(i)=max(t);
	end
end

function [residualAcc,residualP,bestInd] = best_method_matrix(A,p,rho,delta)

	% Compute best
	[bestVal,bestInd] = max(A);
	bestVal=squeeze(bestVal);
	bestInd = squeeze(bestInd);

	% First subplot
	figure;
	%tiledlayout(1,3);
	%nexttile
	plt(bestInd .* squeeze(max(A) > .33),rho,delta,2);
	title("Best method")

	% Complute second best
	B=A; N=size(A,2); for i=1:N
		for j=1:N
			B(bestInd(i,j),i,j)=-Inf;
		end
	end
	[secondVal,secondInd] = max(B);
	secondVal=squeeze(secondVal);
	secondInd=squeeze(secondInd);

	% Second subplot
	residualAcc = bestVal-secondVal;
	%nexttile;
    figure;
	plt(residualAcc,rho,delta,2);
	title("Accuracy difference")
	colorbar('Ticks',[0, .1,.2,.3,.4])

	% Second best p
	N = size(A,2);
	residualP = zeros(N);
	for i=1:N
		for j=1:N
			if max(A(:,i,j)) > 0
				residualP(i,j) = -log10(1-p(bestInd(i,j),i,j)) + log10(1-p(secondInd(i,j),i,j));
			end
		end
	end
	min(residualP(:))

	% Third subplot
	%nexttile;
    figure;
	plt(residualP,rho,delta,2);
	title("p value ratio")
	colorbar('Ticks',[0, 1, 2],'TickLabels',{'1','10','100'})

	set(gcf,'Position',[100 100 1500 500])
end

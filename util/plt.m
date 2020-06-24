function plt(p,rho,delta,ytix)

	N_rho = length(rho);
	N_delta = length(delta);
	imagesc(p(1:N_delta,1:N_rho));
	xtklbl = cell(N_rho);
	ytklbl = cell(N_delta);
	for k1 = 1:length(rho)
		xtklbl{k1} = sprintf('%d', round(100*rho(k1)));
		%xtklbl{k1} = sprintf('%.2f', rho(k1));
	end
	for k1 = 1:length(delta)
		ytklbl{k1} = sprintf('%d', round(100*delta(k1)));
		%ytklbl{k1} = sprintf('%.2f', delta(k1));
	end
	set(gca,'XTick',1:5:N_rho,'XTickLabel',xtklbl(1:5:N_rho));
	set(gca,'YTick',1:ytix:N_delta,'YTickLabel',ytklbl(1:ytix:N_delta));
	xlabel('$100\rho$','Interpreter','latex');
	ylabel('$100\Delta$','Interpreter','latex');

end

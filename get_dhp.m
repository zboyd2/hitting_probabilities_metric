% compute the hitting probability distance
% P is th etransition probability from i to j
function dhp=get_dhp(P,beta)

	addpath('util')

	if nargin < 2
		beta=.5;
	end

	dhp = -log10(get_Ahp(P,beta));
	for i=1:size(dhp,1); dhp(i,i)=0; end;

end

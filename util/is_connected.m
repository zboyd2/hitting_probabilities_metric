function ret = is_connected(A)
	ret = (sum(get_one_component(A),1) == size(A,1)) ...
	   && (sum(get_one_component(A'),1)== size(A,1));
end

% two cycles glued along a path
tt = 7; % number of vertices in each cycle
m = 3; % number of shared vertices
n = 2*tt - m; % total number of vertices
P = zeros(n,n);

addpath('util')

% fill first cycle
for ii = 1:(tt-1)
    P(ii,ii+1) = 1;
end

% last vertex in the first cycle splits
P(tt,1) = 0.5; 
P(tt,tt+1) = 0.5; 

% last vertex of second cycle joins back to the tt-m+1 vertex
P(n,tt-m+1) = 1; 

% second cycle
for ii = (tt+1):(n-1)
    P(ii,ii+1) = 1;
end

% add an equatorial link
i1=ceil( (tt-m)/2);
i2=n - ceil((tt-m)/2);
P(i1,i2)=.5;
P(i2,i1)=.5;
P(i1,i1+1)=.5;
P(i2,i2+1)=.5;

assert(all(abs(sum(P,2)-1) < .001))

% get fiedler vectors and eigenvectors three ways
[Fiedler_vec, lam] = naive_fiedler(P);
[Fiedler_vecFC, lamFC] = fiedler_FC(P);
[Fiedler_vecHT,lamHT] = ht_spectral(P); % hitting probabilities

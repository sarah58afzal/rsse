function trell = poly2bcjrtrellis(poly, xin)

% number of input symbols
[N, M] = size(xin);

% constraint length
k = length(poly)-1;

% number of states
trell.S = M^k;

% states identifier in decimal notation
next_states_dec = zeros(trell.S, M);

outputs_dec = zeros(N, trell.S, M);

for i=1:trell.S, % for all current states
    for j=1:M, % for all possible inputs
        current_mary = de2ma(i-1, M, k);
        current_ip = j-1;
        next_states_mary = [current_ip current_mary(1:end-1)];  
        next_states_dec(i, j) = ma2de(next_states_mary, M);
        outputs_dec(:, i, j) = xin(:, [current_ip current_mary]+1)*poly';
    end
end
trell.nstates = next_states_dec;
trell.outputs = outputs_dec;
trell.ninputs = M;
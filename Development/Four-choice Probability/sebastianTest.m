k = 25;
nLoops = 10000;
permArray = size([65,1]);
simData = size([nLoops,1]);
probabilities = size([nLoops,1]);
for i = 1:nLoops
    for j = 1:65
       permArray(j) = ceil(rand(1)*4);
    end
    simData(i) = (sum(permArray == 3)>=k);
end

p = sum(simData) / nLoops);
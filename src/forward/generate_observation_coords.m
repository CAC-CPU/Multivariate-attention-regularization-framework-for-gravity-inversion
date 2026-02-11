function generate_observation_coords
dataNum = 41;
xObs = linspace(0,2000,dataNum)';
zObs = ones(dataNum,1).*(-50);
forwardObs = [xObs,zObs];
save 'data/forward_obs.txt' 'forwardObs' -ascii;
end
distanceMatrices = cell(1);
distanceMatrices{1} = magic(10);

filtrationValues = cell(1);
filtrationValues{1} = 1:5;

maxDimension = 2;

computePersistenceModules(distanceMatrices, filtrationValues, maxDimension)
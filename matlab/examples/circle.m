
distanceMatrices = cell(1);
filtrationValues = cell(1);

%points = util.Point.getRandomSpherePoints(20, 1);
points = util.Point.circle2D(1, 10);

coords = zeros(points.size, 2);

for i=0:points.size()-1
    p = points.get(i).getX();
    for j=0:p.size()-1
        coords(i+1, j+1) = p.get(j);
    end
end

plot(coords(:, 1), coords(:, 2), 'ro')

%%

distanceMatrices{1} = matrix.distancematrix.DistanceMatrix.computeEuclideanDistanceMatrix(...
    points);
distanceMatrices{2} = matrix.distancematrix.DistanceMatrix.computeKNNMatrix(...
    distanceMatrices{1});

filtrationValues{1} = [0 0.1 2];
filtrationValues{2} = [0 20 100];
    
cpers = computePersistenceModules(distanceMatrices, filtrationValues, 2)



%%

persistenceModule = cpers.get(1);

fcf = standardFCF(persistenceModule)

plotFCF(fcf)


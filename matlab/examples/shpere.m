
distanceMatrices = cell(1);
filtrationValues = cell(1);

s_dim = 2;

points = util.Point.getRandomSpherePoints(50, s_dim);

coords = zeros(points.size, s_dim+1);

for i=0:points.size()-1
    p = points.get(i).getX();
    for j=0:p.size()-1
        coords(i+1, j+1) = p.get(j);
    end
end

plot3(coords(:, 1), coords(:, 2), coords(:, 3), 'ro')

%%

distanceMatrices{1} = matrix.distancematrix.DistanceMatrix.computeEuclideanDistanceMatrix(...
    points);
distanceMatrices{2} = matrix.distancematrix.DistanceMatrix.computeKNNMatrix(...
    distanceMatrices{1});

filtrationValues{1} = [0 0.2 2];
filtrationValues{2} = [0 20 100];

spherePersistenceModules = computePersistenceModules(distanceMatrices, filtrationValues, 3);

%%

persistenceModule = spherePersistenceModules.get(2);

fcf = standardFCF(persistenceModule)

plotFCF(fcf)


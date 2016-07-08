
distanceMatrices = cell(1);
filtrationValues = cell(1);

points = util.Point.circle2D(1, 10);

coords = zeros(points.size, 2);

for i=0:points.size()-1
    p = points.get(i).getX();
    for j=0:p.size()-1
        coords(i+1, j+1) = p.get(j);
    end
end

plot(coords(:, 1), coords(:, 2), 'ro')

distanceMatrices{1} = matrix.distancematrix.DistanceMatrix.computeEuclideanDistanceMatrix(...
    points);

filtrationValues{1} = 0:0.1:2;

circlePersistenceModules = computePersistenceModules(distanceMatrices, filtrationValues, 2)

persistenceModule = circlePersistenceModules.get(1);
direction = [0];

fcf = domainFCF(persistenceModule, direction)

plotFCF(fcf)


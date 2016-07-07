
distanceMatrices = cell(1);
filtrationValues = cell(1);

s_dim = 2;

points = util.Point.getRandomSpherePoints(100, s_dim);

coords = zeros(points.size, s_dim+1);

for i=0:points.size()-1
    p = points.get(i).getX();
    for j=0:p.size()-1
        coords(i+1, j+1) = p.get(j);
    end
end

plot3(coords(:, 1), coords(:, 2), coords(:, 3), 'ro')

distanceMatrices{1} = matrix.distancematrix.DistanceMatrix.computeEuclideanDistanceMatrix(...
    points);

filtrationValues{1} = 0:0.1:1;

spherePersistenceModules = computePersistenceModules(distanceMatrices, filtrationValues, 2)
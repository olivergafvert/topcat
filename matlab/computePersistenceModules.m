% Computes the persistence modules of the filtered space described by the
% 'distanceMatrices' and 'filtrationVlues'. A multifiltered Vietoris-Rips
% complex is contructed and the homology of this is computed.
function persistenceModules = getPersistenceModules(distanceMatrices, filtrationValues, maxDimension)

dataSystem = topcat.util.DataSystem.create();

for i=1:length(distanceMatrices)
    dataSystem.addDistanceMatrix(distanceMatrices{i});
end

for i=1:length(filtrationValues)
    dataSystem.addFiltrationValues(filtrationValues{i});
end

persistenceModules = topcat.persistence.PersistenceModuleCollection.create(dataSystem, ...
                                                                    maxDimension);

% Computes the persistence modules of the filtered space described by the
% 'distanceMatrices' and 'filtrationVlues'. A multifiltered Vietoris-Rips
% complex is contructed and the homology of this is computed.
function persistenceModules = computePersistenceModules(distanceMatrices, filtrationValues, maxDimension)

dataSystem = topcat.util.DataSystem.create();

for i=1:length(distanceMatrices)
    dataSystem.addDistanceMatrix(distanceMatrices{i});
end

for i=1:length(filtrationValues)
    if length(filtrationValues{i}) ~= 3
        error('filtrationValues{i} must have length equal to 3. Format: [start step_size end]');
    end
    f_arr = filtrationValues{i};
    dataSystem.addFiltrationValues(f_arr(1):f_arr(2):f_arr(3));
end

persistenceModules = topcat.persistence.PersistenceModuleCollection.create(dataSystem, ...
                                                                    maxDimension);

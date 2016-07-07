function f = standardFCF(persistenceModule, direction)

fcf = topcat.persistence.noise.StandardNoise.computeFCF(...
    persistenceModule.getFunctor(), persistenceModule.getFiltrationValues(), ...
    direction);

%Turn fcf into matlab array

function f = domainFCF(persistenceModule, direction)

fcf = topcat.persistence.noise.DomainNoise.computeFCF(...
    persistenceModule.getFunctor(), persistenceModule.getFiltrationValues(), ...
    direction);

%Turn fcf into matlab array

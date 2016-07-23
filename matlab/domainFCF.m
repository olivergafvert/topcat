function f = domainFCF(persistenceModule, direction)

fcf = topcat.persistence.noise.DomainNoise.computeFCFStatic(...
    persistenceModule.getFunctor(), persistenceModule.getFiltrationValues());

f = zeros(fcf.size(), 2);

f(:, 1) = fcf.getEpsilons();
f(:, 2) = fcf.getValues();


function f = standardFCF(persistenceModule, direction)

fcf = topcat.persistence.noise.StandardNoise.computeFCF(...
    persistenceModule.getFunctor(), persistenceModule.getFiltrationValues(), ...
    direction);

f = zeros(fcf.size(), 2);

f(:, 1) = fcf.getEpsilons();
f(:, 2) = fcf.getValues();

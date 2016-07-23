function r = rankInv(persistenceModule, u, v)
%RANKINV
%   Computes the rank invarant between the positions u and v.
%   'persistenceModule' - an object of the class PersistenceModule
%   'u' - an array with integer values describing the starting position.
%   'v' - an array with integer values describing the end position.

F = persistenceModule.getFunctor();
start_pos = topcat.util.IntTuple(u);
end_pos = topcat.util.IntTuple(v);
r = topcat.matrix.BMatrix.rank(F.getMap(start_pos, end_pos));

end


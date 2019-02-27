/*
Topcat - a multidimensional persistent homology library.
Copyright (C) 2016 Oliver Gäfvert

This file is part of Topcat.

Topcat is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Topcat is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

package topcat.persistence;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import topcat.matrix.exception.NoSolutionException;
import topcat.persistence.functor.Functor;
import topcat.persistence.functor.exception.MalformedFunctorException;
import topcat.persistence.homology.HomologyUtil;
import topcat.persistence.simplex.SimplexStorageStructure;
import topcat.persistence.simplex.SimplicialComplex;
import topcat.matrix.distancematrix.DistanceMatrix;
import topcat.util.IntTuple;

import java.util.ArrayList;
import java.util.List;

/**
 * Represents a collection of persistence modules over different dimensions.
 */
public class PersistenceModuleCollection extends ArrayList<PersistenceModule>{
    private static final Logger log = LoggerFactory.getLogger(PersistenceModuleCollection.class);

    private PersistenceModuleCollection(){
        super();
    }

    public int getMaxDimension(){
        return get(size()-1).getDimension();
    }

    public static PersistenceModuleCollection create(List<DistanceMatrix> distanceMatrices,
                                                     List<List<Double>> filtrationValues, int maxDimension){
        SimplexStorageStructure simplexStorageStructure = SimplicialComplex.computeSimplexStream(distanceMatrices, filtrationValues, maxDimension);
        return create(simplexStorageStructure, filtrationValues, maxDimension);
    }

    public static PersistenceModuleCollection create(SimplexStorageStructure simplexStorageStructure,
                                                     List<List<Double>> filtrationValues, int maxDimension){
        PersistenceModuleCollection persistenceModuleCollection = new PersistenceModuleCollection();
        try {
            IntTuple size = IntTuple.zeros(filtrationValues.size());
            for(int i=0;i<filtrationValues.size();i++){
                size.set(i, filtrationValues.get(i).size()-1);
            }
            List<Functor> functors = HomologyUtil.computeHomologyFunctors(simplexStorageStructure, size, maxDimension);
            for(int i=0; i<functors.size();i++) {
                PersistenceModule persistenceModule = new PersistenceModule(functors.get(i), i, filtrationValues);
                persistenceModuleCollection.add(persistenceModule);
            }
        }catch (MalformedFunctorException | NoSolutionException ex){
            log.error("Failed to create persistence module.", ex);
        }
        return persistenceModuleCollection;
    }
}

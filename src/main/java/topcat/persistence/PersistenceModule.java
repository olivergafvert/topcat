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

import topcat.persistence.functor.Functor;
import topcat.util.Pair;

import java.util.List;

/**
 * Represents a persistence module, i.e a functor F: Q^r -> Vect_K. It is implemented as a functor F: N^r -> Vect_K
 * where we have a map \alpha : N^r -> Q^r mapping the filtrationvalues (see [1] for details).
 *
 * [1] - Multidimensional Persistence and Noise by Chachólski et al. (arXiv:1505.06929).
 */
public class PersistenceModule {
    protected int dimension;
    protected Functor F;

    //Tracks the indexing of the filtration values in each dimension
    protected List<List<Double>> filtrationIndices;

    PersistenceModule(Functor F, int dimension, List<List<Double>> filtrationIndices){
        this.F = F;
        this.dimension = dimension;
        this.filtrationIndices = filtrationIndices;
    }

    public Functor getFunctor(){
        return F;
    }

    public int getDimension(){
        return dimension;
    }

    public List<List<Double>> getFiltrationValues(){
        return filtrationIndices;
    }

}

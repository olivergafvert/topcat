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

package topcat.matrix.rankminimization;

import topcat.matrix.BMatrix;
import topcat.matrix.BVector;
import topcat.matrix.exception.AffineVectorSpaceDimensionException;
import topcat.matrix.exception.NoSolutionException;
import topcat.matrix.exception.WrongDimensionException;

/**
 * Same as VectorSpaceIterator but for an affine subspace.
 */
public class AffineVectorSpaceIterator extends VectorSpaceIterator {

    protected BVector v;

    protected AffineVectorSpaceIterator(BMatrix ker, BVector v, int cache_size) {
        super(ker, cache_size);
        this.v = v;
    }

    /**
     * Constructor for AffineVectorSpaceIterator. Iterates over all elements
     * in the subspace given by v + <basis>.
     * @param basis
     * @param v
     * @param cache_size
     * @return
     * @throws WrongDimensionException
     */
    public static AffineVectorSpaceIterator create(BMatrix basis, BVector v, int cache_size) throws AffineVectorSpaceDimensionException {
        if(basis.rows > 64) {
            throw new AffineVectorSpaceDimensionException("Does not support vector spaces of dim > 64");
        }
        return new AffineVectorSpaceIterator(basis, v, cache_size);
    }

    @Override
    public BVector next() {
        BVector n = null;
        try{
            n = super.next().plus(v);
        }catch (WrongDimensionException wde){
            wde.printStackTrace();
        }
        return n;
    }
}

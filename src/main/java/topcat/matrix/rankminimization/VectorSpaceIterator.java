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

import gnu.trove.map.hash.TIntObjectHashMap;
import topcat.matrix.BMatrix;
import topcat.matrix.BVector;
import topcat.matrix.exception.WrongDimensionException;

import java.util.Iterator;

/**
 * Iterates over all vectors in a subspace of a vector space over the field Z/2Z given a basis.
 */
public class VectorSpaceIterator implements Iterator<BVector>{

    //A matrix where the rows are a basis for the vector space
    protected BMatrix basis;
    private long max_value;
    private long pos = 0;
    private int cache_size = 0;
    private TIntObjectHashMap<BVector> cache = new TIntObjectHashMap<>();
    private BVector prev;

    protected VectorSpaceIterator(BMatrix basis, int cache_size){
        this.basis = basis;
        this.max_value = 1L << basis.rows;
        this.cache_size = cache_size;
    }

    public static VectorSpaceIterator create(BMatrix basis, int cache_size) throws WrongDimensionException{
        if(basis.rows > 64)
            throw new WrongDimensionException("Does not support vector spaces of dim > 64");
        return new VectorSpaceIterator(basis, cache_size);
    }

    public int getDimension(){
        return basis.rows;
    }

    public int getAmbientDimension(){
        return basis.cols;
    }

    public void reset(){
        prev = null;
        pos = 0;
    }

    @Override
    public boolean hasNext() {
        return pos < max_value;
    }

    @Override
    public BVector next(){
        if(prev==null){
            prev = new BVector(basis.cols);
            pos++;
            return prev.copy();
        }
        BVector v = prev;
        try {
            long diff = pos^(pos-1); //Calc bit diff
            if(diff < this.cache_size){
                if(cache.contains((int)diff)){
                    v = prev.plus(cache.get((int)diff));
                }else{
                    BVector w = new BVector(basis.cols);
                    for (int i = 0; i < Integer.highestOneBit((int) diff); i++) {
                        if (diff >> i == 1) {
                            w = w.plus(basis.getRow(i));
                        }
                    }
                    cache.put((int) diff, w);
                    v = v.plus(w);
                }
            }else {
                v = new BVector(basis.cols);
                for (int i = 0; i < basis.rows; i++) {
                    if (pos >> i == 1) {
                        v = v.plus(basis.getRow(i));
                    }
                }
            }
            pos++;
        }catch (WrongDimensionException wde){
            wde.printStackTrace();
            pos = max_value;
            return null;
        }
        prev = v;
        return v;
    }
}

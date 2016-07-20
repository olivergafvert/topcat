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

package topcat.persistence.functor;

import topcat.matrix.BMatrix;
import topcat.matrix.exception.WrongDimensionException;
import topcat.persistence.functor.exception.MalformedNaturalTransformation;
import topcat.util.Grid;
import topcat.util.GridIterator;
import topcat.util.IntTuple;

import java.util.List;

/**
 * Represents a natural transformation between two functors.
 */
public class Nat {
    private IntTuple size;
    Grid<BMatrix> maps;

    public Nat(IntTuple size){
        this.size = size;
        maps = Grid.create(size);
    }

    public void setMap(IntTuple v, BMatrix A){
        maps.set(v, A);
    }

    public BMatrix getMap(IntTuple v){
        return maps.get(v);
    }

    public IntTuple getSize(){
        return size;
    }

    /**
     * Verify that the natural transformation nat: F -> G respects commutativity.
     * @param F
     * @param G
     * @param nat
     * @return
     */
    public static void verify(Functor F, Functor G, Nat nat) throws WrongDimensionException, MalformedNaturalTransformation{
        for(IntTuple v : GridIterator.getSequence(G.getSize())){
            List<IntTuple> basis = IntTuple.getStandardBasisSequence(v.length());
            for (int i = 0; i < basis.size(); i++) {
                IntTuple w = v.plus(basis.get(i));
                if(w.leq(nat.size)) {
                    assert (nat.getMap(w).mult(F.getMap(v, w)).equals(G.getMap(v, w).mult(nat.getMap(v))));
                }
            }
        }
    }

    @Override
    public String toString(){
        StringBuilder sb = new StringBuilder();
        for(IntTuple v : GridIterator.getSequence(size)){
            sb.append(v).append(":\n");
            sb.append(getMap(v));
        }
        return sb.toString();
    }

}

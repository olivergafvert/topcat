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

package topcat.util;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * An iterator over a n-dimensional grid.
 */
public class GridIterator implements Iterator<IntTuple> {
    private IntTuple bounds;
    private IntTuple current;

    public GridIterator(IntTuple bounds){
        this.bounds = bounds;
        List<Integer> zeros = new ArrayList<>();
        for(int i=0;i<bounds.length();i++){
            zeros.add(0);
        }
        this.current = new IntTuple(zeros);
    }

    @Override
    public boolean hasNext() {
        return current.get(0).compareTo(bounds.get(0)) <= 0;
    }

    @Override
    public IntTuple next() {
        //Copy current
        IntTuple ret = new IntTuple(current);
        increment();
        return ret;
    }

    @Override
    public void remove() {
        //TODO: Implement.
    }

    private void increment(){
        List<Integer> elements = current.getElements();
        int k=elements.size()-1;
        while(elements.get(k).equals(bounds.get(k)) && k>0){
            elements.set(k--, 0);
        }
        elements.set(k, elements.get(k)+1);
    }

    /**
     * Returns a sequence ordered by the leftmost coordinate with non-zero difference.
     * I.e (1, 2) < (2, 1) and (1, 1) < (1, 2).
     * @param bound
     * @return
     */
    public static List<IntTuple> getSequence(IntTuple bound){
        GridIterator iterator = new GridIterator(bound);
        List<IntTuple> sequence = new ArrayList<>();
        while(iterator.hasNext()) sequence.add(iterator.next());
        return sequence;
    }

    /**
     * Returns a path in the grid from position 'from' to position 'to'.
     * @param from
     * @param to
     * @return
     */
    public static List<IntTuple> getSequence(IntTuple from, IntTuple to){
        List<IntTuple> sequence = new ArrayList<>();
        IntTuple current = new IntTuple(from);
        sequence.add(current);
        for(int i=0;i<current.length();i++){
            for(int j=0;j<to.get(i)-current.get(i);j++){
                IntTuple w = IntTuple.getStandardBasisElement(current.length(), i);
                current = current.plus(w);
                sequence.add(new IntTuple(current));
            }
        }
        return sequence;
    }

    /**
     * Returns a path in the grid from position 'from' to position 'to'.
     * @param from
     * @param to
     * @return
     */
    public static List<Integer> getIntegerSequence(IntTuple from, IntTuple to){
        List<Integer> sequence = new ArrayList<>();
        IntTuple current = new IntTuple(from);
        for(int i=0;i<current.length();i++){
            for(int j=0;j<to.get(i)-current.get(i);j++){
                IntTuple w = IntTuple.getStandardBasisElement(current.length(), i);
                current = current.plus(w);
                sequence.add(i);
            }
        }
        return sequence;
    }

    /**
     * Returns a path in the grid from position 'from' to position 'to'.
     * @param from
     * @param to
     * @return
     */
    public static List<Pair<IntTuple, Integer>> getPairSequence(IntTuple from, IntTuple to){
        List<Pair<IntTuple, Integer>> sequence = new ArrayList<>();
        IntTuple current = new IntTuple(from);
        for(int i=0;i<current.length();i++){
            int diff = to.get(i)-current.get(i);
            for (int j = 0; j < diff; j++) {
                IntTuple w = IntTuple.getStandardBasisElement(current.length(), i);
                sequence.add(new Pair<>(current, i));
                current = current.plus(w);
            }
        }

        return sequence;
    }
}

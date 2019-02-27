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
import java.util.Arrays;
import java.util.List;

/**
 * Represents an element in the poset Z^r.
 */
public class IntTuple extends Tuple<Integer>{

    public IntTuple(List<Integer> tuple){
        super(tuple);
    }

    public IntTuple(IntTuple intTuple){
        super(new ArrayList<>(intTuple.getElements()));
    }

    public IntTuple(Integer... val){
        super(Arrays.asList(val));
    }

    /**
     * Take double array as input and return the closest IntTuple.
     * @param vals
     */
    public IntTuple(double[] vals){
        super(new ArrayList<Integer>());
        for(int i=0;i<vals.length;i++){
            this.tuple.add((int)Math.round(vals[i]));
        }
    }

    /**
     * Computes the element-wise addition of this tuple with 'v'.
     * @param v
     * @return
     */
    public IntTuple plus(IntTuple v){
        IntTuple w = IntTuple.zeros(v.length());
        for(int i=0;i<v.length();i++){
            w.set(i, v.get(i)+get(i));
        }
        return w;
    }

    /**
     * Computes the element-wise subtraction of this tuple with 'v'.
     * @param v
     * @return
     */
    public IntTuple minus(IntTuple v){
        IntTuple w = IntTuple.zeros(v.length());
        for(int i=0;i<v.length();i++){
            w.set(i, get(i)-v.get(i));
        }
        return w;
    }

    @Override
    public boolean equals(Object o){
        if(o==null || !o.getClass().equals(this.getClass())) return false;
        IntTuple intTuple = (IntTuple) o;
        if(this.length() != intTuple.length()) return false;
        for(int i=0;i<length();i++){
            if(!get(i).equals(intTuple.get(i))) return false;
        }
        return true;
    }

    /**
     * Returns true if this tuple is less than 'v' in the poset.
     * @param v
     * @return
     */
    public boolean lt(IntTuple v){
        for(int i=0;i<length();i++) {
            if (get(i) >= v.get(i)) return false;
        }
        return true;
    }

    /**
     * Returns true if this tuple is less than or equal to 'v' in the poset.
     * @param v
     * @return
     */
    public boolean leq(IntTuple v){
        for(int i=0;i<length();i++) {
            if (get(i) > v.get(i)) return false;
        }
        return true;
    }

    /**
     * Returns the minimal value of the entries in the tuple.
     * @return
     */
    public Integer min(){
        if(length() == 0) return 0;
        int min = get(0);
        for(int i=1;i<length();i++) {
            if (get(i) < min) min = get(i);
        }
        return min;
    }

    /**
     * Returns the maximal value of the entries in the tuple.
     * @return
     */
    public Integer max(){
        if(length() == 0) return 0;
        int max = get(0);
        for(int i=1;i<length();i++) {
            if (get(i) > max) max = get(i);
        }
        return max;
    }

    /**
     * Computes the poset join of v and w.
     * @param v
     * @param w
     * @return
     */
    public static IntTuple join(IntTuple v, IntTuple w){
        IntTuple z = IntTuple.zeros(v.length());
        for(int i=0;i<z.length();i++) z.set(i, v.get(i)<w.get(i) ? w.get(i) : v.get(i));
        return z;
    }

    /**
     * Computes the poset meet of v and w.
     * @param v
     * @param w
     * @return
     */
    public static IntTuple meet(IntTuple v, IntTuple w){
        IntTuple z = IntTuple.zeros(v.length());
        for(int i=0;i<z.length();i++) z.set(i, v.get(i)<w.get(i) ? v.get(i) : w.get(i));
        return z;
    }

    /**
     * Returns true if 'v' has any negative elements.
     * @param v
     * @return
     */
    public boolean hasNegativeElements(IntTuple v){
        for(int i=0;i<v.length();i++)
            if(v.get(i)<0) return true;
        return false;
    }

    /**
     * Returns a tuple length 'length' with a 1 in position 'i' and zeros everywhere else.
     * @param length
     * @param i
     * @return
     */
    public static IntTuple getStandardBasisElement(int length, int i){
        IntTuple e = IntTuple.zeros(length);
        e.set(i, 1);
        return e;
    }

    /**
     * Returns a sequence of tuples that corresponds to the standard basis in N^d.
     * @param d - the dimension of the ambient space
     * @return
     */
    public static List<IntTuple> getStandardBasisSequence(int d){
        List<IntTuple> sequence = new ArrayList<>();
        for(int i=0;i<d;i++){
            IntTuple tuple = IntTuple.zeros(d);
            tuple.set(i, 1);
            sequence.add(tuple);
        }
        return sequence;
    }

    /**
     * Returns an tuple of length 'length' containing 1's.
     * @param length
     * @return
     */
    public static IntTuple ones(int length){
        return uniform(length, 1);
    }

    /**
     * Returns an tuple of length 'length' containing 0's.
     * @param length
     * @return
     */
    public static IntTuple zeros(int length){
        return uniform(length, 0);
    }

    /**
     * Returns an tuple of length 'length' with 'value' in each position.
     * @param length
     * @return
     */
    public static IntTuple uniform(int length, int value){
        List<Integer> e = new ArrayList<>();
        for(int i=0;i< length;i++){
            e.add(value);
        }
        return new IntTuple(e);
    }


    private static void choose_helper(IntTuple v, int level, int index, List<IntTuple> vs){
        if(level == v.length()-index){
            for(int i=index;i<v.length();i++){
                v.set(i, 1);
            }
            vs.add(new IntTuple(v));
        }
        else if(level==0){
            for(int i=index;i<v.length();i++){
                v.set(i, 0);
            }
            vs.add(new IntTuple(v));
        }
        else {
            v.set(index, 1);
            choose_helper(v, level-1, index+1, vs);
            v.set(index, 0);
            choose_helper(v, level, index+1, vs);
        }
    }

    /**
     * Returns a list of binary vectors with all combinations of tuples with k ones.
     * @param n
     * @param k
     * @return
     */
    public static List<IntTuple> choose(int n, int k){
        List<IntTuple> vs = new ArrayList<>();
        choose_helper(IntTuple.zeros(n), k, 0, vs);
        return vs;
    }
}

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

package topcat.matrix;

import gnu.trove.iterator.TIntIterator;
import gnu.trove.set.hash.TIntHashSet;
import topcat.matrix.exception.WrongDimensionException;

/**
 * Represents a vector with coefficients in the field Z/2Z. It is implemented as a boolean matrix with element-wise
 * XOR as additive operation and the AND operation for matrix multiplication.
 */
public class BVector {
    private TIntHashSet pos; //the non-zero positions
    private int length;

    public BVector(int length){
        this(length, new TIntHashSet());
    }

    public BVector(int length, int[] nonzeroIndices){
        this(length);
        for(int i : nonzeroIndices) pos.add(i);
    }

    public BVector(int length, TIntHashSet pos){
        this.length = length;
        this.pos = pos;
    }

    public BVector(BVector V){
        this.length = V.length;
        this.pos = new TIntHashSet();
        TIntIterator iterator = V.getIndexSetIterator();
        while(iterator.hasNext()){
            set(iterator.next(), true);
        }
    }

    public BVector copy(){
        return new BVector(this);
    }

    public int getNumberOfNonZeroElements(){
        return pos.size();
    }

    public boolean get(int i){
        return pos.contains(i);
    }

    public void set(int i, boolean val){
        if(val){
            pos.add(i);
        }else if(pos.contains(i)){
                pos.remove(i);
        }
    }

    public int getLength(){
        return length;
    }

    /**
     * Computes the result of this + other binary vector
     * @param other
     * @return
     */
    public BVector plus(BVector other) throws WrongDimensionException {
        if(other.length != this.length){
            throw new WrongDimensionException("Length of vectors don't match.");
        }
        TIntHashSet larger;
        TIntHashSet smaller;

        if(pos.size() < other.pos.size()){
            larger = other.pos;
            smaller = pos;
        }else{
            larger = pos;
            smaller = other.pos;
        }
        TIntIterator iterator = smaller.iterator();
        TIntHashSet n = new TIntHashSet(larger);
        while(iterator.hasNext()){
            int val = iterator.next();
            if(n.contains(val)){
                n.remove(val); //XOR operation
            }else{
                n.add(val);
            }
        }
        return new BVector(length, n);
    }

    public TIntIterator getIndexSetIterator(){
        return pos.iterator();
    }

    public static BVector concat(BVector v1, BVector v2){
        BVector v = new BVector(v1.length+v2.length);
        TIntIterator iterator = v1.getIndexSetIterator();
        while(iterator.hasNext()){
            v.set(iterator.next(), true);
        }
        iterator = v2.getIndexSetIterator();
        while(iterator.hasNext()){
            v.set(iterator.next()+v1.length, true);
        }
        return v;
    }

    public BVector subvector(int start, int end){
        BVector v = new BVector(end-start);
        TIntIterator iterator = this.getIndexSetIterator();
        while(iterator.hasNext()){
            int n = iterator.next();
            if(n >= start && n < end){
                v.set(n, true);
            }
        }
        return v;
    }

    @Override
    public boolean equals(Object o){
        if(o == null || !o.getClass().equals(this.getClass())) return false;
        BVector v = (BVector) o;
        if(length != v.length) return false;
        if(getNumberOfNonZeroElements() != v.getNumberOfNonZeroElements()) return false;
        TIntIterator iterator = getIndexSetIterator();
        while(iterator.hasNext()){
            if(!v.get(iterator.next())) return false;
        }
        return true;
    }

    @Override
    public String toString(){
        StringBuilder sb = new StringBuilder();
        sb.append("[ ");
        if(pos.contains(0)) sb.append(1);
        else sb.append(0);
        for(int i=1;i<length;i++){
            sb.append(", ");
            if(pos.contains(i)) sb.append(1);
            else sb.append(0);
        }
        sb.append(" ]");
        return sb.toString();
    }
}

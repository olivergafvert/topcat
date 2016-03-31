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


/**
 * Represents a 2-dimensional integer tuple. I.e an element in Z^2.
 */
public class IntPair extends Pair<Integer, Integer> implements Comparable<IntPair>{

    public IntPair(Integer i1, Integer i2) {
        super(i1, i2);
    }

    public IntPair minus(IntPair a){
        return new IntPair(this._1()-a._1(), this._2()-a._2());
    }

    public IntPair plus(IntPair a){
        return new IntPair(this._1()+a._1(), this._2()+a._2());
    }

    @Override
    public boolean equals(Object o){
        if(o==null || !o.getClass().equals(IntPair.class)) return false;
        IntPair ip = (IntPair) o;
        return this._1()==ip._1() && this._2()==ip._2();
    }

    public int compareTo(IntPair o) {
        int res = this._1().compareTo(o._1());
        if(res != 0) return res;
        return this._2().compareTo(o._2());
    }

    @Override
    public int hashCode(){
        return this._1()<<16 + this._2();
    }

    /**
     * Computes the poset meet of two elements in N^2.
     * @param i1
     * @param i2
     * @return
     */
    public static IntPair meet(IntPair i1, IntPair i2){
        int a = i1._1() < i2._1() ? i1._1() : i2._1();
        int b = i1._2() < i2._2() ? i1._2() : i2._2();
        return new IntPair(a, b);
    }

    /**
     * Computes the poset join of two elements in N^2.
     * @param i1
     * @param i2
     * @return
     */
    public static IntPair join(IntPair i1, IntPair i2){
        int a = i1._1() > i2._1() ? i1._1() : i2._1();
        int b = i1._2() > i2._2() ? i1._2() : i2._2();
        return new IntPair(a, b);
    }
}

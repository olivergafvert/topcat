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
 * Represents a pair of objects.
 */
public class Pair<T, U> {
    private final T i1;
    private final U i2;

    public Pair(T i1, U i2){
        this.i1 = i1;
        this.i2 = i2;
    }

    public T _1(){
        return i1;
    }

    public U _2(){
        return i2;
    }

    @Override
    public boolean equals(Object o){
        if(o==null || !o.getClass().equals(getClass())){
            return false;
        }
        Pair<T, U> other = (Pair<T, U>) o;
        if(_1().equals(other._1()) && _2().equals(other._2())){
            return true;
        }
        return false;
    }

    @Override
    public String toString(){
        return (new StringBuilder())
                .append("[ ").append(this._1()).append(", ").append(this._2()).append(" ]")
                .toString();
    }

}

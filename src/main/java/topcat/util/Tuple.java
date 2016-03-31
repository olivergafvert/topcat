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

import java.util.Arrays;
import java.util.List;

/**
 * Represents an n-dimensional tuple for some n>0.
 */
public class Tuple<T> {
    protected List<T> tuple;

    public Tuple(List<T> tuple){
        this.tuple = tuple;
    }

    public Tuple(T... val){
        this(Arrays.asList(val));
    }

    public List<T> getElements(){
        return tuple;
    }

    public T get(int i){
        return tuple.get(i);
    }

    public void set(int i, T val){
        tuple.set(i, val);
    }

    public int length(){ return tuple.size(); }

    public List<T> toList(){
        return tuple;
    }

    @Override
    public int hashCode(){
        return tuple.hashCode();
    }

    @Override
    public String toString(){
        if(length() == 0){
            return "[]";
        }
        StringBuilder sb = new StringBuilder();
        sb.append("[ ").append(get(0));
        for(int i=1;i<length();i++){
            sb.append(", ").append(get(i));
        }
        sb.append(" ]");
        return sb.toString();
    }
}

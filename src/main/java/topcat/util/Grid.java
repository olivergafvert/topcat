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

import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;

/**
 * Represents an n-dimensional grid for some n>0. It is implemented using
 * a hierarchy of hashmaps.
 */
public class Grid<T> {

    private IntTuple size;
    private SubGrid grid = new SubGrid(0);

    public Grid(IntTuple size){
        this.size = size;
    }

    public static Grid create(IntTuple size){
        if(size.length() > 1){
            return new Grid<>(size);
        }else{
            return new Grid1D<>(size);
        }
    }

    public T get(IntTuple v){
        return grid.get(v);
    }

    public void set(IntTuple v, T t){
        grid.set(v, t);
    }

    public IntTuple size(){
        return size;
    }

    @Override
    public String toString(){
        StringBuilder sb = new StringBuilder();
        for(IntTuple v : GridIterator.getSequence(size)){
            sb.append(v).append(":\n");
            sb.append(get(v)).append("\n");
        }
        return sb.toString();
    }

    class SubGrid{
        protected int level;
        private Int2ObjectOpenHashMap<SubGrid> subgrid = new Int2ObjectOpenHashMap<>();

        SubGrid(int level){
            this.level = level;
        }

        T get(IntTuple v){
            if(v.length()-level > 2){
                if(!subgrid.containsKey(v.get(level))){
                    return null;
                }
                return subgrid.get(v.get(level)).get(v);
            }
            if(!subgrid.containsKey(v.get(level))){
                return null;
            }
            return ((BaseGrid)subgrid.get(v.get(level))).get(v.get(level+1));
        }

        void set(IntTuple v, T t){
            if(v.length()-level > 2) {
                if (!subgrid.containsKey(v.get(level))){
                    subgrid.put(v.get(level), new SubGrid(level+1));
                }
                subgrid.get(v.get(level)).set(v, t);
            }else{
                if (!subgrid.containsKey(v.get(level))){
                    subgrid.put(v.get(level), new BaseGrid(level+1));
                }
                BaseGrid baseGrid = (BaseGrid) subgrid.get(v.get(level));
                baseGrid.set(v.get(level+1), t);
            }
        }
    }

    class BaseGrid extends SubGrid{
        private Int2ObjectOpenHashMap<T> elements = new Int2ObjectOpenHashMap<>();

        BaseGrid(int level) {
            super(level);
        }

        T get(int i){
            if(!elements.containsKey(i)){
                return null;
            }
            return elements.get(i);
        }

        void set(int i, T t){
            elements.put(i, t);
        }
    }
}

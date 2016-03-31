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
 * Represents a 1-dimensional grid, i.e a list.
 */
public class Grid1D<T> extends Grid<T> {
    private BaseGrid baseGrid = new BaseGrid(0);

    Grid1D(IntTuple size) {
        super(size);
    }

    @Override
    public T get(IntTuple v) {
        return baseGrid.get(v.get(0));
    }

    @Override
    public void set(IntTuple v, T val){
        baseGrid.set(v.get(0), val);
    }
}

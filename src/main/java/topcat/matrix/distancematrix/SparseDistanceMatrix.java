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

package topcat.matrix.distancematrix;

import it.unimi.dsi.fastutil.ints.Int2DoubleOpenHashMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;

import java.util.List;

/**
 * Created by oliver on 2016-04-05.
 */
public class SparseDistanceMatrix extends DistanceMatrix {
    private Int2ObjectOpenHashMap<Int2DoubleOpenHashMap> valueMap = new Int2ObjectOpenHashMap<>();

    public SparseDistanceMatrix(int rows, int cols) {
        super(rows, cols);
    }

    @Override
    public double get(int i, int j){
        if(valueMap.containsKey(i) && valueMap.get(i).containsKey(j)){
            return valueMap.get(i).get(j);
        }else{
            return Double.MAX_VALUE;
        }
    }

    @Override
    public void set(int i, int j, double val){
        if(!valueMap.containsKey(i)){
            valueMap.put(i, new Int2DoubleOpenHashMap());
        }
        valueMap.get(i).put(j, val);
    }

    @Override
    public void setRow(int i, List<Double> list){
        for(int j=0;j<list.size();j++){
            set(i, j, list.get(j));
        }
    }

    @Override
    public double[] getRow(int i){
        double[] row = new double[cols];
        for(int j=0;j<cols;j++){
            row[j] = get(i, j);
        }
        return row;
    }

    @Override
    public int[] getNonZeroRows() {
        return valueMap.keySet().toIntArray();
    }

    @Override
    public int[] getNonZeroRowEntries(int i) {
        if(valueMap.containsKey(i)){
            return valueMap.get(i).keySet().toIntArray();
        }
        return new int[0];
    }


}

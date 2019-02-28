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

import java.util.Arrays;
import java.util.List;

/**
 * Created by oliver on 2016-04-05.
 */
public class ArrayDistanceMatrix extends DistanceMatrix {
    private double[][] distanceMatrix;

    public ArrayDistanceMatrix(double[][] distanceMatrix){
        super(distanceMatrix.length, distanceMatrix.length > 0 ? distanceMatrix[0].length : 0);
        this.distanceMatrix = distanceMatrix;
    }

    public ArrayDistanceMatrix(List<List<Double>> distanceMatrix){
        this(distanceMatrix.size(), distanceMatrix.size() > 0 ? distanceMatrix.get(0).size() : 0);
        for(int i=0;i<distanceMatrix.size();i++){
            for(int j=0;j<distanceMatrix.get(i).size();j++){
                set(i, j, distanceMatrix.get(i).get(j));
            }
        }
    }

    public ArrayDistanceMatrix(int rows, int cols){
        super(rows, cols);
        this.distanceMatrix = new double[rows][];
        for(int i=0;i<rows;i++){
            this.distanceMatrix[i] = new double[cols];
        }
    }

    @Override
    public double get(int i, int j){
        return distanceMatrix[i][j];
    }

    @Override
    public void set(int i, int j, double val){
        distanceMatrix[i][j] = val;
    }

    @Override
    public void setRow(int i, List<Double> list){
        for(int j=0;j<list.size();j++){
            distanceMatrix[i][j] = list.get(j);
        }
    }

    @Override
    public double[] getRow(int i){
        return distanceMatrix[i];
    }

    @Override
    public int[] getNonZeroRows() {
        return range(0, rows);
    }

    @Override
    public int[] getNonZeroRowEntries(int i) {
        return range(0, cols);
    }

    private static int[] range(int start, int end){
        int length = end-start;
        int[] arr = new int[length];
        for(int i=0;i<length;i++){
            arr[i] = start+i;
        }
        return arr;
    }

    @Override
    public String toString(){
        StringBuilder sb = new StringBuilder();
        for(int i=0;i<rows;i++){
            sb.append(Arrays.toString(getRow(i)));
        }
        return sb.toString();
    }
}

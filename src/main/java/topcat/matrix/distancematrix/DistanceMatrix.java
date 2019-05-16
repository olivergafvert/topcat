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

import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import topcat.util.Pair;
import topcat.util.Point;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

public abstract class DistanceMatrix {
    public final int cols;
    public final int rows;

    public DistanceMatrix(int rows, int cols){
        this.rows = rows;
        this.cols = cols;
    }

    public abstract double get(int i, int j);

    public abstract void set(int i, int j, double val);

    public abstract void setRow(int i, List<Double> list);

    public abstract double[] getRow(int i);

    /**
     * Returns the indices of the rows which have set values.
     * @return
     */
    public abstract int[] getNonZeroRows();

    /**
     * Returns the indices of row 'i' which have set values.
     * @param i
     * @return
     */
    public abstract int[] getNonZeroRowEntries(int i);

    /**
     * Returns a collection of hashsets where the ith hashset contains the vertices
     * with ordering larger than i such that the value of (i, j) is less than or equal
     * to 'val'.
     * @param val
     * @return
     */
    public Int2ObjectOpenHashMap<IntOpenHashSet> getEdgesLEQThan(double val){
        Int2ObjectOpenHashMap<IntOpenHashSet> edges = new Int2ObjectOpenHashMap<>();
        int[] nonzero_rows = getNonZeroRows();
        for(int i=0;i<nonzero_rows.length;i++){
            int row = nonzero_rows[i];
            IntOpenHashSet vertices = new IntOpenHashSet();
            int[] nonzero_cols = getNonZeroRowEntries(row);
            for(int j : nonzero_cols){
                if(j > row && get(row, j) <= val){
                    vertices.add(j);
                }
            }
            edges.put(row, vertices);
        }
        return edges;
    }

    public static DistanceMatrix computeAxisDistanceMatrix(List<Point> points, int axis){
        DistanceMatrix distanceMatrix = new ArrayDistanceMatrix(points.size(), points.size());
        for(int i=0;i<points.size();i++){
            distanceMatrix.set(i, i, 0);
            for(int j=i+1;j<points.size();j++){
                double d = (points.get(i).getX().get(axis) - points.get(j).getX().get(axis))*(points.get(i).getX().get(axis) - points.get(j).getX().get(axis));
                distanceMatrix.set(i, j, d);
                distanceMatrix.set(j, i, d);
            }
        }
        return distanceMatrix;
    }

    public static DistanceMatrix computeEuclideanDistanceMatrix(List<Point> points){
        DistanceMatrix distanceMatrix = new ArrayDistanceMatrix(points.size(), points.size());
        for(int i=0;i<points.size();i++){
            distanceMatrix.set(i, i, 0);
            for(int j=i+1;j<points.size();j++){
                double d = euclideanDistance(points.get(i).getX(), points.get(j).getX());
                distanceMatrix.set(i, j, d);
                distanceMatrix.set(j, i, d);
            }
        }
        return distanceMatrix;
    }

    public static double euclideanDistance(List<Double> a, List<Double> b){
        double d = 0;
        for(int i=0;i<a.size();i++){
            d += (a.get(i)-b.get(i))*(a.get(i)-b.get(i));
        }
        return Math.sqrt(d);
    }

    /**
     * Computes the p-norm of the elements in the list 'a'.
     * @param a
     * @param p
     * @return
     */
    public static double norm(List<Double> a, double p){
        double n = 0;
        for(int i=0;i<a.size();i++) n += Math.pow(Math.abs(a.get(i)), p);
        return Math.pow(n, 1/p);
    }


    public static DistanceMatrix codensityMatrix(DistanceMatrix distanceMatrix){
        DistanceMatrix densityMatrix = new ArrayDistanceMatrix(distanceMatrix.rows, distanceMatrix.cols);
        for(int r=0; r<distanceMatrix.rows;r++){
            double[] row = distanceMatrix.getRow(r);

            //Elements of the row with positions
            List<Pair<Double, Integer>> elements = new ArrayList<>();
            for(int i=0;i<row.length;i++){
                elements.add(new Pair<>(row[i], i));
            }
            Collections.sort(elements, new Comparator<Pair<Double, Integer>>() {
                @Override
                public int compare(Pair<Double, Integer> o1, Pair<Double, Integer> o2) {
                    return o1._1().compareTo(o2._1());
                }
            });

            densityMatrix.set(r, r, -densityMatrix.rows);
            //Set the nearest neighbors
            for(int i=1;i<elements.size();i++){
                if(elements.get(i)._2() < r){
                    if(densityMatrix.get(elements.get(i)._2(), r) > -i/elements.get(i)._1()){
                        densityMatrix.set(r, elements.get(i)._2(), densityMatrix.get(elements.get(i)._2(), r));
                    }else{
                        densityMatrix.set(r, elements.get(i)._2(), -i/elements.get(i)._1());
                    }
                }else{
                    densityMatrix.set(r, elements.get(i)._2(), -i/elements.get(i)._1());
                }
            }
        }
        return densityMatrix;
    }

    public static DistanceMatrix computeKNNMatrix(DistanceMatrix distanceMatrix){
        DistanceMatrix densityMatrix = new ArrayDistanceMatrix(distanceMatrix.rows, distanceMatrix.cols);
        for(int r=0; r<distanceMatrix.rows;r++){
            double[] row = distanceMatrix.getRow(r);

            //Elements of the row with positions
            List<Pair<Double, Integer>> elements = new ArrayList<>();
            for(int i=0;i<row.length;i++){
                elements.add(new Pair<>(row[i], i));
            }
            Collections.sort(elements, new Comparator<Pair<Double, Integer>>() {
                @Override
                public int compare(Pair<Double, Integer> o1, Pair<Double, Integer> o2) {
                    return o1._1().compareTo(o2._1());
                }
            });

            //Set the nearest neighbors
            for(int i=0;i<elements.size();i++){
                densityMatrix.set(r, elements.get(i)._2(), i);
            }
        }
        return densityMatrix;
    }

    public List<List<Double>> toList(){
        List<List<Double>> dmat = new ArrayList<>();
        for(int i=0;i<this.rows;i++){
            List<Double> row = new ArrayList<>();
            for(int j=0;j<this.cols;j++){
                row.add(get(i, j));
            }
            dmat.add(row);
        }
        return dmat;
    }
}

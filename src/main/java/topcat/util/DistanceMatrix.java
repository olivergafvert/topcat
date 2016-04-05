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

import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.hash.TIntHashSet;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.function.BiFunction;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

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
    public TIntObjectHashMap<TIntHashSet> getEdgesLEQThan(double val){
        TIntObjectHashMap<TIntHashSet> edges = new TIntObjectHashMap<>();
        int[] nonzero_rows = getNonZeroRows();
        IntStream.range(0, nonzero_rows.length).parallel().mapToObj(i -> {
            int row = nonzero_rows[i];
            TIntHashSet vertices = new TIntHashSet();
            int[] nonzero_cols = getNonZeroRowEntries(row);
            for(int j : nonzero_cols){
                if(j > row && get(row, j) <= val){
                    vertices.add(j);
                }
            }
            return new Pair<>(row, vertices);
        }).collect(Collectors.toList())
                .stream()
                .forEach(pair -> edges.put(pair._1(), pair._2()));
        return edges;
    }

    public static DistanceMatrix computeDistanceMatrix(List<Point> points, BiFunction<Point, Point, Double> function){
        List<List<Double>> distanceMatrix = points.parallelStream()
                .map(
                        p1 -> points.stream()
                                .map(
                                        p2 -> function.apply(p1, p2)
                                ).collect(Collectors.toList())
                ).collect(Collectors.toList());
        DistanceMatrix d = new ArrayDistanceMatrix(points.size(), points.size());
        IntStream.range(0, distanceMatrix.size()).forEach(i -> d.setRow(i, distanceMatrix.get(i)));
        return d;
    }

    public static DistanceMatrix computeDensityMatrix(DistanceMatrix distanceMatrix){
        List<List<Double>> densityMatrix = IntStream.range(0, distanceMatrix.rows).parallel().mapToObj(r -> {
            double[] row = distanceMatrix.getRow(r);
            List<Pair<Double, Integer>> elements = new ArrayList<>();
            for(int i=0;i<row.length;i++){
                elements.add(new Pair<>(row[i], i));
            }
            Collections.sort(elements, (o1, o2) -> o1._1()<o2._1() ? -1 : 1);
            List<Pair<Double, Integer>> density = new ArrayList<>();
            for(int i=0;i<elements.size();i++){
                double radius = elements.get(i)._1();
                double V = 4*3.14*radius*radius*radius/3;
                if(radius < 1E-12){
                    density.add(new Pair<>(i+1.0, elements.get(i)._2()));
                }else{
                    density.add(new Pair<>((i+1.0)/V, elements.get(i)._2()));
                }
            }
            Collections.sort(density, (o1, o2) -> o1._2()-o2._2());
            return density.stream().map(e -> e._1()).collect(Collectors.toList());
        }).collect(Collectors.toList());
        DistanceMatrix d = new ArrayDistanceMatrix(distanceMatrix.rows, distanceMatrix.cols);
        IntStream.range(0, densityMatrix.size()).forEach(i -> d.setRow(i, densityMatrix.get(i)));
        return d;
    }

    public static DistanceMatrix computeInverseDensityMatrix(DistanceMatrix distanceMatrix){
        DistanceMatrix inverseDensity = new ArrayDistanceMatrix(distanceMatrix.rows, distanceMatrix.cols);
        DistanceMatrix densityMatrix = computeDensityMatrix(distanceMatrix);
        IntStream.range(0, densityMatrix.rows).forEach(r -> {
            double[] row = densityMatrix.getRow(r);
            List<Double> list = new ArrayList<>();
            for(double d : row){
                if(d == 0.0){
                    list.add(Double.MAX_VALUE);
                }else{
                    list.add(1.0/d);
                }
            }
            inverseDensity.setRow(r, list);
        });
        return inverseDensity;
    }

    public static DistanceMatrix computeKNNMatrix(DistanceMatrix distanceMatrix){
        List<List<Double>> densityMatrix = IntStream.range(0, distanceMatrix.rows).parallel().mapToObj(r -> {
            double[] row = distanceMatrix.getRow(r);
            List<Pair<Double, Integer>> elements = new ArrayList<>();
            for(int i=0;i<row.length;i++){
                elements.add(new Pair<>(row[i], i));
            }
            Collections.sort(elements, (o1, o2) -> o1._1().compareTo(o2._1()));
            List<Pair<Double, Integer>> knn = new ArrayList<>();
            for(int i=0;i<elements.size();i++){
                knn.add(new Pair<>((double)i, elements.get(i)._2()));
            }
            Collections.sort(knn, (o1, o2) -> o1._2()-o2._2());
            return knn.stream().map(e -> e._1()).collect(Collectors.toList());
        }).collect(Collectors.toList());
        DistanceMatrix d = new ArrayDistanceMatrix(distanceMatrix.rows, distanceMatrix.cols);
        IntStream.range(0, densityMatrix.size()).forEach(i -> d.setRow(i, densityMatrix.get(i)));
        return d;
    }
}

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

import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import topcat.matrix.distancematrix.DistanceMatrix;
import topcat.persistence.PersistenceModuleCollection;
import topcat.persistence.contours.PersistenceContour;
import topcat.persistence.contours.StandardContour;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

/**
 * Represents a point in some n-dimensional R-vector space.
 */
public class Point{
    List<Double> x;

    public Point(List<Double> x){
        this.x = x;
    }
    public Point(Double... x){
        this(Arrays.asList(x));
    }

    public List<Double> getX(){
        return this.x;
    }

    /**
     * Returns the Euclidean distance between two points.
     * @param p1
     * @param p2
     * @return
     */
    public static Double euclideanDistance(Point p1, Point p2){
        double d = 0;
        for(int i=0;i<p1.x.size();i++){
            double diff = p1.x.get(i)-p2.x.get(i);
            d += diff*diff;
        }
        return Math.sqrt(d);
    }

    /**
     * Sample points from a list with uniform distribution.
     * @param points
     * @param samples
     * @return
     */
    public static List<Point> samplePoints(List<Point> points, int samples){
        List<Point> sampledPoints = new ArrayList<>();
        IntOpenHashSet chosen = new IntOpenHashSet();
        Random r = new Random();
        for(int i=0;i<samples;i++){
            int next = r.nextInt(points.size());
            while(chosen.contains(next)){
                next = r.nextInt(points.size());
            }
            chosen.add(next);
            sampledPoints.add(points.get(next));
        }
        return sampledPoints;
    }

    /**
     * Returns 'N' number of equidistant points on a circle.
     * @param radius
     * @param N
     * @return
     */
    public static List<Point> circle2D(double radius, int N){
        List<Point> points = new ArrayList<>();
        for(int i=0;i<N;i++){
            List<Double> coordinate = new ArrayList<>();
            double theta = i*2*Math.PI/((double)N);
            coordinate.add(radius*Math.cos(theta));
            coordinate.add(radius*Math.sin(theta));
            points.add(new Point(coordinate));
        }
        return points;
    }

    /**
     * Returns 'n' number of points samples from a 'd'-dimensional shpere.
     * @param n
     * @param d
     * @return
     */
    public static List<Point> getRandomSpherePoints(int n, int d) {
        List<Point> points = new ArrayList<>();
        Random r = new Random();
        for (int i = 0; i < n; i++) {
            List<Double> x = new ArrayList<>();
            double norm = 0;
            for(int j=0;j<d+1;j++){
                double val = r.nextDouble()-0.5;
                x.add(val);
                norm += val*val;
            }
            norm = Math.sqrt(norm);
            for(int j=0;j<d+1;j++){
                x.set(j, x.get(j)/norm);
            }
            points.add(new Point(x));
        }
        return points;
    }

    public static void main(String[] args){
        List<Point> points = Point.getRandomSpherePoints(100, 2);

        //Add distancematrices
        DistanceMatrix distanceMatrix = DistanceMatrix.computeEuclideanDistanceMatrix(points);
        DistanceMatrix densityMatrix = DistanceMatrix.computeKNNMatrix(distanceMatrix);
        List<DistanceMatrix> distanceMatrices = new ArrayList<>();
        distanceMatrices.add(distanceMatrix);
        distanceMatrices.add(densityMatrix);

        //Add filtrationvalues
        List<List<Double>> filtrationValues = new ArrayList<>();
        List<Double> radiusFiltrationValues = new ArrayList<>();
        for(int i=0;i<=5;i++){
            radiusFiltrationValues.add(i*0.2);
        }
        filtrationValues.add(radiusFiltrationValues);
        filtrationValues.add(Arrays.asList(new Double[]{0.0, 10.0, 20.0, 30.0, 40.0}));

        int maxDimension = 3; //Compute 0th and 1 homology.
        PersistenceModuleCollection persistenceModules = PersistenceModuleCollection.create(distanceMatrices, filtrationValues, maxDimension);

        PersistenceContour persistenceContour = new StandardContour(filtrationValues);
        persistenceModules.get(2).computeStableRank(filtrationValues.get(0), persistenceContour);
    }
}

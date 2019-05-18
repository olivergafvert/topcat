/*
Topcat - a multidimensional persistent homology library.
Copyright (C) 2019 Oliver Gäfvert

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

package topcat;

import org.junit.Assert;
import org.junit.Test;
import topcat.mains.PythonInterface;
import topcat.matrix.distancematrix.DistanceMatrix;
import topcat.persistence.PersistenceModuleCollection;
import topcat.persistence.contours.PersistenceContour;
import topcat.persistence.contours.StandardContour;
import topcat.persistence.stablerank.StableRankFunction;
import topcat.util.Pair;
import topcat.util.Point;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class PythonInterfaceTest {
    @Test
    public void computePersitenceModuleTest(){
        List<Point> points = Point.circle2D(1, 10);

        //Add distancematrices
        List<List<Double>> distanceMatrix = DistanceMatrix.computeEuclideanDistanceMatrix(points).toList();

        List<List<List<Double>>> distanceMatrices = new ArrayList<>();
        distanceMatrices.add(distanceMatrix);

        //Add filtrationvalues
        List<List<Double>> filtrationValues = new ArrayList<>();
        List<Double> radiusFiltrationValues = new ArrayList<>();
        for(int i=0;i<10;i++){
            radiusFiltrationValues.add(i*0.1);
        }
        filtrationValues.add(radiusFiltrationValues);

        int maxDimension = 2; //Compute 0th and 1 homology.
        PersistenceModuleCollection persistenceModules = PythonInterface.computePersistenceModules(distanceMatrices, filtrationValues, maxDimension, false);

        PersistenceContour persistenceContour = new StandardContour(filtrationValues);
        StableRankFunction barcode0 = new StableRankFunction();
        for(int i=0;i<7;i++){
            barcode0.add(new Pair<>(radiusFiltrationValues.get(i), 10));
        }
        for(int i=7;i<10;i++){
            barcode0.add(new Pair<>(radiusFiltrationValues.get(i), 1));
        }
        StableRankFunction barcode0_test = persistenceModules.get(0).computeStableRank(persistenceModules.get(0).getFiltrationValues().get(0), persistenceContour);
        Assert.assertEquals(barcode0, barcode0_test);

        StableRankFunction barcode1 = new StableRankFunction();
        for(int i=0;i<10;i++){
            barcode1.add(new Pair<>(radiusFiltrationValues.get(i), 1));
        }
        StableRankFunction barcode1_test = persistenceModules.get(1).computeStableRank(persistenceModules.get(1).getFiltrationValues().get(0), persistenceContour);
        Assert.assertEquals(barcode1, barcode1_test);
    }

    @Test
    public void noisyCircleTest(){
        List<Point> points = Point.getRandomSpherePoints(300, 2);
        List<List<Double>> points_d = new ArrayList<>();
        for(Point p : points) points_d.add(p.getX());

        //Add filtrationvalues
        List<List<Double>> filtrationValues = new ArrayList<>();
        List<Double> radiusFiltrationValues = new ArrayList<>();
        for(int i=0;i<6;i++){
            radiusFiltrationValues.add(i*0.1);
        }
        filtrationValues.add(radiusFiltrationValues);
        filtrationValues.add(Arrays.asList(new Double[]{-100.0, -50.0, -40.0, -30.0, -20.0, -15.0, -10.0, -1.0, 0.0}));

        int maxDimension = 3; //Compute 0th and 1 homology.
        PersistenceModuleCollection persistenceModules = PythonInterface.computePersistenceModules(points_d, Arrays.asList(new String[]{"euclidean_codensity"}), filtrationValues, maxDimension, false);
    }
}

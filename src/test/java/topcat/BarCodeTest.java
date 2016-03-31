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

package topcat;

import org.junit.Assert;
import org.junit.Test;
import topcat.persistence.PersistenceModuleCollection;
import topcat.persistence.barcode.BasicBarcode;
import topcat.persistence.noise.Noise;
import topcat.persistence.noise.StandardNoise;
import topcat.util.DistanceMatrix;
import topcat.util.Pair;
import topcat.util.Point;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by oliver on 2016-02-29.
 */
public class BarCodeTest {

    @Test
    public void circleTest(){
        List<Point> points = Point.circle2D(1, 10);

        //Add distancematrices
        DistanceMatrix distanceMatrix = DistanceMatrix.computeDistanceMatrix(points, Point::euclideanDistance);
        List<DistanceMatrix> distanceMatrices = new ArrayList<>();
        distanceMatrices.add(distanceMatrix);

        //Add filtrationvalues
        List<List<Double>> filtrationValues = new ArrayList<>();
        List<Double> radiusFiltrationValues = new ArrayList<>();
        for(int i=0;i<10;i++){
            radiusFiltrationValues.add(i*0.1);
        }
        filtrationValues.add(radiusFiltrationValues);

        int maxDimension = 2; //Compute 0th and 1 homology.
        PersistenceModuleCollection persistenceModules = PersistenceModuleCollection.create(distanceMatrices, filtrationValues, maxDimension);

        Noise noise = new StandardNoise();
        BasicBarcode barcode0 = new BasicBarcode();
        for(int i=0;i<7;i++){
            barcode0.add(new Pair<>(radiusFiltrationValues.get(i), 10));
        }
        for(int i=7;i<10;i++){
            barcode0.add(new Pair<>(radiusFiltrationValues.get(i), 1));
        }
        BasicBarcode barcode0_test = noise.computeBasicBarcode(persistenceModules.get(0).getFunctor(), persistenceModules.get(0).getFiltrationValues());
        Assert.assertEquals(barcode0, barcode0_test);

        BasicBarcode barcode1 = new BasicBarcode();
        for(int i=0;i<10;i++){
            barcode1.add(new Pair<>(radiusFiltrationValues.get(i), 1));
        }
        BasicBarcode barcode1_test = noise.computeBasicBarcode(persistenceModules.get(1).getFunctor(), persistenceModules.get(1).getFiltrationValues());
        Assert.assertEquals(barcode1, barcode1_test);
    }
}

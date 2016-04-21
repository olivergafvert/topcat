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

package topcat.dataset.caltech101;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import topcat.dataset.util.Image;
import topcat.persistence.PersistenceModuleCollection;
import topcat.persistence.noise.Noise;
import topcat.persistence.noise.StandardNoise;
import topcat.matrix.distancematrix.DistanceMatrix;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Created by oliver on 2016-04-05.
 */
public class CALTECH101 {
    private static Logger log = LoggerFactory.getLogger(CALTECH101.class);

    public static void main(String[] args){
        Image img = Image.loadImageRGB("local/Caltech101/101_ObjectCategories/cougar_body/image_00351.jpg");
        log.info(img.width + " " + img.height);
        DistanceMatrix d1 = Image.getDistanceMatrix(img, 0, 1);
        DistanceMatrix d2 = Image.getDistanceMatrix(img, 1, 1);
        DistanceMatrix d3 = Image.getDistanceMatrix(img, 2, 1);
        List<DistanceMatrix> distanceMatrices = new ArrayList<>();
        distanceMatrices.add(d1);
        distanceMatrices.add(d2);
        //distanceMatrices.add(d3);

        List<List<Double>> filtrationValues = new ArrayList<>();
        List<Double> filtrationValues1 = Arrays.asList(new Double[]{0.0, 10.0, 20.0});
        List<Double> filtrationValues2 = Arrays.asList(new Double[]{0.0, 10.0, 20.0});
        List<Double> filtrationValues3 = Arrays.asList(new Double[]{0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0});
        filtrationValues.add(filtrationValues1);
        filtrationValues.add(filtrationValues2);
        //filtrationValues.add(filtrationValues3);

        PersistenceModuleCollection persistenceModules = PersistenceModuleCollection.create(distanceMatrices,
                filtrationValues, 2);
        Noise noise = new StandardNoise(10000000, 10);
        noise.computeBasicBarcode(persistenceModules.get(1).getFunctor(), persistenceModules.get(1).getFiltrationValues());
    }
}

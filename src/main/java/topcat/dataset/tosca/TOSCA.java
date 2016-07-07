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

package topcat.dataset.tosca;

import topcat.persistence.PersistenceModuleCollection;
import topcat.persistence.noise.Noise;
import topcat.persistence.noise.StandardNoise;
import topcat.matrix.distancematrix.DistanceMatrix;
import topcat.util.Point;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Created by oliver on 2016-03-13.
 */
public class TOSCA {

    public static List<Point> readPointsFromFile(String file) throws IOException{
        BufferedReader reader = new BufferedReader(new FileReader(new File(file)));
        String line;
        List<Point> points = new ArrayList<>();
        while((line = reader.readLine()) != null){
            List<Double> point = new ArrayList<>();
            for(String a : line.trim().split(" ")){
                point.add(Double.parseDouble(a));
            }
            points.add(new Point(point));
        }
        return points;
    }

    public static void main(String[] args){
        try {
            System.out.println("Reading points...");
            List<Point> points = Point.samplePoints(readPointsFromFile("local/TOSCA/cat0.vert"), 200);
            System.out.println("Finished reading points. Computing distance matrix...");
            DistanceMatrix distanceMatrix = DistanceMatrix.computeEuclideanDistanceMatrix(points);
            System.out.println("Finished computing distance matrix. Computing density matrix...");
            DistanceMatrix densityMatrix = DistanceMatrix.computeKNNMatrix(distanceMatrix);
            System.out.println("Finished computing density matrix.");
            List<List<Double>> filtrationValues = new ArrayList<>();
            List<Double> filtrationValues1 = Arrays.asList(new Double[]{0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0});
            List<Double> filtrationValues2 = Arrays.asList(new Double[]{0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0});
            filtrationValues.add(filtrationValues1);
            filtrationValues.add(filtrationValues2);

            List<DistanceMatrix> distanceMatrices = new ArrayList<>();
            distanceMatrices.add(distanceMatrix);
            distanceMatrices.add(densityMatrix);
            PersistenceModuleCollection persistenceModules = PersistenceModuleCollection.create(distanceMatrices,
                    filtrationValues, 2);
            Noise noise = new StandardNoise(10000000, 10);
            noise.computeFCF(persistenceModules.get(1).getFunctor(), persistenceModules.get(1).getFiltrationValues());
        }catch (IOException ioe){
            ioe.printStackTrace();
        }
    }

}

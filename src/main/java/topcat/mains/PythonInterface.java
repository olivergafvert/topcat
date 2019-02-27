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

package topcat.mains;
import py4j.GatewayServer;
import topcat.matrix.distancematrix.ArrayDistanceMatrix;
import topcat.matrix.distancematrix.DistanceMatrix;
import topcat.persistence.PersistenceModule;
import topcat.persistence.PersistenceModuleCollection;
import topcat.persistence.contours.ExponentialContour;
import topcat.persistence.contours.PersistenceContour;
import topcat.persistence.contours.ProductContour;
import topcat.persistence.contours.StandardContour;
import topcat.persistence.contours.kernels.KernelFunction;
import topcat.persistence.contours.kernels.StepKernelFunction;
import topcat.persistence.simplex.SimplexStorageStructure;
import topcat.persistence.simplex.SimplicialComplex;
import topcat.persistence.stablerank.StableRankFunction;
import topcat.util.Point;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * A python interface to access topcat via python.
 */
public class PythonInterface {

    public static List<Point> transformPoints(List<List<Double>> points){
        List<Point> _points = new ArrayList<>();
        for(int i=0;i<points.size();i++){
            _points.add(new Point(points.get(i)));
        }
        return _points;
    }

    public static DistanceMatrix euclideanDistanceMatrix(List<List<Double>> points){
        return DistanceMatrix.computeEuclideanDistanceMatrix(transformPoints(points));
    }

    public static PersistenceModuleCollection computePersistenceModules(List<List<List<Double>>> distanceMatrices, List<List<Double>> filtrationValues, Integer maxDimension){
        List<DistanceMatrix> _distanceMatrices = new ArrayList<>();
        for(List<List<Double>> dmat : distanceMatrices) _distanceMatrices.add(new ArrayDistanceMatrix(dmat));
        SimplexStorageStructure simplexStorageStructure = SimplicialComplex.computeSimplexStream(_distanceMatrices, filtrationValues, maxDimension);
        PersistenceModuleCollection persistenceModules = PersistenceModuleCollection.create(simplexStorageStructure, filtrationValues, maxDimension);
        return persistenceModules;
    }

    public static List<List<List<Double>>> computeStableRank(List<List<List<Double>>> distanceMatrices, List<List<Double>> filtrationValues, Integer maxDimension) {
        return computeStableRank(distanceMatrices, filtrationValues, maxDimension, new StandardContour(filtrationValues));
    }

    public static List<List<List<Double>>> computeStableRank(List<List<List<Double>>> distanceMatrices, List<List<Double>> filtrationValues, Integer maxDimension, List<List<List<Double>>> stepkernels) {
        List<KernelFunction> kernels = new ArrayList<>();
        for(int i=0;i<stepkernels.size();i++){
            kernels.add(new StepKernelFunction(stepkernels.get(i).get(0), stepkernels.get(i).get(1)));
        }
        PersistenceContour contour = new ProductContour(filtrationValues, kernels);
        return computeStableRank(distanceMatrices, filtrationValues, maxDimension, contour);
    }

    public static List<List<List<Double>>> computeStableRank(List<List<List<Double>>> distanceMatrices, List<List<Double>> filtrationValues, Integer maxDimension, PersistenceContour persistenceContour){
        List<DistanceMatrix> _distanceMatrices = new ArrayList<>();
        for(List<List<Double>> dmat : distanceMatrices) _distanceMatrices.add(new ArrayDistanceMatrix(dmat));
        SimplexStorageStructure simplexStorageStructure = SimplicialComplex.computeSimplexStream(_distanceMatrices, filtrationValues, maxDimension);
        PersistenceModuleCollection persistenceModules = PersistenceModuleCollection.create(simplexStorageStructure, filtrationValues, maxDimension);
        List<List<List<Double>>> stableRankFunctions = new ArrayList<>();
        for(int i=0;i<=persistenceModules.getMaxDimension();i++){
            PersistenceModule persistenceModule = persistenceModules.get(i);
            StableRankFunction stableRankFunction = persistenceModule.computeStableRank(filtrationValues.get(0), persistenceContour);
            stableRankFunctions.add(stableRankFunction.toList());
            System.out.println(i+": "+ stableRankFunction);
        }
        return stableRankFunctions;
    }

    public static void main(String[] args) {
        PythonInterface app = new PythonInterface();
        GatewayServer server = new GatewayServer(app);
        server.start();
    }
}

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
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
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
import topcat.persistence.simplex.SparseSimplexStorageStructure;
import topcat.persistence.stablerank.StableRankFunction;
import topcat.util.Point;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * A python interface to access topcat via python.
 */
public class PythonInterface {
    private static Logger log = LoggerFactory.getLogger(PythonInterface.class);

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

    public static List<DistanceMatrix> parseDistances(List<Point> points, List<String> distances){
        List<DistanceMatrix> distanceMatrices = new ArrayList<>();
        for(String distance : distances){
            String dist = distance.trim().toLowerCase();
            if(dist.equals("euclidean")){
                distanceMatrices.add(DistanceMatrix.computeEuclideanDistanceMatrix(points));
            }else if(dist.equals("euclidean_codensity")){
                distanceMatrices.add(DistanceMatrix.computeEuclideanDistanceMatrix(points));
                distanceMatrices.add(DistanceMatrix.codensityMatrix(distanceMatrices.get(distanceMatrices.size()-1)));
            }else if(dist.length() == 1){
                try{
                    int  axis = Integer.parseInt(dist);
                    distanceMatrices.add(DistanceMatrix.computeAxisDistanceMatrix(points, axis));
                }catch(NumberFormatException nfe){
                    log.debug("Could not parse distance: "+dist);
                    log.error(nfe.toString());
                    return null;
                }
            }
        }
        return distanceMatrices;
    }

    public static PersistenceModuleCollection computePersistenceModules(List<List<Double>> points, List<String> distances, List<List<Double>> filtrationValues, Integer maxDimension, boolean sparse){
        List<Point> t_points = transformPoints(points);
        List<DistanceMatrix> _distanceMatrices = parseDistances(t_points, distances);
        SimplexStorageStructure simplexStorageStructure = SimplicialComplex.computeSimplexStream(_distanceMatrices, filtrationValues, maxDimension, sparse);
        PersistenceModuleCollection persistenceModules = PersistenceModuleCollection.create(simplexStorageStructure, filtrationValues, maxDimension);
        return persistenceModules;
    }

    public static PersistenceModuleCollection computePersistenceModules(List<List<List<Double>>> distanceMatrices, List<List<Double>> filtrationValues, Integer maxDimension, boolean sparse){
        List<DistanceMatrix> _distanceMatrices = new ArrayList<>();
        for(List<List<Double>> dmat : distanceMatrices) _distanceMatrices.add(new ArrayDistanceMatrix(dmat));
        SimplexStorageStructure simplexStorageStructure = SimplicialComplex.computeSimplexStream(_distanceMatrices, filtrationValues, maxDimension, sparse);
        PersistenceModuleCollection persistenceModules = PersistenceModuleCollection.create(simplexStorageStructure, filtrationValues, maxDimension);
        return persistenceModules;
    }

    public static List<List<List<Double>>> computeStableRank(List<List<Double>> points, List<String> distances, List<List<Double>> filtrationValues, Integer maxDimension) {
        return computeStableRank(points, distances, filtrationValues, maxDimension, null);
    }

    public static List<List<List<Double>>> computeStableRank(List<List<Double>> points, List<String> distances, List<List<Double>> filtrationValues, Integer maxDimension, List<List<List<Double>>> stepkernels) {
        List<Point> t_points = transformPoints(points);
        List<DistanceMatrix> distanceMatrices = parseDistances(t_points, distances);
        PersistenceContour contour;
        if(stepkernels == null){
            contour = new StandardContour(filtrationValues);
        }else {
            List<KernelFunction> kernels = new ArrayList<>();
            for (int i = 0; i < stepkernels.size(); i++) {
                kernels.add(new StepKernelFunction(stepkernels.get(i).get(0), stepkernels.get(i).get(1)));
            }
            contour = new ProductContour(filtrationValues, kernels);
        }
        return computeStableRank(distanceMatrices, filtrationValues, maxDimension, contour);
    }

    public static List<List<List<Double>>> computeStableRank(List<List<List<Double>>> distanceMatrices, List<List<Double>> filtrationValues, Integer maxDimension) {
        return computeStableRank(distanceMatrices, filtrationValues, maxDimension, null);
    }

    public static List<List<List<Double>>> computeStableRank(List<List<List<Double>>> distanceMatrices, List<List<Double>> filtrationValues, Integer maxDimension, List<List<List<Double>>> stepkernels) {
        List<DistanceMatrix> _distanceMatrices = new ArrayList<>();
        for(List<List<Double>> dmat : distanceMatrices) _distanceMatrices.add(new ArrayDistanceMatrix(dmat));
        PersistenceContour contour;
        if(stepkernels == null){
            contour = new StandardContour(filtrationValues);
        }else {
            List<KernelFunction> kernels = new ArrayList<>();
            for (int i = 0; i < stepkernels.size(); i++) {
                kernels.add(new StepKernelFunction(stepkernels.get(i).get(0), stepkernels.get(i).get(1)));
            }
            contour = new ProductContour(filtrationValues, kernels);
        }
        return computeStableRank(_distanceMatrices, filtrationValues, maxDimension, contour);
    }

    public static List<List<List<Double>>> computeStableRank(List<DistanceMatrix> distanceMatrices, List<List<Double>> filtrationValues, Integer maxDimension, PersistenceContour persistenceContour){
        SimplexStorageStructure simplexStorageStructure = SimplicialComplex.computeSimplexStream(distanceMatrices, filtrationValues, maxDimension);
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

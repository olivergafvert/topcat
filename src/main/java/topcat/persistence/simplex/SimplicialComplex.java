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

package topcat.persistence.simplex;


import gnu.trove.iterator.TIntIterator;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.hash.TIntHashSet;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import topcat.persistence.PersistenceModuleCollection;
import topcat.persistence.noise.Noise;
import topcat.persistence.noise.StandardNoise;
import topcat.util.DistanceMatrix;
import topcat.util.Pair;
import topcat.util.Point;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Utils for constructing a simplicial complex for distance matrices and filtration values.
 */
public class SimplicialComplex {
    private static Logger log = LoggerFactory.getLogger(SimplicialComplex.class);

    /**
     * Computes the Vietoris-Rips for the filtration value given by 'filtrationValue'.
     * @param distanceMatrices
     * @param filtrationValue
     * @param maxDimension
     * @return
     */
    public static List<Simplex> computeVietorisRipsComplex(List<DistanceMatrix> distanceMatrices, List<Double> filtrationValue, int maxDimension){
        DistanceMatrix distanceMatrix = distanceMatrices.get(0);
        List<Simplex> simplices = new ArrayList<>();

        //Add 0-dimensional simplices
        for(int i=0;i<distanceMatrix.rows;i++){
            simplices.add(new Simplex(i));
        }

        //Initialize neighbor sets
        TIntObjectHashMap<TIntHashSet> baseVertices = new TIntObjectHashMap<>();
        for(int i=0;i<distanceMatrix.rows;i++){
            baseVertices.put(i, new TIntHashSet());
        }

        //Compute all edges that will appear in the multifiltration and add them as 1-simplices
        IntStream.range(0, distanceMatrix.rows).parallel().mapToObj(i -> {
            TIntHashSet upperNeighbors = new TIntHashSet();
            List<Simplex> local_simplices = new ArrayList<>();
            for(int j=i+1; j<distanceMatrix.rows;j++){
                //Check if edge will appear in the multifiltration
                boolean isLEQ = true;
                for(int k=0;k<filtrationValue.size();k++){
                    if(distanceMatrices.get(k).get(i, j) > filtrationValue.get(k)){
                        isLEQ = false;
                        break;
                    }
                }
               if(isLEQ){
                   upperNeighbors.add(j);
                   local_simplices.add(new Simplex(i, j));
               }
            }
            return new Pair<>(new Pair<>(i, upperNeighbors), local_simplices);
        }).collect(Collectors.toList()) //Collect first to make it thread safe
                .stream()
                .forEach(p -> {
                    baseVertices.put(p._1()._1(), p._1()._2());
                    simplices.addAll(p._2());
                });


        //Add cofaces to each vertex inductively.
        int[] keys = baseVertices.keys();
        IntStream.range(0, keys.length).parallel().mapToObj(k -> {
            List<Simplex> local_simplices = new ArrayList<>();
            int[] vertices = new int[maxDimension + 1];
            vertices[0] = keys[k];
            addCofaces(local_simplices, baseVertices.get(keys[k]), vertices, 0, baseVertices);
            return local_simplices;
        }).collect(Collectors.toList()) //Collect first to make it thread safe
                .stream()
                .forEach(l -> simplices.addAll(l));


        return simplices;
    }

    /**
     * Inductively adds the cofaces to a simplex.
     * @param simplices
     * @param candidates
     * @param vertices
     * @param index
     * @param local_baseVertices
     */
    private static void addCofaces(List<Simplex> simplices, TIntHashSet candidates, int[] vertices, int index, TIntObjectHashMap<TIntHashSet> local_baseVertices){
        if(index >= 2 && index < vertices.length){
            List<Integer> t_vertices = new ArrayList<>(index+1);
            for(int i=0;i<index+1;i++){
                t_vertices.add(vertices[i]);
            }
            simplices.add(new Simplex(t_vertices));
        }
        if(index < vertices.length-1) {
            TIntIterator iterator = candidates.iterator();
            while(iterator.hasNext()) {
                vertices[index + 1] = iterator.next();
                TIntHashSet lower_candidates = intersect(candidates, local_baseVertices.get(vertices[index+1]));
                addCofaces(simplices, lower_candidates, vertices, index + 1, local_baseVertices);
            }
        }
    }

    /**
     * Computes the intersection of two TIntHashSets.
     * @param h1
     * @param h2
     * @return
     */
    private static TIntHashSet intersect(TIntHashSet h1, TIntHashSet h2){
        TIntHashSet intersection = new TIntHashSet();
        TIntIterator iterator = h2.iterator();
        while(iterator.hasNext()){
            int n = iterator.next();
            if(h1.contains(n)){
                intersection.add(n);
            }
        }
        return intersection;
    }

    /**
     * Computes a multifiltered simplicial complex.
     * @param distanceMatrices
     * @param filtrationValues
     * @param maxDimension
     * @return
     */
    public static SimplexStorageStructure computeSimplexStream(List<DistanceMatrix> distanceMatrices, List<List<Double>> filtrationValues, int maxDimension){
        //Pick out the largest filtration value in each direction
        List<Double> maxFiltrationValues = new ArrayList<>();
        for(int i=0;i<filtrationValues.size();i++){
            maxFiltrationValues.add(filtrationValues.get(i).get(filtrationValues.get(i).size()-1));
        }
        SimplexStorageStructure storageStructure = new SimplexStorageStructure(filtrationValues);
        log.debug("Starting to compute simplicial complex...");
        //Compute the Vietoris-Rips complex for the maximal filtration value
        List<Simplex> simplices = computeVietorisRipsComplex(distanceMatrices, maxFiltrationValues, maxDimension);
        log.debug("Finished computing simplicial complex. (Computed "+simplices.size()+" number of simplices.)");
        log.debug("Starting to compute filtrationIndices for each simplex...");
        //Compute the filtration indices for each simplex.
        for(int i=0;i<simplices.size();i++){
            List<Integer> filtrationIndexes = calcFiltrationIndexes(simplices.get(i), distanceMatrices, filtrationValues);
            if (filtrationIndexes != null) {
                storageStructure.addElement(simplices.get(i), filtrationIndexes);
            }
        }
        log.debug("Finished computing filtrationIndices.");
        return storageStructure;
    }

    /**
     * Calculates the filtrationIndices of a simplex.
     * @param simplex
     * @param distanceMatrices
     * @param filtrationValues
     * @return
     */
    public static List<Integer> calcFiltrationIndexes(Simplex simplex, List<DistanceMatrix> distanceMatrices, List<List<Double>> filtrationValues){
        //Find maximum value of the weights on the edges for each metric
        List<Integer> vertices = simplex.getVertices();


        return IntStream.range(0, filtrationValues.size()).mapToObj(k -> {
            DistanceMatrix distanceMatrix = distanceMatrices.get(k);
            double f_max = -1;
            for(int i=0;i<vertices.size();i++){
                for(int j=i;j<vertices.size();j++){
                    double f = distanceMatrix.get(vertices.get(i), vertices.get(j));
                    if(f > f_max){
                        f_max = f;
                    }
                }
            }

            List<Double> filtrationvalues_k = filtrationValues.get(k);
            int filtrationIndex = 0, i=0;
            while(i<filtrationvalues_k.size() && f_max > filtrationvalues_k.get(i)){
                filtrationIndex = ++i;
            }
            return filtrationIndex;
        }).collect(Collectors.toList());
    }

    public static void main(String[] args){
        List<Point> points = Point.circle2D(1, 50);
        DistanceMatrix distanceMatrix = DistanceMatrix.computeDistanceMatrix(points, topcat.util.Point::euclideanDistance);
        List<List<Double>> filtrationValues = new ArrayList<>();
        filtrationValues.add(Arrays.asList(new Double[]{0.0, 0.1, 0.3, 1.8, 1.85, 1.9, 2.0, 2.1, 50.0}));
        List<DistanceMatrix> distanceMatrices = new ArrayList<>();
        distanceMatrices.add(distanceMatrix);
        PersistenceModuleCollection persistenceModules = PersistenceModuleCollection.create(distanceMatrices, filtrationValues, 2);
        Noise noise = new StandardNoise();
        noise.computeBasicBarcode(persistenceModules.get(1).getFunctor(), persistenceModules.get(1).getFiltrationValues());
    }
}

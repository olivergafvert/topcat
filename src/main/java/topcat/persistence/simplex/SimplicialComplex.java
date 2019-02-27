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

import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.ints.IntIterator;
import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import topcat.persistence.PersistenceModuleCollection;
import topcat.persistence.contours.PersistenceContour;
import topcat.persistence.contours.StandardContour;
import topcat.matrix.distancematrix.DistanceMatrix;
import topcat.util.BinomialCoeffTable;
import topcat.util.IntTuple;
import topcat.util.Point;

import java.util.*;

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
    public static List<Simplex> computeVietorisRipsComplex(List<DistanceMatrix> distanceMatrices, List<Double> filtrationValue, int maxDimension, BinomialCoeffTable binom_coeff){
        DistanceMatrix distanceMatrix = distanceMatrices.get(0);
        List<Simplex> simplices = new ArrayList<>();

        //Add 0-dimensional simplices
        for(int i=0;i<distanceMatrix.rows;i++){
            simplices.add(new Simplex(i, 0));
        }

        //Initialize neighbor sets
        Int2ObjectOpenHashMap<IntOpenHashSet> baseVertices = new Int2ObjectOpenHashMap<>();
        for(int i=0;i<distanceMatrix.rows;i++){
            baseVertices.put(i, new IntOpenHashSet());
        }

        //Compute all edges that will appear in the multifiltration and add them as 1-simplices
        int[] nonzeros_rows = distanceMatrix.getNonZeroRows();

        for(int r=0; r < nonzeros_rows.length; r++) {
            int i = nonzeros_rows[r];
            IntOpenHashSet upperNeighbors = new IntOpenHashSet();
            List<Simplex> local_simplices = new ArrayList<>();
            int[] nonzero_cols = distanceMatrix.getNonZeroRowEntries(i);
            for (int j : nonzero_cols) {
                if (j <= i) {
                    continue;
                }
                //Check if edge will appear in the multifiltration
                boolean isLEQ = true;
                for (int k = 0; k < filtrationValue.size(); k++) {
                    if (distanceMatrices.get(k).get(i, j) > filtrationValue.get(k)) {
                        isLEQ = false;
                        break;
                    }
                }
                if (isLEQ) {
                    upperNeighbors.add(j);
                    local_simplices.add(new Simplex(binom_coeff.computeIndex(i, j), 1));
                }
            }
            baseVertices.put(i, upperNeighbors);
            simplices.addAll(local_simplices);
        }


        //Add cofaces to each vertex inductively.
        int[] keys = baseVertices.keySet().toIntArray();
        for(int k = 0; k < keys.length; k++){
            List<Simplex> local_simplices = new ArrayList<>();
            int[] vertices = new int[maxDimension + 1];
            vertices[0] = keys[k];
            addCofaces(local_simplices, baseVertices.get(keys[k]), vertices, 0, baseVertices, binom_coeff);
            simplices.addAll(local_simplices);
        }


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
    private static void addCofaces(List<Simplex> simplices, IntOpenHashSet candidates, int[] vertices, int index, Int2ObjectOpenHashMap<IntOpenHashSet> local_baseVertices, BinomialCoeffTable binomial_coeff){
        if(index >= 2 && index < vertices.length){
            List<Integer> t_vertices = new ArrayList<>(index+1);
            for(int i=0;i<index+1;i++){
                t_vertices.add(vertices[i]);
            }
            simplices.add(new Simplex(binomial_coeff.computeIndex(t_vertices), t_vertices.size()-1));
        }
        if(index < vertices.length-1) {
            IntIterator iterator = candidates.iterator();
            while(iterator.hasNext()) {
                vertices[index + 1] = iterator.nextInt();
                IntOpenHashSet lower_candidates = intersect(candidates, local_baseVertices.get(vertices[index+1]));
                addCofaces(simplices, lower_candidates, vertices, index + 1, local_baseVertices, binomial_coeff);
            }
        }
    }

    /**
     * Computes the intersection of two TIntHashSets.
     * @param h1
     * @param h2
     * @return
     */
    private static IntOpenHashSet intersect(IntOpenHashSet h1, IntOpenHashSet h2){
        IntOpenHashSet intersection = new IntOpenHashSet();
        IntIterator iterator = h2.iterator();
        while(iterator.hasNext()){
            int n = iterator.nextInt();
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
        IntTuple gridSize = IntTuple.zeros(filtrationValues.size());
        for(int i=0;i<filtrationValues.size();i++){
            gridSize.set(i, filtrationValues.get(i).size()-1);
        }
        SimplexStorageStructure storageStructure = new SimplexStorageStructure(filtrationValues, gridSize, maxDimension, distanceMatrices.get(0).cols);
        log.debug("Starting to compute simplicial complex...");
        //Compute the Vietoris-Rips complex for the maximal filtration value
        List<Simplex> simplices = computeVietorisRipsComplex(distanceMatrices, maxFiltrationValues, maxDimension, storageStructure.binomialCoeffTable);
        log.debug("Finished computing simplicial complex. (Computed "+simplices.size()+" number of simplices.)");
        log.debug("Starting to compute filtrationValues for each simplex...");
        //Compute the filtration indices for each simplex.
        for(int i=0;i<simplices.size();i++){
            List<Integer> filtrationIndexes = calcFiltrationIndexes(simplices.get(i), distanceMatrices, filtrationValues, storageStructure.binomialCoeffTable);
            if (filtrationIndexes != null) {
                storageStructure.addElement(simplices.get(i), new IntTuple(filtrationIndexes));
            }else{
                log.error("Failed to add simplex "+simplices.get(i)+" to simplex storage structure");
            }
        }
        log.debug("Finished computing filtrationValues.");
        return storageStructure;
    }

    /**
     * Calculates the filtrationValues of a simplex.
     * @param simplex
     * @param distanceMatrices
     * @param filtrationValues
     * @return
     */
    public static List<Integer> calcFiltrationIndexes(Simplex simplex, List<DistanceMatrix> distanceMatrices, List<List<Double>> filtrationValues, BinomialCoeffTable binomial_coeff){
        //Find maximum value of the weights on the edges for each metric
        int[] vertices = Simplex.get_simplex_vertices(simplex.getIndex(), simplex.getDimension(), distanceMatrices.get(0).cols-1, binomial_coeff);

        List<Integer> filtrationIndices = new ArrayList<>();
        for(int k=0 ; k < filtrationValues.size(); k++){
            DistanceMatrix distanceMatrix = distanceMatrices.get(k);
            double f_max = -1;
            for(int i=0;i<vertices.length;i++){
                for(int j=i;j<vertices.length;j++){
                    double f = distanceMatrix.get(vertices[i], vertices[j]);
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
            filtrationIndices.add(filtrationIndex);
        }
        return filtrationIndices;
    }

    public static void main(String[] args){
        List<Point> points = Point.circle2D(1, 50);
        DistanceMatrix distanceMatrix = DistanceMatrix.computeEuclideanDistanceMatrix(points);
        List<List<Double>> filtrationValues = new ArrayList<>();
        filtrationValues.add(Arrays.asList(new Double[]{0.0, 0.1, 0.3, 1.8, 1.85, 1.9, 2.0, 2.1, 50.0}));
        List<DistanceMatrix> distanceMatrices = new ArrayList<>();
        distanceMatrices.add(distanceMatrix);
        PersistenceModuleCollection persistenceModules = PersistenceModuleCollection.create(distanceMatrices, filtrationValues, 2);
        PersistenceContour persistenceContour = new StandardContour(filtrationValues);
        persistenceModules.get(1).computeStableRank(filtrationValues.get(0), persistenceContour);
    }
}

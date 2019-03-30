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

package topcat.persistence.homology;

import it.unimi.dsi.fastutil.objects.Object2IntOpenHashMap;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import topcat.matrix.GradedColumn;
import topcat.persistence.functor.Functor;
import topcat.persistence.functor.Nat;
import topcat.persistence.functor.exception.MalformedFunctorException;
import topcat.matrix.BMatrix;
import topcat.matrix.exception.NoSolutionException;
import topcat.matrix.exception.WrongDimensionException;
import topcat.persistence.simplex.Simplex;
import topcat.persistence.simplex.SimplexStorageStructure;
import topcat.util.Grid;
import topcat.util.GridIterator;
import topcat.util.IntTuple;
import topcat.util.Pair;
import topcat.util.paralelliterator.ParalellIntIterator;

import java.util.*;
import java.util.concurrent.*;

/**
 * Tools for computing the homology of a multifiltration.
 */
public class HomologyUtil {
    private static Logger log = LoggerFactory.getLogger(HomologyUtil.class);

    /**
     * Computes the inclusion map C^i_n \hookrightarrow C^{i+j}_n.
     * @param lower the simplices in C^i_n
     * @param current the simplices in C^{i+j}_n
     * @return
     */
    public static BMatrix computeInclusionMap(List<? extends Object> lower, List<? extends Object> current){
        if(lower.size() == 0){
            return new BMatrix(current.size(), lower.size());
        }
        BMatrix A = new BMatrix(current.size(), lower.size());
        Object2IntOpenHashMap<Object> lowerMap = new Object2IntOpenHashMap<>();
        for(int i=0;i<lower.size();i++){
            lowerMap.put(lower.get(i), i);
        }
        for(int i=0;i<current.size();i++){
            if(lowerMap.containsKey(current.get(i))){
                A.set(i, lowerMap.getInt(current.get(i)), true);
            }
        }
        return A;
    }

    /**
     * Computes the chain functors C_n: Q^2 \to Vect_K for each n <= 'maxdimension'
     * @param simplexStorageStructure
     * @param size
     * @param maxDimension
     * @return
     * @throws MalformedFunctorException
     * @throws WrongDimensionException
     */
    private static List<Grid<Integer>> computeChainFunctorDimensions(final SimplexStorageStructure simplexStorageStructure, IntTuple size, int maxDimension) throws MalformedFunctorException {
        log.debug("Starting to compute chain functor dimensions...");
        List<Grid<Integer>> chainFunctors = new ArrayList<>();

        for(int k=0;k<maxDimension+1;k++){
            final Grid<Integer> grid = Grid.create(size);
            //Compute the maps
            for(IntTuple v : GridIterator.getSequence(size)) {
                grid.set(v, simplexStorageStructure.getSimplicesLEQThan(k, v).size());

            }
            chainFunctors.add(grid);
        }
        log.debug("Finished to computing chain functor dimensions.");
        return chainFunctors;
    }

    /**
     * Finds indices where the chain complex does not change so we don't have to compute it's homology.
     * @param v
     * @param chainFunctorDims
     * @return
     */
    private static int findInvariantDimension(IntTuple v, List<Grid<Integer>> chainFunctorDims){
        List<IntTuple> basis = IntTuple.getStandardBasisSequence(v.length());
        for(int i=0;i<basis.size();i++){
            boolean is_invariant = true;
            for(int j=0;j<chainFunctorDims.size();j++){
                if(chainFunctorDims.get(j).get(v) != chainFunctorDims.get(j).get(v.minus(basis.get(i)))){
                    is_invariant = false;
                }
            }
            if(is_invariant) return i;
        }
        return -1;
    }

    static List<GradedColumn<Simplex>> columns_LEQThan(IntTuple v, Grid<List<GradedColumn<Simplex>>>  columns_grid){
        List<GradedColumn<Simplex>> columns = new ArrayList<>();
        for(IntTuple w : GridIterator.getSequence(v)) {
            if(columns_grid.get(w) != null) {
                columns.addAll(columns_grid.get(w));
            }
        }
        return columns;
    }


    /**
     * Computes the homology functors for each dimension less than 'maxdimension'.
     * @param simplexStorageStructure
     * @param size
     * @param maxDimension
     * @throws WrongDimensionException
     * @throws NoSolutionException
     */
    public static List<Functor> computeHomologyFunctors(final SimplexStorageStructure simplexStorageStructure, final IntTuple size, final int maxDimension) throws MalformedFunctorException, NoSolutionException{
        log.debug("Starting to compute homology functors...");
        List<Grid<Integer>> chainFunctorDims = computeChainFunctorDimensions(simplexStorageStructure, size, maxDimension);

        //The natural transformations from the chain functors to a basis change of the chain modules
        final List<Nat> naturalTransformation = new ArrayList<>();
        final List<Nat> naturalTransformation_inverse = new ArrayList<>();

        //The persistence modules, i.e homology of the multifiltration
        List<Functor> homfunctors = new ArrayList<>();

        for(int k=0;k<maxDimension;k++){
            naturalTransformation.add(new Nat(size));
            naturalTransformation_inverse.add(new Nat(size));
            homfunctors.add(new Functor(size));
        }

        //The dimension of the basis for the homology in C(v) for each v in N^r
        final List<Grid<Integer>> homologyDimension = new ArrayList<>();
        for(int i=0;i<maxDimension;i++) {
            homologyDimension.add(Grid.create(size));
        }

        log.debug("Starting to compute basis change in each position...");
        int n_threads = 2*(Runtime.getRuntime().availableProcessors() < 5 ? 5 : Runtime.getRuntime().availableProcessors());
        ExecutorService exec = Executors.newFixedThreadPool(n_threads);
        List<Future> futures = new ArrayList<>();
        List<HomologyWorker2> workers = new ArrayList<>();
        List<Pair<IntTuple, IntTuple>> invariantIndices = new ArrayList<>();

        //Sequential reduction
        Grid<List<HashMap<Simplex, GradedColumn<Simplex>>>> pivot_to_column_grid = Grid.create(size);
        Grid<List<List<GradedColumn<Simplex>>>> kernel_grid = Grid.create(size);
        List<Grid<List<GradedColumn<Simplex>>>> columns_grid = new ArrayList<>();
        for(int i=0;i<maxDimension;i++)
            columns_grid.add(Grid.create(size));
        //for(int s = 0 ; s<=size.max();s++){
        for(IntTuple v : GridIterator.getSequence(size)){
            log.debug("Starting reduction in grade: "+v);
            HomologyWorker2 worker = new HomologyWorker2(simplexStorageStructure, v, maxDimension);
            workers.add(worker);
            List<Pair<Integer, Pair<List<GradedColumn<Simplex>>, HashMap<Simplex, GradedColumn<Simplex>>>>> tt = new ArrayList<>();
            for(int i = 1 ; i<=maxDimension;i++){
                List<GradedColumn<Simplex>> columns = worker.computeBoundaries(simplexStorageStructure.getSimplicesAt(i, v));
                Collections.sort(columns);
                HashMap<Simplex, GradedColumn<Simplex>> pivot_to_columns = new HashMap<>();
                List<IntTuple> basis = IntTuple.getStandardBasisSequence(size.length());
                for(IntTuple b : basis){
                    if(pivot_to_column_grid.get(v.minus(b)) != null)
                        pivot_to_columns.putAll(pivot_to_column_grid.get(v.minus(b)).get(i-1));
                }
                List<GradedColumn<Simplex>> kernel = worker.reduce_matrix(columns, pivot_to_columns);
                //return new Pair<>(kernel, pivot_to_columns);
                tt.add(new Pair<>(i, new Pair<>(kernel, pivot_to_columns)));
            }

            List<HashMap<Simplex, GradedColumn<Simplex>>> pivot_to_columns = new ArrayList<>();
            List<List<GradedColumn<Simplex>>> kernelPartial = new ArrayList<>();
            for(int i=0;i<maxDimension;i++){
                pivot_to_columns.add(new HashMap<>());
                kernelPartial.add(new ArrayList<>());
            }
            for(Pair<Integer, Pair<List<GradedColumn<Simplex>> , HashMap<Simplex, GradedColumn<Simplex>>>> pair : tt){
                pivot_to_columns.get(pair._1()-1).putAll(pair._2()._2());
                kernelPartial.get(pair._1()-1).addAll(pair._2()._1());
            }
            pivot_to_column_grid.set(v, pivot_to_columns);
            kernel_grid.set(v, kernelPartial);
        }

        //Compute homology basis and basis change maps
        for(HomologyWorker2 worker : workers){
            //int invariantDim = findInvariantDimension(worker.v, chainFunctorDims);
            //if(invariantDim == -1) {
            List<List<GradedColumn<Simplex>>> images = new ArrayList<>();
            List<List<GradedColumn<Simplex>>> kernels = new ArrayList<>();
            List<List<Simplex>> basis = new ArrayList<>();
            for(int i=0;i<maxDimension;i++){

                images.add(new ArrayList<>(pivot_to_column_grid.get(worker.v).get(i).values()));
                kernels.add(new ArrayList<>());
                basis.add(new ArrayList<>(simplexStorageStructure.getSimplicesLEQThan(i, worker.v)));
            }
            for(IntTuple w : GridIterator.getSequence(worker.v)){
                for(int i=0;i<maxDimension;i++)
                    kernels.get(i).addAll(kernel_grid.get(w).get(i));
            }
            worker.images = images;
            worker.kernels = kernels;
            worker.basis = basis;
            //worker.run();
            futures.add(exec.submit(worker));
//            }else{
//                log.debug("Invariant index: "+worker.v);
//                invariantIndices.add(new Pair<>(worker.v, worker.v.minus(IntTuple.getStandardBasisElement(worker.v.length(), invariantDim))));
//            }
        }
        exec.shutdown();

        while(!exec.isTerminated()){
            try{
                exec.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
            }catch (InterruptedException ire){
                log.error("Failed to terminate executor service.", ire);
            }
        }
        try {
            for (Future future : futures) {
                future.get();
            }
        }catch (InterruptedException ire){
            log.error("Failed to terminate executor service.", ire);
        }catch (ExecutionException exe){
            log.error("Failed to terminate executor service.", exe);
        }

        for(HomologyWorker2 worker : workers){
            for(int k=0;k<maxDimension;k++){
                homologyDimension.get(k).set(worker.v, worker.homologyDimension[k]);
                naturalTransformation.get(k).setMap(worker.v, worker.naturalTransformation[k]);
                naturalTransformation_inverse.get(k).setMap(worker.v, worker.naturalTransformation_inverse[k]);
            }
        }
//        for(Pair<IntTuple, IntTuple> p : invariantIndices){
//            for(int k=0;k<maxDimension;k++){
//                homologyDimension.get(k).set(p._1(), homologyDimension.get(k).get(p._2()));
//                naturalTransformation.get(k).setMap(p._1(), new BMatrix(naturalTransformation.get(k).getMap(p._2())));
//                naturalTransformation_inverse.get(k).setMap(p._1(), new BMatrix(naturalTransformation_inverse.get(k).getMap(p._2())));
//            }
//        }
        log.debug("Finished computing basis change.");

        log.debug("Starting to apply basis change...");
        for(int k=0;k<homfunctors.size();k++){
            Functor H = homfunctors.get(k);
            for(IntTuple v : GridIterator.getSequence(size)){
                List<Simplex> currentSimplices = simplexStorageStructure.getSimplicesLEQThan(k, v);
                for(int i=0;i<v.length();i++){
                    BMatrix inclusionMap = computeInclusionMap(currentSimplices, simplexStorageStructure.getSimplicesLEQThan(k, v.plus(IntTuple.getStandardBasisElement(v.length(), i))));
                    if(v.get(i).equals(size.get(i))){
                        H.setMap(v, BMatrix.identity(naturalTransformation.get(k).getMap(v).rows), i);
                    }else{
                        H.setMap(v, inclusionMap.mult(naturalTransformation_inverse.get(k).getMap(v)), i);
                    }
                    if(v.get(i) > 0){
                        IntTuple w = v.minus(IntTuple.getStandardBasisElement(v.length(), i));
                        BMatrix hmap = H.getMap(w, i);
                        BMatrix natmap = naturalTransformation.get(k).getMap(v);
                        H.setMap(w, natmap.mult(hmap), i);
                    }
                }
            }
        }
        log.debug("Finished applying basis change.");

        log.debug("Starting verification...");
        for(int k=0;k<maxDimension;k++){
            Functor.verify(homfunctors.get(k));
            log.debug("Basis change dimension "+k+" OK.");
            //Nat.verify(chainFunctors.get(k), homfunctors.get(k), naturalTransformation.get(k));
            log.debug("Basis change map dimension "+k+" OK.");
            Pair<Functor, Nat> gnat = homfunctors.get(k).getSubFunctor(homologyDimension.get(k));
            Functor.verify(gnat._1());
            log.debug("Homology functor dimension "+k+" OK.");
            homfunctors.set(k, gnat._1());
        }
        log.debug("Finished verification.");

        log.debug("Finished computing homology functors.");
        return homfunctors;
    }
}

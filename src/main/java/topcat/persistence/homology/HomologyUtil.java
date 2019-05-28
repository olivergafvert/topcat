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
     * Computes the homology functors for each dimension less than 'maxdimension'.
     * @param simplexStorageStructure
     * @param size
     * @param maxDimension
     * @throws WrongDimensionException
     * @throws NoSolutionException
     */
    public static List<Functor> computeHomologyFunctors(final SimplexStorageStructure simplexStorageStructure, final IntTuple size, final int maxDimension) throws MalformedFunctorException, NoSolutionException{
        log.debug("Starting to compute homology functors...");

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
        List<HomologyWorker> workers = new ArrayList<>();
        for(IntTuple v : GridIterator.getSequence(size)){
            HomologyWorker worker = new HomologyWorker(simplexStorageStructure, v, maxDimension);
            workers.add(worker);
            futures.add(exec.submit(worker));
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

        for(HomologyWorker worker : workers){
            for(int k=0;k<maxDimension;k++){
                homologyDimension.get(k).set(worker.v, worker.homologyDimension[k]);
                naturalTransformation.get(k).setMap(worker.v, worker.naturalTransformation[k]);
                naturalTransformation_inverse.get(k).setMap(worker.v, worker.naturalTransformation_inverse[k]);
            }
        }

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

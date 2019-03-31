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

import it.unimi.dsi.fastutil.longs.Long2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.objects.Object2IntOpenHashMap;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import topcat.matrix.Column;
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
import topcat.util.paralelliterator.ParalellIterator;

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
     * Finds indices where the chain complex does not change so we don't have to compute it's homology
     * @return
     */
    private static int findInvariantDimension(IntTuple v, List<Grid<List<GradedColumn<Simplex>>>> columns_grid){
        List<IntTuple> basis = IntTuple.getStandardBasisSequence(v.length());
        for(int i=0;i<basis.size();i++){
            boolean is_invariant = true;
            for(int j=0;j<columns_grid.size();j++){
                if(columns_grid.get(j).get(v) != null || columns_grid.get(j).get(v.minus(basis.get(i)))!= null){
                    is_invariant = false;
                }
            }
            if(is_invariant) return i;
        }
        return -1;
    }

    private static  Simplex pop_local_pivot(PriorityQueue<Simplex> column){
        if(column.isEmpty()) return null;
        else{
            Simplex pivot = column.poll();
            IntTuple filtration_value = pivot.getValue();
            while(!column.isEmpty() && column.peek().equals(pivot)){
                column.poll();
                if(column.isEmpty()) return null;
                else{
                    pivot = column.poll();
                }
            }
            if(filtration_value.equals(pivot.getValue()))
                return pivot;
            return null;
        }
    }


    public static Simplex get_local_pivot(PriorityQueue<Simplex> column){
        Simplex result = pop_local_pivot(column);
        if(result != null) column.add(result);
        return result;
    }


    static void compute_local_pairs_h(List<Simplex> columns_to_reduce, SimplexStorageStructure simplexStorageStructure, IntTuple filtrationValue){
        Object2IntOpenHashMap<Simplex> pivot_column_index = new Object2IntOpenHashMap<>();
        for(int index_column_to_reduce = 0; index_column_to_reduce<columns_to_reduce.size();index_column_to_reduce++) {
            if(columns_to_reduce.get(index_column_to_reduce).local == 1)
                continue;
            PriorityQueue<Simplex> working_boundary= new PriorityQueue<Simplex>(10, new Comparator<Simplex>() {
                @Override
                public int compare(Simplex o1, Simplex o2) {
                    if(o1.getValue().lt(o2.getValue())) return 1;
                    if(o2.getValue().lt(o1.getValue())) return -1;
                    if(o1.getIndex() < o2.getIndex()) return 1;
                    if(o1.getIndex() > o2.getIndex()) return -1;
                    return 0;
                }
            });
            int index_column_to_add = index_column_to_reduce;
            while(true) {
                working_boundary.addAll(columns_to_reduce.get(index_column_to_add).getBoundary(simplexStorageStructure));
                Simplex pivot = get_local_pivot(working_boundary);
                if (pivot != null && pivot.getValue().equals(filtrationValue)) {
                    if (pivot_column_index.containsKey(pivot)) {
                        index_column_to_add = pivot_column_index.getInt(pivot);
                    } else {
                        pivot_column_index.addTo(pivot, index_column_to_reduce);
                        pivot.local = 1;
                        columns_to_reduce.get(index_column_to_reduce).local = -1;
                        pivot.pair = columns_to_reduce.get(index_column_to_reduce);
                        break;
                    }
                }else{
                    break;
                }
            }
        }
    }

    public static List<Grid<List<GradedColumn<Simplex>>>> chunk_reduction(final SimplexStorageStructure simplexStorageStructure, IntTuple size, final int maxDimension){

        log.debug("Starting to compute local pairs");
        //Compute local pairs
        new ParalellIterator<IntTuple, Void>(GridIterator.getSequence(size)){
                @Override
                public Void method(IntTuple index) {
                    for(int dim = 0; dim<maxDimension;dim++){
                        List<Simplex> simplices = simplexStorageStructure.getSimplicesAt(maxDimension-dim, index);
                        Collections.sort(simplices);
                        compute_local_pairs_h(simplices, simplexStorageStructure, index);
                    }
                    return null;
            }
        }.run();

        log.debug("Finished local pair computation.");
        log.debug("Starting global column red...");

        //Compute global columns
        List<Pair<IntTuple, List<List<GradedColumn<Simplex>>>>> global_columns_grid = new ParalellIterator<IntTuple, List<List<GradedColumn<Simplex>>>>(GridIterator.getSequence(size)){
            @Override
            public List<List<GradedColumn<Simplex>>> method(IntTuple index) {
                List<List<GradedColumn<Simplex>>> g_columns_v = new ArrayList<>();
                for(int dim=maxDimension;dim>0;dim--) {
                    List<Simplex> simplices = simplexStorageStructure.getSimplicesAt(dim, index);
                    Collections.sort(simplices);
                    List<GradedColumn<Simplex>> global_columns = compute_global_columns(simplices, simplexStorageStructure);
                    g_columns_v.add(global_columns);
                }
                Collections.reverse(g_columns_v);
                return g_columns_v;
            }
        }.run();

        List<Grid<List<GradedColumn<Simplex>>>> global_columns = new ArrayList<>();
        for(int i=0;i<maxDimension;i++){
            global_columns.add(Grid.create(size));
        }
        for(Pair<IntTuple, List<List<GradedColumn<Simplex>>>> g_columns : global_columns_grid){
            for(int i=0;i<maxDimension;i++){
                global_columns.get(i).set(g_columns._1(), g_columns._2().get(i));
            }
        }

        return global_columns;
    }

    public static List<GradedColumn<Simplex>> compute_global_columns(List<Simplex> columns_to_reduce, SimplexStorageStructure simplexStorageStructure){
        if(columns_to_reduce == null) return new ArrayList<>();
        List<GradedColumn<Simplex>> global_columns = new ArrayList<>();
        Long2ObjectOpenHashMap<List<Simplex>> cache = new Long2ObjectOpenHashMap<>();
        for(int index_column_to_reduce = 0; index_column_to_reduce<columns_to_reduce.size();index_column_to_reduce++) {
            if(columns_to_reduce.get(index_column_to_reduce).local != 0)
                continue;
            Column<Simplex> working_boundary= new Column<Simplex>();
            GradedColumn<Simplex> global_column= new GradedColumn<Simplex>(columns_to_reduce.get(index_column_to_reduce));
            working_boundary.addAll(columns_to_reduce.get(index_column_to_reduce).getBoundary(simplexStorageStructure));
            while(true) {
                Simplex pivot = working_boundary.get_pivot();
                if (pivot != null) {
                    if (pivot.local == -1) {
                        working_boundary.pop_pivot();
                    }else if(pivot.local == 1) {
                        if(cache.containsKey(pivot.getIndex()))
                            working_boundary.addAll(cache.get(pivot.getIndex()));
                        else {
                            List<Simplex> boundary = pivot.pair.getBoundary(simplexStorageStructure);
                            cache.put(pivot.getIndex(), boundary);
                            working_boundary.addAll(boundary);
                        }

                    }else{
                        global_column.add(working_boundary.pop_pivot());
                    }
                }else{
                    break;
                }
            }
            if(global_column.size() > 0)
                global_columns.add(global_column);
        }
        return global_columns;
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
        List<Pair<IntTuple, IntTuple>> invariantIndices = new ArrayList<>();

        //Chunk reduction
        List<Grid<List<GradedColumn<Simplex>>>> columns_grid = chunk_reduction(simplexStorageStructure, size, maxDimension);
        List<Grid<List<Simplex>>> basis_grid = new ArrayList<>();
        for(int i=0;i<maxDimension;i++){
            basis_grid.add(Grid.create(size));
        }

        HashMap<IntTuple, HomologyWorker> workerMap = new HashMap<>();

        //Compute homology basis and basis change maps
        for(IntTuple v : GridIterator.getSequence(size)){
            int invariantDim = findInvariantDimension(v, columns_grid);
            if(invariantDim == -1) {
                HomologyWorker worker = new HomologyWorker(simplexStorageStructure, columns_grid, v, maxDimension);
                workerMap.put(v, worker);
                workers.add(worker);
                futures.add(exec.submit(worker));
            }else{
                log.debug("Invariant index: "+v);
                invariantIndices.add(new Pair<>(v, v.minus(IntTuple.getStandardBasisElement(v.length(), invariantDim))));
            }
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
//        for(Pair<IntTuple, IntTuple> p : invariantIndices){
//            for(int k=0;k<maxDimension;k++){
//                homologyDimension.get(k).set(p._1(), homologyDimension.get(k).get(p._2()));
//                naturalTransformation.get(k).setMap(p._1(), new BMatrix(naturalTransformation.get(k).getMap(p._2())));
//                naturalTransformation_inverse.get(k).setMap(p._1(), new BMatrix(naturalTransformation_inverse.get(k).getMap(p._2())));
//            }
//        }
        log.debug("Finished computing basis change.");


        List<Functor> chainFunctors = new ArrayList<>();

        log.debug("Starting to apply basis change...");
        for(int k=0;k<homfunctors.size();k++){
            Functor H = homfunctors.get(k);
            Functor C = new Functor(size);
            for(IntTuple v : GridIterator.getSequence(size)){
                List<Simplex> currentSimplices = workerMap.get(v).basis.get(k);
                for(int i=0;i<v.length();i++){
                    BMatrix inclusionMap;
                    if(workerMap.containsKey(v.plus(IntTuple.getStandardBasisElement(v.length(), i))))
                        inclusionMap = computeInclusionMap(currentSimplices, workerMap.get(v.plus(IntTuple.getStandardBasisElement(v.length(), i))).basis.get(k));
                    else
                        inclusionMap = BMatrix.identity(currentSimplices.size());
                    C.setMap(v, inclusionMap, i);
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
            chainFunctors.add(C);
        }
        log.debug("Finished applying basis change.");

        log.debug("Starting verification...");
        for(int k=0;k<maxDimension;k++){
            Functor.verify(homfunctors.get(k));
            log.debug("Basis change dimension "+k+" OK.");
            Nat.verify(chainFunctors.get(k), homfunctors.get(k), naturalTransformation.get(k));
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

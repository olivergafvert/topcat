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
import topcat.matrix.BVector;
import topcat.matrix.Column;
import topcat.matrix.GradedColumn;
import topcat.persistence.functor.Functor;
import topcat.persistence.functor.Nat;
import topcat.persistence.functor.SparseFunctor;
import topcat.persistence.functor.exception.MalformedFunctorException;
import topcat.matrix.BMatrix;
import topcat.matrix.exception.NoSolutionException;
import topcat.matrix.exception.WrongDimensionException;
import topcat.persistence.simplex.Simplex;
import topcat.persistence.simplex.SimplexStorageStructure;
import topcat.persistence.simplex.SparseSimplexStorageStructure;
import topcat.util.*;
import topcat.util.paralelliterator.ParalellIntIterator;
import topcat.util.paralelliterator.ParalellIterator;

import java.util.*;

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

        log.debug("Starting to compute local pairs...");
        //Compute local pairs
        for(int dim = 0; dim<maxDimension;dim++) {
            final int dimension = dim;
            new ParalellIterator<IntTuple, Void>(simplexStorageStructure.gridSequence(dim)) {
                @Override
                public Void method(IntTuple index) {
                    List<Simplex> simplices = simplexStorageStructure.getSimplicesAt(maxDimension - dimension, index);
                    Collections.sort(simplices);
                    compute_local_pairs_h(simplices, simplexStorageStructure, index);
                    return null;
                }
            }.run();
        }

        log.debug("Finished local pair computation.");
        log.debug("Starting global column reduction...");

        List<List<Pair<IntTuple, List<GradedColumn<Simplex>>>>> global_columns_grid = new ArrayList<>();
        //Compute global columns
        for(int dim=maxDimension;dim>0;dim--) {
            final int dimension = dim;
            global_columns_grid.add(new ParalellIterator<IntTuple, List<GradedColumn<Simplex>>>(simplexStorageStructure.gridSequence(dimension)) {
                @Override
                public List<GradedColumn<Simplex>> method(IntTuple index) {
                    List<Simplex> simplices = simplexStorageStructure.getSimplicesAt(dimension, index);
                    Collections.sort(simplices);
                    List<GradedColumn<Simplex>> global_columns = compute_global_columns(simplices, simplexStorageStructure);
                    return global_columns;
                }
            }.run());
        }
        log.debug("Finished computing global columns.");

        List<Grid<List<GradedColumn<Simplex>>>> global_columns = new ArrayList<>();

        for(int i=0;i<global_columns_grid.size();i++){
            global_columns.add(simplexStorageStructure.getClass().equals(SparseSimplexStorageStructure.class) ? new SparseGrid<>(size) : Grid.create(size));
            for(int j=0;j<global_columns_grid.get(i).size();j++){
                global_columns.get(i).set(global_columns_grid.get(i).get(j)._1(), global_columns_grid.get(i).get(j)._2());
            }
        }

        Collections.reverse(global_columns);
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
            global_columns.add(global_column);
        }
        return global_columns;
    }


    public static Pair<List<GradedColumn<Simplex>>, List<GradedColumn<Simplex>>> reduce_matrix(List<GradedColumn<Simplex>> columns_to_reduce, HashMap<Simplex, GradedColumn<Simplex>> pivot_column_index){
        if(columns_to_reduce == null) return new Pair<>(new ArrayList<>(), new ArrayList<>());
        List<GradedColumn<Simplex>> kernel = new ArrayList<>();
        List<GradedColumn<Simplex>> image = new ArrayList<>();
        for(int index_column_to_reduce = 0; index_column_to_reduce<columns_to_reduce.size();index_column_to_reduce++) {
            Column<Simplex> working_boundary= new Column<>();
            GradedColumn<Simplex> working_reduction = new GradedColumn<>(columns_to_reduce.get(index_column_to_reduce).getGrade());
            working_reduction.add(columns_to_reduce.get(index_column_to_reduce).getGrade());
            working_boundary.addAll(columns_to_reduce.get(index_column_to_reduce));
            Simplex pivot = working_boundary.get_pivot();
            while(true) {
                if (pivot != null) {
                    if (pivot_column_index.containsKey(pivot)) {
                        working_reduction.add(pivot_column_index.get(pivot).getGrade());
                        working_boundary.addAll(pivot_column_index.get(pivot));
                        pivot = working_boundary.get_pivot();
                    } else {
                        pivot_column_index.put(pivot, columns_to_reduce.get(index_column_to_reduce));
                        image.add(columns_to_reduce.get(index_column_to_reduce));
                        break;
                    }
                }else{
                    kernel.add(working_reduction);
                    break;
                }
            }
        }
        return new Pair<>(kernel, image);
    }

    public static List<GradedColumn<Simplex>> basisChange(List<GradedColumn<Simplex>> columns_to_reduce, HashMap<Simplex, GradedColumn<Simplex>> pivot_column_index){
        if(columns_to_reduce == null || columns_to_reduce.size() == 0) return new ArrayList<>();

        List<Pair<Integer, GradedColumn<Simplex>>> kernel = new ParalellIntIterator<GradedColumn<Simplex>>(0, columns_to_reduce.size()) {
            @Override
            public GradedColumn<Simplex> method(Integer index_column_to_reduce) {
                Column<Simplex> working_boundary= new Column<>();
                GradedColumn<Simplex> working_reduction = new GradedColumn<>(columns_to_reduce.get(index_column_to_reduce).getGrade());
                working_boundary.addAll(columns_to_reduce.get(index_column_to_reduce));
                Simplex pivot = working_boundary.get_pivot();
                while(true) {
                    if (pivot != null) {
                        if (!pivot_column_index.containsKey(pivot)) throw new AssertionError();
                        working_reduction.add(pivot_column_index.get(pivot).getGrade());
                        working_boundary.addAll(pivot_column_index.get(pivot));
                        pivot = working_boundary.get_pivot();
                    }else{
                        return working_reduction;
                    }
                }
            }
        }.run();

        Collections.sort(kernel, new Comparator<Pair<Integer, GradedColumn<Simplex>>>() {
            @Override
            public int compare(Pair<Integer, GradedColumn<Simplex>> o1, Pair<Integer, GradedColumn<Simplex>> o2) {
                return o1._1().compareTo(o2._1());
            }
        });
        List<GradedColumn<Simplex>> ret = new ArrayList<>(columns_to_reduce.size());
        for(Pair<Integer, GradedColumn<Simplex>> pair : kernel)
            ret.add(pair._2());
        return ret;
    }

    public static HashMap<Simplex, GradedColumn<Simplex>> mergePivots(List<HashMap<Simplex, GradedColumn<Simplex>>> pivot_cols){
        if(pivot_cols == null || pivot_cols.size() == 0) return new HashMap<>();
        HashMap<Simplex, GradedColumn<Simplex>> main_map = pivot_cols.get(0);
        int k = 0;
        for(int i=1;i<pivot_cols.size();i++){
            if(pivot_cols.get(i).size()>main_map.size()){
                main_map = pivot_cols.get(i);
                k = i;
            }
        }
        HashMap<Simplex, GradedColumn<Simplex>> pivot_col = new HashMap<>();
        pivot_col.putAll(main_map);
        List<GradedColumn<Simplex>> columns = new ArrayList<>();
        for(int i=0;i<pivot_cols.size();i++){
            if(i!=k)
                columns.addAll(pivot_cols.get(i).values());
        }

        Collections.sort(columns);
        pivots(columns, pivot_col);
        return pivot_col;
    }

    public static List<Pair<Simplex, Integer>> pivots(List<GradedColumn<Simplex>> columns_to_reduce, HashMap<Simplex, GradedColumn<Simplex>> pivot_column_index){
        if(columns_to_reduce == null) return null;
        List<Pair<Simplex, Integer>> pivot_columns = new ArrayList<>();
        for(int index_column_to_reduce = 0; index_column_to_reduce<columns_to_reduce.size();index_column_to_reduce++) {
            Column<Simplex> working_boundary= new Column<>();
            working_boundary.addAll(columns_to_reduce.get(index_column_to_reduce));
            Simplex pivot = working_boundary.get_pivot();
            while(pivot != null) {
                if (pivot_column_index.containsKey(pivot)) {
                    working_boundary.addAll(pivot_column_index.get(pivot));
                    pivot = working_boundary.get_pivot();
                } else {
                    pivot_column_index.put(pivot, columns_to_reduce.get(index_column_to_reduce));
                    pivot_columns.add(new Pair<>(pivot, index_column_to_reduce));
                    break;
                }
            }
        }
        return pivot_columns;
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

        //The persistence modules, i.e homology of the multifiltration
        List<Functor> homfunctors = new ArrayList<>();

        //Chunk reduction
        List<Grid<List<GradedColumn<Simplex>>>> columns_grid = chunk_reduction(simplexStorageStructure, size, maxDimension);

        //Filter out the local dim 0 simplices
        Grid<List<GradedColumn<Simplex>>> basis_0_grid = new SparseGrid<>(size);
        List<Pair<IntTuple, List<GradedColumn<Simplex>>>> basis_res = new ParalellIterator<IntTuple, List<GradedColumn<Simplex>>>(simplexStorageStructure.gridSequence(0)) {
            @Override
            public List<GradedColumn<Simplex>> method(IntTuple index) {
                List<GradedColumn<Simplex>> basis = new ArrayList<>();
                for(Simplex s : simplexStorageStructure.getSimplicesAt(0, index)) {
                    if (s.local == 0) {
                        GradedColumn<Simplex> column = new GradedColumn<>(s);
                        column.add(s);
                        basis.add(column);
                    }
                }
                return basis;
            }
        }.run();
        for(Pair<IntTuple, List<GradedColumn<Simplex>>> pair : basis_res)
            basis_0_grid.set(pair._1(), pair._2());

        //Compute the support of image and kernel bases
        if(simplexStorageStructure.getClass().equals(SparseSimplexStorageStructure.class))
            ((SparseSimplexStorageStructure)simplexStorageStructure).preparePoset();

        List<IntTuple> indices = simplexStorageStructure.gridSequence();

        //Clear simplex storage structure
        simplexStorageStructure.clearSimplices();

        //Group supported cross-diagonal elements.
        HashMap<Integer, List<IntTuple>> diagonal_sequence_map = new HashMap<>();
        for(int i=0;i<indices.size();i++){
            int sum = indices.get(i).sum();
            if(!diagonal_sequence_map.containsKey(sum)) diagonal_sequence_map.put(sum, new ArrayList<>());
            diagonal_sequence_map.get(sum).add(indices.get(i));
        }
        List<List<IntTuple>> diagonal_sequence = new ArrayList<>();
        List<Integer> diagonal_keyset = new ArrayList<>(diagonal_sequence_map.keySet());
        Collections.sort(diagonal_keyset);
        for(int i=0;i<diagonal_keyset.size();i++)
            diagonal_sequence.add(diagonal_sequence_map.get(diagonal_keyset.get(i)));

        //Sequential reduction
        final List<Grid<HashMap<Simplex, GradedColumn<Simplex>>>> pivot_to_column_grid = new ArrayList<>();
        final List<Grid<List<GradedColumn<Simplex>>>> kernel_grid = new ArrayList<>();
        final List<Grid<HashSet<GradedColumn<Simplex>>>> image_grid = new ArrayList<>();

        for(int k=0;k<maxDimension;k++){
            if(simplexStorageStructure.getClass().equals(SparseSimplexStorageStructure.class)) {
                pivot_to_column_grid.add(new SparseGrid<>(size));
                kernel_grid.add(new SparseGrid<>(size));
                image_grid.add(new SparseGrid<>(size));
            }else{
                pivot_to_column_grid.add(Grid.create(size));
                kernel_grid.add(Grid.create(size));
                image_grid.add(Grid.create(size));
            }
        }

        class Packet{
            List<HashMap<Simplex, GradedColumn<Simplex>>> pivots = new ArrayList<>();
            List<List<GradedColumn<Simplex>>> kernels = new ArrayList<>();
            List<HashSet<GradedColumn<Simplex>>> images = new ArrayList<>();
            List<List<Simplex>> basis = new ArrayList<>();
        }

        log.debug("Starting boundary matrix reduction...");
        for(int d=0;d<diagonal_sequence.size();d++) {
            List<Pair<IntTuple, Packet>> results = new ParalellIterator<IntTuple, Packet>(diagonal_sequence.get(d)) {
                @Override
                public Packet method(IntTuple v) {
                    Packet packet = new Packet();
                    for (int i = 0; i < maxDimension; i++) {
                        //Fetch global columns and sort
                        List<GradedColumn<Simplex>> columns = columns_grid.get(i).get(v) == null ? new ArrayList<>() : new ArrayList<>(columns_grid.get(i).get(v));
                        Collections.sort(columns);

                        //Record the grades of the global columns. This is the basis of the reduced chain complex at this index
                        if (i < maxDimension - 1) {
                            packet.basis.add(new ArrayList<>());
                            for (int j = 0; j < columns.size(); j++) {
                                packet.basis.get(i).add(columns.get(j).getGrade());
                            }
                        }

                        //Prepare the kernel and image of previous index
                        List<HashMap<Simplex, GradedColumn<Simplex>>> prev_pivots = new ArrayList<>(); //The pivots of the prev neighbours
                        List<IntTuple> prev_v = simplexStorageStructure.getAdjacent(v); //The prev neighbours
                        HashSet<GradedColumn<Simplex>> image = new HashSet<>();
                        for (IntTuple w : prev_v) {
                            prev_pivots.add(pivot_to_column_grid.get(i).get(w));
                            image.addAll(image_grid.get(i).get(w));
                        }

                        HashMap<Simplex, GradedColumn<Simplex>> pivot_to_columns = new HashMap<>(); //The pivots we'll use for reduction
                        int k = 0;
                        List<GradedColumn<Simplex>> non_red_image = new ArrayList<>();
                        if (prev_pivots.size() > 0) { //If we have more than one neighbour we need to merge pivots
                            //Choose the biggest pivot set as reference
                            HashMap<Simplex, GradedColumn<Simplex>> main_map = prev_pivots.get(0);
                            for (int j = 1; j < prev_pivots.size(); j++) {
                                if (prev_pivots.get(j).size() > main_map.size()) {
                                    main_map = prev_pivots.get(j);
                                    k = j;
                                }
                            }
                            HashSet<GradedColumn<Simplex>> reduced = new HashSet<>();
                            reduced.addAll(main_map.values());

                            for (int j = 0; j < prev_pivots.size(); j++) {
                                if (j != k) {
                                    for (GradedColumn<Simplex> column : prev_pivots.get(j).values()) {
                                        if (!reduced.contains(column)) {
                                            reduced.add(column);
                                            non_red_image.add(column);
                                        }
                                    }
                                }
                            }
                            pivot_to_columns.putAll(main_map);
                        }

                        Pair<List<GradedColumn<Simplex>>, List<GradedColumn<Simplex>>> off_kernel = reduce_matrix(non_red_image, pivot_to_columns);
                        Pair<List<GradedColumn<Simplex>>, List<GradedColumn<Simplex>>> kerim = reduce_matrix(columns, pivot_to_columns);

                        //Add off kernel elements to kernel set
                        kerim._1().addAll(off_kernel._1());

                        //Add previous kernel elements to kernel
                        HashSet<GradedColumn<Simplex>> p_kernel_cols = new HashSet<>();
                        p_kernel_cols.addAll(kerim._1());
                        for (int j = 0; j < prev_v.size(); j++) {
                            for (GradedColumn<Simplex> column : kernel_grid.get(i).get(prev_v.get(j))) {
                                if (!p_kernel_cols.contains(column)) {
                                    p_kernel_cols.add(column);
                                    kerim._1().add(column);
                                }
                            }
                        }

                        //If i>0 we can express the image in the kernel basis
                        if (i > 0) {
                            HashMap<Simplex, GradedColumn<Simplex>> kernel_pivs = new HashMap<>();
                            List<GradedColumn<Simplex>> kernel = packet.kernels.get(packet.kernels.size() - 1);
                            pivots(kernel, kernel_pivs);
                            List<GradedColumn<Simplex>> image_s = basisChange(kerim._2(), kernel_pivs);
                            image.addAll(image_s);
                        } else {
                            image.addAll(kerim._2());
                        }
                        packet.pivots.add(pivot_to_columns);
                        packet.kernels.add(kerim._1());
                        packet.images.add(image);
                    }
                    return packet;
                }
            }.run();
            for(Pair<IntTuple, Packet> pair : results){
                for(int k=0;k<pair._2().pivots.size();k++) {
                    pivot_to_column_grid.get(k).set(pair._1(), pair._2().pivots.get(k));
                    kernel_grid.get(k).set(pair._1(), pair._2().kernels.get(k));
                    image_grid.get(k).set(pair._1(), pair._2().images.get(k));
                }
            }
        }


        log.debug("Finished reduction.");

        log.debug("Starting homology basis computation...");

        List<Grid<HashMap<Simplex, GradedColumn<Simplex>>>> kernel_pivots_grid = new ArrayList<>();
        List<Grid<List<GradedColumn<Simplex>>>> homology_grid = new ArrayList<>();
        List<Grid<Integer>> homology_dim_grid = new ArrayList<>();

        for(int k=0;k<maxDimension;k++){
            if(simplexStorageStructure.getClass().equals(SparseSimplexStorageStructure.class)){
                homology_grid.add(new SparseGrid<>(size));
                homology_dim_grid.add(new SparseGrid<>(size));
                kernel_pivots_grid.add(new SparseGrid<>(size));
            }else{
                homology_grid.add(Grid.create(size));
                homology_dim_grid.add(Grid.create(size));
                kernel_pivots_grid.add(Grid.create(size));
            }
        }

        class HomologyPacket{
            List<Integer> homology_dim = new ArrayList<>();
            List<List<GradedColumn<Simplex>>> homologybases = new ArrayList<>();
            List<GradedColumn<Simplex>> basis_0 = new ArrayList<>();
            List<HashMap<Simplex, GradedColumn<Simplex>>> kernel_pivots = new ArrayList<>();
        }

        for(int d=0;d<diagonal_sequence.size();d++) {
            List<Pair<IntTuple, HomologyPacket>> results = new ParalellIterator<IntTuple, HomologyPacket>(diagonal_sequence.get(d)) {
                @Override
                public HomologyPacket method(IntTuple v) {
                    HomologyPacket packet = new HomologyPacket();
                    for (int i = 0; i < maxDimension; i++) {
                        List<GradedColumn<Simplex>> image = new ArrayList<>(image_grid.get(i).get(v));
                        Collections.sort(image);

                        List<IntTuple> prev_v = simplexStorageStructure.getAdjacent(v);

                        //Prepare kernel basis
                        List<GradedColumn<Simplex>> kernel_basis;
                        if (i == 0) {
                            HashSet<GradedColumn<Simplex>> kernel_hashset = new HashSet<>();
                            if (basis_0_grid.get(v) != null) kernel_hashset.addAll(basis_0_grid.get(v));
                            for (IntTuple w : prev_v) {
                                kernel_hashset.addAll(basis_0_grid.get(w));
                            }
                            kernel_basis = new ArrayList<>(kernel_hashset);
                            Collections.sort(kernel_basis);
                            packet.basis_0 = kernel_basis;
                        } else {
                            kernel_basis = kernel_grid.get(i - 1).get(v);
                        }
                        Object2IntOpenHashMap<Simplex> index_map = new Object2IntOpenHashMap<>();
                        int j = 0;
                        for (GradedColumn<Simplex> column : kernel_basis) {
                            index_map.put(column.getGrade(), j++);
                        }

                        List<Pair<Simplex, Integer>> image_pivots = pivots(image, new HashMap<>());
                        List<GradedColumn<Simplex>> red_image = new ArrayList<>();
                        for (Pair<Simplex, Integer> pivot_pair : image_pivots) {
                            GradedColumn<Simplex> image_column = image.get(pivot_pair._2());
                            GradedColumn<Simplex> copy = new GradedColumn<>(image_column.getGrade());
                            copy.addAll(image_column);
                            red_image.add(copy);
                        }

                        HashMap<Simplex, GradedColumn<Simplex>> kernel_pivots = new HashMap<>();
                        pivots(kernel_basis, kernel_pivots);

                        int image_size = red_image.size();

                        if (image_size > kernel_basis.size()) throw new AssertionError();

                        //Extend basis
                        for (j = 0; j < kernel_basis.size(); j++) {
                            GradedColumn<Simplex> column = new GradedColumn<>(new Simplex(Long.MAX_VALUE, i));
                            column.add(kernel_basis.get(j).getGrade());
                            red_image.add(column);
                        }

                        HashMap<Simplex, GradedColumn<Simplex>> pivs = new HashMap<>();
                        List<Pair<Simplex, Integer>> pivotlist = pivots(red_image, pivs);

                        List<GradedColumn<Simplex>> homologybasis = new ArrayList<>();

                        int homology_basis_size = pivotlist.size() - image_size;

                        HashMap<Simplex, GradedColumn<Simplex>> basis_pivots = new HashMap<>();
                        int k = 0;
                        for (j = image_size; j < pivotlist.size(); j++) {
                            GradedColumn<Simplex> column = new GradedColumn<>(new Simplex(k++, i));
                            GradedColumn<Simplex> image_column = red_image.get(pivotlist.get(j)._2());
                            while (!image_column.isEmpty()) {
                                column.addAll(kernel_basis.get(index_map.getInt(image_column.pop_pivot())));
                            }
                            basis_pivots.put(pivotlist.get(j)._1(), column);
                            homologybasis.add(column);
                        }
                        for (j = 0; j < image_size; j++) {
                            GradedColumn<Simplex> column = new GradedColumn<>(new Simplex(k++, i));
                            GradedColumn<Simplex> image_column = red_image.get(pivotlist.get(j)._2());
                            while (!image_column.isEmpty()) {
                                column.addAll(kernel_basis.get(index_map.getInt(image_column.pop_pivot())));
                            }
                            basis_pivots.put(pivotlist.get(j)._1(), column);
                            homologybasis.add(column);
                        }
                        packet.homology_dim.add(homology_basis_size);
                        packet.homologybases.add(homologybasis);
                        packet.kernel_pivots.add(basis_pivots);
                    }
                    return packet;
                }
            }.run();
            for(Pair<IntTuple, HomologyPacket> pair : results){
                for(int k=0;k<maxDimension;k++) {
                    homology_dim_grid.get(k).set(pair._1(), pair._2().homology_dim.get(k));
                    homology_grid.get(k).set(pair._1(), pair._2().homologybases.get(k));
                    kernel_pivots_grid.get(k).set(pair._1(), pair._2().kernel_pivots.get(k));
                }
                basis_0_grid.set(pair._1(), pair._2().basis_0);
            }
        }




        if(simplexStorageStructure.getClass().equals(SparseSimplexStorageStructure.class)) {
            HashSet<IntTuple> support = new HashSet<>(indices);
            for(int k=0;k<maxDimension;k++)
                homfunctors.add(new SparseFunctor(size, kernel_pivots_grid.get(k), homology_grid.get(k), homology_dim_grid.get(k), support, ((SparseSimplexStorageStructure)simplexStorageStructure).getAdjacentMap()));
        }
        else{
            for(int k=0;k<maxDimension;k++)
                homfunctors.add(new Functor(size));

            //Glue maps
            List<Pair<IntTuple, List<List<Pair<IntTuple, BMatrix>>>>> h_maps = new ParalellIterator<IntTuple, List<List<Pair<IntTuple, BMatrix>>>>(indices) {
                @Override
                public List<List<Pair<IntTuple, BMatrix>>> method(IntTuple index) {
                    List<List<Pair<IntTuple, BMatrix>>> maps = new ArrayList<>();
                    List<IntTuple> adjacent = simplexStorageStructure.getAdjacent(index);
                    for(int i=0;i<adjacent.size();i++){
                        maps.add(new ArrayList<>());
                    }
                    for(int i=0;i<maxDimension;i++) {
                        int homologyDimension_v = homology_dim_grid.get(i).get(index);
                        HashMap<Simplex, GradedColumn<Simplex>> pivot_to_column_index = kernel_pivots_grid.get(i).get(index);
                        for(int j=0;j<adjacent.size();j++) {
                            List<GradedColumn<Simplex>> homologyBasis = homology_grid.get(i).get(adjacent.get(j)).subList(0, homology_dim_grid.get(i).get(adjacent.get(j)));
                            List<GradedColumn<Simplex>> hom_trans = basisChange(homologyBasis, pivot_to_column_index);
                            BMatrix M = new BMatrix(homologyDimension_v, homologyBasis.size());
                            for(int k=0;k<hom_trans.size();k++){
                                GradedColumn<Simplex> column = hom_trans.get(k);
                                while(!column.isEmpty()){
                                    Simplex grade = column.pop_pivot();
                                    if(grade.getIndex()<homologyDimension_v){
                                        M.set((int)grade.getIndex(), k, true);
                                    }
                                }
                            }
                            maps.get(j).add(new Pair<>(adjacent.get(j), M));
                        }
                    }
                    return maps;
                }
            }.run();

            for(Pair<IntTuple, List<List<Pair<IntTuple, BMatrix>>>> pair : h_maps){
                for(int i=0;i<pair._2().size();i++){
                    List<Pair<IntTuple, BMatrix>> maps = pair._2().get(i);
                    for(int j=0;j<maxDimension;j++){
                        homfunctors.get(j).setMap(maps.get(j)._1(), pair._1(), maps.get(j)._2());
                    }
                }
            }
        }

        log.debug("Finished applying basis change.");

        log.debug("Finished computing homology functors.");
        return homfunctors;
    }
}

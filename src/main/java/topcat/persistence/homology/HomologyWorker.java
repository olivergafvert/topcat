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

package topcat.persistence.homology;

import it.unimi.dsi.fastutil.ints.Int2IntOpenHashMap;
import it.unimi.dsi.fastutil.longs.Long2IntOpenHashMap;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import topcat.matrix.BMatrix;
import topcat.matrix.BVector;
import topcat.matrix.Column;
import topcat.matrix.GradedColumn;
import topcat.persistence.simplex.Simplex;
import topcat.persistence.simplex.SimplexStorageStructure;
import topcat.util.*;

import java.util.*;

/**
 * Worker node in the parallelized computation of the homology.
 */
public class HomologyWorker implements Runnable{
        public static Logger log = LoggerFactory.getLogger(HomologyWorker.class);
        int maxDimension;
        IntTuple v;
        List<List<Simplex>> basis = new ArrayList<>();
        List<Grid<List<GradedColumn<Simplex>>>> global_columns_grid;
        int[] homologyDimension;
        BMatrix[] naturalTransformation;
        BMatrix[] naturalTransformation_inverse;
        SimplexStorageStructure simplexStorageStructure;
        BinomialCoeffTable binomialCoeffTable;

    public HomologyWorker(SimplexStorageStructure simplexStorageStructure, List<Grid<List<GradedColumn<Simplex>>>> global_columns_grid, IntTuple v, int maxDimension) {
        this.maxDimension = maxDimension;
        this.v = v;
        this.homologyDimension = new int[maxDimension];
        this.naturalTransformation = new BMatrix[maxDimension];
        this.naturalTransformation_inverse = new BMatrix[maxDimension];
        this.simplexStorageStructure = simplexStorageStructure;
        if(simplexStorageStructure!=null)
            this.binomialCoeffTable = new BinomialCoeffTable(simplexStorageStructure.getNumberOfVertices(), maxDimension+2);
        this.global_columns_grid = global_columns_grid;
    }


    public List<GradedColumn<Simplex>> reduce_matrix(List<GradedColumn<Simplex>> columns_to_reduce, HashMap<Simplex, GradedColumn<Simplex>> pivot_column_index){
        List<GradedColumn<Simplex>> reductionMatrix = new ArrayList<>();
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
                        break;
                    }
                }else{
                    reductionMatrix.add(working_reduction);
                    break;
                }
            }
        }
        return reductionMatrix;
    }

    /**
     * Computes a set of pivots of the matrix given by the list of columns 'columns_to_reduce'.
     * @param columns_to_reduce - a list of columns of the matrix.
     * @param pivot_column_index
     * @return A list of pivot and column indices.
     */
        public List<Pair<Simplex, Integer>> computePivots(List<GradedColumn<Simplex>> columns_to_reduce, Long2IntOpenHashMap pivot_column_index){
            List<Pair<Simplex, Integer>> pivot_columns = new ArrayList<>();
            for(Integer index_column_to_reduce = 0; index_column_to_reduce<columns_to_reduce.size();index_column_to_reduce++) {
                Column<Simplex> working_column= new Column<>();
                working_column.addAll(columns_to_reduce.get(index_column_to_reduce));
                Simplex pivot;
                while((pivot=working_column.get_pivot()) != null) {
                    if(pivot_column_index.containsKey(pivot.getIndex())){
                        working_column.addAll(columns_to_reduce.get(pivot_column_index.get(pivot.getIndex())));
                    }else{
                        pivot_column_index.addTo(pivot.getIndex(), index_column_to_reduce);
                        pivot_columns.add(new Pair<>(pivot, index_column_to_reduce));
                        break;
                    }
                }
            }
            return pivot_columns;
        }

    public List<Pair<Integer, Integer>> computePivots2(List<Column<Integer>> columns_to_reduce, Int2IntOpenHashMap pivot_column_index){
        List<Pair<Integer, Integer>> pivot_columns = new ArrayList<>();
        for(Integer index_column_to_reduce = 0; index_column_to_reduce<columns_to_reduce.size();index_column_to_reduce++) {
            Column<Integer> working_column= new Column<>();
            working_column.addAll(columns_to_reduce.get(index_column_to_reduce));
            Integer pivot;
            while((pivot=working_column.get_pivot()) != null) {
                if(pivot_column_index.containsKey((int)pivot)){
                    working_column.addAll(columns_to_reduce.get(pivot_column_index.get((int)pivot)));
                }else{
                    pivot_column_index.addTo(pivot, index_column_to_reduce);
                    pivot_columns.add(new Pair<>(pivot, index_column_to_reduce));
                    break;
                }
            }
        }
        return pivot_columns;
    }


    public List<Pair<Long, Integer>> computePivots3(List<Column<Long>> columns_to_reduce, Long2IntOpenHashMap pivot_column_index){
        List<Pair<Long, Integer>> pivot_columns = new ArrayList<>();
        for(Integer index_column_to_reduce = 0; index_column_to_reduce<columns_to_reduce.size();index_column_to_reduce++) {
            Column<Long> working_column= new Column<>();
            working_column.addAll(columns_to_reduce.get(index_column_to_reduce));
            Long pivot;
            while((pivot=working_column.get_pivot()) != null) {
                if(pivot_column_index.containsKey((long)pivot)){
                    working_column.addAll(columns_to_reduce.get(pivot_column_index.get(pivot)));
                }else{
                    pivot_column_index.addTo(pivot, index_column_to_reduce);
                    pivot_columns.add(new Pair<Long, Integer>(pivot, index_column_to_reduce));
                    break;
                }
            }
        }
        return pivot_columns;
    }

    /**
     * Computes the matrix inverse of the matrix whose columns are given by the list 'columns_to_reduce'.
     * @param columns_to_reduce
     * @return the matrix inverse.
     */
    public BMatrix inverse(List<Column<Integer>> columns_to_reduce){
        Int2IntOpenHashMap pivot_column_index = new Int2IntOpenHashMap();
        computePivots2(columns_to_reduce, pivot_column_index);
        assert pivot_column_index.size() == columns_to_reduce.size();
        BMatrix reductionMatrix = new BMatrix(columns_to_reduce.size(), columns_to_reduce.size());//BMatrix.identity(columns_to_reduce.size());

        for(int i = 0; i<columns_to_reduce.size();i++) {
            Integer index_column_to_reduce = pivot_column_index.get(i);
            Column<Integer> working_column= new Column<>();
            working_column.addAll(columns_to_reduce.get(index_column_to_reduce));
            working_column.add(i);
            reductionMatrix.set(index_column_to_reduce, i, true);
            Integer pivot;
            while((pivot=working_column.get_pivot()) != null) {
                working_column.addAll(columns_to_reduce.get(pivot_column_index.get((int)pivot)));
                reductionMatrix.set(pivot_column_index.get((int)pivot), i, !reductionMatrix.get(pivot_column_index.get((int)pivot), i));
            }
        }
        return reductionMatrix;
    }

    /**
     * Computes a basis for the homology in each dimension of the chain complex given by 'chain'.
     *
     */
    protected void computeHomologyBasis() {
        List<List<GradedColumn<Simplex>>> boundaryMatrices = new ArrayList<>();
        for(int i=0;i<maxDimension;i++){
            boundaryMatrices.add(new ArrayList<>());
            basis.add(new ArrayList<>());
        }
        for(IntTuple w : GridIterator.getSequence(v)){
            boundaryMatrices.get(0).addAll(global_columns_grid.get(0).get(w));
            basis.get(0).addAll(simplexStorageStructure.getSimplicesAt(0, w));
            Collections.sort(basis.get(0));
            for(int i=1;i<maxDimension;i++) {
                List<GradedColumn<Simplex>> columns = global_columns_grid.get(i-1).get(w);
                boundaryMatrices.get(i).addAll(global_columns_grid.get(i).get(w));
                for(int j=0;j<columns.size();j++)
                    basis.get(i).add(columns.get(j).getGrade());
            }
        }

        List<GradedColumn<Simplex>> kernel = new ArrayList<>();

        for (int dim = 0; dim < maxDimension; dim++) {
            HashMap<Simplex, GradedColumn<Simplex>> pivot_column_index = new HashMap<>();
            Collections.sort(boundaryMatrices.get(dim));
            List<GradedColumn<Simplex>> kernel_p1 = reduce_matrix(boundaryMatrices.get(dim), pivot_column_index);
            List<GradedColumn<Simplex>> image_basis = new ArrayList<>(pivot_column_index.values());
            Collections.sort(image_basis);
            int image_basis_size = image_basis.size();
            if(dim>0) {
                image_basis.addAll(kernel);
            }else{
                for(Simplex s : basis.get(0)){
                    GradedColumn<Simplex> column = new GradedColumn<>(new Simplex(Long.MAX_VALUE, dim));
                    column.add(s);
                    image_basis.add(column);
                }
            }
            List<Pair<Simplex, Integer>> pivot_column = computePivots(image_basis, new Long2IntOpenHashMap());

            for(int i=0;i<image_basis_size;i++){
                assert pivot_column.get(i)._2() < image_basis_size;
            }

            List<Column<Long>> homology_basis = new ArrayList<>();
            for (int i = image_basis_size; i < pivot_column.size(); i++) {
                Column<Long> column = new Column<>();
                for(Simplex s : image_basis.get(pivot_column.get(i)._2()))
                    column.add(s.getIndex());
                homology_basis.add(column);
            }

            int homology_basis_size = homology_basis.size();
            for (int i = 0; i < basis.get(dim).size(); i++) {
                Column<Long> column = new Column<>();
                column.add(basis.get(dim).get(i).getIndex());
                homology_basis.add(column);
            }

            List<Pair<Long, Integer>> pivot_column2 = computePivots3(homology_basis, new Long2IntOpenHashMap());
            Long2IntOpenHashMap reverse_basis_map = new Long2IntOpenHashMap();
            for(int i=0;i<basis.get(dim).size();i++){
                reverse_basis_map.put(basis.get(dim).get(i).getIndex(), i);
            }
            List<BVector> extended_basis = new ArrayList<>();
            List<Column<Integer>> e_basis = new ArrayList<>();
            for (Pair<Long, Integer> pivot : pivot_column2) {
                BVector row = new BVector(basis.get(dim).size());
                Column<Integer> column = new Column<>();
                for (Long pos : homology_basis.get(pivot._2())) {
                    row.set(reverse_basis_map.get((long)pos), !row.get(reverse_basis_map.get((long)pos)));
                    column.add(reverse_basis_map.get((long)pos));
                }
                extended_basis.add(row);
                e_basis.add(column);
            }

            BMatrix extendedBasis_inv = inverse(e_basis);
            BMatrix extendedBasis = (new BMatrix(extended_basis)).transpose();

            homologyDimension[dim] = homology_basis_size;
            naturalTransformation[dim] = extendedBasis_inv;
            naturalTransformation_inverse[dim] = extendedBasis;

            kernel = new ArrayList<>(kernel_p1);
        }
    }

    @Override
    public void run() {
        log.debug("Starting basis change computation of position: "+v);
        computeHomologyBasis();
        log.debug("Finished basis change computation of position: "+v);
    }
}


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

import it.unimi.dsi.fastutil.ints.IntIterator;
import it.unimi.dsi.fastutil.longs.Long2IntOpenHashMap;
import it.unimi.dsi.fastutil.longs.LongOpenHashSet;
import it.unimi.dsi.fastutil.objects.Object2IntOpenHashMap;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import topcat.matrix.BMatrix;
import topcat.matrix.BVector;
import topcat.matrix.Column;
import topcat.matrix.GradedColumn;
import topcat.persistence.simplex.Simplex;
import topcat.persistence.simplex.SimplexCoboundaryEnumerator;
import topcat.persistence.simplex.SimplexStorageStructure;
import topcat.util.BinomialCoeffTable;
import topcat.util.IntTuple;
import topcat.util.Pair;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

/**
 * Worker node in the parallelized computation of the homology.
 */
public class HomologyWorker2 implements Runnable{
        public static Logger log = LoggerFactory.getLogger(HomologyWorker2.class);
        int maxDimension;
        IntTuple v;
        List<List<GradedColumn<Simplex>>> kernels = new ArrayList<>();
        List<List<GradedColumn<Simplex>>> images = new ArrayList<>();
        List<List<Simplex>> basis = new ArrayList<>();
        int[] homologyDimension;
        BMatrix[] naturalTransformation;
        BMatrix[] naturalTransformation_inverse;
        SimplexStorageStructure simplexStorageStructure;
        BinomialCoeffTable binomialCoeffTable;

    HomologyWorker2(SimplexStorageStructure simplexStorageStructure, IntTuple v, int maxDimension) {
        this.maxDimension = maxDimension;
        this.v = v;
        this.homologyDimension = new int[maxDimension];
        this.naturalTransformation = new BMatrix[maxDimension];
        this.naturalTransformation_inverse = new BMatrix[maxDimension];
        this.simplexStorageStructure = simplexStorageStructure;
        this.binomialCoeffTable = new BinomialCoeffTable(simplexStorageStructure.getNumberOfVertices(), maxDimension+2);
    }

    public List<GradedColumn<Simplex>> computeBoundaries(List<Simplex> simplices){
        if(simplices == null) return new ArrayList<>();
        List<GradedColumn<Simplex>> boundaries = new ArrayList<>();
        for(Simplex s : simplices){
            GradedColumn<Simplex> column = new GradedColumn<>(s);
            column.addAll(s.getBoundary(simplexStorageStructure));
            boundaries.add(column);
        }
        return boundaries;
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
     * Adds the coboundary of 's' to the current column 'working_columns' and returns the largest index where the column
     * is non-zero.
     * @param s
     * @param index_lookup
     * @param working_coboundary
     * @return
     */
    public Long add_coboundary_and_get_pivot(Simplex s, LongOpenHashSet index_lookup, Column<Long> working_coboundary){
        SimplexCoboundaryEnumerator enumerator = new SimplexCoboundaryEnumerator(s, simplexStorageStructure.getNumberOfVertices(), binomialCoeffTable);
        while (enumerator.hasNext()) {
            long index = enumerator.next();
            if (index_lookup.contains(index)) {
                working_coboundary.add(index);
            }
        }
        return working_coboundary.get_pivot();
    }

    /**
     * Performs an implicit matrix reduction of the coboundary matrix.
     * @param columns_to_reduce - a list of the non-zero columns of the coboundary matrix.
     * @param pivot_column_index
     * @param index_lookup - hash set containing the indexes of the non-zero rows in the coboundary matrix.
     * @return
     */
    public BMatrix reduce_coboundary_matrix(List<Simplex> columns_to_reduce, Long2IntOpenHashMap pivot_column_index, LongOpenHashSet index_lookup){
            BMatrix reductionMatrix = BMatrix.identity(columns_to_reduce.size());
            for(int index_column_to_reduce = 0; index_column_to_reduce<columns_to_reduce.size();index_column_to_reduce++) {
                Column<Long> working_coboundary= new Column<>();
                int index_column_to_add = index_column_to_reduce;
                while(true) {
                    Long pivot = add_coboundary_and_get_pivot(columns_to_reduce.get(index_column_to_add), index_lookup, working_coboundary);
                    if (pivot != null) {
                        if (pivot_column_index.containsKey((long)pivot)) {
                            index_column_to_add = pivot_column_index.get((long)pivot);
                            reductionMatrix.set(index_column_to_reduce, index_column_to_add, !reductionMatrix.get(index_column_to_reduce, index_column_to_add));
                        } else {
                            pivot_column_index.addTo(pivot, index_column_to_reduce);
                            break;
                        }
                    }else{
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


    /**
     * Computes the matrix inverse of the matrix whose columns are given by the list 'columns_to_reduce'.
     * @param columns_to_reduce
     * @return the matrix inverse.
     */
    public BMatrix inverse(List<GradedColumn<Simplex>> columns_to_reduce){
        Long2IntOpenHashMap pivot_column_index = new Long2IntOpenHashMap();
        computePivots(columns_to_reduce, pivot_column_index);
        assert pivot_column_index.size() == columns_to_reduce.size();
        BMatrix reductionMatrix = new BMatrix(columns_to_reduce.size(), columns_to_reduce.size());//BMatrix.identity(columns_to_reduce.size());

        long[] pivots = pivot_column_index.keySet().toLongArray();
        Arrays.sort(pivots);

        for(int i = 0; i<columns_to_reduce.size();i++) {
            int index_column_to_reduce = pivot_column_index.get(pivots[i]);
            Column<Simplex> working_column= new Column<>();
            working_column.addAll(columns_to_reduce.get(index_column_to_reduce));
            reductionMatrix.set(index_column_to_reduce, i, true);
            Simplex pivot;
            while((pivot=working_column.get_pivot()) != null) {
                if(pivot_column_index.get(pivot.getIndex()) == index_column_to_reduce){
                    working_column.pop_pivot();
                }
                else{
                    working_column.addAll(columns_to_reduce.get(pivot_column_index.get(pivot.getIndex())));
                    reductionMatrix.set(pivot_column_index.get(pivot.getIndex()), i, !reductionMatrix.get(pivot_column_index.get(pivot.getIndex()), i));
                }
            }
        }
        return reductionMatrix;
    }

    /**
     * Computes a basis for the homology in each dimension of the chain complex given by 'chain'.
     *
     */
    protected void computeHomologyBasis() {
        for (int dim = 0; dim < maxDimension; dim++) {
            List<GradedColumn<Simplex>> image_basis = images.get(dim);
            int image_basis_size = image_basis.size();
            if(dim>0) {
                image_basis.addAll(kernels.get(dim - 1));
            }else{
                for(Simplex s : basis.get(0)){
                    GradedColumn<Simplex> column = new GradedColumn<>(new Simplex(Long.MAX_VALUE, dim));
                    column.add(s);
                    image_basis.add(column);
                }
            }
            List<Pair<Simplex, Integer>> pivot_column = computePivots(image_basis, new Long2IntOpenHashMap());

            List<GradedColumn<Simplex>> homology_basis = new ArrayList<>();
            for (int i = image_basis_size; i < pivot_column.size(); i++) {
                homology_basis.add(image_basis.get(pivot_column.get(i)._2()));
            }

            int homology_basis_size = homology_basis.size();
            for (int i = 0; i < basis.get(dim).size(); i++) {
                GradedColumn<Simplex> column = new GradedColumn<>(new Simplex(Long.MAX_VALUE, dim)); //TODO: Set grade as null?
                column.add(basis.get(dim).get(i));
                homology_basis.add(column);
            }

            pivot_column = computePivots(homology_basis, new Long2IntOpenHashMap());
            Object2IntOpenHashMap<Simplex> reverse_basis_map = new Object2IntOpenHashMap<>();
            for(int i=0;i<basis.get(dim).size();i++){
                reverse_basis_map.put(basis.get(dim).get(i), i);
            }
            List<BVector> extended_basis = new ArrayList<>();
            List<GradedColumn<Simplex>> e_basis = new ArrayList<>();
            for (Pair<Simplex, Integer> pivot : pivot_column) {
                BVector row = new BVector(basis.get(dim).size());
                for (Simplex pos : homology_basis.get(pivot._2()))
                    row.set(reverse_basis_map.getInt(pos), true);
                extended_basis.add(row);
                e_basis.add(homology_basis.get(pivot._2()));
            }

            BMatrix extendedBasis_inv = inverse(e_basis);
            BMatrix extendedBasis = (new BMatrix(extended_basis)).transpose();

            homologyDimension[dim] = homology_basis_size;
            naturalTransformation[dim] = extendedBasis_inv;
            naturalTransformation_inverse[dim] = extendedBasis;
        }
    }

    @Override
    public void run() {
        log.debug("Starting basis change computation of position: "+v);
        computeHomologyBasis();
        log.debug("Finished basis change computation of position: "+v);
    }
}


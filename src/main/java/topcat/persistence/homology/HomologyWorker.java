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
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import topcat.matrix.BMatrix;
import topcat.matrix.BVector;
import topcat.matrix.Column;
import topcat.persistence.simplex.Simplex;
import topcat.persistence.simplex.SimplexStorageStructure;
import topcat.util.BinomialCoeffTable;
import topcat.util.IntTuple;
import topcat.util.Pair;
import topcat.persistence.simplex.SimplexCoboundaryEnumerator;
import java.util.*;

/**
 * Worker node in the parallelized computation of the homology.
 */
public class HomologyWorker implements Runnable{
        public static Logger log = LoggerFactory.getLogger(topcat.persistence.homology.HomologyWorker.class);
        int maxDimension;
        IntTuple v;
        List<List<Simplex>> chain = new ArrayList<>();
        int[] homologyDimension;
        BMatrix[] naturalTransformation;
        BMatrix[] naturalTransformation_inverse;
        SimplexStorageStructure simplexStorageStructure;
        BinomialCoeffTable binomialCoeffTable;

    HomologyWorker(SimplexStorageStructure simplexStorageStructure, IntTuple v, int maxDimension) {
            this.maxDimension = maxDimension;
            this.v = v;
            this.homologyDimension = new int[maxDimension];
            this.naturalTransformation = new BMatrix[maxDimension];
            this.naturalTransformation_inverse = new BMatrix[maxDimension];
            this.simplexStorageStructure = simplexStorageStructure;
            this.binomialCoeffTable = new BinomialCoeffTable(simplexStorageStructure.getNumberOfVertices(), maxDimension+2);

            //For each index in the grid we compute the homology up to dimension 'maxdimension'
            for (int k = 0; k <= maxDimension; k++) {
                chain.add(simplexStorageStructure.getSimplicesLEQThan(k, v));
            }
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
        public List<Pair<Long, Integer>> computePivots(List<List<Long>> columns_to_reduce, Long2IntOpenHashMap pivot_column_index){
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
    public BMatrix inverse(List<List<Long>> columns_to_reduce){
        Long2IntOpenHashMap pivot_column_index = new Long2IntOpenHashMap();
        computePivots(columns_to_reduce, pivot_column_index);
        assert pivot_column_index.size() == columns_to_reduce.size();
        BMatrix reductionMatrix = new BMatrix(columns_to_reduce.size(), columns_to_reduce.size());//BMatrix.identity(columns_to_reduce.size());
        for(int i = 0; i<columns_to_reduce.size();i++) {
            Integer index_column_to_reduce = pivot_column_index.get((long)i);
            Column<Long> working_column= new Column<>();
            working_column.addAll(columns_to_reduce.get(index_column_to_reduce));
            reductionMatrix.set(index_column_to_reduce, i, true);
            Long pivot;
            while((pivot=working_column.get_pivot()) != null) {
                if(pivot_column_index.get(pivot) == index_column_to_reduce){
                    working_column.pop_pivot();
                }
                else{
                    working_column.addAll(columns_to_reduce.get(pivot_column_index.get((long)pivot)));
                    reductionMatrix.set(pivot_column_index.get((long)pivot), i, !reductionMatrix.get(pivot_column_index.get((long)pivot), i));
                }
            }
        }
        return reductionMatrix;
    }

    /**
     * Computes a basis for the homology in each dimension of the chain complex given by 'chain'.
     * @param chain
     */
    protected void computeHomologyBasis(final List<List<Simplex>> chain) {
        List<LongOpenHashSet> index_lookups = new ArrayList<>();
        for (int i = 0; i < chain.size(); i++) {
            List<Simplex> simplices = chain.get(i);
            log.debug("Number of simplices of dim " + i + ": " + chain.get(i).size());
            LongOpenHashSet index_lookup = new LongOpenHashSet();
            for (Simplex s : simplices) index_lookup.add(s.getIndex());
            index_lookups.add(index_lookup);
        }

        List<List<Long>> image_basis = new ArrayList<>();

        for (int dim = 0; dim < maxDimension; dim++) {
            Long2IntOpenHashMap pivot_column_index = new Long2IntOpenHashMap();
            BMatrix reduction_matrix = reduce_coboundary_matrix(chain.get(dim), pivot_column_index, index_lookups.get(dim + 1));

            List<List<Long>> kernel_basis = new ArrayList<>();
            for (int i = 0; i < chain.get(dim).size(); i++) {
                if (!pivot_column_index.containsValue(i)) {
                    BVector row = reduction_matrix.getRow(i);
                    List<Long> column = new ArrayList<>();
                    IntIterator iterator = row.getIndexSetIterator();
                    while (iterator.hasNext()) {
                        column.add(((Integer) iterator.nextInt()).longValue());
                    }
                    kernel_basis.add(column);
                }
            }

            int image_basis_size = image_basis.size();
            image_basis.addAll(kernel_basis);
            List<Pair<Long, Integer>> pivot_column = computePivots(image_basis, new Long2IntOpenHashMap());

            List<List<Long>> homology_basis = new ArrayList<>();
            for (int i = image_basis_size; i < pivot_column.size(); i++) {
                homology_basis.add(image_basis.get(pivot_column.get(i)._2()));
            }

            assert homology_basis.size() == kernel_basis.size() - image_basis_size;

            int homology_basis_size = homology_basis.size();
            for (int i = 0; i < chain.get(dim).size(); i++) {
                List<Long> column = new ArrayList<>();
                column.add((long) i);
                homology_basis.add(column);
            }
            pivot_column = computePivots(homology_basis, new Long2IntOpenHashMap());
            List<BVector> extended_basis = new ArrayList<>();
            List<List<Long>> e_basis = new ArrayList<>();
            for (Pair<Long, Integer> pivot : pivot_column) {
                BVector row = new BVector(chain.get(dim).size());
                for (long pos : homology_basis.get(pivot._2()))
                    row.set((int) pos, true);
                extended_basis.add(row);
                e_basis.add(homology_basis.get(pivot._2()));
            }

            BMatrix extendedBasis_inv = inverse(e_basis);
            BMatrix extendedBasis = (new BMatrix(extended_basis)).transpose();

            homologyDimension[dim] = homology_basis_size;
            naturalTransformation[dim] = extendedBasis_inv;
            naturalTransformation_inverse[dim] = extendedBasis;

            Long2IntOpenHashMap index_column_lookup = new Long2IntOpenHashMap();
            for (int i = 0; i < chain.get(dim + 1).size(); i++)
                index_column_lookup.addTo(chain.get(dim + 1).get(i).getIndex(), i);

            image_basis = new ArrayList<>();
            int[] pivot_columns = pivot_column_index.values().toIntArray();
            for (int i = 0; i < pivot_columns.length; i++) {
                List<Long> column = new ArrayList<>();
                Simplex s = chain.get(dim).get(pivot_columns[i]);
                SimplexCoboundaryEnumerator enumerator = new SimplexCoboundaryEnumerator(s, simplexStorageStructure.getNumberOfVertices(), binomialCoeffTable);
                while (enumerator.hasNext()) {
                    long simplex = enumerator.next();
                    if (index_column_lookup.containsKey(simplex)) {
                        column.add(((long) index_column_lookup.get(simplex)));
                    }
                }
                image_basis.add(column);
            }
        }
    }

    @Override
    public void run() {
        log.debug("Starting basis change computation of position: "+v);
        computeHomologyBasis(chain);
        log.debug("Finished basis change computation of position: "+v);
    }
}


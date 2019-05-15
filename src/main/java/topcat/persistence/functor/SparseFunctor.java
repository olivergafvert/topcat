package topcat.persistence.functor;

import it.unimi.dsi.fastutil.objects.Object2IntOpenHashMap;
import topcat.matrix.BMatrix;
import topcat.matrix.GradedColumn;
import topcat.persistence.simplex.Simplex;
import topcat.util.Grid;
import topcat.util.IntTuple;
import topcat.util.Pair;
import topcat.util.SparseGrid;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

public class SparseFunctor extends Functor {

    Grid<List<GradedColumn<Simplex>>> kernel_grid = new SparseGrid<>(size);
    Grid<List<GradedColumn<Simplex>>> homology_grid = new SparseGrid<>(size);
    Grid<Integer> homology_dim_grid = new SparseGrid<>(size);

    public SparseFunctor(IntTuple size, Grid<List<GradedColumn<Simplex>>> kernel_grid, Grid<List<GradedColumn<Simplex>>> homology_grid, Grid<Integer> homology_dim_grid){
        super(size);
    }

//    /**
//     * Returns the map at position 'v' in direction 'dim'.
//     * @param v
//     * @param dim
//     * @return
//     */
//    @Override
//    public BMatrix getMap(IntTuple v, int dim){
//
//        int homologyDimension_v = homology_dim_grid.get(i).get(index);
//        List<GradedColumn<Simplex>> homologyBasis_v = homology_grid.get(i).get(index);
//        List<GradedColumn<Simplex>> kernelBasis_v = kernel_grid.get(i).get(index);
//        Object2IntOpenHashMap<Simplex> index_map_v = new Object2IntOpenHashMap<>();
//        for(int k=0;k<kernelBasis_v.size();k++)
//            index_map_v.put(kernelBasis_v.get(k).getGrade(), k);
//        List<GradedColumn<Simplex>> working_columns_v = new ArrayList<>();
//        for(int k=0;k<homologyBasis_v.size();k++){
//            GradedColumn<Simplex> column = homologyBasis_v.get(k);
//            GradedColumn<Simplex> working_column = new GradedColumn<>(new Simplex(k, i));
//            while(!column.isEmpty()){
//                working_column.addAll(kernelBasis_v.get(index_map_v.getInt(column.pop_pivot())));
//            }
//            working_columns_v.add(working_column);
//        }
//        HashMap<Simplex, GradedColumn<Simplex>> pivot_to_column_index = new HashMap<>();
//        pivots(working_columns_v, pivot_to_column_index);
//
//        for(int j=0;j<adjacent.size();j++) {
//            List<GradedColumn<Simplex>> homologyBasis = homology_grid.get(i).get(adjacent.get(j)).subList(0, homology_dim_grid.get(i).get(adjacent.get(j)));
//            List<GradedColumn<Simplex>> kernelBasis = kernel_grid.get(i).get(adjacent.get(j));
//            Object2IntOpenHashMap<Simplex> index_map = new Object2IntOpenHashMap<>();
//            for(int k=0;k<kernelBasis.size();k++)
//                index_map.put(kernelBasis.get(k).getGrade(), k);
//            List<GradedColumn<Simplex>> working_columns = new ArrayList<>();
//            for(int k=0;k<homologyBasis.size();k++){
//                GradedColumn<Simplex> column = homologyBasis.get(k);
//                GradedColumn<Simplex> working_column = new GradedColumn<>(new Simplex(k, i));
//                while(!column.isEmpty()){
//                    working_column.addAll(kernelBasis.get(index_map.getInt(column.pop_pivot())));
//                }
//                working_columns.add(working_column);
//            }
//            List<GradedColumn<Simplex>> hom_trans = basisChange(working_columns, pivot_to_column_index);
//            BMatrix M = new BMatrix(homology_dim_grid.get(i).get(index), homologyBasis.size());
//            for(int k=0;k<hom_trans.size();k++){
//                GradedColumn<Simplex> column = hom_trans.get(k);
//                while(!column.isEmpty()){
//                    Simplex grade = column.pop_pivot();
//                    if(index_map_v.getInt(grade)<homologyDimension_v){
//                        M.set(index_map_v.getInt(grade), k, true);
//                    }
//                }
//            }
//            maps.get(j).add(new Pair<>(adjacent.get(j), M));
//        }
//    }
//
//
//
//        v = parseTuple(v);
//        if(v.hasNegativeElements(v)){
//            return new BMatrix(getDimension(v.plus(IntTuple.getStandardBasisElement(v.length(), dim))), 0);
//        }
//        if(!maps.get(v).containsKey(v.plus(IntTuple.getStandardBasisElement(v.length(), dim)))) return BMatrix.identity(getDimension(v));
//        return maps.get(dim).get(v);
//    }
//
//
//    /**
//     * Sets the map 'A' at position 'v' in direction 'dim'.
//     * @param v
//     * @param A
//     * @param dim
//     */
//    @Override
//    public void setMap(IntTuple v, BMatrix A, int dim){
//        if(isOutOfBounds(v)){
//            return;
//        }
//        maps.get(v).put(v.plus(IntTuple.getStandardBasisElement(v.length(), dim)), A);
//    }
//
//    @Override
//    public void setMap(IntTuple v, IntTuple w, BMatrix A){
//        if(isOutOfBounds(v)){
//            return;
//        }
//        if(!maps.containsKey(v)) maps.put(v, new HashMap<>());
//        maps.get(v).put(w, A);
//    }
}

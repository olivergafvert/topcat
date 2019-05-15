package topcat.persistence.functor;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import topcat.matrix.BMatrix;
import topcat.matrix.GradedColumn;
import topcat.matrix.exception.NoSolutionException;
import topcat.persistence.homology.HomologyUtil;
import topcat.persistence.simplex.Simplex;
import topcat.util.Grid;
import topcat.util.GridIterator;
import topcat.util.IntTuple;
import topcat.util.Pair;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

public class SparseFunctor extends Functor {
    private static Logger log = LoggerFactory.getLogger(SparseFunctor.class);

    Grid<HashMap<Simplex, GradedColumn<Simplex>>> basis_pivot_grid;
    Grid<List<GradedColumn<Simplex>>> homology_grid;
    Grid<Integer> homology_dim_grid;
    HashSet<IntTuple> support;
    HashMap<IntTuple, List<IntTuple>> adjacent_map;

    public SparseFunctor(IntTuple size, Grid<HashMap<Simplex, GradedColumn<Simplex>>> basis_pivot_grid, Grid<List<GradedColumn<Simplex>>> homology_grid, Grid<Integer> homology_dim_grid, HashSet<IntTuple> support, HashMap<IntTuple, List<IntTuple>> adjacent_map){
        super(size);
        this.basis_pivot_grid = basis_pivot_grid;
        this.homology_grid = homology_grid;
        this.homology_dim_grid = homology_dim_grid;
        this.support = support;
        this.adjacent_map = adjacent_map;
    }


    @Override
    public BMatrix getMap(IntTuple v, IntTuple w){
        int homologyDimension_v = homology_dim_grid.get(v);
        List<GradedColumn<Simplex>> homologyBasis = homology_grid.get(w).subList(0, homology_dim_grid.get(w));
        List<GradedColumn<Simplex>> hom_trans = HomologyUtil.basisChange(homologyBasis, basis_pivot_grid.get(v));
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
        return M;
    }


    @Override
    public BMatrix getMap(IntTuple v, int dim){
        throw new UnsupportedOperationException("Not supported for sparse functors.");
    }

    @Override
    public void setMap(IntTuple v, BMatrix A, int dim){
        throw new UnsupportedOperationException("Not supported for sparse functors.");
    }

    @Override
    public void setMap(IntTuple v, IntTuple w, BMatrix A){
        throw new UnsupportedOperationException("Not supported for sparse functors.");
    }

    @Override
    public List<Generator> getGenerators(){
        List<Generator> generators = new ArrayList<>();
        for(IntTuple v : support){
            try {
                List<IntTuple> basis = adjacent_map.get(v);
                if(basis.size()>0) {
                    BMatrix A = getMap(basis.get(0), v);
                    for (int i = 1; i < basis.size(); i++) {
                        A = BMatrix.concat(A, getMap(v.minus(basis.get(i)), v));
                    }
                    Pair<BMatrix, BMatrix> kerim = BMatrix.reduction(A);
                    if (kerim._2().rows != kerim._2().cols) {
                        BMatrix generatorBasis = BMatrix.extendBasis(kerim._2());
                        generatorBasis = generatorBasis.subMatrix(kerim._2().rows, -1, 0, -1);
                        for (int k = 0; k < generatorBasis.rows; k++) {
                            generators.add(new Generator(new IntTuple(v), generatorBasis.getRow(k)));
                        }
                    }
                }else{
                    BMatrix A = BMatrix.identity(homology_dim_grid.get(v));
                    for (int k = 0; k < A.rows; k++) {
                        generators.add(new Generator(new IntTuple(v), A.getRow(k)));
                    }
                }
            }catch (NoSolutionException nse){
                log.info("Failed to compute minimal set of generators.", nse);
            }
        }
        return generators;
    }
}

package topcat.persistence.functor;

import topcat.matrix.BMatrix;
import topcat.util.Grid;
import topcat.util.IntTuple;
import topcat.util.SparseGrid;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

public class SparseFunctor extends Functor {

    HashMap<IntTuple, HashMap<IntTuple, BMatrix>> maps;

    public SparseFunctor(IntTuple size){
        super(size);
        this.maps = new HashMap<>();
    }

    /**
     * Returns the map at position 'v' in direction 'dim'.
     * @param v
     * @param dim
     * @return
     */
    @Override
    public BMatrix getMap(IntTuple v, int dim){
        v = parseTuple(v);
        if(v.hasNegativeElements(v)){
            return new BMatrix(getDimension(v.plus(IntTuple.getStandardBasisElement(v.length(), dim))), 0);
        }
        if(!maps.get(v).containsKey(v.plus(IntTuple.getStandardBasisElement(v.length(), dim)))) return BMatrix.identity(getDimension(v));
        return maps.get(dim).get(v);
    }


    /**
     * Sets the map 'A' at position 'v' in direction 'dim'.
     * @param v
     * @param A
     * @param dim
     */
    @Override
    public void setMap(IntTuple v, BMatrix A, int dim){
        if(isOutOfBounds(v)){
            return;
        }
        maps.get(v).put(v.plus(IntTuple.getStandardBasisElement(v.length(), dim)), A);
    }

    @Override
    public void setMap(IntTuple v, IntTuple w, BMatrix A){
        if(isOutOfBounds(v)){
            return;
        }
        if(!maps.containsKey(v)) maps.put(v, new HashMap<>());
        maps.get(v).put(w, A);
    }
}

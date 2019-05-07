package topcat.util;

import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.ints.IntOpenHashSet;

import java.util.*;

public class SparseGrid<T> extends Grid<T> {

    HashMap<IntTuple, T> grid = new HashMap<>();

    public SparseGrid(IntTuple size) {
        super(size);
    }

    @Override
    public T get(IntTuple v){
        if(grid.containsKey(v))
            return grid.get(v);
        return null;
    }

    @Override
    public void set(IntTuple v, T t){
        grid.put(v, t);
    }

    @Override
    public List<IntTuple> indexSet(){
        return new ArrayList<>(grid.keySet());
    }

    @Override
    public List<IntTuple> indexSet(IntTuple v){
        List<IntTuple> indices = new ArrayList<>();
        for(IntTuple w : grid.keySet()){
            if(w.leq(v))
                indices.add(w);
        }
        return indices;
    }

    public boolean containsIndex(IntTuple v){
        return grid.containsKey(v);
    }




    public static SparseGrid create(List<IntTuple> indices){
        Int2ObjectOpenHashMap<IntOpenHashSet> mins = new Int2ObjectOpenHashMap<>();
        Collections.sort(indices, new Comparator<IntTuple>() {
            @Override
            public int compare(IntTuple o1, IntTuple o2) {
                if(o1.lexLt(o2)) return 1;
                if(o2.lexLt(o1)) return -1;
                return 0;
            }
        });
        for(int i=0;i<indices.size();i++){
            for(int j=i+1;j<indices.size();j++){
                if(indices.get(i).leq(indices.get(j))){
                    //Build incidence sets recursively. If v<w then neighbour set is contained
                }
            }
        }
        return null;
    }


}

package topcat.matrix;

import it.unimi.dsi.fastutil.longs.Long2IntOpenHashMap;

import java.util.ArrayList;
import java.util.List;

public class PMatrix {
    int rows, cols;
    List<Column<Integer>> row_vectors;

    public PMatrix(BMatrix M){
        this.rows = M.rows;
        this.cols = M.cols;
        this.row_vectors = new ArrayList<>();
        for(int i=0;i<M.rows;i++){
            if(M.hasCreatedRow(i)){
                if(M.getRow(i).getNumberOfNonZeroElements()>0){
                    Column<Integer> indices = new Column<>();
                    indices.addAll(M.getRow(i).getIndexSet());
                    this.row_vectors.add(indices);
                }
            }
        }
    }


    public static int rank(BMatrix M){
        if(M==null) return 0;
        List<List<Integer>> columns_to_reduce = new ArrayList<>();
        for(int i=0;i<M.rows;i++)
            if(M.hasCreatedRow(i) && M.getRow(i).getNumberOfNonZeroElements()>0)
                columns_to_reduce.add(M.getRow(i).getIndexSet());
        Long2IntOpenHashMap pivot_column_index = new Long2IntOpenHashMap();
        for(Integer index_column_to_reduce = 0; index_column_to_reduce<columns_to_reduce.size();index_column_to_reduce++) {
            Column<Integer> working_column= new Column<>();
            working_column.addAll(columns_to_reduce.get(index_column_to_reduce));
            long pivot;
            while((pivot=working_column.get_pivot()) != -1) {
                if(pivot_column_index.containsKey(pivot)){
                    working_column.addAll(columns_to_reduce.get(pivot_column_index.get(pivot)));
                }else{
                    pivot_column_index.addTo(pivot, index_column_to_reduce);
                    break;
                }
            }
        }
        return pivot_column_index.size();
    }


}

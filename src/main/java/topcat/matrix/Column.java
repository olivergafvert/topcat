package topcat.matrix;

import java.util.Collections;
import java.util.PriorityQueue;

public class Column<T> extends PriorityQueue<T> {

    public Column(){
        super(10, Collections.reverseOrder());
    }


    /**
     * Helper function to @get_pivot.
     * @return
     */
    public T pop_pivot(){
        if(this.isEmpty()) return null;
        else{
            T pivot = this.poll();
            while(!this.isEmpty() && this.peek() == pivot){
                this.poll();
                if(this.isEmpty()) return null;
                else{
                    pivot = this.poll();
                }
            }
            return pivot;
        }
    }

    /**
     * Returns the largest index of the column which is non-zero.
     * @return
     */
    public T get_pivot(){
        T result = pop_pivot();
        if(result != null) this.add(result);
        return result;
    }

    public Column<T> plus(Column<T> v){
        Column<T> w = new Column<>();
        w.addAll(v);
        w.addAll(this);
        return w;
    }

    public static Column<Integer> create(BVector v){
        Column<Integer> c = new Column<>();
        c.addAll(v.getIndexSet());
        return c;
    }
}

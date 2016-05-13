package topcat.util.paralelliterator;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.ArrayList;

/**
 * Created by oliver on 2016-05-13.
 */
public abstract class ParalellIntIterator<T> extends ParalellIterator<Integer, T> {
    private Logger log = LoggerFactory.getLogger(ParalellIterator.class);
    private int from, to;

    public ParalellIntIterator(int from, int to){
        super(null);
        this.indices = new ArrayList<>();
        for(int i=from;i<to;i++){
            this.indices.add(i);
        }
    }

}

package topcat.util.paralelliterator;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import topcat.util.Pair;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.*;

/**
 * Created by oliver on 2016-05-13.
 */
public abstract class ParalellIterator<S, T> {
    private Logger log = LoggerFactory.getLogger(ParalellIterator.class);
    protected List<S> indices;

    public ParalellIterator(List<S> indices){
        this.indices = indices;
    }

    public abstract T method(S index);

    public List<Pair<S, T>> run(){
        ExecutorService exec = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());
        List<Future> futures = new ArrayList<>();
        List<Worker> workers = new ArrayList<>();
        for(S i : indices){
            Worker worker = new Worker(i);
            workers.add(worker);
            futures.add(exec.submit(worker));
        }
        exec.shutdown();

        while(!exec.isTerminated()){
            try{
                exec.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
            }catch (InterruptedException ire){
                log.error("Failed to terminate executor service.", ire);
            }
        }
        try {
            for (Future future : futures) {
                future.get();
            }
        }catch (InterruptedException ire){
            log.error("Failed to terminate executor service.", ire);
        }catch (ExecutionException exe){
            log.error("Failed to terminate executor service.", exe);
        }

        List<Pair<S, T>> results = new ArrayList<>();
        for(Worker worker : workers){
            results.add(new Pair<S, T>(worker.i, worker.result));
        }
        return results;
    }


    private class Worker implements Runnable{
        S i;
        T result;
        private Worker(S i){
            this.i = i;
        }

        @Override
        public void run() {
            result = method(i);
        }
    }

}

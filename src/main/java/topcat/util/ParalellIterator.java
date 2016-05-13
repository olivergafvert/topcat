package topcat.util;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.*;

/**
 * Created by oliver on 2016-05-13.
 */
public abstract class ParalellIterator<T> {
    private Logger log = LoggerFactory.getLogger(ParalellIterator.class);
    private int from, to;

    public ParalellIterator(int from, int to){
        this.from = from;
        this.to = to;
    }

    public abstract T method(int index);

    public List<Pair<Integer, T>> run(){
        ExecutorService exec = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());
        List<Future> futures = new ArrayList<>();
        List<Worker> workers = new ArrayList<>();
        for(int i=from;i<to;i++){
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

        List<Pair<Integer, T>> results = new ArrayList<>();
        for(Worker worker : workers){
            results.add(new Pair<Integer, T>(worker.i, worker.result));
        }
        return results;
    }


    private class Worker implements Runnable{
        int i;
        T result;
        private Worker(int i){
            this.i = i;
        }

        @Override
        public void run() {
            result = method(i);
        }
    }

}

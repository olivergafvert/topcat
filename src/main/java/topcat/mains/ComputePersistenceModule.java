package topcat.mains;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import topcat.persistence.PersistenceModuleCollection;
import topcat.persistence.simplex.SimplexStorageStructure;
import topcat.persistence.simplex.SparseSimplexStorageStructure;

import java.util.List;

public class ComputePersistenceModule {
    private static Logger log = LoggerFactory.getLogger(SimplexStorageStructure.class);
    /**
     * Computes the stable rank function of a given persistence contour for a given simplicial complex.
     *
     * Usage: simplex_file maxDimension contour [domain or standard]
     *
     * The simplex file should be in the following format:
     * number of filtrationdirections
     * filtration values in direction 1
     * ...
     * vertices in simplex : filtration indices
     * ...
     *
     *
     * @param args
     */
    public static void main(String[] args){
        log.debug("Reading input...");
        SimplexStorageStructure simplexStorageStructure;
        try {
            simplexStorageStructure = SparseSimplexStorageStructure.readOFFFile(args[0]);
        }catch(Exception ioe){
            ioe.printStackTrace();
            return;
        }
        log.debug("Finished reading input.");
        log.debug("Computing persistence modules...");
        int maxDimension = 2;
        List<List<Double>> filtrationValues = simplexStorageStructure.getFiltrationValues();
        PersistenceModuleCollection persistenceModules = PersistenceModuleCollection.create(simplexStorageStructure, filtrationValues, maxDimension);
        log.debug("Finished computing persistence modules.");

    }
}

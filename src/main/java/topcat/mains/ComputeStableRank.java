/*
Topcat - a multidimensional persistent homology library.
Copyright (C) 2016 Oliver Gäfvert

This file is part of Topcat.

Topcat is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Topcat is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

package topcat.mains;

import topcat.persistence.PersistenceModule;
import topcat.persistence.PersistenceModuleCollection;
import topcat.persistence.contours.ExponentialContour;
import topcat.persistence.contours.PersistenceContour;
import topcat.persistence.contours.StandardContour;
import topcat.persistence.stablerank.StableRankFunction;
import topcat.persistence.simplex.SimplexStorageStructure;

import java.io.File;
import java.io.IOException;
import java.util.List;

public class ComputeStableRank {

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
        SimplexStorageStructure simplexStorageStructure;
        try {
            simplexStorageStructure = SimplexStorageStructure.readFromFile(new File(args[0]));
        }catch(IOException ioe){
            ioe.printStackTrace();
            return;
        }
        int maxDimension = Integer.parseInt(args[1]);
        List<List<Double>> filtrationValues = simplexStorageStructure.getFiltrationValues();
        PersistenceModuleCollection persistenceModules = PersistenceModuleCollection.create(simplexStorageStructure, filtrationValues, maxDimension);

        PersistenceContour persistenceContour;
        if(args[2].equals("exponentaial")){
            persistenceContour = new ExponentialContour(filtrationValues);
        }else{
            persistenceContour = new StandardContour(filtrationValues);
        }

        for(int i=0;i<persistenceModules.getMaxDimension();i++){
            PersistenceModule persistenceModule = persistenceModules.get(i);
            StableRankFunction stableRankFunction = persistenceModule.computeStableRank(filtrationValues.get(0), persistenceContour);
            System.out.println(i+": "+ stableRankFunction);
        }
    }

}

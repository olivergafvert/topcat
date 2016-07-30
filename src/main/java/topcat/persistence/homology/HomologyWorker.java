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

package topcat.persistence.homology;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import topcat.matrix.BMatrix;
import topcat.matrix.BVector;
import topcat.matrix.exception.NoSolutionException;
import topcat.persistence.simplex.Simplex;
import topcat.persistence.simplex.SimplexStorageStructure;
import topcat.util.IntTuple;
import topcat.util.Pair;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

/**
 * Created by oliver on 2016-07-28.
 */
class HomologyWorker implements Runnable {
    public static Logger log = LoggerFactory.getLogger(HomologyWorker.class);
    int maxDimension;
    IntTuple v;
    List<List<Simplex>> chain = new ArrayList<>();
    int[] homologyDimension;
    BMatrix[] naturalTransformation;
    BMatrix[] naturalTransformation_inverse;


    HomologyWorker(SimplexStorageStructure simplexStorageStructure, IntTuple v, int maxDimension){
        this.maxDimension = maxDimension;
        this.v = v;
        this.homologyDimension = new int[maxDimension];
        this.naturalTransformation = new BMatrix[maxDimension];
        this.naturalTransformation_inverse = new BMatrix[maxDimension];

        //For each index in the grid we compute the homology up to dimension 'maxdimension'
        for (int k = 0; k <= maxDimension; k++) {
            chain.add(simplexStorageStructure.getSimplicesLEQThan(k, v));
        }
    }

    /**
     * Computes the boundary matrix \partial_n: C_n \to C_n-1.
     * @param lower the simplices in C_n-1
     * @param current the simplices in C_n
     * @return
     */
    private BMatrix computeBoundaryMatrix(List<Simplex> lower, List<Simplex> current){
        HashMap<Simplex, BVector> t_matrix = new HashMap<>();
        for(int i=0;i<lower.size();i++){
            t_matrix.put(lower.get(i), new BVector(current.size()));
        }
        for(int i=0;i<current.size();i++){
            Simplex s = current.get(i);
            for(Simplex ss : s.getBoundaryList()){
                t_matrix.get(ss).set(i, true);
            }
        }
        BMatrix boundaryMatrix = new BMatrix(lower.size(), current.size());
        for(int i=0;i<lower.size();i++){
            boundaryMatrix.setRow(i, t_matrix.get(lower.get(i)));
        }
        return boundaryMatrix;
    }

    /**
     * Computes a basis for Z that contains B as a sub-basis.
     * @param Z
     * @param B
     * @return
     * @throws NoSolutionException
     */
    private static Pair<BMatrix, BMatrix> computeNormalizedBasis(BMatrix Z, BMatrix B) throws NoSolutionException {
        BMatrix Hbasis = BMatrix.getBasis(BMatrix.concat(B.transpose(), Z.transpose()).transpose()).subMatrix(B.rows, -1, 0, -1);
        return new Pair<>(Hbasis, B);
    }

    /**
     * Computes a basis for the homology in each position of the chain complex
     *   C_n -> C_n-1 -> ... -> C_1 -> C_0
     * where C_i is given by the simplices at position i in the list 'chain'.
     * @param chain - a list of lists of simplices describing a chain complex.
     * @return - a basis for the homology in each position in the chain complex.
     */
    private List<BMatrix> computeHomologyBasis(List<List<Simplex>> chain){
        BMatrix[] Z = new BMatrix[maxDimension]; //The cycle subspaces
        BMatrix[] B = new BMatrix[maxDimension]; //The boundary subspaces
        int[] dims = new int[maxDimension];

        for (int k = 1; k <= maxDimension; k++) {
            //The chain of the lower dimension
            List<Simplex> lower = chain.get(k - 1);
            dims[k - 1] = lower.size();

            //The chain of the current dimension
            List<Simplex> current = chain.get(k);

            //Compute the boundary matrix
            BMatrix boundaryMatrix = computeBoundaryMatrix(lower, current);
            if (boundaryMatrix != null) {
                Pair<BMatrix, BMatrix> basis = BMatrix.reduction(boundaryMatrix);
                Z[k - 1] = basis._1();
                B[k - 1] = basis._2();
            } else {
                Z[k - 1] = null;
                B[k - 1] = null;
            }
        }
        //Compute the quotients Z/B = H;
        List<BMatrix> homologyBasis = new ArrayList<>(); //The homology subspaces
        if (B[0] == null || B[0].rows == 0) {
            homologyBasis.add(BMatrix.identity(dims[0]));
        } else {
            BMatrix H = BMatrix.extendBasis(B[0]).subMatrix(B[0].rows, -1, 0, -1);
            homologyBasis.add(H);
        }
        for (int k = 1; k < maxDimension; k++) {
            if (B[k] == null || B[k].rows == 0) {
                if (Z[k - 1] == null || Z[k - 1].rows == 0) {
                    homologyBasis.add(new BMatrix(0, dims[k]));
                } else {
                    BMatrix H = Z[k - 1];
                    homologyBasis.add(H);
                }
            } else {
                if(Z[k-1].rows == B[k].rows){
                    homologyBasis.add(new BMatrix(0, dims[k]));
                }else {
                    Pair<BMatrix, BMatrix> basis = computeNormalizedBasis(Z[k - 1], B[k]);
                    BMatrix H = basis._1();
                    homologyBasis.add(H);
                }
            }
        }
        return homologyBasis;
    }

    @Override
    public void run() {
        log.debug("Starting basis change computation of position: "+v);
        List<BMatrix> homologyBasis = computeHomologyBasis(chain);
        for(int k=0;k<maxDimension;k++){
            BMatrix H = homologyBasis.get(k);
            BMatrix M = BMatrix.extendBasis(H).transpose();
            BMatrix Minv = M.inverse();
            homologyDimension[k] = H.rows;
            naturalTransformation[k] = Minv;
            naturalTransformation_inverse[k] = M;
        }
        log.debug("Finished basis change computation of position: "+v);
    }
}

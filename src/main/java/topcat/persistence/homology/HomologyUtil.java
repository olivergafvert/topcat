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

import gnu.trove.map.hash.TObjectIntHashMap;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import topcat.persistence.functor.Functor;
import topcat.persistence.functor.Nat;
import topcat.persistence.functor.exception.MalformedFunctorException;
import topcat.matrix.BMatrix;
import topcat.matrix.BVector;
import topcat.matrix.exception.NoSolutionException;
import topcat.matrix.exception.WrongDimensionException;
import topcat.persistence.simplex.Simplex;
import topcat.persistence.simplex.SimplexStorageStructure;
import topcat.util.Grid;
import topcat.util.GridIterator;
import topcat.util.IntTuple;
import topcat.util.Pair;
import topcat.util.paralelliterator.ParalellIterator;

import java.util.*;

/**
 * Tools for computing the homology of a multifiltration.
 */
public class HomologyUtil {
    private static Logger log = LoggerFactory.getLogger(HomologyUtil.class);

    /**
     * Computes the boundary matrix \partial_n: C_n \to C_n-1.
     * @param lower the simplices in C_n-1
     * @param current the simplices in C_n
     * @return
     */
    private static BMatrix computeBoundaryMatrix(List<Simplex> lower, List<Simplex> current){
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
     * Computes the inclusion map C^i_n \hookrightarrow C^{i+j}_n.
     * @param lower the simplices in C^i_n
     * @param current the simplices in C^{i+j}_n
     * @return
     */
    private static BMatrix computeInclusionMap(List<Simplex> lower, List<Simplex> current){
        if(lower.size() == 0){
            return new BMatrix(current.size(), lower.size());
        }
        BMatrix A = new BMatrix(current.size(), lower.size());
        TObjectIntHashMap<Simplex> lowerMap = new TObjectIntHashMap<>();
        for(int i=0;i<lower.size();i++){
            lowerMap.put(lower.get(i), i);
        }
        for(int i=0;i<current.size();i++){
            if(lowerMap.containsKey(current.get(i))){
                A.set(i, lowerMap.get(current.get(i)), true);
            }
        }
        return A;
    }

    /**
     * Computes the chain functors C_n: Q^2 \to Vect_K for each n <= 'maxdimension'
     * @param simplexStorageStructure
     * @param size
     * @param maxDimension
     * @return
     * @throws MalformedFunctorException
     * @throws WrongDimensionException
     */
    private static List<Functor> computeChainFunctors(final SimplexStorageStructure simplexStorageStructure, IntTuple size, int maxDimension) throws MalformedFunctorException {
        log.debug("Starting to compute chain functors...");
        List<Functor> chainFunctors = new ArrayList<>();

        for(int k=0;k<maxDimension+1;k++){
            final Functor C = new Functor(size);
            //Compute the maps
            for(IntTuple v : GridIterator.getSequence(size)) {
                List<Simplex> currentSimplices = simplexStorageStructure.getSimplicesLEQThan(k, v);
                for(int i=0;i<v.length();i++){
                    C.setMap(v, computeInclusionMap(currentSimplices,
                            simplexStorageStructure.getSimplicesLEQThan(k, v.plus(IntTuple.getStandardBasisElement(v.length(), i)))), i);
                }
            }
            chainFunctors.add(C);
        }

        for(Functor C : chainFunctors){
            Functor.verify(C);
        }
        log.debug("Finished to computing chain functors.");

        return chainFunctors;
    }

    /**
     * Computes a basis for the homology in each position of the chain complex
     *   C_n -> C_n-1 -> ... -> C_1 -> C_0
     * where C_i is given by the simplices at position i in the list 'chain'.
     * @param chain - a list of lists of simplices describing a chain complex.
     * @return - a basis for the homology in each position in the chain complex.
     */
    private static List<BMatrix> computeHomologyBasis(List<List<Simplex>> chain){
        int maxDimension = chain.size()-1;
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

    /**
     * Computes the homology functors for each dimension less than 'maxdimension'.
     * @param simplexStorageStructure
     * @param size
     * @param maxDimension
     * @throws WrongDimensionException
     * @throws NoSolutionException
     */
    public static List<Functor> computeHomologyFunctors(final SimplexStorageStructure simplexStorageStructure, final IntTuple size, final int maxDimension) throws MalformedFunctorException, NoSolutionException{
        log.debug("Starting to compute homology functors...");
        List<Functor> chainFunctors = computeChainFunctors(simplexStorageStructure, size, maxDimension);

        //The natural transformations from the chain functors to a basis change of the chain modules
        final List<Nat> naturalTransformation = new ArrayList<>();
        final List<Nat> naturalTransformation_inverse = new ArrayList<>();

        //The persistence modules, i.e homology of the multifiltration
        List<Functor> homfunctors = new ArrayList<>();

        for(int k=0;k<maxDimension;k++){
            naturalTransformation.add(new Nat(size));
            naturalTransformation_inverse.add(new Nat(size));
            homfunctors.add(new Functor(chainFunctors.get(k)));
        }

        //The dimension of the basis for the homology in C(v) for each v in N^r
        final List<Grid<Integer>> homologyDimension = new ArrayList<>();
        for(int i=0;i<maxDimension;i++) {
            homologyDimension.add(Grid.create(size));
        }

        log.debug("Starting to compute basis change in each position...");
        for(IntTuple v : GridIterator.getSequence(size)) {
            log.debug("Starting to compute position "+v+"...");
            //For each index in the grid we compute the homology up to dimension 'maxdimension'
            List<List<Simplex>> chain = new ArrayList<>();
            for (int k = 0; k <= maxDimension; k++) {
                chain.add(simplexStorageStructure.getSimplicesLEQThan(k, v));
            }
            List<BMatrix> homologyBasis = computeHomologyBasis(chain);
            for(int k=0;k<homologyBasis.size();k++){
                BMatrix H = homologyBasis.get(k);
                BMatrix M = BMatrix.extendBasis(H).transpose();
                BMatrix Minv = M.inverse();
                homologyDimension.get(k).set(v, H.rows);
                naturalTransformation.get(k).setMap(v, Minv);
                naturalTransformation_inverse.get(k).setMap(v, M);
            }
            log.debug("Finished computing position "+v+".");
        }
        log.debug("Finished computing basis change.");

        log.debug("Starting to apply basis change...");
        for(int k=0;k<homfunctors.size();k++){
            Functor H = homfunctors.get(k);
            for(IntTuple v : GridIterator.getSequence(size)){
                for(int i=0;i<v.length();i++){
                    if(v.get(i).equals(size.get(i))){
                        H.setMap(v, BMatrix.identity(naturalTransformation.get(k).getMap(v).rows), i);
                    }else{
                        H.setMap(v, H.getMap(v, i).mult(naturalTransformation_inverse.get(k).getMap(v)), i);
                    }
                    if(v.get(i) > 0){
                        IntTuple w = v.minus(IntTuple.getStandardBasisElement(v.length(), i));
                        BMatrix hmap = H.getMap(w, i);
                        BMatrix natmap = naturalTransformation.get(k).getMap(v);
                        H.setMap(w, natmap.mult(hmap), i);
                    }
                }
            }
        }
        log.debug("Finished applying basis change.");

        log.debug("Starting verification...");
        for(int k=0;k<maxDimension;k++){
            Functor.verify(homfunctors.get(k));
            log.debug("Basis change dimension "+k+" OK.");
            Nat.verify(chainFunctors.get(k), homfunctors.get(k), naturalTransformation.get(k));
            log.debug("Basis change map dimension "+k+" OK.");
            Pair<Functor, Nat> gnat = homfunctors.get(k).getSubFunctor(homologyDimension.get(k));
            Functor.verify(gnat._1());
            log.debug("Homology functor dimension "+k+" OK.");
            homfunctors.set(k, gnat._1());
        }
        log.debug("Finished verification.");

        log.debug("Finished computing homology functors.");
        return homfunctors;
    }
}

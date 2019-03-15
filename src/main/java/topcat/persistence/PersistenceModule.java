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

package topcat.persistence;

import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import topcat.matrix.BMatrix;
import topcat.matrix.BVector;
import topcat.matrix.exception.AffineVectorSpaceDimensionException;
import topcat.matrix.exception.NoSolutionException;
import topcat.matrix.exception.WrongDimensionException;
import topcat.matrix.rankminimization.AffineVectorSpaceIterator;
import topcat.matrix.rankminimization.RankTreeSearch;
import topcat.persistence.contours.PersistenceContour;
import topcat.persistence.contours.StandardContour;
import topcat.persistence.functor.Functor;
import topcat.persistence.landscape.PersistenceLandscape;
import topcat.persistence.stablerank.StableRankFunction;
import topcat.util.IntTuple;
import topcat.util.Pair;

import java.util.ArrayList;
import java.util.List;

/**
 * Represents a persistence module, i.e a functor F: Q^r -> Vect_K. It is implemented as a functor F: N^r -> Vect_K
 * where we have a map \alpha : N^r -> Q^r mapping the filtrationvalues (see [1] for details).
 *
 * [1] - Stable Invariants for Multidimensional Persistence by Gäfvert and Chachólski (arXiv:1703.03632).
 */
public class PersistenceModule {
    private static final Logger log = LoggerFactory.getLogger(PersistenceModule.class);
    protected int dimension;
    protected Functor F;

    //Tracks the indexing of the filtration values in each dimension
    protected List<List<Double>> filtrationValues;

    PersistenceModule(Functor F, int dimension, List<List<Double>> filtrationValues){
        this.F = F;
        this.dimension = dimension;
        this.filtrationValues = filtrationValues;
    }

    public Functor getFunctor(){
        return F;
    }

    public int getDimension(){
        return dimension;
    }

    public List<List<Double>> getFiltrationValues(){
        return filtrationValues;
    }

    /**
     * Computes the rank of the map from index 'from' to index 'to'.
     * @param from
     * @param to
     * @return
     */
    public Integer rank(IntTuple from, IntTuple to){
        return BMatrix.rank(F.getMap(from, to));
    }

    public Integer rank(List<Integer> from, List<Integer> to) { return rank(new IntTuple(from), new IntTuple(to)); }

    /**
     * Computes the cartesian persistence landscape, as described in [1].
     *
     * [1] - Multiparameter Persistence Landscapes, Oliver Vipond (arXiv:1812.09935).
     * @param p - using the p-norm
     * @return
     */
    public PersistenceLandscape persistenceLandscape(Integer p){
        return PersistenceLandscape.cartesian(this, p);
    }

    /**
     * Computes the Stable Rank of the persistence module at shift values 'epsilon' with
     * respect to the standard contour.
     * @param epsilons
     * @return
     */
    public StableRankFunction computeStableRank(List<Double> epsilons){
        PersistenceContour contour = new StandardContour(filtrationValues);
        return computeStableRank(epsilons, contour);
    }

    /**
     * Computes the Stable Rank of a the persistence module at shift values 'epsilon' with
     * respect to the persistence contour 'contour'.
     * @param epsilons
     * @param contour
     * @return
     */
    public StableRankFunction computeStableRank(List<Double> epsilons, PersistenceContour contour){
        StableRankFunction stableRankFunction = new StableRankFunction();

        List<Functor.Generator> f_generators = F.getGenerators();

        stableRankFunction.add(new Pair<>(0.0, f_generators.size()));

        for(int i=1;i<epsilons.size();i++){
            Double epsilon = epsilons.get(i);
            log.debug("Shift value: "+epsilon);
            try {
                Integer bar = computeStableRank(F, f_generators, epsilon, contour)._1();
                stableRankFunction.add(new Pair<>(epsilon, bar));
                log.debug("Bar: " + bar);
            }catch (WrongDimensionException wde){
                log.error("Failed to compute bar", wde);
            }
        }
        log.debug(stableRankFunction.toString());
        return stableRankFunction;
    }

    private Pair<Integer, BMatrix> computeStableRank(Functor F, List<Functor.Generator> f_generators, Double epsilon, PersistenceContour contour)  throws WrongDimensionException {
        return computeStableRank(F, f_generators, epsilon, contour, Long.MAX_VALUE, 10);
    }

    /**
     * Computes the Stable Rank invariant for of a functor F at shift epsilon.
     * @param F
     * @param f_generators
     * @param epsilon
     * @param contour
     * @param threshold
     * @param cache_size
     * @return
     * @throws WrongDimensionException
     */
    private Pair<Integer, BMatrix> computeStableRank(Functor F, List<Functor.Generator> f_generators, Double epsilon, PersistenceContour contour, long threshold, int cache_size)  throws WrongDimensionException{
        //f_generators.forEach(f -> System.out.println(f));
        List<Functor.Generator> g_generators = F.generatorShift(f_generators, epsilon, contour);

        List<AffineVectorSpaceIterator> solution_sets = new ArrayList<>();


        //Remove all f-generators that don't map to any g-generator
        List<Functor.Generator> non_zero_f_generators = new ArrayList<>();

        for(Functor.Generator f : f_generators){
            boolean is_nonzero = false;
            for(Functor.Generator g : g_generators){
                BMatrix M = F.getMap(f.position, g.position);
                if(M==null) continue;
                BVector v = M.mult(f.v);
                if(v.getNumberOfNonZeroElements() > 0) {
                    is_nonzero = true;
                    break;
                }
            }
            if(is_nonzero){
                non_zero_f_generators.add(f);
            }
        }

        f_generators = non_zero_f_generators;


        //Calcutale the solution-sets
        for(Functor.Generator generator : g_generators){
            //The position of the G-generator in F

            List<BVector> vectors = new ArrayList<>();
            IntOpenHashSet incomparables = new IntOpenHashSet();
            for(int i=0;i<f_generators.size();i++){
                BMatrix M = F.getMap(f_generators.get(i).position, generator.position);
                if(M!=null){
                    vectors.add(M.mult(f_generators.get(i).v));
                }else{
                    incomparables.add(i);
                }
            }

            BMatrix A = new BMatrix(vectors);
            BMatrix At = A.transpose();

            BVector v = null;
            BMatrix K = null;
            try {
                v = BMatrix.solve(At, generator.v);
                K = BMatrix.ker(At).transpose();
            }catch (NoSolutionException nse){
                log.info("Failed to solve linear system when computing affine solution set.", nse);
            }

            BMatrix Kn = new BMatrix(K.rows+incomparables.size(), K.cols);
            BVector vn = new BVector(v.getLength()+incomparables.size());
            int k=0;
            for(int i=0;i<Kn.rows;i++){
                if(!incomparables.contains(i)){
                    Kn.setRow(i, K.getRow(k));
                    vn.set(i, v.get(k));
                    k++;
                }
            }

            AffineVectorSpaceIterator iterator;
            try{
                iterator = AffineVectorSpaceIterator.create(Kn.transpose(), vn, cache_size);
                solution_sets.add(iterator);
            }catch (AffineVectorSpaceDimensionException afe){
                log.error("Failed to create affine vector space iterator.", afe);
                return new Pair<>(-1, null);
            }
        }
        RankTreeSearch rankTreeSearch = new RankTreeSearch(solution_sets);
        Pair<Integer, BMatrix> bar = rankTreeSearch.findMinRank(threshold);
        return bar;
    }
}

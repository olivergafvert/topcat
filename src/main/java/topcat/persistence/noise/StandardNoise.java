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

package topcat.persistence.noise;

import gnu.trove.set.hash.TIntHashSet;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import topcat.matrix.BMatrix;
import topcat.matrix.BVector;
import topcat.matrix.exception.AffineVectorSpaceDimensionException;
import topcat.matrix.exception.NoSolutionException;
import topcat.matrix.rankminimization.AffineVectorSpaceIterator;
import topcat.matrix.rankminimization.RankTreeSearch;
import topcat.persistence.barcode.BasicBarcode;
import topcat.persistence.functor.Functor;
import topcat.matrix.exception.WrongDimensionException;
import topcat.util.IntTuple;
import topcat.util.Pair;

import java.util.ArrayList;
import java.util.List;

/**
 * Represents standard noise in the direction of a ray (which is equivalent to
 * standard noise in the direction of a cone). For details see the paper
 * Multidimensional Persistence and Noise by Chachólski et al. (arXiv:1505.06929).
 */
public class StandardNoise extends Noise{
    private static final Logger log = LoggerFactory.getLogger(StandardNoise.class);
    private long threshold; //The maximal size of the search space
    private int cache_size;

    public StandardNoise(){
        this(Long.MAX_VALUE, 10);
    }

    public StandardNoise(long threshold, int cache_size){
        this.threshold = threshold;
        this.cache_size = cache_size;
    }

    /**
     * Computes the basic barcode in the diagonal direction using standard noise in the direction of a cone.
     * @param F
     * @param filtrationValues
     * @return
     */
    public BasicBarcode computeBasicBarcode(Functor F, List<List<Double>> filtrationValues){
        BasicBarcode basicBarcode = new BasicBarcode();
        List<IntTuple> indices = getIndexSequence(filtrationValues);

        List<Functor.Generator> f_generators = F.getGenerators();

        basicBarcode.add(new Pair<>(0.0, f_generators.size()));

        for(int i=1;i<indices.size();i++){
            IntTuple epsilon = indices.get(i);
            log.debug("Filtrationindices: "+epsilon);
            try {
                Integer bar = computeBasicBarcode(F, f_generators, epsilon)._1();
                double max = -Double.MAX_VALUE;
                for(int j=0;j<epsilon.length();j++){
                    if(filtrationValues.get(j).get(epsilon.get(j)) > max){
                        max = filtrationValues.get(j).get(epsilon.get(j));
                    }
                }
                basicBarcode.add(new Pair<>(max, bar));
                log.debug("Bar: " + bar);
            }catch (WrongDimensionException wde){
                log.error("Failed to compute bar", wde);
            }
        }
        log.debug(basicBarcode.toString());
        return basicBarcode;
    }

    private Pair<Integer, BMatrix> computeBasicBarcode(Functor F, List<Functor.Generator> f_generators, IntTuple epsilon)  throws WrongDimensionException{
        //f_generators.forEach(f -> System.out.println(f));
        List<Functor.Generator> g_generators = F.generatorShift(f_generators, epsilon);

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
            TIntHashSet incomparables = new TIntHashSet();
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

            //System.out.println("F-matrix: \n"+At);

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

            //System.out.println("Solution: v -\n"+vn+"\n ker -\n"+K);

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
        //System.out.println("Bar = "+bar._1());
        return bar;
    }

}

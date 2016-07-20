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

import topcat.persistence.fcf.FeatureCountingFunction;
import topcat.persistence.functor.Functor;
import topcat.util.GridIterator;
import topcat.util.IntTuple;
import topcat.util.Pair;
import topcat.util.Tuple;

import java.util.List;

/**
 * Represents a type of domain noise described in the paper Multidimensional Persistence and Noise by Chachólski et al.
 * (arXiv:1505.06929).
 */
public class DomainNoise extends Noise{
    //private static final Logger log = LoggerFactory.getLogger(DomainNoise.class);

    public static FeatureCountingFunction computeFCF(Functor F, List<List<Double>> filtrationValues, double[] weight){
        Noise noise = new DomainNoise();
        //TODO: Turn weight into epsilon.
        return noise.computeFCF(F, filtrationValues);
    }

    public FeatureCountingFunction computeFCF(Functor F, List<List<Double>> filtrationValues, Tuple<Double> weight){
        FeatureCountingFunction featureCountingFunction = new FeatureCountingFunction();
        List<IntTuple> indices = getIndexSequence(filtrationValues, weight);
        List<Functor.Generator> f_generators = F.getGenerators();

        featureCountingFunction.add(new Pair<>(0.0, f_generators.size()));

        for(IntTuple epsilon : indices){
            Functor H = new Functor(F.getSize().minus(epsilon));
            for(IntTuple v : GridIterator.getSequence(H.getSize())){
                for(int i=0;i<v.length();i++){
                    H.setMap(v, F.getMap(v.plus(epsilon), i), i);
                }
            }
            List<Functor.Generator> generators = H.getGenerators();
            Integer bar = generators.size();

            double max = Double.MIN_VALUE;
            for(int i=0;i<epsilon.length();i++){
                if(filtrationValues.get(i).get(epsilon.get(i)) > max){
                    max = filtrationValues.get(i).get(epsilon.get(i));
                }
            }

            featureCountingFunction.add(new Pair<>(max, bar));
        }

        System.out.println(featureCountingFunction);

        return featureCountingFunction;
    }
}

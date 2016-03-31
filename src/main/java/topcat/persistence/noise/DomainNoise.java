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

import topcat.persistence.barcode.BasicBarcode;
import topcat.persistence.functor.Functor;
import topcat.util.GridIterator;
import topcat.util.IntTuple;
import topcat.util.Pair;

import java.util.List;

/**
 * Represents a type of domain noise described in the paper Multidimensional Persistence and Noise by Chachólski et al.
 * (arXiv:1505.06929).
 */
public class DomainNoise extends Noise{
    //private static final Logger log = LoggerFactory.getLogger(DomainNoise.class);

    public BasicBarcode computeBasicBarcode(Functor F, List<List<Double>> filtrationValues){
        BasicBarcode basicBarcode = new BasicBarcode();
        List<IntTuple> indices = getIndexSequence(filtrationValues);
        List<Functor.Generator> f_generators = F.getGenerators();

        basicBarcode.add(new Pair<>(0.0, f_generators.size()));

        for(IntTuple epsilon : indices){
            Functor H = new Functor(F.getSize().minus(epsilon));
            GridIterator.getSequence(H.getSize()).parallelStream().forEach(v ->{
                for(int i=0;i<v.length();i++){
                    H.setMap(v, F.getMap(v.plus(epsilon), i), i);
                }
            });
            List<Functor.Generator> generators = H.getGenerators();
            Integer bar = generators.size();

            double max = Double.MIN_VALUE;
            for(int i=0;i<epsilon.length();i++){
                if(filtrationValues.get(i).get(epsilon.get(i)) > max){
                    max = filtrationValues.get(i).get(epsilon.get(i));
                }
            }

            basicBarcode.add(new Pair<>(max, bar));
        }

        System.out.println(basicBarcode);

        return basicBarcode;
    }
}

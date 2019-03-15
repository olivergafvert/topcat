/*
Topcat - a multidimensional persistent homology library.
Copyright (C) 2019 Oliver Gäfvert

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

package topcat.persistence.landscape;

import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;
import topcat.matrix.PMatrix;
import topcat.matrix.distancematrix.DistanceMatrix;
import topcat.persistence.PersistenceModule;
import topcat.util.Grid;
import topcat.util.GridIterator;
import topcat.util.IntTuple;

import java.util.ArrayList;
import java.util.List;

public class PersistenceLandscape {
    IntTuple size;
    Int2ObjectOpenHashMap<Grid<Double>> lambdas = new Int2ObjectOpenHashMap<>();

    public PersistenceLandscape(IntTuple size) {
        this.size = size;
    }

    public void set(int k, IntTuple v, Double value){
        if(!lambdas.containsKey(k)) lambdas.put(k, Grid.create(size));
        lambdas.get(k).set(v, value);
    }

    public double get(int k, IntTuple v){
        if(!lambdas.containsKey(k) || lambdas.get(k).get(v) == null) return 0;
        return lambdas.get(k).get(v);
    }

    public List<List<Double>> vectorize(){
        List<List<Double>> v_lambdas = new ArrayList<>();
        for(int i=0;i<lambdas.keySet().size();i++) v_lambdas.add(new ArrayList<>());
        for(IntTuple v : GridIterator.getSequence(size)){
            for(int k : lambdas.keySet()){
                v_lambdas.get(k).add(get(k, v));
            }
        }
        return v_lambdas;
    }

    public static PersistenceLandscape cartesian(PersistenceModule P, Integer p){
        PersistenceLandscape landscape = new PersistenceLandscape(P.getFunctor().getSize());
        for(IntTuple v : GridIterator.getSequence(P.getFunctor().getSize())){
            Int2ObjectOpenHashMap<List<Double>> lambdas = new Int2ObjectOpenHashMap<>();
            List<IntTuple> basis = IntTuple.getStandardBasisSequence(v.length());
            for (int i=0;i<basis.size();i++) {
                IntTuple shift = basis.get(i);
                int j=0;
                int k=P.getFunctor().getDimension(v);
                while(k>=0) {
                    while (PMatrix.rank(P.getFunctor().getMap(v, shift)) >= k && shift.lt(P.getFunctor().getSize())) {
                        shift = shift.plus(basis.get(i));
                        j++;
                    }
                    if(!lambdas.containsKey(k)) lambdas.put(k, new ArrayList<>());
                    lambdas.get(k).add(P.getFiltrationValues().get(i).get(v.get(i) + j) - P.getFiltrationValues().get(i).get(v.get(i)));
                    k--;
                }
            }
            for(int k : lambdas.keySet()) {
                landscape.set(k, v, DistanceMatrix.norm(lambdas.get(k), p));
            }
        }
        return landscape;
    }

}

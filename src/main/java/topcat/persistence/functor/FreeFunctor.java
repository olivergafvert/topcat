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

package topcat.persistence.functor;

import topcat.matrix.BMatrix;
import topcat.persistence.homology.HomologyUtil;
import topcat.util.IntTuple;

import java.util.ArrayList;
import java.util.List;

public class FreeFunctor extends Functor {

    public List<Generator> generators;

    FreeFunctor(List<Generator> generators, IntTuple size){
        this.generators = generators;
        this.size = size;
    }

    public List<IntTuple> betti0(){
        List<IntTuple> gens = new ArrayList<>();
        for(int i=0;i<generators.size();i++){
            gens.add(generators.get(i).position);
        }
        return gens;
    }

    @Override
    public List<Generator> getGenerators(){
        return generators;
    }

    @Override
    public BMatrix getMap(IntTuple u, IntTuple v){
        List<Generator> u_gens = Generator.getGensLEQThan(generators, u);
        List<Generator> v_gens = Generator.getGensLEQThan(generators, v);
        if(u_gens.equals(v_gens)) return BMatrix.identity(u_gens.size());
        return HomologyUtil.computeInclusionMap(u_gens, v_gens);
    }

    @Override
    public BMatrix getMap(IntTuple u, int dim){
        IntTuple v = u.plus(IntTuple.getStandardBasisElement(u.length(), dim));
        return getMap(u, v);
    }

    @Override
    public int getDimension(IntTuple v){
        if(v.hasNegativeElements(v)){
            return 0;
        }
        return Generator.getGensLEQThan(generators, v).size();
    }




}

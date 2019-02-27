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

package topcat.persistence.simplex;

import topcat.util.BinomialCoeffTable;


/**
 * Enumerates the simplices in the coboundary of a simplex. This is inspired by the combinatorial
 * number system and coboundary enumeration in Ripser [1].
 *
 * [1] - http://ripser.org
 */
public class SimplexCoboundaryEnumerator {

    long idx_below, idx_above;
    int v, k;
    BinomialCoeffTable binomial_coeff;

    public SimplexCoboundaryEnumerator(Simplex simplex, int n_vertices, BinomialCoeffTable binomial_coeff){
        this.idx_above=0;
        this.idx_below=simplex.getIndex();
        this.binomial_coeff = binomial_coeff;
        this.k = simplex.getDimension()+1;
        this.v = n_vertices;
    }

    public boolean hasNext(){
        while((v!=-1) && binomial_coeff.get(v, k) <= idx_below){
            idx_below -= binomial_coeff.get(v, k);
            idx_above += binomial_coeff.get(v, k+1);
            --v;
            --k;
            assert(k!=-1);
        }
        return v!=-1;
    }

    public long next(){
        long index = idx_above + binomial_coeff.get(v--, k+1) + idx_below;
        return index;
    }
}

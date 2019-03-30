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

package topcat.persistence.simplex;

import topcat.util.BinomialCoeffTable;

import java.util.ArrayList;
import java.util.List;

/**
 * Represents an n-simplex for some n >= 0. The implementation is inspired by the combinatorial number system
 * used in Ripser [1].
 *
 * [1] - http://ripser.org
 */
public class Simplex implements Comparable<Simplex>{
    private final long index;
    private final int dimension;

    public Simplex(long index, int dimension){
        this.index = index;
        this.dimension = dimension;
    }

    public long getIndex(){ return index; }

    public int getDimension() {
        return dimension;
    }

    public List<Simplex> getBoundary(SimplexStorageStructure simplexStorageStructure){
        int[] vertices = get_simplex_vertices(getIndex(), getDimension(), simplexStorageStructure.getNumberOfVertices(), simplexStorageStructure.binomialCoeffTable);
        List<Simplex> boundary = new ArrayList<>();
        for(int i=0;i<vertices.length;i++){
            List<Integer> new_v = new ArrayList<>();
            for(int k=0;k<vertices.length;k++){
                if(k!=i){
                    new_v.add(vertices[k]);
                }
            }
            boundary.add(simplexStorageStructure.indexLookup(simplexStorageStructure.binomialCoeffTable.computeIndex(new_v), getDimension()-1));
        }
        return boundary;
    }

    @Override
    public boolean equals(Object o){
        if(o == null || !o.getClass().equals(this.getClass())){
            return false;
        }
        Simplex simplex = (Simplex) o;
        return this.dimension == simplex.dimension && this.index == simplex.index;
    }

    @Override
    public int hashCode(){
        return (int)index;
    }

    @Override
    public int compareTo(Simplex o) {
        if(this.index < o.index) return -1;
        if(this.index > o.index) return 1;
        return 0;
    }

    /**
     * Enumerates the vertices of a simplex of dimension 'dim' with index 'idx' given a total number
     * of vertices 'v'.
     * @param idx - index of simplex
     * @param dim - dimension of simplex
     * @param v - total number of vertices
     * @return an int array with the vertices of the simplex.
     */
    public static int[] get_simplex_vertices(long idx, int dim, int v, BinomialCoeffTable binomial_coeff){
        int[] vertices = new int[dim+1];
        for (int k = dim + 1; k > 0; --k) {
            v = get_next_vertex(v, idx, k, binomial_coeff);
            vertices[k-1] = v;
            idx -= binomial_coeff.get(v, k);
        }
        return vertices;
    }

    /**
     * Performs a binary search to find the largest vertex in a (k-1)-simplex.
     * @param v - total number of vertices
     * @param idx - the index of the (k-1)-simplex
     * @param k
     * @return
     */
    static int get_next_vertex(int v, long idx, int k, BinomialCoeffTable binomial_coeff) {
        if (binomial_coeff.get(v, k) > idx) {
            int count = v;
            while (count > 0) {
                int i = v;
                int step = count >> 1;
                i -= step;
                if (binomial_coeff.get(i, k) > idx) {
                    v = --i;
                    count -= step + 1;
                } else
                    count = step;
            }
        }
        assert(binomial_coeff.get(v, k) <= idx && binomial_coeff.get(v + 1, k) > idx);
        return v;
    }
}

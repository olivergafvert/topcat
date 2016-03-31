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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Represents an n-simplex for some n >= 0.
 */
public class Simplex implements Comparable<Simplex>{
    private final int dimension;
    private final List<Integer> vertices;
    private final int hash;

    public Simplex(List<Integer> vertices){
        this.dimension = vertices.size()-1;
        this.vertices = vertices;
        this.hash = vertices.hashCode();
    }

    public Simplex(Integer... vertices){
        this(Arrays.asList(vertices));
    }

    public Simplex(int[] vertices_arr){
        this.vertices = new ArrayList<>(vertices_arr.length);
        for(int v : vertices_arr){
            this.vertices.add(v);
        }
        this.dimension = this.vertices.size()-1;
        this.hash = vertices.hashCode();
    }

    public int getDimension() {
        return dimension;
    }

    public List<Integer> getVertices() {
        return vertices;
    }

    public int getHash(){
        return hash;
    }

    /**
     * Returns a list of the boundary simplices (i.e (n-1)-dimensional faces).
     * @return
     */
    public List<Simplex> getBoundaryList(){
        List<Simplex> simplices = new ArrayList<>(dimension*dimension);
        for(int i=0;i<vertices.size();i++){
            List<Integer> boundaryVertices = new ArrayList<>(dimension-1);
            for(int j=0;j<vertices.size();j++){
                if(i!=j){
                    boundaryVertices.add(vertices.get(j));
                }
            }
            simplices.add(new Simplex(boundaryVertices));
        }
        return simplices;
    }

    public String toString(){
        StringBuilder sb = new StringBuilder();
        return sb.append(dimension).append("-Simplex: vertices: ").append(vertices).toString();
    }

    @Override
    public boolean equals(Object o){
        if(o == null || !o.getClass().equals(this.getClass())){
            return false;
        }
        Simplex simplex = (Simplex) o;
        if(this.dimension != simplex.dimension){
            return false;
        }
        if(this.hash != simplex.hash){
            return false;
        }
        if(simplex.vertices.equals(this.vertices)){
            return true;
        }
        return false;
    }

    @Override
    public int hashCode(){
        return hash;
    }

    @Override
    public int compareTo(Simplex o) {
        for(int i=0;i<vertices.size();i++) {
            if(!vertices.get(i).equals(o.vertices.get(i))) {
                return vertices.get(i).compareTo(o.vertices.get(i));
            }
        }
        return 0;
    }
}

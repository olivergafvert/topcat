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

import it.unimi.dsi.fastutil.ints.Int2IntOpenHashMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.ints.IntIterator;
import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import it.unimi.dsi.fastutil.longs.Long2IntOpenHashMap;
import it.unimi.dsi.fastutil.longs.Long2ObjectOpenHashMap;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import topcat.util.*;
import topcat.util.paralelliterator.ParalellIterator;

import java.io.*;
import java.util.*;

/**
 * Represents a grid storage structure based on a hierarchy of hashmaps for
 * storing a multifiltered simplicial complex.
 */
public class SimplexStorageStructure {
    private static Logger log = LoggerFactory.getLogger(SimplexStorageStructure.class);
    Int2ObjectOpenHashMap<Grid<List<Simplex>>> simplexContainer;
    List<Long2ObjectOpenHashMap<Simplex>> index_lookup;
    List<List<Double>> filtrationValues;
    IntTuple gridSize;
    public BinomialCoeffTable binomialCoeffTable;
    Integer n_vertices;
    int maxDimension;

    public SimplexStorageStructure(List<List<Double>> filtrationValues, IntTuple gridSize, Integer maxDimension, Integer n_vertices){
        simplexContainer = new Int2ObjectOpenHashMap<>();
        this.filtrationValues = filtrationValues;
        this.gridSize = gridSize;
        this.n_vertices = n_vertices;
        this.maxDimension = maxDimension;
        this.binomialCoeffTable = new BinomialCoeffTable(n_vertices, maxDimension);
        this.index_lookup = new ArrayList<>();
        for(int i=0;i<=maxDimension;i++){
            index_lookup.add(new Long2ObjectOpenHashMap<>());
        }
    }

    public Simplex indexLookup(Long index, Integer dim){
        return index_lookup.get(dim).get((long)index);
    }

    public List<List<Double>> getFiltrationValues() { return filtrationValues; }

    public void addElement(Simplex simplex, IntTuple filtrationIndex){
        if(!simplexContainer.containsKey(simplex.getDimension())){
            simplexContainer.put(simplex.getDimension(), Grid.create(gridSize));
        }
        List<Simplex> simplices = simplexContainer.get(simplex.getDimension()).get(filtrationIndex);
        if(simplices == null){
            simplices = new ArrayList<>();
            simplexContainer.get(simplex.getDimension()).set(filtrationIndex, simplices);
        }
        index_lookup.get(simplex.getDimension()).put(simplex.getIndex(), simplex);
        simplices.add(simplex);
    }

    public List<Simplex> getSimplicesAt(int dim, IntTuple filtrationIndex){
        if(!simplexContainer.containsKey(dim)){
            return new ArrayList<>();
        }
        List<Simplex> simplices = simplexContainer.get(dim).get(filtrationIndex);
        if(simplices == null) return new ArrayList<>();
        return simplices;
    }

    public List<Simplex> getSimplicesLEQThan(int dim, IntTuple filtrationIndex){
        List<Simplex> simplices = new ArrayList<>();
        for(IntTuple v : GridIterator.getSequence(filtrationIndex)){
            List<Simplex> local_simplices = getSimplicesAt(dim, v);
            if(local_simplices != null) {
                simplices.addAll(local_simplices);
            }
        }
        Collections.sort(simplices);
        return simplices;
    }

    public List<IntTuple> gridSequence(int dim){
        return gridSequence();
    }

    public List<IntTuple> gridSequence(){
        return GridIterator.getSequence(gridSize);
    }

    public List<IntTuple> getAdjacent(IntTuple v){
        List<IntTuple> basis = IntTuple.getStandardBasisSequence(v.length());
        List<IntTuple> ret = new ArrayList<>();
        for(int i=0;i<basis.size();i++){
            IntTuple w = v.minus(basis.get(i));
            if(w.min()>=0)
                ret.add(v.minus(basis.get(i)));
        }
        return ret;
    }

    public int getNumberOfVertices(){
        return n_vertices;
    }

    /**
     * Returns the simplex storage structure and filtration values.
     * @param f
     * @return
     */
    public static SimplexStorageStructure readFromFile(File f) throws IOException{
        BufferedReader reader = new BufferedReader(new FileReader(f));
        String line;
        Integer dimension = Integer.parseInt(reader.readLine());
        List<List<Double>> filtrationValues = new ArrayList<>();
        for(int i=0;i<dimension;i++){
            List<Double> values = new ArrayList<>();
            for(String x : reader.readLine().trim().split(" ")){
                values.add(Double.parseDouble(x));
            }
            filtrationValues.add(values);
        }
        IntTuple gridSize = IntTuple.zeros(filtrationValues.size());
        for(int i=0;i<filtrationValues.size();i++){
            gridSize.set(i, filtrationValues.get(i).size()-1);
        }
        Integer n_vertices = Integer.parseInt(reader.readLine());
        Integer maxDimension = Integer.parseInt(reader.readLine());
        SimplexStorageStructure simplexStorageStructure = new SimplexStorageStructure(filtrationValues, gridSize, maxDimension, n_vertices);
        while((line = reader.readLine()) != null){

            //Parse simplex
            String[] part = line.trim().split(":");
            if(part.length != 2){
                throw new IOException("Malformed Simplex format. Line: '"+line+"' contains more or less than one ':'.");
            }

            //Parse vertices
            List<Integer> vertices = new ArrayList<>();
            for(String v : part[0].trim().split(" ")){
                vertices.add(Integer.parseInt(v));
            }
            long index = simplexStorageStructure.binomialCoeffTable.computeIndex(vertices);

            //Parse filtration index
            List<Integer> filtrationIndex = new ArrayList<>();
            for(String x : part[1].trim().split(" ")){
                filtrationIndex.add(Integer.parseInt(x));
            }

            simplexStorageStructure.addElement(new Simplex(index, vertices.size()-1), new IntTuple(filtrationIndex));
        }
        return simplexStorageStructure;
    }

    public void clearSimplices(){
        simplexContainer.clear();
        index_lookup.clear();
    }

    @Override
    public String toString(){
        IntTuple size = IntTuple.zeros(filtrationValues.size());
        for(int i=0;i<filtrationValues.size();i++){
            size.set(i, filtrationValues.get(i).size());
        }
        StringBuilder sb = new StringBuilder();
        for(IntTuple v : GridIterator.getSequence(size)){
            sb.append(v).append(":\n");
            sb.append("Dimension 0: ").append(getSimplicesAt(0, v)).append("\n");
            sb.append("Dimension 1; ").append(getSimplicesAt(1, v)).append("\n");
            sb.append("Dimension 2; ").append(getSimplicesAt(2, v)).append("\n");
        }
        return sb.toString();
    }
}

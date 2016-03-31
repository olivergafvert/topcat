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

import gnu.trove.map.hash.TIntObjectHashMap;
import topcat.util.GridIterator;
import topcat.util.IntTuple;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * Represents a grid storage structure based on a hierarchy of hashmaps for
 * storing a multifiltered simplicial complex.
 */
public class SimplexStorageStructure {
    TIntObjectHashMap<Container> simplexContainer;
    List<List<Double>> filtrationValues;

    public SimplexStorageStructure(List<List<Double>> filtrationValues){
        simplexContainer = new TIntObjectHashMap<>();
        this.filtrationValues = filtrationValues;
    }

    public List<List<Double>> getFiltrationValues() { return filtrationValues; }

    public void addElement(Simplex simplex, List<Integer> filtrationIndex){
        if(!simplexContainer.containsKey(simplex.getDimension())){
            if(this.filtrationValues.size() > 1)
                simplexContainer.put(simplex.getDimension(), new Container());
            else
                simplexContainer.put(simplex.getDimension(), new SimplexContainer());
        }
        simplexContainer.get(simplex.getDimension()).set(filtrationIndex, simplex);
    }

    public List<Simplex> getSimplicesAt(int dim, List<Integer> filtrationIndex){
        if(!simplexContainer.containsKey(dim)) return new ArrayList<>();
        return simplexContainer.get(dim).get(filtrationIndex);
    }

    private void getSimplicesLEQThanHelper(int dim, List<Integer> filtrationIndex, int pos, List<Integer> current, List<Simplex> simplices){
        if(pos == filtrationIndex.size()){
            List<Simplex> simplexList = getSimplicesAt(dim, current);
            if(simplexList != null)
                simplices.addAll(simplexList);
        }else {
            for (int i = 0; i <= filtrationIndex.get(pos); i++) {
                current.set(pos, i);
                getSimplicesLEQThanHelper(dim, filtrationIndex, pos + 1, current, simplices);
            }
        }
    }

    public List<Simplex> getSimplicesLEQThan(int dim, List<Integer> filtrationIndex){
        List<Simplex> simplices = new ArrayList<>();
        List<Integer> positionVector = new ArrayList<>();
        for(int i=0;i<filtrationIndex.size();i++) positionVector.add(0);
        getSimplicesLEQThanHelper(dim, filtrationIndex, 0, positionVector, simplices);
        return simplices;
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
            Arrays.asList(reader.readLine().trim().split(" "))
                    .stream()
                    .forEach(x -> values.add(Double.parseDouble(x)));
            filtrationValues.add(values);
        }
        SimplexStorageStructure simplexStorageStructure = new SimplexStorageStructure(filtrationValues);
        while((line = reader.readLine()) != null){

            //Parse simplex
            String[] part = line.trim().split(":");
            if(part.length != 2){
                throw new IOException("Malformed Simplex format. Line: '"+line+"' contains more or less than one ':'.");
            }

            //Parse vertices
            List<Integer> vertices = new ArrayList<>();
            Arrays.asList(part[0].trim().split(" "))
                    .stream()
                    .forEach(v -> vertices.add(Integer.parseInt(v)));
            Collections.sort(vertices);

            //Parse filtration index
            List<Integer> filtrationIndex = new ArrayList<>();
            Arrays.asList(part[1].trim().split(" "))
                    .stream()
                    .forEach(x -> filtrationIndex.add(Integer.parseInt(x)));

            simplexStorageStructure.addElement(new Simplex(vertices), filtrationIndex);
        }
        return simplexStorageStructure;
    }

    @Override
    public String toString(){
        IntTuple size = new IntTuple(filtrationValues.size());
        for(int i=0;i<filtrationValues.size();i++){
            size.set(i, filtrationValues.get(i).size());
        }
        StringBuilder sb = new StringBuilder();
        GridIterator.getSequence(size).stream().forEach(v -> {
            sb.append(v).append(":\n");
            sb.append("Dimension 0: ").append(getSimplicesAt(0, v.toList())).append("\n");
            sb.append("Dimension 1; ").append(getSimplicesAt(1, v.toList())).append("\n");
            sb.append("Dimension 2; ").append(getSimplicesAt(2, v.toList())).append("\n");
        });
        return sb.toString();
    }

    class Container{
        private TIntObjectHashMap<Container> container = new TIntObjectHashMap<>();

        List<Simplex> get(List<Integer> vals){
            if(vals.size()>2){
                if(!container.containsKey(vals.get(0))) return null;
                return container.get(vals.get(0))
                        .get(vals.subList(1, vals.size()));
            }
            if(!container.containsKey(vals.get(1))) return null;
            return ((SimplexContainer)container.get(vals.get(1))).get(vals.get(0));
        }

        void set(List<Integer> vals, Simplex simplex){
            set(vals, Arrays.asList(new Simplex[]{simplex}));
        }

        void set(List<Integer> vals, List<Simplex> simplices){
            if(vals.size() > 2) {
                if (!container.containsKey(vals.get(0))) container.put(vals.get(0), new Container());
                container.get(vals.get(0)).set(vals.subList(1, vals.size()), simplices);
            }else{
                if (!container.containsKey(vals.get(0))) container.put(vals.get(0), new SimplexContainer());
                SimplexContainer simplexContainer = (SimplexContainer) container.get(vals.get(0));
                simplexContainer.set(vals.get(1), simplices);
            }
        }
    }

    class SimplexContainer extends Container{
        private TIntObjectHashMap<List<Simplex>> simplexContainer = new TIntObjectHashMap<>();

        List<Simplex> get(int val){
            if(!simplexContainer.containsKey(val)) return null;
            return simplexContainer.get(val);
        }

        @Override
        List<Simplex> get(List<Integer> vals){
            return get(vals.get(0));
        }

        void set(int val, List<Simplex> simplices){
            if(!simplexContainer.containsKey(val)) simplexContainer.put(val, new ArrayList<>());
            simplexContainer.get(val).addAll(simplices);
        }

        @Override
        void set(List<Integer> vals, List<Simplex> simplices){
            set(vals.get(0), simplices);
        }

    }
}

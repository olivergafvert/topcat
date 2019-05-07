package topcat.persistence.simplex;

import it.unimi.dsi.fastutil.ints.Int2IntOpenHashMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.ints.IntIterator;
import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import it.unimi.dsi.fastutil.longs.Long2IntOpenHashMap;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import topcat.util.BinomialCoeffTable;
import topcat.util.Grid;
import topcat.util.IntTuple;
import topcat.util.Pair;
import topcat.util.paralelliterator.ParalellIntIterator;
import topcat.util.paralelliterator.ParalellIterator;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;
import java.util.logging.Level;

public class SparseSimplexStorageStructure extends SimplexStorageStructure{
    private static Logger log = LoggerFactory.getLogger(SparseSimplexStorageStructure.class);

    List<HashMap<IntTuple, List<Simplex>>> simplices = new ArrayList<>();
    HashMap<IntTuple, List<IntTuple>> adjacent_map = new HashMap<>();
    List<IntTuple> grid_sequence = null;

    public SparseSimplexStorageStructure(List<List<Double>> filtrationValues, IntTuple gridSize, Integer max_dimesion, Integer n_vertices) {
        super(filtrationValues, gridSize, max_dimesion, n_vertices);
        for(int i=0;i<=max_dimesion;i++)
            simplices.add(new HashMap<>());
    }

    public void addElement(Simplex simplex, IntTuple filtrationIndex){
        if(!simplices.get(simplex.getDimension()).containsKey(filtrationIndex))
            simplices.get(simplex.getDimension()).put(filtrationIndex, new ArrayList<>());
        simplices.get(simplex.getDimension()).get(filtrationIndex).add(simplex);
    }

    public List<Simplex> getSimplicesAt(int dim, IntTuple filtrationIndex){
        if(dim > this.maxDimension || !simplices.get(dim).containsKey(filtrationIndex)) return new ArrayList<>();
        return simplices.get(dim).get(filtrationIndex);
    }

    public List<IntTuple> gridSequence(int dim){
        return new ArrayList<>(simplices.get(dim).keySet());
    }

    public List<IntTuple> gridSequence(){
        if(this.grid_sequence != null) return grid_sequence;
        HashSet<IntTuple> keys = new HashSet<>();
        for(int i=0;i<=maxDimension;i++)
            keys.addAll(simplices.get(i).keySet());
        return new ArrayList<>(keys);
    }


    public List<Pair<IntTuple, List<IntTuple>>> paretoFrontiers(final List<IntTuple> elements){
        Collections.sort(elements, new Comparator<IntTuple>() {
            @Override
            public int compare(IntTuple o1, IntTuple o2) {
                if(o1.lt(o2)) return -1;
                if(o2.lt(o1)) return 1;
                if(o1.lexLt(o2)) return -1;
                if(o2.lexLt(o1)) return 1;
                return 0;
            }
        });

        final List<Pair<Integer, Pair<IntTuple, List<IntTuple>>>> maximalElements = new ParalellIntIterator<Pair<IntTuple, List<IntTuple>>>(0, elements.size()) {
            @Override
            public Pair<IntTuple, List<IntTuple>> method(Integer index) {
                List<IntTuple> maxElements = new ArrayList<>();
                IntTuple z = elements.get(index);
                for(int i=index-1;i>=0;i--){
                    IntTuple v = elements.get(i);
                    if(v.leq(z)) {
                        boolean shouldAdd = true;
                        for (IntTuple w : maxElements) {
                            if (v.leq(w))
                                shouldAdd = false;
                        }
                        if (shouldAdd)
                            maxElements.add(v);
                    }
                }
                return new Pair<>(z, maxElements);
            }
        }.run();
        List<Pair<IntTuple, List<IntTuple>>> ret = new ArrayList<>();
        for(Pair<Integer, Pair<IntTuple, List<IntTuple>>> pair : maximalElements)
            ret.add(pair._2());
        return ret;
    }


    public void preparePoset1(){
        final List<IntTuple> index_seq = gridSequence();
        log.debug("Number of indices: "+index_seq.size());
        Collections.sort(index_seq, new Comparator<IntTuple>() {
            @Override
            public int compare(IntTuple o1, IntTuple o2) {
                if(o1.lt(o2)) return -1;
                if(o2.lt(o1)) return 1;
                if(o1.lexLt(o2)) return -1;
                if(o2.lexLt(o1)) return 1;
                return 0;
            }
        });

        long time_s = System.nanoTime();

//        Int2ObjectOpenHashMap<IntOpenHashSet> leqmap = new Int2ObjectOpenHashMap<>();
//
//        for(int i=0;i<index_seq.size();i++){
//            IntTuple z = index_seq.get(i);
//            List<IntTuple> maxElements = new ArrayList<>();
//            IntOpenHashSet visited = new IntOpenHashSet();
//            leqmap.put(i, visited);
//            for(int j=i-1;j>=0;j--){
//                if(!visited.contains(j) && index_seq.get(j).leq(z)){
//                    visited.add(j);
//                    visited.addAll(leqmap.get(j));
//                    maxElements.add(index_seq.get(j));
//                }
//            }
//        }


        List<Pair<Integer, List<IntTuple>>> maximalElements = new ParalellIntIterator<List<IntTuple>>(0, index_seq.size()) {
            @Override
            public List<IntTuple> method(Integer index) {
                List<IntTuple> maxElements = new ArrayList<>();
                IntTuple z = index_seq.get(index);
                for(int i=index-1;i>=0;i--){
                    IntTuple v = index_seq.get(i);
                    if(v.leq(z)) {
                        boolean shouldAdd = true;
                        for (IntTuple w : maxElements) {
                            if (v.leq(w))
                                shouldAdd = false;
                        }
                        if (shouldAdd)
                            maxElements.add(v);
                    }
                }
                return maxElements;
            }
        }.run();

        for(Pair<Integer, List<IntTuple>> pair : maximalElements){
            adjacent_map.put(index_seq.get(pair._1()), pair._2());
        }

        log.debug("Time: "+(System.nanoTime()-time_s));
    }


    public void preparePoset(){
        final List<IntTuple> index_seq = gridSequence();
        log.debug("Starting to prepare poset with "+index_seq.size()+" number of elements...");
        final List<Pair<IntTuple, List<IntTuple>>> maximalElements = paretoFrontiers(index_seq);

        //Saturate poset
        List<Pair<Pair<IntTuple, List<IntTuple>>, List<Pair<IntTuple, List<IntTuple>>>>> expandedElements = new ParalellIterator<Pair<IntTuple, List<IntTuple>>, List<Pair<IntTuple, List<IntTuple>>>>(maximalElements) {
            @Override
            public List<Pair<IntTuple, List<IntTuple>>> method(Pair<IntTuple, List<IntTuple>> z) {
                List<Pair<IntTuple, List<IntTuple>>> expanded = new ArrayList<>();
                if(z._2().size()<=2){
                    expanded.add(z);
                    return expanded;
                }
                HashSet<IntTuple> visited = new HashSet<>();
                HashSet<IntTuple> isbasis = new HashSet<>();
                visited.add(z._1());
                visited.addAll(z._2());
                isbasis.addAll(z._2());
                List<IntTuple> baseElements = z._2();
                List<IntTuple> syzygies = new ArrayList<>();
                for(int i=0;i<baseElements.size();i++){
                    for(int j=i+1;j<baseElements.size();j++){
                        IntTuple join = IntTuple.join(baseElements.get(i), baseElements.get(j));
                        if(!visited.contains(join)){
                            visited.add(join);
                            syzygies.add(join);
                        }
                    }
                }
                baseElements.addAll(syzygies);
                baseElements.add(z._1());
                List<Pair<IntTuple, List<IntTuple>>> paretofront = paretoFrontiers(baseElements);
                List<Pair<IntTuple, List<IntTuple>>> ret = new ArrayList<>();
                for(int i=0;i<paretofront.size();i++) {
                    if(!isbasis.contains(paretofront.get(i)._1())) {
                        ret.add(new Pair<>(paretofront.get(i)._1(), paretofront.get(i)._2()));
                    }
                }
                return ret;
            }

        }.run();

        this.grid_sequence = new ArrayList<>();
        for(Pair<Pair<IntTuple, List<IntTuple>>, List<Pair<IntTuple, List<IntTuple>>>> pair : expandedElements){
            for(Pair<IntTuple, List<IntTuple>> element : pair._2()) {
                adjacent_map.put(element._1(), element._2());
                grid_sequence.add(element._1());
            }
        }
        Collections.sort(grid_sequence, new Comparator<IntTuple>() {
            @Override
            public int compare(IntTuple o1, IntTuple o2) {
                if(o1.lt(o2)) return -1;
                if(o2.lt(o1)) return 1;
                if(o1.lexLt(o2)) return -1;
                if(o2.lexLt(o1)) return 1;
                return 0;
            }
        });
        log.debug("Finished preparing poset with a resulting number of elements: "+grid_sequence.size());
    }


    private List<IntTuple> getLt(IntTuple w, List<IntTuple> indices){
        List<IntTuple> leqs = new ArrayList<>();
        for(IntTuple v : indices){
            if(v.lt(w))
                leqs.add(v);
        }
        return leqs;
    }

    private List<IntTuple> findMaximalElements(List<IntTuple> indices){
        IntOpenHashSet visited = new IntOpenHashSet();
        List<IntTuple> maximal = new ArrayList<>();
        for(int i=0;i<indices.size();i++){
            if(!visited.contains(i)) {
                visited.add(i);
                int k = i;
                boolean shouldAdd = true;
                for(int j=i+1;j<indices.size();j++){
                    if(indices.get(k).leq(indices.get(j))){
                        if(visited.contains(j)){
                            shouldAdd=false;
                            break;
                        }
                        visited.add(j);
                        k=j;
                    }
                }
                if(shouldAdd)
                    maximal.add(indices.get(k));
            }
        }
        return maximal;
    }

    private void max_pairs(List<IntTuple> vs, List<IntTuple> leqs){
        for(IntTuple v : vs){
            if(!adjacent_map.containsKey(v)){
                List<IntTuple> v_leqs = getLt(v, leqs);
                List<IntTuple> v_max = findMaximalElements(v_leqs);
                adjacent_map.put(v, v_max);
                max_pairs(v_max, v_leqs);
            }
        }
    }

    public void prepareAdjacentMap(){
        List<IntTuple> indices = gridSequence();
        Collections.sort(indices, new Comparator<IntTuple>() {
            @Override
            public int compare(IntTuple o1, IntTuple o2) {
                if(o1.lt(o2)) return -1;
                if(o2.lt(o1)) return 1;
                if(o1.lexLt(o2)) return -1;
                if(o2.lexLt(o1)) return 1;
                return 0;
            }
        });
        long time_s = System.nanoTime();
        List<IntTuple> maximal = findMaximalElements(indices);
        max_pairs(maximal, indices);
        log.debug("Time: "+(System.nanoTime()-time_s));
    }


    @Override
    public List<IntTuple> getAdjacent(IntTuple v){
        if(adjacent_map.size()==0){
            preparePoset();
        }
        return adjacent_map.get(v);
    }

    public static SimplexStorageStructure readOFFFile(String file){
        // This will reference one line at a time
        String line;
        SimplexStorageStructure simplexStorageStructure = null;
        Double precision = 1e0;

        try {
            BufferedReader bufferedReader = new BufferedReader(new FileReader(file));

            bufferedReader.readLine();//Discard first line
            line = bufferedReader.readLine();
            String[] s = line.split(" ");
            int NV = Integer.parseInt(s[0]);
            int NF = Integer.parseInt(s[1]);

            List<IntOpenHashSet> filtrationValues = new ArrayList<>();
            filtrationValues.add(new IntOpenHashSet());
            filtrationValues.add(new IntOpenHashSet());

            List<Simplex> simplices = new ArrayList<>();

            for(int i=0;i<NV;i++) {
                Simplex simplex = new Simplex(i, 0);
                line = bufferedReader.readLine();
                s = line.split(" ");
                int x = ((Double)(Double.parseDouble(s[0])*precision)).intValue();
                filtrationValues.get(0).add(x);
                int y = ((Double)(Double.parseDouble(s[1])*precision)).intValue();
                filtrationValues.get(1).add(y);
                simplex.setValue(new IntTuple(x, y));
                simplices.add(simplex);
            }

            List<Int2IntOpenHashMap> filtrationIndices = new ArrayList<>();
            List<List<Double>> d_filtrationValues = new ArrayList<>();
            for(IntOpenHashSet set : filtrationValues) {
                List<Integer> values = new ArrayList<>();
                IntIterator iterator = set.iterator();
                while(iterator.hasNext())
                    values.add(iterator.nextInt());
                Collections.sort(values);
                Int2IntOpenHashMap value_index = new Int2IntOpenHashMap();
                List<Double> d_values = new ArrayList<>();
                for(int i=0;i<values.size();i++) {
                    value_index.put((int) values.get(i), i);
                    d_values.add(((double)values.get(i))/precision);
                }
                filtrationIndices.add(value_index);
                d_filtrationValues.add(d_values);
            }

            new ParalellIterator<Simplex,Void>(simplices) {
                @Override
                public Void method(Simplex index) {
                    IntTuple values = index.getValue();
                    List<Integer> indices = new ArrayList<>();
                    for(int i=0;i<values.length();i++)
                        indices.add(filtrationIndices.get(i).get((int)values.get(i)));
                    index.setValue(new IntTuple(indices));
                    return null;
                }
            }.run();

            List<List<Integer>> triangles = new ArrayList<>();

            for(int i=0;i<NF;i++) {
                line = bufferedReader.readLine();
                s = line.split(" ");
                int n = Integer.parseInt(s[0]);
                List<Integer> vertices = new ArrayList<>();
                for(int j=0;j<n;j++) {
                    vertices.add(Integer.parseInt(s[j+1]));
                }
                triangles.add(vertices);
            }
            bufferedReader.close();


            BinomialCoeffTable binomialCoeffTable = new BinomialCoeffTable(NV, 1);
            Long2IntOpenHashMap edge_map = new Long2IntOpenHashMap();

            int index = NV;
            for(int i=0;i<triangles.size();i++) {
                Simplex simplex = new Simplex(index++, 2);
                simplex.setValue(IntTuple.zeros(filtrationIndices.size()));
                simplices.add(simplex);
                simplex.setBoundary(new ArrayList<>());
                List<Integer> triangle = triangles.get(i);
                for(int m=0;m<3;m++){
                    int v1 = triangle.get(m);
                    int v2 = triangle.get((m+1)%3);
                    if(v2<v1){
                        int temp = v1;
                        v1 = v2;
                        v2 = temp;
                    }
                    long e_index = binomialCoeffTable.computeIndex(v1, v2);
                    if(!edge_map.containsKey(e_index)){
                        edge_map.put(e_index, index);
                        Simplex edge = new Simplex(index++, 1);
                        List<Simplex> boundary = new ArrayList<>();
                        boundary.add(simplices.get(v1));
                        boundary.add(simplices.get(v2));
                        edge.setBoundary(boundary);
                        edge.setValue(IntTuple.join(simplices.get(v1).getValue(), simplices.get(v2).getValue()));
                        simplices.add(edge);
                        simplex.getBoundary(null).add(edge);
                        simplex.setValue(IntTuple.join(simplex.getValue(), edge.getValue()));
                    }else{
                        Simplex edge = simplices.get(edge_map.get(e_index));
                        simplex.getBoundary(null).add(edge);
                        simplex.setValue(IntTuple.join(simplex.getValue(), edge.getValue()));
                    }
                }
            }

            Collections.sort(simplices, new Comparator<Simplex>() {
                @Override
                public int compare(Simplex o1, Simplex o2) {
                    if(o1.getValue().equals(o2.getValue())) {
                        if (o1.getDimension() == o2.getDimension())
                            return (int) (o1.getIndex() - o2.getIndex());
                        return o1.getDimension() - o2.getDimension();
                    }
                    if(o1.getValue().lexLt(o2.getValue())) return -1;
                    return 1;
                }
            });


            simplexStorageStructure = new SparseSimplexStorageStructure(d_filtrationValues, new IntTuple(filtrationIndices.get(0).size(), filtrationIndices.get(1).size()), 2, NV);

            for(int i=0;i<simplices.size();i++){
                simplices.get(i).setIndex(i);
                simplexStorageStructure.addElement(simplices.get(i), simplices.get(i).getValue());
            }
        }
        catch(FileNotFoundException ex) {
            log.error(
                    "Unable to open file '" +
                            file + "'");
        }
        catch(IOException ex) {
            log.error(
                    "Error reading file '"
                            + file + "'");
        }
        return simplexStorageStructure;
    }
}

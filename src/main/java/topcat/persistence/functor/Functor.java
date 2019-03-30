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

package topcat.persistence.functor;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import topcat.matrix.BMatrix;
import topcat.matrix.BVector;
import topcat.matrix.exception.NoSolutionException;
import topcat.matrix.exception.WrongDimensionException;
import topcat.persistence.contours.PersistenceContour;
import topcat.persistence.homology.HomologyUtil;
import topcat.util.*;

import java.util.*;

/**
 * Represents a functor F: N^r -> Vect_K, for some r > 0.
 */
public class Functor {
    private static final Logger log = LoggerFactory.getLogger(Functor.class);

    private List<Grid<BMatrix>> maps; //The maps of the functor
    protected IntTuple size; //The size of the grid of the multifiltration on which the functor is defined

    public Functor(IntTuple size){
        this.size = size;
        this.maps = new ArrayList<>();
        for(int i=0;i<size.length();i++){
            this.maps.add(Grid.create(size));
        }
    }

    public Functor(Integer... d){
        this.size = new IntTuple(d);
        this.maps = new ArrayList<>();
        for(int i=0;i<this.size.length();i++){
            this.maps.add(Grid.create(this.size));
        }
    }

    public Functor(Functor F){
        this.size = new IntTuple(F.getSize());
        this.maps = new ArrayList<>();
        for(Grid<BMatrix> f_grid : F.maps){
            Grid<BMatrix> grid = Grid.create(this.size);
            for(IntTuple v : GridIterator.getSequence(this.size)){
                grid.set(v, new BMatrix(f_grid.get(v)));
            }
            this.maps.add(grid);
        }
    }


    //Private methods

    /**
     * Make sure tuple has indices within the defined domain.
     * @param v
     * @return
     */
    private IntTuple parseTuple(IntTuple v){
        IntTuple w = IntTuple.zeros(v.length());
        for(int i=0;i<v.length();i++){
            w.set(i, v.get(i) <= size.get(i) ? v.get(i) : size.get(i));
        }
        return w;
    }

    /**
     * Returns true if tuple have indices outside of defined domain.
     * @param v
     * @return
     */
    private boolean isOutOfBounds(IntTuple v){
        if(v.hasNegativeElements(v)){
            return true;
        }
        if(!v.equals(parseTuple(v))){
            return true;
        }
        return false;
    }


    //Public methods

    public int getDimension(){
        return this.size.length();
    }

    /**
     * Computes the shift of the generators of F by epsilon.
     * @param f_generators
     * @param epsilon
     * @return
     */
    public List<Generator> generatorShift(List<Generator> f_generators,
                                           Double epsilon, PersistenceContour contour){
        List<Generator> generators = new ArrayList<>(); //Stores the shifted generators

        //First shift all f-generators by epsilon
        for(Generator f : f_generators){
            IntTuple pos = contour.shift(f.position, epsilon);
            if(getMap(f.position, pos)==null){
                int ee= 0;
            }
            BVector g_v = getMap(f.position, pos).mult(f.v);
            if(g_v.getNumberOfNonZeroElements() > 0) {
                generators.add(new Generator(pos, g_v));
            }
        }

        //Now determine which shifted generators are independent
        List<Generator> independent_generators = new ArrayList<>();
        for(int i=0;i<generators.size();i++){
            Generator g = generators.get(i);
            List<BVector> vectors = new ArrayList<>();
            for(int j=0;j<independent_generators.size();j++){
                if(independent_generators.get(j).equals(g)){
                    break;
                }
                BMatrix M = getMap(independent_generators.get(j).position, g.position);
                if(M!=null) {
                    BVector image = M.mult(independent_generators.get(j).v);
                    if(image.getNumberOfNonZeroElements()>0)
                        vectors.add(image);
                }
            }
            if(vectors.size() > 0) {
                BMatrix A = new BMatrix(vectors);
                if(!BMatrix.hasSolution(A.transpose(), g.v)){
                    independent_generators.add(g);
                }
            }else {
                independent_generators.add(g);
            }
        }
        return independent_generators;
    }

    /**
     * Returns the map at position 'v' in direction 'dim'.
     * @param v
     * @param dim
     * @return
     */
    public BMatrix getMap(IntTuple v, int dim){
        if(v.hasNegativeElements(v)){
            return new BMatrix(getDimension(v.plus(IntTuple.getStandardBasisElement(v.length(), dim))), 0);
        }
        return maps.get(dim).get(parseTuple(v));
    }


    /**
     * Sets the map 'A' at position 'v' in direction 'dim'.
     * @param v
     * @param A
     * @param dim
     */
    public void setMap(IntTuple v, BMatrix A, int dim){
        if(isOutOfBounds(v)){
            return;
        }
        maps.get(dim).set(v, A);
    }

    /**
     * Returns the size domain of the functor.
     * @return
     */
    public IntTuple getSize(){
        return size;
    }

    /**
     * Returns the dimension of the vector space at position (i, j).
     * @return
     */
    public int getDimension(IntTuple v){
        if(v.hasNegativeElements(v)){
            return 0;
        }
        return maps.get(0).get(parseTuple(v)).cols;
    }

    /**
     * Returns the map from position 'from' to position 'to'.
     * @param from
     * @param to
     * @return
     */
    public BMatrix getMap(IntTuple from, IntTuple to){
        if(from.equals(to)){
            return BMatrix.identity(getDimension(from));
        }
        if(!from.leq(to)){
            return null;
        }
        BMatrix A = null;
        try {
            for (Pair<IntTuple, Integer> p : GridIterator.getPairSequence(from, to)) {
                if (A == null){
                    A = getMap(p._1(), p._2());
                }
                else {
                    A = getMap(p._1(), p._2()).mult(A);
                }
            }
        }catch (WrongDimensionException wde){
            log.error("Failed to compute map between points "+from+ " and "+to, wde);
        }
        return A;
    }

    /**
     * Returns a list of generators that generate the functor.
     * @return
     */
    public List<Generator> getGenerators(){
        List<Generator> generators = new ArrayList<>();
        for(IntTuple v : GridIterator.getSequence(size)){
            try {
                List<IntTuple> basis = IntTuple.getStandardBasisSequence(v.length());
                BMatrix A = getMap(v.minus(basis.get(0)), v);
                for(int i=1;i<basis.size();i++){
                    A = BMatrix.concat(A, getMap(v.minus(basis.get(i)), v));
                }
                Pair<BMatrix, BMatrix> kerim = BMatrix.reduction(A);
                //System.err.println(kerim._2());
                if (kerim._2().rows != kerim._2().cols) {
                    BMatrix generatorBasis = BMatrix.extendBasis(kerim._2());
                    //System.err.println("Generatorbasis: \n"+generatorBasis);
                    generatorBasis = generatorBasis.subMatrix(kerim._2().rows, -1, 0, -1);
                    for (int k = 0; k < generatorBasis.rows; k++) {
                        generators.add(new Generator(new IntTuple(v), generatorBasis.getRow(k)));
                    }
                }
            }catch (NoSolutionException nse){
                log.info("Filed to compute minimal set of generators.", nse);
            }
        }
        return generators;
    }

    /**
     * Computes the projective resolution of the functor. (OBS Currently only computes up to P_1)
     * TODO: Implement computation of P_2 and higher.
     * @return
     */
    public List<FreeFunctor> projectiveResolution(){
        List<Generator> generators = getGenerators();
        List<FreeFunctor> proj = new ArrayList<>();
        proj.add(new FreeFunctor(generators, size));
        for(int d=0;d<getDimension();d++) {
            Nat kernel = new Nat(size);
            List<Functor.Generator> kernelgens = new ArrayList<>();
            for (IntTuple v : GridIterator.getSequence(getSize())) {
                List<Functor.Generator> gens = Generator.getGensLEQThan(generators, v);
                BMatrix A = new BMatrix(gens.size(), getDimension(v));
                for (int i = 0; i < gens.size(); i++) {
                    A.setRow(i, new BVector(getMap(gens.get(i).position, v).mult(gens.get(i).v)));
                }
                Pair<BMatrix, BMatrix> kerim = BMatrix.reduction(A.transpose());
                if(kerim._1().rows == 0) continue;
                kernel.setMap(v, kerim._1().transpose());
                List<IntTuple> basis = IntTuple.getStandardBasisSequence(v.length());
                List<Generator> ugens = Generator.getGensLEQThan(generators, v.minus(basis.get(0)));
                BMatrix I = HomologyUtil.computeInclusionMap(ugens, gens).mult(kernel.getMap(v.minus(basis.get(0))));
                for(int i=1;i<basis.size();i++){
                    ugens = Generator.getGensLEQThan(generators, v.minus(basis.get(i)));
                    I = BMatrix.concat(I, HomologyUtil.computeInclusionMap(ugens, gens).mult(kernel.getMap(v.minus(basis.get(i)))));
                }
                if(I==null){
                    for (int i = 0; i < kerim._1().rows; i++) {
                        kernelgens.add(new Generator(new IntTuple(v), kerim._1().getRow(i)));
                    }
                }
                else if(BMatrix.rank(I) < kernel.getMap(v).cols){
                    kerim = BMatrix.reduction(I);
                    BMatrix Hbasis = BMatrix.getBasis(BMatrix.vconcat(kerim._2(), kernel.getMap(v).transpose())).subMatrix(kerim._2().rows, -1, 0, -1);
                    for (int i = 0; i < Hbasis.rows; i++) {
                        kernelgens.add(new Generator(new IntTuple(v), Hbasis.getRow(i)));
                    }
                }
            }
            proj.add(new FreeFunctor(kernelgens, size));
            generators = new ArrayList<>();
            Collections.copy(kernelgens, generators);
        }
        return proj;
    }

    /**
     * Picks out the functor described by the first dims(v) coordinates at position v, for all v in N^r.
     * @param dims
     * @return
     * @throws WrongDimensionException
     */
    public Pair<Functor, Nat> getSubFunctor(Grid<Integer> dims) throws WrongDimensionException{
        if(!size.equals(dims.size())){
            throw new WrongDimensionException("Size of dimension grid is not equal to size of functor.");
        }
        Functor G = new Functor(size);
        Nat nat = new Nat(size);
        for(IntTuple v : GridIterator.getSequence(size)){
            Integer d = dims.get(v);
            try {
                BMatrix A = BMatrix.concat(BMatrix.identity(d), new BMatrix(d, getDimension(v) - d));
                nat.setMap(v, A.transpose());
            }catch (WrongDimensionException wde){
                log.error("Failed to concatenate matrix at position "+v+" when constructing natural transformation.", wde);
            }
            for(int i=0;i<v.length();i++){
                if(v.get(i).equals(size.get(i))){
                    G.setMap(v, BMatrix.identity(d), i);
                }else{
                    Integer dplus = dims.get(v.plus(IntTuple.getStandardBasisElement(v.length(), i)));
                    G.setMap(v, getMap(v, i).subMatrix(0, dplus, 0, d), i);
                }
            }
        }
        return new Pair<>(G, nat);
    }

    @Override
    public String toString(){
        StringBuilder sb = new StringBuilder();
        for(IntTuple v : GridIterator.getSequence(size)){
            sb.append(v).append(":\n");
            for(int i=0;i<v.length();i++){
                sb.append(i).append(":").append(getMap(v, i));
            }
        }
        return sb.toString();
    }


    //Static methods

    /**
     * Verifies that the functor is properly constructed. I.e that all maps are commutative etc.
     * @param F
     */
    public static void verify(Functor F){
        for(IntTuple v : GridIterator.getSequence(F.size)){
            List<IntTuple> basis = IntTuple.getStandardBasisSequence(v.length());
            for (int i = 0; i < basis.size(); i++) {
                for (int j = i + 1; j < basis.size(); j++) {
                    IntTuple w = v.plus(basis.get(i));
                    IntTuple wp = v.plus(basis.get(j));
                    IntTuple z = w.plus(basis.get(j));
                    if(!F.getMap(w, z).mult(F.getMap(v, w)).equals(F.getMap(wp, z).mult(F.getMap(v, wp)))){
                        System.out.println(w);
                    }
                    assert(F.getMap(w, z).mult(F.getMap(v, w)).equals(F.getMap(wp, z).mult(F.getMap(v, wp))));
                }
            }
        }
    }

    /**
     * Describes a generator of the functor.
     */
    public static class Generator{
        public IntTuple position;
        public BVector v;
        private int hash;

        public Generator(IntTuple position, BVector v) {
            this.position = position;
            this.v = v;
            this.hash = position.hashCode()+v.hashCode();
        }

        public static List<Generator> getGensLEQThan(List<Generator> generators, IntTuple u){
            List<Generator> gens = new ArrayList<>();
            for(int i=0;i<generators.size();i++){
                if(generators.get(i).position.leq(u)) gens.add(generators.get(i));
            }
            return gens;
        }

        @Override
        public String toString(){
            return (new StringBuilder())
                    .append("Position: ").append(position).append('\n')
                    .append(v)
                    .toString();
        }

        @Override
        public boolean equals(Object o){
            if(o == null || !o.getClass().equals(this.getClass())){
                return false;
            }
            Generator generator = (Generator) o;
            return this.position.equals(generator.position) && this.v.equals(generator.v);
        }

        @Override
        public int hashCode(){
            return hash;
        }
    }

}





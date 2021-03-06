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

package topcat.matrix;

import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.ints.IntIterator;
import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import topcat.matrix.exception.NoSolutionException;
import topcat.matrix.exception.WrongDimensionException;
import topcat.util.IntPair;
import topcat.util.Pair;

import java.util.*;

/**
 * Represents a matrix with coefficients in the field Z/2Z. It is implemented as a boolean matrix with element-wise
 * XOR as additive operation and the AND operation for matrix multiplication.
 */
public class BMatrix {
    private Int2ObjectOpenHashMap<BVector> A = new Int2ObjectOpenHashMap<>();
    public int rows, cols;

    public BMatrix(int m, int n){
        this.rows = m;
        this.cols = n;
    }

    public BMatrix(List<BVector> vectors){
        if(vectors == null || vectors.size() == 0){
            this.rows = 0;
            this.cols = 0;
        }else {
            this.rows = vectors.size();
            this.cols = vectors.get(0).getLength();
            for (int i = 0; i < this.rows; i++) {
                this.A.put(i, vectors.get(i));
            }
        }
    }

    public BMatrix(BMatrix B){
        this.rows = B.rows;
        this.cols = B.cols;
        for(int i=0;i<this.rows;i++){
            if(B.hasCreatedRow(i)) {
                this.A.put(i, new BVector(B.getRow(i)));
            }
        }
    }

    public boolean get(int i, int j){
        if(A.containsKey(i)){
            return A.get(i).get(j);
        }
        return false;
    }

    public BVector getRow(int i){
        if(A.containsKey(i)){
            return A.get(i);
        }
        return new BVector(cols);
    }

    public void setRow(int i, BVector row){
        A.put(i, row);
    }

    public void set(int i, int j, boolean val){
        if(!A.containsKey(i)){
            A.put(i, new BVector(cols));
        }
        this.A.get(i).set(j, val);
    }

    public boolean hasCreatedRow(int i){
        return A.containsKey(i);
    }

    public int[] getCreatedRowIndices(){
        return A.keySet().toIntArray();
    }

    public BMatrix mult(BMatrix B) {
        if(B == null) return null;
        if(cols != B.rows) {
            throw new WrongDimensionException("Matrix dimensions don't match. dim(A) = (" + rows + ", " + cols + "), dim(B) = (" + B.rows + ", " + B.cols + ")");
        }
        BMatrix C = new BMatrix(rows, B.cols);
        for(int i=0;i<rows;i++){
            if(A.containsKey(i)) {
                IntIterator iterator = A.get(i).getIndexSetIterator();
                BVector row = new BVector(B.cols);
                while (iterator.hasNext()) {
                    row = row.plus(B.getRow(iterator.nextInt()));
                }
                C.setRow(i, row);
            }
        }
        return C;
    }

    public BVector mult(BVector b){
        if(cols != b.getLength()) {
            throw new WrongDimensionException("Vector dimension don't match");
        }
        BVector x = new BVector(rows);
        for(int i=0;i<rows;i++){
            if(A.containsKey(i)) {
                IntIterator iterator = A.get(i).getIndexSetIterator();
                boolean val = false;
                while (iterator.hasNext()) {
                    val = val ^ b.get(iterator.nextInt());
                }
                x.set(i, val);
            }
        }
        return x;
    }

    public BMatrix plus(BMatrix B) throws WrongDimensionException{
        if(rows!=B.rows || cols!=B.cols) throw new WrongDimensionException("Matrix dimensions don't agree. dim(A)" +
                " = ("+rows+", "+cols+") dim(B) = ("+B.rows+", "+B.cols+")");
        BMatrix C = new BMatrix(rows, cols);
        for(int i=0;i<rows;i++){
            C.setRow(i, getRow(i).plus(B.getRow(i)));
        }
        return C;
    }

    public static BMatrix ker(BMatrix A){
        //Init basis elements
        BMatrix sourceBasis = identity(A.cols);
        BMatrix At = A.transpose();

        reduceRows(At, sourceBasis);

        int ctr = 0;
        for(int i=0;i<At.rows;i++){
            if(At.getRow(i).getNumberOfNonZeroElements() == 0) ctr++;
        }
        BMatrix kernelBasis = new BMatrix(ctr, At.rows);
        ctr=0;
        for(int i=0;i<At.rows;i++){
            if(At.getRow(i).getNumberOfNonZeroElements() == 0){
                kernelBasis.setRow(ctr, sourceBasis.getRow(i));
                ctr++;
            }
        }
        return kernelBasis;
    }

    public BMatrix transpose(){
        BMatrix At = new BMatrix(cols, rows);
        for(int i=0;i<rows;i++){
            if(A.containsKey(i)) {
                IntIterator iterator = A.get(i).getIndexSetIterator();
                while (iterator.hasNext()) {
                    At.set(iterator.nextInt(), i, true);
                }
            }
        }
        return At;
    }

    /**
     *
     * @param A - a matrix where the rows are a basis for a subspace
     * @return a matrix where the rows are the basis vectors
     */
    public static BMatrix extendBasis(BMatrix A) throws NoSolutionException{
        BMatrix M = vconcat(A, identity(A.cols));
        BMatrix Mreduced = new BMatrix(M);

        List<IntPair> pivots = reduceRows(Mreduced);

        BMatrix B = new BMatrix(A.cols, A.cols);
        for(int i=0;i<pivots.size();i++){
            B.setRow(i, M.getRow(pivots.get(i)._1()));
        }
        if(!A.equals(B.subMatrix(0, A.rows, 0, -1))){
            throw new NoSolutionException("The matrix A is not a basis! rank: "+rank(A)+" rows: "+A.rows);
        }
        return B;
    }

    /**
     * Returns a basis for the space spanned by the vectors described by the rows of A.
     * @param A
     * @return a basis (as row vectors) for the subspace spanned by the rows of A.
     */
    public static BMatrix getBasis(BMatrix A) throws WrongDimensionException{
        BMatrix At = A.transpose();
        final int[] column_sizes = new int[At.rows];
        for(int i=0;i<At.rows;i++) column_sizes[i] = At.getRow(i).getNumberOfNonZeroElements();
        List<IntPair> pivots = reduceRows(new BMatrix(A));
        BMatrix B = new BMatrix(pivots.size(), A.cols);
        for(int i=0;i<pivots.size();i++){
            B.setRow(i, A.getRow(pivots.get(i)._1()));
        }
        return B;
    }

    public static int rank(BMatrix A){
        if(A==null) return 0;
        int pivots=-1;
        try{
            pivots = reduceRows(new BMatrix(A), new BMatrix(A.rows, 1)).size();
        }catch (WrongDimensionException wde){
            //This cannot happen...
            wde.printStackTrace();
        }
        return pivots;
    }

    public static List<IntPair> reduceRows(BMatrix A) throws WrongDimensionException{
        List<IntPair> pivots = new ArrayList<>();
        int[] rowIndices = A.getCreatedRowIndices();
        for(int i=0; i<A.rows;i++){
            if(A.hasCreatedRow(i)) {
                BVector row = A.getRow(i);
                if (row.getNumberOfNonZeroElements() == 0){
                    continue;
                }
                int pivot = row.getIndexSetIterator().nextInt();
                pivots.add(new IntPair(i, pivot));
                for (int j : rowIndices) {
                    if (i == j){
                        continue;
                    }
                    if (A.get(j, pivot)) {
                        A.setRow(j, A.getRow(j).plus(A.getRow(i)));
                    }
                }
            }

        }
        return pivots;
    }

    public static List<IntPair> reduceRows(BMatrix A, BMatrix B) throws WrongDimensionException{
        List<IntPair> pivots = new ArrayList<>();
        int[] rowIndices = A.getCreatedRowIndices();
        for(int i=0; i<A.rows;i++){
            if(A.hasCreatedRow(i)) {
                BVector row = A.getRow(i);
                if (row.getNumberOfNonZeroElements() == 0){
                    continue;
                }
                int pivot = row.getIndexSetIterator().nextInt();
                pivots.add(new IntPair(i, pivot));
                for (int j : rowIndices) {
                    if (i == j){
                        continue;
                    }
                    if (A.get(j, pivot)) {
                        A.setRow(j, A.getRow(j).plus(A.getRow(i)));
                        B.setRow(j, B.getRow(j).plus(B.getRow(i)));
                    }
                }
            }

        }
        return pivots;
    }


    public static void reduceRows(BMatrix A, BMatrix B, List<IntPair> pivots) throws WrongDimensionException{
        int[] rowIndices = A.getCreatedRowIndices();
        for(IntPair pivot : pivots){
            int i = pivot._1();
            int _pivot = pivot._2();
            BVector row = A.getRow(i);
            if (row.getNumberOfNonZeroElements() == 0){
                continue;
            }
            for (int j : rowIndices) {
                if (i == j){
                    continue;
                }
                if (A.get(j, _pivot)) {
                    A.setRow(j, A.getRow(j).plus(A.getRow(i)));
                    B.setRow(j, B.getRow(j).plus(B.getRow(i)));
                }
            }

        }
    }


    public static Pair<BMatrix, BMatrix> reduction(final BMatrix A, List<Pair<Integer, Integer>> pivots){
        BMatrix At = A.transpose();
        BMatrix Atred = new BMatrix(At);
        BMatrix sourceBasis = identity(A.cols);
        //Reduce rows

        //reduceRows(Atred, sourceBasis, pivots);
        reduceRows(Atred, sourceBasis);

        List<BVector> imageVectors = new ArrayList<>();
        List<BVector> kernelVectors = new ArrayList<>();

        for(int i=0;i<At.rows;i++){
            if(Atred.hasCreatedRow(i) && Atred.getRow(i).getNumberOfNonZeroElements() > 0){
                imageVectors.add(At.getRow(i));
            }else{
                kernelVectors.add(sourceBasis.getRow(i));
            }
        }

        BMatrix imageBasis;
        if(imageVectors.size() > 0){
            imageBasis = new BMatrix(imageVectors);
        }else{
            imageBasis = new BMatrix(0, At.cols);
        }

        BMatrix kernelBasis;
        if(kernelVectors.size() > 0){
            kernelBasis = new BMatrix(kernelVectors);
        }else{
            kernelBasis = new BMatrix(0, A.cols);
        }
        return new Pair<>(kernelBasis, imageBasis);
    }


    /**
     * Performs a matrix reduction on A that reduces it to Smith-normal form. A basis
     * for the image and kernel of A is returned.
     * @param A
     * @return a pair of matrices (kernelbasis, imagebasis) where the rows of the
     * matrices are the basis vectors.
     */
    public static Pair<BMatrix, BMatrix> reduction(final BMatrix A){
        BMatrix At = A.transpose();
        BMatrix Atred = new BMatrix(At);
        BMatrix sourceBasis = identity(A.cols);
        //Reduce rows

        reduceRows(Atred, sourceBasis);//, column_ordering);

        List<BVector> imageVectors = new ArrayList<>();
        List<BVector> kernelVectors = new ArrayList<>();

        for(int i=0;i<At.rows;i++){
            if(Atred.hasCreatedRow(i) && Atred.getRow(i).getNumberOfNonZeroElements() > 0){
                imageVectors.add(At.getRow(i));
            }else{
                kernelVectors.add(sourceBasis.getRow(i));
            }
        }

        BMatrix imageBasis;
        if(imageVectors.size() > 0){
            imageBasis = new BMatrix(imageVectors);
        }else{
            imageBasis = new BMatrix(0, At.cols);
        }

        BMatrix kernelBasis;
        if(kernelVectors.size() > 0){
            kernelBasis = new BMatrix(kernelVectors);
        }else{
            kernelBasis = new BMatrix(0, A.cols);
        }
        return new Pair<>(kernelBasis, imageBasis);
    }

    /**
     * Returns a identity matrix of dimension 'dim'.
     * @param dim
     * @return
     */
    public static BMatrix identity(int dim){
        BMatrix A = new BMatrix(dim, dim);
        for(int i=0;i<dim;i++){
            A.set(i, i, true);
        }
        return A;
    }

    /**
     * Solves the linear equation AX=B where A, X and B are matrices.
     * @param A
     * @param B
     * @return
     * @throws NoSolutionException
     */
    public static boolean hasSolution(BMatrix A, BMatrix B) throws NoSolutionException{
        BMatrix Ap = new BMatrix(A);
        BMatrix Bp = new BMatrix(B);
        List<IntPair> pivots = reduceRows(Ap, Bp);

        IntOpenHashSet pivotHashSet = new IntOpenHashSet();
        for(IntPair pair : pivots){
            pivotHashSet.add(pair._1());
        }
        for(int i=0; i<Bp.rows; i++){
            if(!pivotHashSet.contains(i) && Bp.getRow(i).getNumberOfNonZeroElements() !=0){
                return false;
            }
        }
        return true;
    }

    public static boolean hasSolution(BMatrix A, BVector b) throws NoSolutionException{
        BMatrix B = new BMatrix(1, b.getLength());
        B.setRow(0, b);
        return hasSolution(A, B.transpose());
    }


    public static BMatrix solve(BMatrix A, BMatrix B){
        return solve(A, B, null);
    }

    /**
     * Solves the linear equation AX=B where A, X and B are matrices.
     * @param A
     * @param B
     * @return
     * @throws NoSolutionException
     */
    public static BMatrix solve(BMatrix A, BMatrix B, List<IntPair> pivots) throws NoSolutionException{
        BMatrix Ap = new BMatrix(A);
        BMatrix Bp = new BMatrix(B);

        pivots = null;

        if(pivots != null) {
            reduceRows(Ap, Bp, pivots);
        }else {
            pivots = reduceRows(Ap, Bp);
            //Sort the pivots by column
            Collections.sort(pivots, new Comparator<IntPair>() {
                @Override
                public int compare(IntPair o1, IntPair o2) {
                    return o1._2() - o2._2();
                }
            });
        }

        IntOpenHashSet pivotHashSet = new IntOpenHashSet();
        for(IntPair pair : pivots){
            pivotHashSet.add(pair._1());
        }
        for(int i=0; i<Bp.rows; i++){
            if(!pivotHashSet.contains(i) && Bp.getRow(i).getNumberOfNonZeroElements() !=0){
                throw new NoSolutionException("Found no solution to system of equations");
            }
        }

        BMatrix C = new BMatrix(Ap.cols, Bp.cols);
        for(IntPair pair : pivots){
            C.setRow(pair._2(), new BVector(Bp.getRow(pair._1())));
        }
        return C;
    }

    /**
     * Solves the linear equation Ax=b where A is a matrix and x and b are vectors.
     * @param A
     * @param b
     * @return
     * @throws NoSolutionException
     */
    public static BVector solve(BMatrix A, BVector b) throws NoSolutionException{
        BMatrix B = new BMatrix(1, b.getLength());
        B.setRow(0, b);
        BMatrix C = solve(A, B.transpose());
        return C.transpose().getRow(0);
    }

    public BMatrix inverse() throws NoSolutionException{
        return inverse(null);
    }

    public BMatrix inverse(List<IntPair> pivots) throws NoSolutionException{
        if(rows != cols){
            throw new WrongDimensionException("Matrix is not square! rows = "+rows+" cols = "+cols);
        }
        return solve(this, BMatrix.identity(rows), pivots);
    }

    /**
     * Returns the submatrix A(start_row:end_row-1, start_col:end_col-1) (in MATLAB notation).
     * @param start_row
     * @param end_row
     * @param start_col
     * @param end_col
     * @return
     */
    public BMatrix subMatrix(int start_row, int end_row, int start_col, int end_col){
        if(end_row==-1) end_row=rows;
        if(end_col==-1) end_col=cols;
        BMatrix A = new BMatrix(end_row-start_row, end_col-start_col);
        for(int i=start_row;i<end_row;i++){
            IntIterator iterator = getRow(i).getIndexSetIterator();
            while(iterator.hasNext()){
                int index = iterator.nextInt();
                if(index >= start_col && index < end_col) A.set(i-start_row, index-start_col, true);
            }
        }
        return A;
    }

    /**
     * Concatenates the matrix A and B horizontally into [A B].
     * @param A
     * @param B
     * @return [A B]
     * @throws WrongDimensionException
     */
    public static BMatrix concat(BMatrix A, BMatrix B) throws WrongDimensionException{
        if(A==null) return B;
        if(B == null) return A;
        if(A.rows != B.rows) {
            throw new WrongDimensionException("Dimensions do not agree: Rows must be equal A.rows = " + A.rows + " B.rows = " + B.rows);
        }
        BMatrix C = new BMatrix(A.rows, A.cols+B.cols);
        for(int i=0;i<C.rows;i++){
            C.setRow(i, BVector.concat(A.getRow(i), B.getRow(i)));
        }
        return C;
    }

    /**
     * Concatenates the matrix A and B vertically into [A; B].
     * @param A
     * @param B
     * @return [A; B]
     * @throws WrongDimensionException
     */
    public static BMatrix vconcat(BMatrix A, BMatrix B) throws WrongDimensionException{
        if(A==null) return B;
        if(B == null) return A;
        if(A.cols != B.cols) {
            throw new WrongDimensionException("Dimensions do not agree: cols must be equal A.cols = " + A.cols + " B.cols = " + B.cols);
        }
        BMatrix C = new BMatrix(A.rows+B.rows, A.cols);
        for(int i=0;i<A.rows;i++){
            C.setRow(i, A.getRow(i));
        }
        for(int i=0;i<B.rows;i++){
            C.setRow(A.rows+i, B.getRow(i));
        }
        return C;
    }

    @Override
    public boolean equals(Object o){
        if(o == null || !o.getClass().equals(this.getClass())) return false;
        BMatrix B = (BMatrix) o;
        if(rows != B.rows || cols != B.cols) return false;
        for(int i=0;i<rows;i++){
            if(!getRow(i).equals(B.getRow(i))) return false;
        }
        return true;
    }

    @Override
    public String toString(){
        if(rows == 0 && cols > 0){
            return "[]";
        }
        StringBuilder sb = new StringBuilder();
        for(int i=0;i<rows;i++){
            sb.append("[ ");
            for(int j=0;j<cols;j++){
                if(get(i, j)){
                    sb.append(1).append(' ');
                }else{
                    sb.append(0).append(' ');
                }
            }
            sb.append(']').append('\n');
        }
        return sb.toString();
    }
}

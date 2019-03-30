package topcat.matrix;

import topcat.persistence.simplex.Simplex;

public class GradedColumn<T> extends Column<T> implements Comparable<GradedColumn<T>>{
    private Simplex grade;

    public GradedColumn(Simplex grade){
        super();
        this.grade = grade;
    }

    public Simplex getGrade(){
        return grade;
    }

    @Override
    public int compareTo(GradedColumn<T> o) {
        return this.grade.compareTo(o.getGrade());
    }
}

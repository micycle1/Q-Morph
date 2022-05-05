package meshditor;

/**
 * This class holds information for lines, and has methods for dealing with
 * line-related issues. The purpose of this class is solely to determine the
 * intersection point between two lines. The length of a line is, of course,
 * infinite.
 */

public class MyLine {
	// n1: A point that the line passes through
	// n2: A point that the line passes through
	public MyLine(Vertex n1, Vertex n2) {
		ref = n1;
		x = n2.x - n1.x;
		y = n2.y - n1.y;
	}

	public MyLine(Vertex n1, double x, double y) {
		ref = n1;
		this.x = x;
		this.y = y;
	}

	public double cross(MyLine l) {
		return x * l.y - l.x * y;
	}

	public Vertex pointIntersectsAt(MyLine d1) {
		Vertex p0 = ref, p1 = d1.ref;
		MyLine delta = new MyLine(p0, p1.x - p0.x, p1.y - p0.y);
		MyLine d0 = this;
		double d0crossd1 = d0.cross(d1);

		if (d0crossd1 == 0) { // Parallel and, alas, no pointintersection
			return null;
		} else {
			double t = delta.cross(d0) / d0crossd1;

			double x = d1.ref.x + t * d1.x;
			double y = d1.ref.y + t * d1.y;
			return new Vertex(x, y); // Intersects at this line point
		}
	}

	public String descr() {
		return ref.descr() + ", (" + (x + ref.x) + ", " + (y + ref.y) + ")";
	}

	Vertex ref;
	double x, y;
}

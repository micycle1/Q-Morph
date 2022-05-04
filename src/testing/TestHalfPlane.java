package testing;

import meshditor.Edge;
import meshditor.Msg;
import meshditor.MyVector;
import meshditor.Vertex;
import meshditor.Quad;
import meshditor.Triangle;

class TestHalfPlane {

	public static void main(String[] args) {
		Vertex pa = new Vertex(0.0, 0.0);
		Vertex pb = new Vertex(0.0, 2.0);
		Vertex pc = new Vertex(-1.0, 0.0);

		Edge e1 = new Edge(pa, pb);
		Edge e2 = new Edge(pb, pc);
		Edge e3 = new Edge(pa, pc);

		Triangle t = new Triangle(e1, e2, e3);

		Msg.debug("*** First set of test Vertexes (not in halfplane) ***");
		Vertex p1 = new Vertex(7.2, 0.0);
		Vertex p2 = new Vertex(0.2, -30.0);
		Vertex p3 = new Vertex(48, 53.2);
		Vertex p4 = new Vertex(0.01, 0.0);

		if (p1.inHalfplane(t, e1) == 1) {
			Msg.error("p1 incorrectly detected to belong to halfplane left of Edge " + e1.descr());
		} else {
			Msg.debug("p1 correctly detected not to belong to halfplane left of Edge " + e1.descr());
		}

		if (p2.inHalfplane(t, e1) == 1) {
			Msg.error("p2 incorrectly detected to belong to halfplane left of Edge " + e1.descr());
		} else {
			Msg.debug("p2 correctly detected not to belong to halfplane left of Edge " + e1.descr());
		}

		if (p3.inHalfplane(t, e1) == 1) {
			Msg.error("p3 incorrectly detected to belong to halfplane left of Edge " + e1.descr());
		} else {
			Msg.debug("p3 correctly detected not to belong to halfplane left of Edge " + e1.descr());
		}

		if (p4.inHalfplane(t, e1) == 1) {
			Msg.error("p4 incorrectly detected to belong to halfplane left of Edge " + e1.descr());
		} else {
			Msg.debug("p4 correctly detected not to belong to halfplane left of Edge " + e1.descr());
		}

		Msg.debug("*** Second set of test Vertexes (in halfplane) ***");
		p1 = new Vertex(-0.1, 0.0);
		p2 = new Vertex(-21.3, -123.2);
		p3 = new Vertex(-0.01, 100.0);
		p4 = new Vertex(-2.0, 12.3);

		if (p1.inHalfplane(t, e1) == 1) {
			Msg.debug("p1 correctly detected to belong to halfplane left of Edge " + e1.descr());
		} else {
			Msg.error("p1 incorrectly detected not to belong to halfplane left of Edge " + e1.descr());
		}

		if (p2.inHalfplane(t, e1) == 1) {
			Msg.debug("p2 correctly detected to belong to halfplane left of Edge " + e1.descr());
		} else {
			Msg.error("p2 incorrectly detected not to belong to halfplane left of Edge " + e1.descr());
		}

		if (p3.inHalfplane(t, e1) == 1) {
			Msg.debug("p3 correctly detected to belong to halfplane left of Edge " + e1.descr());
		} else {
			Msg.error("p3 incorrectly detected not to belong to halfplane left of Edge " + e1.descr());
		}

		if (p4.inHalfplane(t, e1) == 1) {
			Msg.debug("p4 correctly detected to belong to halfplane left of Edge " + e1.descr());
		} else {
			Msg.error("p4 incorrectly detected not to belong to halfplane left of Edge " + e1.descr());
		}

		Msg.debug("*** A 'real life' test (in halfplane) ***");
		Vertex a = new Vertex(-1.3, 4.8);
		Vertex b = new Vertex(-8.8, 0.7);
		Vertex c = new Vertex(-2.3, -2.7);

		Vertex p = new Vertex(-2.4, 0.7);

		Edge ab = new Edge(a, b);
		Edge bc = new Edge(b, c);
		Edge ac = new Edge(a, c);

		Triangle abc = new Triangle(ab, bc, ac);

		if (p.inHalfplane(abc, bc) == 1) {
			Msg.debug("p correctly detected to belong to halfplane right of Edge " + bc.descr());
		} else {
			Msg.error("p incorrectly detected not to belong to halfplane right of Edge " + bc.descr());
		}

		testConvexity();
	}

	public static void testConvexity() {
		Vertex n1 = new Vertex(0.0, 0.0);
		Vertex n2 = new Vertex(0.0, 1.4);
		Vertex n3 = new Vertex(-0.4, 2.9);
		Vertex n4 = new Vertex(3.8, 0.9);

		Edge baseEdge = new Edge(n1, n2);
		Edge leftEdge = new Edge(n2, n3);
		Edge rightEdge = new Edge(n1, n4);
		Edge topEdge = new Edge(n3, n4);

		Quad q = new Quad(baseEdge, leftEdge, rightEdge, topEdge);

		if (q.isStrictlyConvex()) {
			Msg.error("method says Quad " + q.descr() + " is strictly convex, but it isn't");
		} else {
			Msg.debug("method correctly reports that Quad " + q.descr() + " is not strictly convex");
		}

		// Run a test on the intersects method of MyVector:
		Vertex n5 = new Vertex(0.0, 0.2);
		Vertex n6 = new Vertex(0.0, 1.1);

		MyVector v0 = new MyVector(n1, n2);
		MyVector v1 = new MyVector(n5, n6);

		if (v0.intersects(v1)) {
			Msg.debug("Correctly identified intersection");
		} else {
			Msg.error("Incorrectly rejected intersection");
		}
	}

}

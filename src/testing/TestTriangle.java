package testing;

import meshditor.Edge;
import meshditor.GeomBasics;
import meshditor.Msg;
import meshditor.Triangle;
import meshditor.Vertex;

public class TestTriangle extends GeomBasics {

	public static void main(String[] args) {

		Msg.debugMode = true;

		Vertex n1 = new Vertex(0, 0);
		Vertex n2 = new Vertex(3, 0);
		Vertex n3 = new Vertex(0, 3);

		Edge e1 = new Edge(n1, n2);
		Edge e2 = new Edge(n1, n3);
		Edge e3 = new Edge(n2, n3);
		Edge cwEdge, ccwEdge;

		Triangle t = new Triangle(e1, e2, e3);
		Msg.debug("Triangle is " + t.descr());

		Msg.debug("- --- - Test of nextCWEdge and nextCCWEdge - --- -");
		cwEdge = t.nextCWEdge(e1);
		ccwEdge = t.nextCCWEdge(e1);
		Msg.debug("cw to " + e1.descr() + " is " + cwEdge.descr());
		Msg.debug("ccw to " + e1.descr() + " is " + ccwEdge.descr());
	}
}

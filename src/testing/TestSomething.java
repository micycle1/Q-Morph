package testing;

import meshditor.Constants;
import meshditor.Msg;
import meshditor.MyVector;
import meshditor.Vertex;

public class TestSomething extends Constants {

	public static void main(String[] args) {
		double dx = 3.9;
		int ix = (int) dx;

		Msg.debugMode = true;

		Msg.debug("double value is: " + dx);
		Msg.debug("Value of double casted to int is: " + ix);

		Msg.debug("Value of remainder of " + dx + "/3.0==" + Math.IEEEremainder(dx, 3.0));

		Msg.debug("Testing the old and new MyVector.dot(..) methods:");

		Vertex n0 = new Vertex(0, 0);
		Vertex n1 = new Vertex(1, 0);
		Vertex n2 = new Vertex(1, 1);
		Vertex n3 = new Vertex(0, 1);
		Vertex n4 = new Vertex(-1, 1);
		Vertex n5 = new Vertex(-1, 0);

		MyVector v1 = new MyVector(n0, n1);
		MyVector v2 = new MyVector(n0, n2);
		MyVector v3 = new MyVector(n0, n3);
		MyVector v5 = new MyVector(n0, n5);

		double v2x = 1.0, v2y = 1.0;
		double v3x = 0.0, v3y = 1.0;
		double v5x = -1.0, v5y = 0.0;

		Msg.debug("v1== " + v1.descr());
		Msg.debug("v2== " + v2.descr());
		Msg.debug("v3== " + v3.descr());
		Msg.debug("v5== " + v5.descr());

		Msg.debug("v1.dot(v3)== " + v1.dot(v3));
		Msg.debug("v1.dot(v3x, v3y)== " + v1.dot(v3x, v3y));

		Msg.debug("v1.dot(v2)== " + v1.dot(v2));
		Msg.debug("v1.dot(v2x, v2y)== " + v1.dot(v2x, v2y));

		Msg.debug("v1.dot(v5)== " + v1.dot(v5));
		Msg.debug("v1.dot(v5x, v5y)== " + v1.dot(v5x, v5y));
	}
}

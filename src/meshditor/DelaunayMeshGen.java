package meshditor;

import java.util.ArrayList;

/**
 * This class offers methods for incrementally constructing Delaunay triangle
 * meshes.
 */

public class DelaunayMeshGen extends GeomBasics {
	public DelaunayMeshGen() {
	}

	/**
	 * Initialize the class
	 * 
	 * @param delaunayCompliant boolean to indicate whether to create a delaunay
	 *                          compliant mesh or not.
	 */
	boolean delaunayCompliant = false;

	public void init(boolean delaunayCompliant) {
		this.delaunayCompliant = delaunayCompliant;
		setCurMethod(this);

		// Perform the steps necessary before inserting the Vertices in incrDelaunay():
		// Create the two initial Delaunay triangles from the four most extreme Vertices.

		triangleList = new ArrayList();
		edgeList = new ArrayList();
		findExtremeVertices();

		

		// The boundary edges
		Edge edge2 = new Edge(leftmost, uppermost);
		Edge edge3 = new Edge(leftmost, lowermost);
		Edge edge4 = new Edge(lowermost, rightmost);
		Edge edge5 = new Edge(uppermost, rightmost);

		// * Create the two initial triangles
		// * They are made Delaunay by selecting the correct interior edge.
		Triangle del1, del2;
		Edge edge1;
		if (!rightmost.inCircle(uppermost, leftmost, lowermost)) {
			edge1 = new Edge(uppermost, lowermost);
			del1 = new Triangle(edge1, edge2, edge3);
			del2 = new Triangle(edge1, edge4, edge5);
		} else {
			edge1 = new Edge(leftmost, rightmost);
			del1 = new Triangle(edge1, edge2, edge5);
			del2 = new Triangle(edge1, edge3, edge4);
		}

		edge1.connectVertices();
		edge2.connectVertices();
		edge3.connectVertices();
		edge4.connectVertices();
		edge5.connectVertices();

		del1.connectEdges();
		del2.connectEdges();

		// Update "global" edgeList, triangleList, and vertexList
		edgeList.add(edge1);
		edgeList.add(edge2);
		edgeList.add(edge3);
		edgeList.add(edge4);
		edgeList.add(edge5);
		triangleList.add(del1);
		triangleList.add(del2);
	}

	// Run the implementation on the give set of vertices. */
	public void run() {
		
		// Point insertions
		Vertex n;
		for (Object element : vertexList) {
			n = (Vertex) element;

			// The extreme vertices have been inserted already, so we skip them here
			if (n != leftmost && n != rightmost && n != uppermost && n != lowermost) {
				insertVertex(n, delaunayCompliant);
			}
		}
	}

	private int counter = 0;

	/**
	 * Step through the insertions one by one. Upon completion, run
	 * incrDelaunayPost().
	 */
	@Override
	public void step() {
		Vertex n;
		if (counter < vertexList.size()) {
			n = (Vertex) vertexList.get(counter);
			counter++;
			// The extreme vertices have been inserted already, so we skip them here
			while (n == leftmost || n == rightmost || n == uppermost || n == lowermost) {
				if (counter < vertexList.size()) {
					n = (Vertex) vertexList.get(counter);
					counter++;
				} else if (counter == vertexList.size()) {
					counter++;
					return;
				} else {
					return;
				}
			}
			if (n != leftmost && n != rightmost && n != uppermost && n != lowermost) {
				insertVertex(n, delaunayCompliant);
			}
		} else if (counter == vertexList.size()) {
			counter++;
		}
	}

	/**
	 * Find the triangle (in the list) that contains the specified Vertex. A method
	 * highly inspired by the dart based localization procedure.
	 * 
	 * @return the Triangle containing the Vertex, but if the Vertex is located on an
	 *         Edge, this Edge is returned instead. Set inside= true if found. If
	 *         the Vertex is located outside of the current triangulation, then set
	 *         inside= false and return a boundary triangle that can be seen from
	 *         the Vertex.
	 */
	private Object findTriangleContaining(Vertex newVertex, Triangle start) {
		
		// The initial "dart" must be ccw in the initial triangle:
		Triangle ts = start; // d_start
		Edge es; // d_start
		Edge e1 = ts.edgeList[0];
		Edge e2 = ts.edgeList[1];
		Vertex ns = e1.commonVertex(e2); // d_start
		MyVector v1 = new MyVector(ns, e1.otherVertex(ns));
		MyVector v2 = new MyVector(ns, e2.otherVertex(ns));
		Edge online = null;

		if (v1.isCWto(v2)) {
			es = e2;
		} else {
			es = e1;
		}

		Triangle t = ts; // d_i
		Edge e = es; // d_i
		Vertex n = ns; // d_i

		int count = 0;
		int hp;
		while (true) {

			
			/*
			 * Msg.debug("Loop nr. "+count++); Msg.debug("newVertex= "+newVertex.descr());
			 * Msg.debug("Current dart has triangle= "+t.descr());
			 * Msg.debug("Current dart has edge= "+e.descr());
			 * Msg.debug("Current dart has Vertex= "+n.descr());
			 * Msg.debug("Stop dart has triangle= "+ts.descr());
			 * Msg.debug("Stop dart has edge= "+es.descr());
			 * Msg.debug("Stop dart has Vertex= "+ns.descr());
			 */

			hp = newVertex.inHalfplane(t, e);
			if (hp == 1 || hp == 0) { // || (hp== 0 && !newVertex.inBoundedPlane(e))) {
				if (hp == 0) {
					online = e;
				}
				// is newVertex in halfplane defined by (t, e)?
				n = e.otherVertex(n);
				e = t.neighborEdge(n, e);
				if (ts == t && es == e && ns == n) {
					inside = true;

					if (online != null) {
						return online;
					} else {
						return t;
					}
				}
			}
			/*
			 * else if (hp==0 && newVertex.inBoundedPlane(e)) { // is newVertex actually on Edge
			 * e? Msg.debug("Leaving findTriangleContaining(..), returning Edge"); inside=
			 * true; return (Object) e; }
			 */
			else { // try to move to the adjacent triangle
				online = null;
				ts = (Triangle) t.neighbor(e);
				if (ts == null) {
					/*
					 * if (hp== 0) { e= ; t= ; es= ; ts= ; } else {
					 */
					inside = false;
					return e; // outside triangulation
					// }
				} else {
					t = ts; // d_start= alpha_0 o alpha_2(d_i)
					es = e;
					ns = e.otherVertex(n);

					// d_i= alpha_1 o alpha_2(d_i)
					e = t.neighborEdge(n, e);
				}
			}
		}
	}

	/** Simple method to perform single swap. */
	private void swap(Edge e) {
		Triangle t1 = (Triangle) e.element1;
		Triangle t2 = (Triangle) e.element2;

		if (t1 == null | t2 == null) {
			return;
		}

		Vertex na, nb, nc, nd;
		na = e.leftVertex;
		nb = e.rightVertex;
		nc = t1.oppositeOfEdge(e);
		nd = t2.oppositeOfEdge(e);

		double cross1 = cross(nc, na, nd, na); // The cross product ca x da
		double cross2 = cross(nc, nb, nd, nb); // The cross product cb x db

		if (cross1 == 0 || cross2 == 0) {
			// if (!q.isStrictlyConvex()) {
			return;
		}

		// Create the new Edge, do the swap
		Edge ei = new Edge(nc, nd);
		e.swapToAndSetElementsFor(ei);
		Triangle tNew1 = (Triangle) ei.element1;
		Triangle tNew2 = (Triangle) ei.element2;

		// Update "global" lists: remove old triangles and edge e, add new ones
		int ind1 = triangleList.indexOf(t1);
		if (ind1 != -1) {
			triangleList.remove(ind1);
		}
		int ind2 = triangleList.indexOf(t2);
		if (ind2 != -1) {
			triangleList.remove(ind2);
		}

		edgeList.remove(edgeList.indexOf(e));
		edgeList.add(ei);
		triangleList.add(tNew1);
		triangleList.add(tNew2);
	}

	/** Recursive method that swaps Edges in order to maintain Delaunay property. */
	private void recSwapDelaunay(Edge e, Vertex n) {
		Triangle t1 = (Triangle) e.element1;
		Triangle t2 = (Triangle) e.element2;
		Vertex na, nb, nc, nd;

		if (t1 == null || t2 == null) {// Make sure we're dealing with an interior edge
			return;
		}

		Triangle t;
		if (!e.element1.hasVertex(n)) {
			t = (Triangle) e.element1;
		} else {
			t = (Triangle) e.element2;
		}

		Vertex p1, p2, p3, opposite = t.oppositeOfEdge(e);
		Quad q = new Quad(e, opposite, n);

		nc = n;
		nd = e.oppositeVertex(n);
		na = e.leftVertex;
		nb = e.rightVertex;

		double cross1 = cross(nc, na, nd, na); // The cross product ca x da
		double cross2 = cross(nc, nb, nd, nb); // The cross product cb x db

		if (cross1 == 0 || cross2 == 0) {
			// if (!q.isStrictlyConvex()) {
			return;
		}

		p1 = q.nextCCWVertex(n);
		p2 = q.nextCCWVertex(p1);
		p3 = q.nextCCWVertex(p2);

		if (!n.inCircle(p1, p2, p3)) { // If n lies outside the cicrumcircle..
			return;
		}

		// Create the new Edge, do the swap
		Edge ei = new Edge(p2, n);
		e.swapToAndSetElementsFor(ei);
		Triangle tNew1 = (Triangle) ei.element1;
		Triangle tNew2 = (Triangle) ei.element2;

		// Update "global" lists: remove old triangles and edge e, add new ones
		int ind1 = triangleList.indexOf(t1);
		if (ind1 != -1) {
			triangleList.remove(ind1);
		}
		int ind2 = triangleList.indexOf(t2);
		if (ind2 != -1) {
			triangleList.remove(ind2);
		}

		edgeList.remove(edgeList.indexOf(e));
		edgeList.add(ei);
		triangleList.add(tNew1);
		triangleList.add(tNew2);

		// Proceed with recursive calls
		recSwapDelaunay(tNew1.oppositeOfVertex(n), n);
		recSwapDelaunay(tNew2.oppositeOfVertex(n), n);
	}

	/**
	 * Recursive method that deletes the part of the mesh boundary that is no longer
	 * Delaunay compliant (the "influence region") when a new Vertex has been inserted
	 * exterior to the current mesh. Also compiles a list of Vertices in this region.
	 * The method must be called for each boundary triangle that is affected by the
	 * newly inserted Vertex. Some interior triangles are also affected, but you don't
	 * need to worry about them, because they get deleted automagically by this
	 * method.
	 *
	 * @param t a triangle on the boundary that is possibly no longer Delaunay
	 *          compliant
	 * @param e the boundary edge of this triangle
	 * @param n the exterior Vertex that has recently been inserted
	 */
	private void makeDelaunayTriangle(Triangle t, Edge e, Vertex n) {
		int j;
		Edge e1, e2;
		Triangle t1, t2;
		Vertex p1, p2, p3, opposite = t.oppositeOfEdge(e);
		Quad q = new Quad(e, opposite, n);

		if (!q.isStrictlyConvex()) {
			j = irVertices.indexOf(e.leftVertex);
			if (j == -1) {
				irVertices.add(e.leftVertex);
			}
			j = irVertices.indexOf(e.rightVertex);
			if (j == -1) {
				irVertices.add(e.rightVertex);
			}

			return;
		}

		p1 = q.nextCCWVertex(n);
		p2 = q.nextCCWVertex(p1);
		p3 = q.nextCCWVertex(p2);

		if (n.inCircle(p1, p2, p3)) {
			j = irVertices.indexOf(p1);
			if (j == -1) {
				irVertices.add(p1);
			}
			j = irVertices.indexOf(p2);
			if (j == -1) {
				irVertices.add(p2);
			}
			j = irVertices.indexOf(p3);
			if (j == -1) {
				irVertices.add(p3);
			}

			e1 = t.otherEdge(e);
			t1 = (Triangle) t.neighbor(e1);
			e2 = t.otherEdge(e, e1);
			t2 = (Triangle) t.neighbor(e2);

			e.disconnectVertices();
			j = edgeList.indexOf(e);
			if (j != -1) {
				edgeList.remove(j);
			}

			if (t1 == null) {
				e1.disconnectVertices();
				j = edgeList.indexOf(e1);
				if (j != -1) {
					edgeList.remove(j);
				}
			}
			if (t2 == null) {
				e2.disconnectVertices();
				j = edgeList.indexOf(e2);
				if (j != -1) {
					edgeList.remove(j);
				}
			}

			t.disconnectEdges();
			triangleList.remove(triangleList.indexOf(t));

			if (t1 != null) {
				makeDelaunayTriangle(t1, e1, n);
			}
			if (t2 != null) {
				makeDelaunayTriangle(t2, e2, n);
			}
		} else {
			j = irVertices.indexOf(e.leftVertex);
			if (j == -1) {
				irVertices.add(e.leftVertex);
			}
			j = irVertices.indexOf(e.rightVertex);
			if (j == -1) {
				irVertices.add(e.rightVertex);
			}
		}
	}

	private boolean inside = false;
	private ArrayList irVertices = new ArrayList();

	/**
	 * Insert a interior/ exterior Vertex and update mesh to remain Delaunay
	 * compliant.
	 */
	private void insertVertex(Vertex n, boolean remainDelaunay) {
		Triangle t, t1, t2 = null, t3, t4 = null, oldt1, oldt2;
		Object o;
		Edge e, e1 = null, e2, e3, e4 = null, e12, e22 = null, e32, e42 = null, eOld, b0 = null, b1 = null;
		Vertex vertex, nVertex, pVertex, other, other2, n0, n1;
		ArrayList boundaryEdges = new ArrayList();
		int i, j;
		boolean loop = true;

		// Locate the triangle that contains the point
		o = findTriangleContaining(n, (Triangle) triangleList.get(0));
		if (!inside) {
			// --- the Vertex is to be inserted outside the current triangulation --- //
			e = (Edge) o;
			

			// Compile an ordered list of boundary edges that is part of the i. polygon.
			// First find the leftmost edge of these boundary edges...
			pVertex = e.rightVertex;
			vertex = e.leftVertex;

			while (loop) {
				e1 = vertex.anotherBoundaryEdge(e);
				nVertex = e1.otherVertex(vertex);
				j = nVertex.inHalfplane(n, vertex, pVertex);
				if (j == -1) { // pVertex and nVertex must lie on different sides of (n, Vertex)
					pVertex = vertex;
					vertex = nVertex;
					e = e1;
				} else {
					loop = false;
				}
			}
			

			// ... then traverse the boundary edges towards the right, adding one
			// edge at a time until the rightmost edge is encountered
			loop = true;
			pVertex = vertex;
			n0 = vertex; // Most distant Vertex to the left
			vertex = e.otherVertex(vertex);
			boundaryEdges.add(e);

			while (loop) {
				e1 = vertex.anotherBoundaryEdge(e);

				nVertex = e1.otherVertex(vertex);
				j = nVertex.inHalfplane(n, vertex, pVertex);
				if (j == -1) { // pVertex and nVertex must lie on different sides of (n, Vertex)
					pVertex = vertex;
					vertex = nVertex;
					boundaryEdges.add(e1);
					e = e1;
				} else {
					loop = false;
				}
			}
			n1 = vertex; // Most distant Vertex to the right
			

			// From this list, find each triangle in the influence region and delete it.
			// Also create the list of vertices in the influence region
			// makeDelaunayTriangle does the job.
			for (i = 0; i < boundaryEdges.size(); i++) {
				e = (Edge) boundaryEdges.get(i);
				t = e.getTriangleElement();
				if (t != null) {
					makeDelaunayTriangle(t, e, n);
				}
			}
			

			// Create new triangles by connecting the vertices in the
			// influence region to Vertex n:

			// Build initial edgeList for Vertex n
			for (i = 0; i < irVertices.size(); i++) {
				other = (Vertex) irVertices.get(i);
				e = new Edge(n, other);
				e.connectVertices();
				edgeList.add(e);
				if (other == n0) {
					b0 = e;
				} else if (other == n1) {
					b1 = e;
				}
			}
			

			// Sort this list of edges in ccw order
			n.edgeList = n.calcCCWSortedEdgeList(b0, b1);
			

			// Create each triangle
			printEdgeList(n.edgeList);
			for (i = 0; i < n.edgeList.size() - 1; i++) {
				e1 = (Edge) n.edgeList.get(i);
				e2 = (Edge) n.edgeList.get(i + 1);
				other = e1.otherVertex(n);
				other2 = e2.otherVertex(n);
				e = new Edge(other, other2);
				j = other.edgeList.indexOf(e);
				if (j != -1) {
					e = (Edge) other.edgeList.get(j);
				} else {
					e.connectVertices();
					edgeList.add(e);
				}

				t = new Triangle(e, e1, e2);
				triangleList.add(t);
				t.connectEdges();
			}
			irVertices.clear(); // NB! IMPORTANT!
		} else if (o instanceof Triangle) {
			// --- the Vertex is to be inserted inside the current triangulation --- //
			t = (Triangle) o;
			

			// Make some pointers for its edges
			Edge te1 = t.edgeList[1], te2 = t.edgeList[2], te3 = t.edgeList[0];

			// Split the triangle into three new triangles with n as a common Vertex.
			e1 = new Edge(n, te3.commonVertex(te1));
			e2 = new Edge(n, te1.commonVertex(te2));
			e3 = new Edge(n, te2.commonVertex(te3));

			e1.connectVertices();
			e2.connectVertices();
			e3.connectVertices();

			edgeList.add(e1);
			edgeList.add(e2);
			edgeList.add(e3);
			t1 = new Triangle(e1, e2, te1); // This should be correct...
			t2 = new Triangle(e2, e3, te2); //
			t3 = new Triangle(e3, e1, te3); //

			// Disconnect Edges from old Triangle & connect Edges to new triangles.
			t.disconnectEdges();

			t1.connectEdges();
			t2.connectEdges();
			t3.connectEdges();

			// Update triangleList.
			triangleList.remove(triangleList.indexOf(t));
			triangleList.add(t1);
			triangleList.add(t2);
			triangleList.add(t3);

			// Swap edges so that the new triangulation becomes Delauney
			if (remainDelaunay) {
				recSwapDelaunay(te1, n);
				recSwapDelaunay(te2, n);
				recSwapDelaunay(te3, n);
			} else {
				if (!t1.areaLargerThan0()) {
					swap(te1);
				}
				if (!t2.areaLargerThan0()) {
					swap(te2);
				}
				if (!t3.areaLargerThan0()) {
					swap(te3);
				}
			}
		} else if (o instanceof Edge) { // n lies on Edge e:
			// Split the (1 or) 2 Triangles adjacent Edge e into (2 or) 4 new Triangles.
			e = (Edge) o;

			oldt1 = (Triangle) e.element1;
			oldt2 = (Triangle) e.element2;

			// Create the (2 or) 4 new Edges, get ptrs for the (2 or) 4 outer Edges,
			// remove the old Edge e.
			e1 = new Edge(e.leftVertex, n);
			e2 = new Edge(e.rightVertex, n);
			e12 = oldt1.neighborEdge(e.leftVertex, e);
			e3 = new Edge(n, e12.otherVertex(e.leftVertex));
			e32 = oldt1.neighborEdge(e.rightVertex, e);

			e1.connectVertices();
			e2.connectVertices();
			e3.connectVertices();
			edgeList.add(e1);

			edgeList.add(e2);
			edgeList.add(e3);
			if (oldt2 != null) {
				e22 = oldt2.neighborEdge(e.leftVertex, e);
				e4 = new Edge(n, e22.otherVertex(e.leftVertex));
				e4.connectVertices();
				e42 = oldt2.neighborEdge(e.rightVertex, e);

				edgeList.add(e4);
			}

			e.disconnectVertices();
			edgeList.remove(edgeList.indexOf(e));
			

			// Create the (2 or) 4 new triangles
			t1 = new Triangle(e1, e12, e3); // This should be correct...
			t3 = new Triangle(e2, e3, e32);
			if (oldt2 != null) {
				t2 = new Triangle(e1, e22, e4);
				t4 = new Triangle(e2, e4, e42);
			}

			// Disconnect Edges from old Triangles & connect Edges to new triangles.
			oldt1.disconnectEdges();
			if (oldt2 != null) {
				oldt2.disconnectEdges();
			}

			t1.connectEdges();

			t3.connectEdges();
			if (oldt2 != null) {
				t2.connectEdges();
				t4.connectEdges();
			}

			// Update triangleList.
			triangleList.remove(triangleList.indexOf(oldt1));
			triangleList.add(t1);
			triangleList.add(t3);
			if (oldt2 != null) {
				triangleList.remove(triangleList.indexOf(oldt2));
				triangleList.add(t2);
				triangleList.add(t4);
			}

			// Swap edges so that the new triangulation becomes Delauney
			if (remainDelaunay) {
				recSwapDelaunay(e12, n);
				recSwapDelaunay(e32, n);
				if (oldt2 != null) {
					recSwapDelaunay(e22, n);
					recSwapDelaunay(e42, n);
				}
			} else {
				if (!t1.areaLargerThan0()) {
					swap(e12);
				}
				if (!t3.areaLargerThan0()) {
					swap(e32);
				}
				if (oldt2 != null) {
					if (!t2.areaLargerThan0()) {
						swap(e22);
					}
					if (!t4.areaLargerThan0()) {
						swap(e42);
					}
				}
			}
		}
	}

}

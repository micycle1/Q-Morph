package meshditor;

import java.util.ArrayList;
import java.util.List;

/**
 * This class constitutes a simple implementation of the cleanup process as
 * outlined by Paul Kinney: "CleanUp: Improving Quadrilateral Finite Element
 * Meshes" (1997). Please note that this is not a complete and accurate
 * implementation, as it had to be somewhat adapted to work with the Q-Morph
 * implementation. E.g. the size cleanup is not implemented. Neither is cleaning
 * up bowties and the mesh topology inspection. Furthermore, the number of
 * cleanup patterns is limited to those described in Kinney's paper.<br>
 * <br>
 * TODO<br>
 * Draw all cases and test them. The ones running correctly: Connectivity -
 * stdCase1a - stdCase1b - stdCase2a - stdCase2b - stdCase3a - stdCase3b -
 * case1a og case1b - case2 - case3 - case4 - case5
 *
 * Boundary: - diamond - case1 - case2 - case3 - case4 - boundary angle > 150
 * degrees with one & two row transition
 *
 * Shape: - case1 - case2
 *
 * TODO: All good, but: - Why do I now see so little effect from class
 * GlobalSmooth???
 *
 * @author Karl Erik Levik
 *
 */
public class TopoCleanup extends GeomBasics {

	public TopoCleanup() {
	}

	/** A dart used when traversing the mesh in cleanup operations. */
	private Dart d;

	/** Initialize the object */
	public void init() {
		Quad q;
		Triangle tri;

		setCurMethod(this);

		elimChevsFinished = false;
		connCleanupFinished = false;
		boundaryCleanupFinished = false;
		shapeCleanupFinished = false;

		bcaseTriQFin = false;
		bcaseValPat1Fin = false;
		bcaseValPat2Fin = false;
		bcaseValPat3Fin = false;
		bcaseValPat4Fin = false;
		bcaseDiamondFin = false;

		shape1stTypeFin = false;
		passNum = 0;
		count = 0;

		// Ok... Remove fake quads and replace with triangles:
		for (int i = 0; i < elementList.size(); i++) {
			if (elementList.get(i) instanceof Quad) {
				q = (Quad) elementList.get(i);
				if (q.isFake) {
					tri = new Triangle(q.edgeList[base], q.edgeList[left], q.edgeList[right]);
					elementList.set(i, tri);
					q.disconnectEdges();
					tri.connectEdges();
				}
			}
		}

		deleteList = new ArrayList<>();
		vertices = new ArrayList<>(vertexList);
		d = new Dart();
	}

	/** Main loop for global topological clean-up */
	public void run() {
		int i, j;
		Quad q;
		Triangle tri;

		if (!step) {
			// Initial pass to eliminate chevrons
			while (!elimChevsFinished) {
				elimChevsStep();
			}

			// Then the major cleanup processes:
			for (j = 0; j < 3; j++) {

				// Perform connectivity cleanup:
				// Parse the list of vertices looking for cases that match. Fix these.
				connCleanupFinished = false;
				while (!connCleanupFinished) {
					connCleanupStep();
				}
				// Run some kind of global smooth
				globalSmooth();

				// Boundary cleanup:
				// Parse the list of elements looking for cases that match. Fix these.
				boundaryCleanupFinished = false;
				bcaseValPat1Fin = false;
				bcaseValPat2Fin = false;
				bcaseValPat3Fin = false;
				bcaseValPat4Fin = false;
				bcaseTriQFin = false;
				bcaseDiamondFin = false;
				while (!boundaryCleanupFinished) {
					boundaryCleanupStep();
				}
				// Run some kind of global smooth
				globalSmooth();

				// Shape cleanup:
				// Parse the list of elements looking for cases that match. Fix these.
				// Run a local smooth after each action.
				shapeCleanupFinished = false;
				shape1stTypeFin = false;
				while (!shapeCleanupFinished) {
					shapeCleanupStep();
				}

				// Run some kind of global smooth
				globalSmooth();
			}
		}

		if (doSmooth) {
			setCurMethod(globalSmooth);
		} else {
			setCurMethod(null);
		}
	}

	private boolean elimChevsFinished = false;
	private boolean connCleanupFinished = false;
	private boolean boundaryCleanupFinished = false;
	private boolean shapeCleanupFinished = false;

	private int passNum = 0;

	/** Method for stepping through the implementation one step at the time. */
	@Override
	public void step() {
		Element elem;
		int i, j;

		if (!elimChevsFinished) {
			elimChevsStep();
		} else if (!connCleanupFinished) {
			connCleanupStep();
			if (connCleanupFinished) {
				globalSmooth();
			}
		} else if (!boundaryCleanupFinished) {
			boundaryCleanupStep();
			if (boundaryCleanupFinished) {
				globalSmooth();
			}
		} else if (!shapeCleanupFinished) {
			shapeCleanupStep();
			if (shapeCleanupFinished) {
				globalSmooth();
			}
		} else {
			passNum++;
			if (passNum == 3) {
				passNum = 0;
				if (doSmooth) {
					setCurMethod(globalSmooth);
				} else {
					setCurMethod(null);
				}
			} else {
				connCleanupFinished = false;
				boundaryCleanupFinished = false;
				bcaseValPat1Fin = false;
				bcaseValPat2Fin = false;
				bcaseValPat3Fin = false;
				bcaseValPat4Fin = false;
				bcaseTriQFin = false;
				bcaseDiamondFin = false;
				shapeCleanupFinished = false;
				shape1stTypeFin = false;
			}
		}
	}

	int count = 0;
	List<Element> deleteList;
	List<Vertex> vertices;

	/** Initial pass to detect and cleanup chevrons. */
	private void elimChevsStep() {
		Element elem;
		Quad q = null;
		int i, j;

		while (count < elementList.size()) {
			elem = (Element) elementList.get(count);

			if (elem == null || !(elem instanceof Quad)) {
				count++;
			} else {
				q = (Quad) elem;
				if (q.isFake) {
					Msg.error("...Fake quad encountered!!!");
				}

				if (q.isChevron()) {
					eliminateChevron(q);
					count++;
					return;
				} else {
					count++;
				}
			}
		}

		for (i = 0; i < deleteList.size(); i++) {
			elem = deleteList.get(i);
			elementList.remove(elementList.indexOf(elem));
		}
		deleteList.clear();

		elimChevsFinished = true;
		count = 0;
	}

	/**
	 * The chevron is deleted along with one of its neighbors. A new Vertex is created
	 * and the deleted quads are replaced with three new ones surrounding the new
	 * Vertex. The neighbor chosen for deletion is the one that, when replaced by the
	 * new quads, yields the optimal Vertex valences.
	 * 
	 * @param q the chevron to be eliminated
	 */
	private void eliminateChevron(Quad q) {
		Element neighbor1, neighbor2, neighbor3, neighbor4;
		Quad q1, q2, qn1, qn2, qn3;
		Triangle tri;
		Edge e1, e2, e3, e4;
		Vertex n, n1, n2, n3;

		int i, j;
		int[] valenceAlt1 = new int[6];
		int[] valenceAlt2 = new int[6];
		int irrAlt1 = 0, irrAlt2 = 0, badAlt1 = 0, badAlt2 = 0;

		n = q.VertexAtLargestAngle();
		e3 = q.neighborEdge(n);
		n3 = e3.otherVertex(n);
		e4 = q.neighborEdge(n, e3);
		n1 = e4.otherVertex(n);

		e2 = q.neighborEdge(n3, e3);
		n2 = e2.otherVertex(n3);
		e1 = q.neighborEdge(n2, e2);

		neighbor1 = q.neighbor(e1);
		neighbor2 = q.neighbor(e2);

		if (!n.boundaryVertex()) {
			// Then the Vertex can be relocated so that the element is no longer a chevron
			Vertex nOld = new Vertex(n.x, n.y), nNew = n.laplacianSmooth();

			if (!n.equals(nNew)) {
				n.moveTo(nNew);
				inversionCheckAndRepair(n, nOld);
				if (!q.isChevron()) {
					n.update();
					return;
				} else {
					n.setXY(nOld.x, nOld.y);
					n.update();
				}
			}
		}

		if (neighbor1 != null && neighbor1 instanceof Quad && !((Quad) neighbor1).largestAngleGT180()) {
			q1 = (Quad) neighbor1;
		} else {
			q1 = null;
		}
		if (neighbor2 != null && neighbor2 instanceof Quad && !((Quad) neighbor2).largestAngleGT180()) {
			q2 = (Quad) neighbor2;
		} else {
			q2 = null;
		}

		if (q1 != null && q2 != null) {
			if (q.ang[q.angleIndex(n1)] + q1.ang[q1.angleIndex(n1)] < DEG_180) {
				valenceAlt1[0] = n.valence() + 1;
				valenceAlt1[1] = n1.valence() - 1;
				valenceAlt1[2] = n2.valence();
				valenceAlt1[3] = n3.valence();
				valenceAlt1[4] = q1.oppositeVertex(n2).valence() + 1;
				valenceAlt1[5] = q1.oppositeVertex(n1).valence();

				irrAlt1 = 1; // the new central Vertex will have valence 3
			} else { // then we must consider a fill_4(q, e1, n1):
				valenceAlt1[0] = n.valence() + 1;
				valenceAlt1[1] = n1.valence();
				valenceAlt1[2] = n2.valence();
				valenceAlt1[3] = n3.valence();
				valenceAlt1[4] = q1.oppositeVertex(n2).valence();
				valenceAlt1[5] = q1.oppositeVertex(n1).valence() + 1;

				irrAlt1 = 2; // the two new central vertices will each have valence 3
			}

			if (q.ang[q.angleIndex(n3)] + q2.ang[q2.angleIndex(n3)] < DEG_180) {
				valenceAlt2[0] = n.valence() + 1;
				valenceAlt2[1] = n1.valence();
				valenceAlt2[2] = n2.valence();
				valenceAlt2[3] = n3.valence() - 1;
				valenceAlt2[4] = q2.oppositeVertex(n2).valence() + 1;
				valenceAlt2[5] = q2.oppositeVertex(n3).valence();

				irrAlt2 = 1; // the new central Vertex will have valence 3
			} else { // then we must consider a fill_4(q, e2, n3):
				valenceAlt2[0] = n.valence() + 1;
				valenceAlt2[1] = n1.valence();
				valenceAlt2[2] = n2.valence();
				valenceAlt2[3] = n3.valence();
				valenceAlt2[4] = q2.oppositeVertex(n2).valence();
				valenceAlt2[5] = q2.oppositeVertex(n3).valence() + 1;

				irrAlt2 = 2; // the two new central vertices will each have valence 3
			}

			for (j = 0; j <= 5; j++) {
				if (valenceAlt1[j] != 4) {
					irrAlt1++;
				}
				if (valenceAlt1[j] < 3 || valenceAlt1[j] > 5) {
					badAlt1++;
				}

				if (valenceAlt2[j] != 4) {
					irrAlt2++;
				}
				if (valenceAlt2[j] < 3 || valenceAlt2[j] > 5) {
					badAlt2++;
				}
			}
		}

		if ((q1 != null && q2 == null) || (q1 != null && q2 != null && (badAlt1 < badAlt2 || (badAlt1 == badAlt2 && irrAlt1 <= irrAlt2)))) {

			deleteList.add(null);
			elementList.set(elementList.indexOf(q), null);
			deleteList.add(null);
			elementList.set(elementList.indexOf(q1), null);

			if (q.ang[q.angleIndex(n1)] + q1.ang[q1.angleIndex(n1)] < DEG_180) {
				fill3(q, e1, n2, true);
			} else {
				fill4(q, e1, n2); // n1
			}

		} else if (q2 != null) {
			deleteList.add(null); // q don't need any special treatment;
			elementList.set(elementList.indexOf(q), null);
			deleteList.add(null); // but q2 does, because it migth be at a later pos
			elementList.set(elementList.indexOf(q2), null);

			if (q.ang[q.angleIndex(n3)] + q2.ang[q2.angleIndex(n3)] < DEG_180) {
				fill3(q, e2, n2, true);
			} else {
				fill4(q, e2, n2); // n3
			}
		} else {
		}
	}

	/**
	 * Combine with neighbor and fill with "fill_3" as defined in the paper by
	 * P.Kinney Note that the method doesn't remove q and its neighbor from
	 * elementList.
	 * 
	 * @param q    the first quad
	 * @param e    the edge which is adjacent to both q and its neighbor
	 * @param n    a Vertex belonging to q, e, and one of the three new edges to be
	 *             created
	 * @param safe boolean indicating whether a safe pos must be attempted for the
	 *             new Vertex
	 */
	private Dart fill3(Quad q, Edge e, Vertex n, boolean safe) {
		Quad qn = (Quad) q.neighbor(e);
		Quad qn1, qn2, qn3;
		Edge b, l, r, t, ea, eb, ec;
		Vertex nOpp = q.oppositeVertex(n);
		Edge e2 = q.neighborEdge(n, e);
		Vertex eother = e.otherVertex(n), e2other = e2.otherVertex(n);
		Dart d = new Dart();

		Vertex newVertex;
		if (safe) {
			newVertex = e.midPoint();
		} else {
			newVertex = q.centroid();
		}

		newVertex.color = java.awt.Color.red; // creation in tCleanup
		ea = new Edge(nOpp, newVertex);
		eb = new Edge(n, newVertex);
		ec = new Edge(newVertex, qn.oppositeVertex(n));

		// Some minor updating...
		q.disconnectEdges();
		qn.disconnectEdges();

		ea.connectVertices();
		eb.connectVertices();
		ec.connectVertices();

		b = q.neighborEdge(e2other, e2);
		if (nOpp == b.rightVertex) {
			l = q.neighborEdge(b.leftVertex, b);
			r = ea;
		} else {
			l = ea;
			r = q.neighborEdge(b.rightVertex, b);
		}
		t = eb;
		qn1 = new Quad(b, l, r, t); // 1st replacement quad

		b = eb;
		if (ec.hasVertex(b.rightVertex)) { // b.leftVertex== n
			l = qn.neighborEdge(n, e);
			r = ec;
		} else {
			r = qn.neighborEdge(n, e);
			l = ec;
		}
		t = qn.oppositeEdge(e);
		qn2 = new Quad(b, l, r, t); // 2nd replacement quad

		b = q.neighborEdge(eother, e);
		if (b.leftVertex == nOpp) {
			l = ea;
			r = qn.neighborEdge(eother, e);
		} else {
			r = ea;
			l = qn.neighborEdge(eother, e);
		}
		t = ec;
		qn3 = new Quad(b, l, r, t); // 3rd replacement quad

		// remember to update the lists (vertexList, edgeList,
		// elementList, the vertices' edgeLists, ...
		qn1.connectEdges();
		qn2.connectEdges();
		qn3.connectEdges();

		elementList.add(qn1);
		elementList.add(qn2);
		elementList.add(qn3);

		e.disconnectVertices();
		edgeList.remove(edgeList.indexOf(e));

		edgeList.add(ea);
		edgeList.add(eb);
		edgeList.add(ec);

		vertexList.add(newVertex);
		vertices.add(newVertex);

		

		// Try smoothing the pos of newVertex:
		Vertex nOld = new Vertex(newVertex.x, newVertex.y), smoothed = newVertex.laplacianSmooth();
		if (!newVertex.equals(smoothed)) {
			newVertex.moveTo(smoothed);
			inversionCheckAndRepair(newVertex, nOld);
			newVertex.update();
		}

		d.elem = qn1;
		d.e = eb;
		d.n = n;

		return d;
	}

	/**
	 * Combine with neighbor and fill with "fill_4" as defined in paper by P.Kinney
	 * Note that the method doesn't remove q and its neighbor from elementList.
	 * 
	 * @param q  the first of the two quads to be combined
	 * @param e  the common edge of q and the second quad
	 * @param n2 one of the vertices of edge e and whose opposite Vertex in q will not
	 *           get connected to any new edge.
	 */
	private Dart fill4(Quad q, Edge e, Vertex n2) {
		Dart d = new Dart();
		Quad qn = (Quad) q.neighbor(e);
		// First get the vertices and edges in the two quads
		Edge temp;
		Vertex n5 = e.otherVertex(n2);

		Edge e1 = q.neighborEdge(n5, e);
		Vertex n0 = e1.otherVertex(n5);
		Edge e2 = q.neighborEdge(n0, e1);
		Vertex n1 = e2.otherVertex(n0);

		Edge e3 = q.neighborEdge(n1, e2);
		Edge e4 = qn.neighborEdge(n2, e);
		Vertex n3 = e4.otherVertex(n2);
		Edge e5 = qn.neighborEdge(n3, e4);
		Edge e6 = qn.neighborEdge(n5, e);
		Vertex n4 = e6.otherVertex(n5);

		// Create new vertices and edges
		Vertex n6 = e.midPoint();
		// This is new... hope it works...
		Vertex n7 = qn.centroid(); // new Vertex((n5.x + n3.x)*0.5, (n5.y + n3.y)*0.5);

		n6.color = java.awt.Color.red; // creation in tCleanup
		n7.color = java.awt.Color.red; // creation in tCleanup

		Edge eNew1 = new Edge(n0, n6);
		Edge eNew2 = new Edge(n2, n6);
		Edge eNew3 = new Edge(n6, n7);
		Edge eNew4 = new Edge(n5, n7);
		Edge eNew5 = new Edge(n3, n7);

		eNew1.connectVertices();
		eNew2.connectVertices();
		eNew3.connectVertices();
		eNew4.connectVertices();
		eNew5.connectVertices();

		// Create the new quads
		Edge l, r;
		if (eNew1.leftVertex == n0) {
			l = e2;
			r = eNew2;
		} else {
			r = e2;
			l = eNew2;
		}
		Quad qNew1 = new Quad(eNew1, l, r, e3);

		if (e1.leftVertex == n5) {
			l = eNew4;
			r = eNew1;
		} else {
			r = eNew4;
			l = eNew1;
		}
		Quad qNew2 = new Quad(e1, l, r, eNew3);

		if (e4.leftVertex == n3) {
			l = eNew5;
			r = eNew2;
		} else {
			r = eNew5;
			l = eNew2;
		}
		Quad qNew3 = new Quad(e4, l, r, eNew3);

		if (e6.leftVertex == n4) {
			l = e5;
			r = eNew4;
		} else {
			r = e5;
			l = eNew4;
		}
		Quad qNew4 = new Quad(e6, l, r, eNew5);

		// Update lists etc.
		e.disconnectVertices();
		q.disconnectEdges();
		qn.disconnectEdges();

		qNew1.connectEdges();
		qNew2.connectEdges();
		qNew3.connectEdges();
		qNew4.connectEdges();

		edgeList.remove(edgeList.indexOf(e));
		edgeList.add(eNew1);
		edgeList.add(eNew2);
		edgeList.add(eNew3);
		edgeList.add(eNew4);
		edgeList.add(eNew5);

		elementList.add(qNew1);
		elementList.add(qNew2);
		elementList.add(qNew3);
		elementList.add(qNew4);

		vertexList.add(n6);
		vertexList.add(n7);
		vertices.add(n6);
		vertices.add(n7);

		Vertex nOld = new Vertex(n6.x, n6.y), nNew = n6.laplacianSmooth();
		if (!n6.equals(nNew)) {
			n6.moveTo(nNew);
			inversionCheckAndRepair(n6, nOld);
			n6.update();
		}

		nOld = new Vertex(n7.x, n7.y);
		nNew = n7.laplacianSmooth();
		if (!n7.equals(nNew)) {
			n7.moveTo(nNew);
			inversionCheckAndRepair(n7, nOld);
			n7.update();
		}
		d.elem = qNew1;
		d.e = eNew2;
		d.n = n2;

		return d;
	}

	/**
	 * Replace the specified surrounding mesh elements (quads) with some other
	 * specified quads. This is accomplished by applying a composition of alpha
	 * iterators and mesh modify codes.
	 *
	 * @param startDart the dart where the composition is to start.
	 * @param fillPat   a composition of alpha iterators and mesh modification codes
	 *                  that moves around on the existing mesh and modifies it
	 *                  according to the codes in the composition. The format of the
	 *                  byte array is:<br>
	 *                  [Total number of bytes in array], [number of quads to be
	 *                  created], [number of quads to be deleted], [iterators and
	 *                  modify codes]. In addition to the normal alpha iterator
	 *                  codes 0,1,2, we have mesh modify codes:
	 *                  <ul>
	 *                  <li>Code 0 for alpha iterator 0
	 *                  <li>Code 1 for alpha iterator 1
	 *                  <li>Code 2 for alpha iterator 2
	 *                  <li>Code 3 for closing the current quad, new pos of cur.
	 *                  Vertex at the opposite Vertex
	 *                  <li>Code 4 for closing the current quad, new pos of cur.
	 *                  Vertex midway to oppos. Vertex
	 *                  <li>Code 5 for filling cur. quad and neighbour with fill_3
	 *                  <li>Code 6 for filling cur. quad and neighbour with fill_4
	 *                  <li>Code 7 for splitting cur. quad into two new quads along
	 *                  diagonal from cur. Vertex
	 *                  <li>Code 8 for switching cur. edge clockwise
	 *                  <li>Code 9 for switching cur. edge counter-clockwise
	 */
	private void applyComposition(Dart startDart, byte[] fillPat) {
		byte a;
		int qaIndex, qbIndex;
		d = startDart;
		for (int i = 1; i < fillPat[0]; i++) {
			a = fillPat[i];

			

			// Alpha iterators:
			if (a == 0) {
				d.n = d.e.otherVertex(d.n);
			} else if (a == 1) {
				d.e = d.elem.neighborEdge(d.n, d.e);
			} else if (a == 2) {
				d.elem = d.elem.neighbor(d.e);
			} else if (a == 3) {
				d = closeQuad((Quad) d.elem, d.e, d.n, false);
			} else if (a == 4) {
				d = closeQuad((Quad) d.elem, d.e, d.n, true);
			} else if (a == 5) {
				qaIndex = elementList.indexOf(d.elem);
				elementList.remove(qaIndex);
				qbIndex = elementList.indexOf(d.elem.neighbor(d.e));
				elementList.remove(qbIndex);
				d = fill3((Quad) d.elem, d.e, d.n, true);
			} else if (a == 6) {
				qaIndex = elementList.indexOf(d.elem);
				elementList.remove(qaIndex);
				qbIndex = elementList.indexOf(d.elem.neighbor(d.e));
				elementList.remove(qbIndex);
				d = fill4((Quad) d.elem, d.e, d.n);
			} else if (a == 7) {
				d = openQuad((Quad) d.elem, d.e, d.n);
			} else if (a == 8) {
				d = switchDiagonalCW((Quad) d.elem, d.e, d.n);
			} else if (a == 9) {
				d = switchDiagonalCCW((Quad) d.elem, d.e, d.n);
			} else {
				Msg.error("Illegal mesh modification code: " + a);
			}

			if (d == null) {
				return;
			}

			if (d.elem == null) {
				Msg.error("d.elem== null");
			}
			if (d.e == null) {
				Msg.error("d.e== null");
			}
			if (d.n == null) {
				Msg.error("d.n== null");
			}
		}
	}

	/*
	 * The different cases: Their valence patterns, vertex patterns, boundary
	 * patterns and compositions of mesh modification codes.
	 *
	 * In the valence patterns, the following codes apply: 4- (code 14) means a
	 * valence of 4 or less 4+ (code 24) means a valence of 4 or more 5 means a
	 * valence of 5 or more 0 means the valence is ignored and unchanged and is
	 * usually drawn as valence 4 1st value in pattern is total length of pattern
	 * 2nd value is the valence of the central Vertex The rest of the valences then
	 * follow in ccw order around the central Vertex
	 */

	/* The connectivity cases */
	static final byte[] stdCase1 = { 12, 5, 24, 3, 24, 3, 4, 0, 4, 0, 4, 3 };
	static final boolean[] stdVertexCase1 = { true, false, true, false, false, true, false, true, false, false }; // ok...
	static final byte[] stdComp1 = { 4, 5, 1, 9 };
	static final byte[] stdCase2a = { 14, 6, 14, 24, 4, 3, 4, 24, 14, 3, 24, 3, 24, 3 };
	static final boolean[] stdVertexCase2a = { false, true, false, false, false, true, false, false, true, false, true, false }; // ok...
	static final byte[] stdComp2a = { 22, 2, 1, 5, 2, 1, 8, 2, 1, 0, 2, 1, 0, 5, 0, 2, 1, 2, 0, 1, 9, 5 };
	static final byte[] stdCase2b = { 14, 6, 24, 3, 4, 3, 24, 3, 4, 0, 4, 0, 4, 3 };
	static final boolean[] stdVertexCase2b = { true, false, false, false, true, false, false, true, false, true, false, false }; // ok...
	static final byte[] stdComp2b = { 9, 5, 1, 2, 1, 5, 1, 0, 5 };
	static final byte[] stdCase3a = { 12, 5, 24, 3, 4, 0, 4, 0, 4, 0, 4, 3 };
	static final boolean[] stdVertexCase3a = { true, false, false, true, false, true, false, true, false, false }; // ok...
	static final byte[] stdComp3a = { 2, 5 };
	static final byte[] stdCase3b = { 12, 5, 24, 0, 4, 3, 24, 3, 24, 3, 24, 3 };
	static final boolean[] stdVertexCase3b = { false, true, false, false, true, false, true, false, true, false }; // ok...
	static final byte[] stdComp3b = { 18, 1, 2, 1, 0, 8, 5, 0, 1, 2, 0, 1, 9, 2, 1, 0, 1, 9 };
	static final byte[] case1a = { 12, 5, 3, 4, 4, 3, 0, 0, 0, 0, 0, 0 }; // mirror of case1b
	static final byte[] comp1a = { 3, 1, 9 };
	static final byte[] case1b = { 12, 5, 4, 4, 3, 0, 0, 0, 0, 0, 0, 3 }; // mirror of case1a
	static final byte[] comp1b = { 2, 8 };
	static final byte[] case2 = { 12, 5, 3, 4, 24, 4, 3, 0, 0, 0, 0, 0 };
	static final byte[] comp2 = { 4, 1, 0, 5 };
	static final byte[] case3 = { 8, 3, 24, 3, 24, 0, 0, 0 };
	static final boolean[] vertexCase3 = { false, false, false, true, false, true };
	static final boolean[] ipat3 = { true, true, true, false, false, false };
	static final byte[] comp3 = { 2, 4 };
	static final byte[] case4 = { 8, 3, 3, 5, 4, 5, 4, 4 };
	static final boolean[] ipat4 = { true, false, false, false, false, false };
	static final byte[] comp4 = { 5, // Consider 4 instead of 3... nah, 3 is best here
			0, 3, 0, 3 };
	static final byte[] case5 = { 10, 4, 3, 4, 4, 3, 4, 4, 4, 5 };
	static final boolean[] ipat5 = { true, true, true, true, true, false, false, false };
	static final byte[] comp5 = { 13, // Consider 4 instead of 3...
			4, 0, 3, 2, 0, 1, 0, 2, 1, 4, 0, 3 };

	/* The boundary cases */
	static final byte[] bcase1a = { 9, 5, 4, 3, 5, 4, 3, 4, 4 };
	static final boolean[] bpat1a = { true, true, false, true, true, false, false, true };
	static final byte[] bcomp1a = { 3, 1, 8 };
	static final byte[] bcase1b = { 9, 5, 4, 24, 3, 4, 5, 3, 4 };
	static final boolean[] bpat1b = { true, true, false, false, true, true, false, true };
	static final byte[] bcomp1b = { 5, 1, 2, 1, 9 };
	static final byte[] bcase2a = { 7, 4, 4, 3, 5, 4, 4 };
	static final boolean[] bpat2a = { true, true, false, true, true, true };
	static final byte[] bcomp2a = { 3, 1, 5 };
	static final byte[] bcase2b = { 7, 4, 4, 4, 5, 3, 4 };
	static final boolean[] bpat2b = { true, true, true, true, false, true };
	static final byte[] bcomp2b = { 3, 1, 5 };

	static final byte[] bcase3 = { 10, 4, 5, 4, 3, 4, 5, 4, 3, 4 };
	static final boolean[] bpat3 = { false, true, false, false, true, true, false, false, true };
	static final byte[] bcomp3 = { 12, 8, 2, 0, 1, 0, 2, 3, 1, 0, 2, 3 };
	// 3,1,8,4};
	static final byte[] bcase4 = { 9, 5, 4, 3, 5, 3, 5, 3, 4 };
	static final boolean[] bpat4 = { true, true, false, true, true, true, false, true };
	static final byte[] bcomp4 = { 7, 1, 2, 1, 5, 1, 8 };

	/** Perform one more step of connectivity cleanup. */
	private void connCleanupStep() {
		int i, vInd;
		Vertex c = null;
		Element elem;
		Dart d;
		byte[] pattern;
		Vertex[] ccwNeighbors = null;
		double[] angles;

		// First check for the standard patterns:
		for (i = 0; i < vertices.size(); i++) {
			c = (Vertex) vertices.get(i);
			if (c == null || c.boundaryOrTriangleVertex()) {
				continue;
			}
			ccwNeighbors = c.ccwSortedNeighbors();
			c.createValencePattern(ccwNeighbors);
			if (c.irregNeighborVertices() <= 2) {
				continue;
			}

			angles = c.surroundingAngles(ccwNeighbors, c.pattern[0] - 2);

			if ((vInd = c.patternMatch(stdCase1, stdVertexCase1, angles)) != -1) {
				d = getDartAt(c, ccwNeighbors, vInd - 2);
				applyComposition(d, stdComp1);

				vInd = vertices.indexOf(c);
				if (vInd != -1) {
					vertices.remove(vInd);
				}
				addVertices(ccwNeighbors, c.pattern[0] - 2);
			} else if ((vInd = c.patternMatch(stdCase2a, stdVertexCase2a, angles)) != -1) {
				d = getDartAt(c, ccwNeighbors, vInd - 2);
				applyComposition(d, stdComp2a);

				vInd = vertices.indexOf(c);
				if (vInd != -1) {
					vertices.remove(vInd);
				}
				addVertices(ccwNeighbors, c.pattern[0] - 2);
			} else if ((vInd = c.patternMatch(stdCase2b, stdVertexCase2b, angles)) != -1) {
				d = getDartAt(c, ccwNeighbors, vInd - 2);
				applyComposition(d, stdComp2b);

				vInd = vertices.indexOf(c);
				if (vInd != -1) {
					vertices.remove(vInd);
				}
				addVertices(ccwNeighbors, c.pattern[0] - 2);
			} else if ((vInd = c.patternMatch(stdCase3a, stdVertexCase3a, angles)) != -1) {
				d = getDartAt(c, ccwNeighbors, vInd - 2);
				applyComposition(d, stdComp3a);

				vInd = vertices.indexOf(c);
				if (vInd != -1) {
					vertices.remove(vInd);
				}
				addVertices(ccwNeighbors, c.pattern[0] - 2);
			} else if ((vInd = c.patternMatch(stdCase3b, stdVertexCase3b, angles)) != -1) {
				d = getDartAt(c, ccwNeighbors, vInd - 2);
				applyComposition(d, stdComp3b);

				vInd = vertices.indexOf(c);
				if (vInd != -1) {
					vertices.remove(vInd);
				}
				addVertices(ccwNeighbors, c.pattern[0] - 2);
			}
		}

		// Then check for the other patterns:
		for (i = 0; i < vertices.size(); i++) {
			c = (Vertex) vertices.get(i);
			if (c == null || c.boundaryOrTriangleVertex()) {
				continue;
			}
			ccwNeighbors = c.ccwSortedNeighbors();
			c.createValencePattern(ccwNeighbors);
			if (c.irregNeighborVertices() <= 2) {
				continue;
			}

			angles = c.surroundingAngles(ccwNeighbors, c.pattern[0] - 2);

			if ((vInd = c.patternMatch(case1a)) != -1) {
				d = getDartAt(c, ccwNeighbors, vInd - 2);
				if (d != null) {
					applyComposition(d, comp1a);
					vInd = vertices.indexOf(c);
					if (vInd != -1) {
						vertices.remove(vInd);
					}
					addVertices(ccwNeighbors, c.pattern[0] - 2);
				}
			} else if ((vInd = c.patternMatch(case1b)) != -1) {
				d = getDartAt(c, ccwNeighbors, vInd - 2);
				if (d != null) {
					applyComposition(d, comp1b);
					vInd = vertices.indexOf(c);
					if (vInd != -1) {
						vertices.remove(vInd);
					}
					addVertices(ccwNeighbors, c.pattern[0] - 2);
				}
			} else if ((vInd = c.patternMatch(case2)) != -1) {
				d = getDartAt(c, ccwNeighbors, vInd - 2);
				if (d != null) {
					applyComposition(d, comp2);
					vInd = vertices.indexOf(c);
					if (vInd != -1) {
						vertices.remove(vInd);
					}
					addVertices(ccwNeighbors, c.pattern[0] - 2);
				}
			} else if ((vInd = c.patternMatch(case3, vertexCase3, angles)) != -1
					&& internalVertices(ipat3, ccwNeighbors, vInd - 2, c.pattern[0] - 2)) {
				d = getDartAt(c, ccwNeighbors, vInd - 2);
				if (d != null) {
					applyComposition(d, comp3);
					vInd = vertices.indexOf(c);
					if (vInd != -1) {
						vertices.remove(vInd);
					}
					addVertices(ccwNeighbors, c.pattern[0] - 2);
				}
			} else if ((vInd = c.patternMatch(case4)) != -1 && internalVertices(ipat4, ccwNeighbors, vInd - 2, c.pattern[0] - 2)) {
				d = getDartAt(c, ccwNeighbors, vInd - 2);
				if (d != null) {
					applyComposition(d, comp4);
					vInd = vertices.indexOf(c);
					if (vInd != -1) {
						vertices.remove(vInd);
					}
					addVertices(ccwNeighbors, c.pattern[0] - 2);
				}
			} else if ((vInd = c.patternMatch(case5)) != -1 && internalVertices(ipat5, ccwNeighbors, vInd - 2, c.pattern[0] - 2)) {
				d = getDartAt(c, ccwNeighbors, vInd - 2);
				if (d != null) {
					applyComposition(d, comp5);
					vInd = vertices.indexOf(c);
					if (vInd != -1) {
						vertices.remove(vInd);
					}
					addVertices(ccwNeighbors, c.pattern[0] - 2);
				}
			}
		}

		vertices = new ArrayList<>(vertexList);
		connCleanupFinished = true;
		count = 0;
	}

	/** Method to confirm whether marked vertices are truely internal . */
	private boolean internalVertices(boolean[] ipat, Vertex[] ccwVertices, int index, int len) {
		int i = index, j = 0;
		while (j < len) {
			if (ipat[j] && ccwVertices[i].boundaryVertex()) {
				return false;
			}
			j++;
			i++;
			if (i == len) {
				i = 0;
			}
		}
		return true;
	}

	private void addVertices(Vertex[] arrayOfVertices, int len) {
		Vertex n;
		int j;
		for (int i = 0; i < len; i++) {
			n = arrayOfVertices[i];
			if (n != null) {
				j = vertices.indexOf(n);
				if (j == -1) {
					j = vertexList.indexOf(n);
					if (j != -1) {
						vertices.add(n);
					}
				}
			}
		}
	}

	private boolean bcaseTriQFin = false;
	private boolean bcaseValPat1Fin = false;
	private boolean bcaseValPat2Fin = false;
	private boolean bcaseValPat3Fin = false;
	private boolean bcaseValPat4Fin = false;
	private boolean bcaseDiamondFin = false;

	/** Perform one more steps of boundary cleanup. */
	private void boundaryCleanupStep() {
		int i, j, index;
		Element elem;
		Triangle tri;
		Quad q = null, pq = null;
		Vertex n, n0, n1 = null, n2, n3;
		Edge e1 = null, e2 = null, ep = null, e;
		Vertex[] ccwNeighbors;
		Dart d;

		if (!bcaseValPat1Fin) {

			for (i = 0; i < vertices.size(); i++) {
				n1 = (Vertex) vertices.get(i);
				if (n1 == null) {
					continue;
				}
				ccwNeighbors = n1.ccwSortedNeighbors();
				if (ccwNeighbors == null) {
					continue;
				}
				n1.createValencePattern((byte) (n1.edgeList.size() * 2 - 1), ccwNeighbors);

				if (n1.boundaryPatternMatch(bcase1a, bpat1a, ccwNeighbors)) {
					d = getDartAt(n1, ccwNeighbors, 0);
					applyComposition(d, bcomp1a);

					j = vertices.indexOf(n1);
					if (j != -1) {
						vertices.remove(j);
					}
					addVertices(ccwNeighbors, n1.pattern[0] - 2);
				} else if (n1.boundaryPatternMatch(bcase1b, bpat1b, ccwNeighbors)) {
					d = getDartAt(n1, ccwNeighbors, 0);
					applyComposition(d, bcomp1b);

					j = vertices.indexOf(n1);
					if (j != -1) {
						vertices.remove(j);
					}
					addVertices(ccwNeighbors, n1.pattern[0] - 2);
				}
			}
			// vertices= (ArrayList)vertexList.clone();
			bcaseValPat1Fin = true;
			return;
		} else if (!bcaseValPat2Fin) {

			for (i = 0; i < vertices.size(); i++) {
				n1 = (Vertex) vertices.get(i);
				if (n1 == null) {
					continue;
				}
				ccwNeighbors = n1.ccwSortedNeighbors();
				if (ccwNeighbors == null) {
					continue;
				}
				n1.createValencePattern((byte) (n1.edgeList.size() * 2 - 1), ccwNeighbors);

				if (n1.boundaryPatternMatch(bcase2a, bpat2a, ccwNeighbors)) {
					d = getDartAt(n1, ccwNeighbors, 0);
					applyComposition(d, bcomp2a);

					j = vertices.indexOf(n1);
					if (j != -1) {
						vertices.remove(j);
					}
					addVertices(ccwNeighbors, n1.pattern[0] - 2);
				} else if (n1.boundaryPatternMatch(bcase2b, bpat2b, ccwNeighbors)) {
					d = getDartAt(n1, ccwNeighbors, 0);
					applyComposition(d, bcomp2b);

					j = vertices.indexOf(n1);
					if (j != -1) {
						vertices.remove(j);
					}
					addVertices(ccwNeighbors, n1.pattern[0] - 2);
				}
			}
			// vertices= (ArrayList)vertexList.clone();
			bcaseValPat2Fin = true;
			return;
		} else if (!bcaseValPat3Fin) {

			for (i = 0; i < vertices.size(); i++) {
				n1 = (Vertex) vertices.get(i);
				if (n1 == null) {
					continue;
				}
				if (n1.boundaryVertex()) {
					continue;
				}

				ccwNeighbors = n1.ccwSortedNeighbors();
				if (ccwNeighbors == null) {
					continue;
				}
				n1.createValencePattern((byte) (n1.edgeList.size() * 2), ccwNeighbors);

				if ((index = n1.boundaryPatternMatchSpecial(bcase3, bpat3, ccwNeighbors)) != -1) {
					d = getDartAt(n1, ccwNeighbors, index - 2);
					applyComposition(d, bcomp3);

					j = vertices.indexOf(n1);
					if (j != -1) {
						vertices.remove(j);
					}
					addVertices(ccwNeighbors, n1.pattern[0] - 2);
				}
			}
			// vertices= (ArrayList)vertexList.clone();
			bcaseValPat3Fin = true;
			return;
		} else if (!bcaseValPat4Fin) {
			for (i = 0; i < vertices.size(); i++) {
				n1 = (Vertex) vertices.get(i);
				if (n1 == null) {
					continue;
				}
				ccwNeighbors = n1.ccwSortedNeighbors();
				if (ccwNeighbors == null) {
					continue;
				}
				n1.createValencePattern((byte) (n1.edgeList.size() * 2 - 1), ccwNeighbors);

				if (n1.boundaryPatternMatch(bcase4, bpat4, ccwNeighbors)) {
					d = getDartAt(n1, ccwNeighbors, 0);
					applyComposition(d, bcomp4);

					j = vertices.indexOf(n1);
					if (j != -1) {
						vertices.remove(j);
					}
					addVertices(ccwNeighbors, n1.pattern[0] - 2);
				}
			}
			// vertices= (ArrayList)vertexList.clone();
			bcaseValPat4Fin = true;
			return;
		} else if (!bcaseTriQFin) {
			for (i = 0; i < vertices.size(); i++) {
				n1 = (Vertex) vertices.get(i);

				if (n1 == null) {
					continue;
				}

				if (n1.boundaryVertex()) {

					e1 = (Edge) n1.edgeList.get(0);
					elem = e1.element1;
					if (!(elem instanceof Quad)) {
						continue;
					}
					q = (Quad) elem;

					if (n1.edgeList.size() == 2 && q.ang[q.angleIndex(n1)] > DEG_150) {
						e2 = q.neighborEdge(n1, e1);
						if (e1.boundaryEdge() && e2.boundaryEdge()) {
							Quad q3 = null, q4 = null, q33 = null, q44 = null, qNew, qn = null;
							Edge e3, e4, e33, e44;
							n0 = e1.otherVertex(n1);
							n2 = e2.otherVertex(n1);
							e3 = q.neighborEdge(n2, e2);
							n3 = e3.otherVertex(n2);
							e4 = q.neighborEdge(n3, e3);

							elem = q.neighbor(e3);
							if (elem instanceof Quad) {
								q3 = (Quad) elem;
							}

							elem = q.neighbor(e4);
							if (elem instanceof Quad) {
								q4 = (Quad) elem;
							}

							if (q3 != null && q4 != null) {

								e33 = q3.neighborEdge(n3, e3);
								e44 = q4.neighborEdge(n3, e4);

								elem = q3.neighbor(e33);
								if (elem instanceof Quad) {
									q33 = (Quad) elem;
								}

								elem = q4.neighbor(e44);
								if (elem instanceof Quad) {
									q44 = (Quad) elem;
								}

								if (q33 != null) {
									e = q33.neighborEdge(n3, e33);
									elem = q33.neighbor(e);
									if (elem instanceof Quad) {
										qn = (Quad) elem;
									}
								}

								if (q33 == null || q44 == null || qn != q44) {
									// One row transition
									elementList.remove(elementList.indexOf(q));
									elementList.remove(elementList.indexOf(q3));
									elementList.remove(elementList.indexOf(q4));

									fill4(q, e4, n3);
									qNew = (Quad) q3.neighbor(e3);
									fill3(q3, e3, n2, true);
									elementList.remove(elementList.indexOf(qNew));

									j = vertices.indexOf(n1);
									if (j != -1) {
										vertices.remove(j);
									}
								} else if (qn == q44) {
									// Two row transition
									elementList.remove(elementList.indexOf(q));
									elementList.remove(elementList.indexOf(q4));
									elementList.remove(elementList.indexOf(q44));

									fill3(q4, e44, e44.otherVertex(n3), true);
									qNew = (Quad) q.neighbor(e4);
									fill3(q, e4, n3, true);

									elementList.remove(elementList.indexOf(qNew));

									j = vertices.indexOf(n1);
									if (j != -1) {
										vertices.remove(j);
									}
								}
							}
						}
					}
				}
			}
			bcaseTriQFin = true;
			return;
		} else if (!bcaseDiamondFin) {

			for (i = 0; i < vertices.size(); i++) {
				n1 = (Vertex) vertices.get(i);

				if (n1 == null) {
					continue;
				}

				if (n1.boundaryVertex() && n1.valence() > 4) {
					

					// First find a quad having an edge at the boundary at Vertex n1,
					// or if this does not exist, the find the first quad when looking
					// from the boundary an inwards around Vertex n1:
					e1 = n1.anotherBoundaryEdge(null);
					if (e1.element1 instanceof Quad) {
						q = (Quad) e1.element1;
						pq = q;
						ep = e1;
					} else {
						tri = (Triangle) e1.element1;
						e1 = e1.nextQuadEdgeAt(n1, e1.element1);
						if (e1 != null) {
							q = e1.getQuadElement();
							pq = q;
							ep = e1;
						}
					}

					if (q != null) {
					}

					// Then parse each quad from one boundary edge to the other until
					// a match is made:
					int val, n2val, n3val, maxdev;
					while (q != null) {
						if (q.boundaryDiamond()) {

							n2 = e1.otherVertex(n1);
							e2 = q.neighborEdge(n1, e1);
							n3 = e2.otherVertex(n1);

							n2val = n2.valence();
							n3val = n3.valence();
							maxdev = Math.max(Math.abs(4 - n2val), Math.abs(4 - n3val));
							val = n2val + n3val;

							if (Math.abs(4 - (val - 2)) <= maxdev) {
								d = closeQuad(q, e1, e1.otherVertex(n1), true);
								if (d != null) {
									q = pq;
									e1 = ep;
								}
							}
						}
						pq = q;
						ep = e1;
						e1 = e1.nextQuadEdgeAt(n1, q);
						if (e1 != null) {
							q = (Quad) q.neighbor(e1);
						} else {
							break;
						}
					}
				}
			}
			bcaseDiamondFin = true;
			return;
		} else {
			vertices = new ArrayList<>(vertexList);
			boundaryCleanupFinished = true;
			count = 0;
			return;
		}
	}

	private boolean shape1stTypeFin = false;

	/** The shape cleanup */
	private void shapeCleanupStep() {
		Element elem;
		Quad q = null, q2 = null, qo = null, qtemp;
		Edge e1, e2, e3, e4, eo;
		Vertex n, n1, n2, n3, n4, nqOpp, nq2Opp, noOpp;
		int i, j;
		double ang, ang1 = 0, ang2 = 0, ango = 0, angtmp, q2angn3;

		if (!shape1stTypeFin) {
			for (i = 0; i < vertices.size(); i++) {
				n = (Vertex) vertices.get(i);

				if (n == null) {
					continue;
				}

				if (n.boundaryVertex()) {

					if (n.edgeList.size() == 3) {

						e1 = n.anotherBoundaryEdge(null);
						elem = e1.element1;
						if (!(elem instanceof Quad)) {
							continue;
						}
						q = (Quad) elem;
						ang = q.ang[q.angleIndex(n)];

						e2 = q.neighborEdge(n, e1);
						elem = q.neighbor(e2);
						if (elem instanceof Quad) {
							q2 = (Quad) elem;

							ang2 = q2.ang[q2.angleIndex(n)];
							if (ang2 > DEG_160 && ang2 > ang) {
								qtemp = q;
								q = q2;
								q2 = qtemp;
								angtmp = ang;
								ang = ang2;
								ang2 = angtmp;

								e1 = q.neighborEdge(n, e2);
							}
						}

						if (ang < DEG_160) {
							continue;
						}

						n1 = e1.otherVertex(n);
						n2 = e2.otherVertex(n);
						nqOpp = q.oppositeVertex(n);

						ang1 = q.ang[q.angleIndex(n1)];

						if (!n1.boundaryVertex() || !n2.boundaryVertex() || !nqOpp.boundaryVertex()) {
							continue;
						}

						eo = q.neighborEdge(n1, e1);

						elem = q.neighbor(eo);
						if (elem instanceof Quad) {
							qo = (Quad) elem;
							ango = qo.ang[qo.angleIndex(n1)];
						}

						if (q2 == null || qo == null) {
							continue;
						}

						e3 = q2.neighborEdge(n, e2);
						n3 = e3.otherVertex(n);
						q2angn3 = q2.ang[q2.angleIndex(n3)];
						nq2Opp = q2.oppositeVertex(n);

						e4 = qo.neighborEdge(n1, eo);
						n4 = e4.otherVertex(n1);
						noOpp = qo.oppositeVertex(n1);

						if (ang2 != 0 && q2angn3 > ango && n3.boundaryVertex() && nq2Opp.boundaryVertex()) {

							elementList.remove(elementList.indexOf(q));
							elementList.remove(elementList.indexOf(q2));
							fill4(q, e2, n);
						} else if (ango != 0 && ango > q2angn3 && n4.boundaryVertex() && nq2Opp.boundaryVertex()) {

							elementList.remove(elementList.indexOf(q));
							elementList.remove(elementList.indexOf(qo));
							fill4(qo, eo, n1);
						}
					} else if (n.edgeList.size() == 4) {
						e1 = n.anotherBoundaryEdge(null);
						elem = e1.element1;
						if (elem instanceof Quad) {
							q2 = (Quad) elem;
						} else {
							continue;
						}

						e2 = q2.neighborEdge(n, e1);
						elem = q2.neighbor(e2);
						if (elem instanceof Quad) {
							q = (Quad) elem;
						} else {
							continue;
						}

						if (q.ang[q.angleIndex(n)] < DEG_160) {
							continue;
						}

						e3 = q.neighborEdge(n, e2);
						elem = q.neighbor(e3);
						if (elem instanceof Quad) {
							qo = (Quad) elem;
						} else {
							continue;
						}
						e4 = qo.neighborEdge(n, e3);

						n1 = e1.otherVertex(n);
						n2 = e2.otherVertex(n);
						n3 = e3.otherVertex(n);
						n4 = e4.otherVertex(n);

						if (!n1.boundaryVertex() || !n2.boundaryVertex() || !n3.boundaryVertex() || !n4.boundaryVertex()) {
							continue;
						}

						openQuad(q, e2, n);
						elementList.remove(elementList.indexOf(q2.neighbor(e2)));
						elementList.remove(elementList.indexOf(qo.neighbor(e3)));
						fill3(q2, e2, n2, true);
						fill3(qo, e3, n3, true);

						elementList.remove(elementList.indexOf(q2));
						elementList.remove(elementList.indexOf(qo));
					}
				}
			}
			shape1stTypeFin = true;
			return;
		}

		while (count < elementList.size()) {
			elem = (Element) elementList.get(count);

			if (elem == null || !(elem instanceof Quad)) {
				count++;
			} else {
				q = (Quad) elem;
				if (q.isFake) {
					Msg.error("...Fake quad encountered!!!");
				}

				if (q.isChevron()) {
					eliminateChevron(q);
					count++;
				} else {
					count++;
				}
			}
		}

		for (i = 0; i < deleteList.size(); i++) {
			elem = (Element) deleteList.get(i);
			elementList.remove(elementList.indexOf(elem));
		}
		deleteList.clear();

		vertices = new ArrayList<>(vertexList);
		shapeCleanupFinished = true;
		count = 0;
	}

	/**
	 * Return the dart with Vertex c, the edge connecting c and the Vertex at pos. i in
	 * neighbors, and the quad with that edge and Vertex at pos. i+1
	 * 
	 * @param c         the central Vertex
	 * @param neighbors array of neighboring vertices to c
	 * @param i         index into neighbors
	 */
	private Dart getDartAt(Vertex c, Vertex[] neighbors, int i) {
		Edge e = c.commonEdge(neighbors[i]);
		if (e == null) {
			return null;
		}
		Quad q1 = (Quad) e.element1;
		Quad q2 = (Quad) e.element2;

		if (q1.hasVertex(neighbors[i + 1])) {
			return new Dart(c, e, q1);
		} else if (q2.hasVertex(neighbors[i + 1])) {
			return new Dart(c, e, q2);
		} else {
			return null;
		}
	}

	/**
	 * Collapse a quad by joining two and two of its consecutive edges.
	 * 
	 * @param q        the quad to be collapsed
	 * @param e1       an edge of q that has the Vertex nK
	 * @param nK       the Vertex that is to be joined with its opposite Vertex in q
	 * @param centroid boolean indicating whether to look for a new pos for the
	 *                 joined vertices somewhere between the original positions,
	 *                 starting at the centroid of q, or to unconditionally try
	 *                 using the position of the Vertex in q which is opposite to nK.
	 * @return the new current dart.
	 */
	private Dart closeQuad(Quad q, Edge e1, Vertex nK, boolean centroid) {
		Dart d = new Dart();
		Element nElem = q.neighbor(e1); // Save for later...
		Vertex nKOpp = q.oppositeVertex(nK);
		Vertex nKp1 = e1.otherVertex(nK); // , nKm1= eKm1.otherVertex(nK);
		Edge e2 = q.neighborEdge(nKp1, e1), e4 = q.neighborEdge(nK, e1);

		List<Element> lK = nK.adjElements();
		List<Element> lKOpp = nKOpp.adjElements();
		Vertex n = null;
		int i;

		if (centroid) {
			n = safeNewPosWhenCollapsingQuad(q, nK, nKOpp);
			if (n == null) {
				return null;
			}
		} else if (q.anyInvertedElementsWhenCollapsed(nKOpp, nK, nKOpp, lK, lKOpp)) {
			return null;
		} else {
			n = nKOpp;
		}
		elementList.remove(elementList.indexOf(q));

		edgeList.remove(edgeList.indexOf(e1)); // e2
		edgeList.remove(edgeList.indexOf(q.neighborEdge(nK, e1))); // (nKOpp, e2)
		q.disconnectEdges();
		q.closeQuad(e2, e1); // (e1,e2)

		nKOpp.setXY(n); // nK.setXY(n);
		vertexList.remove(vertexList.indexOf(nK)); // nKOpp
		i = vertices.indexOf(nK);
		if (i != -1) {
			vertices.set(i, null); // nKOpp
		}
		nKOpp.update(); // nK.update();

		d.elem = nElem;
		d.e = e2; // e1;
		d.n = nKOpp; // nK;

		return d;
	}

	/**
	 * Create a new quad by "opening" one at the specified Vertex inside the specified
	 * quad. This effectively results in a splitting of the specified quad.
	 * 
	 * @param q  the quad
	 * @param e  an edge in q (not used by method), but contained in returned dart
	 * @param n1 split q along the diagonal between n1 and its opposite Vertex
	 */
	private Dart openQuad(Quad q, Edge e, Vertex n1) {
		Dart d = new Dart();
		Vertex c = q.centroid();
		Edge e1 = q.neighborEdge(n1);
		Edge e2 = q.neighborEdge(n1, e1);

		Vertex n2 = e2.otherVertex(n1);
		Edge e3 = q.neighborEdge(n2, e2);
		Vertex n3 = e3.otherVertex(n2);
		Edge e4 = q.neighborEdge(n3, e3);
		Vertex n4 = e4.otherVertex(n3);

		Edge e1New = new Edge(c, n1);
		Edge e4New = new Edge(c, n3);

		e1New.connectVertices();
		e4New.connectVertices();

		q.replaceEdge(e1, e1New);
		q.replaceEdge(e4, e4New);
		e1.disconnectFromElement(q);
		e4.disconnectFromElement(q);

		e1New.connectToQuad(q);
		e4New.connectToQuad(q);

		Quad qNew;
		if (e1.leftVertex == n1) {
			qNew = new Quad(e1, e1New, e4, e4New);
		} else {
			qNew = new Quad(e1, e4, e1New, e4New);
		}

		qNew.connectEdges();

		edgeList.add(e1New);
		edgeList.add(e4New);

		elementList.add(qNew);
		c.color = java.awt.Color.red; // Indicate it was created during clean-up
		vertexList.add(c);
		vertices.add(c);

		q.updateLR();

		d.elem = q;
		d.e = e;
		d.n = n1;

		return d;
	}

	/**
	 * Create 2 new quads from 2 specified quads, in which the common edge of the
	 * given quads has been rotated one step in the CCW direction. Delete the old
	 * quads.
	 * 
	 * @param qa  one of the two quads adjacent the edge to be switched, e1a
	 * @param e1a the edge to be switched
	 * @param n   one of e1a's vertices
	 * @return a dart representing the input dart after the operation is performed.
	 */
	private Dart switchDiagonalCCW(Quad qa, Edge e1a, Vertex n) {
		Dart d = new Dart();
		Vertex n1a, n2a, n3a, n4a, n1b, n2b, n3b, n4b;
		Edge e2a, e3a, e4a, e1b, e2b, e3b, e4b;
		Edge eNew, l, r;
		Quad q1, q2;

		Quad qb = (Quad) qa.neighbor(e1a);
		int qaIndex = elementList.indexOf(qa), qbIndex = elementList.indexOf(qb);

		// First get the edges of qa in ccw order:
		n2a = qa.nextCCWVertex(e1a.leftVertex);
		if (n2a == e1a.rightVertex) {
			n1a = e1a.leftVertex;
		} else {
			n1a = e1a.rightVertex;
			n2a = e1a.leftVertex;
		}

		e2a = qa.neighborEdge(n2a, e1a);
		n3a = e2a.otherVertex(n2a);
		e3a = qa.neighborEdge(n3a, e2a);
		n4a = e3a.otherVertex(n3a);
		e4a = qa.neighborEdge(n4a, e3a);

		// Now get the edges of qb in ccw order:
		e1b = e1a;
		n2b = qb.nextCCWVertex(e1b.leftVertex);
		if (n2b == e1b.rightVertex) {
			n1b = e1b.leftVertex;
		} else {
			n1b = e1b.rightVertex;
			n2b = e1b.leftVertex;
		}
		e2b = qb.neighborEdge(n2b, e1b);
		n3b = e2b.otherVertex(n2b);
		e3b = qb.neighborEdge(n3b, e2b);
		n4b = e3b.otherVertex(n3b);
		e4b = qb.neighborEdge(n4b, e3b);

		Vertex nOld, smoothed;
		// Check to see if the switch will violate the mesh topology:
		if (e4a.sumAngle(qa, n1a, e2b) >= Math.PI) { // if angle >= 180 degrees...
			if (n1a.boundaryVertex()) { // exit if Vertex on boundary
				return null;
			}

			// ...then try smoothing the pos of the Vertex:
			nOld = new Vertex(n1a.x, n1a.y);
			smoothed = n1a.laplacianSmoothExclude(n2a);
			if (!n1a.equals(smoothed)) {
				n1a.setXY(smoothed.x, smoothed.y);
				inversionCheckAndRepair(n1a, nOld);
				n1a.update();
			}

			if (e4a.sumAngle(qa, n1a, e2b) >= Math.PI) { // Still angle >= 180 degrees?
				return null;
			}
		}

		if (e2a.sumAngle(qa, n2a, e4b) >= Math.PI) { // if angle >= 180 degrees...
			if (n2a.boundaryVertex()) { // exit if Vertex on boundary
				return null;
			}

			// ...then try smoothing the pos of the Vertex:
			nOld = new Vertex(n2a.x, n2a.y);
			smoothed = n2a.laplacianSmoothExclude(n1a);
			if (!n2a.equals(smoothed)) {
				n2a.setXY(smoothed.x, smoothed.y);
				inversionCheckAndRepair(n2a, nOld);
				n2a.update();
			}

			if (e2a.sumAngle(qa, n2a, e4b) >= Math.PI) { // Still angle >= 180 degrees?
				return null;
			}
		}
		// The new diagonal:
		eNew = new Edge(n3a, n3b);

		// Create the new quads:
		l = qa.neighborEdge(e4a.leftVertex, e4a);
		r = qa.neighborEdge(e4a.rightVertex, e4a);
		if (l == e1a) {
			l = e2b;
		} else {
			r = e2b;
		}
		q1 = new Quad(e4a, l, r, eNew);

		l = qb.neighborEdge(e4b.leftVertex, e4b);
		r = qb.neighborEdge(e4b.rightVertex, e4b);
		if (l == e1b) {
			l = e2a;
		} else {
			r = e2a;
		}
		q2 = new Quad(e4b, l, r, eNew);

		qa.disconnectEdges();
		qb.disconnectEdges();
		e1a.disconnectVertices();
		q1.connectEdges();
		q2.connectEdges();
		eNew.connectVertices();

		// Update lists:
		edgeList.set(edgeList.indexOf(e1a), eNew);

		elementList.set(qaIndex, q1);
		elementList.set(qbIndex, q2);

		d.elem = q1;
		d.e = eNew;
		if (n == n1a) {
			d.n = n3b;
		} else {
			d.n = n3a;
		}

		return d;
	}

	/**
	 * Create 2 new quads from 2 specified quads, in which the common edge of the
	 * given quads has been rotated one step in the CW direction. Delete the old
	 * quads. Update the vertices list.
	 * 
	 * @param qa  one of the two quads adjacent the edge to be switched, e1a
	 * @param e1a the edge to be switched
	 * @param n   one of e1a's vertices
	 * @return a dart representing the input dart after the operation is performed.
	 */
	private Dart switchDiagonalCW(Quad qa, Edge e1a, Vertex n) {
		Dart d = new Dart();
		Vertex n1a, n2a, n3a, n4a, n1b, n2b, n3b, n4b;
		Edge e2a, e3a, e4a, e1b, e2b, e3b, e4b;
		Edge eNew, l, r;
		Quad q1, q2;

		Quad qb = (Quad) qa.neighbor(e1a);
		int qaIndex = elementList.indexOf(qa), qbIndex = elementList.indexOf(qb);

		// First get the edges of qa in ccw order:
		n2a = qa.nextCCWVertex(e1a.leftVertex);
		if (n2a == e1a.rightVertex) {
			n1a = e1a.leftVertex;
		} else {
			n1a = e1a.rightVertex;
			n2a = e1a.leftVertex;
		}
		e2a = qa.neighborEdge(n2a, e1a);
		n3a = e2a.otherVertex(n2a);
		e3a = qa.neighborEdge(n3a, e2a);
		n4a = e3a.otherVertex(n3a);
		e4a = qa.neighborEdge(n4a, e3a);

		// Now get the edges of qb in ccw order:
		e1b = e1a;
		n2b = qb.nextCCWVertex(e1b.leftVertex);
		if (n2b == e1b.rightVertex) {
			n1b = e1b.leftVertex;
		} else {
			n1b = e1b.rightVertex;
			n2b = e1b.leftVertex;
		}
		e2b = qb.neighborEdge(n2b, e1b);
		n3b = e2b.otherVertex(n2b);
		e3b = qb.neighborEdge(n3b, e2b);
		n4b = e3b.otherVertex(n3b);
		e4b = qb.neighborEdge(n4b, e3b);

		Vertex nOld, smoothed;
		// Check to see if the switch will violate the mesh topology:
		if (e4a.sumAngle(qa, n1a, e2b) >= Math.PI) { // if angle >= 180 degrees...
			// ...then try smoothing the pos of the Vertex:
			nOld = new Vertex(n1a.x, n1a.y);
			smoothed = n1a.laplacianSmooth();
			if (!n1a.equals(smoothed)) {
				n1a.moveTo(smoothed);
				inversionCheckAndRepair(n1a, nOld);
				n1a.update();
			}

			if (e4a.sumAngle(qa, n1a, e2b) >= Math.PI) { // Still angle >= 180 degrees?
				return null;
			}
		}

		// Check to see if the switch will violate the mesh topology:
		if (e2a.sumAngle(qa, n2a, e4b) >= Math.PI) { // if angle >= 180 degrees...
			// ...then try smoothing the pos of the Vertex:
			nOld = new Vertex(n2a.x, n2a.y);
			smoothed = n2a.laplacianSmooth();
			if (!n2a.equals(smoothed)) {
				n2a.moveTo(smoothed);
				inversionCheckAndRepair(n2a, nOld);
				n2a.update();
			}

			if (e2a.sumAngle(qa, n2a, e4b) >= Math.PI) { // Still angle >= 180 degrees?
				return null;
			}
		}

		// The new diagonal:
		eNew = new Edge(n4a, n4b);

		// Create the new quads:
		l = qa.neighborEdge(e2a.leftVertex, e2a);
		r = qa.neighborEdge(e2a.rightVertex, e2a);
		if (l == e1a) {
			l = e4b;
		} else {
			r = e4b;
		}
		q1 = new Quad(e2a, l, r, eNew);

		l = qb.neighborEdge(e2b.leftVertex, e2b);
		r = qb.neighborEdge(e2b.rightVertex, e2b);
		if (l == e1b) {
			l = e4a;
		} else {
			r = e4a;
		}
		q2 = new Quad(e2b, l, r, eNew);

		qa.disconnectEdges();
		qb.disconnectEdges();
		e1a.disconnectVertices();
		q1.connectEdges();
		q2.connectEdges();
		eNew.connectVertices();

		// Update lists:
		edgeList.set(edgeList.indexOf(e1a), eNew);

		elementList.set(qaIndex, q1);
		elementList.set(qbIndex, q2);

		d.elem = q1;
		d.e = eNew;
		if (n == n1a) {
			d.n = n4a;
		} else {
			d.n = n4b;
		}

		return d;
	}

	private void globalSmooth() {
		Vertex n, nn, nOld;

		for (Object element : vertexList) {
			n = (Vertex) element;

			if (!n.boundaryVertex()) {

				// Try smoothing the pos of the Vertex:
				nOld = new Vertex(n.x, n.y);
				nn = n.laplacianSmooth();
				if (!n.equals(nn)) {
					n.setXY(nn.x, nn.y);
					inversionCheckAndRepair(n, nOld);
					n.update();
				}
			}
		}
	}

}

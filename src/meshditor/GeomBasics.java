package meshditor;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

/**
 * This is a basic geometry class with methods for reading and writing meshes,
 * sorting Vertex lists, printing lists, topology inspection, etc.
 */
public class GeomBasics extends Constants {

	public static List<Element> elementList;
	protected static List<Triangle> triangleList;
	protected static List<Vertex> vertexList;
	protected static List<Edge> edgeList;

	protected static Vertex leftmost = null, rightmost = null, uppermost = null, lowermost = null;

	public static boolean step = false;

	protected static TopoCleanup topoCleanup;
	protected static GlobalSmooth globalSmooth;

	static String meshFilename;
	static String meshDirectory = ".";
	static boolean meshLenOpt;
	static boolean meshAngOpt;

	public static void createNewLists() {
		elementList = new ArrayList<>();
		triangleList = new ArrayList<>();
		edgeList = new ArrayList<>();
		vertexList = new ArrayList<>();
	}

	public static void setParams(String filename, String dir, boolean len, boolean ang) {
		meshFilename = filename;
		meshDirectory = dir;
		meshLenOpt = len;
		meshAngOpt = ang;
	}

	/** Return the edgeList */
	public static List<Edge> getEdgeList() {
		return edgeList;
	}

	/** Return the vertexList */
	public static List<Vertex> getvertexList() {
		return vertexList;
	}

	/** Return the triangleList */
	public static List<Triangle> getTriangleList() {
		return triangleList;
	}

	/** Return the elementList */
	public static List<Element> getElementList() {
		return elementList;
	}

	private static GeomBasics curMethod = null;

	public static void setCurMethod(GeomBasics method) {
		curMethod = method;
	}

	public static GeomBasics getCurMethod() {
		return curMethod;
	}

	/** This method should be implemented in each of the subclasses. */
	public void step() {
	}

	/** Delete all the edges in the mesh. */
	public static void clearEdges() {
		Element curElem;
		for (Object element : elementList) {
			curElem = (Element) element;
			curElem.disconnectEdges();
		}

		for (Edge curEdge : edgeList) {
			curEdge.disconnectVertices();
		}

		Vertex curVertex;
		for (int i = 0; i < vertexList.size(); i++) {
			curVertex = vertexList.get(i);
			curVertex.edgeList.clear();
		}

		elementList.clear();
		if (triangleList != null) {
			triangleList.clear();
		}
		edgeList.clear();
	}

	/** Clear the vertexList, edgeList, triangleList and elementList. */
	public static void clearLists() {
		if (vertexList != null) {
			vertexList.clear();
		} else {
			vertexList = new ArrayList<>();
		}
		if (edgeList != null) {
			edgeList.clear();
		} else {
			edgeList = new ArrayList<Edge>();
		}
		if (triangleList != null) {
			triangleList.clear();
		} else {
			triangleList = new ArrayList<Triangle>();
		}
		if (elementList != null) {
			elementList.clear();
		} else {
			elementList = new ArrayList<Element>();
		}
	}

	/** Update distortion metric for all elements in mesh. */
	public static void updateMeshMetrics() {
		Triangle tri;
		Element elem;
		if (elementList.size() > 0) {
			for (Object element : elementList) {
				elem = (Element) element;
				if (elem != null) {
					elem.updateDistortionMetric();
				}
			}
		} else {
			for (Object element : triangleList) {
				tri = (Triangle) element;
				if (tri != null) {
					tri.updateDistortionMetric();
				}
			}
		}
	}

	/** @return a string containing the average and minimum element metrics. */
	public static String meshMetricsReport() {
		Triangle tri;
		Quad q;
		Element elem;
		double sumMetric = 0.0, minDM = java.lang.Double.MAX_VALUE;
		int size = 0, nTris = 0, nQuads = 0;
		String s;
		byte temp, no2valents = 0, no3valents = 0, no4valents = 0, no5valents = 0, no6valents = 0, noXvalents = 0;

		if (elementList.size() > 0) {

			for (Object element : elementList) {
				elem = (Element) element;

				if (elem != null) {
					size++;
					if (elem instanceof Triangle) {
						nTris++;
					} else if (elem instanceof Quad) {
						q = (Quad) elem;
						if (q.isFake) {
							nTris++;
						} else {
							nQuads++;
						}
					}
					elem.updateDistortionMetric();
					sumMetric += elem.distortionMetric;
					if (elem.distortionMetric < minDM) {
						minDM = elem.distortionMetric;
					}
				}
			}
		}

		if (triangleList.size() > 0) {
			for (Object element : triangleList) {
				tri = (Triangle) element;
				if (tri != null) {
					size++;
					nTris++;
					tri.updateDistortionMetric();
					sumMetric += tri.distortionMetric;
					if (tri.distortionMetric < minDM) {
						minDM = tri.distortionMetric;
					}
				}
			}

		}

		s = "Average distortion metric: " + (sumMetric / size) + "\n" + "Minimum distortion metric: " + minDM + "\n";

		for (Vertex n : vertexList) {

			temp = n.valence();
			if (temp == 2) {
				no2valents++;
			} else if (temp == 3) {
				no3valents++;
			} else if (temp == 4) {
				no4valents++;
			} else if (temp == 5) {
				no5valents++;
			} else if (temp == 6) {
				no6valents++;
			} else if (temp > 6) {
				noXvalents++;
			}
		}

		if (no2valents > 0) {
			s = s + "Number of 2-valent vertices: " + no2valents + "\n";
		}
		if (no3valents > 0) {
			s = s + "Number of 3-valent vertices: " + no3valents + "\n";
		}
		if (no4valents > 0) {
			s = s + "Number of 4-valent vertices: " + no4valents + "\n";
		}
		if (no5valents > 0) {
			s = s + "Number of 5-valent vertices: " + no5valents + "\n";
		}
		if (no6valents > 0) {
			s = s + "Number of 6-valent vertices: " + no6valents + "\n";
		}
		if (noXvalents > 0) {
			s = s + "Number of vertices with valence > 6: " + noXvalents + "\n";
		}

		s = s + "Number of quadrilateral elements: " + nQuads + "\n" + "Number of triangular elements: " + nTris + "\n"
				+ "Number of edges: " + edgeList.size() + "\n" + "Number of vertices: " + vertexList.size();

		return s;
	}

	/** Find inverted elements and paint them with red colour. */
	public static void detectInvertedElements() {
		Element elem;
		Triangle t;
		int i;
		for (i = 0; i < elementList.size(); i++) {
			elem = (Element) elementList.get(i);
			if (elem != null && elem.inverted()) {
				elem.markEdgesIllegal();
				Msg.warning("Element " + elem.descr() + " is inverted.");
				Msg.warning("It has firstVertex " + elem.firstVertex.descr());
			}
		}

		for (i = 0; i < triangleList.size(); i++) {
			t = (Triangle) triangleList.get(i);
			if (t != null && t.inverted()) {
				t.markEdgesIllegal();
				Msg.warning("Triangle " + t.descr() + " is inverted.");
				Msg.warning("It has firstVertex " + t.firstVertex.descr());
			}
		}
	}

	/** Output nr of tris and fake quads in mesh. */
	public static void countTriangles() {
		int fakes = 0, tris = 0;
		Triangle t;
		Element elem;
		for (Object element : elementList) {
			elem = (Element) element;
			if (elem instanceof Quad && ((Quad) elem).isFake) {
				fakes++;
			} else if (elem instanceof Triangle) {
				tris++;
			}
		}

		Msg.debug("Counted # of fake quads: " + fakes);
		Msg.debug("Counted # of triangles: " + tris);
	}

	/** Output warnings if mesh is not consistent. */
	public static void consistencyCheck() {
		Msg.debug("Entering consistencyCheck()");
		Vertex n;
		for (int i = 0; i < vertexList.size(); i++) {
			n = vertexList.get(i);

			Edge e;
			if (n.edgeList.size() == 0) {
				Msg.warning("edgeList.size()== 0 for Vertex " + n.descr());
			}

			for (int j = 0; j < n.edgeList.size(); j++) {
				e = n.edgeList.get(j);
				if (e == null) {
					Msg.warning("Vertex " + n.descr() + " has a null in its edgeList.");
				} else if (edgeList.indexOf(e) == -1) {
					Msg.warning("Edge " + e.descr() + " found in the edgeList of Vertex " + n.descr() + ", but not in global edgeList");
				}
			}
		}

		Edge e;
		for (int i = 0; i < edgeList.size(); i++) {
			e = edgeList.get(i);
			if (e.leftVertex.edgeList.indexOf(e) == -1) {
				Msg.warning("leftVertex of edge " + e.descr() + " has not got that edge in its .edgeList");
			}
			if (e.rightVertex.edgeList.indexOf(e) == -1) {
				Msg.warning("rightVertex of edge " + e.descr() + " has not got that edge in its .edgeList");
			}

			if (e.element1 == null && e.element2 == null) {
				Msg.warning("Edge " + e.descr() + " has null in both element pointers");
			}

			if (e.element1 == null && e.element2 != null) {
				Msg.warning("Edge " + e.descr() + " has null in element1 pointer");
			}

			if (e.element1 != null && !triangleList.contains(e.element1) && !elementList.contains(e.element1)) {
				Msg.warning("element1 of edge " + e.descr() + " is not found in triangleList or elementList");
			}

			if (e.element2 != null && !triangleList.contains(e.element2) && !elementList.contains(e.element2)) {
				Msg.warning("element2 of edge " + e.descr() + " is not found in triangleList or elementList");
			}

			if (vertexList.indexOf(e.leftVertex) == -1) {
				Msg.warning("leftVertex of edge " + e.descr() + " not found in vertexList.");
			}
			if (vertexList.indexOf(e.rightVertex) == -1) {
				Msg.warning("rightVertex of edge " + e.descr() + " not found in vertexList.");
			}
		}

		Triangle t;
		Vertex na, nb, nc;
		double cross1;

		for (Object element : triangleList) {
			t = (Triangle) element;
			if (!t.edgeList[0].hasElement(t)) {
				Msg.warning("edgeList[0] of triangle " + t.descr() + " has not got that triangle as an adjacent element");
			}

			if (!t.edgeList[1].hasElement(t)) {
				Msg.warning("edgeList[1] of triangle " + t.descr() + " has not got that triangle as an adjacent element");
			}

			if (!t.edgeList[2].hasElement(t)) {
				Msg.warning("edgeList[2] of triangle " + t.descr() + " has not got that triangle as an adjacent element");
			}

			if (t.edgeList[0].commonVertex(t.edgeList[1]) == null) {
				Msg.warning("edgeList[0] and edgeList[1] of triangle " + t.descr() + " has no common Vertex");
			}
			if (t.edgeList[1].commonVertex(t.edgeList[2]) == null) {
				Msg.warning("edgeList[1] and edgeList[2] of triangle " + t.descr() + " has no common Vertex");
			}
			if (t.edgeList[2].commonVertex(t.edgeList[0]) == null) {
				Msg.warning("edgeList[2] and edgeList[0] of triangle " + t.descr() + " has no common Vertex");
			}

			na = t.edgeList[0].leftVertex;
			nb = t.edgeList[0].rightVertex;
			nc = t.oppositeOfEdge(t.edgeList[0]);

			cross1 = cross(na, nc, nb, nc); // The cross product nanc x nbnc

			if (cross1 == 0 /* !t.areaLargerThan0() */) {
				Msg.warning("Degenerate triangle in triangleList, t= " + t.descr());
			}
		}

		Element elem;
		Quad q;
		for (Object element : elementList) {
			elem = (Element) element;
			if (elem == null) {
				Msg.debug("elementList has a null-entry.");
			} else if (elem instanceof Quad) {
				q = (Quad) elem;

				if (!q.edgeList[base].hasElement(q)) {
					Msg.warning("edgeList[base] of quad " + q.descr() + " has not got that quad as an adjacent element");
				}

				if (!q.edgeList[left].hasElement(q)) {
					Msg.warning("edgeList[left] of quad " + q.descr() + " has not got that quad as an adjacent element");
				}

				if (!q.edgeList[right].hasElement(q)) {
					Msg.warning("edgeList[right] of quad " + q.descr() + " has not got that quad as an adjacent element");
				}

				if (!q.isFake && !q.edgeList[top].hasElement(q)) {
					Msg.warning("edgeList[top] of quad " + q.descr() + " has not got that quad as an adjacent element");
				}

				if (q.edgeList[base].commonVertex(q.edgeList[left]) == null) {
					Msg.warning("edgeList[base] and edgeList[left] of quad " + q.descr() + " has no common Vertex");
				}
				if (q.edgeList[base].commonVertex(q.edgeList[right]) == null) {
					Msg.warning("edgeList[base] and edgeList[right] of quad " + q.descr() + " has no common Vertex");
				}
				if (!q.isFake && q.edgeList[left].commonVertex(q.edgeList[top]) == null) {
					Msg.warning("edgeList[left] and edgeList[top] of quad " + q.descr() + " has no common Vertex");
				}
				if (!q.isFake && q.edgeList[right].commonVertex(q.edgeList[top]) == null) {
					Msg.warning("edgeList[right] and edgeList[top] of quad " + q.descr() + " has no common Vertex");
				}

				if (q.isFake && q.edgeList[left].commonVertex(q.edgeList[right]) == null) {
					Msg.warning("edgeList[left] and edgeList[right] of fake quad " + q.descr() + " has no common Vertex");
				}
			}
		}

		Msg.debug("Leaving consistencyCheck()");

	}

	/** Load a mesh from a file. */
	public static List<Element> loadMesh() {
		FileInputStream fis;
		Vertex Vertex1, Vertex2, Vertex3, Vertex4;
		Edge edge1, edge2, edge3, edge4;
		Triangle t;
		Quad q;

		elementList = new ArrayList<>();
		triangleList = new ArrayList<>();
		edgeList = new ArrayList<>();
		ArrayList<Vertex> usvertexList = new ArrayList<>();
		ArrayList<double[]> triangles = new ArrayList<>();

		try {
			String file = "C:\\Users\\micyc\\Documents\\Repositories\\Q-Morph\\examples\\thesis-tri\\difficult.mesh";
//			String file = meshDirectory + meshFilename;
			fis = new FileInputStream(file);
			BufferedReader in = new BufferedReader(new InputStreamReader(fis));
			double x1, x2, x3, x4, y1, y2, y3, y4;
			int i = 0;

			try {
				String inputLine;
				inputLine = in.readLine();
				while (inputLine != null) {
					cInd = 0;
					x1 = nextDouble(inputLine);
					y1 = nextDouble(inputLine);
					x2 = nextDouble(inputLine);
					y2 = nextDouble(inputLine);
					x3 = nextDouble(inputLine);
					y3 = nextDouble(inputLine);
					x4 = nextDouble(inputLine); // quad possible
					y4 = nextDouble(inputLine); // quad possible

					double[] triangle = new double[] { x1, y1, x2, y2, x3, y3 };
					triangles.add(triangle);

					Vertex1 = new Vertex(x1, y1);
					if (!usvertexList.contains(Vertex1)) {
						usvertexList.add(Vertex1);
					} else {
						Vertex1 = usvertexList.get(usvertexList.indexOf(Vertex1));
					}
					Vertex2 = new Vertex(x2, y2);
					if (!usvertexList.contains(Vertex2)) {
						usvertexList.add(Vertex2);
					} else {
						Vertex2 = usvertexList.get(usvertexList.indexOf(Vertex2));
					}
					Vertex3 = new Vertex(x3, y3);
					if (!usvertexList.contains(Vertex3)) {
						usvertexList.add(Vertex3);
					} else {
						Vertex3 = usvertexList.get(usvertexList.indexOf(Vertex3));
					}

					edge1 = new Edge(Vertex1, Vertex2);
					if (!edgeList.contains(edge1)) {
						edgeList.add(edge1);
						edge1.connectVertices();
					} else {
						edge1 = edgeList.get(edgeList.indexOf(edge1));
					}

					edge2 = new Edge(Vertex1, Vertex3);
					if (!edgeList.contains(edge2)) {
						edgeList.add(edge2);
						edge2.connectVertices();
					} else {
						edge2 = edgeList.get(edgeList.indexOf(edge2));
					}

					if (!Double.isNaN(x4) && !Double.isNaN(y4)) {
						Vertex4 = new Vertex(x4, y4);
						if (!usvertexList.contains(Vertex4)) {
							usvertexList.add(Vertex4);
						} else {
							Vertex4 = usvertexList.get(usvertexList.indexOf(Vertex4));
						}

						edge3 = new Edge(Vertex2, Vertex4);
						if (!edgeList.contains(edge3)) {
							edgeList.add(edge3);
							edge3.connectVertices();
						} else {
							edge3 = edgeList.get(edgeList.indexOf(edge3));
						}

						edge4 = new Edge(Vertex3, Vertex4);
						if (!edgeList.contains(edge4)) {
							edgeList.add(edge4);
							edge4.connectVertices();
						} else {
							edge4 = edgeList.get(edgeList.indexOf(edge4));
						}

						q = new Quad(edge1, edge2, edge3, edge4);
						q.connectEdges();
						elementList.add(q);
					} else {
						edge3 = new Edge(Vertex2, Vertex3);
						if (!edgeList.contains(edge3)) {
							edgeList.add(edge3);
							edge3.connectVertices();
						} else {
							edge3 = edgeList.get(edgeList.indexOf(edge3));
						}

						t = new Triangle(edge1, edge2, edge3);
						t.connectEdges();
						triangleList.add(t);
						// elementList.add(t);
					}
					inputLine = in.readLine();
				}
			} catch (Exception e) {
				Msg.error("Cannot read triangle-mesh data.");
			}
		} catch (Exception e) {
//			Msg.error("File " + meshFilename + " not found.");
			e.printStackTrace();
		}

		vertexList = usvertexList; // sortVertices(usvertexList);
		double[][] ts = triangles.toArray(new double[triangles.size()][6]);
		loadTriangleMesh(ts);
		return elementList;
	}

	/** Load a triangle mesh from a file. */
	public static List<Triangle> loadTriangleMesh() {
		FileInputStream fis;
		Vertex Vertex1, Vertex2, Vertex3;
		Edge edge1, edge2, edge3;
		Triangle t;

		triangleList = new ArrayList<Triangle>();
		edgeList = new ArrayList<Edge>();
		ArrayList<Vertex> usvertexList = new ArrayList<Vertex>();
		
		ArrayList<double[]> triangles = new ArrayList<>();

		try {
			String file = "C:\\Users\\micyc\\Documents\\Repositories\\Q-Morph\\examples\\thesis-tri\\difficult.mesh";
//			String file = meshDirectory + meshFilename;
			fis = new FileInputStream(file);
			BufferedReader in = new BufferedReader(new InputStreamReader(fis));
			double x1, x2, x3, y1, y2, y3, len1 = 0, len2 = 0, len3 = 0, ang1 = 0, ang2 = 0, ang3 = 0;
			int i = 0;

			try {
				String inputLine;
				inputLine = in.readLine();
				while (inputLine != null) {
					cInd = 0;
					x1 = nextDouble(inputLine);
					y1 = nextDouble(inputLine);
					x2 = nextDouble(inputLine);
					y2 = nextDouble(inputLine);
					x3 = nextDouble(inputLine);
					y3 = nextDouble(inputLine);
					
					double[] triangle = new double[] { x1, y1, x2, y2, x3, y3 };
					triangles.add(triangle);

//					Vertex1 = new Vertex(x1, y1);
//					if (!usvertexList.contains(Vertex1)) {
//						usvertexList.add(Vertex1);
//					} else {
//						Vertex1 = usvertexList.get(usvertexList.indexOf(Vertex1));
//					}
//					Vertex2 = new Vertex(x2, y2);
//					if (!usvertexList.contains(Vertex2)) {
//						usvertexList.add(Vertex2);
//					} else {
//						Vertex2 = usvertexList.get(usvertexList.indexOf(Vertex2));
//					}
//					Vertex3 = new Vertex(x3, y3);
//					if (!usvertexList.contains(Vertex3)) {
//						usvertexList.add(Vertex3);
//					} else {
//						Vertex3 = usvertexList.get(usvertexList.indexOf(Vertex3));
//					}
//
//					edge1 = new Edge(Vertex1, Vertex2);
//					if (!edgeList.contains(edge1)) {
//						edgeList.add(edge1);
//					} else {
//						edge1 = edgeList.get(edgeList.indexOf(edge1));
//					}
//					edge1.leftVertex.connectToEdge(edge1);
//					edge1.rightVertex.connectToEdge(edge1);
//
//					edge2 = new Edge(Vertex2, Vertex3);
//					if (!edgeList.contains(edge2)) {
//						edgeList.add(edge2);
//					} else {
//						edge2 = edgeList.get(edgeList.indexOf(edge2));
//					}
//					edge2.leftVertex.connectToEdge(edge2);
//					edge2.rightVertex.connectToEdge(edge2);
//
//					edge3 = new Edge(Vertex1, Vertex3);
//					if (!edgeList.contains(edge3)) {
//						edgeList.add(edge3);
//					} else {
//						edge3 = edgeList.get(edgeList.indexOf(edge3));
//					}
//					edge3.leftVertex.connectToEdge(edge3);
//					edge3.rightVertex.connectToEdge(edge3);
//
//					if (meshLenOpt) {
//						len1 = nextDouble(inputLine);
//						len2 = nextDouble(inputLine);
//						len3 = nextDouble(inputLine);
//					}
//
//					if (meshAngOpt) {
//						ang1 = nextDouble(inputLine);
//						ang2 = nextDouble(inputLine);
//						ang3 = nextDouble(inputLine);
//					}
//
//					t = new Triangle(edge1, edge2, edge3, len1, len2, len3, ang1, ang2, ang3, meshLenOpt, meshAngOpt);
//					t.connectEdges();
//					triangleList.add(t);
					inputLine = in.readLine();
				}
			} catch (Exception e) {
				Msg.error("Cannot read triangle-mesh data.");
			}
		} catch (Exception e) {
			Msg.error("File " + meshFilename + " not found.");
		}
		vertexList = usvertexList; // sortVertices(usvertexList);
		return loadTriangleMesh(triangles.toArray(new double[triangles.size()][6]));
	}

	/**
		 * 
		 * @param triangles [t1=[x1,y1,x2,y2,x3,y3], t2=[x1,y1,x2,y2,x3,y3]...]
		 * @return
		 */
		public static List<Triangle> loadTriangleMesh(double[][] triangles) {
	
			Vertex Vertex1, Vertex2, Vertex3;
			Edge edge1, edge2, edge3;
			Triangle t;
	
			elementList = new ArrayList<>();
			triangleList = new ArrayList<>();
			edgeList = new ArrayList<>();
			Map<Vertex, Vertex> vertexSet = new HashMap<>();
			Map<Edge, Edge> edgeSet = new HashMap<>();
	
			for (double[] triangle : triangles) {
				double len1 = 0, len2 = 0, len3 = 0, ang1 = 0, ang2 = 0, ang3 = 0;
				double x1 = triangle[0];
				double y1 = triangle[1];
				double x2 = triangle[2];
				double y2 = triangle[3];
				double x3 = triangle[4];
				double y3 = triangle[5];
	
				final Vertex Vertex1F = new Vertex(x1, y1);
				Vertex1 = vertexSet.computeIfAbsent(Vertex1F, v -> Vertex1F);
				final Vertex Vertex2F = new Vertex(x2, y2);
				Vertex2 = vertexSet.computeIfAbsent(Vertex2F, v -> Vertex2F);
				final Vertex Vertex3F = new Vertex(x3, y3);
				Vertex3 = vertexSet.computeIfAbsent(Vertex3F, v -> Vertex3F);
	
				final Edge edge1F = new Edge(Vertex1, Vertex2);
				edge1 = edgeSet.computeIfAbsent(edge1F, e -> edge1F);
				edge1.leftVertex.connectToEdge(edge1);
				edge1.rightVertex.connectToEdge(edge1);
				
				final Edge edge2F = new Edge(Vertex2, Vertex3);
				edge2 = edgeSet.computeIfAbsent(edge2F, e -> edge2F);
				edge2.leftVertex.connectToEdge(edge2);
				edge2.rightVertex.connectToEdge(edge2);
				
				final Edge edge3F = new Edge(Vertex1, Vertex3);
				edge3 = edgeSet.computeIfAbsent(edge3F, e -> edge3F);
				edge3.leftVertex.connectToEdge(edge3);
				edge3.rightVertex.connectToEdge(edge3);
	
				if (meshLenOpt) {
				}
	
				if (meshAngOpt) {
				}
	
				t = new Triangle(edge1, edge2, edge3, len1, len2, len3, ang1, ang2, ang3, meshLenOpt, meshAngOpt);
	//			t = new Triangle(triangle);
				t.connectEdges();
				triangleList.add(t);
	
			}
			vertexList = new ArrayList<>(vertexSet.keySet());
			edgeList = new ArrayList<>(edgeSet.keySet());
	
			return triangleList;
		}

	/** A method to read Vertex files. */
	public static ArrayList<Vertex> loadVertices() {
		FileInputStream fis;
		Vertex Vertex1, Vertex2, Vertex3, Vertex4;
		ArrayList<Vertex> usvertexList = new ArrayList<Vertex>();

		try {
			fis = new FileInputStream(meshDirectory + meshFilename);
			BufferedReader in = new BufferedReader(new InputStreamReader(fis));
			double x1, x2, x3, y1, y2, y3, x4, y4;

			try {
				String inputLine;
				inputLine = in.readLine();
				while (inputLine != null) {
					cInd = 0;
					x1 = nextDouble(inputLine);
					y1 = nextDouble(inputLine);
					x2 = nextDouble(inputLine);
					y2 = nextDouble(inputLine);
					x3 = nextDouble(inputLine);
					y3 = nextDouble(inputLine);
					x4 = nextDouble(inputLine);
					y4 = nextDouble(inputLine);

					if (!Double.isNaN(x1) && !Double.isNaN(y1)) {
						Vertex1 = new Vertex(x1, y1);
						if (!usvertexList.contains(Vertex1)) {
							usvertexList.add(Vertex1);
						}
					}
					if (!Double.isNaN(x2) && !Double.isNaN(y2)) {
						Vertex2 = new Vertex(x2, y2);
						if (!usvertexList.contains(Vertex2)) {
							usvertexList.add(Vertex2);
						}
					}
					if (!Double.isNaN(x3) && !Double.isNaN(y3)) {
						Vertex3 = new Vertex(x3, y3);
						if (!usvertexList.contains(Vertex3)) {
							usvertexList.add(Vertex3);
						}
					}
					if (!Double.isNaN(x4) && !Double.isNaN(y4)) {
						Vertex4 = new Vertex(x4, y4);
						if (!usvertexList.contains(Vertex4)) {
							usvertexList.add(Vertex4);
						}
					}
					inputLine = in.readLine();
				}
			} catch (Exception e) {
				Msg.error("Cannot read Vertex file data.");
			}
		} catch (Exception e) {
			Msg.error("File " + meshFilename + " not found.");
		}

		// vertexList= sortVertices(usvertexList);
		vertexList = usvertexList;
		return usvertexList;
	}

	/**
	 * Method for writing to a LaTeX drawing format (need the epic and eepic
	 * packages).
	 */
	public static boolean exportMeshToLaTeX(String filename, int unitlength, double xcorr, double ycorr, boolean visibleVertices) {
		FileOutputStream fos;
		Edge edge;
		Vertex n;
		int i;
		ArrayList<Edge> boundary = new ArrayList<Edge>();

		findExtremeVertices();

		// Collect boundary edges in a list
		for (i = 0; i < edgeList.size(); i++) {
			edge = edgeList.get(i);
			if (edge.boundaryEdge()) {
				boundary.add(edge);
			}
		}

		try {
			fos = new FileOutputStream(filename);
			BufferedWriter out = new BufferedWriter(new OutputStreamWriter(fos));
			double x1, x2, x3, x4, y1, y2, y3, y4;
			double width = rightmost.x - leftmost.x, height = uppermost.y - lowermost.y;

			try {

				out.write("% Include in the header of your file:");
				out.newLine();
				out.write("% \\usepackage{epic, eepic}");
				out.newLine();
				out.newLine();
				out.write("\\begin{figure}[!Htbp]");
				out.newLine();
				out.write("\\begin{center}");
				out.newLine();
				out.write("\\setlength{\\unitlength}{" + unitlength + "mm}");
				out.newLine();
				out.write("\\begin{picture}(" + width + "," + height + ")");
				out.newLine();
				out.write("\\filltype{black}");
				out.newLine();

				// All boundary edges...
				out.write("\\thicklines");
				out.newLine();
				for (i = 0; i < boundary.size(); i++) {
					edge = boundary.get(i);

					x1 = edge.leftVertex.x + xcorr;
					y1 = edge.leftVertex.y + ycorr;
					x2 = edge.rightVertex.x + xcorr;
					y2 = edge.rightVertex.y + ycorr;

					out.write("\\drawline[1](" + x1 + "," + y1 + ")(" + x2 + "," + y2 + ")");
					out.newLine();
				}

				// All other edges...
				out.write("\\thinlines");
				out.newLine();
				for (i = 0; i < edgeList.size(); i++) {
					edge = edgeList.get(i);

					if (!edge.boundaryEdge()) {
						x1 = edge.leftVertex.x + xcorr;
						y1 = edge.leftVertex.y + ycorr;
						x2 = edge.rightVertex.x + xcorr;
						y2 = edge.rightVertex.y + ycorr;

						out.write("\\drawline[1](" + x1 + "," + y1 + ")(" + x2 + "," + y2 + ")");
						out.newLine();
					}
				}

				// All vertices...
				if (visibleVertices) {
					for (i = 0; i < vertexList.size(); i++) {
						n = vertexList.get(i);
						out.write("\\put(" + (n.x + xcorr) + "," + (n.y + ycorr) + "){\\circle*{0.1}}");
						out.newLine();
					}
				}

				out.write("\\end{picture}");
				out.newLine();
				out.write("\\end{center}");
				out.newLine();
				out.write("\\end{figure}");
				out.newLine();

				out.close();
			} catch (Exception e) {
				Msg.error("Cannot write quad-mesh data export file.");
			}
		} catch (Exception e) {
			Msg.error("File " + filename + " not found.");
		}
		return true;
	}

	/** Write all elements in elementList to a file. */
	public static boolean writeQuadMesh(String filename, List<? extends Element> elementList) {
		FileOutputStream fos;
		Triangle t;
		Quad q;

		try {
			fos = new FileOutputStream(filename);
			BufferedWriter out = new BufferedWriter(new OutputStreamWriter(fos));
			double x1, x2, x3, x4, y1, y2, y3, y4;

			try {
				for (Element elem : elementList) {
					if (elem instanceof Quad) {
						q = (Quad) elem;
						x1 = q.edgeList[base].leftVertex.x;
						y1 = q.edgeList[base].leftVertex.y;
						x2 = q.edgeList[base].rightVertex.x;
						y2 = q.edgeList[base].rightVertex.y;
						x3 = q.edgeList[left].otherVertex(q.edgeList[base].leftVertex).x;
						y3 = q.edgeList[left].otherVertex(q.edgeList[base].leftVertex).y;
						x4 = q.edgeList[right].otherVertex(q.edgeList[base].rightVertex).x;
						y4 = q.edgeList[right].otherVertex(q.edgeList[base].rightVertex).y;

						out.write(x1 + ", " + y1 + ", " + x2 + ", " + y2 + ", " + x3 + ", " + y3 + ", " + x4 + ", " + y4);
					} else {
						t = (Triangle) elem;
						x1 = t.edgeList[0].leftVertex.x;
						y1 = t.edgeList[0].leftVertex.y;
						x2 = t.edgeList[0].rightVertex.x;
						y2 = t.edgeList[0].rightVertex.y;
						if (!t.edgeList[1].leftVertex.equals(t.edgeList[0].leftVertex)
								&& !t.edgeList[1].leftVertex.equals(t.edgeList[0].rightVertex)) {
							x3 = t.edgeList[1].leftVertex.x;
							y3 = t.edgeList[1].leftVertex.y;
						} else {
							x3 = t.edgeList[1].rightVertex.x;
							y3 = t.edgeList[1].rightVertex.y;
						}
						out.write(x1 + ", " + y1 + ", " + x2 + ", " + y2 + ", " + x3 + ", " + y3);
					}
					out.newLine();
				}
				out.close();
			} catch (Exception e) {
				Msg.error("Cannot write quad-mesh data.");
			}
		} catch (Exception e) {
			Msg.error("File " + filename + " not found.");
		}
		return true;
	}

	/** Write all elements in elementList and triangleList to a file. */
	public static boolean writeMesh(String filename) {
		FileOutputStream fos = null;
		Element elem;
		Triangle t;
		Quad q;

		try {
			fos = new FileOutputStream(filename);
		} catch (Exception e) {
			Msg.error("File " + filename + " not found.");
		}

		BufferedWriter out = new BufferedWriter(new OutputStreamWriter(fos));
		double x1, x2, x3, x4, y1, y2, y3, y4;

		if (triangleList != null) {
			for (Object element : triangleList) {
				t = (Triangle) element;
				x1 = t.edgeList[0].leftVertex.x;
				y1 = t.edgeList[0].leftVertex.y;
				x2 = t.edgeList[0].rightVertex.x;
				y2 = t.edgeList[0].rightVertex.y;
				if (!t.edgeList[1].leftVertex.equals(t.edgeList[0].leftVertex)
						&& !t.edgeList[1].leftVertex.equals(t.edgeList[0].rightVertex)) {
					x3 = t.edgeList[1].leftVertex.x;
					y3 = t.edgeList[1].leftVertex.y;
				} else {
					x3 = t.edgeList[1].rightVertex.x;
					y3 = t.edgeList[1].rightVertex.y;
				}
				try {
					out.write(x1 + ", " + y1 + ", " + x2 + ", " + y2 + ", " + x3 + ", " + y3);
					out.newLine();
				} catch (Exception e) {
					Msg.error("Cannot write quad-mesh data.");
				}

			}
		}

		if (elementList != null) {
			for (Object element : elementList) {

				if (element instanceof Quad) {
					q = (Quad) element;

					x1 = q.edgeList[base].leftVertex.x;
					y1 = q.edgeList[base].leftVertex.y;
					x2 = q.edgeList[base].rightVertex.x;
					y2 = q.edgeList[base].rightVertex.y;
					x3 = q.edgeList[left].otherVertex(q.edgeList[base].leftVertex).x;
					y3 = q.edgeList[left].otherVertex(q.edgeList[base].leftVertex).y;
					if (!q.isFake) {
						x4 = q.edgeList[right].otherVertex(q.edgeList[base].rightVertex).x;
						y4 = q.edgeList[right].otherVertex(q.edgeList[base].rightVertex).y;
						try {
							out.write(x1 + ", " + y1 + ", " + x2 + ", " + y2 + ", " + x3 + ", " + y3 + ", " + x4 + ", " + y4);
							out.newLine();
						} catch (Exception e) {
							Msg.error("Cannot write quad-mesh data.");
						}

					} else {
						try {
							out.write(x1 + ", " + y1 + ", " + x2 + ", " + y2 + ", " + x3 + ", " + y3);
						} catch (Exception e) {
							Msg.error("Cannot write quad-mesh data.");
						}
					}
				} else if (element instanceof Triangle) {
					t = (Triangle) element;
					x1 = t.edgeList[0].leftVertex.x;
					y1 = t.edgeList[0].leftVertex.y;
					x2 = t.edgeList[0].rightVertex.x;
					y2 = t.edgeList[0].rightVertex.y;
					if (!t.edgeList[1].leftVertex.equals(t.edgeList[0].leftVertex)
							&& !t.edgeList[1].leftVertex.equals(t.edgeList[0].rightVertex)) {
						x3 = t.edgeList[1].leftVertex.x;
						y3 = t.edgeList[1].leftVertex.y;
					} else {
						x3 = t.edgeList[1].rightVertex.x;
						y3 = t.edgeList[1].rightVertex.y;
					}
					try {
						out.write(x1 + ", " + y1 + ", " + x2 + ", " + y2 + ", " + x3 + ", " + y3);
						out.newLine();
					} catch (Exception e) {
						Msg.error("Cannot write quad-mesh data.");
					}
				}
			}
		}

		try {
			out.close();
		} catch (Exception e) {
			Msg.error("Cannot write quad-mesh data.");
		}
		return true;
	}

	/** Write all vertices in vertexList to a file. */
	public static boolean writeVertices(String filename) {
		FileOutputStream fos;

		try {
			fos = new FileOutputStream(filename);
			BufferedWriter out = new BufferedWriter(new OutputStreamWriter(fos));
			double x, y;

			try {
				if (vertexList != null) {
					for (Vertex n : vertexList) {
						x = n.x;
						y = n.y;
						out.write(x + ", " + y);
						out.newLine();
					}
				}
				out.close();
			} catch (Exception e) {
				Msg.error("Cannot write Vertex data.");
			}
		} catch (Exception e) {
			Msg.error("Could not open file " + filename);
		}
		return true;
	}

	/** Find the leftmost, rightmost, uppermost, and lowermost vertices. */
	public static void findExtremeVertices() {
		// vertexList= sortVertices(vertexList);
		if (vertexList == null || vertexList.size() == 0) {
			leftmost = null;
			rightmost = null;
			uppermost = null;
			lowermost = null;
			return;
		}

		leftmost = vertexList.get(0);
		rightmost = leftmost;
		uppermost = leftmost;
		lowermost = leftmost;

		Vertex curVertex;
		for (int i = 1; i < vertexList.size(); i++) {
			curVertex = vertexList.get(i);

			if ((curVertex.x < leftmost.x) || (curVertex.x == leftmost.x && curVertex.y > leftmost.y)) {
				leftmost = curVertex;
			}
			if ((curVertex.x > rightmost.x) || (curVertex.x == rightmost.x && curVertex.y < rightmost.y)) {
				rightmost = curVertex;
			}

			if ((curVertex.y > uppermost.y) || (curVertex.y == uppermost.y && curVertex.x < uppermost.x)) {
				uppermost = curVertex;
			}
			if ((curVertex.y < lowermost.y) || (curVertex.y == lowermost.y && curVertex.x < lowermost.x)) {
				lowermost = curVertex;
			}
		}

	}

	/** Sort vertices left to right. Higher y-values are preferred to lower ones. */
	public static List<Vertex> sortVertices(List<Vertex> unsortedVertices) {
		List<Vertex> sortedVertices = new ArrayList<>();
		Vertex curVertex, candVertex;
		while (unsortedVertices.size() > 0) {
			curVertex = unsortedVertices.get(0);
			for (int i = 1; i < unsortedVertices.size(); i++) {
				candVertex = unsortedVertices.get(i);
				if (candVertex.x < curVertex.x || (candVertex.x == curVertex.x && candVertex.y < curVertex.y)) {
					curVertex = candVertex;
				}
			}
			sortedVertices.add(curVertex);
			unsortedVertices.remove(unsortedVertices.indexOf(curVertex));
		}

		// Find the leftmost, rightmost, uppermost, and lowermost vertices.
		leftmost = sortedVertices.get(0);
		rightmost = sortedVertices.get(sortedVertices.size() - 1);
		uppermost = leftmost;
		lowermost = leftmost;

		for (int i = 1; i < sortedVertices.size(); i++) {
			curVertex = sortedVertices.get(i);
			if (curVertex.y > uppermost.y) {
				uppermost = curVertex;
			}
			if (curVertex.y < lowermost.y) {
				lowermost = curVertex;
			}
		}

		return sortedVertices;
	}

	private static int cInd = 0;

	/**
	 * Method to assist the different load... methods.
	 * 
	 * @param iline a comma-separated string
	 * @return the next double value from iline. If no more numbers, then return
	 *         NaN.
	 */
	private static double nextDouble(String iline) {
		String ndouble;
		if (cInd > iline.length()) {
			return Double.NaN;
		}
		int nInd = iline.indexOf(",", cInd);
		if (nInd == -1) {
			nInd = iline.length();
		}

		ndouble = iline.substring(cInd, nInd);
		cInd = nInd + 1;
		return java.lang.Double.valueOf(ndouble).doubleValue();
	}

	public static void printVectors(List<?> vectorList) {
		if (Msg.debugMode) {
			MyVector v;
			for (Object element : vectorList) {
				v = (MyVector) element;
				v.printMe();
			}
		}
	}

	public static void printElements(List<?> elemList) {
		if (Msg.debugMode) {
			Element elem;
			for (Object element : elemList) {
				elem = (Element) element;
				elem.printMe();
			}
		}
	}

	public static void printTriangles(List<Triangle> triangleList) {
		Msg.debug("triangleList: (size== " + triangleList.size() + ")");
		printElements(triangleList);
	}

	public static void printQuads(List<?> quadList) {
		Msg.debug("quadList: (size== " + quadList.size() + ")");
		printElements(quadList);
	}

	public static void printEdgeList(List<Edge> edgeList) {
		if (Msg.debugMode) {
			for (Edge edge : edgeList) {
				edge.printMe();
			}
		}
	}

	public static void printVertices(List<Vertex> vertexList) {
		if (Msg.debugMode) {
			Msg.debug("vertexList:");
			for (Vertex vertex : vertexList) {
				vertex.printMe();
			}
		}
	}

	/** */
	public static void printValences() {
		for (Vertex n : vertexList) {
			Msg.debug("Vertex " + n.descr() + " has valence " + n.valence());
		}
	}

	/** */
	public static void printValPatterns() {
		Vertex[] neighbors;
		for (Vertex n : vertexList) {
			if (!n.boundaryVertex()) {
				neighbors = n.ccwSortedNeighbors();
				n.createValencePattern(neighbors);
				Msg.debug("Vertex " + n.descr() + " has valence pattern " + n.valDescr());
			}
		}
	}

	/** */
	public static void printAnglesAtSurrondingVertices() {
		Vertex[] neighbors;
		double[] angles;
		for (Vertex n : vertexList) {
			if (!n.boundaryVertex()) {
				neighbors = n.ccwSortedNeighbors();
				n.createValencePattern(neighbors);
				angles = n.surroundingAngles(neighbors, n.pattern[0] - 2);

				Msg.debug("Angles at the vertices surrounding Vertex " + n.descr() + ":");
				for (int j = 0; j < n.pattern[0] - 2; j++) {
					Msg.debug("angles[" + j + "]== " + Math.toDegrees(angles[j]) + " (in degrees)");
				}
			}
		}
	}

	/**
	 * Do inversion test and repair inversion if requiered
	 * 
	 * @return true if any repairing was neccessary, else return false.
	 */
	public static boolean inversionCheckAndRepair(Vertex newN, Vertex oldPos) {
		Msg.debug("Entering inversionCheckAndRepair(..), Vertex oldPos: " + oldPos.descr());
		List<Element> elements = newN.adjElements();
		if (newN.invertedOrZeroAreaElements(elements)) {
			if (!newN.incrAdjustUntilNotInvertedOrZeroArea(oldPos, elements)) {
				for (Element elem : elements) {
					if (elem.invertedOrZeroArea()) {
						break;
					}
				}
				Msg.error("It seems that an element was inverted initially.");
				return false;
			}
			Msg.debug("Leaving inversionCheckAndRepair(..)");
			return true;
		} else {
			Msg.debug("Leaving inversionCheckAndRepair(..)");
			return false;
		}
	}

	/**
	 * Quad q is to be collapsed. Vertices n1 and n2 are two opposite vertices in q.
	 * This method tries to find a location inside the current q to which n1 and n2
	 * can safely be relocated and joined without causing any adjacent elements to
	 * become inverted. The first candidate location is the centroid of the quad. If
	 * that location is not suitable, the method tries locations on the vectors from
	 * the centroid towards n1 and from the centroid towards n2. The first suitable
	 * location found is returned.
	 * 
	 * @param q  the quad to be collapsed
	 * @param n1 the Vertex in quad q that is to be joined with opposite Vertex n2
	 * @param n2 the Vertex in quad q that is to be joined with opposite Vertex n1
	 * @return a position inside quad q to which both n1 and n2 can be relocated
	 *         without inverting any of their adjacent elements.
	 */
	public static Vertex safeNewPosWhenCollapsingQuad(Quad q, Vertex n1, Vertex n2) {
		Msg.debug("Entering safeNewPosWhenCollapsingQuad(..)");

		Vertex n = q.centroid();
		MyVector back2n1 = new MyVector(n, n1), back2n2 = new MyVector(n, n2);
		double startX = n.x, startY = n.y;
		double xstepn1 = back2n1.x / 50.0, ystepn1 = back2n1.y / 50.0, xstepn2 = back2n2.x / 50.0, ystepn2 = back2n2.y / 50.0;
		double xincn1, yincn1, xincn2, yincn2;
		int steps2n1, steps2n2, i;
		List<Element> l1 = n1.adjElements(), l2 = n2.adjElements();

		if (!q.anyInvertedElementsWhenCollapsed(n, n1, n2, l1, l2)) {
			Msg.debug("Leaving safeNewPosWhenCollapsingQuad(..): found");
			return n;
		}

		// Calculate the parameters for direction n to n1
		if (Math.abs(xstepn1) < COINCTOL || Math.abs(ystepn1) < COINCTOL) {
			Msg.debug("...ok, resorting to use of minimum increment");
			if (Math.abs(back2n1.x) < Math.abs(back2n1.y)) {
				if (back2n1.x < 0) {
					xstepn1 = -COINCTOL;
				} else {
					xstepn1 = COINCTOL;
				}

				// abs(ystepn1/xstepn1) = abs(n1.y/n1.x)
				ystepn1 = Math.abs(n1.y) * COINCTOL / Math.abs(n1.x);
				if (back2n1.y < 0) {
					ystepn1 = -ystepn1;
				}

				steps2n1 = (int) (back2n1.x / xstepn1);
			} else {
				if (back2n1.y < 0) {
					ystepn1 = -COINCTOL;
				} else {
					ystepn1 = COINCTOL;
				}

				// abs(xstepn1/ystepn1) = abs(n1.x/n1.y)
				xstepn1 = Math.abs(n1.x) * COINCTOL / Math.abs(n1.y);
				if (back2n1.x < 0) {
					xstepn1 = -xstepn1;
				}

				steps2n1 = (int) (back2n1.y / ystepn1);
			}
		} else {
			xstepn1 = back2n1.x / 50.0;
			ystepn1 = back2n1.x / 50.0;
			steps2n1 = 50;
		}

		// Calculate the parameters for direction n to n2
		if (Math.abs(xstepn2) < COINCTOL || Math.abs(ystepn2) < COINCTOL) {
			Msg.debug("...ok, resorting to use of minimum increment");
			if (Math.abs(back2n2.x) < Math.abs(back2n2.y)) {
				if (back2n2.x < 0) {
					xstepn2 = -COINCTOL;
				} else {
					xstepn2 = COINCTOL;
				}

				// abs(ystepn2/xstepn2) = abs(n2.y/n2.x)
				ystepn2 = Math.abs(n2.y) * COINCTOL / Math.abs(n2.x);
				if (back2n2.y < 0) {
					ystepn2 = -ystepn2;
				}

				steps2n2 = (int) (back2n2.x / xstepn2);
			} else {
				if (back2n2.y < 0) {
					ystepn2 = -COINCTOL;
				} else {
					ystepn2 = COINCTOL;
				}

				// abs(xstepn2/ystepn2) = abs(n2.x/n2.y)
				xstepn2 = Math.abs(n2.x) * COINCTOL / Math.abs(n2.y);
				if (back2n2.x < 0) {
					xstepn2 = -xstepn2;
				}

				steps2n2 = (int) (back2n2.y / ystepn2);
			}
		} else {
			xstepn2 = back2n2.x / 50.0;
			ystepn2 = back2n2.x / 50.0;
			steps2n2 = 50;
		}

		Msg.debug("...back2n1.x is: " + back2n1.x);
		Msg.debug("...back2n1.y is: " + back2n1.y);
		Msg.debug("...xstepn1 is: " + xstepn1);
		Msg.debug("...ystepn1 is: " + ystepn1);

		Msg.debug("...back2n2.x is: " + back2n2.x);
		Msg.debug("...back2n2.y is: " + back2n2.y);
		Msg.debug("...xstepn2 is: " + xstepn2);
		Msg.debug("...ystepn2 is: " + ystepn2);

		// Try to find a location
		for (i = 1; i <= steps2n1 || i <= steps2n2; i++) {
			if (i <= steps2n1) {
				n.x = startX + xstepn1 * i;
				n.y = startY + ystepn1 * i;
				if (!q.anyInvertedElementsWhenCollapsed(n, n1, n2, l1, l2)) {
					Msg.debug("Leaving safeNewPosWhenCollapsingQuad(..): found");
					return n;
				}
			}
			if (i <= steps2n2) {
				n.x = startX + xstepn2 * i;
				n.y = startY + ystepn2 * i;
				if (!q.anyInvertedElementsWhenCollapsed(n, n1, n2, l1, l2)) {
					Msg.debug("Leaving safeNewPosWhenCollapsingQuad(..): found");
					return n;
				}
			}
		}

		Msg.debug("Leaving safeNewPosWhenCollapsingQuad(..): not found");
		return null;
	}

	/**
	 * To be used only with all-triangle meshes.
	 * 
	 * @return true if any zero area triangles were removed.
	 */
	boolean repairZeroAreaTriangles() {
		Msg.debug("Entering GeomBasics.repairZeroAreaTriangles()");
		boolean res = false;
		int j;
		Triangle t, old1, old2;
		Edge e, eS, e1, e2;

		for (int i = 0; i < triangleList.size(); i++) {
			if (!(triangleList.get(i) instanceof Triangle)) {
				continue;
			}
			t = (Triangle) triangleList.get(i);
			if (t.zeroArea()) {
				e = t.longestEdge();
				e1 = t.otherEdge(e);
				e2 = t.otherEdge(e, e1);
				res = true;

				Msg.debug("...longest edge is " + e.descr());
				if (!e.boundaryEdge()) {
					Msg.debug("...longest edge not on boundary!");
					old1 = (Triangle) e.element1;
					old2 = (Triangle) e.element2;
					eS = e.getSwappedEdge();
					e.swapToAndSetElementsFor(eS);

					triangleList.set(triangleList.indexOf(old1), null);
					triangleList.set(triangleList.indexOf(old2), null);

					triangleList.add((Triangle) eS.element1);
					triangleList.add((Triangle) eS.element2);

					edgeList.remove(edgeList.indexOf(e));
					edgeList.add(eS);
				} else {
					// The zero area triangle has its longest edge on the boundary...
					// Then we can just remove the triangle and the long edge!
					// Note that we now get a new boundary Vertex...
					Msg.debug("...longest edge is on boundary!");
					triangleList.set(triangleList.indexOf(t), null);
					t.disconnectEdges();
					edgeList.remove(edgeList.indexOf(e));
					e.disconnectVertices();
				}

			}
		}

		// Remove those entries that were set to null above.
		int i = 0;
		do {
			t = (Triangle) triangleList.get(i);
			if (t == null) {
				triangleList.remove(i);
			} else {
				i++;
			}
		} while (i < triangleList.size());

		Msg.debug("Leaving GeomBasics.repairZeroAreaTriangles()");
		return res;
	}

	/**
	 * A method for fast computation of the cross product of two vectors.
	 * 
	 * @param o1 origin of first vector
	 * @param p1 endpoint of first vector
	 * @param o2 origin of second vector
	 * @param p2 endpoint of second vector
	 * @return the cross product of the two vectors
	 */
	public static double cross(Vertex o1, Vertex p1, Vertex o2, Vertex p2) {
		double x1 = p1.x - o1.x;
		double x2 = p2.x - o2.x;
		double y1 = p1.y - o1.y;
		double y2 = p2.y - o2.y;
		return x1 * y2 - x2 * y1;
	}
}

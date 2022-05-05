package meshditor;

import java.util.ArrayList;
import java.util.List;

// ==== ---- ==== ---- ==== ---- ==== ---- ==== ---- ==== ---- ==== ---- 
/**
 * This class is an implementation of the algorithm described in the paper "An
 * approach to Combined Laplacian and Optimization-Based Smoothing for
 * Triangular, Quadrilateral and Quad-Dominant Meshes" (1998) by Cannan,
 * Tristano, and Staten.
 * 
 * The meshes produced by Q-Morph are indeed highly quad-dominant, with at most
 * one single triangle, so Q-Morph should work well with this algorithm.
 * 
 * Note that the boundary layer smoothing is not implemented.
 * 
 * @author Karl Erik Levik
 *
 */
// ==== ---- ==== ---- ==== ---- ==== ---- ==== ---- ==== ---- ==== ---- 

public class GlobalSmooth extends GeomBasics {
	
	public GlobalSmooth() {
	}

	/**
	 * Compute the constrained Laplacian smoothed position of a Vertex.
	 * 
	 * @param n the Vertex which is to be subjected to the smooth.
	 * @return the smoothed position of Vertex n.
	 */
	private Vertex constrainedLaplacianSmooth(Vertex n) {
		List<Element> elements = n.adjElements();
		Element oElem, sElem;
		MyVector vL = n.laplacianMoveVector();
		double deltaMy = 0, theta = 0, temp;
		int N = elements.size(), Nminus = 0, Nplus = 0, Nup = 0, Ndown = 0, Ninverted = 0;

		Vertex nLPos = new Vertex(n.x + vL.x, n.y + vL.y);
		nLPos.edgeList = n.edgeList;

		for (int count = 0; count < 20; count++) {
			deltaMy = 0.0;
			theta = 0.0;
			Nminus = 0;
			Nplus = 0;
			Nup = 0;
			Ndown = 0;
			Ninverted = 0;

			// Update the adjacent Elements' distortion metrics
			for (int i = 0; i < N; i++) {
				oElem = (Element) elements.get(i);

				sElem = oElem.elementWithExchangedVertices(n, nLPos);
				sElem.updateDistortionMetric();

				if (oElem.distortionMetric > sElem.distortionMetric) {
					Nminus++;
				} else if (oElem.distortionMetric < sElem.distortionMetric) {
					Nplus++;
				}

				deltaMy += sElem.distortionMetric - oElem.distortionMetric;

				// # elements whose metric improves significantly
				if ((oElem.distortionMetric < 0 && sElem.distortionMetric >= 0)
						|| (oElem.distortionMetric < 0 && sElem.distortionMetric > oElem.distortionMetric)
						|| (oElem.distortionMetric < MYMIN && sElem.distortionMetric >= MYMIN)) {
					Nup++;
				} else if ((oElem.distortionMetric >= 0 && sElem.distortionMetric < 0)
						|| (oElem.distortionMetric < 0 && sElem.distortionMetric < oElem.distortionMetric)
						|| (oElem.distortionMetric >= MYMIN && sElem.distortionMetric < MYMIN)) {
					Ndown++;
				}

				if (sElem.inverted()) {
					Ninverted++;
				}

				// For efficiency i could break out here... test if theta> THETAMAX here
				temp = sElem.largestAngle();
				if (temp > theta) {
					theta = temp;
				}
			}
			deltaMy = deltaMy / N;

			if (acceptable(N, Nminus, Nplus, Nup, Ndown, Ninverted, deltaMy, theta)) {
				// Return the proposed new position for Vertex n
				return nLPos;
			} else {
				vL = vL.div(2.0);
				nLPos = new Vertex(n.x + vL.x, n.y + vL.y);
			}
		}

		// Return the old position
		return n;
	}

	/**
	 * @return true if the new constrained-smoothed position is acceptable according
	 *         to the criteria given in section 4.2 of the article.
	 */
	private boolean acceptable(int N, int Nminus, int Nplus, int Nup, int Ndown, int Ninverted, double deltaMy, double theta) {

		if (Nminus == N || Ninverted > 0 || Ndown > Nup || deltaMy < MYMIN || theta > THETAMAX) {
			return false;
		} else if (Nplus == N || (Nup > 0 && Ndown == 0) || (Nup >= Ndown && deltaMy > MYMIN)) {
			return true;
		} else {
			return false;
		}
	}

	/**
	 * Compute the optimization-based smoothed position of a Vertex. As described in
	 * section 5 in the paper. Warning: The fields of the argument Vertex will be
	 * altered.
	 * 
	 * @return a Vertex with a position that is the optimaization-based smoothed
	 *         position of Vertex n.
	 */
	private Vertex optBasedSmooth(Vertex x, List<Element> elements) {
		Element sElem;
		double delta = Constants.DELTAFACTOR * maxModDim;
		Vertex xPX = new Vertex(x.x, x.y), xPY = new Vertex(x.x, x.y), xNew = new Vertex(x.x, x.y);
		double gX, gY;
		double minDM, newMinDM = java.lang.Double.MAX_VALUE;
		int iterations = 0;

		do { // Iterate until the min. dist. metric is acceptable
			minDM = java.lang.Double.MAX_VALUE;
			gX = 0.0;
			gY = 0.0;
			xPX.setXY(x.x + delta, x.y);
			xPY.setXY(x.x, x.y + delta);

			// Estimate the gradient vector g for each element:
			for (Element oElem : elements) {
				sElem = oElem.elementWithExchangedVertices(x, xPX);
				sElem.updateDistortionMetric();
				oElem.gX = (sElem.distortionMetric - oElem.distortionMetric) / delta;

				sElem = oElem.elementWithExchangedVertices(x, xPY);
				sElem.updateDistortionMetric();
				oElem.gY = (sElem.distortionMetric - oElem.distortionMetric) / delta;

				// Find the minimal DM and its gvec, but skip those with small gvecs
				if ((Math.abs(oElem.gX) > 0.00001 || Math.abs(oElem.gY) > 0.00001) && oElem.distortionMetric < minDM) {
					minDM = oElem.distortionMetric;
					gX = oElem.gX;
					gY = oElem.gY;
				}
			}
			

			// Which yields the final gradient vector g:
			MyVector g = new MyVector(x, gX, gY);

			// Determine gamma:
			// The dot product of two vectors u and v:
			// u dot v = sqrt(u.x^2 + u.y^2) * sqrt(v.x^2 + v.y^2) * cos(angle(u,v))
			boolean flag = false;
			double gamma = java.lang.Double.MAX_VALUE, gammaI = 0, gdotgi, gdotg = gX * gX + gY * gY;

			for (Element oElem : elements) {
				gdotgi = g.dot(oElem.gX, oElem.gY); // gX * oElem.gX + gY * oElem.gY;
				if (gdotgi < 0) {
					flag = true;
					gammaI = (oElem.distortionMetric - minDM) / (gdotg - gdotgi);
					if (gammaI < gamma) {
						gamma = gammaI;
					}
				}
			}

			if (!flag) {
				gamma = GAMMA; // I suppose something in the range (0, 1] so 0.8 is ok?
			}

			

			// Attempt the move x= x + gamma * g:
			xNew.setXY(x.x + gamma * gX, x.y + gamma * gY);

			for (int j = 0; j < 4; j++) {
				newMinDM = java.lang.Double.MAX_VALUE;

				for (Element oElem : elements) {
					sElem = oElem.elementWithExchangedVertices(x, xNew);
					sElem.updateDistortionMetric();
					oElem.newDistortionMetric = sElem.distortionMetric;

					if (sElem.distortionMetric < newMinDM) {
						newMinDM = sElem.distortionMetric;
					}
				}

				if (newMinDM > minDM + TOL) {
					break; // Escape from for-loop
				} else {
					gamma = gamma / 2.0;
					xNew.setXY(x.x + gamma * gX, x.y + gamma * gY);
				}
			}

			if (newMinDM > minDM + TOL) {
				x.setXY(xNew);
				x.update();
				minDM = newMinDM;

				// Update the adjacent Elements' distortion metrics
				for (Element oElem : elements) {
					oElem.distortionMetric = oElem.newDistortionMetric; // TODO ???
				}
			} else {
				return x;
			}
		} while (minDM <= OBSTOL && iterations++ <= 3); // Set max # of iterations

		return x;
	}

	private double maxModDim = 0.0;

	/** Initialize the object. */
	public void init() {
	}

	/** Perform the smoothing of the vertices in a step-wise manner. */
	@Override
	public void step() {
		run();
		setCurMethod(null);
	}

	/** The overall smoothing algorithm from section 3 in the paper. */
	public void run() {
		
		// Variables
		int i, j;
		List<Vertex> vertices = new ArrayList<>();
		List<Element> elements = new ArrayList<>();
		Element elem;
		Triangle t;
		double curLen, oldX, oldY;
		Vertex v, v_moved, n;

		// Get the internal vertices from vertexList.
		for (i = 0; i < vertexList.size(); i++) {
			v = (Vertex) vertexList.get(i);
			if (!v.boundaryVertex()) {
				vertices.add(v);
			}
		}

		for (i = 0; i < triangleList.size(); i++) {
			t = (Triangle) triangleList.get(i);
			t.updateDistortionMetric();
			// Find the largest edge length in the mesh
			curLen = t.longestEdgeLength();
			if (curLen > maxModDim) {
				maxModDim = curLen;
			}
		}

		for (i = 0; i < elementList.size(); i++) {
			elem = (Element) elementList.get(i);
			elem.updateDistortionMetric();
			// Find the largest edge length in the mesh
			curLen = elem.longestEdgeLength();
			if (curLen > maxModDim) {
				maxModDim = curLen;
			}
		}

		Edge e;

		double distance, maxMoveDistance = 0.0;
		int niter = 1;
		boolean VertexMoved;
		do {
			VertexMoved = false;
			for (i = 0; i < vertices.size(); i++) {
				v = (Vertex) vertices.get(i);

				if (v == null) {
					continue;
				}

				if (!v.movedByOBS) {
					v_moved = constrainedLaplacianSmooth(v);
					distance = v.length(v_moved);
					if (distance < Constants.MOVETOLERANCE) {
						vertices.set(i, null);
					} else {
						// Allow the move
						v.setXY(v_moved);
						v.update();
						VertexMoved = true;
						
						// Put neighbor vertices back in list, if they are not already there
						for (j = 0; j < v.edgeList.size(); j++) {
							e = (Edge) v.edgeList.get(j);
							n = e.otherVertex(v);
							if (!n.boundaryVertex() && !vertices.contains(n)) {
								vertices.add(n);
							}
						}
						// Keep track of the largest distance moved
						if (distance > maxMoveDistance) {
							maxMoveDistance = distance;
						}
						// Update the adjacent Elements' distortion metrics
						elements = v.adjElements();
						for (j = 0; j < elements.size(); j++) {
							elem = (Element) elements.get(j);
							elem.updateDistortionMetric();
						}
					}
				}
				if (niter >= 2) {
					
					// Find minimum distortion metric for the elements adjacent Vertex v
					elements = v.adjElements();
					elem = (Element) elements.get(0);
					double minDistMetric = elem.distortionMetric;
					for (j = 1; j < elements.size(); j++) {
						elem = (Element) elements.get(j);
						if (elem.distortionMetric < minDistMetric) {
							minDistMetric = elem.distortionMetric;
						}
					}
					if (minDistMetric <= OBSTOL) {
						oldX = v.x;
						oldY = v.y;
						v_moved = optBasedSmooth(v, elements);
						if (v_moved.x != oldX || v_moved.y != oldY) {
							// Put neighbor vertices back in list, if they're not there
							for (j = 0; j < v.edgeList.size(); j++) {
								e = (Edge) v.edgeList.get(j);
								n = e.otherVertex(v);
								if (!n.boundaryVertex() && !vertices.contains(n)) {
									vertices.add(n);
								}
							}
							// Keep track of the largest distance moved
							distance = v_moved.length(oldX, oldY);
							if (distance > maxMoveDistance) {
								maxMoveDistance = distance;
							}
							// Do the move
							v.setXY(v_moved);
							v.update();
							VertexMoved = true;
							v.movedByOBS = true; // Mark v as recently moved by OptBS
						} else {
							v.setXY(oldX, oldY);
							v.update();
						}
					} else {
						v.movedByOBS = false; // Mark v as not recently moved by OptBS
					}
				}
			}
			niter++;
		} while (VertexMoved && maxMoveDistance >= 1.75 * MOVETOLERANCE && niter < MAXITER);
	}

} // End of class GlobalSmooth

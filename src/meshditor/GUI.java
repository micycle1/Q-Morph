package meshditor;

import java.awt.CheckboxMenuItem;
import java.awt.Color;
import java.awt.FileDialog;
import java.awt.Font;
import java.awt.Frame;
import java.awt.Menu;
import java.awt.MenuBar;
import java.awt.MenuItem;
import java.awt.MenuShortcut;
import java.awt.ScrollPane;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.event.KeyEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.util.List;

/** This class implements the graphical user interface. */
public class GUI extends Constants implements ActionListener, ItemListener {

	/** Create frame, set font */
	public GUI() {
		f = new Frame("MeshDitor");
		// Font font= new Font("SansSerif", Font.PLAIN, 12);
		// f.setFont(font);
		f.setFont(new Font("Calibri", Font.PLAIN, 12));
		f.setIconImage(null); // MyIconImage.makeIconImage()
		GeomBasics.createNewLists();
	}

	/** Create frame, set font, instantiate QMorph */
	public GUI(String dir, String filename) {
		f = new Frame("MeshDitor: " + filename);
		// Font font= new Font("SansSerif", Font.PLAIN, 12);
		// f.setFont(font);
		f.setFont(new Font("Calibri", Font.PLAIN, 12));
		f.setIconImage(null); // MyIconImage.makeIconImage()

		this.filename = filename;
		// Parameters params= new Parameters(filename, false, false);
		GeomBasics.setParams(filename, dir, false, false);
		// qm= new QuadMorph(params);

		GeomBasics.loadMesh();
		GeomBasics.findExtremeVertices();
		// edgeList= GeomBasics.getEdgeList();
		// vertexList= GeomBasics.getvertexList();
	}

	/** The filename of the current mesh. */
	public String filename;
	/** Boolean indicating that we are currently defining vertices. */
	public boolean VertexMode = false;
	/** Boolean indicating that we are currently defining triangles. */
	public boolean triangleMode = true;
	/** Boolean indicating that we are currently defining quads. */
	public boolean quadMode = false;
	/** Boolean indicating visibility for the grid */
	public boolean grid = true;
	/** Boolean indicating visibility for the axis */
	public boolean axis = true;

	/** Pointer to an instance of the QMorph class. */
	public QMorph qm = null;
	/** Pointer to an instance of the DelaunayMeshGen class. */
	public DelaunayMeshGen tri = null;

	Frame f;
	private GCanvas cvas;
	private GControls gctrls;
	private ScrollPane sp;
	private MenuBar mb;

	MenuItem mi;
	Menu fileMenu, editMenu, modeMenu, debugMenu, runMenu, helpMenu;
	int width = 640, height = 480;
	int scale = 100;
	MyMouseListener myMouseListener;

	MenuItem newItem, loadMeshItem, loadVerticesItem, saveItem, saveAsItem, saveVerticesItem, saveVerticesAsItem, saveTriAsItem, exportItem,
			exitItem;
	MenuItem undoItem, clearEdgesItem;
	CheckboxMenuItem VertexModeItem, triModeItem, quadModeItem, debugModeItem, stepModeItem;
	MenuItem consistencyItem, detectInversionItem, printElementsItem, printTrianglesItem, reportMetricsItem, printValencesItem,
			printValPatItem, printAngAtSurVerticesItem, centroidItem, triCountItem, delauneyItem, qmorphItem, globalCleanUpItem,
			globalSmoothItem, helpItem, aboutItem;

	MenuShortcut qkey;

	/** Start up the GUI. */
	public void startGUI() {
		f.setSize(width, height);
		fileMenu = new Menu("File");
		newItem = new MenuItem("New");
		loadMeshItem = new MenuItem("Load mesh");
		loadVerticesItem = new MenuItem("Load vertices");
		saveItem = new MenuItem("Save mesh");
		saveAsItem = new MenuItem("Save mesh as...");
		saveVerticesItem = new MenuItem("Save vertices");
		saveVerticesAsItem = new MenuItem("Save vertices as...");
		saveTriAsItem = new MenuItem("Save triangle mesh as...");
		exportItem = new MenuItem("Export mesh to LaTeX file");
		exitItem = new MenuItem("Exit");
		qkey = new MenuShortcut(KeyEvent.VK_Q, false);
		exitItem.setShortcut(qkey);

		newItem.addActionListener(this);
		loadMeshItem.addActionListener(this);
		loadVerticesItem.addActionListener(this);
		saveItem.addActionListener(this);
		saveVerticesItem.addActionListener(this);
		saveVerticesAsItem.addActionListener(this);
		saveAsItem.addActionListener(this);
		saveTriAsItem.addActionListener(this);
		exportItem.addActionListener(this);
		exitItem.addActionListener(this);

		fileMenu.add(newItem);
		fileMenu.add(loadMeshItem);
		fileMenu.add(loadVerticesItem);
		fileMenu.add(saveItem);
		fileMenu.add(saveAsItem);
		fileMenu.add(saveVerticesItem);
		fileMenu.add(saveVerticesAsItem);
		fileMenu.add(saveTriAsItem);
		fileMenu.add(exportItem);
		fileMenu.addSeparator();
		fileMenu.add(exitItem);

		editMenu = new Menu("Edit");
		undoItem = new MenuItem("Undo last Vertex/edge creation or move");
		undoItem.addActionListener(this);
		editMenu.add(undoItem);
		clearEdgesItem = new MenuItem("Clear all edges");
		clearEdgesItem.addActionListener(this);
		editMenu.add(clearEdgesItem);

		modeMenu = new Menu("Mode");

		VertexModeItem = new CheckboxMenuItem("Plot vertices");
		VertexModeItem.setState(false);
		VertexModeItem.addItemListener(this);
		modeMenu.add(VertexModeItem);

		triModeItem = new CheckboxMenuItem("Construct triangles");
		triModeItem.setState(true);
		triModeItem.addItemListener(this);
		modeMenu.add(triModeItem);
		quadModeItem = new CheckboxMenuItem("Construct quads");
		quadModeItem.setState(false);
		quadModeItem.addItemListener(this);
		modeMenu.add(quadModeItem);
		modeMenu.addSeparator();
		debugModeItem = new CheckboxMenuItem("Debug mode");
		debugModeItem.setState(Msg.debugMode);
		debugModeItem.addItemListener(this);
		modeMenu.add(debugModeItem);
		stepModeItem = new CheckboxMenuItem("Step mode");
		stepModeItem.setState(false);
		stepModeItem.addItemListener(this);
		modeMenu.add(stepModeItem);

		debugMenu = new Menu("Debug");
		consistencyItem = new MenuItem("Test consistency of mesh");
		consistencyItem.addActionListener(this);
		detectInversionItem = new MenuItem("Detect inverted elements");
		detectInversionItem.addActionListener(this);
		printTrianglesItem = new MenuItem("Print triangleList");
		printTrianglesItem.addActionListener(this);
		printElementsItem = new MenuItem("Print elementList");
		printElementsItem.addActionListener(this);
		reportMetricsItem = new MenuItem("Report mesh metrics");
		reportMetricsItem.addActionListener(this);
		printValencesItem = new MenuItem("Print valences of all vertices");
		printValencesItem.addActionListener(this);
		printValPatItem = new MenuItem("Print valence patterns of all vertices");
		printValPatItem.addActionListener(this);
		printAngAtSurVerticesItem = new MenuItem("Print angles at surrounding vertices");
		printAngAtSurVerticesItem.addActionListener(this);

		centroidItem = new MenuItem("Create centroid for last quad");
		centroidItem.addActionListener(this);
		triCountItem = new MenuItem("Count triangles");
		triCountItem.addActionListener(this);

		debugMenu.add(consistencyItem);
		debugMenu.add(detectInversionItem);
		debugMenu.add(printTrianglesItem);
		debugMenu.add(printElementsItem);
		debugMenu.add(reportMetricsItem);
		debugMenu.add(printValencesItem);
		debugMenu.add(printValPatItem);
		debugMenu.add(printAngAtSurVerticesItem);
		debugMenu.add(centroidItem);
		debugMenu.add(triCountItem);

		runMenu = new Menu("Run");
		qmorphItem = new MenuItem("Run QMorph");
		delauneyItem = new MenuItem("Run Delauney generator");
		// globalCleanUpItem= new MenuItem("Run topological cleanup");
		// globalSmoothItem= new MenuItem("Run smooth");
		delauneyItem.addActionListener(this);
		qmorphItem.addActionListener(this);
		// globalCleanUpItem.addActionListener(this);
		// globalSmoothItem.addActionListener(this);

		runMenu.add(qmorphItem);
		runMenu.add(delauneyItem);
		// runMenu.add(globalCleanUpItem);
		// runMenu.add(globalSmoothItem);

		helpMenu = new Menu("Help");
		helpItem = new MenuItem("Help");
		aboutItem = new MenuItem("About");
		helpItem.addActionListener(this);
		aboutItem.addActionListener(this);

		helpMenu.add(helpItem);
		helpMenu.add(aboutItem);

		mb = new MenuBar();
		mb.add(fileMenu);
		mb.add(editMenu);
		mb.add(modeMenu);
		mb.add(debugMenu);
		mb.add(runMenu);
		mb.setHelpMenu(helpMenu);

		f.setMenuBar(mb);

		f.setBackground(Color.lightGray);
		f.setForeground(Color.black);

		if (GeomBasics.leftmost == null) {
			cvas = new GCanvas(this, scale);
		} else {
			cvas = new GCanvas(this, GeomBasics.leftmost.x, GeomBasics.lowermost.y, GeomBasics.rightmost.x, GeomBasics.uppermost.y, scale);
		}

		myMouseListener = new MyMouseListener();
		cvas.addMouseListener(myMouseListener);
		cvas.repaint();

		f.add("South", gctrls = new GControls(this, cvas));
		f.add("Center", sp = new ScrollPane(ScrollPane.SCROLLBARS_ALWAYS));
		sp.setForeground(Color.darkGray);
		sp.setBackground(Color.lightGray);

		sp.add(cvas);

		cvas.setForeground(Color.black);
		cvas.setBackground(Color.black);

		f.setVisible(true);
	}

	void commandNew() {
		GeomBasics.clearLists();
		GeomBasics.setParams(null, ".", false, false);
		f.setTitle("MeshDitor:");
		cvas.resize(-2, -2, 2, 2, 100);
		cvas.clear();
		qm = null;
	}

	void commandLoadMesh() {
		FileDialog fd = new FileDialog(f, "Load mesh from file", FileDialog.LOAD);
		fd.setDirectory(GeomBasics.meshDirectory);
		fd.show();
		String dir = fd.getDirectory();
		String loadName = fd.getFile();
		if (dir != null && dir != "" && loadName != null && loadName != "") {
			GeomBasics.clearLists();
			cvas.clear();

			GeomBasics.setParams(loadName, dir, false, false);
			GeomBasics.loadMesh();
			f.setTitle("MeshDitor: " + loadName);

			GeomBasics.findExtremeVertices();
			cvas.resize(GeomBasics.leftmost.x, GeomBasics.lowermost.y, GeomBasics.rightmost.x, GeomBasics.uppermost.y, scale);
		}
	}

	void commandLoadVertices() {
		FileDialog fd = new FileDialog(f, "Load vertices from file", FileDialog.LOAD);
		fd.setDirectory(GeomBasics.meshDirectory);
		fd.show();
		String dir = fd.getDirectory();
		String loadName = fd.getFile();
		if (dir != null && dir != "" && loadName != null && loadName != "") {
			GeomBasics.clearLists();
			cvas.clear();

			GeomBasics.setParams(loadName, dir, false, false);
			GeomBasics.loadVertices();
			f.setTitle("MeshDitor: " + loadName);

			GeomBasics.findExtremeVertices();
			cvas.resize(GeomBasics.leftmost.x, GeomBasics.lowermost.y, GeomBasics.rightmost.x, GeomBasics.uppermost.y, scale);
		}
	}

	void commandSaveMesh() {
		if (filename == null || filename == "") {
			commandSaveMeshAs();
		} else {
			GeomBasics.writeMesh(GeomBasics.meshDirectory + GeomBasics.meshFilename);
		}
	}

	void commandSaveVertices() {
		if (filename == null || filename == "") {
			commandSaveVerticesAs();
		} else {
			GeomBasics.writeMesh(GeomBasics.meshDirectory + GeomBasics.meshFilename);
		}
	}

	void commandSaveVerticesAs() {
		FileDialog fd = new FileDialog(f, "Save vertices to file", FileDialog.SAVE);
		fd.setDirectory(GeomBasics.meshDirectory);
		fd.show();
		String dir = fd.getDirectory();
		String saveName = fd.getFile();
		if (dir != null && dir != "" && saveName != null && saveName != "") {
			GeomBasics.writeVertices(dir + saveName);
			f.setTitle("MeshDitor: " + saveName);
			filename = dir + saveName;
			GeomBasics.setParams(saveName, dir, false, false);
		}
	}

	void commandSaveMeshAs() {
		FileDialog fd = new FileDialog(f, "Save mesh to file", FileDialog.SAVE);
		fd.setDirectory(GeomBasics.meshDirectory);
		fd.show();
		String dir = fd.getDirectory();
		String saveName = fd.getFile();
		if (dir != null && dir != "" && saveName != null && saveName != "") {
			GeomBasics.writeMesh(dir + saveName);
			f.setTitle("MeshDitor: " + saveName);
			filename = dir + saveName;
			GeomBasics.setParams(saveName, dir, false, false);
		}
	}

	void commandSaveTriangleMeshAs() {
		FileDialog fd = new FileDialog(f, "Save triangle mesh to file", FileDialog.SAVE);
		fd.setDirectory(GeomBasics.meshDirectory);
		fd.show();
		String dir = fd.getDirectory();
		String saveName = fd.getFile();
		if (dir != null && dir != "" && saveName != null && saveName != "") {
			GeomBasics.writeQuadMesh(saveName, GeomBasics.triangleList);
			f.setTitle("MeshDitor: " + saveName);
			filename = saveName;
		}
	}

	void commandExportMeshToLaTeX() {
		int ul;
		double xcorr, ycorr;
		boolean visibleVertices;

		ExportToLaTeXOptionsDialog ed = new ExportToLaTeXOptionsDialog(f, "Set export parameters", true);
		ed.show();

		if (ed.okPressed()) {
			ul = ed.getUnitlength();
			xcorr = ed.getXCorr();
			ycorr = ed.getYCorr();
			visibleVertices = ed.getVisibleVertices();

			FileDialog fd = new FileDialog(f, "Export mesh to file", FileDialog.SAVE);
			fd.setVisible(true);
			String dir = fd.getDirectory();
			String saveName = fd.getFile();
			if (dir != null && dir != "" && saveName != null && saveName != "") {
				GeomBasics.exportMeshToLaTeX(dir + saveName, ul, xcorr, ycorr, visibleVertices);
			}
		}
	}

	void commandUndo() {
		myMouseListener.undo();
		cvas.repaint();
	}

	void commandClearEdges() {
		GeomBasics.clearEdges();
		cvas.repaint();
	}

	void commandVertexMode() {
		VertexMode = true;
		triangleMode = false;
		quadMode = false;
		VertexModeItem.setState(true);
		triModeItem.setState(false);
		quadModeItem.setState(false);
		gctrls.clickStatus.setText("1");
	}

	void commandTriMode() {
		VertexMode = false;
		triangleMode = true;
		quadMode = false;
		VertexModeItem.setState(false);
		triModeItem.setState(true);
		quadModeItem.setState(false);
		gctrls.clickStatus.setText("3");
	}

	void commandQuadMode() {
		VertexMode = false;
		triangleMode = false;
		quadMode = true;
		VertexModeItem.setState(false);
		triModeItem.setState(false);
		quadModeItem.setState(true);
		gctrls.clickStatus.setText("4");
	}

	void commandToggleDebugMode() {
		if (Msg.debugMode) {
			Msg.debugMode = false;
			debugModeItem.setState(false);
		} else {
			Msg.debugMode = true;
			debugModeItem.setState(true);
		}
	}

	void commandToggleStepMode() {
		if (GeomBasics.step) {
			GeomBasics.step = false;
			stepModeItem.setState(false);
		} else {
			GeomBasics.step = true;
			stepModeItem.setState(true);
		}
	}

	void commandQMorph() {
		QMorphOptionsDialog qmod = new QMorphOptionsDialog(f);
		qmod.setSize(qmod.getPreferredSize());
		qmod.show();

		if (qmod.runPressed()) {
			qmod.copyToProgramParameters();
			qm = new QMorph();
			qm.init();

			if (!GeomBasics.step) {
				qm.run();
				GeomBasics.findExtremeVertices();
				cvas.resize(GeomBasics.leftmost.x, GeomBasics.lowermost.y, GeomBasics.rightmost.x, GeomBasics.uppermost.y, scale);
			} else {
			}
		}
	}

	void commandDelaunay() {
		tri = new DelaunayMeshGen();
		tri.init(true); // false

		if (!GeomBasics.step) {
			/* Straightforwardly run the method */
			tri.run();
			// elementList= tri.getTriangleList(); // tri.incrDelauney(vertexList);
			// edgeList= tri.getEdgeList();
			// vertexList= tri.getvertexList();
			cvas.resize(GeomBasics.leftmost.x, GeomBasics.lowermost.y, GeomBasics.rightmost.x, GeomBasics.uppermost.y, scale);
		} else {
			/* Run method in step mode */
			cvas.resize(GeomBasics.leftmost.x, GeomBasics.lowermost.y, GeomBasics.rightmost.x, GeomBasics.uppermost.y, scale);
		}
	}

	void commandTopoCleanup() {
		if (GeomBasics.topoCleanup == null) {
			GeomBasics.topoCleanup = new TopoCleanup();
		}

		GeomBasics.topoCleanup.init();
		if (!GeomBasics.step) {
			GeomBasics.topoCleanup.run();
			cvas.repaint();
		}
	}

	void commandSmooth() {
		if (GeomBasics.globalSmooth == null) {
			GeomBasics.globalSmooth = new GlobalSmooth();
		}

		GeomBasics.globalSmooth.init();
		GeomBasics.globalSmooth.run();
		cvas.repaint();
	}

	/**
	 * Invoked when a registered item change occurs. The item is identified, and the
	 * corresponding action is invoked.
	 */
	@Override
	public void itemStateChanged(ItemEvent e) {
		String command = (String) e.getItem();
		if (command.equals("Plot vertices")) {
			commandVertexMode();
		} else if (command.equals("Construct triangles")) {
			commandTriMode();
		} else if (command.equals("Construct quads")) {
			commandQuadMode();
		} else if (command.equals("Debug mode")) {
			commandToggleDebugMode();
		} else if (command.equals("Step mode")) {
			commandToggleStepMode();
		}
	}

	/**
	 * Invoked when a registered action command occurs. The command is identified,
	 * and the corresponding action is invoked.
	 */
	@Override
	public void actionPerformed(ActionEvent e) {
		HelpDialog hd;
		AboutDialog ad;
		MsgDialog rd;

		String command = e.getActionCommand();

		if (command.equals("New")) {
			commandNew();
		} else if (command.equals("Load mesh")) {
			commandLoadMesh();
		} else if (command.equals("Load vertices")) {
			commandLoadVertices();
		} else if (command.equals("Save mesh")) {
			commandSaveMesh();
		} else if (command.equals("Save vertices")) {
			commandSaveVertices();
		} else if (command.equals("Save vertices as...")) {
			commandSaveVerticesAs();
		} else if (command.equals("Save mesh as...")) {
			commandSaveMeshAs();
		} else if (command.equals("Save triangle mesh as...")) {
			commandSaveTriangleMeshAs();
		} else if (command.equals("Export mesh to LaTeX file")) {
			commandExportMeshToLaTeX();
		} else if (command.equals("Exit")) {
			System.exit(0);
		} else if (command.equals("Undo last Vertex/edge creation or move")) {
			commandUndo();
		} else if (command.equals("Clear all edges")) {
			commandClearEdges();
		}

		else if (command.equals("Test consistency of mesh")) {
			GeomBasics.consistencyCheck();
		} else if (command.equals("Detect inverted elements")) {
			GeomBasics.detectInvertedElements();
			cvas.repaint();
		} else if (command.equals("Print triangleList")) {
			GeomBasics.printTriangles(GeomBasics.getTriangleList());
		} else if (command.equals("Print elementList")) {
			GeomBasics.printQuads(GeomBasics.getElementList());
		}
		// else if (command.equals("Update mesh metrics")) {
		// GeomBasics.updateMeshMetrics();
		// }
		else if (command.equals("Report mesh metrics")) {
			GeomBasics.updateMeshMetrics();
			rd = new MsgDialog(f, "Mesh Metrics Report", GeomBasics.meshMetricsReport(), 80, 18);
			rd.show();
		} else if (command.equals("Print valences of all vertices")) {
			GeomBasics.printValences();
		} else if (command.equals("Print valence patterns of all vertices")) {
			GeomBasics.printValPatterns();
		} else if (command.equals("Print angles at surrounding vertices")) {
			GeomBasics.printAnglesAtSurrondingVertices();
		} else if (command.equals("Create centroid for last quad")) {

			Vertex n;
			Quad q;
			Element elem;
			int size = GeomBasics.elementList.size();
			if (size > 0) {
				elem = GeomBasics.elementList.get(size - 1);
				if (elem instanceof Quad) {
					q = (Quad) elem;
					n = q.centroid();
					GeomBasics.vertexList.add(n);
					cvas.repaint();
				}
			}
		} else if (command.equals("Count triangles")) {
			GeomBasics.countTriangles();
		}

		else if (command.equals("Run QMorph")) {
			commandQMorph();
		} else if (command.equals("Run Delauney generator")) {
			commandDelaunay();
		} else if (command.equals("Run topological cleanup")) {
			commandTopoCleanup();
		} else if (command.equals("Run smooth")) {
			commandSmooth();
		} else if (command.equals("Help")) {
			hd = new HelpDialog(f);
			hd.show();
		} else if (command.equals("About")) {
			ad = new AboutDialog(f);
			ad.show();
		}
	}

	/** A class for handling mouse actions. */
	class MyMouseListener extends MouseAdapter {
		Vertex movingVertex = null, oldMovingVertex = null;
		int VertexCnt = 0;
		Edge edge1, edge2, edge3, edge4;
		Vertex[] myvertexList = new Vertex[4];
		Triangle tri;
		Quad q;
		boolean lastActionMergeVertices = false;
		boolean lastActionMoveVertex = false;
		boolean lastActionNewVertex = false;
		boolean lastActionNewEdge = false;
		boolean lastActionTwoNewEdges = false;
		boolean lastActionNewTriangle = false;
		boolean lastActionNewQuad = false;

		double oldX = 0, oldY = 0;
		int nONewEdges = 0;

		public MyMouseListener() {
		}

		/** Invoked when the mouse has been clicked on a component. */
		@Override
		public void mouseClicked(MouseEvent e) {
			Edge b, l, r, t;
			lastActionMoveVertex = false;
			lastActionMergeVertices = false;
			lastActionNewVertex = false;
			lastActionNewEdge = false;
			lastActionTwoNewEdges = false;
			lastActionNewTriangle = false;
			lastActionNewQuad = false;
			nONewEdges = 0;

			double x = Math.rint(e.getX() / 10.0) * 10;
			double y = Math.rint(e.getY() / 10.0) * 10;
			x = (x - cvas.getYAxisXPos()) / scale;
			y = (y - cvas.getXAxisYPos()) / -scale;
			Vertex n = new Vertex(x, y);
			if (!GeomBasics.vertexList.contains(n)) {
				GeomBasics.vertexList.add(n);
				lastActionNewVertex = true;
			} else {
				n = GeomBasics.vertexList.get(GeomBasics.vertexList.indexOf(n));
			}

			if (!VertexMode) {
				myvertexList[VertexCnt] = n;
				VertexCnt++;
				if (VertexCnt == 2) {
					if (myvertexList[0] == myvertexList[1]) {
						VertexCnt = 1;
						return;
					}

					edge1 = new Edge(myvertexList[0], myvertexList[1]);
					if (!GeomBasics.edgeList.contains(edge1)) {
						GeomBasics.edgeList.add(edge1);
						edge1.connectVertices();
						nONewEdges++;
						lastActionNewEdge = true;
					} else {
						edge1 = GeomBasics.edgeList.get(GeomBasics.edgeList.indexOf(edge1));
					}
				} else if (VertexCnt == 3 && triangleMode) {

					if (myvertexList[2] == myvertexList[0] || myvertexList[2] == myvertexList[1]) {
						VertexCnt = 2;
						return;
					}

					edge2 = new Edge(myvertexList[1], myvertexList[2]);
					if (!GeomBasics.edgeList.contains(edge2)) {
						GeomBasics.edgeList.add(edge2);
						edge2.connectVertices();
						nONewEdges++;
						if (nONewEdges >= 2) {
							lastActionTwoNewEdges = true;
							lastActionNewEdge = false;
						} else {
							lastActionNewEdge = true;
						}
					} else {
						edge2 = GeomBasics.edgeList.get(GeomBasics.edgeList.indexOf(edge2));
					}

					edge3 = new Edge(myvertexList[0], myvertexList[2]);
					if (!GeomBasics.edgeList.contains(edge3)) {
						GeomBasics.edgeList.add(edge3);
						edge3.connectVertices();
						nONewEdges++;
						if (nONewEdges >= 2) {
							lastActionTwoNewEdges = true;
							lastActionNewEdge = false;
						} else {
							lastActionNewEdge = true;
						}
					} else {
						edge3 = GeomBasics.edgeList.get(GeomBasics.edgeList.indexOf(edge3));
					}

					tri = new Triangle(edge1, edge2, edge3);

					if (!GeomBasics.triangleList.contains(tri)) {
						GeomBasics.triangleList.add(tri);
						tri.connectEdges();
						lastActionNewTriangle = true;
					}
					VertexCnt = 0;
				} else if (VertexCnt == 3 && quadMode) {
					if (myvertexList[2] == myvertexList[0] || myvertexList[2] == myvertexList[1]) {
						VertexCnt = 2;
						return;
					}

					edge2 = new Edge(myvertexList[1], myvertexList[2]);
					if (!GeomBasics.edgeList.contains(edge2)) {
						GeomBasics.edgeList.add(edge2);
						edge2.connectVertices();
						nONewEdges++;
						lastActionNewEdge = true;
					} else {
						edge2 = GeomBasics.edgeList.get(GeomBasics.edgeList.indexOf(edge2));
					}
				} else if (VertexCnt == 4 && quadMode) {
					if (myvertexList[3] == myvertexList[0] || myvertexList[3] == myvertexList[1] || myvertexList[3] == myvertexList[2]) {
						VertexCnt = 3;
						return;
					}

					edge3 = new Edge(myvertexList[2], myvertexList[3]);
					if (!GeomBasics.edgeList.contains(edge3)) {
						GeomBasics.edgeList.add(edge3);
						edge3.connectVertices();
						nONewEdges++;
						if (nONewEdges >= 2) {
							lastActionTwoNewEdges = true;
							lastActionNewEdge = false;
						} else {
							lastActionNewEdge = true;
						}
					} else {
						edge3 = GeomBasics.edgeList.get(GeomBasics.edgeList.indexOf(edge3));
					}

					edge4 = new Edge(myvertexList[0], myvertexList[3]);
					if (!GeomBasics.edgeList.contains(edge4)) {
						GeomBasics.edgeList.add(edge4);
						edge4.connectVertices();
						nONewEdges++;
						if (nONewEdges >= 2) {
							lastActionTwoNewEdges = true;
							lastActionNewEdge = false;
						} else {
							lastActionNewEdge = true;
						}
					} else {
						edge4 = GeomBasics.edgeList.get(GeomBasics.edgeList.indexOf(edge4));
					}

					// Decide which is base, left, right, and top:
					if (edge2.hasVertex(edge1.leftVertex)) {
						b = edge1;
						l = edge2;
						t = edge3;
						r = edge4;
					} else if (edge2.hasVertex(edge1.rightVertex)) {
						b = edge1;
						r = edge2;
						t = edge3;
						l = edge4;
					} else {
						Msg.error("Weird stuff while creating new quad...");
						return;
					}
					q = new Quad(b, l, r, t);

					if (!GeomBasics.elementList.contains(q)) {
						GeomBasics.elementList.add(q);
						q.connectEdges();
						/*
						 * b.connectVertices(); l.connectVertices(); r.connectVertices();
						 * t.connectVertices();
						 */
						lastActionNewQuad = true;
					}
					VertexCnt = 0;
				}

				if (VertexMode) {
					gctrls.clickStatus.setText("1");
				} else if (triangleMode) {
					gctrls.clickStatus.setText(Integer.toString(3 - VertexCnt));
				} else if (quadMode) {
					gctrls.clickStatus.setText(Integer.toString(4 - VertexCnt));
				}
			}
			cvas.repaint();
		}

		/**
		 * Invoked when a mouse button is pressed (but not yet released). If it is
		 * pressed on a particular Vertex, then remember which.
		 */
		@Override
		public void mousePressed(MouseEvent e) {
			double x = Math.rint(e.getX() / 10.0) * 10;
			double y = Math.rint(e.getY() / 10.0) * 10;

			x = (x - cvas.getYAxisXPos()) / scale;
			y = (y - cvas.getXAxisYPos()) / -scale;

			movingVertex = new Vertex(x, y);
			int j = GeomBasics.vertexList.indexOf(movingVertex);
			if (j != -1) {
				movingVertex = GeomBasics.vertexList.get(j);
				oldX = movingVertex.x;
				oldY = movingVertex.y;
			} else {
				movingVertex = null;
			}
		}

		/** Invoked when a mouse button is released (after being pressed). */
		@Override
		public void mouseReleased(MouseEvent e) {
			Edge ei, ej, oldE;
			Vertex n, other;
			int ind, j, k;

			double x = Math.rint(e.getX() / 10.0) * 10;
			double y = Math.rint(e.getY() / 10.0) * 10;
			x = (x - cvas.getYAxisXPos()) / scale;
			y = (y - cvas.getXAxisYPos()) / -scale;

			if (movingVertex != null && (x != movingVertex.x || y != movingVertex.y)) {
				oldMovingVertex = movingVertex.copy();
				ind = GeomBasics.vertexList.indexOf(new Vertex(x, y));
				movingVertex.setXY(x, y);
				movingVertex.update();

				if (ind != -1) { // We have to merge the vertices
					n = GeomBasics.vertexList.get(ind);

					for (int i = 0; i < n.edgeList.size(); i++) {
						ei = n.edgeList.get(i);
						j = movingVertex.edgeList.indexOf(ei);
						if (j == -1) {
							ei.replaceVertex(n, movingVertex);
							movingVertex.edgeList.add(ei);
							// ei.connectTo()
						} else { // keep only one copy of each edge
							// if (ei.leftVertex== ei.rightVertex) {
							ej = movingVertex.edgeList.get(j);

							if (ej.element1 != null) {
								ej.element1.replaceEdge(ej, ei);
								ei.connectToTriangle((Triangle) ej.element1);
							}
							if (ej.element2 != null) {
								ej.element2.replaceEdge(ej, ei);
								ei.connectToTriangle((Triangle) ej.element2);
							}

							if (ei.element1 != null) {
							}
							if (ei.element2 != null) {
							}

							movingVertex.edgeList.set(j, ei);
							ei.replaceVertex(n, movingVertex);
							other = ei.otherVertex(movingVertex);

							// Remove the correct edge from "global" edgelist
							k = GeomBasics.edgeList.indexOf(ej);
							oldE = GeomBasics.edgeList.get(k);
							if (oldE == ei) {
								k = GeomBasics.edgeList.lastIndexOf(oldE);
							}
							GeomBasics.edgeList.remove(k);

							// Remove the correct edge from other's edgelist
							k = other.edgeList.indexOf(ej);
							oldE = other.edgeList.get(k);
							if (oldE == ei) {
								k = other.edgeList.lastIndexOf(oldE);
							}
							other.edgeList.remove(k);
							// }
						}
					}

					List<Element> aeList = movingVertex.adjElements();
					Element elem;
					for (Object element : aeList) {
						elem = (Element) element;
						elem.updateAngles();
					}

					GeomBasics.vertexList.remove(ind);
					lastActionMergeVertices = true;
				} else {
					lastActionMergeVertices = false;
				}

				lastActionMoveVertex = true;
				lastActionNewVertex = false;
				lastActionNewEdge = false;
				lastActionTwoNewEdges = false;
				lastActionNewTriangle = false;

				cvas.repaint();
			}
		}

		/** Undo the last mouse action. */
		public void undo() {
			Edge e;
			int j;
			if (lastActionMoveVertex) {
				if (lastActionMergeVertices) {
					for (int i = 0; i < oldMovingVertex.edgeList.size(); i++) {
						e = oldMovingVertex.edgeList.get(i);
						j = movingVertex.edgeList.indexOf(e);
						if (j != -1) {
							if (!(e.hasVertex(movingVertex) && e.hasVertex(oldMovingVertex))) {
								movingVertex.edgeList.remove(j);
								e.replaceVertex(movingVertex, oldMovingVertex);
							}
						} else { // common edge that collapsed because of the merge
							movingVertex.edgeList.add(e);
							GeomBasics.edgeList.add(e);
						}
					}
					GeomBasics.vertexList.add(oldMovingVertex);
				} else {
					movingVertex.setXY(oldX, oldY);
					movingVertex.update();
				}
			}

			else if (triangleMode && (lastActionNewVertex || lastActionNewEdge)) {
				if (VertexCnt == 0) {
					VertexCnt = 2;
				} else {
					VertexCnt--;
				}
			}

			else if (quadMode && (lastActionNewVertex || lastActionNewEdge)) {
				if (VertexCnt == 0) {
					VertexCnt = 3;
				} else {
					VertexCnt--;
				}
			}

			if (lastActionNewVertex) {
				GeomBasics.vertexList.remove(GeomBasics.vertexList.size() - 1);
			}
			if (lastActionNewEdge) {
				e = GeomBasics.edgeList.get(GeomBasics.edgeList.size() - 1);
				e.disconnectVertices();
				GeomBasics.edgeList.remove(GeomBasics.edgeList.size() - 1);
			} else if (lastActionTwoNewEdges) {
				e = GeomBasics.edgeList.get(GeomBasics.edgeList.size() - 1);
				e.disconnectVertices();
				GeomBasics.edgeList.remove(GeomBasics.edgeList.size() - 1);
				e = GeomBasics.edgeList.get(GeomBasics.edgeList.size() - 1);
				e.disconnectVertices();
				GeomBasics.edgeList.remove(GeomBasics.edgeList.size() - 1);
			}

			if (lastActionNewTriangle) {
				tri.disconnectEdges();
				GeomBasics.triangleList.remove(GeomBasics.triangleList.size() - 1);
			} else if (lastActionNewQuad) {
				q.disconnectEdges();
				GeomBasics.elementList.remove(GeomBasics.elementList.size() - 1);
			}

			lastActionMoveVertex = false;
			lastActionNewVertex = false;
			lastActionNewEdge = false;
			lastActionTwoNewEdges = false;
			lastActionNewTriangle = false;
			lastActionNewQuad = false;
		}

	}

}

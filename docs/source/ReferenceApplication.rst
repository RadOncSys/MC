Reference application
=====================

Transport modules coding in XML geometry file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


.. list-table:: List of transport module types
   :name: _module_types-table
   :widths: 4, 6, 6
   :width: 100%
   :header-rows: 1

   * - Type name
     - Transport class
     - Description
   * - *cylinder*
     - *mcTransportCylinder*
     - Cylinder
   * - *cone*
     - *mcTransportCone*
     - 
   * - *prism*
     - *mcTransportPrism*
     - 
   * - *wedge*
     - *mcTransportWedge*
     - 
   * - *ring*
     - *mcTransportRing*
     - 
   * - *conicalring*
     - *mcTransportConicalRing*
     - 
   * - *conicalhole*
     - *mcTransportConicalHole*
     - 
   * - *rectanglering*
     - *mcTransportRectangleRing*
     - 
   * - *jaw_pair*
     - *mcTransportJawPairRounded*
     - 
   * - *jaw_pair_focused*
     - *mcTransportJawPairFocused*
     - 
   * - *rectanglepolygonsidehole*
     - *mcTransportRectanglePolygonSideHole*
     - 
   * - *mlc*
     - *mcTransportMLC*
     - 
   * - *gridcylinder*
     - *mcTransportCylinderStack*
     - 
   * - *planefilter*
     - *mcTransportPlaneFilter*
     - 
   * - *etrap*
     - *mcETransportTrap*
     - 
   * - *rectangletrap*
     - *mcTransportRectangleTrap*
     - 
   * - *spheretrap*
     - *mcTransportSphereTrap*
     - 
   * - *axial_splitter*
     - *mcTransportAxialSymmetricSplitter*
     - 
   * - *simple_splitter*
     - *mcTransportSimpleSplitter*
     - 
   * - *slab*
     - *mcTransportSlab*
     - 
   * - *group*
     - *mcTransport*
     - Search and parse embedded modules
   * - *embedding*
     - *mcTransport*
     - 
   * - *embedded_group*
     - *mcTransportEmbeddedGroup*
     - 
   * - *linear_chain*
     - *mcTransportLinearChain*
     - 
   * - *ptlasvegas*
     - *mcPTLasVegas*
     - Las Vegas phantom for testing portal imaging systems in radiation therapy units
   * - *esphere*
     - *mcETransportSphere*
     - 
   * - *e_convex_polygon_circle*
     - *mcETransportConvexPolygonCircle*
     - 

.. code-block:: CPP

	if (_wcsicmp(geomType.c_str(), L"ptlasvegas") == 0)
	{
		t = new mcPTLasVegas(origin, normal, xaxis);
	}
	else if (_wcsicmp(geomType.c_str(), L"esphere") == 0)
	{
		t = new mcETransportSphere(origin, normal, xaxis, radius);
	}
	else if (_wcsicmp(geomType.c_str(), L"e_convex_polygon_circle") == 0)
	{
		t = new mcETransportConvexPolygonCircle(origin, normal, xaxis, poly_z, poly_x);
	}


Module "*cylinder*"
-------------------

.. list-table::
   :name: _module_cylinder-table
   :widths: 3, 15
   :width: 100%
   :header-rows: 0

   * - **Class:**
     - *mcTransportCylinder*
   * - **Description:**
     - Cylinder object
   * - **Example:**
     - 

Module "*rectanglepolygonsidehole*"
-----------------------------------

.. list-table::
   :name: _module_rectanglepolygonsidehole-table
   :widths: 3, 15
   :width: 100%
   :header-rows: 0

   * - **Class:**
     - *mcTransportRectanglePolygonSideHole*
   * - **Description:**
     - Cylindrically symmetric object with shaped by polygon hole 
   * - **Example:**
     - 

Module "*group*"
----------------

.. list-table::
   :name: _module_group-table
   :widths: 3, 15
   :width: 100%
   :header-rows: 0

   * - **Class:**
     - *mcTransport*
   * - **Description:**
     - Not a real object. It is a special construction to group other object.
       Group is treated by shower as a single object. 
       I.e. it looks like a module in z-ordered chain of modules.  
   * - **Example:**
     - .. code-block:: XML

            <!--CyberKnife primary collimator with aluminum filter layer-->
            <module type="group" name="PRI_1">
              <position unit="cm" x="0" y="0" z="-78.5" />
              <normal x="0" y="0" z="1" />
              <xaxis x="1" y="0" z="0" />
              <module type="ring" name="PRI_1W" medium="W700ICRU" density="1">
                <Color r="0" g="0.5" b="0" t="0.5" />
                <position unit="cm" x="0" y="0" z="0" />
                <normal x="0" y="0" z="1" />
                <xaxis x="1" y="0" z="0" />
                <size unit="cm" r0="0.24" r1="8.25" height="3.1"/>
              </module>
              <module type="cylinder" name="PRI_1A" medium="AL700ICRU" density="1">
                <Color r="1" g="1" b="1" t="0.5" />
                <position unit="cm" x="0" y="0" z="0" />
                <normal x="0" y="0" z="1" />
                <xaxis x="1" y="0" z="0" />
                <size unit="cm" radius="0.24" height="0.85"/>
              </module>
            </module>

Next is a parser code fragment, that demonstrates how content of *"group"* module description is interpreted.
Grouping is implemented through the concept of *Regions*.

.. code-block:: CPP

	else if (_wcsicmp(geomType.c_str(), L"group") == 0)
	{
		t = new mcTransport(origin, normal, xaxis);
		// Search and parse embedded modules
		for (auto node : geometry.Nodes)
		{
			if (_wcsicmp(node.Name.c_str(), L"module") == 0)
			{
				mcTransport* region = GeometryParser::ParseTransport(node, media, nThreads);
				region->setDefaultScore(nThreads);
				t->addRegion(region);
			}
		}
	}

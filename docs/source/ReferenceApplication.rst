.. _ch_reference_application:

Reference application
=====================

Transport modules coding in XML geometry file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table:: List of simple transport module types
   :name: _module_simple_types-table
   :widths: 4, 6, 6
   :width: 100%
   :header-rows: 1

   * - Type name
     - Transport class
     - Description
   * - *cylinder*
     - *mcTransportCylinder*
     - Solid cylinder
   * - *cone*
     - *mcTransportCone*
     - Solid cone
   * - *prism*
     - *mcTransportPrism*
     - Solid prism
   * - *wedge*
     - *mcTransportWedge*
     - Solid rectangle cylinder
   * - *ring*
     - *mcTransportRing*
     - Solid ring
   * - *conicalring*
     - *mcTransportConicalRing*
     - Solid ring with focused edges
   * - *conicalhole*
     - *mcTransportConicalHole*
     - Solid ring with hole and all edges focused
   * - *rectanglering*
     - *mcTransportRectangleRing*
     - Solid rectangle ring
   * - *jaw_pair*
     - *mcTransportJawPairRounded*
     - Pair of solid prism, representing collimator jaws
   * - *jaw_pair_focused*
     - *mcTransportJawPairFocused*
     - Pair of solid prism with focused internal edges, representing collimator jaws
   * - *rectanglepolygonsidehole*
     - *mcTransportRectanglePolygonSideHole*
     - Solid prism with rectangle hole which internal edge is defined by arbitrary polygon
   * - *mlc*
     - *mcTransportMLC*
     - *MLC* modelling module
   * - *gridcylinder*
     - *mcTransportCylinderStack*
     - Stack of nested cylinders with space between them
   * - *planefilter*
     - *mcTransportPlaneFilter*
     - Transparent infinite plane used for hosting objects like scoring
   * - *etrap*
     - *mcETransportTrap*
     - An other infinite plane, which kills particles after processing in scoring
   * - *rectangletrap*
     - *mcTransportRectangleTrap*
     - Rectangular trap, which kills particles except those that go through the internal hole or onside the transport
   * - *spheretrap*
     - *mcTransportSphereTrap*
     - Spherical trap that kills particle after scoring
   * - *axial_splitter*
     - *mcTransportAxialSymmetricSplitter*
     - Splits particles in cylindrical symmetric cases with eventual rotation 
       around *Z* axis and preserving total weight
   * - *simple_splitter*
     - *mcTransportSimpleSplitter*
     - Generates copies of particles while keeping total weight
   * - *slab*
     - *mcTransportSlab*
     - Simple infinite slab
   * - *ptlasvegas*
     - *mcPTLasVegas*
     - Las Vegas phantom for testing portal imaging systems in radiation therapy units
   * - *esphere*
     - *mcETransportSphere*
     - Solid sphere
   * - *e_convex_polygon_circle*
     - *mcETransportConvexPolygonCircle*
     - Cylindrically symmetric solid object formed by polygon rotation around *Z* axis 

Next is an example of code, which parses geometry configuration file and creates network of transport objects.

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

Composite transport modules coding in XML geometry file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table:: List of composite transport module types
   :name: _module_composite_types-table
   :widths: 4, 6, 6
   :width: 100%
   :header-rows: 1

   * - Type name
     - Transport class
     - Description
   * - *group*
     - *mcTransport*
     - Embedding other transport as regions
   * - *embedding*
     - *mcTransport*
     - Embedding nested transport modules
   * - *embedded_group*
     - *mcTransportEmbeddedGroup*
     - Group of transport modules with arbitrary order
   * - *linear_chain*
     - *mcTransportLinearChain*
     - *Z* ordered chain of transports that pass particles to neighbors


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

In this code example variable *t* represents transport module which will be inserted in to the current transport chain.
Nested modules will be simulated as regular modules, but under the multiregion transport *t* roof.

.. note:: 
   Nested modules coordinate system is defined relatively to the parent module.

In this particular module *"group"* example two modules are nested to parent module.
Their registration in parent as regions means, that after each event of exiting 
from nested module particles will be tested on possibility of hitting any other region.
If so, transport will be delegated to the nearest hitting region.
Otherwise particle will exit parent module.

Module "*embedding*"
--------------------

.. list-table::
   :name: _module_embedding-table
   :widths: 3, 15
   :width: 100%
   :header-rows: 0

   * - **Class:**
     - Can be any transport class, that implements nested transport
   * - **Description:**
     - Nested set of transport modules. The most external one represents 
       the whole module from the point of view of consuming linear chain of transports. 
   * - **Example:**
     - .. code-block:: XML

        <?xml version="1.0" encoding="utf-8"?>
        <!-- Faddegon experiment geometry with spherical detector -->
        <accelerator>
        <!--
            Objects are sequence of nested objects.
            First one is experimental unit as whole (wich can be inserted in other chains).
            All its childs are nested as a Russian dole.
            First child is most internal object (target in this case).
            All child coordinates are defined relatively to the parent.
        -->
        <module type="embedding" name="ExperimentalUnit">
            <position unit="cm" x="0" y="0" z="0" />
            <normal x="0" y="0" z="1" />
            <xaxis x="1" y="0" z="0" />
            <!--Target-->
            <!-- Pb (9.13 g/cm2 thick, 17.95 g/cm2 radius RHO= 1.1340E+01 -->
            <module type="cylinder" name="Target" medium="PB700ICRU" density="1">
            <Color r="1" g="1" b="1" t="0" />
            <position unit="cm" x="0" y="0" z="0" />
            <normal x="0" y="0" z="1" />
            <xaxis x="1" y="0" z="0" />
            <size unit="cm" radius="1.58" height="0.805"/>
            </module>

            <!-- Assume transport in Vacuum (i.e. next will be detector) -->

            <!--Spherical detector is requred to atach scoring-->
            <module type="spheretrap" name="Detector" medium="H2O700ICRU" density="1">
            <Color r="1.0" g="0.5" b="0.5" t="0.8" />
            <position unit="cm" x="0" y="0" z="0" />
            <normal x="0" y="0" z="1" />
            <xaxis x="1" y="0" z="0" />
            <size unit="cm" radius="100.0"/>
            </module>
        </module>
        </accelerator>

*XML* configuration file contains nested transports of *"embedding"* module in the order of internal transports first.

Next is a parser code fragment, that demonstrates how content of *"embedding"* module description is interpreted.
Grouping is implemented through the concept of nested transports,
which means *russian doll* type construction, i.e. nested transport is completely covered by parent transport.
Each transport may have only one reference to embedded transport.

.. code-block:: CPP

	else if (_wcsicmp(geomType.c_str(), L"embedding") == 0)
	{
		// Temporal transport needed for group coordinate system
		mcTransport ttmp(origin, normal, xaxis);
		mcTransport* eprev = nullptr;

		// Search and parse nested transports
		for (auto node : geometry.Nodes)
		{
			if (_wcsicmp(node.Name.c_str(), L"module") == 0)
			{
				mcTransport* eobj = GeometryParser::ParseTransport(node, media, nThreads);
				if (eobj == nullptr)
					throw exception("Embedding group should contain only embedding objects");

				// Transfer coordinate system from group to world
				eobj->MoveToCoordinateSystem(ttmp.MT2W());

				eobj->setDefaultScore(nThreads);
				if (eprev != nullptr)
				{
					eprev->setExternalTransport(eobj);
					eobj->setInternalTransport(eprev);
				}
				eprev = eobj;
			}
		}
		t = eprev;
		skipInit = true;
	}

In this code example variable *t* represents transport module 
which will be inserted in to the current transport chain.
Nested modules will be simulated by moving particles between internal / external 
transport until particle leave the most external transport.

.. note:: 
   Nested modules coordinate system is defined relatively to the parent module.

In this particular module *"embedding"* example there are two modules.
Оne (*Target*) is inside the other (*spherical detector*).
In contrast to linear chain of transports ordered by increasing *Z* position
nested transport is implemented by setting references in each transport
*"embedding"* module references to external and internal ones.
Registration transports in a such a way is provided by 
*setExternalTransport() / setInternalTransport()* functions.

Module "*embedded_group*"
-------------------------

.. list-table::
   :name: _module_embedded_group-table
   :widths: 3, 15
   :width: 100%
   :header-rows: 0

   * - **Class:**
     - *mcTransportEmbeddedGroup*
   * - **Description:**
     - Special class with overridden transport, which at transport inside itself 
       detects nearest transport from internal collection which can be hit by particles and pass transport to it.
       If no internal transport can be hit, particle exits from this module. 
   * - **Example:**
     - .. code-block:: XML

        <!--
        Group of modules consisting of target and 
        shielding elements and primary collimator hole
        -->
        <module type="embedded_group" name="ShieldInternalGroup">
          <position unit="cm" x="0" y="0" z="0" />
          <normal x="0" y="0" z="1" />
          <xaxis x="1" y="0" z="0" />

          <!--Accelerator part-->
          <module type="embedding" name="Accelerator">
            <position unit="cm" x="0" y="0" z="0" />
            <normal x="0" y="0" z="1" />
            <xaxis x="1" y="0" z="0" />

            <!--......................-->

          </module>

          <!--Primary collimator hole-->
          <module type="e_convex_polygon_circle" name="PriCollAir" 
                  medium="AIR700ICRU" density="1">
            <Color r="0.2" g="0.2" b="1.0" t="0.8" />
            <position unit="cm" x="0" y="0" z="0" />
            <normal x="0" y="0" z="1" />
            <xaxis x="1" y="0" z="0" />
            <polygon>
              <point z="1.801" x="0.65"/>
              <point z="9.799" x="2.6"/>
            </polygon>
          </module>
        </module>

*XML* configuration file contains nested transports of *"embedded_group"* module in the order of internal transports first.

Next is a parser code fragment, that demonstrates how content of *"embedded_group"* module description is interpreted.
Grouping is implemented through the concept of internal collection of transport modules, 
which can be hit by the particles in any order.

.. code-block:: CPP

	else if (_wcsicmp(geomType.c_str(), L"embedded_group") == 0)
	{
		t = new mcTransportEmbeddedGroup(origin, normal, xaxis);

		// Search and parse nested transports
		for (auto node : geometry.Nodes)
		{
			if (_wcsicmp(node.Name.c_str(), L"module") == 0)
			{
				mcTransport* eobj = GeometryParser::ParseTransport(node, media, nThreads);
				if (eobj == nullptr)
					throw exception("embedded_group internal module parse error");

				// Transfer coordinate system from group to world
				eobj->MoveToCoordinateSystem(t->MT2W());
				((mcTransportEmbeddedGroup*)t)->addTransport(eobj);
			}
		}
	}

In this code example variable *t* represents transport module 
which will be inserted in to the current transport chain.

.. note:: 
   Nested modules coordinate system is defined relatively to the parent module.

In this particular module *"embedded_group"* example demonstrates 
how primary collimator hole can be implemented inside the whole shield 
side by side with electron accelerator inside the same shied object.

Module "*linear_chain*"
-----------------------

.. list-table::
   :name: _module_linear_chain-table
   :widths: 3, 15
   :width: 100%
   :header-rows: 0

   * - **Class:**
     - *mcTransportLinearChain*
   * - **Description:**
     - Special class with overridden transport, which at transport inside itself 
       detects nearest transport from internal collection which can be hit by particles and pass transport to it.
       If no internal transport can be hit, particle exits from this module. 
   * - **Example:**
     - .. code-block:: XML

        <module type="linear_chain" name="RadiationHead">
            <position unit="cm" x="0" y="0" z="0" />
            <normal x="0" y="0" z="1" />
            <xaxis x="1" y="0" z="0" />

            <module type="embedding" name="ShieldGroup">
                <position unit="cm" x="0" y="0" z="0" />
                <normal x="0" y="0" z="1" />
                <xaxis x="1" y="0" z="0" />

                <!--......................-->

            </module>

            <!--Flattening filter-->
            <module type="cone" name="Flattenig Filter" medium="W700ICRU" density="1">
                <Color r="0.7" g="0.7" b="0.7" t="0.2" />
                <position unit="cm" x="0" y="0" z="11" />
                <normal x="0" y="0" z="-1" />
                <xaxis x="-1" y="0" z="0" />
                <size unit="cm" radius="2.4" height="1.0"/>
            </module>

            <!--Ionization chamber-->
            <module type="cylinder" name="IonChamber" medium="AIR700ICRU" density="1">
                <Color r="1.0" g="0.5" b="0.0" t="0.7" />
                <position unit="cm" x="0" y="0" z="11.5" />
                <normal x="0" y="0" z="1" />
                <xaxis x="1" y="0" z="0" />
                <size unit="cm" radius="6.0" height="2"/>
            </module>
        
            <!--......................-->

        </module>

*XML* configuration file contains nested transports of *"linear_chain"* module in the order of internal transports first.

Next is a parser code fragment, that demonstrates how content of *"linear_chain"* module description is interpreted.
Grouping is implemented through the concept of internal linear Z ordered chain of transport modules.

.. code-block:: CPP

	else if (_wcsicmp(geomType.c_str(), L"linear_chain") == 0)
	{
		t = new mcTransportLinearChain(origin, normal, xaxis);

		for (auto node : geometry.Nodes)
		{
			if (_wcsicmp(node.Name.c_str(), L"module") == 0)
			{
				mcTransport* eobj = GeometryParser::ParseTransport(node, media, nThreads);
				if (eobj == nullptr)
					throw exception("Embedding group should contain only embedding objects");
				eobj->MoveToCoordinateSystem(t->MT2W());
				((mcTransportLinearChain*)t)->addTransport(eobj);
			}
		}
		((mcTransportLinearChain*)t)->completeInit();
	}

In this code example variable *t* represents transport module 
which will be inserted in to the current transport chain.

.. note:: 
   Nested modules coordinate system is defined relatively to the parent module.

In this particular module *"linear_chain"* example demonstrates 
implementation accelerator head collimating structures sequence,
where elements can be separated by planes, orthogonal to *Z* axis.

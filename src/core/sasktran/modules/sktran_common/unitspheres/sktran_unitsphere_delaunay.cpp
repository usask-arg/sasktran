#include "../sktran_common.h"
#include <fstream>

SKTRAN_UnitSphere_Delaunay::SKTRAN_UnitSphere_Delaunay()
{
    m_hasDummyPoint = false;
    m_faces         = NULL;
    m_neigs         = NULL;
}

SKTRAN_UnitSphere_Delaunay::~SKTRAN_UnitSphere_Delaunay()
{
    ReleaseResources();
}

void SKTRAN_UnitSphere_Delaunay::ReleaseResources( )
{
    delete[] m_faces; m_faces=NULL;
    delete[] m_neigs; m_neigs=NULL;

    return;
}

bool SKTRAN_UnitSphere_Delaunay::CreateTriangulation( const nxVector* unitVecs_xyz, size_t numTriplets, const nxVector* openAxis )
{
    bool   ok              = true;
//    size_t numWorkingVerts = 0;

    // Need to check whether SKTRAN_UnitSphere_V2 release needs to be called
    ReleaseResources();

    ok = ok && (4<numTriplets) || (3<numTriplets && NULL!=openAxis);
    if(!ok) nxLog::Record(NXLOG_ERROR,"SKTRAN_UnitSphere_Delaunay::CreateTriangulation, Cannot use %d vertices to construct %s triangulation.", numTriplets, NULL==openAxis?"closed":"open");
	ok = ok && CopyVerticesToInternal( unitVecs_xyz, numTriplets, openAxis );
    if(!ok) nxLog::Record(NXLOG_ERROR,"SKTRAN_UnitSphere_Delaunay::CreateTriangulation, Cannot store geometry.");
	ok = ok && ConstructTriangulation( );
    if(!ok) nxLog::Record(NXLOG_ERROR,"SKTRAN_UnitSphere_Delaunay::CreateTriangulation, Cannot create triangulation.");
	ok = ok && ConstructLookupObjects( );
    if(!ok) nxLog::Record(NXLOG_ERROR,"SKTRAN_UnitSphere_Delaunay::CreateTriangulation, Cannot create lookup objects.");
	
    if(!ok) nxLog::Record(NXLOG_ERROR,"SKTRAN_UnitSphere_Delaunay::CreateTriangulation, Cannot remove end vector, this means something went wrong earlier.");
	
	return ok;
}

bool SKTRAN_UnitSphere_Delaunay::CopyVerticesToInternal( const nxVector* unitVecs, size_t numTriplets, const nxVector* openAxis )
{
    bool ok = true;
 //   double* inputPtr      = NULL;
 //   nxVector* internalPtr = NULL;

    m_hasDummyPoint       = NULL!=openAxis;
    if(m_hasDummyPoint) m_openAxis = *openAxis;
    
    ok = ok && AllocateVertices( m_hasDummyPoint ? numTriplets+1 : numTriplets );

    if(ok) for(size_t vidx=0; vidx<numTriplets; vidx++) UnitVectorAtVar(vidx) = unitVecs[vidx];
    if(m_hasDummyPoint) ok = ok && CreateDummyPoint();
	
    return ok;
}

bool SKTRAN_UnitSphere_Delaunay::CreateDummyPoint( )
{
    bool ok = true;

    UnitVectorAtVar(NumUnitVectors()-1) = m_openAxis;

	return ok;
}

/* Returns true if tetrahedron containing the origin was constructed */
bool SKTRAN_UnitSphere_Delaunay::ConstructSimplex( tuple4<size_t>& simplexVerts )
{
    bool ok = true;
    
    ok = ok && ConstructSimplex_bruteForce( simplexVerts );
    return ok;
}

bool SKTRAN_UnitSphere_Delaunay::ConstructSimplex_bruteForce( tuple4<size_t>& simplexVerts )
{
    bool ok                = true;
	bool originIsContained = false;
	bool originIsWellContained = false;
	bool keepTryingToContainOrigin = true;
	bool originCannotBeWellContained = false;

    size_t   n0, n1, n2, n3;
    nxVector v0, v1, v2, v3;
    double   d0, d1, d2, d3;

    // Choose the first set of vertices that contain the origin
	while(keepTryingToContainOrigin){
		for(n0=NumUnitVectors(); keepTryingToContainOrigin && 0<n0--;)
			for(n1=n0; keepTryingToContainOrigin && 0<n1--;)
				for(n2=n1; keepTryingToContainOrigin && 0<n2--;)
					for(n3=n2; keepTryingToContainOrigin && 0<n3--;)
					{
						v0 = UnitVectorAt(n0); v1=UnitVectorAt(n1); v2=UnitVectorAt(n2); v3=UnitVectorAt(n3);
						d0 = -( v1.X()*(v2.Y()*v3.Z()-v2.Z()*v3.Y()) - v1.Y()*(v2.X()*v3.Z()-v2.Z()*v3.X()) + v1.Z()*(v2.X()*v3.Y()-v2.Y()*v3.X()) );   // These determinants give a sort of "signed volume"
						d1 =  ( v0.X()*(v2.Y()*v3.Z()-v2.Z()*v3.Y()) - v0.Y()*(v2.X()*v3.Z()-v2.Z()*v3.X()) + v0.Z()*(v2.X()*v3.Y()-v2.Y()*v3.X()) );   // of the tetrahedron formed between each face and
						d2 = -( v0.X()*(v1.Y()*v3.Z()-v1.Z()*v3.Y()) - v0.Y()*(v1.X()*v3.Z()-v1.Z()*v3.X()) + v0.Z()*(v1.X()*v3.Y()-v1.Y()*v3.X()) );   // the origin. If they all have the same sign as 
						d3 =  ( v0.X()*(v1.Y()*v2.Z()-v1.Z()*v2.Y()) - v0.Y()*(v1.X()*v2.Z()-v1.Z()*v2.X()) + v0.Z()*(v1.X()*v2.Y()-v1.Y()*v2.X()) );   // d0+d1+d2+d3, the origin is contained within v0..v3
                   
						originIsContained = (d0>0 && d1>0 && d2>0 && d3>0) || (d0<0 && d1<0 && d2<0 && d3<0);   // This is equivalent to testing against sgn(d0+d1+d2+d3)
						originIsWellContained = originIsContained && (1e-10<abs(d0) && 1e-10<abs(d1) && 1e-10<abs(d2) && 1e-10<abs(d3));
						keepTryingToContainOrigin = !(originIsWellContained || (originCannotBeWellContained&&originIsContained));
					}

		if( !originIsWellContained ){
			if( !originCannotBeWellContained ){
				// Tried once with "tight" containment criteria; loosen that and try again.
				originCannotBeWellContained=true;
				nxLog::Record(NXLOG_WARNING, "SKTRAN_UnitSphere_Delaunay::ConstructSimplex_bruteForce, The origin should be ``well within'' at least one 4-tuple of points.");
			} else if( !originIsContained ){
				// Couldn't contain the origin at all
				keepTryingToContainOrigin=false;
				nxLog::Record(NXLOG_WARNING, "SKTRAN_UnitSphere_Delaunay::ConstructSimplex_bruteForce, The origin MUST be at least ``barely within'' at least one 4-tuple of points.");
			} else{
				keepTryingToContainOrigin=false;
				nxLog::Record(NXLOG_WARNING, "SKTRAN_UnitSphere_Delaunay::ConstructSimplex_bruteForce, Unexpected logical path followed -- this shouldn't happen.");
			}
		}
	}

    ok = ok && originIsContained;

    if(ok)
	{
        // Create faces from simplex vertices
        m_faces[0].d[0]=n0; m_faces[0].d[1]=n1; m_faces[0].d[2]=n2;
        m_faces[1].d[0]=n0; m_faces[1].d[1]=n1; m_faces[1].d[2]=n3;
        m_faces[2].d[0]=n0; m_faces[2].d[1]=n2; m_faces[2].d[2]=n3;
        m_faces[3].d[0]=n1; m_faces[3].d[1]=n2; m_faces[3].d[2]=n3;
        
        // Set neighbours for the simplex vertices
        m_neigs[0].d[0]=3;  m_neigs[0].d[1]=2;  m_neigs[0].d[2]=1;
        m_neigs[1].d[0]=3;  m_neigs[1].d[1]=2;  m_neigs[1].d[2]=0;
        m_neigs[2].d[0]=3;  m_neigs[2].d[1]=1;  m_neigs[2].d[2]=0;
        m_neigs[3].d[0]=2;  m_neigs[3].d[1]=1;  m_neigs[3].d[2]=0;

        // Tell caller which vertices are in simplex
        simplexVerts.d[0] = n0;
        simplexVerts.d[1] = n1;
        simplexVerts.d[2] = n2;
        simplexVerts.d[3] = n3;
	} else{
        nxLog::Record(NXLOG_ERROR, "SKTRAN_UnitSphere_Delaunay::ConstructSimplex_slow, Could not find simplex to contain origin -- it probably doesn't exist.");
	}

    return ok;
}


/* This function needs to be split into several subfunctions. Need to through
 * and see which code blocks can be isolated
 */
bool SKTRAN_UnitSphere_Delaunay::ConstructTriangulation( )
{
    bool ok = true;

    const size_t numExpectedFaces = 2*(NumUnitVectors() - 2);
    const size_t badFlagVal = numExpectedFaces+1;

    m_faces = new tuple3<size_t>[numExpectedFaces];
    m_neigs = new tuple3<size_t>[numExpectedFaces];
    ok = ok && NULL!=m_faces && NULL!=m_neigs;

    tuple2<size_t>* uncreatedNeigTracker = new tuple2<size_t>[numExpectedFaces];
	tuple2<size_t>* edges = new tuple2<size_t>[3*numExpectedFaces];           // This is actually a fairly arbitrary guess at how large #edges should be; it's way too large in most cases
    tuple4<double>* n4s   = new tuple4<double>[numExpectedFaces];
    size_t* pointsAssignedToFace = new size_t[numExpectedFaces*NumUnitVectors()];   // Should find a more space-efficient way to do this, currently it's n^2
	double* furthestPointDistance = new double[numExpectedFaces];
    ok = ok && NULL!=uncreatedNeigTracker && NULL!=edges && NULL!=n4s && NULL!=pointsAssignedToFace && NULL!=furthestPointDistance;
    
    size_t* vertTemp = new size_t[NumUnitVectors()];
    size_t* faceToWhichPointIsAssigned = new size_t[NumUnitVectors()];
    size_t* furthestPointIndex = new size_t[numExpectedFaces];
    size_t* numPointsAssignedToFace = new size_t[numExpectedFaces];
    size_t* newFaceIndices         = new size_t[numExpectedFaces];
    std::vector<bool> isassigned; isassigned.resize(NumUnitVectors());  // std::vector<bool> is special type that stores bool in a single bit then bit-masks; this thing is really probably a waste of space
    ok = ok && NULL!= vertTemp && NULL!= faceToWhichPointIsAssigned && NULL!=furthestPointIndex && NULL!=numPointsAssignedToFace && NULL!=newFaceIndices && NumUnitVectors()==isassigned.size();
    
	size_t numFaces = 0, numEdges = 0, temp, nim, numFacesToDelete, nned, numTempPoints, memSpId, nfidx;
    bool faceIsFound, isInsideHorizon;
    double valDiff;
    size_t fptidx, numNeedVisit, nid, prevEdgeId;
    std::vector<bool> hasBeenVisited; hasBeenVisited.resize(3*numExpectedFaces);	// Should have try/catch on this and others; ok check does nothing
    size_t* needVisitList = new size_t[numExpectedFaces];
    size_t* fromPreviousFace = new size_t[numExpectedFaces];            // Came from this face index to get to neighbour
    size_t* fromPreviousEdge = new size_t[numExpectedFaces];            // Crossed fromPreviousFace[x]'s edge fromPreviousEdge[x] to get to neighbour
    size_t* screwedUpNeighbourList = new size_t[3*numExpectedFaces];  // This is the absolute max buffer size -- I don't think it's actually possible for this to happen
    size_t* screwedUpNeighbourEdge = new size_t[3*numExpectedFaces];  // The edges for which the neighbour entry has to be changed
    size_t* facesToDelete = new size_t[numExpectedFaces];
    std::vector<bool> markedForDeletion; markedForDeletion.assign(numExpectedFaces,false);
    ok = ok && (3*numExpectedFaces)==hasBeenVisited.size() && NULL!=needVisitList && NULL!=fromPreviousFace && NULL!=fromPreviousEdge && NULL!=screwedUpNeighbourList && NULL!=screwedUpNeighbourEdge && NULL!=facesToDelete && markedForDeletion.size()==numExpectedFaces;

    tuple4<size_t> simplexVerts;

    ok = ok && ConstructSimplex( simplexVerts );  // Assume later that ConstructSimplex returns tetrahedron -- different from in MATLAB version
    if(ok)
    {
        numFaces = 4;
        isassigned[simplexVerts.d[0]] = isassigned[simplexVerts.d[1]] = isassigned[simplexVerts.d[2]] = isassigned[simplexVerts.d[3]] = true;
    } else{
        nxLog::Record(NXLOG_ERROR, "SKTRAN_UnitSphere_Delaunay::ConstructTriangulation, Could not construct simplex.");
	}

	if(ok)
	{
		// Compute hyperplane coefficients for simplex m_faces
		for(size_t fidx=0; fidx<4; ++fidx)
		{
			ok = ok && GetHyperplane(m_faces[fidx], n4s[fidx]);
			furthestPointDistance[fidx] = -1.0;                 // No vertex has been assigned to face
			numPointsAssignedToFace[fidx] = 0;
		}
    
		// assign points to m_faces
		for(size_t vidx=0; vidx<NumUnitVectors(); ++vidx)
		{
			if(!isassigned[vidx])   // Simplex vertices are already assigned
			{
				for(size_t fidx=0; fidx<numFaces; ++fidx)   // Find face
				{
					faceIsFound = TestPointUnderPlane(n4s[fidx],vidx,valDiff);
					if(faceIsFound)
					{
						// Assign face
						faceToWhichPointIsAssigned[vidx] = fidx;
						pointsAssignedToFace[fidx*NumUnitVectors() + numPointsAssignedToFace[fidx]] = vidx;
						numPointsAssignedToFace[fidx] = numPointsAssignedToFace[fidx] + 1;
						isassigned[vidx] = true;	
						if(valDiff > furthestPointDistance[fidx])
						{
							furthestPointDistance[fidx] = valDiff;
							furthestPointIndex[fidx] = vidx;
						}
						break;
					}
				}
				if(!isassigned[vidx]){
					nxLog::Record(NXLOG_WARNING, "SKTRAN_UnitSphere_Delaunay::ConstructTriangulation, Could not associate vertex with face; simplex is broken.");
					ok = false;
				}
			}
		} 


		// Perform iteration step (breadth-first)
    
		size_t itidx = 0;
		faceIsFound = true;
		while(ok && numFaces<numExpectedFaces && faceIsFound)  // Should add a check to see if we go through a full iteration without creating new m_faces (error)
		{
			++itidx;
			faceIsFound = false;    // A new face has not yet been created in this iteration through list of m_faces
			for(size_t fidx=0; ok && fidx<numFaces; ++fidx)
			{
				if(ok && 0!=numPointsAssignedToFace[fidx])        // If there are points to remove from this face
				{
					faceIsFound = true;
					fptidx = furthestPointIndex[fidx];      // Get point furthest from this face

					hasBeenVisited.assign(hasBeenVisited.size(), false);
					hasBeenVisited[3*fidx] = hasBeenVisited[3*fidx+1] = hasBeenVisited[3*fidx+2] = true;   // Don't want neighbours coming back to us for HE
					numNeedVisit = 0;                       // We don't know if we'll have to visit neighbours to find horizon edge
					numEdges     = 0;                       // No new edges have yet been found
               
					//fromPreviousEdge          = zeros(size(faces,1),1);   // I don't need to reset this
					//fromPreviousFace          = zeros(size(faces,1),1);  // same here

					// Add fid's neighbours to visitation queue
					for(size_t fidxNeigIdx=0; fidxNeigIdx<3; ++fidxNeigIdx)
					{
						temp = m_neigs[fidx].d[fidxNeigIdx];
	//                //if( (~isnan(temp)) && (0<temp) )    % If there is a neighbour here
						// Add neighbor to neighbour visit list
						needVisitList[numNeedVisit] = temp;
						// Choose this face's index in its neighbour
						//nim=1; while fidx~=m_neigs(needVisitList(numNeedVisit),nim); nim=nim+1; end
						if(fidx==m_neigs[needVisitList[numNeedVisit]].d[0]) nim = 0;
						else if(fidx==m_neigs[needVisitList[numNeedVisit]].d[1]) nim = 1;
						else nim=2;
						hasBeenVisited[3*needVisitList[numNeedVisit] + nim] = true;
						fromPreviousFace[fidxNeigIdx] = fidx;
						fromPreviousEdge[fidxNeigIdx] = fidxNeigIdx;
						//else
						//    % Add new edge to horizon
						//    numEdges = numEdges+1;
						//    edges(numEdges,:) = m_faces(fidx,(1:3)~=fidxNeigIdx);   % We crossed a boundary when we went across edge nlistum in face #1
						//    screwedUpNeighbourEdge(numEdges) = 0;  % The neighbour does not exist
						//    screwedUpNeighbourList(numEdges) = 0;  % The neighbour does not exist
						//end
						++numNeedVisit;
					}

					numFacesToDelete       = 1;
					facesToDelete[numFacesToDelete-1] = fidx;
					markedForDeletion[fidx] = true;

					// Find horizon edge
					for(size_t visidx=0; visidx<numNeedVisit; ++visidx) // Visit all neighbours that have an edge on horizon
					{
						nid = needVisitList[visidx];
                
						// Check if neighbour is inside horizon
						isInsideHorizon = TestPointUnderPlane(n4s[nid],fptidx, valDiff);
						if(isInsideHorizon)
						{
							// Face is inside horizon

							// It needs to be deleted
							if(!markedForDeletion[nid])
							{
								facesToDelete[numFacesToDelete] = nid;  // Face is inside horizon
								markedForDeletion[nid] = true;          // will be deleted
								numFacesToDelete = numFacesToDelete+1;
							}
							// Check neighbour's neighbours
							for(size_t nnid=0; nnid<3; ++nnid)
							{
								// Neighbour's neighbour's index to the neighbour
	//                            nim=1; while nid~=m_neigs(neigs(nid,nnid),nim); nim=nim+1; end
								if(nid==m_neigs[m_neigs[nid].d[nnid]].d[0]) nim = 0;
								else if(nid==m_neigs[m_neigs[nid].d[nnid]].d[1]) nim = 1;
								else nim=2;
                        
								// See about adding the neighbour's neighbours
								if(!hasBeenVisited[3*m_neigs[nid].d[nnid] + nim])          // If the edge from the neighbour to the neighbour's neighbour has not been crossed
								{
									hasBeenVisited[3*m_neigs[nid].d[nnid] + nim] = true;   // The edge has been covered from the neigbour's side
									hasBeenVisited[3*nid + nnid] = true;                 // As well as my side (the crossings are equivalent)
									needVisitList[numNeedVisit]    = m_neigs[nid].d[nnid];                // The neighbour's neighbour needs to be visited
									fromPreviousEdge[numNeedVisit] = nnid;                              // The visit comes from this edge of the neighbour currently being examined
									fromPreviousFace[numNeedVisit] = nid;                               // The visit comes from this neighbour
									++numNeedVisit; 
								}
							}
						} else{
							// Add edge to horizon
                        
							prevEdgeId = fromPreviousEdge[visidx];
                        
							switch(prevEdgeId)
							{
							case 0:
								edges[numEdges].d[0] = m_faces[fromPreviousFace[visidx]].d[1];
								edges[numEdges].d[1] = m_faces[fromPreviousFace[visidx]].d[2];
								break;
							case 1:
								edges[numEdges].d[0] = m_faces[fromPreviousFace[visidx]].d[0];
								edges[numEdges].d[1] = m_faces[fromPreviousFace[visidx]].d[2];
								break;
							case 2:
								edges[numEdges].d[0] = m_faces[fromPreviousFace[visidx]].d[0];
								edges[numEdges].d[1] = m_faces[fromPreviousFace[visidx]].d[1];
								break;
							default:
								nxLog::Record(NXLOG_ERROR, "SKTRAN_UnitSphere_Delaunay::ConstructTriangulation, Cannot construct horizon edge.");
								ok = false;
							}
                   
                    
							// Find the index of the screwed up edge in the screwed up face
							//  nim=1; while fromPreviousFace(visidx)~=m_neigs(nid,nim); nim=nim+1; end
							if(fromPreviousFace[visidx]==m_neigs[nid].d[0]) nim = 0;
							else if(fromPreviousFace[visidx]==m_neigs[nid].d[1]) nim = 1;
							else nim=2;
							screwedUpNeighbourEdge[numEdges] = nim;
							screwedUpNeighbourList[numEdges] = nid;  // Keep track of which neighbour needs to have *its* neighbour list updated; since this face is just outside the horizon, its neighbours will be screwed up
							uncreatedNeigTracker[numEdges].d[0] = uncreatedNeigTracker[numEdges].d[1] = badFlagVal;
							++numEdges;
						}
					}

					// Put edges into a cyclic order
					for(size_t ned=0; ned<numEdges; ++ned)
					{
						nned = ned;
						while(badFlagVal==uncreatedNeigTracker[ned].d[1])
						{
							nned = nned+1;
							if( edges[ned].d[0]==edges[nned].d[0] )
							{
								uncreatedNeigTracker[ned] .d[1] = nned;
								uncreatedNeigTracker[nned].d[1] = ned;
							}else if( edges[ned].d[0]==edges[nned].d[1] )
							{
								uncreatedNeigTracker[ned] .d[1] = nned; 
								uncreatedNeigTracker[nned].d[0] = ned;
							}
						}
						nned = ned;
						while(badFlagVal==uncreatedNeigTracker[ned].d[0])
						{
							nned = nned+1;
							if( edges[ned].d[1]==edges[nned].d[0] )
							{
								uncreatedNeigTracker[ned] .d[0] = nned;
								uncreatedNeigTracker[nned].d[1] = ned;
							}else if( edges[ned].d[1]==edges[nned].d[1] )
							{
								uncreatedNeigTracker[ned] .d[0] = nned;
								uncreatedNeigTracker[nned].d[0] = ned;
							}
							if( nned >= numExpectedFaces )
							{
								return false;
							}
						}
					}

					// Grab vertices from old m_faces
					numTempPoints = 0;
					for(size_t delid=0; delid<numFacesToDelete; ++delid)
					{
						for(size_t tpid=0; tpid<numPointsAssignedToFace[facesToDelete[delid]]; ++tpid)
						{
							if(fptidx != pointsAssignedToFace[facesToDelete[delid]*NumUnitVectors() + tpid] )   // Don't want to put our new vertex in the temporary points stack
							{
								vertTemp[numTempPoints] = pointsAssignedToFace[facesToDelete[delid]*NumUnitVectors() + tpid];
								++numTempPoints;
							}
						}
					}

				
					// Add new m_faces
					for(size_t nftdid=0; nftdid<numFacesToDelete; ++nftdid)
					{
						newFaceIndices[nftdid] = facesToDelete[nftdid];
						markedForDeletion[facesToDelete[nftdid]] = false;
						//facesToDelete[nftdid] = 0;    // Could clear memory, but it's unnecessary
					}
					for(int nfid=0; nfid<(int(numEdges)-int(numFacesToDelete)); ++nfid)
					{
						newFaceIndices[numFacesToDelete+nfid] = numFaces;
						++numFaces;
					}
					for(size_t ned=0; ned<numEdges; ++ned)
					{
						memSpId = newFaceIndices[ned];
						numPointsAssignedToFace[memSpId] = 0;
						//if(0!=screwedUpNeighbourList[ned])    // This line is from the MATLAB version where tetrahedra could give more than one face to CH
						m_neigs[screwedUpNeighbourList[ned]].d[screwedUpNeighbourEdge[ned]] = memSpId;    // Update neighbour's list of neighbours
						// Create this face
						m_faces[memSpId].d[0] = edges[ned].d[0];
						m_faces[memSpId].d[1] = edges[ned].d[1];
						m_faces[memSpId].d[2] = fptidx; 
						m_neigs[memSpId].d[0] = newFaceIndices[uncreatedNeigTracker[ned].d[0]];
						uncreatedNeigTracker[ned].d[0] = badFlagVal;
						m_neigs[memSpId].d[1] = newFaceIndices[uncreatedNeigTracker[ned].d[1]];
						uncreatedNeigTracker[ned].d[1] = badFlagVal;
						m_neigs[memSpId].d[2] = screwedUpNeighbourList[ned];
						// Create hyperplane for new face
						ok = ok && GetHyperplane(m_faces[memSpId], n4s[memSpId]);
						furthestPointDistance[memSpId] = -1.0;  // No vertex has been assigned to face
					}

					// Assign vertices to new m_faces
					for(size_t pidx=0; pidx<numTempPoints; ++pidx)
					{
						temp  = vertTemp[pidx];
						for(size_t nfed=0; nfed<numEdges; ++nfed)
						{
							//Search for a face
							nfidx = newFaceIndices[nfed];
	//                    %faceDist = (-(n4s(nfidx,1)*verts(vidx,1) + n4s(nfidx,2)*verts(vidx,2) + n4s(nfidx,3)*verts(vidx,3)) / n4s(nfidx,4));
							faceIsFound = TestPointUnderPlane(n4s[nfidx], temp, valDiff);
							if(faceIsFound)
							{
								// Assign face
								faceToWhichPointIsAssigned[temp] = nfidx;
								pointsAssignedToFace[nfidx*NumUnitVectors() + numPointsAssignedToFace[nfidx]] = temp;
								++numPointsAssignedToFace[nfidx];
								if(furthestPointDistance[nfidx]<valDiff)
								{
									furthestPointDistance[nfidx] = valDiff;
									furthestPointIndex[nfidx] = temp;
			            		}
								break;
							}
						}
					}

					// Clear buffer memory to zeros
					hasBeenVisited.assign(hasBeenVisited.size(),false);
					for(size_t nfed=0; nfed<numEdges; ++nfed) newFaceIndices[nfed] = 0;
				}
			}
		}

		if(numExpectedFaces!=numFaces)
		{
			ok = false;
			nxLog::Record(NXLOG_WARNING,"SKTRAN_UnitSphere_Delaunay::ConstructTriangulation, Did not get expected number of faces. Some points may be too close together or are not on unit sphere.");
		}
		m_numFaces = numFaces;
    }

	delete[] uncreatedNeigTracker; delete[] edges; delete[] n4s; delete[] pointsAssignedToFace; delete[] furthestPointDistance;
    delete[] vertTemp; delete[] faceToWhichPointIsAssigned; delete[] furthestPointIndex; delete[] numPointsAssignedToFace; delete[] newFaceIndices;
    delete[] needVisitList; delete[] fromPreviousFace; delete[] fromPreviousEdge; delete[] screwedUpNeighbourList; delete[] screwedUpNeighbourEdge; delete[] facesToDelete;

    return ok;

}

bool SKTRAN_UnitSphere_Delaunay::GetHyperplane(tuple3<size_t>& verts, tuple4<double>& target) const
{
    bool ok = true;

    nxVector v1,v2,v3;
    v1 = UnitVectorAt(verts.d[0]);  // Should make sure this isn't doing something stupid
    v2 = UnitVectorAt(verts.d[1]);
    v3 = UnitVectorAt(verts.d[2]);

	target.d[0] =  v1.Y()*(v2.Z()-v3.Z()) + v2.Y()*(v3.Z()-v1.Z()) + v3.Y()*(v1.Z()-v2.Z());
    target.d[1] =  v1.X()*(v3.Z()-v2.Z()) + v2.X()*(v1.Z()-v3.Z()) + v3.X()*(v2.Z()-v1.Z());
    target.d[2] =  v1.X()*(v2.Y()-v3.Y()) + v2.X()*(v3.Y()-v1.Y()) + v3.X()*(v1.Y()-v2.Y());
    target.d[3] =  v1.X()*(v2.Z()*v3.Y()-v2.Y()*v3.Z()) + 
                    v2.X()*(v1.Y()*v3.Z()-v1.Z()*v3.Y()) +
                    v3.X()*(v1.Z()*v2.Y()-v1.Y()*v2.Z()) ;

	return ok;
}

bool SKTRAN_UnitSphere_Delaunay::TestPointUnderPlane(tuple4<double> n4, size_t vertIndex, double& distanceUnder) const
{
    distanceUnder = -1*((n4.d[0]*UnitVectorAt(vertIndex).X() + n4.d[1]*UnitVectorAt(vertIndex).Y() + n4.d[2]*UnitVectorAt(vertIndex).Z()) / n4.d[3]) - 1.0;

    return -1E-14 < distanceUnder;
}

bool SKTRAN_UnitSphere_Delaunay::Triangulate( const nxVector& unit, size_t* unit_indexptr, double* unit_weightptr, size_t maxvertices) const
{
    bool ok = true;
	size_t dummy=0;

	ok = ok && FindThreeNearestVertices(unit, unit_indexptr, maxvertices, dummy);
	ok = ok && FindInterpolationWeights(unit, unit_indexptr, unit_weightptr, dummy);
    return ok;
}


void SKTRAN_UnitSphere_Delaunay::PrintTriangulation(std::string s) const
{
	std::ofstream f;
	f.open(s.c_str());

	for(size_t n=0; n<m_numFaces; ++n)
	{
		f << m_faces[n].d[0] << "\t" <<  m_faces[n].d[1] << "\t" << m_faces[n].d[2] << std::endl;
	}
	f.close();

}


bool SKTRAN_UnitSphere_Delaunay::Triangulate( const nxVector& unit, size_t* unit_indexptr, double* unit_weightptr, size_t maxvertices, size_t& speedHelper) const
{
	bool ok = true;
	ok = ok && FindThreeNearestVertices(unit, unit_indexptr, maxvertices, speedHelper);
	ok = ok && FindInterpolationWeights(unit, unit_indexptr, unit_weightptr, speedHelper);
    return ok;

}

bool SKTRAN_UnitSphere_Delaunay::InterpolateTriangle( const nxVector& P_sphere, size_t* unit_indexptr, double* unit_weightptr, size_t faceidx ) const
{
	nxVector	ab;
	nxVector	ac;
	nxVector	ap;
	nxVector	xunit;
    nxVector	yunit;
	nxVector	A;
	nxVector	B;
	nxVector	C;
	nxVector    P;
	double		xb;
	double		xc;
	double		yc;
	double		xp;
	double		yp;
	double		u;
	double		t;

	A = UnitVectorAt(unit_indexptr[0]);
	B = UnitVectorAt(unit_indexptr[1]);
	C = UnitVectorAt(unit_indexptr[2]);
	P = (m_faceNormals[faceidx].Dot(A)/m_faceNormals[faceidx].Dot(P_sphere)) * P_sphere;				// Query point needs to be projected onto the plane along its radius -- projecting normal to the face causes negative weights

	ab  = B-A;
	ac  = C-A;
	ap  = P-A;													// Vector from point A to point P. Note P is above plane of triangle but analysis is still good

	xunit = ab.UnitVector();									// x is parallel to AB
    yunit = ac.ComponentPerpendicularTo(xunit).UnitVector();	// y is perpendicular to AB

	xb = ab & xunit;											// X Coordinate of B (by definition B is along X axis
	xc = ac & xunit;											// X coordinate of C
//  yb = 0.0;													// Y coordinate of B is implicitly 0.0
	yc = ac & yunit;											// Y coordinate of C
	xp = ap & xunit;											// X coordinate of P
	yp = ap & yunit;											// Y coordinate of P
	u  = yp/yc;
	t  = (xp - u*xc)/xb;

	// Point is inside triangle if u and t are between 0 and 1

	unit_weightptr[0] = 1.0 - u - t;
	unit_weightptr[1] = t;
	unit_weightptr[2] = u;

#if defined(NXDEBUG)
	nxVector	PTest;
	nxVector	diff;
	nxVector	zunit;
	nxVector	aps;

	zunit  = xunit ^ yunit;
	aps    = ap.ComponentPerpendicularTo(zunit);				// Get the point P in the plane of the triangle
	PTest  = t*ab + u*ac;
	diff   = (PTest-aps);
	NXASSERT(( diff.Magnitude() < 0.000001));
#endif

	return true;
};

/*
bool SKTRAN_UnitSphere_Delaunay_InterpolateTriangle_Spherical
{
	double dividefactor;
	double totalweight = 0;
	std::vector<double> alpha;
	std::vector<double> theta;
	nxVector temp1;
	nxVector temp2;

	alpha.resize(3);
	theta.resize(3);
	for( int weightidx = 0; weightidx < 3; weightidx++ )
	{
		theta[weightidx] = P_sphere.AngleTo( UnitVectorAt( unit_indexptr[weightidx] ) ) * nxmath::Pi/180;
		temp1 = P_sphere.Cross( UnitVectorAt( unit_indexptr[weightidx] ) );
		temp2 = P_sphere.Cross( UnitVectorAt( unit_indexptr[(weightidx+1) %3] ) );

		alpha[weightidx] = temp1.AngleTo( temp2 ) * nxmath::Pi/180;

	}

	dividefactor = 0;
	for( size_t vertexidx = 0; vertexidx < 3; vertexidx++ )
	{
		dividefactor += 1.0/(tan(theta[vertexidx])) * (tan(alpha[vertexidx == 0 ? 2 : vertexidx-1]/2) + tan(alpha[vertexidx]/2));
	}
	for( size_t weightidx = 0; weightidx < 3; weightidx++ )
	{
		unit_weightptr[weightidx] = (tan(alpha[weightidx == 0 ? 2 : weightidx-1]/2) + tan(alpha[weightidx]/2)) / sin(theta[weightidx]);
		unit_weightptr[weightidx] /= dividefactor;
		if( !( unit_weightptr[weightidx] > -1 ) )
		{
			// when theta is 0 we are directly on a vertex
			unit_weightptr[weightidx] = 1.0;
		}
		totalweight += unit_weightptr[weightidx];
	}

	nxVector check;
	check.SetCoords(0,0,0);
	for( size_t weightidx = 0; weightidx < 3; weightidx++ )
	{
		check += UnitVectorAt( unit_indexptr[weightidx] ) * unit_weightptr[weightidx];
	}
	if( (check - P_sphere).Magnitude() > 1E-5 )
	{
		printf("Bad Interp\n");
	}

	for( size_t weightidx = 0; weightidx < 3; weightidx++ )
	{
		// spherical barycentric coordinates are not necessary normalized
		unit_weightptr[weightidx] /= totalweight;
	}


	return true;
}
*/

/********************
 * Non-tabled lookup
 ********************/
SKTRAN_UnitSphere_Delaunay_nonTabledLookup::SKTRAN_UnitSphere_Delaunay_nonTabledLookup()
{
	m_startFace = 0;
}

SKTRAN_UnitSphere_Delaunay_nonTabledLookup::~SKTRAN_UnitSphere_Delaunay_nonTabledLookup()
{
	ReleaseResources();
}

void SKTRAN_UnitSphere_Delaunay_nonTabledLookup::ReleaseResources()
{
	std::vector< tuple3<nxVector> > nodata;

	m_edgenorms.clear();
	m_edgenorms.swap(nodata);	// Force reallocation in case user wants to clear memory (but this function is private, so this is useless)
}

bool SKTRAN_UnitSphere_Delaunay_nonTabledLookup::ConstructLookupObjects()
{
	bool ok = true;

	double tripProd;
	nxVector a, b, c;

	m_edgenorms.resize(m_numFaces);	// Should catch
	m_faceNormals.resize(m_numFaces);
	for(size_t fidx=0; fidx<m_numFaces; ++fidx)
	{
        a = UnitVectorAt(m_faces[fidx].d[0]); b = UnitVectorAt(m_faces[fidx].d[1]); c = UnitVectorAt(m_faces[fidx].d[2]);
		tripProd = a.Dot(b.Cross(c));

		ok = ok && 1e-12 < fabs(tripProd);
		if(ok)
			{
			// a.(bxc), b.(cxa), c.(axb) are positive
			m_edgenorms[fidx].d[0] = b.Cross(c);// m_edgenorms[fidx].d[0] *= (1.0/m_edgenorms[fidx].d[0].Magnitude());
			m_edgenorms[fidx].d[1] = c.Cross(a);// m_edgenorms[fidx].d[1] *= (1.0/m_edgenorms[fidx].d[1].Magnitude());
			m_edgenorms[fidx].d[2] = a.Cross(b);// m_edgenorms[fidx].d[2] *= (1.0/m_edgenorms[fidx].d[2].Magnitude());
			m_faceNormals[fidx] = (b-a).Cross(c-a);	// normal of ABC goes outwards
			m_faceNormals[fidx] *= 1.0/m_faceNormals[fidx].Magnitude();
		
			if( tripProd < 0.0 )
			{
				// a.(cxb), b.(axc), c.(bxa) are positive
				m_edgenorms[fidx].d[0] *= -1.0;
				m_edgenorms[fidx].d[1] *= -1.0;
				m_edgenorms[fidx].d[2] *= -1.0;
				m_faceNormals[fidx]    *= -1.0;	// normal of ACB goes outwards
			}
		} else{
			nxLog::Record(NXLOG_ERROR, "SKTRAN_UnitSphere_Delaunay_nonTabledLookup::ConstructLookupObjects, Unit vectors are too close together or coplanar -- this should've been caught in triangulation stage.");
		}
	}

	return ok;
}


bool SKTRAN_UnitSphere_Delaunay_nonTabledLookup::FindThreeNearestVertices( const nxVector& unit, size_t* unit_indexptr, size_t maxvertices, size_t& faceIndex) const
{
	return FindThreeNearestVertices_directed(unit, unit_indexptr, maxvertices, faceIndex);
}


bool SKTRAN_UnitSphere_Delaunay_nonTabledLookup::FindThreeNearestVertices_directed( const nxVector& unit, size_t* unit_indexptr, size_t maxvertices, size_t& startfidx) const
{
	bool ok = true;
//	bool keepSearching = true;

	size_t maxAttempts=2*m_numFaces;
	size_t fidx = startfidx;
	size_t tryidx = 0;

	ok = ok && 2<maxvertices;
	if(ok)
	{
		for(; tryidx<maxAttempts; ++tryidx)
		{
			// If query point is on wrong side of edge i, see if point is in neighbour
			// on other side of edge i
			if(      m_edgenorms[fidx].d[0].Dot(unit) < -1E-10) fidx = m_neigs[fidx].d[0];
			else if( m_edgenorms[fidx].d[1].Dot(unit) < -1E-10) fidx = m_neigs[fidx].d[1];
			else if( m_edgenorms[fidx].d[2].Dot(unit) < -1E-10) fidx = m_neigs[fidx].d[2];
			else{
				// Point is inside this face
				*(unit_indexptr+0)  = m_faces[fidx].d[0];	// Return the indices of unit vectors defining this face
				*(unit_indexptr+1)  = m_faces[fidx].d[1];
				*(unit_indexptr+2)  = m_faces[fidx].d[2];					
				startfidx = fidx;	// Try this face first next time
				break;
			}
		}
	}

	ok = ok && tryidx<maxAttempts;
	if( !ok )
	{
		ok = FindThreeNearestVertices_bruteforce( unit, unit_indexptr, maxvertices );
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_UnitSphere_Delaunay_nonTabledLookup::FindNearestVertex		2015-02-09*/
/** Called after FindThreeNearestVertices_directed fails, returns back a single
 *  vertex in which the point is closest too
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_UnitSphere_Delaunay_nonTabledLookup::FindThreeNearestVertices_bruteforce( const nxVector& unit, size_t* unit_indexptr, size_t maxvertices ) const
{
	std::array<double,3> dotprods;
	std::fill( std::begin(dotprods), std::end(dotprods), 0.0 );
	std::fill( unit_indexptr, unit_indexptr+3, 0 ); 
	for( size_t idx = 0; idx < NumUnitVectors(); idx++ )
	{
		double dotprod = unit.Dot( UnitVectorAt( idx ) );

		if( dotprod > dotprods[0] )
		{
			dotprods[0] = dotprod;
			unit_indexptr[0] = idx;
			
			// resort the values
			if( dotprods[0] > dotprods[1] )
			{
				std::swap( dotprods[0], dotprods[1] );
				std::swap( unit_indexptr[0], unit_indexptr[1] );
			}
			if( dotprods[0] > dotprods[2] )
			{
				std::swap( dotprods[0], dotprods[2] );
				std::swap( unit_indexptr[0], unit_indexptr[2] );
			}
			if( dotprods[1] > dotprods[2] )
			{
				std::swap( dotprods[1], dotprods[2] );
				std::swap( unit_indexptr[1], unit_indexptr[2] );
			}
		}
	}
	return true;
}

bool SKTRAN_UnitSphere_Delaunay_nonTabledLookup::FindThreeNearestVertices_undirected( const nxVector& unit, size_t* unit_indexptr, size_t maxvertices) const
{
	bool ok = true;
//	bool keepSearching = true;

	size_t fidx=0;

	ok = ok && 2<maxvertices;
	if(ok)
	{
		for(; fidx<m_numFaces; ++fidx)
		{
			if( m_edgenorms[fidx].d[0].Dot(unit) >= 0.0 && 
				m_edgenorms[fidx].d[1].Dot(unit) >= 0.0 && 
				m_edgenorms[fidx].d[2].Dot(unit) >= 0.0)
			{
				break;	// unit is inside all faces -- it's inside this face
			}
		}

		ok = ok && fidx<m_numFaces;
		if(ok)
		{
			*(unit_indexptr+0)  = m_faces[fidx].d[0];
			*(unit_indexptr+1)  = m_faces[fidx].d[1];
			*(unit_indexptr+2)  = m_faces[fidx].d[2];
		} 
	}

	return ok;
}


bool SKTRAN_UnitSphere_Delaunay_nonTabledLookup::FindInterpolationWeights( const nxVector& unit, size_t* unit_indexptr, double* unit_weightptr, size_t faceidx ) const
{
	bool ok = true;

	double delfactor;
	// Should really do a better job of interpolation since faces may be quite large
	ok = ok && InterpolateTriangle( unit, unit_indexptr, unit_weightptr, faceidx );

	// Do not interpolate the dummy point
	if(ok && m_hasDummyPoint)
	{
		for(size_t delid=0; delid<3; ++delid)
		{
			if(NumUnitVectors()==unit_indexptr[delid] && 1e-12<unit_weightptr[delid])
			{
				delfactor = (1.0 / (1.0 - unit_weightptr[delid]));
				unit_weightptr[0] *= delfactor;
				unit_weightptr[1] *= delfactor;
				unit_weightptr[2] *= delfactor;
				unit_weightptr[delid] = 0.0;
			}
		} 
	}

	return ok;
}

bool SKTRAN_UnitSphere_Delaunay_nonTabledLookup::OptimizeForLookupInNeighbourhoodOf (const nxVector& unit)
{
	bool ok = true;

	size_t fidx = 0;

	// Find face containing #unit using undirected search
	for(; fidx<m_numFaces; ++fidx)
	{
		if( m_edgenorms[fidx].d[0].Dot(unit) >= 0.0 && 
			m_edgenorms[fidx].d[1].Dot(unit) >= 0.0 && 
			m_edgenorms[fidx].d[2].Dot(unit) >= 0.0)
		{
			break;	// unit is inside all faces -- it's inside this face
		}
	}
	ok = ok && fidx<m_numFaces;
	if(ok) m_startFace = fidx;	// Store index of "popular" face

	return ok;
}

/**********************
 * Binary lookup tables
 **********************/

SKTRAN_UnitSphere_Delaunay_binaryLookup::SKTRAN_UnitSphere_Delaunay_binaryLookup()
{

}

SKTRAN_UnitSphere_Delaunay_binaryLookup::~SKTRAN_UnitSphere_Delaunay_binaryLookup()
{

}

void SKTRAN_UnitSphere_Delaunay_binaryLookup::GetFaceIndicesOrder(size_t fidx, tuple3<size_t>& order) const
{
	if(UnitVectorAt(m_faces[fidx].d[0]).Z()<UnitVectorAt(m_faces[fidx].d[1]).Z()){
		if(UnitVectorAt(m_faces[fidx].d[2]).Z()<UnitVectorAt(m_faces[fidx].d[0]).Z()){
			order.d[0]=2; order.d[1]=0; order.d[2]=1; // 2 < 0 < 1
		} else if(UnitVectorAt(m_faces[fidx].d[2]).Z() < UnitVectorAt(m_faces[fidx].d[1]).Z()){
			order.d[0]=0; order.d[1]=2; order.d[2]=1; // 0 <= 2 < 1
		} else{
			order.d[0]=0; order.d[1]=1; order.d[2]=2; // 0 < 1 <= 2
		}
	} else{
		if(UnitVectorAt(m_faces[fidx].d[2]).Z()<UnitVectorAt(m_faces[fidx].d[1]).Z()){
			order.d[0]=2; order.d[1]=1; order.d[2]=0; // 2 < 1 <= 0
		} else if(UnitVectorAt(m_faces[fidx].d[2]).Z() < UnitVectorAt(m_faces[fidx].d[0]).Z()){
			order.d[0]=1; order.d[1]=2; order.d[2]=0; // 1 <= 2 < 0
		} else{
			order.d[0]=1; order.d[1]=0; order.d[2]=2; // 1 <= 0 <= 2
		}
	}
}

double SKTRAN_UnitSphere_Delaunay_binaryLookup::GetFaceSignedness( size_t fidx, const tuple3<size_t>& order) const
{
	double tripleProduct;

    tripleProduct =   UnitVectorAt(m_faces[fidx].d[order.d[1]]).X() * ( UnitVectorAt(m_faces[fidx].d[order.d[0]]).Y()*UnitVectorAt(m_faces[fidx].d[order.d[2]]).Z() - UnitVectorAt(m_faces[fidx].d[order.d[0]]).Z()*UnitVectorAt(m_faces[fidx].d[order.d[2]]).Y())
                    - UnitVectorAt(m_faces[fidx].d[order.d[1]]).Y() * ( UnitVectorAt(m_faces[fidx].d[order.d[0]]).X()*UnitVectorAt(m_faces[fidx].d[order.d[2]]).Z() - UnitVectorAt(m_faces[fidx].d[order.d[0]]).Z()*UnitVectorAt(m_faces[fidx].d[order.d[2]]).X())
                    + UnitVectorAt(m_faces[fidx].d[order.d[1]]).Z() * ( UnitVectorAt(m_faces[fidx].d[order.d[0]]).X()*UnitVectorAt(m_faces[fidx].d[order.d[2]]).Y() - UnitVectorAt(m_faces[fidx].d[order.d[0]]).Y()*UnitVectorAt(m_faces[fidx].d[order.d[2]]).X());

	return tripleProduct;
}

bool SKTRAN_UnitSphere_Delaunay_binaryLookup::ConstructLookupObjects( )
{
    bool            ok = true;
	bool            pushThisEdge;
	double          tripleProduct;
	tuple3<size_t>  order;
	std::vector<double>::iterator it;
	std::vector<double> zTableTemp;
	std::vector< std::vector<lutType> > lut_temp;
	nxVector		normal;

	// Get z table
	zTableTemp.reserve(NumUnitVectors());
	m_zTable.reserve(NumUnitVectors());
	for(size_t pidx=0; pidx<NumUnitVectors(); ++pidx) zTableTemp.push_back(UnitVectorAt(pidx).Z());
	std::sort(zTableTemp.begin(), zTableTemp.end());
	m_zTable.push_back(zTableTemp.front());
	for(size_t pidx=1; pidx<zTableTemp.size(); ++pidx) 
		if( 1e-12 < (zTableTemp[pidx]-m_zTable.back()))
		{
			m_zTable.push_back(zTableTemp[pidx]);
		}
	lut_temp.resize(m_zTable.size());
	
	// Get entries to put into intercept lists
	for(size_t fidx=0; fidx<m_numFaces; ++fidx)
	{
		// Find the edges where the third point in the face is below the face formed
		// by [vertex with low z]--[vertex with high z]--[origin]. That is, find edges
		// that bound faces from the counter-clockwise direction (low-\theta direction)
		// when looking at the hull in the negative-z direction. 
		// Another way of saying this: If you're looking at the hull with the z-axis
		// pointed upwards, put edges on the list such that they're associated with
		// the face that has area to the right of the edge.
		GetFaceIndicesOrder(fidx, order);				// figure out how the vertices are ordered along the z-axis
		tripleProduct = GetFaceSignedness(fidx, order);	// ([low-z-vert]x[high-z-vert]).[middle-z-vert] is positive iff edge is to the right of area
		pushThisEdge = tripleProduct < 0.0;             // should we add edge 0-2, or edges 0-1 and 1-2? Note vertices are colinear if tripleProduct==0 (error)
		
		it = std::lower_bound(m_zTable.begin(), m_zTable.end(), UnitVectorAt(m_faces[fidx].d[order.d[0]]).Z());
		if(!pushThisEdge)
		{
			// Push edge 0-1
			normal = UnitVectorAt(m_faces[fidx].d[order.d[0]]).Cross(UnitVectorAt(m_faces[fidx].d[order.d[1]]));
			ok = ok && 1e-12<normal.Magnitude();
			if(ok)
			{
				normal /= normal.Magnitude();
				for(; *it < UnitVectorAt(m_faces[fidx].d[order.d[1]]).Z(); ++it)   lut_temp[it-m_zTable.begin()].push_back(lutType(normal, fidx));
			} else{
				nxLog::Record(NXLOG_ERROR, "SKTRAN_UnitSphere_Delaunay::ConstructLookupObjects, Edge as normal zero -- vertices are too close together.");
			}
		}
		// Push edge 0-2, or 1-2 if the iterator has been moved forward in the above if-statement
		normal = UnitVectorAt(m_faces[fidx].d[order.d[pushThisEdge?0:1]]).Cross(UnitVectorAt(m_faces[fidx].d[order.d[2]]));
		ok = ok && 1e-12<normal.Magnitude();
		if(ok)
		{
			normal /= normal.Magnitude();
			for(; *it < UnitVectorAt(m_faces[fidx].d[order.d[2]]).Z(); ++it)   lut_temp[it-m_zTable.begin()].push_back(lutType(normal, fidx));
		} else{
			nxLog::Record(NXLOG_ERROR, "SKTRAN_UnitSphere_Delaunay::ConstructLookupObjects, Edge as normal zero -- vertices are too close together.");
		}
	}


    return ok;
}




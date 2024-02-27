#include "student_code.h"
#include "mutablePriorityQueue.h"

using namespace std;

namespace CGL
{

  /**
   * Evaluates one step of the de Casteljau's algorithm using the given points and
   * the scalar parameter t (class member).
   *
   * @param points A vector of points in 2D
   * @return A vector containing intermediate points or the final interpolated vector
   */
  std::vector<Vector2D> BezierCurve::evaluateStep(std::vector<Vector2D> const &points)
  { 
    // TODO Part 1.
      
    std::vector<Vector2D> nextPoints;
    for (size_t i = 0; i < points.size() - 1; i++) {
        Vector2D p1 = points[i];
        Vector2D p2 = points[i + 1];
        Vector2D intermediatePoint = p1 * (1 - t) + p2 * t; 
        nextPoints.push_back(intermediatePoint);
    }
    return nextPoints;
      
    //return std::vector<Vector2D>();
  }

  /**
   * Evaluates one step of the de Casteljau's algorithm using the given points and
   * the scalar parameter t (function parameter).
   *
   * @param points    A vector of points in 3D
   * @param t         Scalar interpolation parameter
   * @return A vector containing intermediate points or the final interpolated vector
   */
  std::vector<Vector3D> BezierPatch::evaluateStep(std::vector<Vector3D> const &points, double t) const
  {
    // TODO Part 2.
    std::vector<Vector3D> nextPoints;
    for (size_t i = 0; i < points.size() - 1; i++) {
        Vector3D p1 = points[i];
        Vector3D p2 = points[i + 1];
        Vector3D intermediatePoint = p1 * (1 - t) + p2 * t;
        nextPoints.push_back(intermediatePoint);
    }
    return nextPoints;
    //return std::vector<Vector3D>();
  }

  /**
   * Fully evaluates de Casteljau's algorithm for a vector of points at scalar parameter t
   *
   * @param points    A vector of points in 3D
   * @param t         Scalar interpolation parameter
   * @return Final interpolated vector
   */
  Vector3D BezierPatch::evaluate1D(std::vector<Vector3D> const &points, double t) const
  {
    // TODO Part 2.
   
    //This array stores the intermediate points 
    std::vector<Vector3D> nextPoints = points;

    //This loops continues to run eval steps until we are left with a single point
    while (nextPoints.size() > 1) {
        nextPoints = evaluateStep(nextPoints, t);
    }

    return nextPoints.front();
  }

  /**
   * Evaluates the Bezier patch at parameter (u, v)
   *
   * @param u         Scalar interpolation parameter
   * @param v         Scalar interpolation parameter (along the other axis)
   * @return Final interpolated vector
   */
  Vector3D BezierPatch::evaluate(double u, double v) const 
  {  
    // TODO Part 2.
    



    //Declare an empty array of 3D vectors to accumulate the column points
    std::vector<Vector3D> columnPoints;

    //loop over each of the rows in the grid and run evaluate1D on the row with parameter u as the second argument
    for (const auto& row : controlPoints) {
        // Evaluate the row at u to get a point on the Bezier curve defined by this row
        Vector3D pointOnRow = evaluate1D(row, u);
        columnPoints.push_back(pointOnRow);
    }

    //return the result of calling evaluate1D on the resulting column points with v as the second argument.
    Vector3D finalPoint = evaluate1D(columnPoints, v);

    return finalPoint;

  }

  Vector3D Vertex::normal( void ) const
  {
    // TODO Part 3.
    // Returns an approximate unit normal at this vertex, computed by
    // taking the area-weighted average of the normals of neighboring
    // triangles, then normalizing.

    //Credit: I went back and forth with chatGPT on this function, and also used the code from the 184 website primer to start with. 


    //Declare an empty array of 3D vectors to represent the normals of each face
    std::vector<Vector3D> normals;
    //Declare an empty array of scalars to represent the areas of each face
    std::vector<float> areas;


    //Iterate over each face associated with the vertex
    HalfedgeCIter h = this->halfedge(); // get the first half-edge of the vertex
    do {
        if (!h->face()->isBoundary()) { // Make sure the face is not a boundary face
            FaceCIter f = h->face(); //get the face of the current half-edge 

            // getting the verticies of this trinagle
            Vector3D p0 = h->vertex()->position;
            Vector3D p1 = h->next()->vertex()->position;
            Vector3D p2 = h->next()->next()->vertex()->position;

            //calculate the normal using the cross product
            Vector3D normal = cross(p1 - p0, p2 - p0);
            normals.push_back(normal); // Store the normal

            // Calculate the area as half the magnitude of the cross product
            float area = 0.5f * normal.norm();
            areas.push_back(area); // Store the area
        }

        h = h->twin()->next();;               // move to the next half-edge around the vertex
    } while (h != this->halfedge());    // keep going until we are back where we were

    Vector3D weighted_sum; 
    weighted_sum.x = 0; 
    weighted_sum.y = 0; 
    weighted_sum.z = 0; 

    for (int i = 0; i < areas.size(); ++i) {
        weighted_sum += normals[i] * areas[i];
    }
    
    Vector3D normalized_weighted_sum = weighted_sum.unit();
    //cout << normalized_weighted_sum << endl;

    return weighted_sum.unit();
  }

  EdgeIter HalfedgeMesh::flipEdge(EdgeIter e0)
  {
      // TODO Part 4.
      // This method should flip the given edge and return an iterator to the flipped edge.

      HalfedgeIter h0 = e0->halfedge();
      HalfedgeIter h1 = h0->twin();
      HalfedgeIter h2 = h0->next();
      HalfedgeIter h3 = h2->next();
      HalfedgeIter h4 = h1->next();
      HalfedgeIter h5 = h4->next();

      // Early exit if the edge cannot be flipped
      if (h0->isBoundary() || h1->isBoundary()) {
          return e0;
      }

      VertexIter v0 = h0->vertex();
      VertexIter v1 = h1->vertex();
      VertexIter v2 = h3->vertex();
      VertexIter v3 = h5->vertex();

      FaceIter f0 = h0->face();
      FaceIter f1 = h1->face();

      // Reassign half-edge's vertex and face pointers
      h0->vertex() = v3;
      h1->vertex() = v2;

      // Update the next pointers to reroute the half-edge cycles around the flipped edge
      h0->next() = h3;
      h3->next() = h4;
      h4->next() = h0;

      h1->next() = h5;
      h5->next() = h2;
      h2->next() = h1;

      // Update vertex's half-edge pointers
      v0->halfedge() = h2;
      v1->halfedge() = h4;
      v2->halfedge() = h1; // Point to one of the half-edges emanating from v2
      v3->halfedge() = h0; // Point to one of the half-edges emanating from v3

      // Update face's half-edge pointers
      f0->halfedge() = h0;
      f1->halfedge() = h1;

      return e0; // The edge has now been flipped
  }


  VertexIter HalfedgeMesh::splitEdge( EdgeIter e0 )
  {
    // TODO Part 5.
    // This method should split the given edge and return an iterator to the newly inserted vertex.
    // The halfedge of this vertex should point along the edge that was split, rather than the new edges.
    // 
    // Citation: I used chatGPT to help me with debugging this (especially with checking the consistency of the statements printed out by the check during the debugging process!)
      // Check if the edge is a boundary edge and should not be split
      HalfedgeIter h0 = e0->halfedge();
      if (h0->isBoundary() || h0->twin()->isBoundary()) {
          return VertexIter();
      }

      // Original halfedges, vertices, and faces around the edge
      HalfedgeIter h1 = h0->twin();
      HalfedgeIter h2 = h0->next();
      HalfedgeIter h3 = h2->next();
      HalfedgeIter h4 = h1->next();
      HalfedgeIter h5 = h4->next();

      EdgeIter missing_edge = h2->edge();

      VertexIter v0 = h0->vertex();
      VertexIter v1 = h1->vertex();
      VertexIter v2 = h3->vertex();
      VertexIter v3 = h5->vertex();

      FaceIter f0 = h0->face();
      FaceIter f1 = h1->face();

      // New geometric elements
      VertexIter v4 = newVertex();
      EdgeIter e1 = newEdge(), e2 = newEdge(), e3 = newEdge();
      HalfedgeIter h6 = newHalfedge(), h7 = newHalfedge(); // For e1
      HalfedgeIter h8 = newHalfedge(), h9 = newHalfedge(); // For e2
      HalfedgeIter h10 = newHalfedge(), h11 = newHalfedge(); // For e3
      FaceIter f2 = newFace(), f3 = newFace();

      // Set the position of the new vertex at the midpoint of the original edge
      v4->position = (v0->position + v1->position) / 2.0;

      // Save outer edges before modifications: 
      HalfedgeIter h2_outer = h2->twin();
      HalfedgeIter h3_outer = h3->twin();
      HalfedgeIter h4_outer = h4->twin();
      HalfedgeIter h5_outer = h5->twin();

      // Reassign existing halfedges to new faces and vertices
      h0->setNeighbors(h6, h1, v0, e0, f0);
      h1->setNeighbors(h4, h0, v4, e0, f1);
      h2->setNeighbors(h7, h2_outer, v1, h2->edge(), f2);
      h3->setNeighbors(h0, h3_outer, v2, h3->edge(), f0);
      h4->setNeighbors(h8, h4_outer, v0, h4->edge(), f1);
      h5->setNeighbors(h10, h5_outer, v3, h5->edge(), f3);


      // Setup new halfedges
      h6->setNeighbors(h3, h7, v4, e1, f0);
      h7->setNeighbors(h11, h6, v2, e1, f2);
      h8->setNeighbors(h1, h9, v3, e2, f1);
      h9->setNeighbors(h5, h8, v4, e2, f3);
      h10->setNeighbors(h9, h11, v1, e3, f3);
      h11->setNeighbors(h2, h10, v4, e3, f2);

      // Assign new edges to halfedges
      e0->halfedge() = h0;
      e1->halfedge() = h7;
      e2->halfedge() = h9;
      e3->halfedge() = h11;
      

      // Update faces to point to one of their halfedges
      f0->halfedge() = h3;
      f1->halfedge() = h1;
      f2->halfedge() = h11;
      f3->halfedge() = h9;

      // Update the new vertex to point to one of its outgoing halfedges
      v4->halfedge() = h8;


      return v4;
  }



  void MeshResampler::upsample( HalfedgeMesh& mesh )
  {
    // TODO Part 6.
    // This routine should increase the number of triangles in the mesh using Loop subdivision.
    // One possible solution is to break up the method as listed below.

    // 1. Compute new positions for all the vertices in the input mesh, using the Loop subdivision rule,
    // and store them in Vertex::newPosition. At this point, we also want to mark each vertex as being
    // a vertex of the original mesh.
    
    // 2. Compute the updated vertex positions associated with edges, and store it in Edge::newPosition.
    
    // 3. Split every edge in the mesh, in any order. For future reference, we're also going to store some
    // information about which subdivide edges come from splitting an edge in the original mesh, and which edges
    // are new, by setting the flat Edge::isNew. Note that in this loop, we only want to iterate over edges of
    // the original mesh---otherwise, we'll end up splitting edges that we just split (and the loop will never end!)
    
    // 4. Flip any new edge that connects an old and new vertex.

    // 5. Copy the new vertex positions into final Vertex::position.

  }
}

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

      //outer 
      HalfedgeIter h6 = h2->twin();
      HalfedgeIter h7 = h3->twin();
      HalfedgeIter h8 = h4->twin();
      HalfedgeIter h9 = h5->twin();
      
      EdgeIter e1 = h6->edge();
      EdgeIter e2 = h7->edge();
      EdgeIter e3 = h8->edge();
      EdgeIter e4 = h9->edge();


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



      // Reassign half-edge's
      h0->setNeighbors(h3, h1, v3, e0, f0);
      h3->setNeighbors(h4, h7, v2, e2, f0);
      h4->setNeighbors(h0, h8, v0, e3, f0);

      h1-> setNeighbors(h5, h0, v2, e0, f1);
      h5-> setNeighbors(h2, h9, v3, e4, f1);
      h2->setNeighbors(h1, h6, v1, e1, f1);

      h6->setNeighbors(h6->next(), h2, v2, e1, h6->face());
      h7->setNeighbors(h7->next(), h3, v0, e2, h7->face());
      h8->setNeighbors(h8->next(), h4, v3, e3, h8->face());
      h9->setNeighbors(h9->next(), h5, v1, e4, h9->face());


      // Update vertex's half-edge pointers
      v0->halfedge() = h4;
      v1->halfedge() = h2;
      v2->halfedge() = h3; // Point to one of the half-edges emanating from v2
      v3->halfedge() = h5; // Point to one of the half-edges emanating from v3

      // Update face's half-edge pointers
      f0->halfedge() = h0;
      f1->halfedge() = h1;

      //update edge
      e0->halfedge() = h0;
      e1->halfedge() = h2;
      e2->halfedge() = h3;
      e3->halfedge() = h4;
      e4->halfedge() = h5;

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

      EdgeIter h2_outer_e = h2->edge();
      EdgeIter h3_outer_e = h3->edge();
      EdgeIter h4_outer_e = h4->edge();
      EdgeIter h5_outer_e = h5->edge();

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

      //ensure consistency of outer halfedges
      h2_outer->setNeighbors(h2_outer->next(), h2, v2, h2_outer_e, h2_outer->face());
      h3_outer->setNeighbors(h3_outer->next(), h3, v0, h3_outer_e, h3_outer->face());
      h4_outer->setNeighbors(h4_outer->next(), h4, v3, h4_outer_e, h4_outer->face());
      h5_outer->setNeighbors(h5_outer->next(), h5, v1, h5_outer_e, h5_outer->face());


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

      //Set the isNew values

      e0->isNew = 0; 
      e3->isNew = 0; 
      e1->isNew = 1;
      e2->isNew = 1;

      v4->isNew = 1;
      return v4;
      

  }



  void MeshResampler::upsample( HalfedgeMesh& mesh )
  {
    // TODO Part 6.
    // This routine should increase the number of triangles in the mesh using Loop subdivision.
    // One possible solution is to break up the method as listed below.
    //Citaton: I used chatGPT pretty heavily for this function, but it took many back and forths to with it.

    //Citaton: I used chatGPT pretty heavily for this function, but it took many back and forths to with it.

    // 1. Compute new positions for all the vertices in the input mesh, using the Loop subdivision rule,
    // and store them in Vertex::newPosition. At this point, we also want to mark each vertex as being
    // a vertex of the original mesh.
      for (VertexIter v = mesh.verticesBegin(); v != mesh.verticesEnd(); v++) {
          v->isNew = 0; // Mark vertex as original


          // Compute new position for existing vertex using Loop subdivision rule
          float n = static_cast<float>(v->degree());
          float u = (n == 3) ? 3.0f / 16.0f : 3.0f / (8.0f * n);
          Vector3D original_position = v->position;
          Vector3D sum_neighbor_position(0, 0, 0);

          HalfedgeIter h = v->halfedge();
          do {
              sum_neighbor_position += h->twin()->vertex()->position;
              h = h->twin()->next();
          } while (h != v->halfedge());

          v->newPosition = (1 - n * u) * original_position + u * sum_neighbor_position;

      }

      // 2. Compute the updated vertex positions associated with edges, and store it in Edge::newPosition.
      for (EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++) {

          Vector3D A = e->halfedge()->vertex()->position;
          Vector3D B = e->halfedge()->twin()->vertex()->position;
          Vector3D C = e->halfedge()->next()->next()->vertex()->position;
          Vector3D D = e->halfedge()->twin()->next()->next()->vertex()->position;

          e->newPosition = 3.0f / 8.0f * (A + B) + 1.0f / 8.0f * (C + D);

          e->isNew = 0; // Mark edge as original 
      }

      // 3. Split every edge in the mesh, in any order. For future reference, we're also going to store some
      // information about which subdivide edges come from splitting an edge in the original mesh, and which edges
      // are new, by setting the flat Edge::isNew. Note that in this loop, we only want to iterate over edges of
      // the original mesh---otherwise, we'll end up splitting edges that we just split (and the loop will never end!)

      std::vector<EdgeIter> originalEdges;
      for (EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++) {
          originalEdges.push_back(e);
      }

      for (EdgeIter e : originalEdges) {
          //first make sure we have an old edge 
          if (!(e->halfedge()->vertex()->isNew || e->halfedge()->twin()->vertex()->isNew)) {
              VertexIter newV = mesh.splitEdge(e);
              newV->newPosition = e->newPosition;
          }
      }


      // 4. Flip any new edge that connects an old and new vertex.
      for (EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++) {
          bool vertex_A_new = e->halfedge()->vertex()->isNew;
          bool vertex_B_new = e->halfedge()->twin()->vertex()->isNew;
          if (e->isNew && ((vertex_A_new && !vertex_B_new) ||
              (!vertex_A_new && vertex_B_new))) {
              mesh.flipEdge(e);
              e->isNew = false;
          }
      }


      // 5. Copy the new vertex positions into final Vertex::position.
      for (VertexIter v = mesh.verticesBegin(); v != mesh.verticesEnd(); v++) {
          if (true) { // Only update original vertices: Condition to insert for debugging purposes: !v->isNew
              v->position = v->newPosition;
              v->isNew = 0;
          }
      }

  }
}


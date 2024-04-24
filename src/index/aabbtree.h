/*
Date: 03/04/2024
Author: Lu Chen
*/

#pragma once

#include <vector>
#include <stack>
#include <math.h>
#include <algorithm>
#include <numeric>
#include <iostream>

// #include "corridor.cuh"
#include "cuda_runtime.h"

namespace myspatial
{

class MBB_Tri
{
    public:
        MBB_Tri() = default;
        MBB_Tri(float x_min, float y_min, float z_min, float x_max, float y_max, float z_max):
        xmin(x_min), ymin(y_min), zmin(z_min), xmax(x_max), ymax(y_max), zmax(z_max) {}

    public:
        float xmin, ymin, zmin, xmax, ymax, zmax;
    
    public:
        void update(MBB_Tri &other_mbb)
        {
            xmin = std::min(xmin, other_mbb.xmin);
            ymin = std::min(ymin, other_mbb.ymin);
            zmin = std::min(zmin, other_mbb.zmin);
            xmax = std::max(xmax, other_mbb.xmax);
            ymax = std::max(ymax, other_mbb.ymax);
            zmax = std::max(zmax, other_mbb.zmax);
        }

        bool __device__ __host__ ray_intersection(float3 ray_origin, float3 ray_direction)
        {

            if (ray_origin.y > ymax || ray_origin.y < ymin || ray_origin.z > zmax || ray_origin.z < zmin) return false;
            if ((ray_direction.x > 0 && ray_origin.x > xmax) || (ray_direction.x < 0 && ray_origin.x < xmin)) return false;
            return true;

        }

};


class Node
{
    public:
        // idx in traingle array. If it is an internal node, idx is set to -1.  
        int idx;
        // pointer/idx to the left and right child
        int left;
        int right;
        // division axis
        int axis;
        // MBB of the node
        MBB_Tri mbb;

        
        // The start index of the triangle in indices array
        int start;
        // The end index of the triangle in indices array
        int end;
        // depth of the node
        int depth;

        // constructor
        Node() : idx(-1), axis(-1), left(-1), right(-1), start(-1), end(-1), depth(-1) {}
    
};


class Triangle
{
    public: 
        float3 p1;
        float3 p2;
        float3 p3;
        float3 center;
        MBB_Tri mbb;
    
    public:
        Triangle() = default;
        Triangle(float3 p_1, float3 p_2, float3 p_3):
        p1(p_1), p2(p_2), p3(p_3) {set_mbb(); set_center();}
    
    private:
        void set_mbb()
        {
            mbb.xmin = std::min({p1.x, p2.x, p3.x});
            mbb.ymin = std::min({p1.y, p2.y, p3.y});
            mbb.zmin = std::min({p1.z, p2.z, p3.z});
            
            mbb.xmax = std::max({p1.x, p2.x, p3.x});
            mbb.ymax = std::max({p1.y, p2.y, p3.y});
            mbb.zmax = std::max({p1.z, p2.z, p3.z});
        }

        void set_center()
        {
            center.x = (p1.x + p2.x + p3.x) / 3.0;
            center.y = (p1.y + p2.y + p3.y) / 3.0;
            center.z = (p1.z + p2.z + p3.z) / 3.0;           
        }


};


class AABBTree
{

    public:
        int root_;
        std::vector<Triangle> triangles_;
        // Node* node_pool = new Node();
        std::vector<Node> node_pool;
        std::vector<int> indices;
    
    public:
        int total_node = 0;

    public:
        // The constructors
        AABBTree(): root_(-1) {};
        AABBTree(const std::vector<Triangle> &triangles): root_(-1) {build(triangles); }

        // The destructor
        ~AABBTree () {clear(); }

        // rebuild static AABB tree
        void build(const std::vector<Triangle> &triangles)
        {
            clear();

            triangles_ = triangles;

            indices.resize(triangles.size());
            std::iota(std::begin(indices), std::end(indices), 0);

            root_ = buildIterative(indices.data(), (int) triangles.size());
        }

        bool point_inside(float3 point)
        {
            // odd number of intersections means the point inside.
            int n = ray_intersection_with_aabbtree(point);
            // std::cout << "the number of intersections: " << n << std::endl;
            return n % 2;
        }

        // clean AABB tree
        void clear()
        {
            // clearIterative(root_);
            root_ = -1;
            triangles_.clear();
        }

    
    private:

        class Exception: public std::exception {using std::exception::exception; };


        MBB_Tri compute_mbb(Node& node, int* indices)
        {
            int start = node.start;
            int end = node.end;
            int idx;

            MBB_Tri node_mbb = triangles_[start].mbb;
            
            for (int i = start; i <= end; i++)
            {
                idx = indices[i];
                MBB_Tri cur_mbb = triangles_[idx].mbb;
                node_mbb.update(cur_mbb);
            }

            float span_x = node_mbb.xmax - node_mbb.xmin;
            float span_y = node_mbb.ymax - node_mbb.ymin;
            float span_z = node_mbb.zmax - node_mbb.zmin;

            if (span_x > span_y && span_x > span_z)
                node.axis = 0;
            else if (span_y > span_z)
                node.axis = 1;
            else
                node.axis = 2;
            
            node.mbb = node_mbb;

            return node_mbb;
        }


        int buildIterative(int* indices, int npoints)
        {
            
            if (npoints <= 0) return -1;

            std::stack<int> stack;
            node_pool.resize(2 * npoints);

            int k = 0;
            stack.push(k);
            Node &node = node_pool[k++];
            node.start = 0;
            node.end = npoints - 1;
            node.depth = 0;

            while (!stack.empty())
            {
                int node_idx = stack.top();
                stack.pop();
                Node &node = node_pool[node_idx];
                int start = node.start;
                int end = node.end;
                int depth = node.depth;
                int axis = node.axis;

                if (start == end) 
                {
                    // Leaf Node, which triangle the node contains
                    // node.idx = start;
                    node.idx = indices[start];
                    node.mbb = triangles_[indices[start]].mbb;
                    continue;
                }

                // can be optimized using postorder traversal. So the MBB can be the union of the two child nodes. 
                // compute mbb, set mbb and division axis
                auto mbb = compute_mbb(node, indices);
                
                // partial sort 
                int mid = (start + end) / 2;

                std::nth_element(indices + start, indices + mid, indices + end + 1, [&] (int lhs, int rhs)
                {
                    float3 lcenter = triangles_[lhs].center;
                    float3 rcenter = triangles_[rhs].center;

                    if (axis == 0)
                        return lcenter.x < rcenter.x; 
                    else if (axis == 1)
                        return lcenter.y < rcenter.y;
                    return lcenter.z < rcenter.z;
                });
                
                node.left = k;
                Node& left_node = node_pool[k++];
                node.right = k;
                Node& right_node = node_pool[k++];

                left_node.start = start;
                left_node.end = mid;
                left_node.depth = depth + 1;

                right_node.start = mid + 1;
                right_node.end = end;
                right_node.depth = depth + 1;

                stack.push(node.right);
                stack.push(node.left);

            }

            // root_ always 0 if not empty 

            this->total_node = k;

            // for (int i = 0; i < k; i++)
            // {
            //     Node &node = node_pool[i];
            //     MBB_Tri &mbb = node.mbb;
            //     std::cout << "start: " << node.start << " end: " << node.end << std::endl;
            //     std::cout << mbb.xmin << ", " << mbb.ymin << ", " << mbb.zmin << ", " << mbb.xmax << ", " << mbb.ymax << ", " << mbb.zmax << std::endl;
            // }
            return 0;            
        }


        int ray_intersection_with_aabbtree(float3 ray_origin)
        {
            if (triangles_.size() <= 0) return 0;

            // bfs search
            std::stack<int> stack;
            stack.push(0);
            float intersections = 0.0;

            // mbb of all the triangles of the mesh
            float3 ray_direction;
            MBB_Tri mbb = node_pool[0].mbb;
            if (mbb.xmax + mbb.xmin < 2 * ray_origin.x) ray_direction = make_float3(1, 0, 0);
            else ray_direction = make_float3(-1, 0, 0);

            while (!stack.empty())
            {
                int node_idx = stack.top();
                stack.pop();
                
                Node &node = node_pool[node_idx];
                MBB_Tri &mbb = node.mbb;
                
                if (mbb.ray_intersection(ray_origin, ray_direction))
                {
                    // leave node
                    if (node.start == node.end)
                    {
                        // Triangle &triangle = triangles_[indices[node.start]];
                        Triangle &triangle = triangles_[node.idx];
                        float3 triangle_f[3] = {triangle.p1, triangle.p2, triangle.p3};
                        intersections += ray_triangle_intersection(triangle_f, ray_origin, ray_direction);

                    }
                    else
                    {
                        stack.push(node.left);
                        stack.push(node.right);
                    }
                }

            }
            // std::cout << "float intersections: " << intersections << std::endl;
            return (int)intersections;

        }


};


class AABBTreeCUDA
{
    public:
        Node* nodes_;
        Triangle* triangles_;
    
    public:
        // constructor
        AABBTreeCUDA(const AABBTree &aabbtree)
        {
            int n_nodes = aabbtree.node_pool.size();
            int n_triangles = aabbtree.triangles_.size();

            // allocate memory
            nodes_ = new Node[n_nodes];
            triangles_ = new Triangle[n_triangles];

            auto &node_pool = aabbtree.node_pool;
            auto &triangles = aabbtree.triangles_;
            std::copy(node_pool.begin(), node_pool.end(), nodes_);
            std::copy(triangles.begin(), triangles.end(), triangles_);

        }

        AABBTreeCUDA(const AABBTreeCUDA& copySource)
        {
            int n_nodes = sizeof(copySource.nodes_) / sizeof(Node);
            int n_triangles = sizeof(copySource.triangles_) / sizeof(Triangle);

            nodes_ = new Node[n_nodes];
            triangles_ = new Triangle[n_triangles];

            std::memcpy(nodes_, copySource.nodes_, n_nodes * sizeof(Node));
            std::memcpy(triangles_, copySource.triangles_, n_triangles * sizeof(Triangle));

        }

        int __device__ __host__ ray_intersection_with_aabbtree(float3 ray_origin)
        {

            int n_nodes = sizeof(nodes_) / sizeof(Node);
            int n_triangles = sizeof(triangles_) / sizeof(Triangle);

            if (n_triangles <= 0) return 0;

            // bfs search
            int* stack = new int[n_nodes];
            int k = 0;
            stack[k++] = 0;
            float intersections = 0.0;

            // mbb of all the triangles of the mesh
            float3 ray_direction;
            MBB_Tri &mbb = nodes_[0].mbb;
            if (mbb.xmax + mbb.xmin < 2 * ray_origin.x) ray_direction = make_float3(1, 0, 0);
            else ray_direction = make_float3(-1, 0, 0);

            while (k)
            {
                int node_idx = stack[k--];
                Node &node = nodes_[node_idx];
                MBB_Tri &mbb = node.mbb;
                
                // leave node
                if (node.start == node.end)
                {
                    // Triangle &triangle = triangles_[indices[node.start]];
                    Triangle &triangle = triangles_[node.idx];
                    float3 triangle_f[3] = {triangle.p1, triangle.p2, triangle.p3};
                    intersections += ray_triangle_intersection(triangle_f, ray_origin, ray_direction);

                }
                else if (mbb.ray_intersection(ray_origin, ray_direction))
                {
                    stack[k++] = node.left;
                    stack[k++] = node.right;
                }

            }
            // std::cout << "float intersections: " << intersections << std::endl;
            return (int)intersections;

        }
    

};

}
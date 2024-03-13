/*
Date: 03/04/2024
Author: Lu Chen
*/

#pragma once

#include <vector>

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
        MBB mbb;

        // constructor
        Node() : idx(-1), axis(-1), left(-1), right(-1) {}
};

class MBB
{
    public:
        MBB() = default;
        MBB(float x_min, float y_min, float z_min, float x_max, float y_max, float z_max):
        xmin(x_min), ymin(y_min), zmin(z_min), xmax(x_max), ymax(y_max), zmax(z_max) {}
    
    public:
        float xmin, ymin, zmin, xmax, ymax, zmax;

};

class Triangle
{
    public: 
        float3 p1;
        float3 p2;
        float3 p3;
        float3 center;
        MBB mbb;
    
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
        // The constructors
        AABBTree: root_(-1) {};
        AABBTree(const std::vector<Triangle> &triangles): root_(-1) {build(triangles); }

        // The destructor
        ~AABBTree () {clear(); }

        // rebuild static AABB tree
        void build(const std::vector<Triangle> &triangles)
        {
            clear();

            triangles_ = triangles;

            std::vector<int> indices(triangles.size());
            std::iota(std::begin(indices), std::end(indices), 0);

            root_ = buildIterative(indices.data(), (int) points.size(), 0);
        }

        // clean AABB tree
        void clear()
        {
            clearIterative(root_);
            root_ = -1;
            triangles_.clear();
        }

    
    private:

        class Exception: public std::exception {using std::exception::exception; };


        MBB compute_mbb(Node* node, int* indices)
        {
            int start = node.start;
            int end = node.end;
            int idx;

            MBB node_mbb = triangles_[start].mbb;
            
            for (int i = start; i <= end; i++)
            {
                idx = indices[i];
                MBB cur_mbb = triangles_[idx].mbb;
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

            Node* node_pool = new Node();
            std::vector<Node> node_pool(3 * npoints);
            
            int k = 0;
            stack.push(k);
            Node &node = node_pool[k++];
            node.start = 0;
            node.end = npoints - 1;
            node.depth = 0;

            while (!stack.empty())
            {
                int node_idx = stack.pop();
                Node &node = node_pool[node_idx];
                int start = node.start;
                int end = node.end;
                int depth = node.depth;

                if (start == end) 
                {
                    // Leaf Node, which triangle the node contains
                    node.idx = start;
                    continue;
                }

                // compute mbb, set mbb and division axis
                compute_mbb(node, *indices);
                
                // partial sort 
                int mid = (start + end) / 2;

                std::nth_element(indices + start, indices + mid, imdices + end + 1, [&] (int lhs, int rhs)
                {
                    return triangles_[lhs][center][axis] < triangles_[rhs][center][axis];
                });
                
                node.left = k;
                left_node = node_pool[k++];
                node.right = k;
                right_node = node_pool[k++];

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
            return 0;            
        }

    
    private:
        int root_;
        std::vector<Triangle> triangles_;

}
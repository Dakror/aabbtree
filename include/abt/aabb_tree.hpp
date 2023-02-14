/*
  Copyright (c) 2009 Erin Catto http://www.box2d.org
  Copyright (c) 2016-2018 Lester Hedges <lester.hedges+aabbcc@gmail.com>

  This software is provided 'as-is', without any express or implied
  warranty. In no event will the authors be held liable for any damages
  arising from the use of this software.

  Permission is granted to anyone to use this software for any purpose,
  including commercial applications, and to alter it and redistribute it
  freely, subject to the following restrictions:

  1. The origin of this software must not be misrepresented; you must not
     claim that you wrote the original software. If you use this software
     in a product, an acknowledgment in the product documentation would be
     appreciated but is not required.

  2. Altered source versions must be plainly marked as such, and must not be
     misrepresented as being the original software.

  3. This notice may not be removed or altered from any source distribution.

  This code was adapted from parts of the Box2D Physics Engine,
  http://www.box2d.org
*/

#ifndef _AABB_H
#define _AABB_H

#include <algorithm>
#include <array>
#include <cassert>
#include <stdexcept>
#include <vector>

#include <limits>

#include <unordered_map>
#include <numeric>

//#define vector moo

namespace abt {
template <unsigned Dim, typename ValTy = double>
struct point {
  using value_type = ValTy;

  template <typename... Values>
  point(Values... values) : values{value_type(values)...}
  {
  }
  
  
  template <typename Cont>
  constexpr point(Cont const& c)
  {
      for(unsigned i = 0; i < Dim; ++i) values[i] = c[i];
  }
  
  point() = default;

  point(const point &other) : values(other.values){}
  point(point &&other) : values(std::move(other).values) {}

  point& operator=(const point& other) {
    values = other.values;
    return *this;
  }
  point& operator=(point&& other) {
    values = std::move(other).values;
    return *this;
  }

  bool operator==(const point &other) const { return values == other.values; }

  const value_type &x() const { return values[0]; }
  const value_type &y() const { return values[1]; }
  const value_type &z() const { return values[2]; }

  value_type &x() { return values[0]; }
  value_type &y() { return values[1]; }
  value_type &z() { return values[2]; }

  friend point<Dim, ValTy> operator+(const point<Dim, ValTy> &a,
                                     const std::array<ValTy, Dim> &b)
  {
    point<Dim, ValTy> res;
    for (auto d = 0; d < Dim; ++d) {
      res[d] = a[d] + b[d];
    }
    return res;
  }

  friend point<Dim, ValTy> operator-(const point<Dim, ValTy> &a,
                                     const std::array<ValTy, Dim> &b)
  {
    point<Dim, ValTy> res;
    for (auto d = 0; d < Dim; ++d) {
      res[d] = a[d] - b[d];
    }
    return res;
  }

  point &operator+=(const std::array<ValTy, Dim> &rhs)
  {
    *this = *this + rhs;
    return *this;
  }

  point& operator-=(const std::array<ValTy, Dim>& rhs) {
    *this = *this - rhs;
    return *this;
  }

  constexpr size_t size() const { return Dim; }

  value_type &operator[](size_t idx) { return values[idx]; }
  const value_type &operator[](size_t idx) const { return values[idx]; }

  std::array<value_type, Dim> values = {};
};
/*! \brief The axis-aligned bounding box object.

    Axis-aligned bounding boxes (AABBs) store information for the minimum
    orthorhombic bounding-box for an object. Support is provided for
    dimensions >= 2. (In 2D the bounding box is either a rectangle,
    in 3D it is a rectangular prism.)

    Class member functions provide functionality for merging AABB objects
    and testing overlap with other AABBs.
 */
template <unsigned Dim, typename ValTy = double>
class aabb {
 public:
  using point = abt::point<Dim, ValTy>;
  using value_type = ValTy;

  /// Constructor.
  aabb() = default;

  //! Constructor.
  /*! \param lowerBound_
          The lower bound in each dimension.

      \param upperBound_
          The upper bound in each dimension.
   */
  aabb(const point &lower_bound, const point &upper_bound)
      : lowerBound(lower_bound), upperBound(upper_bound)
  {
    surfaceArea = compute_surface_area();
    centre = compute_center();
  }

  bool operator==(const aabb& other) const {
    return lowerBound == other.lowerBound && upperBound == other.upperBound;
  }

  static aabb of_sphere(const point &center, value_type radius)
  {
    point lb, ub;
    for (unsigned int i = 0; i < Dim; i++) {
      lb[i] = center[i] - radius;
      ub[i] = center[i] + radius;
    }
    return {lb, ub};
  }

  /// Compute the surface area of the box.
  value_type compute_surface_area() const
  {
    // Sum of "area" of all the sides.
    value_type sum = 0;

    // General formula for one side: hold one dimension constant
    // and multiply by all the other ones.
    for (unsigned int d1 = 0; d1 < Dim; d1++) {
      // "Area" of current side.
      value_type product = 1;

      for (unsigned int d2 = 0; d2 < Dim; d2++) {
        if (d1 == d2)
          continue;

        value_type dx = upperBound[d2] - lowerBound[d2];
        product *= dx;
      }

      // Update the sum.
      sum += product;
    }

    return 2 * sum;
  }

  /// Get the surface area of the box.
  value_type get_surface_area() const { return surfaceArea; }

  //! Merge two AABBs into this one.
  /*! \param aabb1
          A reference to the first AABB.

      \param aabb2
          A reference to the second AABB.
   */
  void merge(const aabb &aabb1, const aabb &aabb2)
  {
    for (unsigned int i = 0; i < Dim; i++) {
      lowerBound[i] = std::min(aabb1.lowerBound[i], aabb2.lowerBound[i]);
      upperBound[i] = std::max(aabb1.upperBound[i], aabb2.upperBound[i]);
    }

    surfaceArea = compute_surface_area();
    centre = compute_center();
  }

  //! Test whether the AABB is contained within this one.
  /*! \param aabb
          A reference to the AABB.

      \return
          Whether the AABB is fully contained.
   */
  bool contains(const aabb &aabb) const
  {
    for (unsigned int i = 0; i < Dim; i++) {
      if (aabb.lowerBound[i] < lowerBound[i])
        return false;
      if (aabb.upperBound[i] > upperBound[i])
        return false;
    }

    return true;
  }

  //! Test whether the AABB overlaps this one.
  /*! \param aabb
          A reference to the AABB.

      \param touchIsOverlap
          Does touching constitute an overlap?

      \return
          Whether the AABB overlaps.
   */
  bool overlaps(const aabb &aabb, bool touchIsOverlap) const
  {
    if (touchIsOverlap) {
      for (unsigned int i = 0; i < Dim; ++i) {
        if (aabb.upperBound[i] < lowerBound[i] ||
            aabb.lowerBound[i] > upperBound[i]) {
          return false;
        }
      }
    }
    else {
      for (unsigned int i = 0; i < Dim; ++i) {
        if (aabb.upperBound[i] <= lowerBound[i] ||
            aabb.lowerBound[i] >= upperBound[i]) {
          return false;
        }
      }
    }

    return true;
  }

  //! Test whether the point overlaps this one.
  /*! \param point
          A reference to the point.

      \param touchIsOverlap
          Does touching constitute an overlap?

      \return
          Whether the AABB overlaps.
   */
  bool overlaps(const point &pt, bool touchIsOverlap) const
  {
    if (touchIsOverlap) {
      for (unsigned int i = 0; i < Dim; ++i) {
        if (pt[i] < lowerBound[i] || pt[i] > upperBound[i]) {
          return false;
        }
      }
    }
    else {
      for (unsigned int i = 0; i < Dim; ++i) {
        if (pt[i] <= lowerBound[i] || pt[i] >= upperBound[i]) {
          return false;
        }
      }
    }

    return true;
  }

  //! Compute the centre of the AABB.
  /*! \returns
          The position vector of the AABB centre.
   */
  point compute_center()
  {
    point center;

    for (unsigned int i = 0; i < Dim; i++)
      center[i] = 0.5 * (lowerBound[i] + upperBound[i]);

    return center;
  }

  /// Lower bound of AABB in each dimension.
  point lowerBound;

  /// Upper bound of AABB in each dimension.
  point upperBound;

  /// The position of the AABB centre.
  point centre;

  /// The AABB's surface area.
  value_type surfaceArea = 0;

  friend aabb<Dim, ValTy> operator-(const aabb<Dim, ValTy>& lhs, const std::array<ValTy, Dim>& rhs) {
    auto res = lhs;
    res.lowerBound -= rhs;
    res.upperBound -= rhs;
    return res;
  }

  friend aabb<Dim, ValTy> operator+(const aabb<Dim, ValTy> &lhs,
                                    const std::array<ValTy, Dim> &rhs)
  {
    auto res = lhs;
    res.lowerBound += rhs;
    res.upperBound += rhs;
    return res;
  }
};

template <unsigned Dim, typename ValTy = double>
aabb<Dim, ValTy> fattened(aabb<Dim, ValTy> bb, double skin_thickness)
{
  // Fatten the AABB.
  for (unsigned int i = 0; i < Dim; i++) {
    auto sz = bb.upperBound[i] - bb.lowerBound[i];
    bb.lowerBound[i] -= skin_thickness * sz;
    bb.upperBound[i] += skin_thickness * sz;
  }
  bb.surfaceArea = bb.compute_surface_area();
  bb.centre = bb.compute_center();
  return bb;
}


enum visit_action : char { visit_stop, visit_continue };

/*! \brief The dynamic AABB tree.

    The dynamic AABB tree is a hierarchical data structure that can be used
    to efficiently query overlaps between objects of arbitrary shape and
    size that lie inside of a simulation box. Support is provided for
    periodic and non-periodic boxes, as well as boxes with partial
    periodicity, e.g. periodic along specific axes.
 */
template <unsigned Dim, typename ValTy = double>
class tree {
  static_assert(Dim > 0, "0-dimensional tree is not supported");
 public:
  using value_type = ValTy;
  using aabb = abt::aabb<Dim, value_type>;
  using point = abt::point<Dim, value_type>;
  template <typename Ty>
  using vec = std::array<Ty, Dim>;

  enum class node_id : uint32_t {};
  
  float skin_width = 0.01f;

 private:
  static vec<ValTy> minimum_image_shift(const vec<ValTy> &bounds,
                                 const vec<ValTy> &delta)
  {
    vec<ValTy> shift = {};

    for (unsigned int i = 0; i < Dim; i++) {
      if (delta[i] < -(bounds[i] / 2)) {
        shift[i] = -bounds[i];
      }
      else if (delta[i] >= bounds[i] / 2) {
        shift[i] = bounds[i];
      }
    }
    return shift;
  }
  /*! \brief A node of the AABB tree.

   Each node of the tree contains an AABB object which corresponds to a
   entry, or a group of entrys, in the simulation box. The AABB
   objects of individual entrys are "fattened" before they are stored
   to avoid having to continually update and rebalance the tree when
   displacements are small.

   Nodes are aware of their position within in the tree. The isLeaf member
   function allows the tree to query whether the node is a leaf, i.e. to
   determine whether it holds a single entry.
  */

  static constexpr unsigned int NULL_NODE = 0xffffffff;

  struct node {
    /// Constructor.
    node() = default;

    /// The fattened axis-aligned bounding box.
    aabb bb;

    /// Index of the parent node.
    unsigned int parent = NULL_NODE;

    /// Index of the next node.
    unsigned int next = 0;

    /// Index of the left-hand child.
    unsigned int left = NULL_NODE;

    /// Index of the right-hand child.
    unsigned int right = NULL_NODE;

    /// Height of the node. This is 0 for a leaf and -1 for a free node.
    int height = 0;
    
    void* user_data = nullptr;

    //! Test whether the node is a leaf.
    /*! \return
            Whether the node is a leaf node.
     */
    bool isLeaf() const { return (height == 0); }
  };

 public:
  //! Constructor (non-periodic).
  /*! \param skin_thickness
          The skin thickness for fattened AABBs, as a fraction
          of the AABB base length.

      \param initial_size
          The number of entries (for fixed entry number systems).

   */
  tree(unsigned int initial_size = 16)
  {
    // Initialise the tree.
    m_root = NULL_NODE;
    m_node_count = 0;
    m_node_capacity = initial_size;
    m_nodes.resize(m_node_capacity);

    // Build a linked list for the list of free nodes.
    for (unsigned int i = 0; i < m_node_capacity - 1; i++) {
      m_nodes[i].next = i + 1;
      m_nodes[i].height = -1;
    }
    m_nodes[m_node_capacity - 1].next = NULL_NODE;
    m_nodes[m_node_capacity - 1].height = -1;

    // Assign the index of the first free node.
    m_free_list = 0;
  }

  /// Build an optimal tree.
  tree(std::vector<aabb> bbs)
  {
    auto count = bbs.size();

    m_root = NULL_NODE;
    m_node_count = count;
    m_leaf_count = count;
    m_node_capacity = count * 2;
    m_nodes.resize(m_node_capacity);

    // Build a linked list for the list of free nodes.
    for (unsigned int i = count; i < m_node_capacity - 1; i++) {
      m_nodes[i].next = i + 1;
      m_nodes[i].height = -1;
    }
    m_nodes[m_node_capacity - 1].next = NULL_NODE;
    m_nodes[m_node_capacity - 1].height = -1;

    // Assign the index of the first free node.
    m_free_list = count;

    std::vector<unsigned int> node_indices(count);
    std::iota(node_indices.begin(), node_indices.end(), 0);

    for (unsigned i = 0; i < count; ++i) {
      m_nodes[i].bb = bbs[i];
    }

    while (count > 1) {
      double minCost = std::numeric_limits<double>::max();
      int iMin = -1, jMin = -1;

      for (unsigned int i = 0; i < count; i++) {
        aabb aabbi = m_nodes[node_indices[i]].bb;

        for (unsigned int j = i + 1; j < count; j++) {
          aabb aabbj = m_nodes[node_indices[j]].bb;
          aabb aabb;
          aabb.merge(aabbi, aabbj);
          double cost = aabb.get_surface_area();

          if (cost < minCost) {
            iMin = i;
            jMin = j;
            minCost = cost;
          }
        }
      }

      unsigned int index1 = node_indices[iMin];
      unsigned int index2 = node_indices[jMin];

      unsigned int parent = allocate_node();
      m_nodes[parent].left = index1;
      m_nodes[parent].right = index2;
      m_nodes[parent].height =
          1 + std::max(m_nodes[index1].height, m_nodes[index2].height);
      m_nodes[parent].bb.merge(m_nodes[index1].bb, m_nodes[index2].bb);
      m_nodes[parent].parent = NULL_NODE;

      m_nodes[index1].parent = parent;
      m_nodes[index2].parent = parent;

      node_indices[jMin] = node_indices[count - 1];
      node_indices[iMin] = parent;
      count--;
    }

    m_root = node_indices[0];

    validate();
  }

  //! Insert a entry into the tree (arbitrary shape with bounding box).
  /*! \param index
          The index of the entry.

      \param lowerBound
          The lower bound in each dimension.

      \param upperBound
          The upper bound in each dimension.
   */
  node_id insert(const aabb &bb, void* user_data = nullptr)
  {
    // Allocate a new node for the entry.
    unsigned int node_idx = allocate_node();
    auto &node = m_nodes[node_idx];
    node.bb = fattened(bb, skin_width);
    node.user_data = user_data;

    // Zero the height.
    node.height = 0;

    // Insert a new leaf into the tree.
    insert_leaf(node_idx);
    return to_id(node_idx);
  }
  
  void* data(node_id id) const {
    auto node = to_unsigned(id);
    auto &n = m_nodes[node];
    
    assert(node < m_node_capacity);
    assert(n.isLeaf());
    
    return n.user_data;
  }

  /// Return the number of entrys in the tree.
  unsigned int size() const { return m_leaf_count; }

  //! Remove a entry from the tree.
  /*! \param entry
          The entry index (entryMap will be used to map the node).
   */

  static unsigned to_unsigned(node_id id) { return static_cast<unsigned>(id); }
  static node_id to_id(unsigned node) { return static_cast<node_id>(node); }

  void remove(node_id node_id)
  {
    auto node = to_unsigned(node_id);
    remove_leaf(node);
    free_node(node);
  }

  /// Remove all entrys from the tree.
  void clear()
  {
    for_each([this](node_id id, const auto &) { remove(id); });
  }

  //! Update the tree if a entry moves outside its fattened AABB.
  /*! \param entry
          The entry index (entryMap will be used to map the node).

      \param lowerBound
          The lower bound in each dimension.

      \param upperBound
          The upper bound in each dimension.

      \param alwaysReinsert
          Always reinsert the entry, even if it's within its old AABB
     (default: false)
   */
  bool update(node_id id, aabb bb, bool always_reinsert = false)
  {
    auto node = to_unsigned(id);
    auto &n = m_nodes[node];

    assert(node < m_node_capacity);
    assert(n.isLeaf());

    // No need to update if the entry is still within its fattened AABB.
    if (!always_reinsert && n.bb.contains(bb))
      return false;

    // Remove the current leaf.
    remove_leaf(node);

    // Assign the new AABB.
    n.bb = fattened(bb, skin_width);

    // Insert a new leaf node.
    insert_leaf(node);

    return true;
  }

  //! Query the tree to find candidate interactions for an AABB.
  /*! \param aabb
          The AABB.
      \param out
          An output iterator

      \return entrys
          A vector of entry indices.
   */
  template <class Query, class OutputIterator>
  void get_overlaps(const Query &query,
                    OutputIterator out,
                    bool include_touch = true,
                    const vec<ValTy> &bounds = {}) const
  {
    visit_overlaps(
        query, [&](node_id id) { *out++ = id; }, include_touch, bounds);
  }

  //! Query the tree to find candidate interactions for an AABB.
  /*! \param aabb
          The AABB.

      \return entrys
          A vector of entry indices.
   */
  template <class Query>
  std::vector<node_id> get_overlaps(const Query &query,
                                         bool include_touch = true,
                                         const vec<ValTy> &bounds = {}) const
  {
    std::vector<node_id> overlaps;
    visit_overlaps(
        query, [&](node_id id) { overlaps.push_back(id); }, include_touch, bounds);
    return overlaps;
  }

  template <class Query>
  bool any_overlap(const Query &query,
                   bool include_touch = true,
                   const vec<ValTy> &bounds = {}) const
  {
    return any_overlap(
        query, [] { return true; }, include_touch, bounds);
  }

  template <class Query, class Fn>
  bool any_overlap(const Query &query, Fn &&fn, bool include_touch = true,
                   const vec<ValTy> &bounds = {}) const
  {
    bool overlap = false;
    auto wrap_fn = [&overlap, &fn](node_id id, const aabb &bb) {
      bool success = call_with_args(std::forward<Fn>(fn), id, bb);
      overlap |= success;
      return success ? visit_stop : visit_continue;
    };
    visit_overlaps(query, wrap_fn, include_touch, bounds);
    return overlap;
  }

  template <class Fn>
  void for_each(Fn &&fn) const
  {
    for (auto idx = 0ull; idx < m_nodes.size(); ++idx) {
      const auto &node = m_nodes[idx];
      if (node.isLeaf()) {
        call_with_args(std::forward<Fn>(fn), to_id(idx), node.bb);
      }
    }
  }

  template <class Query, class Fn>
  void visit_overlaps(const Query &query,
                      Fn &&fn,
                      bool include_touch = true,
                      const vec<ValTy> &bounds = {}) const
  {
    static thread_local std::vector<unsigned int> stack(64);
    return visit_overlaps(query, std::forward<Fn>(fn), include_touch, bounds, stack);
  }

  template <class Query, class Fn>
  void visit_overlaps(const Query &query,
                      Fn &&fn,
                      bool include_touch,
                      const vec<ValTy> &bounds, 
                      std::vector<unsigned> &stack) const
  {
    constexpr bool query_is_point = std::is_same_v<Query, point>;
    constexpr bool query_is_aabb = std::is_same_v<Query, aabb>;
    static_assert(query_is_point || query_is_aabb,
                  "Only point or aabb queries are supported");


    using rt = decltype(call_with_args(std::forward<Fn>(fn), to_id(0), aabb{}));
    constexpr bool fn_returns_action = std::is_convertible_v<rt, visit_action>;
    static_assert(fn_returns_action || std::is_same_v<rt, void>,
                  "Only void or visit_action return types are allowed");

    // Make sure the tree isn't empty.
    if (size() == 0) {
      return;
    }

    stack.clear();

    stack.push_back(m_root);

    while (!stack.empty()) {
      unsigned int node = stack.back();
      stack.pop_back();

      // Copy the AABB.
      const auto &nodeAABB = m_nodes[node].bb;

      auto shifted_query = query;
      if (node == NULL_NODE)
        continue;

      vec<ValTy> separation = {};
      point center;
      if constexpr (query_is_point) {
        center = query;
      }
      else {
        center = query.centre;
      }
      for (unsigned int i = 0; i < Dim; i++) {
        separation[i] = nodeAABB.centre[i] - center[i];
      }

      shifted_query = shifted_query + minimum_image_shift(bounds, separation);
      

      // Test for overlap between the AABBs.
      if (nodeAABB.overlaps(shifted_query, include_touch)) {
        // Check that we're at a leaf node.
        if (m_nodes[node].isLeaf()) {
          const auto &n = m_nodes[node];
          if constexpr (fn_returns_action) {
            if (call_with_args(std::forward<Fn>(fn), to_id(node), n.bb) == visit_stop) {
              return;
            }
          }
          else {
            call_with_args(std::forward<Fn>(fn), to_id(node), n.bb);
          }
        }
        else {
          stack.push_back(m_nodes[node].left);
          stack.push_back(m_nodes[node].right);
        }
      }
    }
  }

  //! Get a entry AABB.
  /*! \param entry
          The entry index.
   */
  const aabb &get_aabb(node_id node) const { return m_nodes[to_unsigned(node)].bb; }

  //! Get the height of the tree.
  /*! \return
          The height of the binary tree.
   */
  unsigned int get_height() const
  {
    if (m_root == NULL_NODE)
      return 0;
    return m_nodes[m_root].height;
  }

  /// Validate the tree.
  void validate() const
  {
#ifndef NDEBUG
    validate_structure(m_root);
    validate_metrics(m_root);

    unsigned int freeCount = 0;
    unsigned int freeIndex = m_free_list;

    while (freeIndex != NULL_NODE) {
      assert(freeIndex < m_node_capacity);
      freeIndex = m_nodes[freeIndex].next;
      freeCount++;
    }

    assert(get_height() == compute_height());
    assert((m_node_count + freeCount) == m_node_capacity);
#endif
  }

  /// Rebuild an optimal tree.
  void rebuild()
  {
    std::vector<unsigned int> nodeIndices(m_node_count);
    unsigned int count = 0;

    for (unsigned int i = 0; i < m_node_capacity; i++) {
      // Free node.
      if (m_nodes[i].height < 0)
        continue;

      if (m_nodes[i].isLeaf()) {
        m_nodes[i].parent = NULL_NODE;
        nodeIndices[count] = i;
        count++;
      }
      else
        free_node(i);
    }

    while (count > 1) {
      double minCost = std::numeric_limits<double>::max();
      int iMin = -1, jMin = -1;

      for (unsigned int i = 0; i < count; i++) {
        aabb aabbi = m_nodes[nodeIndices[i]].bb;

        for (unsigned int j = i + 1; j < count; j++) {
          aabb aabbj = m_nodes[nodeIndices[j]].bb;
          aabb aabb;
          aabb.merge(aabbi, aabbj);
          double cost = aabb.get_surface_area();

          if (cost < minCost) {
            iMin = i;
            jMin = j;
            minCost = cost;
          }
        }
      }

      unsigned int index1 = nodeIndices[iMin];
      unsigned int index2 = nodeIndices[jMin];

      unsigned int parent = allocate_node();
      m_nodes[parent].left = index1;
      m_nodes[parent].right = index2;
      m_nodes[parent].height =
          1 + std::max(m_nodes[index1].height, m_nodes[index2].height);
      m_nodes[parent].bb.merge(m_nodes[index1].bb, m_nodes[index2].bb);
      m_nodes[parent].parent = NULL_NODE;

      m_nodes[index1].parent = parent;
      m_nodes[index2].parent = parent;

      nodeIndices[jMin] = nodeIndices[count - 1];
      nodeIndices[iMin] = parent;
      count--;
    }

    m_root = nodeIndices[0];

    validate();
  }

 private:
  /// The index of the root node.
  unsigned int m_root;

  /// The dynamic tree.
  std::vector<node> m_nodes;

  /// The current number of nodes in the tree.
  unsigned int m_node_count;

  unsigned int m_leaf_count = 0;

  /// The current node capacity.
  unsigned int m_node_capacity;

  /// The position of node at the top of the free list.
  unsigned int m_free_list;
  
 private:
  template <class Fn>
  static decltype(auto) call_with_args(Fn &&fn, node_id id, const aabb &bb)
  {
    constexpr bool call_id_bb = std::is_invocable_v<Fn, node_id, aabb>;
    constexpr bool call_id = std::is_invocable_v<Fn, node_id>;
    constexpr bool call_bb = std::is_invocable_v<Fn, aabb>;
    constexpr bool call_none = std::is_invocable_v<Fn>;

    static_assert(call_id_bb || call_id || call_bb || call_none,
                  "Callback has unsupported signature");
    if constexpr (call_id_bb) {
      return std::forward<Fn>(fn)(id, bb);
    }
    else if constexpr (call_id) {
      return std::forward<Fn>(fn)(id);
    }
    else if constexpr (call_bb) {
      return std::forward<Fn>(fn)(bb);
    }
    else {
      return std::forward<Fn>(fn)();
    }
  }
  //! Allocate a new node.
  /*! \return
          The index of the allocated node.
   */
  unsigned int allocate_node()
  {
    // Exand the node pool as needed.
    if (m_free_list == NULL_NODE) {
      assert(m_node_count == m_node_capacity);

      // The free list is empty. Rebuild a bigger pool.
      m_node_capacity *= 2;
      m_nodes.resize(m_node_capacity);

      // Build a linked list for the list of free nodes.
      for (unsigned int i = m_node_count; i < m_node_capacity - 1; i++) {
        m_nodes[i].next = i + 1;
        m_nodes[i].height = -1;
      }
      m_nodes[m_node_capacity - 1].next = NULL_NODE;
      m_nodes[m_node_capacity - 1].height = -1;

      // Assign the index of the first free node.
      m_free_list = m_node_count;
    }

    // Peel a node off the free list.
    unsigned int node = m_free_list;
    m_free_list = m_nodes[node].next;
    m_nodes[node].parent = NULL_NODE;
    m_nodes[node].left = NULL_NODE;
    m_nodes[node].right = NULL_NODE;
    m_nodes[node].height = 0;
    m_node_count++;

    return node;
  }

  //! Free an existing node.
  /*! \param node
          The index of the node to be freed.
   */
  void free_node(unsigned int node)
  {
    assert(node < m_node_capacity);
    assert(0 < m_node_count);

    m_nodes[node].next = m_free_list;
    m_nodes[node].height = -1;
    m_free_list = node;
    m_node_count--;
  }

  //! Insert a leaf into the tree.
  /*! \param leaf
          The index of the leaf node.
   */
  void insert_leaf(unsigned int leaf)
  {
    ++m_leaf_count;
    if (m_root == NULL_NODE) {
      m_root = leaf;
      m_nodes[m_root].parent = NULL_NODE;
      return;
    }

    // Find the best sibling for the node.

    aabb leafAABB = m_nodes[leaf].bb;
    unsigned int index = m_root;

    while (!m_nodes[index].isLeaf()) {
      // Extract the children of the node.
      unsigned int left = m_nodes[index].left;
      unsigned int right = m_nodes[index].right;

      double surfaceArea = m_nodes[index].bb.get_surface_area();

      aabb combinedAABB;
      combinedAABB.merge(m_nodes[index].bb, leafAABB);
      double combinedSurfaceArea = combinedAABB.get_surface_area();

      // Cost of creating a new parent for this node and the new leaf.
      double cost = 2.0 * combinedSurfaceArea;

      // Minimum cost of pushing the leaf further down the tree.
      double inheritanceCost = 2.0 * (combinedSurfaceArea - surfaceArea);

      // Cost of descending to the left.
      double costLeft;
      if (m_nodes[left].isLeaf()) {
        aabb aabb;
        aabb.merge(leafAABB, m_nodes[left].bb);
        costLeft = aabb.get_surface_area() + inheritanceCost;
      }
      else {
        aabb aabb;
        aabb.merge(leafAABB, m_nodes[left].bb);
        double oldArea = m_nodes[left].bb.get_surface_area();
        double newArea = aabb.get_surface_area();
        costLeft = (newArea - oldArea) + inheritanceCost;
      }

      // Cost of descending to the right.
      double costRight;
      if (m_nodes[right].isLeaf()) {
        aabb aabb;
        aabb.merge(leafAABB, m_nodes[right].bb);
        costRight = aabb.get_surface_area() + inheritanceCost;
      }
      else {
        aabb aabb;
        aabb.merge(leafAABB, m_nodes[right].bb);
        double oldArea = m_nodes[right].bb.get_surface_area();
        double newArea = aabb.get_surface_area();
        costRight = (newArea - oldArea) + inheritanceCost;
      }

      // Descend according to the minimum cost.
      if ((cost < costLeft) && (cost < costRight))
        break;

      // Descend.
      if (costLeft < costRight)
        index = left;
      else
        index = right;
    }

    unsigned int sibling = index;

    // Create a new parent.
    unsigned int oldParent = m_nodes[sibling].parent;
    unsigned int newParent = allocate_node();
    m_nodes[newParent].parent = oldParent;
    m_nodes[newParent].bb.merge(leafAABB, m_nodes[sibling].bb);
    m_nodes[newParent].height = m_nodes[sibling].height + 1;

    // The sibling was not the root.
    if (oldParent != NULL_NODE) {
      if (m_nodes[oldParent].left == sibling)
        m_nodes[oldParent].left = newParent;
      else
        m_nodes[oldParent].right = newParent;

      m_nodes[newParent].left = sibling;
      m_nodes[newParent].right = leaf;
      m_nodes[sibling].parent = newParent;
      m_nodes[leaf].parent = newParent;
    }
    // The sibling was the root.
    else {
      m_nodes[newParent].left = sibling;
      m_nodes[newParent].right = leaf;
      m_nodes[sibling].parent = newParent;
      m_nodes[leaf].parent = newParent;
      m_root = newParent;
    }

    // Walk back up the tree fixing heights and AABBs.
    index = m_nodes[leaf].parent;
    while (index != NULL_NODE) {
      index = balance(index);

      unsigned int left = m_nodes[index].left;
      unsigned int right = m_nodes[index].right;

      assert(left != NULL_NODE);
      assert(right != NULL_NODE);

      m_nodes[index].height =
          1 + std::max(m_nodes[left].height, m_nodes[right].height);
      m_nodes[index].bb.merge(m_nodes[left].bb, m_nodes[right].bb);

      index = m_nodes[index].parent;
    }
  }

  //! Remove a leaf from the tree.
  /*! \param leaf
          The index of the leaf node.
   */
  void remove_leaf(unsigned int leaf)
  {
    --m_leaf_count;
    if (leaf == m_root) {
      m_root = NULL_NODE;
      return;
    }

    unsigned int parent = m_nodes[leaf].parent;
    unsigned int grandParent = m_nodes[parent].parent;
    unsigned int sibling;

    if (m_nodes[parent].left == leaf)
      sibling = m_nodes[parent].right;
    else
      sibling = m_nodes[parent].left;

    // Destroy the parent and connect the sibling to the grandparent.
    if (grandParent != NULL_NODE) {
      if (m_nodes[grandParent].left == parent)
        m_nodes[grandParent].left = sibling;
      else
        m_nodes[grandParent].right = sibling;

      m_nodes[sibling].parent = grandParent;
      free_node(parent);

      // Adjust ancestor bounds.
      unsigned int index = grandParent;
      while (index != NULL_NODE) {
        index = balance(index);

        unsigned int left = m_nodes[index].left;
        unsigned int right = m_nodes[index].right;

        m_nodes[index].bb.merge(m_nodes[left].bb, m_nodes[right].bb);
        m_nodes[index].height =
            1 + std::max(m_nodes[left].height, m_nodes[right].height);

        index = m_nodes[index].parent;
      }
    }
    else {
      m_root = sibling;
      m_nodes[sibling].parent = NULL_NODE;
      free_node(parent);
    }
  }

  //! Balance the tree.
  /*! \param leaf
          The index of the node.
   */
  unsigned int balance(unsigned int node)
  {
    assert(node != NULL_NODE);

    if (m_nodes[node].isLeaf() || (m_nodes[node].height < 2))
      return node;

    unsigned int left = m_nodes[node].left;
    unsigned int right = m_nodes[node].right;

    assert(left < m_node_capacity);
    assert(right < m_node_capacity);

    int currentBalance = m_nodes[right].height - m_nodes[left].height;

    // Rotate right branch up.
    if (currentBalance > 1) {
      unsigned int rightLeft = m_nodes[right].left;
      unsigned int rightRight = m_nodes[right].right;

      assert(rightLeft < m_node_capacity);
      assert(rightRight < m_node_capacity);

      // Swap node and its right-hand child.
      m_nodes[right].left = node;
      m_nodes[right].parent = m_nodes[node].parent;
      m_nodes[node].parent = right;

      // The node's old parent should now point to its right-hand child.
      if (m_nodes[right].parent != NULL_NODE) {
        if (m_nodes[m_nodes[right].parent].left == node)
          m_nodes[m_nodes[right].parent].left = right;
        else {
          assert(m_nodes[m_nodes[right].parent].right == node);
          m_nodes[m_nodes[right].parent].right = right;
        }
      }
      else
        m_root = right;

      // Rotate.
      if (m_nodes[rightLeft].height > m_nodes[rightRight].height) {
        m_nodes[right].right = rightLeft;
        m_nodes[node].right = rightRight;
        m_nodes[rightRight].parent = node;
        m_nodes[node].bb.merge(m_nodes[left].bb, m_nodes[rightRight].bb);
        m_nodes[right].bb.merge(m_nodes[node].bb, m_nodes[rightLeft].bb);

        m_nodes[node].height =
            1 + std::max(m_nodes[left].height, m_nodes[rightRight].height);
        m_nodes[right].height =
            1 + std::max(m_nodes[node].height, m_nodes[rightLeft].height);
      }
      else {
        m_nodes[right].right = rightRight;
        m_nodes[node].right = rightLeft;
        m_nodes[rightLeft].parent = node;
        m_nodes[node].bb.merge(m_nodes[left].bb, m_nodes[rightLeft].bb);
        m_nodes[right].bb.merge(m_nodes[node].bb, m_nodes[rightRight].bb);

        m_nodes[node].height =
            1 + std::max(m_nodes[left].height, m_nodes[rightLeft].height);
        m_nodes[right].height =
            1 + std::max(m_nodes[node].height, m_nodes[rightRight].height);
      }

      return right;
    }

    // Rotate left branch up.
    if (currentBalance < -1) {
      unsigned int leftLeft = m_nodes[left].left;
      unsigned int leftRight = m_nodes[left].right;

      assert(leftLeft < m_node_capacity);
      assert(leftRight < m_node_capacity);

      // Swap node and its left-hand child.
      m_nodes[left].left = node;
      m_nodes[left].parent = m_nodes[node].parent;
      m_nodes[node].parent = left;

      // The node's old parent should now point to its left-hand child.
      if (m_nodes[left].parent != NULL_NODE) {
        if (m_nodes[m_nodes[left].parent].left == node)
          m_nodes[m_nodes[left].parent].left = left;
        else {
          assert(m_nodes[m_nodes[left].parent].right == node);
          m_nodes[m_nodes[left].parent].right = left;
        }
      }
      else
        m_root = left;

      // Rotate.
      if (m_nodes[leftLeft].height > m_nodes[leftRight].height) {
        m_nodes[left].right = leftLeft;
        m_nodes[node].left = leftRight;
        m_nodes[leftRight].parent = node;
        m_nodes[node].bb.merge(m_nodes[right].bb, m_nodes[leftRight].bb);
        m_nodes[left].bb.merge(m_nodes[node].bb, m_nodes[leftLeft].bb);

        m_nodes[node].height =
            1 + std::max(m_nodes[right].height, m_nodes[leftRight].height);
        m_nodes[left].height =
            1 + std::max(m_nodes[node].height, m_nodes[leftLeft].height);
      }
      else {
        m_nodes[left].right = leftRight;
        m_nodes[node].left = leftLeft;
        m_nodes[leftLeft].parent = node;
        m_nodes[node].bb.merge(m_nodes[right].bb, m_nodes[leftLeft].bb);
        m_nodes[left].bb.merge(m_nodes[node].bb, m_nodes[leftRight].bb);

        m_nodes[node].height =
            1 + std::max(m_nodes[right].height, m_nodes[leftLeft].height);
        m_nodes[left].height =
            1 + std::max(m_nodes[node].height, m_nodes[leftRight].height);
      }

      return left;
    }

    return node;
  }

  //! Compute the height of the tree.
  /*! \return
          The height of the entire tree.
   */
  unsigned int compute_height() const { return compute_height(m_root); }

  //! Compute the height of a sub-tree.
  /*! \param node
          The index of the root node.

      \return
          The height of the sub-tree.
   */
  unsigned int compute_height(unsigned int node) const
  {
    assert(node < m_node_capacity);

    if (m_nodes[node].isLeaf())
      return 0;

    unsigned int height1 = compute_height(m_nodes[node].left);
    unsigned int height2 = compute_height(m_nodes[node].right);

    return 1 + std::max(height1, height2);
  }

  //! Assert that the sub-tree has a valid structure.
  /*! \param node
          The index of the root node.
   */
  void validate_structure(unsigned int node) const
  {
    if (node == NULL_NODE)
      return;

    if (node == m_root)
      assert(m_nodes[node].parent == NULL_NODE);

    unsigned int left = m_nodes[node].left;
    unsigned int right = m_nodes[node].right;

    if (m_nodes[node].isLeaf()) {
      assert(left == NULL_NODE);
      assert(right == NULL_NODE);
      assert(m_nodes[node].height == 0);
      return;
    }

    assert(left < m_node_capacity);
    assert(right < m_node_capacity);

    assert(m_nodes[left].parent == node);
    assert(m_nodes[right].parent == node);

    validate_structure(left);
    validate_structure(right);
  }

  //! Assert that the sub-tree has valid metrics.
  /*! \param node
          The index of the root node.
   */
  void validate_metrics(unsigned int node) const
  {
    if (node == NULL_NODE)
      return;

    unsigned int left = m_nodes[node].left;
    unsigned int right = m_nodes[node].right;

    if (m_nodes[node].isLeaf()) {
      assert(left == NULL_NODE);
      assert(right == NULL_NODE);
      assert(m_nodes[node].height == 0);
      return;
    }

    assert(left < m_node_capacity);
    assert(right < m_node_capacity);

    int height1 = m_nodes[left].height;
    int height2 = m_nodes[right].height;
    int height = 1 + std::max(height1, height2);
    (void)height;  // Unused variable in Release build
    assert(m_nodes[node].height == height);

    aabb aabb;
    aabb.merge(m_nodes[left].bb, m_nodes[right].bb);

    for (unsigned int i = 0; i < Dim; i++) {
      assert(aabb.lowerBound[i] == m_nodes[node].bb.lowerBound[i]);
      assert(aabb.upperBound[i] == m_nodes[node].bb.upperBound[i]);
    }

    validate_metrics(left);
    validate_metrics(right);
  }

  //! Compute minimum image separation.
  /*! \param separation
          The separation vector.

      \param shift
          The shift vector.

      \return
          Whether a periodic shift has been applied.
   */
};

// Utility typedefs
#define TYPEDEFS(suffix, dim, type)       \
  using point##suffix = point<dim, type>; \
  using aabb##suffix = aabb<dim, type>;   \
  using tree##suffix = tree<dim, type>

TYPEDEFS(2d, 2, double);
TYPEDEFS(2f, 2, float);
TYPEDEFS(2i, 2, int);
TYPEDEFS(3d, 3, double);
TYPEDEFS(3f, 3, float);
TYPEDEFS(3i, 3, int);

#undef TYPEDEFS

}  // namespace abt

#endif /* _AABB_H */

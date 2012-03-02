/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   std
 *
 ************************************************************************/

#ifndef included_tbox_vector
#define included_tbox_vector

#include "SAMRAI/SAMRAI_config.h"
#include <vector>

namespace SAMRAI {
namespace tbox {

/*!
 * Drop-in replacement for std::vector, with some SAMRAI enhancements.
 *
 * When compiled with debug, this class adds array bounds checking to
 * std::vector.  This feature is intended to help catch memory bugs.
 *
 * This class introduces duplicate interface to every std::vector
 * interface.  Hopefully, this trivial wrapper will get optimized
 * out when you use compiler optimization.  However, to be absolutely
 * sure, you may want to switch back to std::vector just in case there
 * is a performance penalty for using this class.
 *
 */
template<class TYPE>
class vector
{

private:
   //! @brief Non-essential shorthand.
   typedef std::vector<TYPE> SVec;

public:
   typedef typename SVec::iterator iterator;
   typedef typename SVec::const_iterator const_iterator;
   typedef typename SVec::reverse_iterator reverse_iterator;
   typedef typename SVec::const_reverse_iterator const_reverse_iterator;
   typedef typename SVec::value_type value_type;
   typedef typename SVec::size_type size_type;
   typedef typename SVec::reference reference;
   typedef typename SVec::const_reference const_reference;

   vector():d_vec() {
   }
   vector(
      const vector& r):d_vec(r.d_vec) {
   }
   explicit vector(
      size_type n):d_vec(n) {
   }
   ~vector() {
   }

   //@{ @name Replication of std::vector methods.

   iterator
   begin() {
      return d_vec.begin();
   }
   iterator
   end() {
      return d_vec.end();
   }
   const_iterator
   begin() const {
      return d_vec.begin();
   }
   const_iterator
   end() const {
      return d_vec.end();
   }
   reverse_iterator
   rbegin() {
      return d_vec.rbegin();
   }
   reverse_iterator
   rend() {
      return d_vec.rend();
   }
   const_reverse_iterator
   rbegin() const {
      return d_vec.rbegin();
   }
   const_reverse_iterator
   rend() const {
      return d_vec.rend();
   }

   iterator
   insert(
      iterator i,
      const value_type& v) {
      return d_vec.insert(i, v);
   }
   void
   insert(
      iterator p,
      size_type n,
      const TYPE& x) {
      d_vec.insert(p, n, x);
   }
   template<class InputIterator>
   void
   insert(
      InputIterator i,
      InputIterator j) {
      d_vec.insert(i, j);
   }
   void
   erase(
      iterator i) {
      d_vec.erase(i);
   }
   size_type
   erase(
      iterator& p) {
      return d_vec.erase(p);
   }
   size_type
   erase(
      iterator& p,
      iterator& q) {
      return d_vec.erase(p, q);
   }

   size_t
   size() const {
      return d_vec.size();
   }
   bool
   empty() const {
      return d_vec.empty();
   }
   void
   clear() {
      d_vec.clear();
   }
   void
   resize(
      size_type n,
      const TYPE& t) {
      d_vec.resize(n, t);
   }
   void
   reserve(
      size_type n) {
      d_vec.reserve(n);
   }

   vector&
   operator = (
      const SVec& r) {
      d_vec = r;
   }
   bool
   operator == (
      const SVec& r) const {
      return d_vec == r;
   }
   bool
   operator != (
      const SVec& r) const {
      return d_vec != r;
   }
   void
   swap(
      SVec& r) {
      d_vec.swap(r.d_vec);
   }

   operator const SVec& () const {
      return d_vec;
   }
   operator SVec& () {
      return d_vec;
   }

   //@}

   //@{ @name Enhancements of std::vector methods.

   //! @brief Vector element accessor with array bound checking.
   TYPE&
   operator [] (
      size_type i) {
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT(i <= d_vec.size() - 1);
#endif
      return d_vec[i];
   }
   //! @brief Vector element accessor with array bound checking.
   const TYPE&
   operator [] (
      size_type i) const {
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT(i <= d_vec.size() - 1);
#endif
      return d_vec[i];
   }

   //@}

private:
   SVec d_vec;
};

}
}

#endif
